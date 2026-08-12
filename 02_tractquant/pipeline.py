"""
TractQuant — Anterograde tracing quantification pipeline
=========================================================
Phase 1: Per-section intensity extraction in construct space

Pipeline:
  1. Parse VisualAlign JSON → anchoring vectors per section
  2. Load Allen CCF annotation volume (25um)
  3. For each section: warp atlas labels onto full-res tif space
  4. Extract per-region integrated intensity (both channels)
  5. Compute modal background per section per channel
  6. Resolve channel → construct assignment per animal
  7. Normalize to injection site mean (per channel independently)
  8. Output CSV in construct space with full compound keys

Requirements:
    pip install numpy scipy tifffile nrrd pandas tqdm

Allen CCF annotation volume (25um):
    Download annotation_25.nrrd from:
    https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/annotation/ccf_2017/
"""

import json
import gzip
import struct
import numpy as np
import tifffile
import nrrd
import pandas as pd
from PIL import Image, ImageDraw
from pathlib import Path
from scipy import ndimage, stats
from scipy.interpolate import LinearNDInterpolator
from tqdm import tqdm
from dataclasses import dataclass, field, asdict
from typing import Optional


BACKGROUND_METHOD = "whole_image_modal"


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class ChannelMeta:
    """Stress metadata for one channel — the compound experimental key."""
    construct: str          # "TRAP-Cre" or "TRAP-tTA"
    stress_type: str        # free text e.g. "Restraint", "Social defeat", "None"
    stress_duration: str    # "Acute" | "Sub-chronic" | "Chronic" | "None"
    timepoint: str          # free text e.g. "Day 1", "Home cage"

    @property
    def compound_key(self) -> str:
        parts = [self.construct, self.stress_type, self.stress_duration, self.timepoint]
        return " · ".join(p for p in parts if p and p != "None")


@dataclass
class AnimalConfig:
    """All metadata for one animal."""
    animal_id: str
    sex: str                          # "Male" | "Female"
    ch1: ChannelMeta                  # Ch1 construct + stress
    # Ch2 is always the opposite construct — stress fields independent
    ch2_stress_type: str = ""
    ch2_stress_duration: str = "None"
    ch2_timepoint: str = ""
    notes: str = ""

    # File paths
    visual_align_json: Path = None
    ch1_image_dir: Path = None        # folder containing ch1_001.tif … ch1_040.tif
    ch2_image_dir: Path = None

    # Derived
    @property
    def ch2_construct(self) -> str:
        return "TRAP-tTA" if self.ch1.construct == "TRAP-Cre" else "TRAP-Cre"

    @property
    def ch2(self) -> ChannelMeta:
        return ChannelMeta(
            construct=self.ch2_construct,
            stress_type=self.ch2_stress_type,
            stress_duration=self.ch2_stress_duration,
            timepoint=self.ch2_timepoint,
        )

    def resolve_construct(self, channel: str) -> ChannelMeta:
        """Return ChannelMeta for 'ch1' or 'ch2'."""
        return self.ch1 if channel == "ch1" else self.ch2


# ---------------------------------------------------------------------------
# Atlas
# ---------------------------------------------------------------------------

class CCFAtlas:
    """
    Loads and queries the Allen CCF annotation volume.

    The annotation volume maps every voxel (x, y, z) to a region ID.
    Resolution: 25 µm/voxel.
    Axes: [ML, DV, AP] — matching VisualAlign convention.
    """

    def __init__(self, annotation_path: str, region_csv_path: Optional[str] = None):
        print(f"Loading CCF annotation: {annotation_path}")
        # Support both .nrrd (raw Allen) and .tiff/.tif (BrainGlobe format)
        ann_path = str(annotation_path)
        if ann_path.endswith((".nii.gz", ".nii")):
            self.volume = self._load_nifti_annotation(ann_path)
            self.axis_order = "VISUALIGN_ML_AP_DV"
            print(f"  Loaded NIfTI/Cutlas format (VisualAlign)")
        elif ann_path.endswith((".tiff", ".tif")):
            import tifffile
            self.volume = tifffile.imread(ann_path)
            self.axis_order = "AP_DV_ML"
            print(f"  Loaded TIFF format (BrainGlobe)")
        else:
            self.volume, header = nrrd.read(ann_path)
            self.axis_order = "ML_DV_AP"
            print(f"  Loaded NRRD format")
        self.shape = self.volume.shape
        print(f"  Atlas shape ({self.axis_order}): {self.shape}")

        # Region lookup: id → (acronym, name, parent_id, ...)
        # If a CSV isn't supplied we build a minimal lookup from unique IDs
        if region_csv_path:
            df = pd.read_csv(region_csv_path)
            self.region_lookup = {
                row["id"]: {"acronym": row["acronym"], "name": row["name"],
                             "parent_id": row.get("parent_id", None),
                             "division": row.get("division", "")}
                for _, row in df.iterrows()
            }
        else:
            # Minimal fallback: just unique IDs present in the volume
            unique_ids = np.unique(self.volume)
            self.region_lookup = {int(i): {"acronym": str(i), "name": str(i),
                                            "parent_id": None, "division": ""}
                                   for i in unique_ids}

    @staticmethod
    def _load_nifti_annotation(annotation_path: str) -> np.ndarray:
        """
        Load VisualAlign/Cutlas NIfTI labels without requiring nibabel.

        VisualAlign's ABA_Mouse_CCFv3_2017_25um.cutlas/labels.nii.gz stores
        labels in [ML, AP, DV] order. Keeping that order is important because
        the JSON anchoring vectors use the same coordinate convention.
        """
        opener = gzip.open if annotation_path.endswith(".gz") else open
        with opener(annotation_path, "rb") as f:
            header = f.read(352)
            if len(header) < 352:
                raise ValueError(f"NIfTI header too short: {annotation_path}")

            sizeof_hdr = struct.unpack("<i", header[0:4])[0]
            if sizeof_hdr != 348:
                raise ValueError(f"Unsupported NIfTI header in {annotation_path}")

            dims = struct.unpack("<8h", header[40:56])
            datatype = struct.unpack("<h", header[70:72])[0]
            vox_offset = int(struct.unpack("<f", header[108:112])[0])

            dtype_map = {
                2: np.uint8,
                4: np.int16,
                8: np.int32,
                16: np.float32,
                512: np.uint16,
                768: np.uint32,
            }
            if datatype not in dtype_map:
                raise ValueError(f"Unsupported NIfTI datatype {datatype} in {annotation_path}")

            shape = tuple(int(d) for d in dims[1:1 + dims[0]])
            if len(shape) != 3:
                raise ValueError(f"Expected 3D NIfTI atlas, got shape {shape}")

            f.seek(vox_offset)
            data = f.read()

        dtype = np.dtype(dtype_map[datatype]).newbyteorder("<")
        expected = int(np.prod(shape))
        volume = np.frombuffer(data, dtype=dtype, count=expected)
        if volume.size != expected:
            raise ValueError(
                f"NIfTI data size mismatch in {annotation_path}: "
                f"expected {expected} voxels, got {volume.size}"
            )
        return volume.reshape(shape, order="F")

    def get_region_info(self, region_id: int) -> dict:
        return self.region_lookup.get(int(region_id),
               {"acronym": str(region_id), "name": "Unknown",
                "parent_id": None, "division": ""})

    def query(self, ml: np.ndarray, dv: np.ndarray, ap: np.ndarray) -> np.ndarray:
        """
        Return region IDs for arrays of atlas coordinates.
        Coordinates outside the volume return 0 (background).
        """
        ml_i = np.round(ml).astype(int)
        dv_i = np.round(dv).astype(int)
        ap_i = np.round(ap).astype(int)

        if self.axis_order == "VISUALIGN_ML_AP_DV":
            a_i, b_i, c_i = ml_i, ap_i, dv_i
        elif self.axis_order == "AP_DV_ML":
            a_i, b_i, c_i = ap_i, dv_i, ml_i
        else:
            a_i, b_i, c_i = ml_i, dv_i, ap_i

        valid = (
            (a_i >= 0) & (a_i < self.shape[0]) &
            (b_i >= 0) & (b_i < self.shape[1]) &
            (c_i >= 0) & (c_i < self.shape[2])
        )

        result = np.zeros(ml_i.shape, dtype=np.int32)
        result[valid] = self.volume[a_i[valid], b_i[valid], c_i[valid]]
        return result


# ---------------------------------------------------------------------------
# VisualAlign JSON parsing
# ---------------------------------------------------------------------------

def parse_visual_align(json_path: Path) -> dict:
    """
    Parse a VisualAlign JSON file.

    Returns dict keyed by section filename with:
        nr          : section number (1-based)
        width       : image width in pixels (of the 25% PNG)
        height      : image height in pixels (of the 25% PNG)
        o           : atlas origin vector [ml, ap, dv]
        u           : atlas step vector in x [ml, ap, dv]
        v           : atlas step vector in y [ml, ap, dv]
    """
    with open(json_path) as f:
        raw = f.read()

    # Fix common VisualAlign export issue: leading zeros in nr field
    # e.g. "nr":000 or "nr":001 which are invalid JSON
    import re
    raw = re.sub(r'"nr"\s*:\s*0+(\d)', r'"nr":\1', raw)  # 001->1, 023->23
    raw = re.sub(r'"nr"\s*:\s*0+,',    r'"nr":0,',  raw)  # 000->0

    try:
        data = json.loads(raw)
    except json.JSONDecodeError as e:
        raise ValueError(
            f"Could not parse VisualAlign JSON: {json_path}\n"
            f"Error at line {e.lineno}, column {e.colno}: {e.msg}\n"
            f"Please check the file is a valid VisualAlign export."
        ) from e

    sections = {}
    for s in data["slices"]:
        anch = s["anchoring"]
        # VisualAlign anchoring: [o_ml, o_ap, o_dv, u_ml, u_ap, u_dv, v_ml, v_ap, v_dv]
        sections[s["filename"]] = {
            "nr":     int(s["nr"]),
            "width":  s["width"],
            "height": s["height"],
            "o": np.array(anch[0:3]),   # atlas coord of top-left pixel
            "u": np.array(anch[3:6]),   # atlas step per pixel rightward
            "v": np.array(anch[6:9]),   # atlas step per pixel downward
            "has_markers": "markers" in s,
            "markers": s.get("markers", []),
            "alignment_dir": str(Path(json_path).parent),
        }

    print(f"  Parsed {len(sections)} sections from {json_path.name}")
    return sections


# ---------------------------------------------------------------------------
# Core: warp atlas labels onto full-res image space
# ---------------------------------------------------------------------------

def apply_marker_warp(section: dict, x_va: np.ndarray, y_va: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Apply VisualAlign marker displacement in the 2D alignment image space.

    VisualAlign stores markers as paired 2D points. In the exports used here,
    the first point is the pre-warp/atlas-plane coordinate and the second point
    is the matching histology image coordinate. For every histology pixel we
    therefore interpolate the inverse displacement back into the atlas-plane
    coordinate system before querying the CCF volume.
    """
    markers = np.asarray(section.get("markers", []), dtype=np.float32)
    if markers.ndim != 2 or markers.shape[0] < 3 or markers.shape[1] < 4:
        return x_va, y_va

    src = markers[:, 0:2]
    dst = markers[:, 2:4]

    width = float(section.get("width", np.nanmax(x_va) + 1))
    height = float(section.get("height", np.nanmax(y_va) + 1))
    border = np.array([
        [0, 0], [width - 1, 0], [0, height - 1], [width - 1, height - 1],
        [width / 2, 0], [width / 2, height - 1], [0, height / 2], [width - 1, height / 2],
    ], dtype=np.float32)
    src = np.vstack([src, border])
    dst = np.vstack([dst, border])

    displacement = src - dst
    dx_interp = LinearNDInterpolator(dst, displacement[:, 0], fill_value=0.0)
    dy_interp = LinearNDInterpolator(dst, displacement[:, 1], fill_value=0.0)

    dx = dx_interp(x_va, y_va)
    dy = dy_interp(x_va, y_va)
    dx = np.nan_to_num(dx, nan=0.0)
    dy = np.nan_to_num(dy, nan=0.0)
    return x_va + dx, y_va + dy


def build_label_image(section: dict, full_res_width: int, full_res_height: int,
                      atlas: CCFAtlas, resize_factor: float = 0.25,
                      use_markers: bool = True) -> np.ndarray:
    """
    For each pixel in the full-resolution image, compute its atlas coordinate
    and return the CCF region label.

    The VisualAlign anchoring is calibrated to the resized PNG dimensions.
    We scale pixel coordinates by resize_factor before applying the transform,
    which is equivalent to scaling u and v vectors by resize_factor.

    Args:
        section         : parsed section dict from parse_visual_align()
        full_res_width  : width of the full-res tif
        full_res_height : height of the full-res tif
        atlas           : CCFAtlas instance
        resize_factor   : the Nutil resize fraction (default 0.25)
        use_markers     : apply VisualAlign marker deformation when available

    Returns:
        label_image : (full_res_height, full_res_width) int32 array of CCF region IDs
    """
    # Build pixel coordinate grids for full-res image
    px = np.arange(full_res_width,  dtype=np.float32)   # (W,)
    py = np.arange(full_res_height, dtype=np.float32)   # (H,)
    px_grid, py_grid = np.meshgrid(px, py)              # both (H, W)

    o = np.array(section["o"])
    u = np.array(section["u"])
    v = np.array(section["v"])

    # VisualAlign stores the anchoring in the coordinate system of the image
    # used for alignment. If the measured channel image is a scaled PNG/JPG/TIF
    # of the same field of view, first map its pixel grid back into that
    # VisualAlign grid. The fluorescence image itself is never rescaled.
    va_width = float(section.get("width", full_res_width))
    va_height = float(section.get("height", full_res_height))
    x_va = px_grid * (va_width / full_res_width)
    y_va = py_grid * (va_height / full_res_height)
    if use_markers:
        x_va, y_va = apply_marker_warp(section, x_va, y_va)

    u_per_va_px = u / va_width
    v_per_va_px = v / va_height

    # Atlas coordinates for every pixel in QuickNII/VisualAlign space:
    # c0 = mediolateral, c1 = anterior-posterior, c2 = dorsoventral.
    c_ml = o[0] + x_va * u_per_va_px[0] + y_va * v_per_va_px[0]
    c_ap = o[1] + x_va * u_per_va_px[1] + y_va * v_per_va_px[1]
    c_dv = o[2] + x_va * u_per_va_px[2] + y_va * v_per_va_px[2]

    label_image = atlas.query(c_ml, c_dv, c_ap)   # (H, W) int32
    return label_image


def save_alignment_overlay(
    section: dict,
    ch1_image: np.ndarray,
    ch2_image: np.ndarray,
    atlas: CCFAtlas,
    output_dir: Path,
    animal_id: str,
    section_nr: int,
    section_file: str | None = None,
    resize_factor: float = 0.25,
    use_markers: bool = True,
    alpha: float = 0.42,
    max_size: int = 1600,
) -> dict:
    """
    Save QC images showing the exact atlas mask used for measurement.

    Output:
        alignment_qc/sXXX_<filename>_ch1_overlay.png
        alignment_qc/sXXX_<filename>_ch2_overlay.png
        alignment_qc/sXXX_<filename>_side_by_side.png

    The colored overlay is generated from build_label_image(...), so it shows
    the same aligned CCF labels that measure_section uses for quantification.
    If VisualAlign markers are present, the marker deformation is included.
    """
    out_dir = Path(output_dir) / animal_id / "alignment_qc"
    out_dir.mkdir(parents=True, exist_ok=True)

    def safe_stem(value: str | None) -> str:
        stem = Path(value or f"s{int(section_nr):03d}").stem
        keep = []
        for ch in stem:
            keep.append(ch if ch.isalnum() or ch in ("-", "_") else "_")
        return "".join(keep).strip("_") or f"s{int(section_nr):03d}"

    stem = f"s{int(section_nr):03d}_{safe_stem(section_file)}"

    def normalize_to_uint8(image: np.ndarray) -> np.ndarray:
        arr = np.asarray(image, dtype=np.float32)
        arr = np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)
        lo, hi = np.percentile(arr, [1, 99.5])
        if hi <= lo:
            lo, hi = float(arr.min()), float(arr.max())
        if hi <= lo:
            return np.zeros(arr.shape, dtype=np.uint8)
        arr = np.clip((arr - lo) / (hi - lo), 0, 1)
        return (arr * 255).astype(np.uint8)

    def label_colors(label_img: np.ndarray) -> np.ndarray:
        labels = label_img.astype(np.uint32)
        colors = np.zeros((*labels.shape, 3), dtype=np.uint8)
        colors[..., 0] = ((labels * 37 + 61) % 255).astype(np.uint8)
        colors[..., 1] = ((labels * 67 + 103) % 255).astype(np.uint8)
        colors[..., 2] = ((labels * 97 + 151) % 255).astype(np.uint8)
        return colors

    def label_boundary(label_img: np.ndarray) -> np.ndarray:
        mask = label_img > 0
        boundary = np.zeros(label_img.shape, dtype=bool)
        boundary[1:, :] |= label_img[1:, :] != label_img[:-1, :]
        boundary[:, 1:] |= label_img[:, 1:] != label_img[:, :-1]
        return boundary & mask

    def make_overlay(image: np.ndarray, channel_name: str) -> tuple[Image.Image, float]:
        height, width = image.shape
        labels = build_label_image(section, width, height, atlas, resize_factor, use_markers=use_markers)
        mask = labels > 0
        atlas_percent = float(mask.mean() * 100.0) if mask.size else 0.0
        reference_iou = None
        ref_path = None
        if section_file and section.get("alignment_dir"):
            candidate = Path(section["alignment_dir"]) / f"{Path(section_file).stem}_nl.png"
            if candidate.exists():
                ref = Image.open(candidate).resize((width, height), Image.Resampling.NEAREST).convert("RGB")
                ref_arr = np.asarray(ref)
                ref_mask = ref_arr.max(axis=2) > 8
                union = np.logical_or(mask, ref_mask).sum()
                if union:
                    reference_iou = float(np.logical_and(mask, ref_mask).sum() / union)
                    ref_path = str(candidate)

        base = np.dstack([normalize_to_uint8(image)] * 3).astype(np.float32)
        colors = label_colors(labels).astype(np.float32)
        overlay = base.copy()
        overlay[mask] = (1.0 - alpha) * overlay[mask] + alpha * colors[mask]
        overlay[label_boundary(labels)] = np.array([255, 220, 0], dtype=np.float32)
        overlay = np.clip(overlay, 0, 255).astype(np.uint8)

        pil = Image.fromarray(overlay, mode="RGB")
        scale = min(1.0, max_size / max(pil.size))
        if scale < 1.0:
            new_size = (max(1, int(pil.width * scale)), max(1, int(pil.height * scale)))
            pil = pil.resize(new_size, Image.Resampling.LANCZOS)

        draw = ImageDraw.Draw(pil)
        marker_text = "markers" if use_markers and section.get("markers") else "anchoring"
        text = f"{animal_id} s{int(section_nr):03d} {channel_name} | {marker_text} mask {atlas_percent:.1f}%"
        if reference_iou is not None:
            text += f" | ref IoU {reference_iou:.2f}"
        draw.rectangle((8, 8, min(pil.width - 8, 8 + len(text) * 7), 30), fill=(0, 0, 0))
        draw.text((12, 12), text, fill=(255, 255, 255))
        pil.info["reference_iou"] = reference_iou
        pil.info["reference_path"] = ref_path
        return pil, atlas_percent

    ch1_overlay, ch1_percent = make_overlay(ch1_image, "Ch1")
    ch2_overlay, ch2_percent = make_overlay(ch2_image, "Ch2")

    ch1_path = out_dir / f"{stem}_ch1_overlay.png"
    ch2_path = out_dir / f"{stem}_ch2_overlay.png"
    side_path = out_dir / f"{stem}_side_by_side.png"
    ch1_overlay.save(ch1_path)
    ch2_overlay.save(ch2_path)

    pad = 12
    side_w = ch1_overlay.width + ch2_overlay.width + pad
    side_h = max(ch1_overlay.height, ch2_overlay.height)
    side = Image.new("RGB", (side_w, side_h), (20, 20, 24))
    side.paste(ch1_overlay, (0, 0))
    side.paste(ch2_overlay, (ch1_overlay.width + pad, 0))
    side.save(side_path)

    return {
        "ch1_overlay": str(ch1_path),
        "ch2_overlay": str(ch2_path),
        "side_by_side": str(side_path),
        "ch1_atlas_percent": ch1_percent,
        "ch2_atlas_percent": ch2_percent,
        "ch1_reference_iou": ch1_overlay.info.get("reference_iou"),
        "ch2_reference_iou": ch2_overlay.info.get("reference_iou"),
        "reference_path": ch1_overlay.info.get("reference_path") or ch2_overlay.info.get("reference_path"),
        "use_markers": bool(use_markers and section.get("markers")),
        "marker_count": len(section.get("markers", [])),
    }


# ---------------------------------------------------------------------------
# Background estimation
# ---------------------------------------------------------------------------

def modal_background(image: np.ndarray, n_bins: int = 256,
                     tissue_mask: np.ndarray = None) -> float:
    """
    Estimate background as the modal intensity of a section.

    If tissue_mask is provided, black pixels outside the atlas-labelled tissue
    are excluded. If no mask is provided, the whole image is used.
    """
    if tissue_mask is not None and tissue_mask.sum() > 100:
        pixels = image[tissue_mask]
    else:
        pixels = image.ravel()

    # Use lower 75% of intensity range to find background mode
    upper = np.percentile(pixels, 75)
    lower_pixels = pixels[pixels <= upper]
    if lower_pixels.size < 100:
        return float(np.median(pixels))
    hist, edges = np.histogram(lower_pixels, bins=n_bins)
    mode_bin = np.argmax(hist)
    return float((edges[mode_bin] + edges[mode_bin + 1]) / 2)


# ---------------------------------------------------------------------------
# Per-section measurement
# ---------------------------------------------------------------------------

def measure_section(
    section: dict,
    ch1_image: np.ndarray,
    ch2_image: np.ndarray,
    atlas: CCFAtlas,
    resize_factor: float = 0.25,
    min_region_pixels: int = 10,
    use_markers: bool = True,
) -> pd.DataFrame:
    """
    For one section: build label image, compute per-region intensity stats.

    Returns DataFrame with columns:
        region_id, region_acronym, region_name, division,
        ch1_sum, ch1_mean, ch1_n_pixels,
        ch2_sum, ch2_mean, ch2_n_pixels,
        ch1_background, ch2_background,
        ch1_corrected_sum, ch2_corrected_sum
    """
    H1, W1 = ch1_image.shape
    H2, W2 = ch2_image.shape
    marker_count = len(section.get("markers", []))
    markers_applied = bool(use_markers and marker_count)

    # Build atlas masks in each channel's native pixel grid. This allows PNG,
    # JPG, and TIFF inputs, even mixed across channels, without rescaling image
    # intensities. Only the atlas mask is sampled at the loaded image size.
    label_img_ch1 = build_label_image(section, W1, H1, atlas, resize_factor, use_markers=use_markers)
    label_img_ch2 = label_img_ch1 if ch2_image.shape == ch1_image.shape else build_label_image(
        section, W2, H2, atlas, resize_factor, use_markers=use_markers
    )

    # Diagnostic: check what values the label image contains
    n_nonzero_ch1 = np.sum(label_img_ch1 > 0)
    n_nonzero_ch2 = np.sum(label_img_ch2 > 0)
    if n_nonzero_ch1 == 0 or n_nonzero_ch2 == 0:
        # Log useful debug info
        import sys
        print(f"  DEBUG: label image all zeros. "
              f"Ch1 shape={label_img_ch1.shape}, Ch2 shape={label_img_ch2.shape}. "
              f"Resize={resize_factor}. "
              f"o={section['o'][:3]}, "
              f"Atlas shape={atlas.shape}",
              file=sys.stderr)

    # Whole-section comparison method:
    # one modal background value per section and channel, measured over the
    # whole section image.
    ch1_bg_whole = modal_background(ch1_image)
    ch2_bg_whole = modal_background(ch2_image)
    ch1_bg_tissue = modal_background(ch1_image, tissue_mask=label_img_ch1 > 0)
    ch2_bg_tissue = modal_background(ch2_image, tissue_mask=label_img_ch2 > 0)
    ch1_bg = ch1_bg_whole
    ch2_bg = ch2_bg_whole

    unique_labels = np.union1d(np.unique(label_img_ch1), np.unique(label_img_ch2))
    unique_labels = unique_labels[unique_labels > 0]

    rows = []
    for rid in unique_labels:
        mask_ch1 = label_img_ch1 == rid
        mask_ch2 = label_img_ch2 == rid
        ch1_n_pix = int(mask_ch1.sum())
        ch2_n_pix = int(mask_ch2.sum())
        if ch1_n_pix < min_region_pixels and ch2_n_pix < min_region_pixels:
            continue

        info = atlas.get_region_info(int(rid))

        ch1_pixels = ch1_image[mask_ch1].astype(np.float64) if ch1_n_pix else np.array([], dtype=np.float64)
        ch2_pixels = ch2_image[mask_ch2].astype(np.float64) if ch2_n_pix else np.array([], dtype=np.float64)

        ch1_sum = float(ch1_pixels.sum())
        ch2_sum = float(ch2_pixels.sum())

        # CTCF per section per region = Integrated Density - (Area x background)
        ch1_corr = float(ch1_sum - ch1_bg * ch1_n_pix)
        ch2_corr = float(ch2_sum - ch2_bg * ch2_n_pix)

        rows.append({
            "region_id":      int(rid),
            "region_acronym": info["acronym"],
            "region_name":    info["name"],
            "division":       info["division"],
            "alignment_uses_markers": markers_applied,
            "alignment_marker_count": marker_count,
            "n_pixels":       ch1_n_pix,
            "ch1_n_pixels":   ch1_n_pix,
            "ch2_n_pixels":   ch2_n_pix,
            "ch1_sum":        ch1_sum,
            "ch2_sum":        ch2_sum,
            "background_method": BACKGROUND_METHOD,
            "ch1_background": ch1_bg,
            "ch2_background": ch2_bg,
            "ch1_background_whole_image": ch1_bg_whole,
            "ch2_background_whole_image": ch2_bg_whole,
            "ch1_background_tissue_only": ch1_bg_tissue,
            "ch2_background_tissue_only": ch2_bg_tissue,
            "ch1_mean":            ch1_sum / max(ch1_n_pix, 1),
            "ch2_mean":            ch2_sum / max(ch2_n_pix, 1),
            "ch1_corrected_sum":   max(ch1_corr, 0.0),
            "ch2_corrected_sum":   max(ch2_corr, 0.0),
        })

    df = pd.DataFrame(rows)
    df.attrs["use_markers"] = markers_applied
    df.attrs["marker_count"] = marker_count
    df.attrs["atlas_match_percent"] = float(n_nonzero_ch1 / label_img_ch1.size * 100.0) if label_img_ch1.size else 0.0
    return df


# ---------------------------------------------------------------------------
# Image loading utilities
# ---------------------------------------------------------------------------

def load_section_images(
    ch1_dir: Path,
    ch2_dir: Path,
    section_nr: int,
    filename: str | None = None,
):
    """
    Load Ch1 and Ch2 image files for a given section number.

    Flexible filename matching — handles any naming convention:
      - ch1_001.png, ch1_010.tif   (zero-padded with prefix)
      - IPAC2s010.jpg, F5s001.tif  (VisualAlign-style with animal prefix)
      - s010.tif, s01.tif          (section number only)

    Matches image files whose name ends with the section number
    (with any zero-padding).
    """
    nr = int(section_nr)
    section_tokens = []
    import re
    if filename:
        stem = Path(filename).stem
        section_tokens.append(stem)
        match = re.search(r"(?:^|[^A-Za-z0-9])s?0*(\d+)(?:[^0-9]*$)", stem, re.IGNORECASE)
        if match:
            section_tokens.extend([
                f"s{int(match.group(1)):03d}",
                f"{int(match.group(1)):03d}",
                str(int(match.group(1))),
            ])

    def find_image(folder: Path) -> Path:
        folder = Path(folder)
        extensions = ["*.png", "*.jpg", "*.jpeg", "*.tif", "*.tiff"]
        candidates = [
            f
            for ext in extensions
            for f in folder.glob(ext)
            if not f.name.startswith("._")
        ]
        if not candidates:
            raise FileNotFoundError(f"No image files found in {folder}")

        # Prefer the VisualAlign filename/section token. Some VisualAlign
        # exports use nr as slice order while the real section id lives in the
        # filename, e.g. nr=1 but filename="...s003.png".
        def token_match_score(path: Path, token: str) -> int | None:
            stem_lower = path.stem.lower()
            token_lower = token.lower()

            if stem_lower == token_lower:
                return 0
            if stem_lower.endswith(token_lower):
                return 1

            if token_lower.isdigit():
                n = int(token_lower)
                if re.search(rf"(^|[^A-Za-z0-9])s0*{n}([^0-9]|$)", stem_lower):
                    return 2
                if re.search(rf"(^|[^A-Za-z0-9])serie_0*{n}([^0-9]|$)", stem_lower):
                    return 3
                if re.search(rf"(^|[^0-9])0*{n}([^0-9]|$)", stem_lower):
                    return 9
                return None

            s_match = re.fullmatch(r"s0*(\d+)", token_lower)
            if s_match:
                n = int(s_match.group(1))
                if re.search(rf"(^|[^A-Za-z0-9])s0*{n}([^0-9]|$)", stem_lower):
                    return 2
                return None

            if token_lower in stem_lower:
                return 8
            return None

        for token in section_tokens:
            token_lower = token.lower()
            scored = [
                (token_match_score(f, token_lower), f)
                for f in candidates
            ]
            matches = [(score, f) for score, f in scored if score is not None]
            if matches:
                return sorted(matches, key=lambda item: (item[0], item[1].name))[0][1]

        # Match files ending with the section number (any zero-padding)
        for pad in [3, 2, 1, 4]:
            pattern = f"{nr:0{pad}d}"
            matches = [f for f in candidates
                      if re.search(rf"[^0-9]{pattern}\.tiff?$|^{pattern}\.tiff?$",
                                   f.name, re.IGNORECASE)
                      or f.stem.endswith(pattern)]
            if matches:
                return matches[0]

        # Last resort: sort all files and pick by bounded numeric match
        sorted_files = sorted(candidates)
        matches = [
            f for f in sorted_files
            if re.search(rf"(^|[^A-Za-z0-9])serie_0*{nr}([^0-9]|$)", f.stem.lower())
            or re.search(rf"(^|[^A-Za-z0-9])s0*{nr}([^0-9]|$)", f.stem.lower())
            or re.search(rf"(^|[^0-9])0*{nr}([^0-9]|$)", f.stem.lower())
        ]
        if matches:
            return matches[0]

        raise FileNotFoundError(
            f"No image found for section {nr} in {folder}\n"
            f"Available files: {[f.name for f in sorted_files[:5]]}"
        )

    ch1_path = find_image(Path(ch1_dir))
    ch2_path = find_image(Path(ch2_dir))

    def read_image(path: Path) -> np.ndarray:
        if path.suffix.lower() in {".tif", ".tiff"}:
            return tifffile.imread(str(path)).astype(np.float32)
        return np.asarray(Image.open(path).convert("F"), dtype=np.float32)

    ch1_img = read_image(ch1_path)
    ch2_img = read_image(ch2_path)

    # Handle multichannel tifs — take first channel
    if ch1_img.ndim == 3:
        ch1_img = ch1_img[0] if ch1_img.shape[0] < ch1_img.shape[-1] else ch1_img[:, :, 0]
    if ch2_img.ndim == 3:
        ch2_img = ch2_img[0] if ch2_img.shape[0] < ch2_img.shape[-1] else ch2_img[:, :, 0]

    return ch1_img, ch2_img


# ---------------------------------------------------------------------------
# Animal-level pipeline
# ---------------------------------------------------------------------------

def process_animal(
    animal: AnimalConfig,
    atlas: CCFAtlas,
    injection_regions: list[str],
    resize_factor: float = 0.25,
    min_sections_per_region: int = 1,
    sparse_flag_threshold: int = 3,
    output_dir: Path = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Full pipeline for one animal.

    Returns:
        section_df  : raw per-section measurements
        animal_df   : aggregated animal-level projections in construct space
    """
    print(f"\nProcessing animal: {animal.animal_id}")

    # Parse VisualAlign JSON
    sections = parse_visual_align(animal.visual_align_json)

    # Match JSON filenames to section numbers
    # JSON filenames are PNGs; we find the corresponding tif by section number
    all_section_rows = []

    for filename, sec in tqdm(sections.items(), desc=f"  Sections", leave=False):
        nr = sec["nr"]
        if not sec.get("has_markers", True):
            print(f"  ⚠ s{nr:03d}: no markers — anchoring only, included")

        try:
            ch1_img, ch2_img = load_section_images(
                animal.ch1_image_dir, animal.ch2_image_dir, nr, filename
            )
        except FileNotFoundError as e:
            print(f"  ⚠ Skipping s{nr:03d}: {e}")
            continue

        sec_df = measure_section(sec, ch1_img, ch2_img, atlas, resize_factor)
        if output_dir:
            qc = save_alignment_overlay(
                sec, ch1_img, ch2_img, atlas, output_dir,
                animal.animal_id, nr, filename, resize_factor
            )
            marker_status = (
                f"{qc['marker_count']} markers" if qc.get("use_markers")
                else "anchoring only"
            )
            print(
                f"  QC overlay s{nr:03d}: "
                f"{qc['ch1_atlas_percent']:.1f}%/{qc['ch2_atlas_percent']:.1f}% atlas mask, "
                f"{marker_status}"
            )
        sec_df["section_nr"]  = nr
        sec_df["section_file"] = filename
        all_section_rows.append(sec_df)

    if not all_section_rows:
        raise RuntimeError(f"No sections processed for {animal.animal_id}")

    section_df = pd.concat(all_section_rows, ignore_index=True)

    # -----------------------------------------------------------------------
    # Injection site baseline — mean corrected intensity per channel
    # across all sections containing the injection regions
    # -----------------------------------------------------------------------
    inj_mask = section_df["region_acronym"].isin(injection_regions)
    inj_sections = section_df[inj_mask]

    if inj_sections.empty:
        raise ValueError(
            f"No injection site regions found for {animal.animal_id}. "
            f"Check injection_regions: {injection_regions}"
        )

    # Mean CTCF per section — matches ImageJ workflow:
    #   Per section: CTCF = Integrated Density - (Area x section background)
    #   Across sections: mean CTCF per section the injection region appears in
    # This is scale-invariant — large/small injection regions treated equally.
    ch1_baseline = inj_sections.groupby("section_nr")["ch1_corrected_sum"].sum().mean()
    ch2_baseline = inj_sections.groupby("section_nr")["ch2_corrected_sum"].sum().mean()

    print(f"  Injection baseline — ch1: {ch1_baseline:.1f}, ch2: {ch2_baseline:.1f}")
    print(f"  Ch1={animal.ch1.construct} [{animal.ch1.compound_key}]")
    print(f"  Ch2={animal.ch2.construct} [{animal.ch2.compound_key}]")

    # -----------------------------------------------------------------------
    # Aggregate across sections per region
    # -----------------------------------------------------------------------
    agg = section_df.groupby(
        ["region_id", "region_acronym", "region_name", "division"]
    ).agg(
        ch1_corrected_sum=("ch1_corrected_sum", "sum"),
        ch2_corrected_sum=("ch2_corrected_sum", "sum"),
        total_pixels=("n_pixels", "sum"),
        n_sections=("section_nr", "nunique"),
    ).reset_index()

    # Projection fractions
    # Mean CTCF per section / injection mean CTCF per section
    # Per region: average the CTCF across sections it appears in
    agg["ch1_mean_ctcf"] = (section_df[section_df["region_id"].isin(agg["region_id"])]
        .groupby(["region_id","section_nr"])["ch1_corrected_sum"].sum()
        .groupby("region_id").mean()
        .reindex(agg["region_id"].values)
        .values)
    agg["ch2_mean_ctcf"] = (section_df[section_df["region_id"].isin(agg["region_id"])]
        .groupby(["region_id","section_nr"])["ch2_corrected_sum"].sum()
        .groupby("region_id").mean()
        .reindex(agg["region_id"].values)
        .values)
    agg["ch1_fraction"] = agg["ch1_mean_ctcf"] / ch1_baseline
    agg["ch2_fraction"] = agg["ch2_mean_ctcf"] / ch2_baseline
    agg["ch1_ch2_ratio"] = np.where(
        agg["ch2_fraction"] > 0,
        agg["ch1_fraction"] / agg["ch2_fraction"],
        np.nan
    )

    # Sparse flag
    agg["sparse_flag"] = agg["n_sections"] < sparse_flag_threshold

    # -----------------------------------------------------------------------
    # Resolve to construct space
    # Ch1 might be TRAP-Cre or TRAP-tTA depending on counterbalancing
    # -----------------------------------------------------------------------
    ch1_meta = animal.ch1
    ch2_meta = animal.ch2

    if ch1_meta.construct == "TRAP-Cre":
        trapcre_frac = agg["ch1_fraction"]
        traptta_frac = agg["ch2_fraction"]
        trapcre_baseline = ch1_baseline
        traptta_baseline = ch2_baseline
    else:
        trapcre_frac = agg["ch2_fraction"]
        traptta_frac = agg["ch1_fraction"]
        trapcre_baseline = ch2_baseline
        traptta_baseline = ch1_baseline
        ch1_meta, ch2_meta = ch2_meta, ch1_meta  # swap for metadata output

    cre_meta  = animal.ch1 if animal.ch1.construct == "TRAP-Cre" else animal.ch2
    tta_meta  = animal.ch1 if animal.ch1.construct == "TRAP-tTA" else animal.ch2

    animal_df = agg[["region_id","region_acronym","region_name","division",
                      "total_pixels","n_sections","sparse_flag"]].copy()

    animal_df["animal_id"]             = animal.animal_id
    animal_df["sex"]                   = animal.sex
    animal_df["ch1_construct"]         = animal.ch1.construct
    # TRAP-Cre construct metadata
    animal_df["trapcre_stress_type"]   = cre_meta.stress_type
    animal_df["trapcre_stress_dur"]    = cre_meta.stress_duration
    animal_df["trapcre_timepoint"]     = cre_meta.timepoint
    animal_df["trapcre_compound_key"]  = cre_meta.compound_key
    animal_df["trapcre_baseline"]      = trapcre_baseline
    animal_df["trapcre_fraction"]      = trapcre_frac.values
    # TRAP-tTA construct metadata
    animal_df["traptta_stress_type"]   = tta_meta.stress_type
    animal_df["traptta_stress_dur"]    = tta_meta.stress_duration
    animal_df["traptta_timepoint"]     = tta_meta.timepoint
    animal_df["traptta_compound_key"]  = tta_meta.compound_key
    animal_df["traptta_baseline"]      = traptta_baseline
    animal_df["traptta_fraction"]      = traptta_frac.values
    # Ratio always in construct space
    animal_df["trapcre_traptta_ratio"] = np.where(
        animal_df["traptta_fraction"] > 0,
        animal_df["trapcre_fraction"] / animal_df["traptta_fraction"],
        np.nan
    )

    # Save per-animal section-level data
    if output_dir:
        out = Path(output_dir) / animal.animal_id
        out.mkdir(parents=True, exist_ok=True)
        section_df.to_csv(out / "section_measurements.csv", index=False)
        animal_df.to_csv(out / "animal_projections.csv", index=False)
        print(f"  Saved to {out}/")

    return section_df, animal_df


# ---------------------------------------------------------------------------
# Cohort aggregation
# ---------------------------------------------------------------------------

def aggregate_cohort(
    animal_dfs: list[pd.DataFrame],
    group_by: list[str] = None,
) -> pd.DataFrame:
    """
    Aggregate projection fractions across animals.

    Default grouping: region × trapcre_compound_key × traptta_compound_key × sex
    Returns mean, SEM, z-score, and n per group per region.
    """
    if group_by is None:
        group_by = [
            "region_acronym", "region_name", "division",
            "trapcre_compound_key", "traptta_compound_key", "sex"
        ]

    all_df = pd.concat(animal_dfs, ignore_index=True)

    # Handle both column name variants
    ratio_col = "trapcre_traptta_ratio" if "trapcre_traptta_ratio" in all_df.columns else "ratio"

    def zscore_col(s):
        if s.std() == 0 or len(s) < 2:
            return pd.Series(np.zeros(len(s)), index=s.index)
        return (s - s.mean()) / s.std()

    cohort = all_df.groupby(group_by).agg(
        trapcre_mean=("trapcre_fraction", "mean"),
        trapcre_sem=("trapcre_fraction", lambda x: x.sem()),
        traptta_mean=("traptta_fraction", "mean"),
        traptta_sem=("traptta_fraction", lambda x: x.sem()),
        ratio_mean=(ratio_col, "mean"),
        n_animals=("animal_id", "nunique"),
        n_sections_mean=("n_sections", "mean"),
    ).reset_index()

    # Z-scores within each region across groups
    for col, zcol in [("trapcre_mean","trapcre_zscore"),("traptta_mean","traptta_zscore")]:
        cohort[zcol] = cohort.groupby("region_acronym")[col].transform(zscore_col)

    # Co-projection flag: both constructs above threshold in same region
    copro_threshold = 0.05
    cohort["coprojection"] = (
        (cohort["trapcre_mean"] > copro_threshold) &
        (cohort["traptta_mean"] > copro_threshold)
    )

    return cohort


# ---------------------------------------------------------------------------
# BrainGlobe export
# ---------------------------------------------------------------------------

def export_brainglobe(
    animal_df: pd.DataFrame,
    output_path: Path,
    value_column: str = "trapcre_fraction",
    min_fraction: float = 0.0,
):
    """
    Export a {region_acronym: value} JSON dict compatible with
    brainglobe-heatmap and brainrender.

    Usage in brainglobe:
        from brainglobe_heatmap import Heatmap
        import json
        values = json.load(open("trapcre_heatmap.json"))
        h = Heatmap(values, position=5000, orientation="frontal")
        h.show()
    """
    data = (
        animal_df[animal_df[value_column] > min_fraction]
        .set_index("region_acronym")[value_column]
        .to_dict()
    )
    import json
    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)
    print(f"BrainGlobe export: {output_path} ({len(data)} regions)")


# ---------------------------------------------------------------------------
# QC: atlas overlay image
# ---------------------------------------------------------------------------

def save_qc_overlay(
    ch1_image: np.ndarray,
    label_image: np.ndarray,
    output_path: Path,
    alpha: float = 0.35,
):
    """
    Save a QC image with atlas region boundaries overlaid on Ch1 fluorescence.
    Requires matplotlib.
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors
    except ImportError:
        print("matplotlib not installed — skipping QC overlay")
        return

    # Edge detection on label image to get region boundaries
    from scipy.ndimage import sobel
    sx = sobel(label_image.astype(float), axis=0)
    sy = sobel(label_image.astype(float), axis=1)
    edges = np.hypot(sx, sy) > 0

    fig, ax = plt.subplots(figsize=(8, 8))
    # Normalize fluorescence for display
    vmin, vmax = np.percentile(ch1_image, [1, 99])
    ax.imshow(ch1_image, cmap="gray", vmin=vmin, vmax=vmax)
    # Overlay edges in red
    overlay = np.zeros((*label_image.shape, 4), dtype=np.float32)
    overlay[edges, 0] = 1.0   # red
    overlay[edges, 3] = alpha
    ax.imshow(overlay)
    ax.axis("off")
    ax.set_title("Atlas overlay — Ch1 fluorescence + CCF boundaries", fontsize=10)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()


# ---------------------------------------------------------------------------
# QC: background diagnostic
# ---------------------------------------------------------------------------

def diagnose_background_method(
    image_path: str | Path,
    visual_align_json: str | Path,
    section_nr: int,
    atlas_path: str = None,
    region_csv_path: str = None,
    atlas_cache_dir: Path = None,
    resize_factor: float = 0.25,
    use_markers: bool = True,
) -> dict:
    """
    Compare legacy whole-image background with atlas tissue-only background.

    Use this when checking whether black border pixels are causing
    under-correction. A high tissue-only value compared with the whole-image
    value means the old modal background was probably pulled toward zero.
    """
    image_path = Path(image_path)
    visual_align_json = Path(visual_align_json)

    from atlas_manager import get_atlas_volume, get_region_csv
    if atlas_path is None:
        atlas_path = get_atlas_volume(atlas_cache_dir)
    if region_csv_path is None:
        region_csv_path = get_region_csv(atlas_cache_dir)

    def read_image(path: Path) -> np.ndarray:
        if path.suffix.lower() in {".tif", ".tiff"}:
            return tifffile.imread(str(path)).astype(np.float64)
        return np.asarray(Image.open(path).convert("F"), dtype=np.float64)

    image = read_image(image_path)
    if image.ndim == 3:
        image = image[0] if image.shape[0] < image.shape[-1] else image[:, :, 0]

    atlas = CCFAtlas(str(atlas_path), str(region_csv_path))
    sections = parse_visual_align(visual_align_json)

    section = None
    section_file = None
    for filename, candidate in sections.items():
        if candidate["nr"] == int(section_nr):
            section = candidate
            section_file = filename
            break
    if section is None:
        raise ValueError(f"Section {section_nr} not found in {visual_align_json}")

    height, width = image.shape
    label_img = build_label_image(
        section, width, height, atlas, resize_factor, use_markers=use_markers
    )
    tissue_mask = label_img > 0

    bg_whole = modal_background(image)
    bg_tissue = modal_background(image, tissue_mask=tissue_mask)

    result = {
        "image_path": str(image_path),
        "visual_align_json": str(visual_align_json),
        "section_nr": int(section_nr),
        "section_file": section_file,
        "width": int(width),
        "height": int(height),
        "total_pixels": int(width * height),
        "tissue_pixels": int(tissue_mask.sum()),
        "black_border_pixels": int((~tissue_mask).sum()),
        "tissue_percent": float(100 * tissue_mask.sum() / (width * height)),
        "background_whole_image": float(bg_whole),
        "background_tissue_only": float(bg_tissue),
        "background_difference": float(bg_tissue - bg_whole),
        "use_markers": bool(use_markers and section.get("markers")),
        "marker_count": len(section.get("markers", [])),
    }

    print(f"Image: {width} x {height} = {width * height:,} total pixels")
    print(
        "Tissue pixels (atlas label > 0): "
        f"{result['tissue_pixels']:,} ({result['tissue_percent']:.1f}%)"
    )
    print(f"Black border pixels: {result['black_border_pixels']:,}")
    print(f"Markers used: {result['use_markers']} ({result['marker_count']})")
    print()
    print(f"Modal background (whole image): {bg_whole:.2f}")
    print(f"Modal background (tissue only): {bg_tissue:.2f}")
    print(f"Difference: {bg_tissue - bg_whole:.2f}")

    return result


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run_pipeline(
    animals: list[AnimalConfig],
    injection_regions: list[str],
    output_dir: str,
    atlas_path: str = None,
    region_csv_path: str = None,
    atlas_cache_dir: Path = None,
    resize_factor: float = 0.25,
    sparse_flag_threshold: int = 3,
    save_qc: bool = True,
):
    """
    Run TractQuant for a list of animals.

    The Allen CCF atlas and region CSV are downloaded automatically
    on first run and cached in ~/.tractquant/atlas/ for future runs.
    Pass atlas_path and region_csv_path to override with local files.

    Args:
        animals             : list of AnimalConfig objects
        injection_regions   : list of CCF acronyms defining injection site
                              e.g. ["MOs", "MOp"]
        atlas_path          : (optional) local path to annotation_25.nrrd
                              if None, auto-downloaded on first run
        atlas_cache_dir     : (optional) override cache dir
                              default: ~/.tractquant/atlas/
        output_dir          : root output directory
        region_csv_path     : optional CSV with columns id, acronym, name, division
        resize_factor       : Nutil resize factor used (default 0.25)
        sparse_flag_threshold : regions in fewer sections than this are flagged
        save_qc             : whether to save atlas overlay images

    Output structure:
        output_dir/
            {animal_id}/
                section_measurements.csv
                animal_projections.csv
                qc_overlays/
                    s001_overlay.png
                    ...
            cohort_summary.csv
            brainglobe/
                trapcre_projection.json
                traptta_projection.json
                ratio.json
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Auto-download atlas if not provided locally
    from atlas_manager import get_atlas_volume, get_region_csv
    if atlas_path is None:
        atlas_path = get_atlas_volume(atlas_cache_dir)
    if region_csv_path is None:
        region_csv_path = get_region_csv(atlas_cache_dir)

    atlas = CCFAtlas(atlas_path, region_csv_path)

    animal_dfs = []
    for animal in animals:
        _, animal_df = process_animal(
            animal=animal,
            atlas=atlas,
            injection_regions=injection_regions,
            resize_factor=resize_factor,
            sparse_flag_threshold=sparse_flag_threshold,
            output_dir=output_dir,
        )
        animal_dfs.append(animal_df)

    # Cohort aggregation
    print("\nAggregating cohort...")
    cohort_df = aggregate_cohort(animal_dfs)
    cohort_df.to_csv(output_dir / "cohort_summary.csv", index=False)
    print(f"Cohort summary: {len(cohort_df)} group × region combinations")

    # BrainGlobe exports from first animal as example
    # (in practice average across animals per group)
    bg_dir = output_dir / "brainglobe"
    bg_dir.mkdir(exist_ok=True)
    if animal_dfs:
        avg = pd.concat(animal_dfs).groupby("region_acronym").agg(
            trapcre_fraction=("trapcre_fraction","mean"),
            traptta_fraction=("traptta_fraction","mean"),
            ratio=("trapcre_traptta_ratio","mean"),
        ).reset_index()
        export_brainglobe(avg, bg_dir/"trapcre_projection.json", "trapcre_fraction")
        export_brainglobe(avg, bg_dir/"traptta_projection.json", "traptta_fraction")
        export_brainglobe(avg, bg_dir/"ratio.json", "ratio")

    print(f"\nPipeline complete. Results in: {output_dir}")
    return cohort_df


def _main():
    import argparse

    parser = argparse.ArgumentParser(description="TractQuant pipeline utilities")
    subparsers = parser.add_subparsers(dest="command")

    bg_parser = subparsers.add_parser(
        "diagnose-background",
        help="Compare whole-image and tissue-only modal background",
    )
    bg_parser.add_argument("image", help="Fluorescence image for one section")
    bg_parser.add_argument("visual_align_json", help="Matching VisualAlign JSON")
    bg_parser.add_argument("section_nr", type=int, help="Section number to test")
    bg_parser.add_argument("--atlas-path", default=None, help="Local annotation volume")
    bg_parser.add_argument("--region-csv-path", default=None, help="Local region CSV")
    bg_parser.add_argument(
        "--atlas-cache-dir",
        type=Path,
        default=None,
        help="Atlas cache directory",
    )
    bg_parser.add_argument(
        "--resize-factor",
        type=float,
        default=0.25,
        help="Resize factor for compatibility; default: 0.25",
    )
    bg_parser.add_argument(
        "--no-markers",
        action="store_true",
        help="Disable VisualAlign marker warp for this diagnostic",
    )

    args = parser.parse_args()
    if args.command == "diagnose-background":
        diagnose_background_method(
            image_path=args.image,
            visual_align_json=args.visual_align_json,
            section_nr=args.section_nr,
            atlas_path=args.atlas_path,
            region_csv_path=args.region_csv_path,
            atlas_cache_dir=args.atlas_cache_dir,
            resize_factor=args.resize_factor,
            use_markers=not args.no_markers,
        )
    else:
        parser.print_help()


if __name__ == "__main__":
    _main()
