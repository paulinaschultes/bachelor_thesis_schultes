"""
TractQuant — atlas manager
==========================
Finds the Allen CCF atlas in this order:

  1. Explicit VisualAlign/Cutlas atlas used by the alignment JSON
  2. BrainGlobe cache (~/.brainglobe/allen_mouse_25um_v1.2/)  ← fallback
  3. Local folder next to this script (annotation_25.nrrd)
  4. atlas/ subfolder next to this script
  5. User cache (~/.tractquant/atlas/)
  6. Download via BrainGlobe API
  7. Direct download from Allen Institute (fallback)

RECOMMENDED SETUP:
    conda activate tractquant
    pip install brainglobe-atlasapi
    python -c "from brainglobe_atlasapi import BrainGlobeAtlas; BrainGlobeAtlas('allen_mouse_25um')"

Then TractQuant finds it automatically — nothing else needed.
"""

import csv
import json
import shutil
import time
import urllib.error
import urllib.request
from pathlib import Path

SCRIPT_DIR        = Path(__file__).parent.resolve()
DEFAULT_CACHE     = Path.home() / ".tractquant" / "atlas"
ATLAS_FILENAME    = "annotation_25.nrrd"
REGIONS_FILENAME  = "ccf_regions.csv"
VISUALIGN_ATLAS_NAME = "ABA_Mouse_CCFv3_2017_25um.cutlas"
VISUALIGN_ATLAS_PATH = (
    Path.home() / "Desktop" / "VisuAlign.app" /
    VISUALIGN_ATLAS_NAME / "labels.nii.gz"
)
ATLAS_SIZE_MB_MIN_NRRD = 400   # uncompressed nrrd
ATLAS_SIZE_MB_MIN_TIFF = 100   # compressed tiff (BrainGlobe uses ~294 MB)
ATLAS_SIZE_MB_MIN = 100        # general minimum

ATLAS_URL = (
    "http://download.alleninstitute.org/informatics-archive/"
    "current-release/mouse_ccf/annotation/ccf_2017/annotation_25.nrrd"
)

ALLEN_API_URL = (
    "https://api.brain-map.org/api/v2/data/query.json"
    "?criteria=model::Structure"
    ",rma::criteria,[graph_id$eq1]"
    ",rma::options[num_rows$eq'all'][order$eq'structures.id']"
    "[include$eq'structure_id_path']"
)


# ── BrainGlobe search ─────────────────────────────────────────────────

def find_visualign_atlas():
    """
    Return the exact atlas shipped with VisualAlign/Cutlas.

    VisualAlign stores the anchoring vectors in the coordinate system of this
    atlas. Using BrainGlobe's TIFF fallback can be geometrically correct in
    content but wrong in axis order for these JSON files.
    """
    if VISUALIGN_ATLAS_PATH.exists() and VISUALIGN_ATLAS_PATH.stat().st_size > 1000:
        print(f"VisualAlign atlas selected: {VISUALIGN_ATLAS_PATH}")
        return str(VISUALIGN_ATLAS_PATH)
    return None

def find_brainglobe_atlas():
    """
    Look for the atlas in BrainGlobe's cache.
    BrainGlobe stores atlases in ~/.brainglobe/{atlas_name}_v{version}/
    The annotation file is called 'annotation.tiff' or 'annotation.nrrd'
    depending on version.
    """
    bg_dir = Path.home() / ".brainglobe"
    if not bg_dir.exists():
        return None, None

    # Search all allen_mouse_25um versions
    for atlas_dir in sorted(bg_dir.glob("allen_mouse_25um*"), reverse=True):
        if not atlas_dir.is_dir():
            continue

        # Try annotation file names BrainGlobe uses
        for fname in ["annotation.tiff", "annotation.tif",
                       "annotation.nrrd", "annotation_25.nrrd"]:
            ann = atlas_dir / fname
            if ann.exists() and ann.stat().st_size > 1000:
                # Also look for structures CSV
                struct_file = None
                for sname in ["structures.csv","structures.json",
                               "structures.txt"]:
                    sf = atlas_dir / sname
                    if sf.exists():
                        struct_file = sf
                        break
                print(f"BrainGlobe atlas found: {ann}")
                return str(ann), struct_file

    return None, None


def load_brainglobe_via_api():
    """
    Use BrainGlobe Python API to get atlas paths.
    Returns (annotation_path, structures_df) or (None, None).
    """
    try:
        from brainglobe_atlasapi import BrainGlobeAtlas
        atlas = BrainGlobeAtlas("allen_mouse_25um", check_latest=False)
        ann_path = atlas.root_dir / "annotation.tiff"
        if not ann_path.exists():
            # older versions use different name
            for fname in ["annotation.nrrd","annotation.tif"]:
                p = atlas.root_dir / fname
                if p.exists():
                    ann_path = p
                    break
        struct_path = atlas.root_dir / "structures.csv"
        print(f"BrainGlobe API: {ann_path}")
        return str(ann_path), str(struct_path) if struct_path.exists() else None
    except Exception as e:
        print(f"  BrainGlobe API: {e}")
        return None, None


# ── Local file search ─────────────────────────────────────────────────

def _find_local(filename):
    candidates = [
        SCRIPT_DIR / filename,
        SCRIPT_DIR / "atlas" / filename,
        DEFAULT_CACHE / filename,
    ]
    for p in candidates:
        if p.exists() and p.stat().st_size > 1000:
            return p
    return None


def _validate_nrrd(path):
    if not path or not Path(path).exists():
        return False
    size_mb = Path(path).stat().st_size / 1024 / 1024
    if size_mb < ATLAS_SIZE_MB_MIN:
        print(f"  File too small ({size_mb:.0f} MB, expected >{ATLAS_SIZE_MB_MIN} MB)")
        return False
    try:
        with open(path, "rb") as f:
            magic = f.read(4)
        # Accept NRRD or TIFF (BrainGlobe uses tiff)
        if not (magic.startswith(b"NRRD") or
                magic[:2] in (b"II", b"MM") or  # TIFF little/big endian
                magic == b"\x49\x49\x2a\x00"):
            pass  # don't fail — some valid files have other headers
    except Exception:
        pass
    return True


# ── Region CSV ────────────────────────────────────────────────────────

def _build_regions_from_brainglobe(struct_file):
    """Convert BrainGlobe structures file to TractQuant CSV format."""
    dest = DEFAULT_CACHE / REGIONS_FILENAME
    DEFAULT_CACHE.mkdir(parents=True, exist_ok=True)

    try:
        import pandas as pd
        df = pd.read_csv(struct_file) if str(struct_file).endswith(".csv") else None
        if df is None:
            return None

        # BrainGlobe CSV has columns: id, name, acronym, structure_id_path, rgb_triplet
        rows = []
        id_to_name = dict(zip(df["id"], df["name"])) if "id" in df.columns else {}

        def get_division(path_str):
            if not isinstance(path_str, str):
                return ""
            parts = [p for p in path_str.strip("/").split("/") if p]
            if len(parts) >= 2:
                return id_to_name.get(int(parts[1]), "")
            return ""

        for _, row in df.iterrows():
            rows.append({
                "id":                row.get("id",""),
                "acronym":           row.get("acronym", str(row.get("id",""))),
                "name":              row.get("name",""),
                "parent_id":         row.get("parent_structure_id",""),
                "depth":             row.get("depth",""),
                "structure_id_path": row.get("structure_id_path",""),
                "division":          get_division(row.get("structure_id_path","")),
            })

        with open(dest, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
        print(f"  Region CSV built from BrainGlobe: {dest}")
        return dest

    except Exception as e:
        print(f"  Could not build regions from BrainGlobe: {e}")
        return None


# ── Public API ────────────────────────────────────────────────────────

def get_atlas_volume(cache_dir=None):
    """
    Return path to atlas annotation volume.

    The VisualAlign/Cutlas atlas is selected first because the alignment JSON
    was generated against that volume and axis order.
    """
    # 1. Explicit VisualAlign atlas matching the JSON alignment.
    ann = find_visualign_atlas()
    if ann:
        return ann

    # 2. BrainGlobe API fallback
    ann, _ = load_brainglobe_via_api()
    if ann and Path(ann).exists():
        size_mb = Path(ann).stat().st_size / 1024 / 1024
        if size_mb > ATLAS_SIZE_MB_MIN:
            return ann

    # 3. BrainGlobe cache folder (manual search)
    ann, _ = find_brainglobe_atlas()
    if ann and Path(ann).exists():
        size_mb = Path(ann).stat().st_size / 1024 / 1024
        if size_mb > ATLAS_SIZE_MB_MIN:
            return ann

    # 4. Local nrrd file next to script
    found = _find_local(ATLAS_FILENAME)
    if found and _validate_nrrd(found):
        return str(found)

    # 5. Download via BrainGlobe API
    print("\nAtlas not found — attempting download via BrainGlobe...")
    try:
        from brainglobe_atlasapi import BrainGlobeAtlas
        print("Downloading allen_mouse_25um via BrainGlobe (this may take a few minutes)...")
        atlas = BrainGlobeAtlas("allen_mouse_25um")
        ann_path = atlas.root_dir / "annotation.tiff"
        if not ann_path.exists():
            for fname in ["annotation.nrrd","annotation.tif"]:
                p = atlas.root_dir / fname
                if p.exists():
                    ann_path = p
                    break
        if ann_path.exists():
            print(f"✓ Downloaded: {ann_path}")
            return str(ann_path)
    except ImportError:
        print("  brainglobe-atlasapi not installed.")
        print("  Run: pip install brainglobe-atlasapi")
    except Exception as e:
        print(f"  BrainGlobe download failed: {e}")

    # 5. Direct Allen download
    dest_dir = Path(cache_dir) if cache_dir else DEFAULT_CACHE
    dest = dest_dir / ATLAS_FILENAME
    dest_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*60}")
    print("Attempting direct download from Allen Institute...")
    print(f"Saving to: {dest}")
    print(f"\nTIP: For reliable download run:")
    print(f"  pip install brainglobe-atlasapi")
    print(f"  python -c \"from brainglobe_atlasapi import BrainGlobeAtlas; BrainGlobeAtlas('allen_mouse_25um')\"")
    print(f"{'='*60}\n")

    try:
        class Progress:
            def __init__(self): self.last=0; self.start=time.time()
            def __call__(self,b,bs,ts):
                done=b*bs; pct=min(int(done/ts*100),100) if ts>0 else 0
                if pct!=self.last and pct%5==0:
                    bar="█"*(pct//5)+"░"*(20-pct//5)
                    mb=done/1024/1024; spd=mb/(time.time()-self.start+0.001)
                    print(f"\r  [{bar}] {pct}%  {mb:.0f}MB  {spd:.1f}MB/s",end="",flush=True)
                    self.last=pct
                if pct>=100: print()

        tmp = dest.with_suffix(".tmp")
        urllib.request.urlretrieve(ATLAS_URL, tmp, reporthook=Progress())
        if tmp.stat().st_size / 1024/1024 < ATLAS_SIZE_MB_MIN:
            tmp.unlink(missing_ok=True)
            raise RuntimeError("Download incomplete")
        shutil.move(str(tmp), dest)
        print(f"✓ Downloaded: {dest}")
        return str(dest)
    except Exception as e:
        raise RuntimeError(
            f"Could not obtain atlas: {e}\n\n"
            f"Please run:\n"
            f"  pip install brainglobe-atlasapi\n"
            f"  python -c \"from brainglobe_atlasapi import BrainGlobeAtlas; "
            f"BrainGlobeAtlas('allen_mouse_25um')\"\n"
            f"Then restart TractQuant."
        )


def get_region_csv(cache_dir=None):
    """
    Return path to region CSV.
    Tries BrainGlobe structures first, then Allen API.
    """
    # Check local first
    found = _find_local(REGIONS_FILENAME)
    if found:
        print(f"Region CSV found: {found}")
        return str(found)

    # Try BrainGlobe structures
    _, struct_file = load_brainglobe_via_api()
    if struct_file and Path(struct_file).exists():
        dest = _build_regions_from_brainglobe(struct_file)
        if dest:
            return str(dest)

    _, struct_file = find_brainglobe_atlas()
    if struct_file and Path(struct_file).exists():
        dest = _build_regions_from_brainglobe(struct_file)
        if dest:
            return str(dest)

    # Fall back to Allen API
    dest_dir = Path(cache_dir) if cache_dir else DEFAULT_CACHE
    dest = dest_dir / REGIONS_FILENAME
    dest_dir.mkdir(parents=True, exist_ok=True)

    print("Fetching region list from Allen Brain Atlas API...")
    try:
        with urllib.request.urlopen(ALLEN_API_URL, timeout=30) as resp:
            data = json.loads(resp.read())
    except Exception as e:
        raise RuntimeError(
            f"Failed to fetch region list: {e}\n"
            f"Check your internet connection."
        ) from e

    structures = data.get("msg", [])
    print(f"  Retrieved {len(structures)} regions")
    id_to_name = {s["id"]: s["name"] for s in structures}

    def get_division(path):
        parts = [p for p in path.strip("/").split("/") if p]
        return id_to_name.get(int(parts[1]), "") if len(parts) >= 2 else ""

    rows = [{
        "id":                s["id"],
        "acronym":           s["acronym"],
        "name":              s["name"],
        "parent_id":         s.get("parent_structure_id", ""),
        "depth":             s.get("depth", ""),
        "structure_id_path": s.get("structure_id_path", ""),
        "division":          get_division(s.get("structure_id_path", "")),
    } for s in structures]

    with open(dest, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    print(f"  Region CSV saved: {dest}")
    return str(dest)


def get_atlas(cache_dir=None):
    return get_atlas_volume(cache_dir), get_region_csv(cache_dir)


def print_status():
    print(f"\nTractQuant Atlas Status")
    print(f"{'='*40}")

    visualign = find_visualign_atlas()
    if visualign:
        mb = Path(visualign).stat().st_size/1024/1024
        print(f"✓ VisualAlign atlas: {visualign} ({mb:.1f} MB)")
    else:
        print(f"✗ VisualAlign atlas: not found")

    # BrainGlobe
    ann, struct = load_brainglobe_via_api()
    if ann and Path(ann).exists():
        mb = Path(ann).stat().st_size/1024/1024
        print(f"✓ BrainGlobe atlas:  {ann} ({mb:.0f} MB)")
    else:
        ann, struct = find_brainglobe_atlas()
        if ann:
            mb = Path(ann).stat().st_size/1024/1024
            print(f"✓ BrainGlobe cache:  {ann} ({mb:.0f} MB)")
        else:
            print(f"✗ BrainGlobe atlas:  not found")
            print(f"  → Run: pip install brainglobe-atlasapi")
            print(f"  → Run: python -c \"from brainglobe_atlasapi import BrainGlobeAtlas; BrainGlobeAtlas('allen_mouse_25um')\"")

    # Local
    local = _find_local(ATLAS_FILENAME)
    if local:
        mb = local.stat().st_size/1024/1024
        print(f"✓ Local nrrd:        {local} ({mb:.0f} MB)")
    else:
        print(f"  Local nrrd:        not found")

    # Region CSV
    csv_f = _find_local(REGIONS_FILENAME)
    if csv_f:
        print(f"✓ Region CSV:        {csv_f}")
    else:
        print(f"  Region CSV:        not found (will download from Allen API)")
    print()


if __name__ == "__main__":
    import sys
    cmd = sys.argv[1] if len(sys.argv) > 1 else "status"
    if cmd == "status":
        print_status()
    elif cmd == "download":
        get_atlas()
    else:
        print("Usage: python atlas_manager.py [status|download]")
