"""
section_matcher.py — collision-safe matching of VisualAlign JSON sections to
image files, plus a self-contained HTML "contact sheet" for manual validation
*before* running the analysis.

Why this exists
---------------
The fuzzy matcher in pipeline.load_section_images could silently pair an
alignment with the wrong image when (a) the animal ID contained digits that
coincided with a section number (e.g. animal "F12", section 12), or (b) the
TIFFs and the JSON used different numbering conventions, in which case one
image could be assigned to two sections with no error raised.

This module matches strictly by a *canonical section id* (the trailing number
of the filename, after stripping a known animal prefix), reports a confidence
level and a status flag per section, detects duplicate assignments, and lets a
human eyeball the pairing in a browser before committing compute.

Public API
----------
extract_section_id(stem, animal_id=None) -> int | None
match_sections_to_images(sections, ch1_dir, ch2_dir, animal_id=None) -> list[dict]
build_contact_sheet(records, out_path, animal_id="", thumb_px=260) -> Path
"""

from __future__ import annotations

import base64
import io
import re
from pathlib import Path

import numpy as np

IMAGE_EXTS = (".png", ".jpg", ".jpeg", ".tif", ".tiff")


# ---------------------------------------------------------------------------
# Canonical section id
# ---------------------------------------------------------------------------
def extract_section_id(stem: str, animal_id: str | None = None) -> int | None:
    """
    Return the canonical section id for a filename stem.

    Strategy (strict, prefix-aware):
      1. Strip a known animal_id token if supplied (case-insensitive). This is
         what kills the "F12 contains 12" class of bug.
      2. Prefer an explicit section token  's<NN>' / 's_<NN>' / 'sec<NN>'.
      3. Otherwise fall back to the *last* run of digits in the stem, which is
         the section index in essentially every microscopy export convention.
    """
    s = stem.lower()

    if animal_id:
        aid = re.escape(animal_id.strip().lower())
        # remove the animal id wherever it appears, plus an adjacent separator
        s = re.sub(rf"{aid}", " ", s)

    # 2) explicit section marker wins
    m = re.search(r"s(?:ec(?:tion)?)?[_\-\s]*0*(\d+)", s)
    if m:
        return int(m.group(1))

    # 3) trailing number
    m = re.search(r"(\d+)\D*$", s)
    if m:
        return int(m.group(1))

    return None


def _list_images(folder: Path) -> list[Path]:
    folder = Path(folder)
    if not folder.is_dir():
        return []
    out = []
    for f in sorted(folder.iterdir()):
        if f.suffix.lower() in IMAGE_EXTS and not f.name.startswith("._"):
            out.append(f)
    return out


# ---------------------------------------------------------------------------
# Matching
# ---------------------------------------------------------------------------
def _match_one_channel(sec_targets, files, animal_id):
    """
    sec_targets : list of (key, target_id)   key = JSON filename
    files       : list[Path]
    Returns: dict key -> {"path": Path|None, "confidence": str, "candidates": [Path]}
             and a dict file_path -> [keys it was assigned to]  (for dup detection)
    """
    # precompute each file's id
    file_ids = {f: extract_section_id(f.stem, animal_id) for f in files}

    result = {}
    usage = {}  # Path -> [keys]
    for key, target in sec_targets:
        cands = [f for f in files if file_ids[f] == target] if target is not None else []
        if len(cands) == 1:
            chosen, conf = cands[0], "exact"
        elif len(cands) > 1:
            # deterministic but flagged — pick shortest/alpha-first, mark ambiguous
            chosen = sorted(cands, key=lambda p: (len(p.name), p.name))[0]
            conf = "ambiguous"
        else:
            chosen, conf = None, "missing"
        result[key] = {"path": chosen, "confidence": conf, "candidates": cands,
                       "target_id": target}
        if chosen is not None:
            usage.setdefault(chosen, []).append(key)
    return result, usage


def match_sections_to_images(sections: dict, ch1_dir, ch2_dir,
                             animal_id: str | None = None) -> list[dict]:
    """
    sections : dict keyed by JSON filename (output of parse_visual_align),
               each value must contain at least {"nr": int}.
    Returns a list of record dicts, one per JSON section, in JSON order:
        {
          nr, json_file, target_id,
          ch1_path, ch1_id, ch2_path, ch2_id,
          status, notes
        }
    status is one of: ok | ambiguous | missing | duplicate | channel_mismatch
    """
    ch1_files = _list_images(ch1_dir)
    ch2_files = _list_images(ch2_dir)

    # canonical target id for each JSON section: from the JSON *filename*,
    # NOT from nr (nr can be slice order, the comment in pipeline.py says so)
    sec_targets = []
    for json_file, sec in sections.items():
        tid = extract_section_id(Path(json_file).stem, animal_id)
        sec_targets.append((json_file, tid))

    r1, use1 = _match_one_channel(sec_targets, ch1_files, animal_id)
    r2, use2 = _match_one_channel(sec_targets, ch2_files, animal_id)

    # files used by more than one section
    dup1 = {f for f, keys in use1.items() if len(keys) > 1}
    dup2 = {f for f, keys in use2.items() if len(keys) > 1}

    records = []
    for json_file, sec in sections.items():
        a, b = r1[json_file], r2[json_file]
        notes = []
        status = "ok"

        if a["confidence"] == "missing" or b["confidence"] == "missing":
            status = "missing"
            if a["confidence"] == "missing":
                notes.append("no Ch1 image with this section id")
            if b["confidence"] == "missing":
                notes.append("no Ch2 image with this section id")
        elif a["path"] in dup1 or b["path"] in dup2:
            status = "duplicate"
            notes.append("this image is also assigned to another section")
        elif a["confidence"] == "ambiguous" or b["confidence"] == "ambiguous":
            status = "ambiguous"
            notes.append("multiple images share this section id")

        # channel consistency: both channels should resolve to the same id
        if a["path"] is not None and b["path"] is not None:
            if a["target_id"] != b["target_id"]:
                pass  # same target by construction
            if extract_section_id(a["path"].stem, animal_id) != \
               extract_section_id(b["path"].stem, animal_id):
                status = "channel_mismatch"
                notes.append("Ch1 and Ch2 resolved to different section ids")

        records.append({
            "nr": sec.get("nr"),
            "json_file": json_file,
            "target_id": a["target_id"],
            "ch1_path": a["path"],
            "ch1_id": extract_section_id(a["path"].stem, animal_id) if a["path"] else None,
            "ch2_path": b["path"],
            "ch2_id": extract_section_id(b["path"].stem, animal_id) if b["path"] else None,
            "status": status,
            "notes": "; ".join(notes),
        })
    return records


def summarize(records: list[dict]) -> dict:
    counts = {"ok": 0, "ambiguous": 0, "missing": 0, "duplicate": 0,
              "channel_mismatch": 0, "low_coverage": 0}
    for r in records:
        counts[r["status"]] = counts.get(r["status"], 0) + 1
    counts["total"] = len(records)
    counts["problems"] = counts["total"] - counts["ok"]
    return counts


# ---------------------------------------------------------------------------
# Analysis-ready loader (mirrors pipeline.load_section_images' reader exactly)
# ---------------------------------------------------------------------------
def read_section_array(path) -> np.ndarray:
    """Read one image as a 2-D float32 array, identical to the pipeline reader."""
    path = Path(path)
    if path.suffix.lower() in (".tif", ".tiff"):
        import tifffile
        img = tifffile.imread(str(path)).astype(np.float32)
    else:
        from PIL import Image
        img = np.asarray(Image.open(path).convert("F"), dtype=np.float32)
    if img.ndim == 3:  # multichannel tif -> first plane
        img = img[0] if img.shape[0] < img.shape[-1] else img[:, :, 0]
    return img


def load_validated_pair(record: dict) -> tuple[np.ndarray, np.ndarray]:
    """Load the Ch1/Ch2 arrays for a matched record. Raises if not loadable."""
    if record.get("ch1_path") is None or record.get("ch2_path") is None:
        raise FileNotFoundError(f"section nr={record.get('nr')} has no matched image "
                                f"(status={record.get('status')})")
    return read_section_array(record["ch1_path"]), read_section_array(record["ch2_path"])


# ---------------------------------------------------------------------------
# Thumbnails
# ---------------------------------------------------------------------------
def _read_image_any(path: Path) -> np.ndarray | None:
    try:
        if path.suffix.lower() in (".tif", ".tiff"):
            import tifffile
            arr = tifffile.imread(str(path))
        else:
            from PIL import Image
            arr = np.asarray(Image.open(path))
        arr = np.asarray(arr)
        if arr.ndim == 3:
            # collapse to single plane: take max-projection of first 3 bands
            if arr.shape[-1] in (3, 4):
                arr = arr[..., :3].max(axis=-1)
            else:
                arr = arr[0]
        return arr.astype(np.float32)
    except Exception:
        return None


def _thumb_data_uri_from_rgb(path: Path, thumb_px: int) -> str | None:
    """Thumbnail an already-rendered RGB image (e.g. an overlay PNG)."""
    try:
        from PIL import Image
        img = Image.open(path).convert("RGB")
        img.thumbnail((thumb_px, thumb_px))
        buf = io.BytesIO()
        img.save(buf, format="PNG")
        return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()
    except Exception:
        return None


def build_overlay_records(records, sections, atlas, output_dir, animal_id,
                          pipeline_module, resize_factor=0.25,
                          thumb_px=420, progress=None):
    """
    For each matched section, render the atlas-on-image overlay using the
    pipeline's own save_alignment_overlay (single source of truth with the
    analysis), and attach a thumbnail + QC numbers to the record.

    `pipeline_module` is your imported `pipeline` (passed in to keep this
    module free of a hard dependency on it).
    Returns a new list of enriched records ready for build_contact_sheet(...).
    """
    pl = pipeline_module
    enriched = []
    n = len(records)
    for i, r in enumerate(records):
        info = dict(r)
        if r["ch1_path"] and r["ch2_path"] and r["status"] != "missing":
            try:
                ch1, ch2 = load_validated_pair(r)
                sec = sections[r["json_file"]]
                qc = pl.save_alignment_overlay(
                    sec, ch1, ch2, atlas, output_dir, animal_id,
                    r["nr"], r["json_file"], resize_factor,
                )
                info["overlay_uri"] = _thumb_data_uri_from_rgb(
                    Path(qc["side_by_side"]), thumb_px)
                info["atlas_percent"] = (qc.get("ch1_atlas_percent"),
                                         qc.get("ch2_atlas_percent"))
                info["reference_iou"] = qc.get("ch1_reference_iou")
                info["use_markers"] = qc.get("use_markers")
                info["marker_count"] = qc.get("marker_count")
                # soft alignment-quality flag: atlas barely lands on tissue
                ap = max(v for v in info["atlas_percent"] if v is not None) \
                    if any(v is not None for v in info["atlas_percent"]) else 0
                if r["status"] == "ok" and ap < 1.0:
                    info["status"] = "low_coverage"
                    info["notes"] = (info["notes"] + "; " if info["notes"] else "") + \
                        f"atlas covers only {ap:.1f}% of the image — check anchoring"
            except Exception as e:
                info["overlay_error"] = str(e)
        enriched.append(info)
        if progress:
            progress(i + 1, n)
    return enriched


def _thumb_data_uri(path: Path, thumb_px: int) -> str | None:
    arr = _read_image_any(path)
    if arr is None:
        return None
    from PIL import Image
    # percentile stretch so faint fluorescence is visible
    lo, hi = np.percentile(arr, [1, 99.5])
    if hi <= lo:
        hi = lo + 1
    arr = np.clip((arr - lo) / (hi - lo), 0, 1)
    img = Image.fromarray((arr * 255).astype(np.uint8))
    img.thumbnail((thumb_px, thumb_px))
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


# ---------------------------------------------------------------------------
# Contact sheet (self-contained HTML)
# ---------------------------------------------------------------------------
_STATUS_COLOR = {
    "ok": "#1D9E75",
    "ambiguous": "#EF9F27",
    "missing": "#E24B4A",
    "duplicate": "#E24B4A",
    "channel_mismatch": "#E24B4A",
    "low_coverage": "#EF9F27",
}
_STATUS_LABEL = {
    "ok": "matched",
    "ambiguous": "ambiguous",
    "missing": "missing",
    "duplicate": "duplicate image",
    "channel_mismatch": "channel mismatch",
    "low_coverage": "low coverage",
}


def build_contact_sheet(records: list[dict], out_path, animal_id: str = "",
                        thumb_px: int = 260) -> Path:
    """Write a self-contained HTML contact sheet for manual validation."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    s = summarize(records)

    cards = []
    for r in records:
        sc = _STATUS_COLOR.get(r["status"], "#888")
        sl = _STATUS_LABEL.get(r["status"], r["status"])

        def cell(path, cid, label):
            if path is None:
                inner = (f'<div class="missing">no image<br>'
                         f'<span>(section id {r["target_id"]})</span></div>')
            else:
                uri = _thumb_data_uri(Path(path), thumb_px)
                if uri:
                    inner = f'<img src="{uri}" alt="{label}">'
                else:
                    inner = f'<div class="missing">unreadable<br><span>{Path(path).name}</span></div>'
                idtag = (f'<span class="idtag" style="color:{sc}">id {cid}</span>'
                         if cid != r["target_id"] else f'<span class="idtag">id {cid}</span>')
                inner += (f'<div class="fname">{Path(path).name} {idtag}</div>')
            return f'<div class="cell"><div class="celllabel">{label}</div>{inner}</div>'

        # overlay mode: a single atlas-on-image picture + QC line
        if r.get("overlay_uri"):
            ap = r.get("atlas_percent") or (None, None)
            iou = r.get("reference_iou")
            qc_bits = []
            if ap[0] is not None:
                qc_bits.append(f'atlas mask {ap[0]:.1f}% / {ap[1]:.1f}%')
            if iou is not None:
                qc_bits.append(f'ref IoU {iou:.2f}')
            qc_bits.append('markers' if r.get("use_markers") else 'anchoring only')
            body = (f'<div class="overlaywrap"><img class="overlay" src="{r["overlay_uri"]}" '
                    f'alt="atlas overlay"></div>'
                    f'<div class="qc">{" · ".join(qc_bits)}</div>'
                    f'<div class="fnames">Ch1 {Path(r["ch1_path"]).name} &nbsp;·&nbsp; '
                    f'Ch2 {Path(r["ch2_path"]).name}</div>')
        elif r.get("overlay_error"):
            body = f'<div class="cells"><div class="missing" style="grid-column:1/3">overlay failed<br><span>{r["overlay_error"]}</span></div></div>'
        else:
            body = f'<div class="cells">{cell(r["ch1_path"], r["ch1_id"], "Ch1")}{cell(r["ch2_path"], r["ch2_id"], "Ch2")}</div>'

        note = f'<div class="note">{r["notes"]}</div>' if r["notes"] else ""
        cards.append(f"""
        <div class="card" data-status="{r['status']}">
          <div class="cardhead">
            <div class="secid">§ {r['target_id'] if r['target_id'] is not None else '?'}</div>
            <div class="meta">
              <div class="jname">{r['json_file']}</div>
              <div class="sub">JSON nr {r['nr']}</div>
            </div>
            <div class="badge" style="background:{sc}">{sl}</div>
          </div>
          {body}
          {note}
        </div>""")

    chips = "".join(
        f'<span class="chip" style="--c:{_STATUS_COLOR.get(k,"#888")}">'
        f'{s.get(k,0)} {_STATUS_LABEL.get(k,k)}</span>'
        for k in ["ok", "ambiguous", "low_coverage", "missing", "duplicate", "channel_mismatch"]
        if s.get(k, 0)
    )
    headline = ("All sections matched cleanly — safe to run."
                if s["problems"] == 0 else
                f"{s['problems']} of {s['total']} sections need a look before running.")
    head_color = "#1D9E75" if s["problems"] == 0 else "#EF9F27"

    html = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>TractQuant — section validation{(' · ' + animal_id) if animal_id else ''}</title>
<style>
  :root {{
    --bg:#0f1117; --bg2:#1a1d27; --bg3:#252836; --fg:#e8e8e8; --fg2:#8a8fa3;
    --blue:#378ADD; --border:#2e3147;
  }}
  * {{ box-sizing:border-box; }}
  body {{ margin:0; background:var(--bg); color:var(--fg);
    font-family:"Segoe UI",-apple-system,BlinkMacSystemFont,sans-serif; }}
  header {{ position:sticky; top:0; z-index:5; background:rgba(15,17,23,.92);
    backdrop-filter:blur(8px); border-bottom:1px solid var(--border);
    padding:18px 24px; }}
  h1 {{ margin:0; font-size:15px; letter-spacing:.14em; text-transform:uppercase;
    color:var(--blue); font-weight:700; }}
  .headline {{ margin-top:6px; font-size:20px; font-weight:700; color:{head_color}; }}
  .chips {{ margin-top:10px; display:flex; gap:8px; flex-wrap:wrap; }}
  .chip {{ font-size:12px; font-weight:700; padding:4px 10px; border-radius:999px;
    color:var(--c); border:1px solid var(--c); background:color-mix(in srgb,var(--c) 12%,transparent); }}
  .filters {{ margin-top:12px; display:flex; gap:6px; flex-wrap:wrap; }}
  .filters button {{ font:inherit; font-size:12px; cursor:pointer; color:var(--fg2);
    background:var(--bg3); border:1px solid var(--border); border-radius:6px;
    padding:5px 12px; }}
  .filters button.active {{ color:var(--blue); border-color:var(--blue); }}
  main {{ padding:20px 24px 60px; display:grid;
    grid-template-columns:repeat(auto-fill,minmax(340px,1fr)); gap:16px; }}
  .card {{ background:var(--bg2); border:1px solid var(--border); border-radius:12px;
    overflow:hidden; }}
  .card[data-status="missing"],.card[data-status="duplicate"],
  .card[data-status="channel_mismatch"] {{ border-color:#E24B4A; }}
  .card[data-status="ambiguous"] {{ border-color:#EF9F27; }}
  .cardhead {{ display:flex; align-items:center; gap:12px; padding:12px 14px;
    border-bottom:1px solid var(--border); }}
  .secid {{ font-size:22px; font-weight:800; color:var(--blue);
    font-variant-numeric:tabular-nums; min-width:42px; }}
  .meta {{ flex:1; min-width:0; }}
  .jname {{ font-size:12.5px; font-weight:600; white-space:nowrap; overflow:hidden;
    text-overflow:ellipsis; }}
  .sub {{ font-size:11px; color:var(--fg2); margin-top:1px; }}
  .badge {{ font-size:11px; font-weight:700; color:#0f1117; padding:3px 9px;
    border-radius:999px; white-space:nowrap; }}
  .cells {{ display:grid; grid-template-columns:1fr 1fr; gap:1px; background:var(--border); }}
  .overlaywrap {{ background:#000; display:flex; justify-content:center; }}
  .overlay {{ max-width:100%; display:block; }}
  .qc {{ font-family:ui-monospace,Consolas,monospace; font-size:11.5px; color:var(--blue);
    padding:8px 14px 2px; }}
  .fnames {{ font-size:10.5px; color:var(--fg2); padding:0 14px 10px; word-break:break-all; }}
  .cell {{ background:var(--bg2); padding:10px; min-height:120px;
    display:flex; flex-direction:column; align-items:center; }}
  .celllabel {{ font-size:10px; letter-spacing:.12em; text-transform:uppercase;
    color:var(--fg2); align-self:flex-start; margin-bottom:6px; }}
  .cell img {{ max-width:100%; border-radius:6px; background:#000; }}
  .fname {{ font-size:10.5px; color:var(--fg2); margin-top:6px; text-align:center;
    word-break:break-all; }}
  .idtag {{ font-family:ui-monospace,Consolas,monospace; color:var(--fg2); }}
  .missing {{ display:flex; flex-direction:column; justify-content:center;
    flex:1; color:#E24B4A; font-size:13px; font-weight:700; text-align:center; }}
  .missing span {{ color:var(--fg2); font-weight:400; font-size:11px; }}
  .note {{ padding:9px 14px; font-size:12px; color:#EF9F27;
    border-top:1px solid var(--border); }}
</style></head><body>
<header>
  <h1>TractQuant · section ↔ image validation{(' · ' + animal_id) if animal_id else ''}</h1>
  <div class="headline">{headline}</div>
  <div class="chips">{chips}</div>
  <div class="filters">
    <button data-f="all" class="active">All ({s['total']})</button>
    <button data-f="problem">Problems only ({s['problems']})</button>
  </div>
</header>
<main id="grid">{''.join(cards)}</main>
<script>
  const btns=document.querySelectorAll('.filters button');
  const cards=[...document.querySelectorAll('.card')];
  btns.forEach(b=>b.onclick=()=>{{
    btns.forEach(x=>x.classList.remove('active')); b.classList.add('active');
    const f=b.dataset.f;
    cards.forEach(c=>{{
      const ok=c.dataset.status==='ok';
      c.style.display=(f==='all'||!ok)?'':'none';
    }});
  }});
</script>
</body></html>"""
    out_path.write_text(html, encoding="utf-8")
    return out_path
