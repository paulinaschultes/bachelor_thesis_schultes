"""
TractQuant — Streamlit App
===========================
Student-friendly GUI for the anterograde tracing quantification pipeline.
Run with:  streamlit run app.py
"""

import streamlit as st
import pandas as pd
import numpy as np
import json
import shutil
import sys
import subprocess
from pathlib import Path
from io import StringIO, BytesIO
import time

# ---------------------------------------------------------------------------
# Page config — must be first Streamlit call
# ---------------------------------------------------------------------------
st.set_page_config(
    page_title="TractQuant",
    page_icon="🧠",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ---------------------------------------------------------------------------
# Custom CSS
# ---------------------------------------------------------------------------
st.markdown("""
<style>
/* Import a clean scientific font */
@import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@300;400;500;600&family=IBM+Plex+Mono:wght@400&display=swap');

html, body, [class*="css"] {
    font-family: 'IBM Plex Sans', sans-serif;
}

/* Sidebar */
[data-testid="stSidebar"] {
    background: #0f1117;
}
[data-testid="stSidebar"] * {
    color: #e0e0e0 !important;
}
[data-testid="stSidebar"] .sidebar-title {
    font-size: 20px;
    font-weight: 600;
    color: #ffffff !important;
    letter-spacing: -0.3px;
    margin-bottom: 2px;
}
[data-testid="stSidebar"] .sidebar-sub {
    font-size: 11px;
    color: #666 !important;
    letter-spacing: 0.05em;
    text-transform: uppercase;
    margin-bottom: 20px;
}

/* Hide default streamlit elements */
#MainMenu { visibility: hidden; }
footer { visibility: hidden; }
header { visibility: hidden; }

/* Step indicator in sidebar */
.step-item {
    display: flex;
    align-items: center;
    gap: 10px;
    padding: 7px 0;
    cursor: pointer;
}
.step-dot {
    width: 8px; height: 8px;
    border-radius: 50%;
    background: #444;
    flex-shrink: 0;
}
.step-dot.done { background: #1D9E75; }
.step-dot.active { background: #7C6FE0; }
.step-label { font-size: 13px; }

/* Cards */
.tq-card {
    background: #ffffff;
    border: 1px solid #f0f0f0;
    border-radius: 12px;
    padding: 20px 24px;
    margin-bottom: 16px;
    box-shadow: 0 1px 3px rgba(0,0,0,0.04);
}
.tq-card h4 {
    font-size: 13px;
    font-weight: 600;
    color: #111;
    margin: 0 0 14px 0;
    letter-spacing: -0.1px;
}

/* Channel blocks */
.ch-block-cre {
    border-left: 3px solid #534AB7;
    background: #fafafe;
    border-radius: 0 8px 8px 0;
    padding: 14px 16px;
    margin-bottom: 10px;
}
.ch-block-tta {
    border-left: 3px solid #1D9E75;
    background: #fafffe;
    border-radius: 0 8px 8px 0;
    padding: 14px 16px;
    margin-bottom: 10px;
}
.ch-label-cre {
    display: inline-block;
    background: #EEEDFE;
    color: #3C3489;
    font-size: 11px;
    font-weight: 600;
    padding: 2px 10px;
    border-radius: 99px;
    margin-bottom: 10px;
}
.ch-label-tta {
    display: inline-block;
    background: #E1F5EE;
    color: #0F6E56;
    font-size: 11px;
    font-weight: 600;
    padding: 2px 10px;
    border-radius: 99px;
    margin-bottom: 10px;
}

/* Animal row in list */
.animal-row {
    background: #fafafa;
    border: 1px solid #efefef;
    border-radius: 8px;
    padding: 12px 16px;
    margin-bottom: 8px;
}
.animal-id {
    font-weight: 600;
    font-size: 14px;
    color: #111;
}
.animal-meta {
    font-size: 12px;
    color: #777;
    margin-top: 2px;
}
.compound-key {
    font-family: 'IBM Plex Mono', monospace;
    font-size: 11px;
    color: #555;
    background: #f5f5f5;
    padding: 2px 8px;
    border-radius: 4px;
    display: inline-block;
    margin-top: 3px;
}

/* Metric cards */
.metric-row {
    display: flex;
    gap: 12px;
    margin-bottom: 16px;
}
.metric-card {
    flex: 1;
    background: #f8f8f8;
    border-radius: 10px;
    padding: 14px 16px;
    text-align: center;
}
.metric-val {
    font-size: 28px;
    font-weight: 600;
    color: #111;
    line-height: 1;
}
.metric-lbl {
    font-size: 11px;
    color: #888;
    margin-top: 4px;
    text-transform: uppercase;
    letter-spacing: 0.04em;
}

/* Tags */
.tag {
    display: inline-block;
    padding: 2px 8px;
    border-radius: 99px;
    font-size: 11px;
    font-weight: 500;
    margin: 1px;
}
.tag-purple { background: #EEEDFE; color: #3C3489; }
.tag-green  { background: #E1F5EE; color: #0F6E56; }
.tag-amber  { background: #FAEEDA; color: #633806; }
.tag-gray   { background: #f0f0f0; color: #555; }
.tag-pink   { background: #FBEAF0; color: #72243E; }
.tag-blue   { background: #E6F1FB; color: #0C447C; }

/* Banner */
.info-banner {
    background: #E6F1FB;
    border-left: 3px solid #378ADD;
    border-radius: 0 6px 6px 0;
    padding: 10px 14px;
    font-size: 13px;
    color: #0C447C;
    margin-bottom: 14px;
}
.warn-banner {
    background: #FAEEDA;
    border-left: 3px solid #EF9F27;
    border-radius: 0 6px 6px 0;
    padding: 10px 14px;
    font-size: 13px;
    color: #633806;
    margin-bottom: 14px;
}
.ok-banner {
    background: #E1F5EE;
    border-left: 3px solid #1D9E75;
    border-radius: 0 6px 6px 0;
    padding: 10px 14px;
    font-size: 13px;
    color: #085041;
    margin-bottom: 14px;
}

/* Log box */
.log-box {
    background: #0f1117;
    border-radius: 8px;
    padding: 14px 16px;
    font-family: 'IBM Plex Mono', monospace;
    font-size: 12px;
    color: #aaa;
    max-height: 280px;
    overflow-y: auto;
    line-height: 1.7;
}
.log-ok   { color: #1D9E75; }
.log-warn { color: #EF9F27; }
.log-err  { color: #E24B4A; }

/* Page header */
.page-header {
    margin-bottom: 24px;
    padding-bottom: 16px;
    border-bottom: 1px solid #f0f0f0;
}
.page-header h2 {
    font-size: 22px;
    font-weight: 600;
    color: #111;
    margin: 0 0 4px 0;
    letter-spacing: -0.3px;
}
.page-header p {
    font-size: 13px;
    color: #888;
    margin: 0;
}

/* Streamlit overrides */
.stButton button {
    border-radius: 8px;
    font-family: 'IBM Plex Sans', sans-serif;
    font-weight: 500;
}
.stTextInput input, .stSelectbox select {
    border-radius: 8px;
    font-family: 'IBM Plex Sans', sans-serif;
}
div[data-testid="stFileUploader"] {
    border-radius: 8px;
}

/* Section divider */
.section-divider {
    height: 1px;
    background: #f0f0f0;
    margin: 20px 0;
}
</style>
""", unsafe_allow_html=True)

# ---------------------------------------------------------------------------
# Session state initialisation
# ---------------------------------------------------------------------------
def init_state():
    defaults = {
        "page":             "setup",
        "project_name":     "",
        "project_notes":    "",
        "atlas_resize":     "25% — QuickNII standard",
        "region_level":     "Leaf level (~1300 regions)",
        "background_method":"Modal intensity per section",
        "animals":          [],
        "injection_regions":["MOs", "MOp"],
        "norm_method":      "Mean across injection site sections (recommended)",
        "sparse_threshold": 3,
        "run_done":         False,
        "section_df":       None,
        "animal_df":        None,
        "cohort_df":        None,
        "run_log":          [],
    }
    for k, v in defaults.items():
        if k not in st.session_state:
            st.session_state[k] = v

init_state()

# ---------------------------------------------------------------------------
# Sidebar navigation
# ---------------------------------------------------------------------------
PAGES = [
    ("setup",   "Project setup"),
    ("data",    "Load data"),
    ("config",  "Configure"),
    ("run",     "Run pipeline"),
    ("results", "Results"),
    ("viz",     "Visualize"),
    ("export",  "Export"),
]

with st.sidebar:
    st.markdown('<div class="sidebar-title">🧠 TractQuant</div>', unsafe_allow_html=True)
    st.markdown('<div class="sidebar-sub">Anterograde tracing pipeline</div>', unsafe_allow_html=True)

    st.markdown("**Navigation**")
    for page_id, page_label in PAGES:
        is_active = st.session_state.page == page_id
        prefix = "→ " if is_active else "   "
        if st.button(
            f"{prefix}{page_label}",
            key=f"nav_{page_id}",
            use_container_width=True,
            type="primary" if is_active else "secondary",
        ):
            st.session_state.page = page_id
            st.rerun()

    st.divider()
    st.caption(f"Animals loaded: **{len(st.session_state.animals)}**")
    st.caption("Atlas: ABA CCFv3 2017 · 25µm")
    st.caption("Cache: ~/.tractquant/atlas/")

# ---------------------------------------------------------------------------
# Helper: render animal card HTML
# ---------------------------------------------------------------------------
def animal_card_html(an: dict) -> str:
    c2 = "TRAP-tTA" if an["ch1_construct"] == "TRAP-Cre" else "TRAP-Cre"
    c1cls = "tag-purple" if an["ch1_construct"] == "TRAP-Cre" else "tag-green"
    c2cls = "tag-green"  if c2 == "TRAP-tTA" else "tag-purple"
    sex_cls = "tag-pink" if an["sex"] == "Female" else "tag-blue"

    ch1_key = " · ".join(filter(None, [
        an["ch1_construct"], an["ch1_stress"], an["ch1_duration"], an["ch1_timepoint"]
    ]))
    ch2_key = " · ".join(filter(None, [
        c2, an["ch2_stress"], an["ch2_duration"], an["ch2_timepoint"]
    ]))

    same = (an["ch1_stress"] == an["ch2_stress"] and
            an["ch1_duration"] == an["ch2_duration"] and
            an["ch1_timepoint"] == an["ch2_timepoint"])
    stress_tag = (
        '<span class="tag tag-green">same stress both channels</span>'
        if same else
        '<span class="tag tag-amber">different stress per channel</span>'
    )

    return f"""
    <div class="animal-row">
        <div style="display:flex;align-items:center;gap:8px;margin-bottom:6px">
            <span class="animal-id">{an['id']}</span>
            <span class="tag {sex_cls}">{an['sex']}</span>
            {stress_tag}
        </div>
        <div style="margin-bottom:3px">
            <span class="tag {c1cls}">Ch1 · {an['ch1_construct']}</span>
            <span class="compound-key">{ch1_key}</span>
        </div>
        <div>
            <span class="tag {c2cls}">Ch2 · {c2}</span>
            <span class="compound-key">{ch2_key}</span>
        </div>
    </div>"""


# ---------------------------------------------------------------------------
# Page: Setup
# ---------------------------------------------------------------------------
def page_setup():
    st.markdown("""
    <div class="page-header">
        <h2>Project setup</h2>
        <p>Name your experiment and configure global analysis settings</p>
    </div>""", unsafe_allow_html=True)

    with st.form("setup_form"):
        st.markdown("#### Project details")
        col1, col2 = st.columns(2)
        with col1:
            name = st.text_input("Project name",
                value=st.session_state.project_name,
                placeholder="e.g. Motor cortex stress series 2025")
        with col2:
            notes = st.text_input("Notes (optional)",
                value=st.session_state.project_notes,
                placeholder="Any additional context")

        st.markdown('<div class="section-divider"></div>', unsafe_allow_html=True)
        st.markdown("#### Atlas and imaging")
        col1, col2 = st.columns(2)
        with col1:
            atlas_resize = st.selectbox("Image resize factor (Nutil)",
                ["25% — QuickNII standard", "50%", "100% full res"],
                index=["25% — QuickNII standard","50%","100% full res"].index(
                    st.session_state.atlas_resize))
        with col2:
            region_level = st.selectbox("Region hierarchy level",
                ["Leaf level (~1300 regions)",
                 "Summary structures (~300 regions)",
                 "Major divisions (12 regions)"],
                index=["Leaf level (~1300 regions)",
                       "Summary structures (~300 regions)",
                       "Major divisions (12 regions)"].index(
                    st.session_state.region_level))

        st.markdown('<div class="section-divider"></div>', unsafe_allow_html=True)
        st.markdown("#### Background estimation")
        bg_method = st.radio("Method",
            ["Modal intensity per section (recommended)",
             "Designated background atlas regions",
             "Border ring per region"],
            index=0)

        st.markdown('<div class="section-divider"></div>', unsafe_allow_html=True)

        st.markdown("""
        <div class="info-banner">
        <strong>Construct space:</strong> All outputs use biological labels
        (TRAP-Cre / TRAP-tTA), not physical channel labels.
        The pipeline resolves which channel carries which construct per animal —
        counterbalancing is handled automatically.
        </div>""", unsafe_allow_html=True)

        submitted = st.form_submit_button("Save and continue →",
            type="primary", use_container_width=True)

        if submitted:
            st.session_state.project_name  = name
            st.session_state.project_notes = notes
            st.session_state.atlas_resize  = atlas_resize
            st.session_state.region_level  = region_level
            st.session_state.page          = "data"
            st.rerun()


# ---------------------------------------------------------------------------
# Page: Load data
# ---------------------------------------------------------------------------
def page_data():
    st.markdown("""
    <div class="page-header">
        <h2>Load data</h2>
        <p>Add animals — stress conditions and construct assignment set independently per channel</p>
    </div>""", unsafe_allow_html=True)

    # Counterbalance check
    if len(st.session_state.animals) >= 2:
        cre_count = sum(1 for a in st.session_state.animals if a["ch1_construct"] == "TRAP-Cre")
        tta_count = len(st.session_state.animals) - cre_count
        if cre_count == 0 or tta_count == 0:
            st.markdown("""<div class="warn-banner">
            ⚠ All animals have the same Ch1 assignment — counterbalancing not yet present.
            </div>""", unsafe_allow_html=True)
        else:
            st.markdown(f"""<div class="ok-banner">
            ✓ Counterbalancing detected — {cre_count} animal(s) with Ch1=TRAP-Cre,
            {tta_count} with Ch1=TRAP-tTA
            </div>""", unsafe_allow_html=True)

    # ---------------------------------------------------------------------------
    # Add animal form
    # ---------------------------------------------------------------------------
    with st.expander("➕  Add animal", expanded=len(st.session_state.animals) == 0):
        with st.form("add_animal_form", clear_on_submit=True):
            st.markdown("**Animal identity**")
            col1, col2 = st.columns(2)
            with col1:
                animal_id = st.text_input("Animal ID *",
                    placeholder="e.g. F5ex")
            with col2:
                sex = st.selectbox("Sex", ["Male", "Female"])

            st.markdown('<div class="section-divider"></div>', unsafe_allow_html=True)

            # CH1 block
            st.markdown('<div class="ch-block-cre">', unsafe_allow_html=True)
            st.markdown('<span class="ch-label-cre">Channel 1</span>', unsafe_allow_html=True)
            col1, col2 = st.columns([1, 3])
            with col1:
                ch1_construct = st.selectbox("Construct", ["TRAP-Cre", "TRAP-tTA"],
                    key="form_ch1_construct")
            ch2_construct = "TRAP-tTA" if ch1_construct == "TRAP-Cre" else "TRAP-Cre"

            col1, col2, col3 = st.columns(3)
            with col1:
                ch1_stress = st.text_input("Stress type",
                    placeholder="e.g. Restraint, CRS, None",
                    key="ch1_stress")
            with col2:
                ch1_duration = st.selectbox("Duration",
                    ["Acute", "Sub-chronic", "Chronic", "None"],
                    key="ch1_duration")
            with col3:
                ch1_timepoint = st.text_input("TRAP timepoint",
                    placeholder="e.g. Day 1, Home cage",
                    key="ch1_timepoint")
            st.markdown('</div>', unsafe_allow_html=True)

            # CH2 block
            st.markdown(f"""
            <div class="ch-block-tta">
                <span class="ch-label-tta">Channel 2 — auto: {ch2_construct}</span>
            </div>""", unsafe_allow_html=True)

            copy_stress = st.checkbox("Same stress as Ch1 (copy across)",
                value=False, key="copy_stress")

            if not copy_stress:
                col1, col2, col3 = st.columns(3)
                with col1:
                    ch2_stress = st.text_input("Stress type",
                        placeholder="e.g. Social defeat, CRS, None",
                        key="ch2_stress")
                with col2:
                    ch2_duration = st.selectbox("Duration",
                        ["Acute", "Sub-chronic", "Chronic", "None"],
                        key="ch2_duration")
                with col3:
                    ch2_timepoint = st.text_input("TRAP timepoint",
                        placeholder="e.g. Day 14, Home cage",
                        key="ch2_timepoint")
            else:
                ch2_stress    = ch1_stress
                ch2_duration  = ch1_duration
                ch2_timepoint = ch1_timepoint

            st.markdown('<div class="section-divider"></div>', unsafe_allow_html=True)
            st.markdown("**Files**")

            col1, col2, col3 = st.columns(3)
            with col1:
                json_file = st.file_uploader("VisualAlign JSON *",
                    type=["json"], key="upload_json",
                    help="The .json output from VisualAlign")
            with col2:
                ch1_dir = st.text_input("Ch1 image folder path *",
                    placeholder="/path/to/ch1/",
                    help="Folder containing ch1_001.tif … ch1_040.tif")
            with col3:
                ch2_dir = st.text_input("Ch2 image folder path *",
                    placeholder="/path/to/ch2/",
                    help="Folder containing ch2_001.tif … ch2_040.tif")

            st.markdown("*All fields marked * are required*")

            submitted = st.form_submit_button("Add animal",
                type="primary", use_container_width=True)

            if submitted:
                errors = []
                if not animal_id.strip():
                    errors.append("Animal ID is required")
                if not json_file:
                    errors.append("VisualAlign JSON is required")
                if not ch1_dir.strip():
                    errors.append("Ch1 folder path is required")
                if not ch2_dir.strip():
                    errors.append("Ch2 folder path is required")
                if animal_id in [a["id"] for a in st.session_state.animals]:
                    errors.append(f"Animal ID '{animal_id}' already exists")

                if errors:
                    for e in errors:
                        st.error(e)
                else:
                    # Save uploaded JSON to temp location
                    json_dir = Path(".tractquant_tmp") / animal_id.strip()
                    json_dir.mkdir(parents=True, exist_ok=True)
                    json_path = json_dir / json_file.name
                    json_path.write_bytes(json_file.getvalue())

                    # Count sections in JSON
                    try:
                        j = json.loads(json_file.getvalue())
                        n_sections = len(j.get("slices", []))
                    except Exception:
                        n_sections = "?"

                    st.session_state.animals.append({
                        "id":             animal_id.strip(),
                        "sex":            sex,
                        "ch1_construct":  ch1_construct,
                        "ch1_stress":     ch1_stress,
                        "ch1_duration":   ch1_duration,
                        "ch1_timepoint":  ch1_timepoint,
                        "ch2_stress":     ch2_stress,
                        "ch2_duration":   ch2_duration,
                        "ch2_timepoint":  ch2_timepoint,
                        "json_path":      str(json_path),
                        "ch1_dir":        ch1_dir.strip(),
                        "ch2_dir":        ch2_dir.strip(),
                        "n_sections":     n_sections,
                    })
                    st.success(f"✓ Added {animal_id.strip()} — {n_sections} sections detected")
                    st.rerun()

    st.markdown('<div class="section-divider"></div>', unsafe_allow_html=True)

    # ---------------------------------------------------------------------------
    # Animal list
    # ---------------------------------------------------------------------------
    n = len(st.session_state.animals)
    st.markdown(f"#### Animals in project ({n})")

    if n == 0:
        st.info("No animals added yet — use the form above to add your first animal.")
    else:
        for i, an in enumerate(st.session_state.animals):
            col1, col2 = st.columns([9, 1])
            with col1:
                st.markdown(animal_card_html(an), unsafe_allow_html=True)
            with col2:
                st.markdown("<div style='margin-top:18px'></div>", unsafe_allow_html=True)
                if st.button("✕", key=f"remove_{i}",
                             help=f"Remove {an['id']}"):
                    st.session_state.animals.pop(i)
                    st.rerun()

        st.markdown('<div class="section-divider"></div>', unsafe_allow_html=True)
        if st.button("Next: Configure →", type="primary"):
            st.session_state.page = "config"
            st.rerun()


# ---------------------------------------------------------------------------
# Page: Configure
# ---------------------------------------------------------------------------
def page_config():
    st.markdown("""
    <div class="page-header">
        <h2>Configure analysis</h2>
        <p>Set injection site regions, normalization method, and region filtering</p>
    </div>""", unsafe_allow_html=True)

    with st.form("config_form"):
        st.markdown("#### Injection site regions")
        st.caption(
            "These CCF regions define your normalization baseline. "
            "Mean background-corrected intensity across all sections "
            "containing these regions — computed independently per channel."
        )

        inj_text = st.text_input(
            "CCF acronyms (comma-separated)",
            value=", ".join(st.session_state.injection_regions),
            placeholder="e.g. MOs, MOp, VISp",
            help="Use official Allen CCF acronyms. Case-sensitive."
        )

        st.markdown('<div class="section-divider"></div>', unsafe_allow_html=True)
        st.markdown("#### Normalization")
        norm = st.radio("Method", [
            "Mean across injection site sections (recommended)",
            "Sum across injection site sections",
        ])
        st.session_state.norm_method = norm
        st.markdown("""
        <div class="info-banner" style="margin-top:8px">
        Each channel is normalized independently — TRAP-Cre to its own injection 
        intensity, TRAP-tTA to its own. This ensures the ratio reflects true 
        differential targeting, not labeling efficiency differences.
        </div>""", unsafe_allow_html=True)

        st.markdown('<div class="section-divider"></div>', unsafe_allow_html=True)
        st.markdown("#### Grouping variables for statistics")
        st.caption("Used when computing group means and z-scores across animals.")
        col1, col2 = st.columns(2)
        with col1:
            st.checkbox("Viral construct (TRAP-Cre vs TRAP-tTA)", value=True, disabled=True)
            st.checkbox("Stress type", value=True)
            st.checkbox("Stress duration", value=True)
        with col2:
            st.checkbox("TRAP timepoint", value=True)
            st.checkbox("Sex", value=True)

        st.markdown('<div class="section-divider"></div>', unsafe_allow_html=True)
        st.markdown("#### Sparse region handling")
        col1, col2 = st.columns(2)
        with col1:
            min_sections = st.number_input(
                "Minimum sections to include region", min_value=1, value=1)
        with col2:
            sparse_threshold = st.number_input(
                "Flag regions sampled by fewer than N sections",
                min_value=1, value=3)

        st.markdown('<div class="section-divider"></div>', unsafe_allow_html=True)
        submitted = st.form_submit_button("Save and continue →",
            type="primary", use_container_width=True)

        if submitted:
            regions = [r.strip() for r in inj_text.split(",") if r.strip()]
            if not regions:
                st.error("Please enter at least one injection region")
            else:
                st.session_state.injection_regions = regions
                st.session_state.norm_method       = norm
                st.session_state.sparse_threshold  = int(sparse_threshold)
                st.session_state.page              = "run"
                st.rerun()


# ---------------------------------------------------------------------------
# Page: Run
# ---------------------------------------------------------------------------
def page_run():
    st.markdown("""
    <div class="page-header">
        <h2>Run pipeline</h2>
        <p>Process all animals — atlas downloads automatically if not already cached</p>
    </div>""", unsafe_allow_html=True)

    n_animals  = len(st.session_state.animals)
    n_sections = sum(a.get("n_sections", 40) for a in st.session_state.animals
                     if isinstance(a.get("n_sections"), int))

    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Animals",  n_animals)
    col2.metric("Sections", n_sections)
    col3.metric("Channels", 2)
    col4.metric("Status", "✓ Done" if st.session_state.run_done else "Ready")

    st.markdown('<div class="section-divider"></div>', unsafe_allow_html=True)

    if n_animals == 0:
        st.warning("No animals loaded. Go back to 'Load data' and add animals first.")
        return

    if not st.session_state.run_done:
        st.markdown("#### Pre-flight check")
        checks = []
        checks.append(("Project name set",     bool(st.session_state.project_name)))
        checks.append(("Animals loaded",        n_animals > 0))
        checks.append(("Injection regions set", len(st.session_state.injection_regions) > 0))

        # Counterbalance check
        cre_n = sum(1 for a in st.session_state.animals if a["ch1_construct"] == "TRAP-Cre")
        tta_n = n_animals - cre_n
        checks.append(("Counterbalancing present",
                        cre_n > 0 and tta_n > 0 or n_animals == 1))

        all_ok = True
        for label, ok in checks:
            icon = "✅" if ok else "❌"
            st.markdown(f"{icon} {label}")
            if not ok:
                all_ok = False

        if not all_ok:
            st.error("Fix the issues above before running.")
            return

        st.markdown('<div class="section-divider"></div>', unsafe_allow_html=True)

        if st.button("🚀  Run TractQuant pipeline",
                     type="primary", use_container_width=True):
            _run_pipeline_ui()
    else:
        st.markdown("""<div class="ok-banner">
        ✓ Pipeline complete — navigate to Results or Visualize to explore your data.
        </div>""", unsafe_allow_html=True)

        # Show log
        with st.expander("View pipeline log"):
            log_html = "<br>".join(st.session_state.run_log)
            st.markdown(f'<div class="log-box">{log_html}</div>',
                        unsafe_allow_html=True)

        col1, col2 = st.columns(2)
        with col1:
            if st.button("View Results →", type="primary", use_container_width=True):
                st.session_state.page = "results"
                st.rerun()
        with col2:
            if st.button("Visualize →", use_container_width=True):
                st.session_state.page = "viz"
                st.rerun()


def _run_pipeline_ui():
    """Run the actual pipeline with live progress display."""
    log_placeholder  = st.empty()
    prog_placeholder = st.progress(0)
    status_placeholder = st.empty()

    log_lines = []

    def log(msg: str, kind: str = "ok"):
        css = {"ok": "log-ok", "warn": "log-warn", "err": "log-err"}.get(kind, "")
        prefix = {"ok": "✓", "warn": "⚠", "err": "✗"}.get(kind, "")
        line = f'<span class="{css}">{prefix} {msg}</span>'
        log_lines.append(line)
        log_html = "<br>".join(log_lines)
        log_placeholder.markdown(
            f'<div class="log-box">{log_html}</div>',
            unsafe_allow_html=True
        )

    def update_progress(pct: int, msg: str = ""):
        prog_placeholder.progress(pct / 100)
        if msg:
            status_placeholder.caption(msg)

    try:
        # Import pipeline modules
        log("Importing pipeline modules...")
        update_progress(2)

        sys.path.insert(0, str(Path(__file__).parent))
        from atlas_manager import get_atlas_volume, get_region_csv
        from pipeline import (
            CCFAtlas, AnimalConfig, ChannelMeta,
            parse_visual_align, measure_section,
            load_section_images, expand_region_acronyms, aggregate_cohort,
            export_brainglobe, save_alignment_overlay,
        )

        # Atlas
        log("Checking atlas cache...")
        update_progress(5, "Loading atlas...")
        atlas_path   = get_atlas_volume()
        region_csv   = get_region_csv()
        log(f"Atlas ready: {Path(atlas_path).name}")
        update_progress(15)

        atlas = CCFAtlas(str(atlas_path), str(region_csv))
        log(f"Atlas loaded — shape: {atlas.shape}")
        update_progress(20)

        resize_map = {
            "25% — QuickNII standard": 0.25,
            "50%": 0.5, "100% full res": 1.0,
        }
        resize_factor = resize_map.get(st.session_state.atlas_resize, 0.25)
        log(f"Resize factor: {resize_factor} (image coords × {1/resize_factor:.0f} → atlas)")

        output_dir = Path("results") / st.session_state.project_name.replace(" ", "_")
        output_dir.mkdir(parents=True, exist_ok=True)

        animal_dfs   = []
        n_animals    = len(st.session_state.animals)
        inj_regions  = st.session_state.injection_regions

        for a_idx, an in enumerate(st.session_state.animals):
            base_pct = 20 + int(a_idx / n_animals * 65)
            update_progress(base_pct, f"Processing {an['id']}...")
            log(f"Animal: {an['id']} ({an['sex']})")

            # Build AnimalConfig
            ch1_meta = ChannelMeta(
                construct=an["ch1_construct"],
                stress_type=an["ch1_stress"],
                stress_duration=an["ch1_duration"],
                timepoint=an["ch1_timepoint"],
            )
            c2 = "TRAP-tTA" if an["ch1_construct"] == "TRAP-Cre" else "TRAP-Cre"
            ch2_meta = ChannelMeta(
                construct=c2,
                stress_type=an["ch2_stress"],
                stress_duration=an["ch2_duration"],
                timepoint=an["ch2_timepoint"],
            )

            log(f"  Ch1: {ch1_meta.construct} — {ch1_meta.compound_key}")
            log(f"  Ch2: {ch2_meta.construct} — {ch2_meta.compound_key}")

            # Parse JSON
            sections = parse_visual_align(Path(an["json_path"]))
            no_marker = [s["nr"] for s in sections.values() if not s.get("has_markers")]
            if no_marker:
                log(f"  Sections without markers: {no_marker} — anchoring only, included", "warn")

            # Process sections
            section_rows = []
            ch1_dir = Path(an["ch1_dir"])
            ch2_dir = Path(an["ch2_dir"])

            for s_idx, (filename, sec) in enumerate(sections.items()):
                nr = sec["nr"]
                try:
                    ch1_img, ch2_img = load_section_images(ch1_dir, ch2_dir, nr, filename)
                    sec_df = measure_section(
                        sec, ch1_img, ch2_img, atlas, resize_factor
                    )
                    qc = save_alignment_overlay(
                        sec, ch1_img, ch2_img, atlas, output_dir, an["id"], nr, filename, resize_factor
                    )
                    if sec_df.empty:
                        atlas_match = sec_df.attrs.get("atlas_match_percent", 0.0)
                        log(f"  s{nr:03d}: 0 regions, atlas match {atlas_match:.2f}% — skipped", "warn")
                        continue
                    sec_df["section_nr"]   = nr
                    sec_df["section_file"] = filename
                    section_rows.append(sec_df)
                    ref_txt = ""
                    if qc.get("ch1_reference_iou") is not None:
                        ref_txt = f", ref IoU {qc['ch1_reference_iou']:.2f}"
                    marker_txt = (
                        f", markers {qc.get('marker_count', 0)}"
                        if qc.get("use_markers") else ", anchoring only"
                    )
                    log(
                        f"  s{nr:03d}: QC overlay saved, mask {qc['ch1_atlas_percent']:.1f}%/{qc['ch2_atlas_percent']:.1f}%{marker_txt}{ref_txt}"
                    )
                except Exception as e:
                    log(f"  s{nr:03d}: skipped — {e}", "warn")

                section_pct = base_pct + int((s_idx+1) / len(sections) * (65 // n_animals))
                update_progress(min(section_pct, 84))

            if not section_rows:
                log(f"  No sections processed for {an['id']}", "err")
                continue

            import pandas as pd
            section_df = pd.concat(section_rows, ignore_index=True)

            # Injection baseline
            inj_regions_expanded = expand_region_acronyms(
                section_df["region_acronym"], inj_regions
            )
            inj_mask = section_df["region_acronym"].isin(inj_regions_expanded)
            if inj_mask.sum() == 0:
                log(f"  Injection regions not found: {inj_regions}", "err")
                continue

            inj_df = section_df[inj_mask]
            if st.session_state.norm_method.lower().startswith("sum"):
                ch1_baseline = inj_df["ch1_corrected_sum"].sum() / inj_df["n_pixels"].sum()
                ch2_baseline = inj_df["ch2_corrected_sum"].sum() / inj_df["n_pixels"].sum()
                norm_label = "sum across injection sections, mean intensity per pixel"
            else:
                inj_by_section = inj_df.groupby("section_nr")
                ch1_baseline = (
                    inj_by_section["ch1_corrected_sum"].sum()
                    / inj_by_section["n_pixels"].sum()
                ).mean()
                ch2_baseline = (
                    inj_by_section["ch2_corrected_sum"].sum()
                    / inj_by_section["n_pixels"].sum()
                ).mean()
                norm_label = "mean across injection sections, mean intensity per pixel"
            if ch1_baseline <= 0 or ch2_baseline <= 0:
                log("  Injection normalization baseline is zero", "err")
                continue
            log(f"  Injection regions: {', '.join(inj_regions_expanded)}")
            log(f"  Normalization: {norm_label}")
            log(f"  Injection baseline — ch1: {ch1_baseline:.0f}, ch2: {ch2_baseline:.0f}")

            # Aggregate
            agg = section_df.groupby(
                ["region_id","region_acronym","region_name","division"]
            ).agg(
                ch1_corrected_sum=("ch1_corrected_sum","sum"),
                ch2_corrected_sum=("ch2_corrected_sum","sum"),
                total_pixels=("n_pixels","sum"),
                n_sections=("section_nr","nunique"),
            ).reset_index()

            agg["ch1_corrected_mean"] = agg["ch1_corrected_sum"] / agg["total_pixels"]
            agg["ch2_corrected_mean"] = agg["ch2_corrected_sum"] / agg["total_pixels"]
            agg["ch1_fraction"] = agg["ch1_corrected_mean"] / ch1_baseline
            agg["ch2_fraction"] = agg["ch2_corrected_mean"] / ch2_baseline
            agg["sparse_flag"]  = agg["n_sections"] < st.session_state.sparse_threshold

            # Resolve to construct space
            if an["ch1_construct"] == "TRAP-Cre":
                cre_frac, tta_frac = agg["ch1_fraction"], agg["ch2_fraction"]
                cre_meta, tta_meta = ch1_meta, ch2_meta
            else:
                cre_frac, tta_frac = agg["ch2_fraction"], agg["ch1_fraction"]
                cre_meta, tta_meta = ch2_meta, ch1_meta

            animal_df = agg[["region_id","region_acronym","region_name",
                              "division","total_pixels","n_sections","sparse_flag"]].copy()
            animal_df["animal_id"]            = an["id"]
            animal_df["sex"]                  = an["sex"]
            animal_df["ch1_construct"]        = an["ch1_construct"]
            animal_df["trapcre_stress_type"]  = cre_meta.stress_type
            animal_df["trapcre_stress_dur"]   = cre_meta.stress_duration
            animal_df["trapcre_timepoint"]    = cre_meta.timepoint
            animal_df["trapcre_compound_key"] = cre_meta.compound_key
            animal_df["traptta_stress_type"]  = tta_meta.stress_type
            animal_df["traptta_stress_dur"]   = tta_meta.stress_duration
            animal_df["traptta_timepoint"]    = tta_meta.timepoint
            animal_df["traptta_compound_key"] = tta_meta.compound_key
            animal_df["trapcre_fraction"]     = cre_frac.values
            animal_df["traptta_fraction"]     = tta_frac.values
            animal_df["ratio"] = np.where(
                animal_df["traptta_fraction"] > 0,
                animal_df["trapcre_fraction"] / animal_df["traptta_fraction"], np.nan)

            # Save
            out = output_dir / an["id"]
            out.mkdir(parents=True, exist_ok=True)
            section_df.to_csv(out / "section_measurements.csv", index=False)
            animal_df.to_csv(out / "animal_projections.csv", index=False)
            log(f"  Saved → {out}/")
            animal_dfs.append(animal_df)

        update_progress(88, "Aggregating cohort...")

        if not animal_dfs:
            log("No animals processed successfully", "err")
            return

        all_df = pd.concat(animal_dfs, ignore_index=True)
        cohort_df = aggregate_cohort(animal_dfs)
        cohort_df.to_csv(output_dir / "cohort_summary.csv", index=False)
        log(f"Cohort summary: {len(cohort_df)} group × region combinations")

        update_progress(95, "Exporting BrainGlobe files...")
        bg_dir = output_dir / "brainglobe"
        bg_dir.mkdir(exist_ok=True)
        avg = all_df.groupby("region_acronym").agg(
            trapcre_fraction=("trapcre_fraction","mean"),
            traptta_fraction=("traptta_fraction","mean"),
            ratio=("ratio","mean"),
        ).reset_index()
        for col, fname in [
            ("trapcre_fraction","trapcre_projection.json"),
            ("traptta_fraction","traptta_projection.json"),
            ("ratio","ratio.json"),
        ]:
            d = avg.dropna(subset=[col]).set_index("region_acronym")[col].to_dict()
            with open(bg_dir / fname, "w") as f:
                json.dump(d, f, indent=2)
        log(f"BrainGlobe exports written to {bg_dir}/")

        update_progress(100, "Complete!")
        log("Pipeline complete ✓")

        st.session_state.run_done  = True
        st.session_state.run_log   = log_lines
        st.session_state.animal_df = all_df
        st.session_state.cohort_df = cohort_df

        time.sleep(0.5)
        st.rerun()

    except ImportError as e:
        log(f"Import error: {e} — make sure you're running inside the tractquant conda env", "err")
    except Exception as e:
        log(f"Pipeline error: {e}", "err")
        import traceback
        log(traceback.format_exc(), "err")


# ---------------------------------------------------------------------------
# Page: Results
# ---------------------------------------------------------------------------
def page_results():
    st.markdown("""
    <div class="page-header">
        <h2>Projection data</h2>
        <p>Per-region projection fractions in construct space</p>
    </div>""", unsafe_allow_html=True)

    if st.session_state.animal_df is None:
        st.warning("No results yet — run the pipeline first.")
        return

    df = st.session_state.animal_df.copy()

    # Summary metrics
    n_regions = df["region_acronym"].nunique()
    copro = ((df["trapcre_fraction"] > 0.05) &
             (df["traptta_fraction"] > 0.05)).sum()
    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Regions quantified", n_regions)
    col2.metric("Co-projection regions", copro)
    col3.metric("Animals", df["animal_id"].nunique())
    col4.metric("Sections (total)", df["n_sections"].sum())

    st.markdown("""<div class="info-banner" style="margin-top:16px">
    Results shown in <strong>construct space</strong> — TRAP-Cre and TRAP-tTA 
    columns correctly assigned regardless of which physical channel carried each 
    construct for each animal.
    </div>""", unsafe_allow_html=True)

    # Filters
    col1, col2, col3 = st.columns(3)
    with col1:
        animal_filter = st.multiselect(
            "Filter by animal",
            options=sorted(df["animal_id"].unique()),
            default=sorted(df["animal_id"].unique())
        )
    with col2:
        division_filter = st.multiselect(
            "Filter by division",
            options=sorted(df["division"].dropna().unique()),
            default=sorted(df["division"].dropna().unique())
        )
    with col3:
        min_frac = st.slider("Min projection fraction (%)", 0, 20, 0) / 100

    filtered = df[
        df["animal_id"].isin(animal_filter) &
        df["division"].isin(division_filter) &
        ((df["trapcre_fraction"] >= min_frac) | (df["traptta_fraction"] >= min_frac))
    ].copy()

    # Display table
    display = filtered[[
        "region_acronym","region_name","division",
        "trapcre_fraction","traptta_fraction","ratio",
        "n_sections","sparse_flag","animal_id","sex",
        "trapcre_compound_key","traptta_compound_key",
    ]].copy()
    display["trapcre_%"] = (display["trapcre_fraction"] * 100).round(2)
    display["traptta_%"] = (display["traptta_fraction"] * 100).round(2)
    display["ratio"]     = display["ratio"].round(3)

    show_cols = ["region_acronym","region_name","division",
                 "trapcre_%","traptta_%","ratio",
                 "n_sections","sparse_flag","animal_id","sex"]

    st.dataframe(
        display[show_cols].sort_values("trapcre_%", ascending=False),
        use_container_width=True,
        height=420,
    )
    st.caption(f"Showing {len(filtered)} region × animal rows")


# ---------------------------------------------------------------------------
# Page: Visualize
# ---------------------------------------------------------------------------
def page_viz():
    st.markdown("""
    <div class="page-header">
        <h2>Visualize</h2>
        <p>Brain-wide projection patterns and TRAP-Cre vs TRAP-tTA comparison</p>
    </div>""", unsafe_allow_html=True)

    if st.session_state.animal_df is None:
        st.warning("No results yet — run the pipeline first.")
        return

    try:
        import matplotlib.pyplot as plt
        import matplotlib
        matplotlib.use("Agg")
    except ImportError:
        st.error("matplotlib not installed in this environment.")
        return

    df = st.session_state.animal_df.copy()

    tab1, tab2, tab3 = st.tabs(["📊 Bar chart", "⚬ Scatter", "🟥 Heatmap"])

    # ---- Bar chart ----
    with tab1:
        col1, col2 = st.columns([2, 1])
        with col2:
            n_top     = st.slider("Top N regions", 5, 30, 15)
            sort_by   = st.radio("Sort by", ["TRAP-Cre", "TRAP-tTA"])
            group_by_animal = st.checkbox("Average across animals", value=True)

        sort_col = "trapcre_fraction" if sort_by == "TRAP-Cre" else "traptta_fraction"

        if group_by_animal:
            plot_df = df.groupby("region_acronym").agg(
                trapcre_fraction=("trapcre_fraction","mean"),
                traptta_fraction=("traptta_fraction","mean"),
                division=("division","first"),
            ).reset_index()
        else:
            plot_df = df.copy()

        plot_df = plot_df.nlargest(n_top, sort_col).sort_values(sort_col, ascending=True)

        with col1:
            fig, ax = plt.subplots(figsize=(7, max(4, n_top * 0.42)))
            y = np.arange(len(plot_df)); h = 0.35
            ax.barh(y + h/2, plot_df["trapcre_fraction"]*100, h,
                    color="#534AB7", alpha=0.85, label="TRAP-Cre")
            ax.barh(y - h/2, plot_df["traptta_fraction"]*100, h,
                    color="#1D9E75", alpha=0.85, label="TRAP-tTA")
            ax.set_yticks(y)
            ax.set_yticklabels(plot_df["region_acronym"], fontsize=9)
            ax.set_xlabel("Projection fraction (% of injection site)")
            ax.legend(frameon=False)
            ax.spines[["top","right"]].set_visible(False)
            plt.tight_layout()
            st.pyplot(fig)
            plt.close()

    # ---- Scatter ----
    with tab2:
        col1, col2 = st.columns([2, 1])
        with col2:
            min_frac = st.slider("Min fraction (%)", 0, 10, 1) / 100
            label_n  = st.slider("Label top N regions", 0, 20, 10)
            group_scatter = st.checkbox("Average across animals", value=True, key="sc_avg")

        if group_scatter:
            sc_df = df.groupby(["region_acronym","division"]).agg(
                trapcre_fraction=("trapcre_fraction","mean"),
                traptta_fraction=("traptta_fraction","mean"),
            ).reset_index()
        else:
            sc_df = df.copy()

        sc_df = sc_df[
            (sc_df["trapcre_fraction"] > min_frac) |
            (sc_df["traptta_fraction"] > min_frac)
        ]

        DIV_COLORS = {
            "Isocortex":"#534AB7","Striatum":"#1D9E75","Thalamus":"#BA7517",
            "Midbrain":"#D85A30","Pallidum":"#888780","Brainstem":"#3C3489",
            "Hypothalamus":"#639922","Hindbrain":"#D4537E",
        }

        with col1:
            fig, ax = plt.subplots(figsize=(6, 6))
            lim = max(sc_df["trapcre_fraction"].max(),
                      sc_df["traptta_fraction"].max()) * 1.1 * 100
            ax.plot([0,lim],[0,lim], color="#ccc", lw=0.8, ls="--", zorder=0)

            for div, grp in sc_df.groupby("division"):
                ax.scatter(grp["trapcre_fraction"]*100,
                           grp["traptta_fraction"]*100,
                           color=DIV_COLORS.get(div,"#aaa"),
                           alpha=0.75, s=40, label=div, zorder=2)

            top = sc_df.nlargest(label_n, "trapcre_fraction")
            for _, row in top.iterrows():
                ax.annotate(row["region_acronym"],
                            xy=(row["trapcre_fraction"]*100, row["traptta_fraction"]*100),
                            xytext=(3,3), textcoords="offset points",
                            fontsize=7, color="#444")

            ax.set_xlabel("TRAP-Cre projection fraction (%)")
            ax.set_ylabel("TRAP-tTA projection fraction (%)")
            ax.set_xlim(0, lim); ax.set_ylim(0, lim)
            ax.text(0.97, 0.03, "Cre dominant", transform=ax.transAxes,
                    ha="right", va="bottom", fontsize=8, color="#999")
            ax.text(0.03, 0.97, "tTA dominant", transform=ax.transAxes,
                    ha="left", va="top", fontsize=8, color="#999")
            ax.legend(frameon=False, fontsize=8, bbox_to_anchor=(1.01,1), loc="upper left")
            ax.spines[["top","right"]].set_visible(False)
            plt.tight_layout()
            st.pyplot(fig)
            plt.close()

    # ---- Heatmap ----
    with tab3:
        col1, col2 = st.columns([2,1])
        with col2:
            construct  = st.radio("Construct", ["TRAP-Cre","TRAP-tTA"])
            n_heat     = st.slider("Top N regions", 10, 50, 25, key="heat_n")

        col_key = "trapcre_fraction" if construct == "TRAP-Cre" else "traptta_fraction"
        color   = "#534AB7" if construct == "TRAP-Cre" else "#1D9E75"

        pivot = df.pivot_table(
            index="region_acronym", columns="animal_id",
            values=col_key, aggfunc="mean"
        )
        top_regs = pivot.mean(axis=1).nlargest(n_heat).index
        pivot    = pivot.loc[top_regs]

        try:
            import seaborn as sns
            with col1:
                fig, ax = plt.subplots(figsize=(max(5, len(pivot.columns)*0.7),
                                                n_heat * 0.35 + 1))
                cmap = sns.light_palette(color, as_cmap=True)
                sns.heatmap(pivot*100, ax=ax, cmap=cmap,
                            linewidths=0.3, linecolor="#eee",
                            cbar_kws={"label":"Projection fraction (%)","shrink":0.6},
                            fmt=".0f",
                            annot=len(pivot.columns) <= 15,
                            annot_kws={"size":7})
                ax.set_title(f"{construct} — region × animal")
                ax.tick_params(axis="x", labelrotation=45, labelsize=8)
                ax.tick_params(axis="y", labelsize=8)
                plt.tight_layout()
                st.pyplot(fig)
                plt.close()
        except ImportError:
            st.warning("Install seaborn for heatmap: pip install seaborn")


# ---------------------------------------------------------------------------
# Page: Export
# ---------------------------------------------------------------------------
def page_export():
    st.markdown("""
    <div class="page-header">
        <h2>Export</h2>
        <p>Download results, BrainGlobe files, and QC reports</p>
    </div>""", unsafe_allow_html=True)

    if st.session_state.animal_df is None:
        st.warning("No results yet — run the pipeline first.")
        return

    df       = st.session_state.animal_df
    cohort   = st.session_state.cohort_df

    st.markdown("#### Data tables")

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("**Animal-level projections**")
        st.caption("Per-region, per-animal — in construct space with full compound keys")
        csv_animal = df.to_csv(index=False).encode()
        st.download_button(
            "⬇ Download animal_projections.csv",
            data=csv_animal,
            file_name="animal_projections.csv",
            mime="text/csv",
            use_container_width=True,
        )

    with col2:
        st.markdown("**Cohort summary**")
        st.caption("Group means, SEMs, z-scores, co-projection flags")
        csv_cohort = cohort.to_csv(index=False).encode()
        st.download_button(
            "⬇ Download cohort_summary.csv",
            data=csv_cohort,
            file_name="cohort_summary.csv",
            mime="text/csv",
            use_container_width=True,
        )

    st.markdown('<div class="section-divider"></div>', unsafe_allow_html=True)
    st.markdown("#### BrainGlobe heatmap exports")
    st.caption(
        "Region acronym → value dictionaries, ready for "
        "`brainglobe-heatmap` and `brainrender`"
    )

    avg = df.groupby("region_acronym").agg(
        trapcre_fraction=("trapcre_fraction","mean"),
        traptta_fraction=("traptta_fraction","mean"),
        ratio=("ratio","mean"),
    ).reset_index()

    col1, col2, col3 = st.columns(3)
    exports = [
        ("TRAP-Cre projection",  "trapcre_fraction",  "trapcre_projection.json",  col1),
        ("TRAP-tTA projection",  "traptta_fraction",  "traptta_projection.json",  col2),
        ("Cre / tTA ratio",      "ratio",             "ratio.json",               col3),
    ]
    for label, col, fname, container in exports:
        d    = avg.dropna(subset=[col]).set_index("region_acronym")[col].to_dict()
        jstr = json.dumps(d, indent=2).encode()
        with container:
            st.markdown(f"**{label}**")
            st.download_button(
                f"⬇ {fname}",
                data=jstr,
                file_name=fname,
                mime="application/json",
                use_container_width=True,
                key=f"dl_{fname}",
            )

    st.markdown('<div class="section-divider"></div>', unsafe_allow_html=True)
    st.markdown("#### QC")

    col1, col2 = st.columns(2)
    with col1:
        # Counterbalance report
        cb_report = df[["animal_id","sex","ch1_construct",
                        "trapcre_compound_key","traptta_compound_key"]].drop_duplicates()
        st.markdown("**Counterbalance report**")
        st.download_button(
            "⬇ counterbalance_report.csv",
            data=cb_report.to_csv(index=False).encode(),
            file_name="counterbalance_report.csv",
            mime="text/csv",
            use_container_width=True,
        )
    with col2:
        # Sparse region flags
        sparse = df[df["sparse_flag"]][
            ["animal_id","region_acronym","region_name","n_sections"]
        ].drop_duplicates().sort_values("n_sections")
        st.markdown("**Sparse region flags**")
        st.download_button(
            "⬇ sparse_regions.csv",
            data=sparse.to_csv(index=False).encode(),
            file_name="sparse_regions.csv",
            mime="text/csv",
            use_container_width=True,
        )

    st.markdown('<div class="section-divider"></div>', unsafe_allow_html=True)
    st.markdown("#### BrainGlobe usage example")
    st.code("""
from brainglobe_heatmap import Heatmap
import json

# Load TractQuant export
values = json.load(open("trapcre_projection.json"))

# Generate flatmap
h = Heatmap(
    values,
    position=5000,
    orientation="frontal",
    title="TRAP-Cre projection",
    cmap="Purples",
)
h.show()
    """, language="python")


# ---------------------------------------------------------------------------
# Router
# ---------------------------------------------------------------------------
PAGE_FUNCS = {
    "setup":   page_setup,
    "data":    page_data,
    "config":  page_config,
    "run":     page_run,
    "results": page_results,
    "viz":     page_viz,
    "export":  page_export,
}

PAGE_FUNCS[st.session_state.page]()
