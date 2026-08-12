"""
TractQuant — example usage
==========================
The Allen CCF atlas (~480 MB) is downloaded automatically on first run
and cached in ~/.tractquant/atlas/ — no manual setup required.
"""

from pathlib import Path
from pipeline import AnimalConfig, ChannelMeta, run_pipeline

# ---------------------------------------------------------------------------
# Define animals
# ---------------------------------------------------------------------------

# Scenario A: different stress duration per channel, same stressor
animal_F5 = AnimalConfig(
    animal_id="F5ex",
    sex="Male",
    ch1=ChannelMeta(
        construct="TRAP-Cre",
        stress_type="Restraint",
        stress_duration="Acute",
        timepoint="Day 1",
    ),
    ch2_stress_type="Restraint",
    ch2_stress_duration="Chronic",
    ch2_timepoint="Day 14",
    visual_align_json=Path("data/F5ex/F5visu.json"),
    ch1_image_dir=Path("data/F5ex/ch1"),
    ch2_image_dir=Path("data/F5ex/ch2"),
)

# Counterbalanced: Ch1 is now TRAP-tTA
animal_F6 = AnimalConfig(
    animal_id="F6ex",
    sex="Female",
    ch1=ChannelMeta(
        construct="TRAP-tTA",
        stress_type="Restraint",
        stress_duration="Acute",
        timepoint="Day 1",
    ),
    ch2_stress_type="Restraint",
    ch2_stress_duration="Chronic",
    ch2_timepoint="Day 14",
    visual_align_json=Path("data/F6ex/F6visu.json"),
    ch1_image_dir=Path("data/F6ex/ch1"),
    ch2_image_dir=Path("data/F6ex/ch2"),
)

animals = [animal_F5, animal_F6]

# ---------------------------------------------------------------------------
# Verify compound keys before running
# ---------------------------------------------------------------------------
print("=== Animal compound keys ===")
for an in animals:
    print(f"\n{an.animal_id} ({an.sex})")
    print(f"  Ch1: {an.ch1.construct:<10} {an.ch1.compound_key}")
    print(f"  Ch2: {an.ch2.construct:<10} {an.ch2.compound_key}")
print()

# ---------------------------------------------------------------------------
# Run pipeline
# ---------------------------------------------------------------------------
# Atlas downloads automatically on first run (~480 MB, one-time only).
# Cached in ~/.tractquant/atlas/ — subsequent runs load from cache instantly.
#
# Override options:
#   atlas_path="path/to/annotation_25.nrrd"   use local file
#   atlas_cache_dir=Path("/your/cache/dir")    change cache location

cohort = run_pipeline(
    animals=animals,
    injection_regions=["MOs", "MOp"],
    output_dir="results/",
    resize_factor=0.25,
    sparse_flag_threshold=3,
    save_qc=True,
)

print("\n=== Top 10 TRAP-Cre projection regions ===")
top = (
    cohort.groupby("region_acronym")["trapcre_mean"]
    .mean().sort_values(ascending=False).head(10)
)
print(top.to_string())
