"""
TractQuant — visualisation module
==================================
All publication-quality plots from pipeline output.

Requires: matplotlib, seaborn, pandas, numpy
Optional: brainglobe-heatmap (for flatmaps)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from pathlib import Path


# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------

CRE_COLOR  = "#534AB7"   # purple — always TRAP-Cre
TTA_COLOR  = "#1D9E75"   # green  — always TRAP-tTA
SEX_COLORS = {"Male": "#378ADD", "Female": "#D4537E"}

DIVISION_COLORS = {
    "Isocortex":  "#534AB7",
    "Striatum":   "#1D9E75",
    "Thalamus":   "#BA7517",
    "Midbrain":   "#D85A30",
    "Pallidum":   "#888780",
    "Brainstem":  "#3C3489",
    "Hypothalamus":"#639922",
    "Hindbrain":  "#D4537E",
}

def set_style():
    plt.rcParams.update({
        "font.family":      "sans-serif",
        "font.size":        10,
        "axes.spines.top":  False,
        "axes.spines.right":False,
        "axes.linewidth":   0.5,
        "xtick.major.width":0.5,
        "ytick.major.width":0.5,
        "figure.dpi":       150,
    })


# ---------------------------------------------------------------------------
# 1. Bar chart — top regions, Ch1 vs Ch2 side by side
# ---------------------------------------------------------------------------

def plot_top_regions_bar(
    animal_df: pd.DataFrame,
    n_top: int = 15,
    sort_by: str = "trapcre_fraction",
    output_path: Path = None,
    title: str = None,
):
    """
    Horizontal bar chart of top N regions, TRAP-Cre and TRAP-tTA side by side.
    Each animal or cohort mean can be passed.
    """
    set_style()

    df = animal_df.copy()
    df = df[~df["region_acronym"].isin(["root","void"])].copy()
    df = df.nlargest(n_top, sort_by).sort_values(sort_by, ascending=True)

    fig, ax = plt.subplots(figsize=(7, n_top * 0.45 + 1))

    y     = np.arange(len(df))
    h     = 0.35

    ax.barh(y + h/2, df["trapcre_fraction"] * 100, h,
            color=CRE_COLOR, alpha=0.85, label="TRAP-Cre")
    ax.barh(y - h/2, df["traptta_fraction"] * 100, h,
            color=TTA_COLOR, alpha=0.85, label="TRAP-tTA")

    ax.set_yticks(y)
    ax.set_yticklabels(df["region_acronym"], fontsize=9)
    ax.set_xlabel("Projection fraction (% of injection site)")
    ax.legend(frameon=False, fontsize=9)
    if title:
        ax.set_title(title, fontsize=11, fontweight="normal")

    plt.tight_layout()
    if output_path:
        plt.savefig(output_path, bbox_inches="tight")
    return fig


# ---------------------------------------------------------------------------
# 2. Scatter — TRAP-Cre vs TRAP-tTA per region
# ---------------------------------------------------------------------------

def plot_cre_vs_tta_scatter(
    animal_df: pd.DataFrame,
    min_fraction: float = 0.01,
    label_top_n: int = 12,
    output_path: Path = None,
    title: str = None,
):
    """
    Scatter plot: x = TRAP-Cre fraction, y = TRAP-tTA fraction.
    Each dot = one brain region, colored by major division.
    Diagonal = equal projection.
    Points above diagonal = tTA dominant.
    Points below = Cre dominant.
    """
    set_style()

    df = animal_df.copy()
    df = df[
        (df["trapcre_fraction"] > min_fraction) |
        (df["traptta_fraction"] > min_fraction)
    ].copy()

    fig, ax = plt.subplots(figsize=(6, 6))

    # Diagonal (equal projection line)
    lim = max(df["trapcre_fraction"].max(), df["traptta_fraction"].max()) * 1.1
    ax.plot([0, lim], [0, lim], color="#BBBBBB", linewidth=0.8,
            linestyle="--", zorder=0, label="Equal projection")

    # Points per division
    for div, grp in df.groupby("division"):
        color = DIVISION_COLORS.get(div, "#AAAAAA")
        ax.scatter(
            grp["trapcre_fraction"] * 100,
            grp["traptta_fraction"] * 100,
            color=color, alpha=0.75, s=40, label=div, zorder=2
        )

    # Label top regions
    top = df.nlargest(label_top_n, "trapcre_fraction")
    for _, row in top.iterrows():
        ax.annotate(
            row["region_acronym"],
            xy=(row["trapcre_fraction"]*100, row["traptta_fraction"]*100),
            xytext=(3, 3), textcoords="offset points",
            fontsize=7, color="#444444"
        )

    ax.set_xlabel("TRAP-Cre projection fraction (%)", fontsize=10)
    ax.set_ylabel("TRAP-tTA projection fraction (%)", fontsize=10)
    ax.set_xlim(0, lim * 100)
    ax.set_ylim(0, lim * 100)

    # Annotations
    ax.text(0.98, 0.02, "Cre dominant", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=8, color="#999999")
    ax.text(0.02, 0.98, "tTA dominant", transform=ax.transAxes,
            ha="left", va="top", fontsize=8, color="#999999")

    ax.legend(frameon=False, fontsize=8, loc="upper left",
              bbox_to_anchor=(1.01, 1))
    if title:
        ax.set_title(title, fontsize=11, fontweight="normal")

    plt.tight_layout()
    if output_path:
        plt.savefig(output_path, bbox_inches="tight")
    return fig


# ---------------------------------------------------------------------------
# 3. Heatmap — region × animal matrix
# ---------------------------------------------------------------------------

def plot_cohort_heatmap(
    animal_dfs: list[pd.DataFrame],
    construct: str = "trapcre",   # "trapcre" or "traptta"
    n_top_regions: int = 30,
    group_col: str = None,        # optional column to sort/annotate animals by
    output_path: Path = None,
):
    """
    Heatmap: rows = brain regions, columns = animals.
    Sorted by mean projection fraction.
    """
    set_style()

    col = f"{construct}_fraction"
    id_col = "animal_id"

    # Pivot: region × animal
    all_df = pd.concat(animal_dfs)
    pivot = all_df.pivot_table(
        index="region_acronym", columns=id_col, values=col, aggfunc="mean"
    )
    # Top N regions by mean
    top_regions = pivot.mean(axis=1).nlargest(n_top_regions).index
    pivot = pivot.loc[top_regions]

    fig, ax = plt.subplots(figsize=(max(6, len(pivot.columns)*0.5), n_top_regions*0.35+1))

    color = CRE_COLOR if construct == "trapcre" else TTA_COLOR
    cmap = sns.light_palette(color, as_cmap=True)

    sns.heatmap(
        pivot * 100,
        ax=ax,
        cmap=cmap,
        linewidths=0.3,
        linecolor="#EEEEEE",
        cbar_kws={"label": "Projection fraction (%)", "shrink": 0.6},
        fmt=".0f",
        annot=len(pivot.columns) <= 20,
        annot_kws={"size": 7},
    )
    construct_label = "TRAP-Cre" if construct == "trapcre" else "TRAP-tTA"
    ax.set_title(f"{construct_label} — region × animal heatmap", fontsize=11)
    ax.set_ylabel("")
    ax.set_xlabel("")
    ax.tick_params(axis="x", labelrotation=45, labelsize=8)
    ax.tick_params(axis="y", labelsize=8)

    plt.tight_layout()
    if output_path:
        plt.savefig(output_path, bbox_inches="tight")
    return fig


# ---------------------------------------------------------------------------
# 4. Z-score map across conditions
# ---------------------------------------------------------------------------

def plot_zscore_comparison(
    cohort_df: pd.DataFrame,
    group_a_key: str,           # compound key for group A e.g. "TRAP-Cre · Restraint · Acute · Day 1"
    group_b_key: str,           # compound key for group B
    n_regions: int = 20,
    output_path: Path = None,
):
    """
    Bar chart of per-region z-scores, highlighting regions
    that differ between two compound groups.
    """
    set_style()

    col_a = cohort_df[cohort_df["trapcre_compound_key"] == group_a_key].set_index("region_acronym")
    col_b = cohort_df[cohort_df["trapcre_compound_key"] == group_b_key].set_index("region_acronym")

    shared = col_a.index.intersection(col_b.index)
    diff = (col_a.loc[shared, "trapcre_mean"] - col_b.loc[shared, "trapcre_mean"]).sort_values()

    # Top + bottom N
    n = n_regions // 2
    display = pd.concat([diff.head(n), diff.tail(n)]).drop_duplicates()

    colors = [CRE_COLOR if v > 0 else TTA_COLOR for v in display.values]

    fig, ax = plt.subplots(figsize=(6, len(display) * 0.4 + 1))
    ax.barh(range(len(display)), display.values * 100, color=colors, alpha=0.85)
    ax.set_yticks(range(len(display)))
    ax.set_yticklabels(display.index, fontsize=8)
    ax.axvline(0, color="#AAAAAA", linewidth=0.8)
    ax.set_xlabel("Projection difference (% A − B)")
    ax.set_title(
        f"Differential projection\n{group_a_key}\nvs {group_b_key}",
        fontsize=9, fontweight="normal"
    )

    plt.tight_layout()
    if output_path:
        plt.savefig(output_path, bbox_inches="tight")
    return fig


# ---------------------------------------------------------------------------
# 5. BrainGlobe flatmap
# ---------------------------------------------------------------------------

def plot_brainglobe_flatmap(
    animal_df: pd.DataFrame,
    output_dir: Path,
    atlas_name: str = "allen_mouse_25um",
):
    """
    Generate BrainGlobe flatmaps for TRAP-Cre, TRAP-tTA, and ratio.
    Requires: pip install brainglobe-heatmap
    """
    try:
        from brainglobe_heatmap import Heatmap
    except ImportError:
        print("brainglobe-heatmap not installed.")
        print("Install with: pip install brainglobe-heatmap")
        return

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for col, label, color in [
        ("trapcre_fraction", "TRAP-Cre projection", CRE_COLOR),
        ("traptta_fraction", "TRAP-tTA projection", TTA_COLOR),
        ("trapcre_traptta_ratio", "Cre/tTA ratio",  "#D85A30"),
    ]:
        values = (
            animal_df.dropna(subset=[col])
            .set_index("region_acronym")[col]
            .to_dict()
        )

        try:
            h = Heatmap(
                values,
                position=5000,
                orientation="frontal",
                atlas_name=atlas_name,
                cmap=sns.light_palette(color, as_cmap=True),
                label_regions=False,
            )
            fig = h.show(show=False)
            safe_label = label.replace("/","_").replace(" ","_")
            fig.savefig(output_dir / f"{safe_label}_flatmap.png",
                        dpi=150, bbox_inches="tight")
            plt.close(fig)
            print(f"Saved flatmap: {label}")
        except Exception as e:
            print(f"Flatmap failed for {label}: {e}")


# ---------------------------------------------------------------------------
# Convenience: generate all figures for one animal
# ---------------------------------------------------------------------------

def generate_all_figures(
    animal_df: pd.DataFrame,
    output_dir: Path,
    animal_id: str = "",
):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    pfx = f"{animal_id}_" if animal_id else ""

    plot_top_regions_bar(
        animal_df, n_top=20,
        title=f"{animal_id} — top projection regions",
        output_path=output_dir / f"{pfx}top_regions.png"
    )
    plot_cre_vs_tta_scatter(
        animal_df,
        title=f"{animal_id} — TRAP-Cre vs TRAP-tTA",
        output_path=output_dir / f"{pfx}scatter.png"
    )
    plot_brainglobe_flatmap(animal_df, output_dir / f"{pfx}flatmaps")

    plt.close("all")
    print(f"Figures saved to {output_dir}/")
