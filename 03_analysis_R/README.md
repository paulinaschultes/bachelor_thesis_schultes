# 03_analysis_R — statistical analysis and figure scripts

R scripts used for the statistical analysis and the figures of the bachelor thesis
*"Stressor-specific whole-brain projection connectivity of IPACL neurons using
activity-dependent viral tracing in mice"* (University of Regensburg, 2026).

Each script carries a header block listing the thesis figures it produces.

## Figure ↔ script mapping

| Figure | Thesis section | Script |
|---|---|---|
| 1–9 | 2.2–2.9 | not script-generated (schematics, timelines, microscopy) |
| 10 | 3.1 | `01_validation.R` |
| 11 | 3.1 | `01_validation.R` |
| 12 | 3.2.1 | `02_wholebrain_heatmaps.R` |
| 13 | 3.3.2 | `03_condition_statistics.R` |
| 14 | 3.3.3 | `04_clustering.R` |
| 15 | 3.3.3 | `04_clustering.R` |
| 16 A / B | 3.4 | `06_top_regions.R` / `05_networks.R` |
| 17 | 3.4 | histological validation, no script |
| 18 | 3.4 | `09_volcano_plots.R` |
| 19 | 3.4 | manual verification, no script |
| 20 | 3.5 | `02_wholebrain_heatmaps.R` |
| 21 | 3.5 | `09_volcano_plots.R` |
| 22 | 3.5 | `04_clustering.R` |
| 23 | 3.5 | `04_clustering.R` |
| 24 A / B | 3.5 | `06_top_regions.R` / `05_networks.R` |
| 25 | 3.5 | `11_venn_diagrams.R` |
| 26 | 3.6 | `08_channel_analyses.R` |
| 27 | 3.7 | `12_retrograde.R` |
| 28 | 3.7 | `11_venn_diagrams.R` |
| 29 | 3.8 | `12_retrograde.R` |
| 30 | 3.9 | `12_retrograde.R` |

Scripts without a main figure:

| Script | Purpose |
|---|---|
| `07_sex_analyses.R` | region-wise Kruskal-Wallis with BH-FDR for section 3.5; source tables for Figures 21 and 24 |
| `10_power_and_tables.R` | achieved power and minimum detectable effect (Kruskal-Wallis based), formatted result tables |
| `10b_posthoc_power_ttest.R` | post-hoc power analysis reported in section 2.10 (two-sample t-test, d = 1.2, Bonferroni over 588 regions → alpha = 8.50e-05) |

**Note on figure numbers inside the code.** Comments within the scripts sometimes
refer to figure numbers from earlier drafts of the thesis. The table above and the
header block of each script reflect the final numbering.

## Requirements

R 4.5.1. Packages used across the scripts:

```
ComplexHeatmap, circlize, dabestr, dendextend, dplyr, ggdendro, ggplot2, ggraph,
ggrepel, grid, gt, igraph, magick, mclust, patchwork, pheatmap, pwr, readr,
readxl, reshape2, tidygraph, tidyr, writexl
```

`10b_posthoc_power_ttest.R` needs base R only.

## Data and paths

The scripts were written to run on the author's machine and read from and write to
absolute paths under `~/Desktop/Bachelor`. To re-run them, adjust these paths.

Only the two validation datasets (`validation_TRAP-Cre.csv`, `validation_TRAP-tTA.csv`
in the repository root) are included here. The whole-brain projection intensity tables,
the retrograde cell-count tables and the Allen CCF structure file are not part of this
repository and are available on request.

## AI use

The code in this repository was written with the assistance of Claude Code (Anthropic),
as declared in the thesis.
