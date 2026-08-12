# =============================================================================
# Whole-brain and correlation heatmaps
#
# Consolidated from the following original scripts:
#   - prepare_heatmap_data.R
#   - plot_heatmap_ALL_DT.R
#   - plot_heatmap_ALL_DT_v2.R
#   - plot_heatmap_ALL_DT_by_sex.R
#   - plot_heatmap_ALL_DT_by_sex_no_chronicFS.R
#   - results_summary_heatmap.R
#   - correlation_heatmap_sex.R
#   - correlation_heatmap_split_sex.R
#   - nullmodel_corrected_correlation_heatmaps.R
# =============================================================================


# ---------------------------------------------------------------------------
# prepare_heatmap_data
# ---------------------------------------------------------------------------
# Shared data prep for the ALL_DT heatmaps (compact overview + full labeled version).
# Builds: mat, condition_per_col, major_division_per_row,
#         row_anno, top_anno, col_fun

library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)

## --- Load data ---
df <- read.csv("~/Desktop/ALL_DT.csv", stringsAsFactors = FALSE)
structures <- read.csv("~/Desktop/tractquant/allen_mouse_25um_v1.2/structures.csv", stringsAsFactors = FALSE)

## --- Fixed condition order (numbered 1-5; stressor identity given in the
## thesis legend, not on the plot itself) ---
condition_order <- c(
  "forced swim · Acute",
  "forced swim · Chronic",
  "restraint · Acute",
  "social defeat · Acute",
  "tail suspension · Acute"
)

## --- Melt Cre/tTA into a single (region x mouse x channel) table with a
## per-row condition. Each animal_group's Cre/tTA -> condition mapping is now
## fixed (verified against the raw per-mouse condition tags), so pooling
## across animal_group and channel by condition is well-defined. ---
strip_condition <- function(x) sub("^TRAP-(Cre|tTA) · ", "", x)

long <- bind_rows(
  df %>% transmute(region_id, region_acronym, division,
                    animal_group, animal_id_original, sex, channel = "Cre",
                    intensity = trapcre_intensity,
                    condition = strip_condition(trapcre_condition)),
  df %>% transmute(region_id, region_acronym, division,
                    animal_group, animal_id_original, sex, channel = "tTA",
                    intensity = traptta_intensity,
                    condition = strip_condition(traptta_condition))
) %>%
  mutate(condition = factor(condition, levels = condition_order))

## --- Aggregate: median intensity per region x condition, pooled across all
## animal groups/channels sharing that condition. ---
agg <- long %>%
  group_by(region_id, region_acronym, division, condition) %>%
  summarise(intensity = median(intensity, na.rm = TRUE), .groups = "drop")

## --- Pivot to wide matrix: rows = regions, columns = condition ---
wide <- agg %>%
  select(region_acronym, condition, intensity) %>%
  pivot_wider(names_from = condition, values_from = intensity)

mat <- as.matrix(wide[, -1])
rownames(mat) <- wide$region_acronym
mat <- mat[, condition_order]   # enforce exact column order

## --- Column grouping/labels: numbered 1-5, in condition_order ---
colnames(mat) <- as.character(seq_along(condition_order))
condition_per_col <- factor(colnames(mat), levels = colnames(mat))

## --- Major anatomical division per region (Allen CCF ontology) ---
## Walk each region's structure_id_path and keep the deepest ID that matches
## one of the standard top-level brain divisions.
major_ids <- c(
  "315"  = "Isocortex",
  "698"  = "Olfactory areas",
  "1089" = "Hippocampal formation",
  "703"  = "Cortical subplate",
  "477"  = "Striatum",
  "803"  = "Pallidum",
  "549"  = "Thalamus",
  "1097" = "Hypothalamus",
  "313"  = "Midbrain",
  "771"  = "Pons",
  "354"  = "Medulla",
  "512"  = "Cerebellum",
  "1009" = "fiber tracts",
  "73"   = "ventricular systems"
)
major_order <- unname(major_ids)  # rostral -> caudal, then fiber tracts, then ventricles

get_major_division <- function(path) {
  ids <- strsplit(gsub("^/|/$", "", path), "/")[[1]]
  hit <- ids[ids %in% names(major_ids)]
  if (length(hit) == 0) return(NA_character_)
  major_ids[[hit[length(hit)]]]
}

region_lookup <- df %>%
  distinct(region_id, region_acronym) %>%
  left_join(structures %>% select(id, structure_id_path), by = c("region_id" = "id")) %>%
  rowwise() %>%
  mutate(major_division = get_major_division(structure_id_path)) %>%
  ungroup()

major_division_per_row <- region_lookup$major_division[match(rownames(mat), region_lookup$region_acronym)]
major_division_per_row <- factor(major_division_per_row, levels = major_order)

## --- Exclude fiber tracts (white matter; not a projection-bearing gray
## matter region, and reviewer comment #7 requested it be excluded from the
## results throughout) from every downstream analysis. Filtered here, at
## the shared data-prep source, so every script that sources this file
## inherits the exclusion automatically without needing its own filter. ---
fiber_tract_regions <- region_lookup$region_acronym[region_lookup$major_division == "fiber tracts"]
fiber_tract_regions <- fiber_tract_regions[!is.na(fiber_tract_regions)]

keep_rows <- !(rownames(mat) %in% fiber_tract_regions)
mat <- mat[keep_rows, , drop = FALSE]
major_division_per_row <- major_division_per_row[keep_rows]

major_order <- setdiff(major_order, "fiber tracts")
major_division_per_row <- factor(as.character(major_division_per_row), levels = major_order)

long <- long %>% filter(!region_acronym %in% fiber_tract_regions)

## --- Row annotation: major anatomical division (muted, print-friendly) ---
major_cols <- c(
  "Isocortex"             = "#5B7C99",
  "Olfactory areas"       = "#8C6D46",
  "Hippocampal formation" = "#6B8F71",
  "Cortical subplate"     = "#9C6B85",
  "Striatum"              = "#4F8F8B",
  "Pallidum"              = "#7A6B9C",
  "Thalamus"              = "#B08B5A",
  "Hypothalamus"          = "#A85C5C",
  "Midbrain"              = "#6B6B99",
  "Pons"                  = "#8FA05E",
  "Medulla"               = "#5E8FA0",
  "Cerebellum"            = "#9C7A5B",
  "ventricular systems"   = "#C4B7A6"
)

row_anno <- rowAnnotation(
  division = major_division_per_row,
  col = list(division = major_cols),
  annotation_name_gp = gpar(fontsize = 20),
  annotation_legend_param = list(division = list(
    ncol = 1, title = "Major division",
    title_gp = gpar(fontsize = 20, fontface = "bold"),
    labels_gp = gpar(fontsize = 18),
    grid_height = unit(4.5, "mm"), grid_width = unit(4.5, "mm")
  ))
)

## --- Column annotation: condition (numbered 1-5; stressor names given in
## the thesis legend, not here) ---
## Conditions 1-2 (forced swim, acute/chronic) share a hue family since they
## are the same stressor at different durations; 3-5 get distinct hues.
condition_cols <- setNames(
  c("#4C72B0", "#2E4459", "#C44E52", "#55A868", "#8172B2"),
  colnames(mat)
)

top_anno <- HeatmapAnnotation(
  condition = condition_per_col,
  col = list(condition = condition_cols),
  annotation_name_gp = gpar(fontsize = 20),
  annotation_legend_param = list(
    condition = list(
      title = "Condition",
      title_gp = gpar(fontsize = 20, fontface = "bold"),
      labels_gp = gpar(fontsize = 18),
      grid_height = unit(4.5, "mm"), grid_width = unit(4.5, "mm")
    )
  ),
  show_legend = TRUE,
  simple_anno_size = unit(4, "mm")
)

## --- Pink-to-violet color scale, matching the viral construct diagram (Fig. 1) ---
## c-Fos/pA = crimson-pink, Cre-ERT2/EYFP/rtTA = violet. Breakpoints at quantiles
## (data is heavily right-skewed) so small differences among the bulk of low
## values stay visible; luminance still decreases monotonically for print.
qs <- quantile(mat, probs = c(0, 0.10, 0.25, 0.50, 0.75, 0.90, 0.99), na.rm = TRUE)
col_fun <- colorRamp2(
  as.numeric(qs),
  c("#FDF0F5", "#FBC9DE", "#F48FB1", "#E05780", "#C2185B", "#8E24AA", "#4A148C")
)


# ---------------------------------------------------------------------------
# plot_heatmap_ALL_DT
# ---------------------------------------------------------------------------
# Compact overview heatmap: fits on a single thesis page.
# Rows grouped into a sidebar by major anatomical division (Allen CCF ontology),
# clustered by similarity within each division. No per-region row labels --
# 670 regions can't stay legible in this footprint (see
# plot_heatmap_ALL_DT_full_labeled.R for the version with all region names).

library(ComplexHeatmap)
library(circlize)

setwd("~/Desktop/Bachelor")
source("prepare_heatmap_data.R")

## --- Plot heatmap ---
ht <- Heatmap(
  mat,
  name = "intensity",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_row_slices = FALSE,     # keep divisions in anatomical order, not reordered by similarity
  row_split = major_division_per_row,
  row_gap = unit(0.6, "mm"),
  row_title = NULL,               # grouping is conveyed by the sidebar, not per-block text
  row_dend_width = unit(6, "cm"),   # extra-wide dendrogram so even small divisions' branches are resolvable
  cluster_columns = FALSE,        # fixed 1-5 condition order, not reordered by similarity
  show_row_names = FALSE,         # 670 regions -> no room for labels on a single page
  show_column_names = TRUE,       # just "1".."5" -- stressor identity is in the thesis legend
  column_names_gp = gpar(fontsize = 6.5),
  column_names_rot = 0,
  top_annotation = top_anno,
  right_annotation = row_anno,
  row_title_gp = gpar(fontsize = 7),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 7, fontface = "bold"),
    labels_gp = gpar(fontsize = 6),
    grid_height = unit(2.8, "mm"), grid_width = unit(2.8, "mm")
  )
)

draw(ht, merge_legend = TRUE)

## --- Save to file: real DIN A4 pages (portrait) ---
## Page 1 = the heatmap. Page 2+ = a legend listing, per major division, the
## region acronyms in the exact top-to-bottom order they appear in that
## division's block (i.e. the row dendrogram's leaf order) -- since the
## heatmap itself has no room for 670 individual row labels.
a4_width  <- 8.27   # DIN A4 portrait, inches
a4_height <- 11.69
pdf("~/Desktop/Bachelor/ALL_DT_heatmap.pdf", width = a4_width, height = a4_height)
ht_drawn <- draw(ht, merge_legend = TRUE)

ro <- row_order(ht_drawn)   # list of row indices per division, in display order
division_legend <- lapply(major_order, function(d) rownames(mat)[ro[[d]]])
names(division_legend) <- major_order

line_height   <- 0.021
top_margin    <- 0.97
bottom_margin <- 0.03
wrap_width    <- 78

new_legend_page <- function() {
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
}

new_legend_page()
text(0.02, top_margin, "Region order within each heatmap division (top to bottom, as clustered)",
     adj = c(0, 1), font = 2, cex = 0.85)
y <- top_margin - 0.045

for (d in major_order) {
  regions <- division_legend[[d]]
  if (length(regions) == 0) next
  header  <- sprintf("%s (n = %d)", d, length(regions))
  body    <- paste(seq_along(regions), regions, sep = ". ", collapse = "   ")
  wrapped <- strwrap(body, width = wrap_width)
  needed  <- (2 + length(wrapped)) * line_height

  if (y - needed < bottom_margin) {
    new_legend_page()
    y <- top_margin
  }

  rect(0.02, y - 0.010, 0.032, y + 0.006, col = major_cols[[d]], border = NA)
  text(0.045, y, header, adj = c(0, 1), font = 2, cex = 0.75)
  y <- y - line_height * 1.4

  for (ln in wrapped) {
    text(0.02, y, ln, adj = c(0, 1), family = "mono", cex = 0.55)
    y <- y - line_height
  }
  y <- y - line_height * 0.7
}

dev.off()


# ---------------------------------------------------------------------------
# plot_heatmap_ALL_DT_v2
# ---------------------------------------------------------------------------
# v2: same single-heatmap overview as plot_heatmap_ALL_DT.R, fixing
# reviewer comment "Need to unify the cluster trees (same height)".
#
# ComplexHeatmap's cluster_rows only accepts TRUE/FALSE / a single hclust /
# a single dendrogram / OR A FUNCTION. When row_split is active and
# cluster_rows is a function, ComplexHeatmap calls that function
# SEPARATELY on each division's own sub-matrix and uses its return value
# (an hclust or dendrogram) for that slice -- this is the hook used here to
# rescale each division's own merge heights to the same 0-1 range before
# they're drawn, so every division's dendrogram fills the same width
# regardless of its absolute distance scale (previously: ventricular
# systems max height 0.26 vs. Striatum max height 12.17, a 48x difference,
# made small divisions' trees an unreadable sliver).
#
# Everything else is unchanged from plot_heatmap_ALL_DT.R: still one single
# Heatmap object, same row_split, same annotations, same layout.

library(ComplexHeatmap)
library(circlize)

setwd("~/Desktop/Bachelor")
source("prepare_heatmap_data.R")

## Custom row-clustering function: called once per row_split slice with
## that slice's own sub-matrix; rescales heights to [0,1] so every
## division's dendrogram is drawn at the same visual scale.
cluster_rows_unified_height <- function(m) {
  if (nrow(m) < 2) return(NULL)
  hc <- hclust(dist(m), method = "complete")
  if (max(hc$height) > 0) hc$height <- hc$height / max(hc$height)
  as.dendrogram(hc)
}

FONT_PT <- 12

## Re-built locally (rather than editing the shared prepare_heatmap_data.R,
## which plot_heatmap_ALL_DT.R and the by-sex variants also source) at the
## requested 12pt and with the condition legend switched off -- the color
## strip stays, it's just the "Condition" legend key that's redundant with
## the numbered column labels below the heatmap.
top_anno <- HeatmapAnnotation(
  condition = condition_per_col,
  col = list(condition = condition_cols),
  annotation_name_gp = gpar(fontsize = FONT_PT),
  show_legend = FALSE,
  simple_anno_size = unit(4, "mm")
)
row_anno <- rowAnnotation(
  division = major_division_per_row,
  col = list(division = major_cols),
  annotation_name_gp = gpar(fontsize = FONT_PT),
  annotation_legend_param = list(division = list(
    ncol = 1, title = "Major division",
    title_gp = gpar(fontsize = FONT_PT, fontface = "bold"),
    labels_gp = gpar(fontsize = FONT_PT),
    grid_height = unit(4.5, "mm"), grid_width = unit(4.5, "mm")
  ))
)

## --- Plot heatmap ---
ht <- Heatmap(
  mat,
  name = "intensity",
  col = col_fun,
  cluster_rows = cluster_rows_unified_height,  # per-slice, height-normalized
  cluster_row_slices = FALSE,     # keep divisions in anatomical order, not reordered by similarity
  row_split = major_division_per_row,
  row_gap = unit(0.6, "mm"),
  row_title = NULL,               # grouping is conveyed by the sidebar, not per-block text
  row_dend_width = unit(4.8, "cm"), # 20% shorter than before (was 6cm) -- more room for the heatmap
  cluster_columns = FALSE,        # fixed 1-5 condition order, not reordered by similarity
  show_row_names = FALSE,         # 670 regions -> no room for labels on a single page
  show_column_names = TRUE,       # just "1".."5" -- stressor identity is in the thesis legend
  column_names_gp = gpar(fontsize = FONT_PT),
  column_names_rot = 0,
  column_title_gp = gpar(fontsize = FONT_PT, fontface = "bold"),
  top_annotation = top_anno,
  right_annotation = row_anno,
  row_title_gp = gpar(fontsize = FONT_PT),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = FONT_PT, fontface = "bold"),
    labels_gp = gpar(fontsize = FONT_PT),
    grid_height = unit(4.5, "mm"), grid_width = unit(4.5, "mm")
  )
)

## --- Save to file: real DIN A4 pages (portrait) ---
a4_width  <- 8.27   # DIN A4 portrait, inches
a4_height <- 11.69
pdf("~/Desktop/Bachelor/ALL_DT_heatmap_v2_unified_dendrograms.pdf", width = a4_width, height = a4_height)
ht_drawn <- draw(ht, merge_legend = TRUE)

ro <- row_order(ht_drawn)   # list of row indices per division, in display order
division_legend <- lapply(major_order, function(d) rownames(mat)[ro[[d]]])
names(division_legend) <- major_order

line_height   <- 0.021
top_margin    <- 0.97
bottom_margin <- 0.03
wrap_width    <- 78

new_legend_page <- function() {
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
}

new_legend_page()
text(0.02, top_margin, "Region order within each heatmap division (top to bottom, as clustered)",
     adj = c(0, 1), font = 2, cex = 0.85)
y <- top_margin - 0.045

for (d in major_order) {
  regions <- division_legend[[d]]
  if (length(regions) == 0) next
  header  <- sprintf("%s (n = %d)", d, length(regions))
  body    <- paste(seq_along(regions), regions, sep = ". ", collapse = "   ")
  wrapped <- strwrap(body, width = wrap_width)
  needed  <- (2 + length(wrapped)) * line_height

  if (y - needed < bottom_margin) {
    new_legend_page()
    y <- top_margin
  }

  rect(0.02, y - 0.010, 0.032, y + 0.006, col = major_cols[[d]], border = NA)
  text(0.045, y, header, adj = c(0, 1), font = 2, cex = 0.75)
  y <- y - line_height * 1.4

  for (ln in wrapped) {
    text(0.02, y, ln, adj = c(0, 1), family = "mono", cex = 0.55)
    y <- y - line_height
  }
  y <- y - line_height * 0.7
}

dev.off()
cat("Saved: ALL_DT_heatmap_v2_unified_dendrograms.pdf\n")

png("~/Desktop/Bachelor/ALL_DT_heatmap_v2_unified_dendrograms.png", width = a4_width, height = a4_height, units = "in", res = 200)
draw(ht, merge_legend = TRUE)
dev.off()
cat("Saved: ALL_DT_heatmap_v2_unified_dendrograms.png\n")


# ---------------------------------------------------------------------------
# plot_heatmap_ALL_DT_by_sex
# ---------------------------------------------------------------------------
# Overview heatmap split by BOTH condition and sex: 670 regions x 8 columns
# (4 conditions x 2 sexes). Extends prepare_heatmap_data.R's `mat`
# (condition-only, sexes pooled) to additionally reveal sex as a source of
# variation within each condition.
#
# "social defeat - Acute" is EXCLUDED entirely from this heatmap: it is a
# single-sex condition (n=6 Male, 0 Female), so a sex-split comparison is not
# meaningful for it (there would be no Female column to compare against, only
# an empty/grey one). This mirrors how "forced swim - Chronic" used to be
# excluded from this figure while it was (mistakenly) all-female; now that
# its sex labels are corrected (6 Female, 3 Male) it is included as a normal
# mixed-sex condition, and social defeat -- now the single-sex condition --
# takes its place as the excluded one.

library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)

setwd("~/Desktop/Bachelor")
source("prepare_heatmap_data.R")   # long (region x mouse x channel x condition x sex), major_division_per_row, condition_order

## Custom row-clustering function (same fix as plot_heatmap_ALL_DT_v2.R):
## called once per row_split slice with that slice's own sub-matrix,
## rescales merge heights to [0,1] so every division's dendrogram is drawn
## at the same visual scale regardless of its absolute distance range.
cluster_rows_unified_height <- function(m) {
  if (nrow(m) < 2) return(NULL)
  hc <- hclust(dist(m), method = "complete")
  if (max(hc$height) > 0) hc$height <- hc$height / max(hc$height)
  as.dendrogram(hc)
}

sex_order <- c("Female", "Male")

## Social defeat is single-sex (Male only) post-correction -- excluded
## entirely from this sex-split heatmap, see header note above.
condition_order_sex_heatmap <- setdiff(condition_order, "social defeat · Acute")
long_sex_heatmap <- long %>% filter(condition %in% condition_order_sex_heatmap)

## --- Aggregate: median intensity per region x condition x sex, pooled
## across all animals/channels sharing that condition+sex combination ---
agg_sex <- long_sex_heatmap %>%
  group_by(region_id, region_acronym, division, condition, sex) %>%
  summarise(intensity = median(intensity, na.rm = TRUE), .groups = "drop") %>%
  mutate(col_id = paste(condition, sex, sep = "_"))

col_levels <- unlist(lapply(condition_order_sex_heatmap, function(c) paste(c, sex_order, sep = "_")))

wide_sex <- agg_sex %>%
  select(region_acronym, col_id, intensity) %>%
  pivot_wider(names_from = col_id, values_from = intensity)

## With social defeat excluded, every remaining condition x sex combination
## has at least one animal, so no explicit NA-fill column should be needed;
## kept as a safety net in case a future data change reintroduces a gap.
missing_cols <- setdiff(col_levels, names(wide_sex))
for (mc in missing_cols) wide_sex[[mc]] <- NA_real_

mat_sex <- as.matrix(wide_sex[, col_levels])
rownames(mat_sex) <- wide_sex$region_acronym
mat_sex <- mat_sex[rownames(mat), ]   # same row order/set as the main heatmap (includes FS)

cat("Condition x sex combinations with no animals (shown as empty/grey columns):\n")
print(missing_cols)

## --- Column annotation: condition (reusing condition_cols) + sex ---
condition_per_col_sex <- factor(rep(seq_along(condition_order_sex_heatmap), each = length(sex_order)),
                                 levels = seq_along(condition_order_sex_heatmap))
sex_per_col <- factor(rep(sex_order, times = length(condition_order_sex_heatmap)), levels = sex_order)

condition_cols_sex <- setNames(condition_cols[condition_order_sex_heatmap],
                                as.character(seq_along(condition_order_sex_heatmap)))
sex_cols <- c("Female" = "#D46A8C", "Male" = "#4C72B0")

FONT_PT <- 12

top_anno_sex <- HeatmapAnnotation(
  condition = condition_per_col_sex,
  sex = sex_per_col,
  col = list(condition = condition_cols_sex, sex = sex_cols),
  annotation_name_gp = gpar(fontsize = FONT_PT),
  annotation_legend_param = list(
    condition = list(title = "Condition", title_gp = gpar(fontsize = FONT_PT, fontface = "bold"),
                      labels_gp = gpar(fontsize = FONT_PT), grid_height = unit(4.5, "mm"), grid_width = unit(4.5, "mm")),
    sex = list(title = "Sex", title_gp = gpar(fontsize = FONT_PT, fontface = "bold"),
               labels_gp = gpar(fontsize = FONT_PT), grid_height = unit(4.5, "mm"), grid_width = unit(4.5, "mm"))
  ),
  show_legend = TRUE,
  simple_anno_size = unit(4, "mm")
)

## --- Plot ---
ht_sex <- Heatmap(
  mat_sex, name = "intensity", col = col_fun,
  cluster_rows = cluster_rows_unified_height, cluster_row_slices = FALSE,
  row_split = major_division_per_row, row_gap = unit(0.6, "mm"), row_title = NULL,
  row_dend_width = unit(6, "cm"),
  cluster_columns = FALSE,
  column_split = condition_per_col_sex, column_gap = unit(1.2, "mm"),
  column_title = as.character(seq_along(condition_order_sex_heatmap)),
  column_labels = rep(c("F", "M"), length(condition_order_sex_heatmap)),
  show_row_names = FALSE, show_column_names = TRUE,
  column_names_gp = gpar(fontsize = FONT_PT), column_names_rot = 0,
  column_title_gp = gpar(fontsize = FONT_PT, fontface = "bold"),
  top_annotation = top_anno_sex,
  right_annotation = row_anno,
  row_title_gp = gpar(fontsize = FONT_PT),
  na_col = "grey85",
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = FONT_PT, fontface = "bold"),
    labels_gp = gpar(fontsize = FONT_PT),
    grid_height = unit(4.5, "mm"), grid_width = unit(4.5, "mm")
  )
)

## Condition-number legend kept in the console log only (not drawn on the
## figure -- provided separately by hand elsewhere).
condition_number_legend <- paste(
  sprintf("%d = %s", seq_along(condition_order_sex_heatmap), condition_order_sex_heatmap),
  collapse = "     "
)

pdf("ALL_DT_heatmap_by_sex.pdf", width = 8.27, height = 11.69)
draw(ht_sex, merge_legend = TRUE)
dev.off()

png("ALL_DT_heatmap_by_sex.png", width = 8.27, height = 11.69, units = "in", res = 300)
draw(ht_sex, merge_legend = TRUE)
dev.off()

cat("\nSaved: ALL_DT_heatmap_by_sex.pdf, ALL_DT_heatmap_by_sex.png\n")
cat("Column order: 1F,1M, 2F,2M, 3F,3M, 4F,4M\n")
cat(sprintf("(%s)\n", condition_number_legend))
cat("social defeat - Acute excluded entirely (single-sex, Male only, n=6)\n")


# ---------------------------------------------------------------------------
# plot_heatmap_ALL_DT_by_sex_no_chronicFS
# ---------------------------------------------------------------------------
# Overview heatmap split by condition and sex, like plot_heatmap_ALL_DT_by_sex.R,
# but with "forced swim · Chronic" dropped entirely (all-Female, 0 Male --
# see sex_kruskal_by_condition.R), so every remaining column has real data
# and no grey/NA columns appear. 670 regions x 8 columns (4 conditions x 2
# sexes).

library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)

setwd("~/Desktop/Bachelor")
source("prepare_heatmap_data.R")   # long, major_division_per_row, condition_cols, col_fun

## Custom row-clustering function (same fix as plot_heatmap_ALL_DT_v2.R):
## called once per row_split slice with that slice's own sub-matrix,
## rescales merge heights to [0,1] so every division's dendrogram is drawn
## at the same visual scale regardless of its absolute distance range.
cluster_rows_unified_height <- function(m) {
  if (nrow(m) < 2) return(NULL)
  hc <- hclust(dist(m), method = "complete")
  if (max(hc$height) > 0) hc$height <- hc$height / max(hc$height)
  as.dendrogram(hc)
}

sex_order <- c("Female", "Male")
condition_order_no_chronicFS <- setdiff(condition_order, "social defeat · Acute")  # social defeat is now single-sex (Male only); excluded here so only conditions comparable across both sexes remain. Chronic FS is included since its sex labels were corrected.

## --- Aggregate: median intensity per region x condition x sex ---
agg_sex <- long %>%
  filter(condition %in% condition_order_no_chronicFS) %>%
  group_by(region_id, region_acronym, division, condition, sex) %>%
  summarise(intensity = median(intensity, na.rm = TRUE), .groups = "drop") %>%
  mutate(col_id = paste(condition, sex, sep = "_"))

col_levels <- unlist(lapply(condition_order_no_chronicFS, function(c) paste(c, sex_order, sep = "_")))

wide_sex <- agg_sex %>%
  select(region_acronym, col_id, intensity) %>%
  pivot_wider(names_from = col_id, values_from = intensity)

missing_cols <- setdiff(col_levels, names(wide_sex))
cat("Condition x sex combinations with no animals (shown as empty/grey columns):\n")
print(missing_cols)
# e.g. "social defeat · Acute_Female" is expected to be empty after the DT4
# sex correction (social defeat is now all-Male); fill with NA instead of
# asserting completeness, consistent with plot_heatmap_ALL_DT_by_sex.R.
for (mc in missing_cols) wide_sex[[mc]] <- NA_real_

mat_sex <- as.matrix(wide_sex[, col_levels])
rownames(mat_sex) <- wide_sex$region_acronym
mat_sex <- mat_sex[rownames(mat), ]   # same row order/set as the main heatmap (includes FS)

## --- Column annotation: condition (reusing condition_cols for the 4 kept
## conditions) + sex ---
condition_per_col_sex <- factor(rep(seq_along(condition_order_no_chronicFS), each = length(sex_order)),
                                 levels = seq_along(condition_order_no_chronicFS))
sex_per_col <- factor(rep(sex_order, times = length(condition_order_no_chronicFS)), levels = sex_order)

## condition_cols is keyed "1".."5" (original 5-condition order); re-key to
## "1".."4" for the kept conditions, preserving each condition's original color.
kept_idx <- match(condition_order_no_chronicFS, condition_order)
condition_cols_sex <- setNames(condition_cols[as.character(kept_idx)], as.character(seq_along(condition_order_no_chronicFS)))
sex_cols <- c("Female" = "#D46A8C", "Male" = "#4C72B0")

top_anno_sex <- HeatmapAnnotation(
  condition = condition_per_col_sex,
  sex = sex_per_col,
  col = list(condition = condition_cols_sex, sex = sex_cols),
  annotation_name_gp = gpar(fontsize = 7),
  annotation_legend_param = list(
    condition = list(title = "Condition", title_gp = gpar(fontsize = 7, fontface = "bold"),
                      labels_gp = gpar(fontsize = 6), grid_height = unit(2.8, "mm"), grid_width = unit(2.8, "mm")),
    sex = list(title = "Sex", title_gp = gpar(fontsize = 7, fontface = "bold"),
               labels_gp = gpar(fontsize = 6), grid_height = unit(2.8, "mm"), grid_width = unit(2.8, "mm"))
  ),
  show_legend = TRUE,
  simple_anno_size = unit(2.5, "mm")
)

## --- Plot ---
ht_sex <- Heatmap(
  mat_sex, name = "intensity", col = col_fun,
  cluster_rows = cluster_rows_unified_height, cluster_row_slices = FALSE,
  row_split = major_division_per_row, row_gap = unit(0.6, "mm"), row_title = NULL,
  row_dend_width = unit(6, "cm"),
  cluster_columns = FALSE,
  column_split = condition_per_col_sex, column_gap = unit(1.2, "mm"),
  column_title = as.character(seq_along(condition_order_no_chronicFS)),
  column_labels = rep(c("F", "M"), length(condition_order_no_chronicFS)),
  show_row_names = FALSE, show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 6.5), column_names_rot = 0,
  top_annotation = top_anno_sex,
  right_annotation = row_anno,
  row_title_gp = gpar(fontsize = 7),
  na_col = "grey85",
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 7, fontface = "bold"),
    labels_gp = gpar(fontsize = 6),
    grid_height = unit(2.8, "mm"), grid_width = unit(2.8, "mm")
  )
)

pdf("ALL_DT_heatmap_by_sex_no_chronicFS.pdf", width = 7.5, height = 11.69)
draw(ht_sex, merge_legend = TRUE)
dev.off()

png("ALL_DT_heatmap_by_sex_no_chronicFS.png", width = 7.5, height = 11.69, units = "in", res = 300)
draw(ht_sex, merge_legend = TRUE)
dev.off()

cat("\nSaved: ALL_DT_heatmap_by_sex_no_chronicFS.pdf / .png\n")
cat("Column order: 1F,1M, 2F,2M, 3F,3M, 4F,4M\n")
cat("(1=forced swim acute, 2=restraint, 3=social defeat, 4=tail suspension; chronic forced swim excluded, 0 Male animals)\n")


# ---------------------------------------------------------------------------
# results_summary_heatmap
# ---------------------------------------------------------------------------
# Whole-results summary heatmap: one compact figure collecting the two
# headline division-level omnibus tests (condition effect, sex effect)
# that are otherwise scattered across separate tables/figures in the
# results section. Rows = the 13 major anatomical divisions (fiber tracts
# excluded, consistent with every other figure); columns = the two test
# families; color = -log10(BH-adjusted p), so darker = stronger effect;
# asterisks mark BH-FDR significance. Meant as a single at-a-glance
# "where are the significant effects" figure for the end of the results
# section, complementing (not replacing) the detailed per-test figures.

library(dplyr)
library(ggplot2)

setwd("~/Desktop/Bachelor")

cond_df <- read.csv("condition_kruskal_by_division.csv", stringsAsFactors = FALSE) %>%
  transmute(major_division, test = "Condition\n(5 groups)", p_adj_BH)
sex_df <- read.csv("sex_kruskal_by_division.csv", stringsAsFactors = FALSE) %>%
  transmute(major_division, test = "Sex\n(F vs M)", p_adj_BH)

summary_df <- bind_rows(cond_df, sex_df) %>%
  mutate(
    neg_log10_p = -log10(p_adj_BH),
    sig_label = case_when(
      p_adj_BH < 0.001 ~ "***",
      p_adj_BH < 0.01  ~ "**",
      p_adj_BH < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

## Row order: divisions ranked by their strongest (smallest) p-value across
## either test, most significant at the top.
division_order <- summary_df %>%
  group_by(major_division) %>%
  summarise(min_p = min(p_adj_BH), .groups = "drop") %>%
  arrange(min_p) %>%
  pull(major_division)

summary_df <- summary_df %>%
  mutate(major_division = factor(major_division, levels = rev(division_order)))

p <- ggplot(summary_df, aes(x = test, y = major_division, fill = neg_log10_p)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sig_label), color = "black", size = 5, vjust = 0.75) +
  scale_fill_gradient(low = "#FDF0F5", high = "#8E24AA", name = "-log10\n(BH-adj p)") +
  labs(
    title = "Summary: division-level effects\nacross the results section",
    subtitle = "Kruskal-Wallis omnibus tests, BH-FDR corrected\nwithin each test family (13 divisions each)\n* p<0.05  ** p<0.01  *** p<0.001",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, color = "grey30"),
    plot.margin = margin(10, 20, 10, 10)
  )

ggsave("results_summary_heatmap.png", p, width = 6.5, height = 7, dpi = 300)
ggsave("results_summary_heatmap.pdf", p, width = 6.5, height = 7)
cat("Saved: results_summary_heatmap.png / .pdf\n")

cat("\n=== Significant (BH-FDR<0.05) cells ===\n")
print(summary_df %>% filter(p_adj_BH < 0.05) %>%
        select(major_division, test, p_adj_BH) %>% arrange(test, p_adj_BH))


# ---------------------------------------------------------------------------
# correlation_heatmap_sex
# ---------------------------------------------------------------------------
# Full pairwise Spearman correlation heatmap (all values shown, not just
# thresholded edges), separately for Female and Male -- the sex-specific
# counterpart to nullmodel_corrected_correlation_heatmaps.R ("Figure C" in
# Results_Networks_Sex_Channel.docx), which did this per condition.
# Regions sorted by descending mean absolute correlation (the ordering
# rule), so the most interconnected regions cluster visually together.
#
# Curated region set (same ~80-region "regions of interest" set used by the
# per-condition correlation heatmaps and coactivation networks):
# top_regions_low_deviation.csv (high intensity, low CV) union
# most_variable_regions.csv (high CV among reasonably-signaled regions).

library(dplyr)
library(tidyr)
library(pheatmap)

setwd("~/Desktop/Bachelor")
source("prepare_heatmap_data.R")   # long (includes sex)
long <- long %>% filter(region_acronym != "FS")

top_regions_df <- read.csv("top_regions_low_deviation.csv", stringsAsFactors = FALSE)
most_var_df <- read.csv("most_variable_regions.csv", stringsAsFactors = FALSE)
node_regions <- setdiff(unique(c(top_regions_df$region_acronym, most_var_df$region_acronym)), fiber_tract_regions)
cat("n curated regions:", length(node_regions), "\n")

sex_order <- c("Female", "Male")

make_corr_matrix <- function(long_data, sx, min_obs = 3) {
  d <- long_data %>% filter(sex == sx, region_acronym %in% node_regions) %>%
    mutate(sample_id = paste(animal_id_original, channel, sep = "_"))
  mat <- d %>% select(sample_id, region_acronym, intensity) %>%
    pivot_wider(names_from = region_acronym, values_from = intensity) %>%
    tibble::column_to_rownames("sample_id") %>% as.matrix()
  keep <- colSums(!is.na(mat)) >= min_obs
  mat <- mat[, keep, drop = FALSE]
  cat(sprintf("[%s] n samples = %d, n regions (n>=%d) = %d\n", sx, nrow(mat), min_obs, ncol(mat)))

  r <- cor(mat, method = "spearman", use = "pairwise.complete.obs")
  mean_abs_r <- rowMeans(abs(r), na.rm = TRUE)
  ord <- order(mean_abs_r, decreasing = TRUE)
  r[ord, ord]
}

plot_corr_heatmap <- function(r, title, fname, max_n = 80) {
  if (nrow(r) > max_n) r <- r[1:max_n, 1:max_n]
  breaks <- seq(-1, 1, length.out = 101)
  colors <- colorRampPalette(c("#4A148C", "#FDF0F5", "#C2185B"))(100)
  pheatmap::pheatmap(r, cluster_rows = FALSE, cluster_cols = FALSE,
           color = colors, breaks = breaks,
           main = title, fontsize = 6, fontsize_row = 5, fontsize_col = 5,
           border_color = NA, filename = fname, width = 9, height = 8)
}

for (sx in sex_order) {
  r <- make_corr_matrix(long, sx)
  fname <- sprintf("correlation_heatmap_sex_%s.png", sx)
  plot_corr_heatmap(r, sprintf("Correlation heatmap (Spearman): %s", sx), fname)
  cat(sprintf("[%s] saved: %s (%d regions)\n", sx, fname, nrow(r)))
}

## --- Combined side-by-side PNG for easy comparison ---
library(magick)
imgs <- lapply(sex_order, function(sx) image_read(sprintf("correlation_heatmap_sex_%s.png", sx)))
combined <- image_append(image_scale(image_join(imgs), "1200"), stack = FALSE)
image_write(combined, "correlation_heatmap_sex_combined.png")

cat("\nSaved: correlation_heatmap_sex_Female.png, correlation_heatmap_sex_Male.png, correlation_heatmap_sex_combined.png\n")


# ---------------------------------------------------------------------------
# correlation_heatmap_split_sex
# ---------------------------------------------------------------------------
# Correlation DIFFERENCE heatmap: Male rho minus Female rho, on a shared
# region axis -- replaces the earlier split-triangle design (upper=Male,
# lower=Female) with a single matrix that directly plots the divergence
# itself, so sex-dependent reversals (e.g. POST/PRE/VISpor6b flipping sign
# between sexes during restraint) are immediately visible as strongly
# colored cells instead of requiring a mental "compare two triangles" step.
#
# Cells are masked (NA, grey) where either sex has fewer than
# min_pairwise_n valid paired observations for that region pair -- Spearman
# rho computed on <5 pairs is dominated by a handful of discrete achievable
# values (see e.g. n=3: rho in {-1,-0.5,0.5,1}) and produced spurious
# "extreme difference" cells in an earlier, unfiltered version of this
# analysis (regions with missing data silently dropping well below the
# condition's nominal n).
#
# Restricted to conditions with n>=7 in BOTH sexes (only restraint Acute
# qualifies; forced swim Acute F=3, forced swim Chronic M=0, social defeat
# 3/3, tail suspension 4/4).

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

setwd("~/Desktop/Bachelor")
source("prepare_heatmap_data.R")
long <- long %>% filter(region_acronym != "FS")

top_regions_df <- read.csv("top_regions_low_deviation.csv", stringsAsFactors = FALSE)
most_var_df <- read.csv("most_variable_regions.csv", stringsAsFactors = FALSE)
node_regions <- setdiff(unique(c(top_regions_df$region_acronym, most_var_df$region_acronym)), fiber_tract_regions)

min_pairwise_n <- 5

## --- Build a sample x region intensity matrix for one sex x condition ---
get_mat <- function(sx, cond, min_obs = 3) {
  d <- long %>% filter(sex == sx, condition == cond, region_acronym %in% node_regions) %>%
    mutate(sample_id = paste(animal_id_original, channel, sep = "_"))
  mat <- d %>% select(sample_id, region_acronym, intensity) %>%
    pivot_wider(names_from = region_acronym, values_from = intensity) %>%
    tibble::column_to_rownames("sample_id") %>% as.matrix()
  keep <- colSums(!is.na(mat)) >= min_obs
  mat[, keep, drop = FALSE]
}

## --- Build the Male-minus-Female difference matrix, masking pairs with
## insufficient pairwise n in either sex ---
build_diff_df <- function(mat_f, mat_m, min_n = min_pairwise_n) {
  shared <- intersect(colnames(mat_f), colnames(mat_m))
  mat_f <- mat_f[, shared, drop = FALSE]; mat_m <- mat_m[, shared, drop = FALSE]

  r_f <- cor(mat_f, method = "spearman", use = "pairwise.complete.obs")
  r_m <- cor(mat_m, method = "spearman", use = "pairwise.complete.obs")
  n_f <- crossprod(!is.na(mat_f))
  n_m <- crossprod(!is.na(mat_m))

  diff <- r_m - r_f
  reliable <- (n_f >= min_n) & (n_m >= min_n)
  diff[!reliable] <- NA
  diag(diff) <- NA

  mean_abs_diff <- rowMeans(abs(diff), na.rm = TRUE)
  ord <- shared[order(mean_abs_diff, decreasing = TRUE)]
  diff <- diff[ord, ord]

  df <- as.data.frame(as.table(diff)) %>%
    rename(region_row = Var1, region_col = Var2, value = Freq) %>%
    mutate(region_row = factor(region_row, levels = ord),
           region_col = factor(region_col, levels = ord))
  list(df = df, order = ord, n_reliable = sum(reliable[upper.tri(reliable)]),
       n_total = sum(upper.tri(reliable)))
}

plot_diff <- function(res, title, fname, max_n = 60, width = 9.5, height = 8, label_top = 10) {
  df <- res$df
  ord <- res$order
  if (length(ord) > max_n) {
    keep_regions <- ord[1:max_n]
    df <- df %>% filter(region_row %in% keep_regions, region_col %in% keep_regions) %>%
      mutate(region_row = droplevels(region_row), region_col = droplevels(region_col))
  }

  ## Label the top N most divergent pairs (upper triangle only, to avoid duplicates)
  df_upper <- df %>% filter(as.integer(region_row) < as.integer(region_col), !is.na(value))
  top_pairs <- df_upper %>% arrange(desc(abs(value))) %>% slice_head(n = label_top)

  p <- ggplot(df, aes(x = region_col, y = region_row, fill = value)) +
    geom_tile() +
    geom_text(data = top_pairs, aes(label = "×"), size = 1.6, color = "black") +
    scale_fill_gradientn(colours = c("#1B5E20", "#FDF0F5", "#B71C1C"), limits = c(-2, 2),
                          na.value = "grey85", name = "Male − Female\n(Spearman rho)") +
    scale_x_discrete(position = "top") +
    scale_y_discrete(limits = rev) +
    labs(title = title,
         subtitle = sprintf("Red = higher in Male | Green = higher in Female | grey = masked (n<%d in either sex) | × = top %d most divergent pairs",
                             min_pairwise_n, label_top),
         x = NULL, y = NULL) +
    coord_fixed() +
    theme_minimal(base_size = 6) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 4.2),
          axis.text.y = element_text(size = 4.2),
          plot.title = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 7),
          panel.grid = element_blank())
  ggsave(fname, p, width = width, height = height, units = "in", dpi = 300)
}

## ============================================================
## Difference heatmap(s) -- restricted to conditions with n>=7 in BOTH sexes.
## Checked all 5 conditions: forced swim Acute (F=3, M=10), forced swim
## Chronic (F=9, M=0), restraint Acute (F=7, M=11), social defeat Acute
## (F=3, M=3), tail suspension Acute (F=4, M=4) -- only restraint qualifies.
## ============================================================
condition_order_n_ok <- c("restraint · Acute")
cat("Difference heatmaps restricted to (n>=7 both sexes):", paste(condition_order_n_ok, collapse=", "), "\n")
cat("Excluded: forced swim Acute (F=3), forced swim Chronic (M=0), social defeat (3/3), tail suspension (4/4)\n\n")

for (cond in condition_order_n_ok) {
  mat_f <- get_mat("Female", cond)
  mat_m <- get_mat("Male", cond)
  res <- build_diff_df(mat_f, mat_m)
  fname <- sprintf("correlation_heatmap_diff_sex_%s.png", gsub("[^A-Za-z0-9]+", "_", cond))
  plot_diff(res, sprintf("Correlation difference (Male − Female) -- %s", cond), fname)
  cat(sprintf("Saved: %s (%d regions, %d/%d pairs reliable [n>=%d both sexes])\n",
              fname, length(res$order), res$n_reliable, res$n_total, min_pairwise_n))

  ## Print the top divergent pairs actually used for the labels
  df_upper <- res$df %>% filter(as.integer(region_row) < as.integer(region_col), !is.na(value))
  top10 <- df_upper %>% arrange(desc(abs(value))) %>% slice_head(n = 10)
  cat("\nTop 10 most divergent (reliable) region pairs:\n")
  print(as.data.frame(top10 %>% select(region_row, region_col, value)), row.names = FALSE)
}

cat("\nOnly one condition qualifies (restraint), so no multi-panel combination is produced.\n")
cat("Done.\n")


# ---------------------------------------------------------------------------
# nullmodel_corrected_correlation_heatmaps
# ---------------------------------------------------------------------------
# Correlation heatmaps, using a null-model-corrected approach (2026) "Estimation statistics and
# visualization": pairwise Spearman correlation matrices for the regions of
# interest, sorted by descending mean absolute correlation, one heatmap per
# stressor condition (her "Level 2: Stressor Specificity" comparison).
#
# Uses pheatmap (same package the reference implementation used). Built for BOTH the raw and the
# per-animal-normalized data, so the before/after contrast from the
# co-activation networks is visible here too (same node/region set as
# nullmodel_corrected_coactivation_networks(.R/_normalized.R) for consistency).

library(dplyr)
library(tidyr)
library(pheatmap)

df <- read.csv("~/Desktop/ALL_DT.csv", stringsAsFactors = FALSE)
cond_map <- read.csv("~/Desktop/group_condition_corrected.csv", stringsAsFactors = FALSE)

## --- Same node/region set as the network analysis ---
pairwise_stats <- read.csv("~/Desktop/nullmodel_corrected_pairwise_stats.csv", stringsAsFactors = FALSE)
top_by_pvalue <- pairwise_stats %>% group_by(region_acronym) %>%
  summarise(min_p = min(p_value), .groups = "drop") %>%
  arrange(min_p) %>% slice_head(n = 80) %>% pull(region_acronym)
top_regions_df <- read.csv("~/Desktop/top_regions_low_deviation.csv", stringsAsFactors = FALSE)
most_var_df <- read.csv("~/Desktop/most_variable_regions.csv", stringsAsFactors = FALSE)
node_regions <- unique(c(top_by_pvalue, top_regions_df$region_acronym, most_var_df$region_acronym))

base_long <- df %>%
  filter(region_acronym != "FS", region_acronym %in% node_regions) %>%
  select(animal_group, animal_id_original, region_acronym,
         trapcre_intensity, traptta_intensity) %>%
  pivot_longer(cols = c(trapcre_intensity, traptta_intensity),
               names_to = "channel", values_to = "value") %>%
  mutate(channel = ifelse(channel == "trapcre_intensity", "Cre", "tTA"),
         sample_id = paste(animal_id_original, channel, sep = "_")) %>%
  rename(Group = animal_group, Channel = channel) %>%
  left_join(cond_map, by = c("Group", "Channel")) %>%
  filter(!is.na(Condition))

animal_overall <- df %>%
  select(animal_id_original, trapcre_intensity, traptta_intensity) %>%
  pivot_longer(cols = c(trapcre_intensity, traptta_intensity), names_to = "channel", values_to = "value") %>%
  mutate(channel = ifelse(channel == "trapcre_intensity", "Cre", "tTA")) %>%
  group_by(animal_id_original, channel) %>%
  summarise(overall_median = median(value, na.rm = TRUE), .groups = "drop")

long_raw <- base_long
long_norm <- base_long %>%
  left_join(animal_overall, by = c("animal_id_original", "Channel" = "channel")) %>%
  mutate(value = value / overall_median)

conditions <- sort(unique(base_long$Condition))

## --- Build region x region Spearman correlation matrix, sorted by
## descending mean absolute correlation (the ordering rule) ---
make_corr_matrix <- function(long_data, cond, min_obs = 3) {
  d <- long_data %>% filter(Condition == cond)
  mat <- d %>% select(sample_id, region_acronym, value) %>%
    pivot_wider(names_from = region_acronym, values_from = value) %>%
    tibble::column_to_rownames("sample_id") %>% as.matrix()
  keep <- colSums(!is.na(mat)) >= min_obs
  mat <- mat[, keep, drop = FALSE]

  r <- cor(mat, method = "spearman", use = "pairwise.complete.obs")
  mean_abs_r <- rowMeans(abs(r), na.rm = TRUE)
  ord <- order(mean_abs_r, decreasing = TRUE)
  r[ord, ord]
}

## --- Plot heatmap (diverging pink-violet scale, consistent with the thesis) ---
plot_corr_heatmap <- function(r, title, fname, max_n = 80) {
  if (nrow(r) > max_n) r <- r[1:max_n, 1:max_n]   # cap size for legibility
  breaks <- seq(-1, 1, length.out = 101)
  colors <- colorRampPalette(c("#4A148C", "#FDF0F5", "#C2185B"))(100)

  pheatmap(r, cluster_rows = FALSE, cluster_cols = FALSE,
           color = colors, breaks = breaks,
           main = title, fontsize = 6, fontsize_row = 4, fontsize_col = 4,
           border_color = NA, filename = fname, width = 9, height = 8)
}

for (cond in conditions) {
  fname_safe <- gsub("[^A-Za-z0-9]+", "_", cond)

  r_raw <- make_corr_matrix(long_raw, cond)
  plot_corr_heatmap(r_raw, sprintf("Correlation heatmap (RAW): %s", cond),
                     sprintf("~/Desktop/nullmodel_corrected_corrheatmap_raw_%s.png", fname_safe))

  r_norm <- make_corr_matrix(long_norm, cond)
  plot_corr_heatmap(r_norm, sprintf("Correlation heatmap (normalized): %s", cond),
                     sprintf("~/Desktop/nullmodel_corrected_corrheatmap_normalized_%s.png", fname_safe))

  cat(sprintf("[%s] raw: %d regions | normalized: %d regions -- saved\n", cond, nrow(r_raw), nrow(r_norm)))
}

