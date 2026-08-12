# =============================================================================
# Region- and division-level statistics between stressor conditions
#
# Consolidated from the following original scripts:
#   - condition_kruskal_top18_table.R
#   - condition_matrix_stats.R
#   - division_level_anova.R
#   - effect_sizes_pairwise.R
#   - curated_regions_all_channels.R
#   - curated_regions_all_channels_raw.R
#   - curated_regions_pairwise_wilcoxon.R
#   - curated_regions_pairwise_wilcoxon_raw.R
#   - nullmodel_corrected_pairwise_stats.R
#   - condition_specific_regions.R
#   - condition_specific_regions_normalized.R
#   - test_condition_pattern_similarity.R
#   - fs_acute_vs_chronic_pooled.R
#   - fs_acute_vs_chronic_female_only.R
# =============================================================================


# ---------------------------------------------------------------------------
# condition_kruskal_top18_table
# ---------------------------------------------------------------------------
# Pretty table (gt) of the 18 regions with raw p < 0.05 in the region-level
# condition Kruskal-Wallis test, styled with the same pink-to-violet palette
# used throughout the thesis heatmaps (prepare_heatmap_data.R col_fun).

library(dplyr)
library(gt)

setwd("~/Desktop/Bachelor")
suppressPackageStartupMessages(source("prepare_heatmap_data.R"))

kw <- read.csv("condition_kruskal_by_region.csv")
top18 <- kw %>%
  filter(p_value < 0.05) %>%
  arrange(p_value) %>%
  head(18) %>%
  mutate(major_division = major_division_per_row[match(region_acronym, rownames(mat))]) %>%
  select(region_acronym, major_division, chi_sq, p_value)

tbl <- gt(top18) %>%
  cols_label(region_acronym = "Region", major_division = "Major division",
             chi_sq = md("&chi;&sup2;"), p_value = "p-value") %>%
  fmt_number(columns = chi_sq, decimals = 2) %>%
  fmt_number(columns = p_value, decimals = 4) %>%
  data_color(columns = chi_sq, colors = scales::col_numeric(
    palette = c("#FDF0F5", "#FBC9DE", "#F48FB1", "#E05780", "#C2185B", "#8E24AA", "#4A148C"),
    domain = range(top18$chi_sq))) %>%
  tab_style(style = cell_text(color = "white"),
            locations = cells_body(columns = chi_sq, rows = chi_sq > 10.7)) %>%
  tab_header(title = "Region-level Kruskal-Wallis: regions with raw p < 0.05",
              subtitle = "Condition effect across 670 regions (5 conditions); none survive BH-FDR correction") %>%
  tab_source_note(source_note = "Kruskal-Wallis rank-sum test, df = 4. None of the 18 regions survive Benjamini-Hochberg correction (min BH q = 0.42).") %>%
  tab_options(table.font.size = 12, heading.title.font.size = 14, heading.subtitle.font.size = 11,
              column_labels.font.weight = "bold", table.width = pct(85))

gtsave(tbl, "condition_kruskal_top18_table.html")
cat("saved html\n")


# ---------------------------------------------------------------------------
# condition_matrix_stats
# ---------------------------------------------------------------------------
# Statistics for the condition-level heatmap matrix (prepare_heatmap_data.R):
#   1. Kruskal-Wallis rank-sum test across all 5 conditions at once -- the
#      multi-group generalization of Wilcoxon (Wilcoxon itself only compares
#      two groups; Kruskal-Wallis is what you get testing "any difference
#      among the 5 conditions" in one test instead of 10 pairwise Wilcoxon
#      comparisons). Run at three resolutions, mirroring the
#      whole-brain/division/region hierarchy used elsewhere in this project
#      (division_level_anova.R, channel_bias_check.R):
#      1a. whole-brain level (1 test; each animal x channel's own
#          whole-brain median intensity vs Condition)
#      1b. division level (14 tests, one per major anatomical division;
#          whole-brain-median-normalized per animal x channel, same
#          normalization as division_level_anova.R, to remove global
#          scaling before testing division-specific effects)
#      1c. region level (669 tests, one per region), pooling all
#          animals/channels sharing that condition -- same pooling used to
#          build `mat`.
#      BH-FDR applied within each resolution separately (1, 14, and 669
#      tests respectively).
#   2. Row-wise z-scores of `mat` (per region, across the 5 conditions) --
#      puts every region on the same scale regardless of its absolute
#      intensity, so cross-condition patterns become comparable.
#   3. Hierarchical clustering (Ward D2, Euclidean) on the z-scored profiles:
#      of the 5 conditions (which stressors evoke similar whole-brain
#      patterns) and of the 670 regions (which regions respond similarly
#      across conditions; cut into k region modules).
#   4. Channel (Cre vs tTA) difference tests:
#      4a. condition-independent, paired by animal, whole-brain + per-region
#          -- same design as channel_bias_check.R (all 27 animals; this is
#          NOT re-run here since it never touched the condition columns and
#          so is unaffected by today's ALL_DT.csv condition-mapping fix --
#          see channel_bias_wholebrain.csv / channel_bias_per_region.csv).
#      4b. within-condition, unpaired, per region -- only for the 3 of 5
#          conditions where BOTH channels are actually represented (forced
#          swim chronic is Cre-only, social defeat is tTA-only, so no
#          channel contrast is possible there).

library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)

setwd("~/Desktop/Bachelor")
source("prepare_heatmap_data.R")   # mat, long, condition_order, major_division_per_row

## "FS" (injection site) is trivially constant across conditions (median
## ratio = 1.0 by construction) -- excluded here, consistent with
## division_level_anova.R and friends, and because a zero-variance row
## breaks z-scoring/clustering.
keep <- rownames(mat) != "FS"
major_division_per_row <- major_division_per_row[keep]
mat <- mat[keep, , drop = FALSE]
long <- long %>% filter(region_acronym != "FS")
regions <- rownames(mat)

## ============================================================
## 1a. Whole-brain-level Kruskal-Wallis (1 test)
## ============================================================
animal_overall <- long %>%
  group_by(animal_id_original, channel, condition) %>%
  summarise(overall_median = median(intensity, na.rm = TRUE), .groups = "drop")

kt_wb <- kruskal.test(overall_median ~ condition, data = animal_overall)
cat("=== 1a. Whole-brain-level Kruskal-Wallis (all 5 conditions) ===\n")
cat(sprintf("n animal x channel samples: %d | chi-sq = %.3f, df = %d, p = %.4f\n",
            nrow(animal_overall), kt_wb$statistic, kt_wb$parameter, kt_wb$p.value))
print(animal_overall %>% group_by(condition) %>%
        summarise(n = n(), median_overall = median(overall_median), .groups = "drop"))

write.csv(
  animal_overall %>% group_by(condition) %>%
    summarise(n = n(), median_overall = median(overall_median), .groups = "drop") %>%
    mutate(chi_sq = unname(kt_wb$statistic), df = unname(kt_wb$parameter), p_value = kt_wb$p.value),
  "condition_kruskal_wholebrain.csv", row.names = FALSE
)

## ============================================================
## 1b. Division-level Kruskal-Wallis (14 tests)
## ============================================================
## Whole-brain-median normalization per animal x channel (removes global
## scaling before testing division-specific effects), same approach as
## division_level_anova.R.
long_norm <- long %>%
  left_join(animal_overall, by = c("animal_id_original", "channel", "condition")) %>%
  mutate(norm_value = intensity / overall_median)

division_scores <- long_norm %>%
  group_by(animal_id_original, channel, condition, major_division = major_division_per_row[match(region_acronym, rownames(mat))]) %>%
  summarise(div_value = median(norm_value, na.rm = TRUE), n_regions = n(), .groups = "drop") %>%
  filter(!is.na(major_division)) %>%
  mutate(log_div_value = log10(div_value + 0.001))

divisions <- unique(division_scores$major_division)
division_results <- vector("list", length(divisions))

for (i in seq_along(divisions)) {
  dv <- divisions[i]
  d <- division_scores %>% filter(major_division == dv)
  if (n_distinct(d$condition) < 2) next
  kt <- kruskal.test(log_div_value ~ condition, data = d)
  medians <- d %>% group_by(condition) %>% summarise(m = median(log_div_value), .groups = "drop")
  division_results[[i]] <- data.frame(
    major_division = dv, chi_sq = unname(kt$statistic), df = unname(kt$parameter),
    n_regions_pooled = round(mean(d$n_regions)), p_value = kt$p.value,
    t(setNames(medians$m, paste0("median_log_", medians$condition)))
  )
}

division_kruskal_df <- bind_rows(division_results) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  arrange(p_value)

write.csv(division_kruskal_df, "condition_kruskal_by_division.csv", row.names = FALSE)

cat("\n=== 1b. Division-level Kruskal-Wallis (all 5 conditions), n =", nrow(division_kruskal_df), "divisions ===\n")
print(as.data.frame(division_kruskal_df[, c("major_division", "chi_sq", "df", "p_value", "p_adj_BH")]))
cat(sprintf("\nDivisions with raw p < 0.05: %d / %d\n", sum(division_kruskal_df$p_value < 0.05), nrow(division_kruskal_df)))
cat(sprintf("Divisions with BH-FDR p < 0.05: %d / %d\n", sum(division_kruskal_df$p_adj_BH < 0.05), nrow(division_kruskal_df)))

## ============================================================
## 1c. Region-level Kruskal-Wallis (669 tests)
## ============================================================
kruskal_results <- vector("list", length(regions))

for (i in seq_along(regions)) {
  reg <- regions[i]
  d <- long %>% filter(region_acronym == reg)
  if (n_distinct(d$condition) < 2 || nrow(d) < 5) next
  kt <- suppressWarnings(kruskal.test(intensity ~ condition, data = d))
  medians <- d %>% group_by(condition) %>% summarise(m = median(intensity, na.rm = TRUE), .groups = "drop")
  kruskal_results[[i]] <- data.frame(
    region_acronym = reg, chi_sq = unname(kt$statistic), df = unname(kt$parameter),
    n_total = nrow(d), p_value = kt$p.value,
    t(setNames(medians$m, paste0("median_", medians$condition)))
  )
}

kruskal_df <- bind_rows(kruskal_results) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj_BH, p_value)

write.csv(kruskal_df, "condition_kruskal_by_region.csv", row.names = FALSE)

cat("\n=== 1c. Region-level Kruskal-Wallis (all 5 conditions at once), per region ===\n")
cat("Total tests:", nrow(kruskal_df), "\n")
cat("Raw p < 0.05:", sum(kruskal_df$p_value < 0.05, na.rm = TRUE), "\n")
cat("BH-FDR p < 0.05:", sum(kruskal_df$p_adj_BH < 0.05, na.rm = TRUE), "\n")
cat("\nTop 15 by BH-adjusted p-value:\n")
print(head(as.data.frame(kruskal_df), 15))

## ============================================================
## 1d. Pairwise Wilcoxon rank-sum tests between all 10 condition pairs,
##     per region -- post-hoc detail complementing the omnibus
##     Kruskal-Wallis test above (WHICH specific condition pairs drive a
##     region's difference). BH-FDR applied GLOBALLY across all tests.
## ============================================================
cond_pairs <- t(combn(condition_order, 2))
pairwise_results <- vector("list", length(regions) * nrow(cond_pairs))
idx <- 1

for (reg in regions) {
  d <- long %>% filter(region_acronym == reg)
  for (p in seq_len(nrow(cond_pairs))) {
    c1 <- cond_pairs[p, 1]; c2 <- cond_pairs[p, 2]
    v1 <- d %>% filter(condition == c1) %>% pull(intensity)
    v2 <- d %>% filter(condition == c2) %>% pull(intensity)
    if (length(v1) < 2 || length(v2) < 2) { idx <- idx + 1; next }
    wt <- suppressWarnings(wilcox.test(v1, v2))
    pairwise_results[[idx]] <- data.frame(
      region_acronym = reg, condition1 = c1, condition2 = c2,
      n1 = length(v1), n2 = length(v2),
      median1 = median(v1), median2 = median(v2),
      p_value = wt$p.value
    )
    idx <- idx + 1
  }
}

pairwise_df <- bind_rows(pairwise_results) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj_BH, p_value)

write.csv(pairwise_df, "condition_pairwise_wilcoxon_by_region.csv", row.names = FALSE)

cat("\n=== 1d. Pairwise Wilcoxon (10 condition pairs), per region ===\n")
cat("Total tests:", nrow(pairwise_df), "\n")
cat("Raw p < 0.05:", sum(pairwise_df$p_value < 0.05, na.rm = TRUE), "\n")
cat("BH-FDR p < 0.05:", sum(pairwise_df$p_adj_BH < 0.05, na.rm = TRUE), "\n")
cat("\nTop 15 by BH-adjusted p-value:\n")
print(head(as.data.frame(pairwise_df), 15))
cat("\nBreakdown of significant pairs (raw p<0.05) by condition pair:\n")
print(pairwise_df %>% filter(p_value < 0.05) %>% count(condition1, condition2, sort = TRUE))

## ============================================================
## 2. Row-wise z-scores of the condition matrix
## ============================================================
mat_z <- t(scale(t(mat)))
write.csv(data.frame(region_acronym = rownames(mat_z), mat_z, check.names = FALSE),
          "condition_zscores_per_region.csv", row.names = FALSE)
cat(sprintf("\n=== 2. Row-wise z-scores computed for %d regions x %d conditions ===\n",
            nrow(mat_z), ncol(mat_z)))

## ============================================================
## 3. Hierarchical clustering (Ward D2, Euclidean) on z-scored profiles
## ============================================================
hc_conditions <- hclust(dist(t(mat_z), method = "euclidean"), method = "ward.D2")
hc_regions <- hclust(dist(mat_z, method = "euclidean"), method = "ward.D2")

k_regions <- 6   # arbitrary, adjustable -- inspect condition_region_clustering.pdf's dendrogram to reconsider
region_clusters <- cutree(hc_regions, k = k_regions)
region_cluster_df <- data.frame(
  region_acronym = names(region_clusters),
  cluster = region_clusters,
  major_division = major_division_per_row[match(names(region_clusters), rownames(mat))]
)
write.csv(region_cluster_df, "condition_region_clusters.csv", row.names = FALSE)

cat("\n=== 3. Condition dendrogram order (Ward D2 on z-scored whole-brain profiles) ===\n")
print(hc_conditions$labels[hc_conditions$order])

cat(sprintf("\nRegion clusters (k = %d), sizes:\n", k_regions))
print(table(region_clusters))

pdf("condition_dendrogram.pdf", width = 5, height = 4)
plot(hc_conditions, main = "Condition clustering (Ward D2, z-scored whole-brain profiles)",
     xlab = "", sub = "", ylab = "Height")
dev.off()

## Same pink-to-violet family as the rest of the thesis (Figure 7 etc.):
## pink = below-average (negative z), white = average, violet = above-average
col_fun_z <- colorRamp2(c(-2, 0, 2), c("#C2185B", "white", "#4A148C"))
## Rebuild the division row annotation scoped to mat_z's 669 rows (FS
## excluded) -- the row_anno from prepare_heatmap_data.R still has 670 rows
## and would mismatch. fiber tracts already excluded upstream (mat comes
## from prepare_heatmap_data.R), so drop the now-unused legend entry.
major_cols_no_ft <- major_cols[names(major_cols) != "fiber tracts"]
row_anno_z <- rowAnnotation(
  division = major_division_per_row,
  col = list(division = major_cols_no_ft),
  annotation_name_gp = gpar(fontsize = 12),
  annotation_legend_param = list(division = list(
    ncol = 1, title = "Major division",
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 12),
    grid_height = unit(4.5, "mm"), grid_width = unit(4.5, "mm")
  ))
)

## Unify cluster tree heights (reviewer comment, same fix as
## plot_heatmap_ALL_DT_v2.R): with row_split cutting ONE shared dendrogram
## into k_regions groups, ComplexHeatmap by default draws each group's
## sub-dendrogram scaled to THAT GROUP's own local max height, filling the
## same panel width regardless of the group's true (dis)similarity to the
## others -- verified empirically (per-slice max heights ranged 2.95-7.08,
## none matching the full tree's max of ~17), so the panel widths as
## rendered do NOT reflect true relative cluster heterogeneity.
##
## Unlike ALL_DT (where each division was an independently-run clustering
## with no shared scale to begin with), these 6 groups come from cutting
## ONE clustering, so their relative heights ARE scientifically meaningful
## -- rescaling each group independently to [0,1] would erase that. Instead
## rescale all groups by ONE shared constant (the largest per-group max
## height), preserving true relative comparability while still using the
## full plot width for the most heterogeneous group.
## group assignment stays IDENTICAL to before (same cutree call), only the
## rendered height axis changes.
region_clusters_f <- factor(region_clusters)
group_max_heights <- sapply(split(rownames(mat_z), region_clusters_f), function(regs) {
  if (length(regs) < 2) return(0)
  max(hclust(dist(mat_z[regs, , drop = FALSE]), method = "ward.D2")$height)
})
shared_max_height <- max(group_max_heights)

cluster_rows_shared_scale <- function(m) {
  if (nrow(m) < 2) return(NULL)
  hc <- hclust(dist(m), method = "ward.D2")
  hc$height <- hc$height / shared_max_height
  as.dendrogram(hc)
}

ht_z <- Heatmap(
  mat_z, name = "z-score", col = col_fun_z,
  cluster_rows = cluster_rows_shared_scale, row_split = region_clusters_f,
  cluster_row_slices = FALSE, show_row_names = FALSE,
  cluster_columns = hc_conditions, show_column_names = TRUE, column_names_rot = 0,
  column_names_gp = gpar(fontsize = 12),
  row_title = NULL, right_annotation = row_anno_z,
  heatmap_legend_param = list(title_gp = gpar(fontsize = 12, fontface = "bold"), labels_gp = gpar(fontsize = 12),
                               grid_height = unit(4.5, "mm"), grid_width = unit(4.5, "mm"))
)
pdf("condition_zscore_heatmap.pdf", width = 10, height = 16)
draw(ht_z)
dev.off()
png("condition_zscore_heatmap.png", width = 10, height = 16, units = "in", res = 200)
draw(ht_z)
dev.off()

## ============================================================
## 3b. Hierarchical clustering of regions WITHIN each condition separately
## ============================================================
## Unlike section 3 (one clustering of regions using their profile ACROSS
## all 5 conditions), this clusters regions independently per condition,
## using that condition's own animals as the clustering dimension (rows =
## regions, columns = the n_mice animals tagged with that condition). Each
## region's row is z-scored across just those animals before clustering, so
## results reflect within-condition covariance structure only, not
## cross-condition differences.
##
## Restricted to the top_n_regions highest-intensity regions per condition
## (by median across that condition's animals) so leaves can be labeled with
## region acronyms and stay legible on a normal page -- with ~600 regions,
## a fully labeled tree isn't readable at any reasonable font size.
top_n_regions <- 40      # adjustable
k_regions_per_cond <- 5  # adjustable; fewer than section 3's k=6 since far fewer leaves

per_condition_clusters <- vector("list", length(condition_order))
hc_by_condition <- vector("list", length(condition_order))
names(per_condition_clusters) <- condition_order
names(hc_by_condition) <- condition_order

for (cond in condition_order) {
  d <- long %>% filter(condition == cond) %>%
    select(region_acronym, animal_id_original, intensity)
  wide_cond <- d %>% pivot_wider(names_from = animal_id_original, values_from = intensity)
  mat_cond <- as.matrix(wide_cond[, -1])
  rownames(mat_cond) <- wide_cond$region_acronym
  ## Keep all animals; drop only regions missing a value for at least one of
  ## this condition's animals (most animals are missing scattered regions,
  ## not each other -- dropping incomplete animals instead would leave 0).
  mat_cond <- mat_cond[rowSums(is.na(mat_cond)) == 0, , drop = FALSE]

  ## Keep only the top_n_regions highest-intensity regions (median across
  ## this condition's animals), consistent with "highest intensity" as the
  ## selection criterion.
  med_per_region <- apply(mat_cond, 1, median, na.rm = TRUE)
  top_regions <- names(sort(med_per_region, decreasing = TRUE))[seq_len(min(top_n_regions, length(med_per_region)))]
  mat_cond_top <- mat_cond[top_regions, , drop = FALSE]

  mat_cond_z <- t(scale(t(mat_cond_top)))
  mat_cond_z <- mat_cond_z[rowSums(is.na(mat_cond_z)) == 0, , drop = FALSE]  # drop zero-variance regions

  hc_cond <- hclust(dist(mat_cond_z, method = "euclidean"), method = "ward.D2")
  clusters_cond <- cutree(hc_cond, k = min(k_regions_per_cond, nrow(mat_cond_z) - 1))

  hc_by_condition[[cond]] <- hc_cond
  per_condition_clusters[[cond]] <- data.frame(
    condition = cond, region_acronym = names(clusters_cond), cluster = clusters_cond,
    median_intensity = med_per_region[names(clusters_cond)],
    major_division = major_division_per_row[match(names(clusters_cond), rownames(mat))]
  )
  cat(sprintf("\n%s: %d animals, top %d regions clustered, cluster sizes: %s\n",
              cond, ncol(mat_cond), nrow(mat_cond_z), paste(table(clusters_cond), collapse = ", ")))
}

per_condition_clusters_df <- bind_rows(per_condition_clusters) %>% arrange(condition, desc(median_intensity))
write.csv(per_condition_clusters_df, "condition_region_clusters_within_condition.csv", row.names = FALSE)

pdf("condition_region_dendrograms_by_condition.pdf", width = 7, height = 9)
for (cond in condition_order) {
  dend <- as.dendrogram(hc_by_condition[[cond]])
  par(mar = c(4, 2, 3, 10))  # room on the right for horizontal leaf labels
  plot(dend, horiz = TRUE, main = sprintf("%s\n(top %d regions by median intensity)", cond, top_n_regions),
       xlab = "Height", cex.main = 0.9)
}
dev.off()

## Combined single-page version: all 5 conditions in a 2x3 grid (one cell
## left empty), sized to fit a standard portrait thesis page.
## Leaf labels colored by major anatomical division (same palette as
## Figures 7/8), so division membership no longer needs to be spelled out
## in parentheses in the text.
library(dendextend)

## Bootstrap stability (ARI, mean +/- SD) per condition, from
## condition_per_condition_cluster_stability.R -- annotated onto each panel
## if that script has already been run; skipped gracefully otherwise.
stability_path <- "condition_per_condition_cluster_stability.csv"
stability_df <- if (file.exists(stability_path)) read.csv(stability_path) else NULL

panel_letters <- setNames(LETTERS[seq_along(condition_order)], condition_order)

## Leaf cex=1.3 needs ~40*1.3*12/72 = 8.7in of plot height per panel (40
## leaves/panel); with 2 rows and ~2.4in of margin per panel, a 34in-tall
## canvas gives ~14.6in per row -- comfortably enough headroom to avoid
## label collisions while still ~37% larger than the original cex=0.95.
png("condition_region_dendrograms_combined.png", width = 16, height = 34, units = "in", res = 300)
par(mfrow = c(2, 3))
for (cond in condition_order) {
  dend <- as.dendrogram(hc_by_condition[[cond]])
  leaf_divs <- major_division_per_row[match(labels(dend), rownames(mat))]
  labels_colors(dend) <- major_cols[as.character(leaf_divs)]
  dend <- dendextend::set(dend, "labels_cex", 1.3)

  stab_line <- ""
  if (!is.null(stability_df) && cond %in% stability_df$condition) {
    s <- stability_df[stability_df$condition == cond, ]
    stab_line <- sprintf("\nARI = %.2f +/- %.2f", s$mean_ari, s$sd_ari)
  }

  par(mar = c(4, 2.5, 7, 4), cex = 1.3)
  plot(dend, horiz = TRUE,
       main = sprintf("(%s) %s\n(top %d regions)%s", panel_letters[cond], cond, top_n_regions, stab_line),
       xlab = "Height", cex.main = 1.4, cex.axis = 1.2, cex.lab = 1.3)
}
## 6th grid cell: legend for the division color code
par(mar = c(2, 0, 0, 0), cex = 1.3)
plot.new()
legend("center", legend = names(major_cols), fill = major_cols, cex = 1.4,
       bty = "n", title = sprintf("(%s) Major division", LETTERS[length(condition_order) + 1]), ncol = 1)
dev.off()
par(mfrow = c(1, 1), cex = 1)

cat("\n=== 3b. Per-condition region clustering (top", top_n_regions, "regions) saved ===\n")
cat("  condition_region_clusters_within_condition.csv\n")
cat("  condition_region_dendrograms_by_condition.pdf\n")

## ============================================================
## 4. Channel (Cre vs tTA) difference tests
## ============================================================
## 4a. Whole-brain, paired by animal, condition-independent -- reported from
##     the existing channel_bias_check.R output (unaffected by the condition
##     column fix; not recomputed here).
if (file.exists("channel_bias_wholebrain.csv")) {
  wb <- read.csv("channel_bias_wholebrain.csv")
  wt_wb <- wilcox.test(wb$Cre, wb$tTA, paired = TRUE)
  cat(sprintf("\n=== 4a. Whole-brain Cre vs tTA (paired by animal, n=%d, condition-independent) ===\n", nrow(wb)))
  cat(sprintf("Wilcoxon signed-rank: V = %.1f, p = %.4f (from channel_bias_check.R design)\n",
              wt_wb$statistic, wt_wb$p.value))
} else {
  cat("\n=== 4a. Run channel_bias_check.R first for the condition-independent whole-brain/per-region channel test ===\n")
}

## 4b. Within-condition, unpaired Cre-vs-tTA per region -- only where both
##     channels are represented for that condition.
channel_counts <- long %>% distinct(condition, channel) %>% count(condition)
conditions_both_channels <- as.character(channel_counts$condition[channel_counts$n == 2])
cat("\nConditions with both channels represented:", paste(conditions_both_channels, collapse = ", "), "\n")
cat("Conditions with only one channel (no contrast possible):",
    paste(setdiff(as.character(condition_order), conditions_both_channels), collapse = ", "), "\n")

channel_diff_results <- vector("list", length(conditions_both_channels) * length(regions))
idx <- 1
for (cond in conditions_both_channels) {
  for (reg in regions) {
    d <- long %>% filter(region_acronym == reg, condition == cond)
    v1 <- d %>% filter(channel == "Cre") %>% pull(intensity)
    v2 <- d %>% filter(channel == "tTA") %>% pull(intensity)
    if (length(v1) < 2 || length(v2) < 2) { idx <- idx + 1; next }
    wt <- suppressWarnings(wilcox.test(v1, v2))
    channel_diff_results[[idx]] <- data.frame(
      condition = cond, region_acronym = reg,
      n_Cre = length(v1), n_tTA = length(v2),
      median_Cre = median(v1), median_tTA = median(v2),
      p_value = wt$p.value
    )
    idx <- idx + 1
  }
}

channel_diff_df <- bind_rows(channel_diff_results) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj_BH, p_value)
write.csv(channel_diff_df, "condition_channel_diff_by_region.csv", row.names = FALSE)

cat("\n=== 4b. Within-condition Cre vs tTA per region (unpaired) ===\n")
cat("Total tests:", nrow(channel_diff_df), "\n")
cat("Raw p < 0.05:", sum(channel_diff_df$p_value < 0.05, na.rm = TRUE), "\n")
cat("BH-FDR p < 0.05:", sum(channel_diff_df$p_adj_BH < 0.05, na.rm = TRUE), "\n")
cat("\nTop 10 by BH-adjusted p-value:\n")
print(head(as.data.frame(channel_diff_df), 10))

cat("\nSaved:\n")
cat("  condition_kruskal_wholebrain.csv, condition_kruskal_by_division.csv, condition_kruskal_by_region.csv\n")
cat("  condition_pairwise_wilcoxon_by_region.csv\n")
cat("  condition_zscores_per_region.csv\n")
cat("  condition_region_clusters.csv\n")
cat("  condition_dendrogram.pdf, condition_zscore_heatmap.pdf\n")
cat("  condition_region_clusters_within_condition.csv, condition_region_dendrograms_by_condition.pdf\n")
cat("  condition_channel_diff_by_region.csv\n")


# ---------------------------------------------------------------------------
# division_level_anova
# ---------------------------------------------------------------------------
# Division-level ANOVA: aggregate the ~660 individual brain regions up to
# the 14 major anatomical divisions (Allen CCF top-level categories, same
# derivation as prepare_heatmap_data.R / coactive_projection_regions.R),
# then run one-way ANOVA (value ~ Condition) at this much coarser
# resolution. Rationale: with n = 660 regions, BH-FDR correction is very
# harsh (0 regions survived at region-level, see
# condition_specific_regions_normalized_anova.csv). Collapsing to 14
# divisions cuts the multiple-testing burden by ~47x and pools many
# region-level samples into one division-level value per animal x channel,
# which increases the effective within-group N and could reveal genuine
# division-wide effects that are invisible at single-region resolution.
#
# Aggregation: for each animal x channel, the division score = median of
# that animal's whole-brain-median-normalized values across all regions
# belonging to that division (median chosen for robustness, consistent with
# the median-first philosophy used throughout this project).
# 'FS' excluded (injection region, trivially 1.0).

library(dplyr)
library(tidyr)
library(ggplot2)

df <- read.csv("~/Desktop/ALL_DT.csv", stringsAsFactors = FALSE)
cond_map <- read.csv("~/Desktop/group_condition_corrected.csv", stringsAsFactors = FALSE)
structures <- read.csv("~/Desktop/tractquant/allen_mouse_25um_v1.2/structures.csv", stringsAsFactors = FALSE)

## --- Major anatomical division (Allen CCF), same derivation used
## throughout the project ---
major_ids <- c(
  "315"  = "Isocortex", "698"  = "Olfactory areas", "1089" = "Hippocampal formation",
  "703"  = "Cortical subplate", "477" = "Striatum", "803" = "Pallidum",
  "549"  = "Thalamus", "1097" = "Hypothalamus", "313" = "Midbrain",
  "771"  = "Pons", "354" = "Medulla", "512" = "Cerebellum",
  "1009" = "fiber tracts", "73" = "ventricular systems"
)
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

long <- df %>%
  filter(region_acronym != "FS") %>%
  select(animal_group, animal_id_original, region_id, region_acronym,
         trapcre_intensity, traptta_intensity) %>%
  left_join(region_lookup %>% select(region_id, major_division), by = "region_id") %>%
  filter(!is.na(major_division)) %>%
  pivot_longer(cols = c(trapcre_intensity, traptta_intensity),
               names_to = "channel", values_to = "value") %>%
  mutate(channel = ifelse(channel == "trapcre_intensity", "Cre", "tTA"),
         sample_id = paste(animal_id_original, channel, sep = "_")) %>%
  rename(Group = animal_group, Channel = channel) %>%
  left_join(cond_map, by = c("Group", "Channel")) %>%
  filter(!is.na(Condition))

## --- Whole-brain-median normalization (per animal x channel) ---
animal_overall <- df %>%
  select(animal_id_original, trapcre_intensity, traptta_intensity) %>%
  pivot_longer(cols = c(trapcre_intensity, traptta_intensity), names_to = "channel", values_to = "value") %>%
  mutate(channel = ifelse(channel == "trapcre_intensity", "Cre", "tTA")) %>%
  group_by(animal_id_original, channel) %>%
  summarise(overall_median = median(value, na.rm = TRUE), .groups = "drop")

long <- long %>%
  left_join(animal_overall, by = c("animal_id_original", "Channel" = "channel")) %>%
  mutate(norm_value = value / overall_median)

## --- Aggregate: one division-level score per animal x channel x division
## (median across all regions belonging to that division) ---
division_scores <- long %>%
  group_by(sample_id, Condition, major_division) %>%
  summarise(div_value = median(norm_value, na.rm = TRUE), n_regions = n(), .groups = "drop") %>%
  mutate(log_div_value = log10(div_value + 0.001))

cat("n animal x channel samples:", n_distinct(division_scores$sample_id), "\n")
cat("n divisions tested:", n_distinct(division_scores$major_division), "\n")
print(division_scores %>% distinct(sample_id, Condition) %>% count(Condition))

## --- One-way ANOVA per division: log_div_value ~ Condition ---
divisions <- unique(division_scores$major_division)
results <- vector("list", length(divisions))

for (i in seq_along(divisions)) {
  dv <- divisions[i]
  d <- division_scores %>% filter(major_division == dv)
  if (n_distinct(d$Condition) < 2) next
  fit <- aov(log_div_value ~ Condition, data = d)
  ss <- summary(fit)[[1]]
  ss_cond <- ss["Condition", "Sum Sq"]; ss_resid <- ss["Residuals", "Sum Sq"]
  eta2 <- ss_cond / (ss_cond + ss_resid)
  pval <- ss["Condition", "Pr(>F)"]

  ## Tukey post-hoc, to see WHICH condition pairs drive any division-level effect
  tuk <- TukeyHSD(fit)
  tuk_df <- as.data.frame(tuk$Condition)
  tuk_df$comparison <- rownames(tuk_df)
  best_pair <- tuk_df$comparison[which.min(tuk_df$`p adj`)]
  best_pair_p <- min(tuk_df$`p adj`)

  results[[i]] <- data.frame(division = dv, eta_sq = eta2, p_value = pval,
                              n_regions_pooled = round(mean(d$n_regions)),
                              best_pair = best_pair, best_pair_p_adj = best_pair_p)
}

anova_results <- bind_rows(results) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  arrange(p_value)

write.csv(anova_results, "~/Desktop/division_level_anova.csv", row.names = FALSE)

cat("\n=== Division-level one-way ANOVA (Condition effect), n =", nrow(anova_results), "divisions ===\n")
print(as.data.frame(anova_results))
cat(sprintf("\nDivisions with raw p < 0.05: %d / %d\n", sum(anova_results$p_value < 0.05), nrow(anova_results)))
cat(sprintf("Divisions with BH-FDR p < 0.05: %d / %d\n", sum(anova_results$p_adj_BH < 0.05), nrow(anova_results)))

## --- Plot: division-level values by condition, ordered by ANOVA p-value ---
div_order <- anova_results$division
plot_data <- division_scores %>%
  mutate(major_division = factor(major_division, levels = div_order))

cond_cols <- c(
  "forced swim · Acute"     = "#FBC9DE",
  "forced swim · Chronic"   = "#F48FB1",
  "restraint · Acute"       = "#E05780",
  "social defeat · Acute"   = "#C2185B",
  "tail suspension · Acute" = "#8E24AA"
)

p <- ggplot(plot_data, aes(x = Condition, y = div_value, color = Condition)) +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.85) +
  stat_summary(fun = median, geom = "crossbar", width = 0.4, color = "grey20", linewidth = 0.3) +
  facet_wrap(~major_division, scales = "free_y", ncol = 4) +
  scale_color_manual(values = cond_cols) +
  labs(title = "Division-level normalized TRAP intensity by stressor condition",
       subtitle = "14 major anatomical divisions (Allen CCF), ordered by ANOVA p-value (best top-left)",
       x = NULL, y = "Normalized intensity (division median)") +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "bottom", strip.text = element_text(face = "bold", size = 7))

ggsave("~/Desktop/division_level_anova.png", p, width = 12, height = 8, units = "in", dpi = 300)
cat("\nSaved: ~/Desktop/division_level_anova.csv, ~/Desktop/division_level_anova.png\n")


# ---------------------------------------------------------------------------
# effect_sizes_pairwise
# ---------------------------------------------------------------------------
# Standardized effect sizes for the null-model-corrected pairwise comparisons
# (nullmodel_corrected_pairwise_stats.R) and the DT1-vs-TRDT1 sex comparison, since
# FDR-significance is essentially unreachable at n = 3-4 per group (see
# earlier mathematically-proven Wilcoxon p-value floor). Effect sizes let us
# report DIRECTION and MAGNITUDE even where classical significance testing
# is underpowered.
#
#   - Wilcoxon-Mann-Whitney U comparisons -> rank-biserial correlation
#     r_rb = 1 - 2U / (n1*n2). Ranges -1..+1; |r_rb| ~ 0.1/0.3/0.5 are
#     conventionally small/medium/large (Cureton 1956 correspondence to
#     Cohen's d benchmarks).
#   - Welch t-test comparisons -> Cohen's d_s (mean difference / pooled SD
#     using the average-variance definition, appropriate for unequal
#     variances). |d| ~ 0.2/0.5/0.8 = small/medium/large (Cohen 1988).
#
# Uses pre-split lookup lists (not per-row dplyr::filter) for speed across
# the ~6,600 region x condition-pair tests.

library(dplyr)
library(tidyr)

df <- read.csv("~/Desktop/ALL_DT.csv", stringsAsFactors = FALSE)
cond_map <- read.csv("~/Desktop/group_condition_corrected.csv", stringsAsFactors = FALSE)
stats <- read.csv("~/Desktop/nullmodel_corrected_pairwise_stats.csv", stringsAsFactors = FALSE)

long <- df %>%
  filter(region_acronym != "FS") %>%
  select(animal_group, animal_id_original, region_acronym,
         trapcre_intensity, traptta_intensity) %>%
  pivot_longer(cols = c(trapcre_intensity, traptta_intensity),
               names_to = "channel", values_to = "value") %>%
  mutate(channel = ifelse(channel == "trapcre_intensity", "Cre", "tTA")) %>%
  rename(Group = animal_group, Channel = channel) %>%
  left_join(cond_map, by = c("Group", "Channel")) %>%
  filter(!is.na(Condition))

## --- Pre-split into a nested lookup: key = "Condition||region" -> values ---
long$key <- paste(long$Condition, long$region_acronym, sep = "||")
value_lookup <- split(long$value, long$key)

rank_biserial <- function(x, y) {
  n1 <- length(x); n2 <- length(y)
  U <- suppressWarnings(wilcox.test(x, y, exact = FALSE)$statistic)
  1 - (2 * U) / (n1 * n2)
}
cohens_d_s <- function(x, y) {
  n1 <- length(x); n2 <- length(y)
  sd1 <- sd(x); sd2 <- sd(y)
  pooled_sd <- sqrt((sd1^2 + sd2^2) / 2)
  if (pooled_sd == 0) return(0)
  (mean(x) - mean(y)) / pooled_sd
}

n_rows <- nrow(stats)
effect <- numeric(n_rows)
effect_type <- character(n_rows)

cat("Computing effect sizes for", n_rows, "region-comparisons...\n")
for (i in seq_len(n_rows)) {
  row <- stats[i, ]
  x <- value_lookup[[paste(row$condition1, row$region_acronym, sep = "||")]]
  y <- value_lookup[[paste(row$condition2, row$region_acronym, sep = "||")]]
  if (row$test_used == "Welch t-test") {
    effect[i] <- cohens_d_s(x, y)
    effect_type[i] <- "Cohen's d_s"
  } else {
    effect[i] <- rank_biserial(x, y)
    effect_type[i] <- "rank-biserial r"
  }
  if (i %% 1000 == 0) cat("  ", i, "/", n_rows, "\n")
}

stats$effect_size <- effect
stats$effect_type <- effect_type
stats$effect_magnitude <- cut(abs(effect), breaks = c(-Inf, 0.1, 0.3, 0.5, Inf),
                               labels = c("negligible", "small", "medium", "large"))

write.csv(stats, "~/Desktop/nullmodel_corrected_pairwise_stats_with_effectsizes.csv", row.names = FALSE)

cat("\n=== Effect size magnitude distribution (all", n_rows, "tests) ===\n")
print(table(stats$effect_magnitude))

cat("\n=== Top 15 largest effect sizes overall (regardless of FDR) ===\n")
top15 <- stats %>% arrange(desc(abs(effect_size))) %>% head(15) %>%
  select(condition1, condition2, region_acronym, division, test_used, effect_size, effect_type, p_value, p_adj_BH)
print(as.data.frame(top15))

cat("\n=== 'Large' effects (|effect| > 0.5) per comparison ===\n")
large_by_pair <- stats %>% filter(effect_magnitude == "large") %>%
  count(condition1, condition2, name = "n_large_effect") %>% arrange(desc(n_large_effect))
print(as.data.frame(large_by_pair))
write.csv(large_by_pair, "~/Desktop/effect_sizes_large_by_comparison.csv", row.names = FALSE)

## --- Same for the sex comparison (DT1 vs TRDT1), top regions ---
cat("\n\n=== Sex comparison (DT1 vs TRDT1), top regions -- effect sizes ===\n")
top_regions_df <- read.csv("~/Desktop/top_regions_low_deviation.csv", stringsAsFactors = FALSE)
top_regions <- unique(top_regions_df$region_acronym)

sex_long <- df %>% filter(animal_group %in% c("DT1", "TRDT1"), region_acronym %in% top_regions) %>%
  select(animal_group, region_acronym, trapcre_intensity, traptta_intensity) %>%
  pivot_longer(cols = c(trapcre_intensity, traptta_intensity), names_to = "channel", values_to = "value") %>%
  mutate(channel = ifelse(channel == "trapcre_intensity", "Cre", "tTA"))

sex_results <- list()
idx <- 1
for (ch in c("Cre", "tTA")) {
  for (reg in top_regions) {
    x <- sex_long %>% filter(channel == ch, region_acronym == reg, animal_group == "DT1") %>% pull(value)
    y <- sex_long %>% filter(channel == ch, region_acronym == reg, animal_group == "TRDT1") %>% pull(value)
    if (length(x) < 2 || length(y) < 2) next
    r_rb <- rank_biserial(x, y)
    d_s <- cohens_d_s(x, y)
    sex_results[[idx]] <- data.frame(channel = ch, region_acronym = reg,
                                      rank_biserial_r = r_rb, cohens_d = d_s)
    idx <- idx + 1
  }
}
sex_effect <- bind_rows(sex_results) %>% arrange(desc(abs(rank_biserial_r)))
write.csv(sex_effect, "~/Desktop/sex_comparison_top_regions_effectsizes.csv", row.names = FALSE)
print(as.data.frame(head(sex_effect, 15)))

cat("\nSaved:\n  ~/Desktop/nullmodel_corrected_pairwise_stats_with_effectsizes.csv\n")
cat("  ~/Desktop/effect_sizes_large_by_comparison.csv\n")
cat("  ~/Desktop/sex_comparison_top_regions_effectsizes.csv\n")


# ---------------------------------------------------------------------------
# curated_regions_all_channels
# ---------------------------------------------------------------------------
# Extend the DT1-vs-TRDT1 curated-region comparison to ALL 16 group x channel
# combinations (DT1...DT8, Cre + tTA -- "my 8 groups" excluding TRDT1).
#
# Same curated region set as sex_comparison_curated_regions.R: top regions
# (high median, low CV) + most variable regions (high CV, signal-filtered),
# from top_regions_low_deviation.csv / most_variable_regions.csv.
#
# Each animal's per-region value is normalized by that animal's own
# whole-brain median (removes individual global scaling). Per region, per
# channel: Kruskal-Wallis test across the 8 groups (non-parametric one-way
# ANOVA), since group sizes are small (n=2-4) and normality can't be assumed.
# BH-FDR applied across the curated regions.

library(dplyr)
library(tidyr)
library(ggplot2)

df <- read.csv("~/Desktop/ALL_DT.csv", stringsAsFactors = FALSE)
group_order <- paste0("DT", 1:8)

## --- Curated region set: top regions + most variable regions ---
top_regions_df <- read.csv("~/Desktop/top_regions_low_deviation.csv", stringsAsFactors = FALSE)
top_regions <- unique(top_regions_df$region_acronym)
most_var_df <- read.csv("~/Desktop/most_variable_regions.csv", stringsAsFactors = FALSE)
most_variable <- unique(most_var_df$region_acronym)
curated_regions <- unique(c(top_regions, most_variable))
region_set_label <- setNames(
  ifelse(curated_regions %in% top_regions & curated_regions %in% most_variable, "Top & most variable",
         ifelse(curated_regions %in% top_regions, "Top region", "Most variable")),
  curated_regions
)
cat("n curated regions:", length(curated_regions), "\n")

## --- Per-animal data, normalized by that animal's own whole-brain median ---
sub <- df %>% filter(animal_group %in% group_order, region_acronym %in% curated_regions) %>%
  select(animal_group, animal_id_original, region_acronym, division,
         trapcre_intensity, traptta_intensity) %>%
  pivot_longer(cols = c(trapcre_intensity, traptta_intensity),
               names_to = "channel", values_to = "value") %>%
  mutate(channel = ifelse(channel == "trapcre_intensity", "Cre", "tTA"),
         animal_group = factor(animal_group, levels = group_order),
         sample_id = paste(animal_group, channel, sep = "_"))

animal_overall <- df %>% filter(animal_group %in% group_order) %>%
  select(animal_id_original, trapcre_intensity, traptta_intensity) %>%
  pivot_longer(cols = c(trapcre_intensity, traptta_intensity), names_to = "channel", values_to = "value") %>%
  mutate(channel = ifelse(channel == "trapcre_intensity", "Cre", "tTA")) %>%
  group_by(animal_id_original, channel) %>%
  summarise(overall_median = median(value, na.rm = TRUE), .groups = "drop")

sub <- sub %>% left_join(animal_overall, by = c("animal_id_original", "channel")) %>%
  mutate(norm_value = value / overall_median)

cat("n unique group x channel samples:", n_distinct(sub$sample_id), "(expect 16)\n")

## --- Per-region, per-channel Kruskal-Wallis test across the 8 groups ---
run_region_tests <- function(ch) {
  d <- sub %>% filter(channel == ch)
  regions <- unique(d$region_acronym)
  out <- vector("list", length(regions))
  for (i in seq_along(regions)) {
    reg <- regions[i]
    dr <- d %>% filter(region_acronym == reg)
    if (n_distinct(dr$animal_group) < 3 || nrow(dr) < 8) { out[[i]] <- NULL; next }
    kt <- suppressWarnings(kruskal.test(norm_value ~ animal_group, data = dr))
    out[[i]] <- data.frame(region_acronym = reg, channel = ch,
                            n_obs = nrow(dr), chi_sq = kt$statistic, p_value = kt$p.value)
  }
  bind_rows(out)
}

results <- bind_rows(run_region_tests("Cre"), run_region_tests("tTA")) %>%
  group_by(channel) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(region_set = region_set_label[region_acronym]) %>%
  arrange(channel, p_adj_BH, p_value)

write.csv(results, "~/Desktop/curated_regions_all_channels_kruskal.csv", row.names = FALSE)

for (ch in c("Cre", "tTA")) {
  r <- results %>% filter(channel == ch)
  cat(sprintf("\n========== %s (n curated regions tested = %d, across 8 groups) ==========\n", ch, nrow(r)))
  cat(sprintf("n regions raw p < 0.05: %d | BH-adj p < 0.05: %d | BH-adj p < 0.20: %d\n",
              sum(r$p_value < 0.05, na.rm = TRUE), sum(r$p_adj_BH < 0.05, na.rm = TRUE),
              sum(r$p_adj_BH < 0.20, na.rm = TRUE)))
  print(as.data.frame(head(r %>% select(region_acronym, region_set, chi_sq, p_value, p_adj_BH), 15)))
}

## --- Plot: all 16 samples per curated region (grid of small multiples) ---
group_cols <- setNames(
  c("#FBC9DE", "#F6A6C8", "#F48FB1", "#E96D9B", "#E05780", "#C2185B", "#8E24AA", "#4A148C"),
  group_order
)

make_plot <- function(ch, top_n = 20) {
  ord <- results %>% filter(channel == ch) %>% slice_head(n = top_n) %>% pull(region_acronym)
  d <- sub %>% filter(channel == ch, region_acronym %in% ord) %>%
    mutate(region_acronym = factor(region_acronym, levels = ord))

  ggplot(d, aes(x = animal_group, y = norm_value, color = animal_group)) +
    geom_jitter(width = 0.15, size = 1.6, alpha = 0.85) +
    stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "grey20", linewidth = 0.25) +
    facet_wrap(~region_acronym, scales = "free_y", ncol = 4) +
    scale_color_manual(values = group_cols, guide = guide_legend(nrow = 1)) +
    labs(title = sprintf("%s: top %d curated regions across all 8 groups (16 group x channel samples total)", ch, top_n),
         subtitle = "Normalized per animal. Regions ordered by Kruskal-Wallis BH-adjusted p-value (best first)",
         x = NULL, y = "Normalized intensity", color = "Group") +
    theme_minimal(base_size = 8) +
    theme(strip.text = element_text(face = "bold", size = 7.5),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "bottom")
}

p_cre <- make_plot("Cre")
p_tta <- make_plot("tTA")

ggsave("~/Desktop/curated_regions_all_channels_Cre.png", p_cre, width = 12, height = 11, units = "in", dpi = 300)
ggsave("~/Desktop/curated_regions_all_channels_tTA.png", p_tta, width = 12, height = 11, units = "in", dpi = 300)


# ---------------------------------------------------------------------------
# curated_regions_all_channels_raw
# ---------------------------------------------------------------------------
# RAW version of curated_regions_all_channels.R -- no global/batch
# normalization. Curated regions (top + most variable), across all 16
# group x channel samples (DT1...DT8, Cre + tTA). Kruskal-Wallis per region
# per channel across the 8 groups, on raw intensity values. BH-FDR across
# the curated regions.

library(dplyr)
library(tidyr)
library(ggplot2)

df <- read.csv("~/Desktop/ALL_DT.csv", stringsAsFactors = FALSE)
group_order <- paste0("DT", 1:8)

top_regions_df <- read.csv("~/Desktop/top_regions_low_deviation.csv", stringsAsFactors = FALSE)
top_regions <- unique(top_regions_df$region_acronym)
most_var_df <- read.csv("~/Desktop/most_variable_regions.csv", stringsAsFactors = FALSE)
most_variable <- unique(most_var_df$region_acronym)
curated_regions <- unique(c(top_regions, most_variable))
region_set_label <- setNames(
  ifelse(curated_regions %in% top_regions & curated_regions %in% most_variable, "Top & most variable",
         ifelse(curated_regions %in% top_regions, "Top region", "Most variable")),
  curated_regions
)
cat("n curated regions:", length(curated_regions), "\n")

sub <- df %>% filter(animal_group %in% group_order, region_acronym %in% curated_regions) %>%
  select(animal_group, animal_id_original, region_acronym, division,
         trapcre_intensity, traptta_intensity) %>%
  pivot_longer(cols = c(trapcre_intensity, traptta_intensity),
               names_to = "channel", values_to = "value") %>%
  mutate(channel = ifelse(channel == "trapcre_intensity", "Cre", "tTA"),
         animal_group = factor(animal_group, levels = group_order),
         sample_id = paste(animal_group, channel, sep = "_"))

cat("n unique group x channel samples:", n_distinct(sub$sample_id), "(expect 16)\n")

run_region_tests <- function(ch) {
  d <- sub %>% filter(channel == ch)
  regions <- unique(d$region_acronym)
  out <- vector("list", length(regions))
  for (i in seq_along(regions)) {
    reg <- regions[i]
    dr <- d %>% filter(region_acronym == reg)
    if (n_distinct(dr$animal_group) < 3 || nrow(dr) < 8) { out[[i]] <- NULL; next }
    kt <- suppressWarnings(kruskal.test(value ~ animal_group, data = dr))
    out[[i]] <- data.frame(region_acronym = reg, channel = ch,
                            n_obs = nrow(dr), chi_sq = kt$statistic, p_value = kt$p.value)
  }
  bind_rows(out)
}

results <- bind_rows(run_region_tests("Cre"), run_region_tests("tTA")) %>%
  group_by(channel) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(region_set = region_set_label[region_acronym]) %>%
  arrange(channel, p_adj_BH, p_value)

write.csv(results, "~/Desktop/curated_regions_all_channels_kruskal_raw.csv", row.names = FALSE)

for (ch in c("Cre", "tTA")) {
  r <- results %>% filter(channel == ch)
  cat(sprintf("\n========== %s RAW (n curated regions tested = %d, across 8 groups) ==========\n", ch, nrow(r)))
  cat(sprintf("n regions raw p < 0.05: %d | BH-adj p < 0.05: %d | BH-adj p < 0.20: %d\n",
              sum(r$p_value < 0.05, na.rm = TRUE), sum(r$p_adj_BH < 0.05, na.rm = TRUE),
              sum(r$p_adj_BH < 0.20, na.rm = TRUE)))
  print(as.data.frame(head(r %>% select(region_acronym, region_set, chi_sq, p_value, p_adj_BH), 15)))
}

## --- Plot ---
group_cols <- setNames(
  c("#FBC9DE", "#F6A6C8", "#F48FB1", "#E96D9B", "#E05780", "#C2185B", "#8E24AA", "#4A148C"),
  group_order
)

make_plot <- function(ch, top_n = 20) {
  ord <- results %>% filter(channel == ch) %>% slice_head(n = top_n) %>% pull(region_acronym)
  d <- sub %>% filter(channel == ch, region_acronym %in% ord) %>%
    mutate(region_acronym = factor(region_acronym, levels = ord))

  ggplot(d, aes(x = animal_group, y = value, color = animal_group)) +
    geom_jitter(width = 0.15, size = 1.6, alpha = 0.85) +
    stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "grey20", linewidth = 0.25) +
    facet_wrap(~region_acronym, scales = "free_y", ncol = 4) +
    scale_color_manual(values = group_cols, guide = guide_legend(nrow = 1)) +
    labs(title = sprintf("%s (RAW, no normalization): top %d curated regions across all 8 groups", ch, top_n),
         subtitle = "Ordered by Kruskal-Wallis BH-adjusted p-value (best first)",
         x = NULL, y = "Raw intensity", color = "Group") +
    theme_minimal(base_size = 8) +
    theme(strip.text = element_text(face = "bold", size = 7.5),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "bottom")
}

p_cre <- make_plot("Cre")
p_tta <- make_plot("tTA")

ggsave("~/Desktop/curated_regions_all_channels_raw_Cre.png", p_cre, width = 12, height = 11, units = "in", dpi = 300)
ggsave("~/Desktop/curated_regions_all_channels_raw_tTA.png", p_tta, width = 12, height = 11, units = "in", dpi = 300)


# ---------------------------------------------------------------------------
# curated_regions_pairwise_wilcoxon
# ---------------------------------------------------------------------------
# Post-hoc pairwise comparisons following the Kruskal-Wallis screen: for each
# curated region x channel, pairwise Wilcoxon rank-sum tests between all
# pairs of the 8 groups (28 pairs), with Benjamini-Hochberg FDR correction
# applied GLOBALLY across every region x channel x pair combination (the
# most conservative/correct option, consistent with the rest of this
# analysis -- not just corrected within each region separately).
#
# Same curated region set + per-animal normalization as
# curated_regions_all_channels.R.

library(dplyr)
library(tidyr)
library(ggplot2)

df <- read.csv("~/Desktop/ALL_DT.csv", stringsAsFactors = FALSE)
group_order <- paste0("DT", 1:8)

top_regions_df <- read.csv("~/Desktop/top_regions_low_deviation.csv", stringsAsFactors = FALSE)
top_regions <- unique(top_regions_df$region_acronym)
most_var_df <- read.csv("~/Desktop/most_variable_regions.csv", stringsAsFactors = FALSE)
most_variable <- unique(most_var_df$region_acronym)
curated_regions <- unique(c(top_regions, most_variable))
region_set_label <- setNames(
  ifelse(curated_regions %in% top_regions & curated_regions %in% most_variable, "Top & most variable",
         ifelse(curated_regions %in% top_regions, "Top region", "Most variable")),
  curated_regions
)

sub <- df %>% filter(animal_group %in% group_order, region_acronym %in% curated_regions) %>%
  select(animal_group, animal_id_original, region_acronym, division,
         trapcre_intensity, traptta_intensity) %>%
  pivot_longer(cols = c(trapcre_intensity, traptta_intensity),
               names_to = "channel", values_to = "value") %>%
  mutate(channel = ifelse(channel == "trapcre_intensity", "Cre", "tTA"),
         animal_group = factor(animal_group, levels = group_order))

animal_overall <- df %>% filter(animal_group %in% group_order) %>%
  select(animal_id_original, trapcre_intensity, traptta_intensity) %>%
  pivot_longer(cols = c(trapcre_intensity, traptta_intensity), names_to = "channel", values_to = "value") %>%
  mutate(channel = ifelse(channel == "trapcre_intensity", "Cre", "tTA")) %>%
  group_by(animal_id_original, channel) %>%
  summarise(overall_median = median(value, na.rm = TRUE), .groups = "drop")

sub <- sub %>% left_join(animal_overall, by = c("animal_id_original", "channel")) %>%
  mutate(norm_value = value / overall_median)

## --- All pairwise Wilcoxon rank-sum tests, per region x channel ---
group_pairs <- t(combn(group_order, 2))
all_results <- list()
idx <- 1

for (reg in curated_regions) {
  for (ch in c("Cre", "tTA")) {
    d <- sub %>% filter(region_acronym == reg, channel == ch)
    if (n_distinct(d$animal_group) < 2) next
    for (p in seq_len(nrow(group_pairs))) {
      g1 <- group_pairs[p, 1]; g2 <- group_pairs[p, 2]
      v1 <- d %>% filter(animal_group == g1) %>% pull(norm_value)
      v2 <- d %>% filter(animal_group == g2) %>% pull(norm_value)
      if (length(v1) < 2 || length(v2) < 2) next
      wt <- suppressWarnings(wilcox.test(v1, v2))
      all_results[[idx]] <- data.frame(
        region_acronym = reg, channel = ch, group1 = g1, group2 = g2,
        median1 = median(v1, na.rm = TRUE), median2 = median(v2, na.rm = TRUE),
        p_value = wt$p.value
      )
      idx <- idx + 1
    }
  }
}

results <- bind_rows(all_results) %>%
  mutate(p_adj_BH_global = p.adjust(p_value, method = "BH"),
         region_set = region_set_label[region_acronym]) %>%
  arrange(p_adj_BH_global, p_value)

write.csv(results, "~/Desktop/curated_regions_pairwise_wilcoxon.csv", row.names = FALSE)

cat("Total pairwise tests run:", nrow(results), "\n")
cat("n raw p < 0.05:", sum(results$p_value < 0.05, na.rm = TRUE), "\n")
cat("n BH-adjusted (global) p < 0.05:", sum(results$p_adj_BH_global < 0.05, na.rm = TRUE), "\n")
cat("n BH-adjusted (global) p < 0.10:", sum(results$p_adj_BH_global < 0.10, na.rm = TRUE), "\n")
cat("n BH-adjusted (global) p < 0.20:", sum(results$p_adj_BH_global < 0.20, na.rm = TRUE), "\n")
cat("Smallest possible p-value at n>=2 per group (mathematical floor for the smallest groups): check min p_value below\n")
cat("Minimum raw p-value observed:", min(results$p_value, na.rm = TRUE), "\n")

cat("\n=== Top 20 pairwise comparisons by global BH-adjusted p-value ===\n")
print(as.data.frame(head(results %>%
  select(region_acronym, region_set, channel, group1, group2, median1, median2, p_value, p_adj_BH_global), 20)))


# ---------------------------------------------------------------------------
# curated_regions_pairwise_wilcoxon_raw
# ---------------------------------------------------------------------------
# RAW version of curated_regions_pairwise_wilcoxon.R -- no global/batch
# normalization. Pairwise Wilcoxon rank-sum tests between all 8 groups, per
# curated region x channel, on raw intensity values. BH-FDR applied globally
# across all region x channel x pair combinations.

library(dplyr)
library(tidyr)

df <- read.csv("~/Desktop/ALL_DT.csv", stringsAsFactors = FALSE)
group_order <- paste0("DT", 1:8)

top_regions_df <- read.csv("~/Desktop/top_regions_low_deviation.csv", stringsAsFactors = FALSE)
top_regions <- unique(top_regions_df$region_acronym)
most_var_df <- read.csv("~/Desktop/most_variable_regions.csv", stringsAsFactors = FALSE)
most_variable <- unique(most_var_df$region_acronym)
curated_regions <- unique(c(top_regions, most_variable))
region_set_label <- setNames(
  ifelse(curated_regions %in% top_regions & curated_regions %in% most_variable, "Top & most variable",
         ifelse(curated_regions %in% top_regions, "Top region", "Most variable")),
  curated_regions
)

sub <- df %>% filter(animal_group %in% group_order, region_acronym %in% curated_regions) %>%
  select(animal_group, animal_id_original, region_acronym, division,
         trapcre_intensity, traptta_intensity) %>%
  pivot_longer(cols = c(trapcre_intensity, traptta_intensity),
               names_to = "channel", values_to = "value") %>%
  mutate(channel = ifelse(channel == "trapcre_intensity", "Cre", "tTA"),
         animal_group = factor(animal_group, levels = group_order))

group_pairs <- t(combn(group_order, 2))
all_results <- list()
idx <- 1

for (reg in curated_regions) {
  for (ch in c("Cre", "tTA")) {
    d <- sub %>% filter(region_acronym == reg, channel == ch)
    if (n_distinct(d$animal_group) < 2) next
    for (p in seq_len(nrow(group_pairs))) {
      g1 <- group_pairs[p, 1]; g2 <- group_pairs[p, 2]
      v1 <- d %>% filter(animal_group == g1) %>% pull(value)
      v2 <- d %>% filter(animal_group == g2) %>% pull(value)
      if (length(v1) < 2 || length(v2) < 2) next
      wt <- suppressWarnings(wilcox.test(v1, v2))
      all_results[[idx]] <- data.frame(
        region_acronym = reg, channel = ch, group1 = g1, group2 = g2,
        median1 = median(v1, na.rm = TRUE), median2 = median(v2, na.rm = TRUE),
        p_value = wt$p.value
      )
      idx <- idx + 1
    }
  }
}

results <- bind_rows(all_results) %>%
  mutate(p_adj_BH_global = p.adjust(p_value, method = "BH"),
         region_set = region_set_label[region_acronym]) %>%
  arrange(p_adj_BH_global, p_value)

write.csv(results, "~/Desktop/curated_regions_pairwise_wilcoxon_raw.csv", row.names = FALSE)

cat("Total pairwise tests run:", nrow(results), "\n")
cat("n raw p < 0.05:", sum(results$p_value < 0.05, na.rm = TRUE), "\n")
cat("n BH-adjusted (global) p < 0.05:", sum(results$p_adj_BH_global < 0.05, na.rm = TRUE), "\n")
cat("n BH-adjusted (global) p < 0.10:", sum(results$p_adj_BH_global < 0.10, na.rm = TRUE), "\n")
cat("n BH-adjusted (global) p < 0.20:", sum(results$p_adj_BH_global < 0.20, na.rm = TRUE), "\n")
cat("Minimum raw p-value observed:", min(results$p_value, na.rm = TRUE), "\n")

cat("\n=== Top 20 pairwise comparisons by global BH-adjusted p-value ===\n")
print(as.data.frame(head(results %>%
  select(region_acronym, region_set, channel, group1, group2, median1, median2, p_value, p_adj_BH_global), 20)))


# ---------------------------------------------------------------------------
# nullmodel_corrected_pairwise_stats
# ---------------------------------------------------------------------------
# Pairwise stressor comparison, following the statistical pipeline of
# the reference implementation (2026) "Stressor- and Time-Specific Modulation of Brain Activity and
# Connectivity in Mice":
#
#   1. Shapiro-Wilk normality test per region, per condition. A region is
#      "normal" only if BOTH conditions being compared pass (p >= 0.05).
#   2. Welch's two-sample t-test for normal regions; two-sided Wilcoxon-
#      Mann-Whitney U test (continuity-corrected) otherwise.
#   3. Benjamini-Hochberg FDR correction (threshold FDR < 0.05), applied
#      within each pairwise comparison across all tested regions.
#   4. Only regions with n >= 3 per condition are tested.
#   5. Effect sizes: mean difference with 95% BCa bootstrap CI (DABEST /
#      Gardner-Altman), computed for the FDR-significant regions (or, if
#      none survive FDR, for the top candidates by raw p-value -- clearly
#      labeled as exploratory).
#
# Adapted to this dataset: there is no untreated control group, so "Stressor
# Specificity" (the Level 2) becomes ALL PAIRWISE comparisons among the
# 5 real stressor conditions (forced swim acute/chronic, restraint acute,
# social defeat acute, tail suspension acute). Each condition pools the
# individual ANIMALS (not group medians) whose Cre or tTA channel carries
# that condition, per group_condition_corrected.csv. 'FS' excluded
# (injection region, trivially 1.0).

library(dplyr)
library(tidyr)
library(ggplot2)
library(dabestr)

df <- read.csv("~/Desktop/ALL_DT.csv", stringsAsFactors = FALSE)
cond_map <- read.csv("~/Desktop/group_condition_corrected.csv", stringsAsFactors = FALSE)

## --- Animal-level long format: one row per animal x region x channel,
## each with its OWN real condition (not the group-level average) ---
long <- df %>%
  filter(region_acronym != "FS") %>%
  select(animal_group, animal_id_original, region_acronym, division,
         trapcre_intensity, traptta_intensity) %>%
  pivot_longer(cols = c(trapcre_intensity, traptta_intensity),
               names_to = "channel", values_to = "value") %>%
  mutate(channel = ifelse(channel == "trapcre_intensity", "Cre", "tTA")) %>%
  rename(Group = animal_group, Channel = channel) %>%
  left_join(cond_map, by = c("Group", "Channel")) %>%
  filter(!is.na(Condition))

conditions <- sort(unique(long$Condition))
cat("Conditions:\n"); print(conditions)
cat("\nn animal-samples per condition:\n")
print(long %>% distinct(animal_id_original, Channel, Condition) %>% count(Condition))

division_lookup <- long %>% distinct(region_acronym, division)
cond_pairs <- t(combn(conditions, 2))

## --- Per-region test: Shapiro-Wilk -> Welch t-test or Wilcoxon ---
test_region <- function(x, y) {
  n1 <- length(x); n2 <- length(y)
  if (n1 < 3 || n2 < 3) return(NULL)

  sw1_p <- tryCatch(shapiro.test(x)$p.value, error = function(e) NA)
  sw2_p <- tryCatch(shapiro.test(y)$p.value, error = function(e) NA)
  both_normal <- !is.na(sw1_p) && !is.na(sw2_p) && sw1_p >= 0.05 && sw2_p >= 0.05

  if (both_normal) {
    tt <- t.test(x, y, var.equal = FALSE)   # Welch
    test_used <- "Welch t-test"
    p_val <- tt$p.value
  } else {
    wt <- suppressWarnings(wilcox.test(x, y, correct = TRUE))
    test_used <- "Wilcoxon-Mann-Whitney U"
    p_val <- wt$p.value
  }
  data.frame(n1 = n1, n2 = n2, mean1 = mean(x), mean2 = mean(y),
             median1 = median(x), median2 = median(y),
             shapiro_p1 = sw1_p, shapiro_p2 = sw2_p,
             test_used = test_used, p_value = p_val)
}

all_results <- list()
idx <- 1

for (p in seq_len(nrow(cond_pairs))) {
  c1 <- cond_pairs[p, 1]; c2 <- cond_pairs[p, 2]
  d <- long %>% filter(Condition %in% c(c1, c2))
  regions <- unique(d$region_acronym)

  for (reg in regions) {
    dr <- d %>% filter(region_acronym == reg)
    x <- dr %>% filter(Condition == c1) %>% pull(value)
    y <- dr %>% filter(Condition == c2) %>% pull(value)
    res <- test_region(x, y)
    if (is.null(res)) next
    all_results[[idx]] <- cbind(condition1 = c1, condition2 = c2, region_acronym = reg, res)
    idx <- idx + 1
  }
  cat(sprintf("Comparison done: %s vs %s (%d regions tested)\n", c1, c2, length(regions)))
}

results <- bind_rows(all_results) %>%
  group_by(condition1, condition2) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  left_join(division_lookup, by = "region_acronym") %>%
  arrange(condition1, condition2, p_adj_BH, p_value)

write.csv(results, "~/Desktop/nullmodel_corrected_pairwise_stats.csv", row.names = FALSE)

cat("\n\n========== SUMMARY ACROSS ALL 10 PAIRWISE COMPARISONS ==========\n")
summary_tbl <- results %>% group_by(condition1, condition2) %>%
  summarise(n_regions_tested = n(),
            n_welch = sum(test_used == "Welch t-test"),
            n_wilcoxon = sum(test_used == "Wilcoxon-Mann-Whitney U"),
            n_raw_p_lt_05 = sum(p_value < 0.05),
            n_FDR_lt_05 = sum(p_adj_BH < 0.05),
            .groups = "drop")
print(as.data.frame(summary_tbl))

write.csv(summary_tbl, "~/Desktop/nullmodel_corrected_pairwise_stats_summary.csv", row.names = FALSE)


# ---------------------------------------------------------------------------
# condition_specific_regions
# ---------------------------------------------------------------------------
# For EVERY region: are the samples sharing the same real stressor condition
# similar to each other, AND clearly different from the other conditions?
# I.e. does this region show a condition-specific signature, on top of the
# dominant "same animal/group" similarity found previously?
#
# Method: per region, one-way ANOVA (log10 intensity ~ condition) across the
# 18 (group,channel) samples. Effect size = eta-squared (share of variance
# explained by condition). Ranked list = candidate condition-specific
# regions. Given n = 2-6 samples per condition, p-values have low power --
# treat eta-squared as a screening signal, not a confirmed finding.
#
# Condition assignment: group_condition_corrected.csv (verified constant per
# group via animalid.xlsx).

library(dplyr)
library(tidyr)
library(ggplot2)

df <- read.csv("~/Desktop/ALL_DT.csv", stringsAsFactors = FALSE)
cond_map <- read.csv("~/Desktop/group_condition_corrected.csv", stringsAsFactors = FALSE)

group_med <- df %>%
  filter(region_acronym != "FS") %>%
  group_by(region_acronym, division, animal_group) %>%
  summarise(Cre = median(trapcre_intensity, na.rm = TRUE),
            tTA = median(traptta_intensity, na.rm = TRUE), .groups = "drop")

long <- group_med %>%
  pivot_longer(cols = c(Cre, tTA), names_to = "Channel", values_to = "value") %>%
  rename(Group = animal_group) %>%
  left_join(cond_map, by = c("Group", "Channel")) %>%
  mutate(log_value = log10(value + 0.001))

division_lookup <- long %>% distinct(region_acronym, division)

## --- Per-region one-way ANOVA: log_value ~ condition ---
regions <- unique(long$region_acronym)
results <- vector("list", length(regions))

for (i in seq_along(regions)) {
  reg <- regions[i]
  d <- long %>% filter(region_acronym == reg)
  if (n_distinct(d$Condition) < 2 || nrow(d) < 5) { results[[i]] <- NULL; next }
  fit <- aov(log_value ~ Condition, data = d)
  ss <- summary(fit)[[1]]
  ss_cond <- ss["Condition", "Sum Sq"]
  ss_resid <- ss["Residuals", "Sum Sq"]
  eta2 <- ss_cond / (ss_cond + ss_resid)
  pval <- ss["Condition", "Pr(>F)"]
  results[[i]] <- data.frame(region_acronym = reg, eta_sq = eta2, p_value = pval)
}

anova_results <- bind_rows(results) %>%
  left_join(division_lookup, by = "region_acronym") %>%
  arrange(desc(eta_sq))

write.csv(anova_results, "~/Desktop/condition_specific_regions_anova.csv", row.names = FALSE)

cat("n regions tested:", nrow(anova_results), "\n")
cat("\n=== Top 20 most condition-discriminating regions (highest eta-squared) ===\n")
print(as.data.frame(head(anova_results, 20)))

cat(sprintf("\nRegions with eta_sq > 0.5: %d\n", sum(anova_results$eta_sq > 0.5, na.rm = TRUE)))
cat(sprintf("Regions with p < 0.05 (uncorrected, low power given n): %d of %d\n",
            sum(anova_results$p_value < 0.05, na.rm = TRUE), nrow(anova_results)))

## --- Plot: top 12 regions, values by condition ---
top12 <- head(anova_results$region_acronym, 12)
plot_data <- long %>% filter(region_acronym %in% top12) %>%
  mutate(region_acronym = factor(region_acronym, levels = top12))

cond_cols <- c(
  "forced swim · Acute"     = "#FBC9DE",
  "forced swim · Chronic"   = "#F48FB1",
  "restraint · Acute"       = "#E05780",
  "social defeat · Acute"   = "#C2185B",
  "tail suspension · Acute" = "#8E24AA"
)

p <- ggplot(plot_data, aes(x = Condition, y = value, color = Condition)) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.85) +
  stat_summary(fun = median, geom = "crossbar", width = 0.4, color = "grey20", linewidth = 0.3) +
  facet_wrap(~region_acronym, scales = "free_y", ncol = 3) +
  scale_color_manual(values = cond_cols) +
  labs(title = "Top 12 most condition-discriminating regions (highest eta-squared)",
       subtitle = "Each point = one group/channel sample, colored by real stressor condition",
       x = NULL, y = "Intensity (group median)") +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "bottom", strip.text = element_text(face = "bold"))

ggsave("~/Desktop/condition_specific_regions_top12.png", p, width = 11, height = 9, units = "in", dpi = 300)
ggsave("~/Desktop/condition_specific_regions_top12.pdf", p, width = 11, height = 9, units = "in")


# ---------------------------------------------------------------------------
# condition_specific_regions_normalized
# ---------------------------------------------------------------------------
# Same condition-discrimination screen as condition_specific_regions.R, but
# each sample is first normalized by its OWN global median (across all
# regions) before testing. This removes the global group/batch-level scaling
# effect confirmed in condition_specific_regions.R (top "condition-specific"
# regions there correlated r = 0.68-0.95 with each sample's overall
# intensity) -- what's left, if anything, is a genuine REGION-SPECIFIC
# condition effect, not just "this animal group had more signal everywhere."
#
# Condition assignment: group_condition_corrected.csv (verified constant per
# group via animalid.xlsx).

library(dplyr)
library(tidyr)
library(ggplot2)

df <- read.csv("~/Desktop/ALL_DT.csv", stringsAsFactors = FALSE)
cond_map <- read.csv("~/Desktop/group_condition_corrected.csv", stringsAsFactors = FALSE)

group_med <- df %>%
  filter(region_acronym != "FS") %>%
  group_by(region_acronym, division, animal_group) %>%
  summarise(Cre = median(trapcre_intensity, na.rm = TRUE),
            tTA = median(traptta_intensity, na.rm = TRUE), .groups = "drop")

long <- group_med %>%
  pivot_longer(cols = c(Cre, tTA), names_to = "Channel", values_to = "value") %>%
  rename(Group = animal_group) %>%
  left_join(cond_map, by = c("Group", "Channel"))

## --- Normalize: divide each region's value by that sample's own global
## median (across all regions), so every (Group, Channel) sample is put on
## the same overall scale before comparing conditions ---
overall <- long %>% group_by(Group, Channel) %>%
  summarise(overall_median = median(value, na.rm = TRUE), .groups = "drop")

long <- long %>%
  left_join(overall, by = c("Group", "Channel")) %>%
  mutate(norm_value = value / overall_median,
         log_norm_value = log10(norm_value + 0.001))

division_lookup <- long %>% distinct(region_acronym, division)

## --- Per-region one-way ANOVA on NORMALIZED values: log_norm_value ~ condition ---
regions <- unique(long$region_acronym)
results <- vector("list", length(regions))

for (i in seq_along(regions)) {
  reg <- regions[i]
  d <- long %>% filter(region_acronym == reg)
  if (n_distinct(d$Condition) < 2 || nrow(d) < 5) { results[[i]] <- NULL; next }
  fit <- aov(log_norm_value ~ Condition, data = d)
  ss <- summary(fit)[[1]]
  ss_cond <- ss["Condition", "Sum Sq"]
  ss_resid <- ss["Residuals", "Sum Sq"]
  eta2 <- ss_cond / (ss_cond + ss_resid)
  pval <- ss["Condition", "Pr(>F)"]
  results[[i]] <- data.frame(region_acronym = reg, eta_sq = eta2, p_value = pval)
}

anova_results <- bind_rows(results) %>%
  left_join(division_lookup, by = "region_acronym") %>%
  arrange(desc(eta_sq))

## --- Sanity check: does eta_sq still track the (now removed) global scale? ---
anova_results$r_with_global <- sapply(anova_results$region_acronym, function(reg) {
  d <- long %>% filter(region_acronym == reg)
  suppressWarnings(cor(log10(d$value + 0.001), log10(d$overall_median + 0.001)))
})

write.csv(anova_results, "~/Desktop/condition_specific_regions_normalized_anova.csv", row.names = FALSE)

cat("n regions tested:", nrow(anova_results), "\n")
cat("\n=== Top 20 most condition-discriminating regions AFTER normalization ===\n")
print(as.data.frame(head(anova_results, 20)))

cat(sprintf("\nRegions with eta_sq > 0.5: %d\n", sum(anova_results$eta_sq > 0.5, na.rm = TRUE)))
cat(sprintf("Regions with p < 0.05 (uncorrected, low power given n): %d of %d\n",
            sum(anova_results$p_value < 0.05, na.rm = TRUE), nrow(anova_results)))
cat(sprintf("\nCorrelation between eta_sq and |r_with_global_scale| (should now be much lower than the raw 0.61): %.3f\n",
            cor(anova_results$eta_sq, abs(anova_results$r_with_global), use = "complete.obs")))

## --- Plot: top 12 regions, NORMALIZED values by condition ---
top12 <- head(anova_results$region_acronym, 12)
plot_data <- long %>% filter(region_acronym %in% top12) %>%
  mutate(region_acronym = factor(region_acronym, levels = top12))

cond_cols <- c(
  "forced swim · Acute"     = "#FBC9DE",
  "forced swim · Chronic"   = "#F48FB1",
  "restraint · Acute"       = "#E05780",
  "social defeat · Acute"   = "#C2185B",
  "tail suspension · Acute" = "#8E24AA"
)

p <- ggplot(plot_data, aes(x = Condition, y = norm_value, color = Condition)) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.85) +
  stat_summary(fun = median, geom = "crossbar", width = 0.4, color = "grey20", linewidth = 0.3) +
  facet_wrap(~region_acronym, scales = "free_y", ncol = 3) +
  scale_color_manual(values = cond_cols) +
  labs(title = "Top 12 most condition-discriminating regions AFTER removing global scale",
       subtitle = "Value = region intensity / sample's own overall median (across all regions)",
       x = NULL, y = "Normalized intensity") +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "bottom", strip.text = element_text(face = "bold"))

ggsave("~/Desktop/condition_specific_regions_normalized_top12.png", p, width = 11, height = 9, units = "in", dpi = 300)
ggsave("~/Desktop/condition_specific_regions_normalized_top12.pdf", p, width = 11, height = 9, units = "in")


# ---------------------------------------------------------------------------
# test_condition_pattern_similarity
# ---------------------------------------------------------------------------
# Test: do (group, channel) samples that share the same real stressor
# condition show MORE similar expression patterns (across ALL regions) than
# samples from different conditions?
#
# Condition assignment: group_condition_corrected.csv (one Cre + one tTA
# condition per GROUP, confirmed constant across all animals within a group
# via animalid.xlsx -- ALL_DT.csv's own per-animal condition columns contain
# a data error and are not used here).
#
# Method: build a region x (group,channel) matrix across ALL regions (not
# just the pre-selected "top" ones -- 'FS' still excluded, injection region,
# trivially 1.0), compute all pairwise Pearson correlations between the 18
# (group,channel) columns, then compare "same condition" pairwise
# correlations against "different condition" pairs (Wilcoxon rank-sum test).

library(dplyr)
library(tidyr)
library(ggplot2)

df <- read.csv("~/Desktop/ALL_DT.csv", stringsAsFactors = FALSE)
cond_map <- read.csv("~/Desktop/group_condition_corrected.csv", stringsAsFactors = FALSE)

## --- Median per region x group x channel, across ALL regions ---
group_med <- df %>%
  filter(region_acronym != "FS") %>%
  group_by(region_acronym, animal_group) %>%
  summarise(Cre = median(trapcre_intensity, na.rm = TRUE),
            tTA = median(traptta_intensity, na.rm = TRUE), .groups = "drop")
cat("n regions used:", n_distinct(group_med$region_acronym), "\n")

long <- group_med %>%
  pivot_longer(cols = c(Cre, tTA), names_to = "channel", values_to = "value") %>%
  mutate(sample_id = paste(animal_group, channel, sep = "_"),
         log_value = log10(value + 0.001))

mat <- long %>% select(region_acronym, sample_id, log_value) %>%
  pivot_wider(names_from = sample_id, values_from = log_value) %>%
  tibble::column_to_rownames("region_acronym") %>%
  as.matrix()

cat("n samples (groups x channels):", ncol(mat), "\n")

sample_cond <- cond_map %>% transmute(sample_id = paste(Group, Channel, sep = "_"), condition = Condition)
stopifnot(!any(duplicated(sample_cond$sample_id)))

## --- Pairwise Pearson correlation between all 18 samples (on top regions only) ---
cor_mat <- cor(mat, use = "pairwise.complete.obs", method = "pearson")

samples <- colnames(mat)
pairs <- t(combn(samples, 2))
cond_lookup <- setNames(sample_cond$condition, sample_cond$sample_id)

pair_df <- data.frame(sample_a = pairs[, 1], sample_b = pairs[, 2], stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(r = cor_mat[sample_a, sample_b],
         condition_a = cond_lookup[[sample_a]],
         condition_b = cond_lookup[[sample_b]]) %>%
  ungroup() %>%
  mutate(same_condition = condition_a == condition_b)

cat("\nTotal pairs:", nrow(pair_df), " (expect choose(", length(samples), ",2) =", choose(length(samples), 2), ")\n")

## --- Compare same-condition vs different-condition pairwise correlations ---
test_result <- wilcox.test(r ~ same_condition, data = pair_df)
summary_stats <- pair_df %>% group_by(same_condition) %>%
  summarise(n_pairs = n(), mean_r = mean(r), median_r = median(r), sd_r = sd(r))

cat("\n=== Pairwise correlation summary ===\n")
print(as.data.frame(summary_stats))
cat("\nWilcoxon rank-sum test (same vs. different condition):\n")
print(test_result)

write.csv(pair_df, "~/Desktop/condition_pattern_similarity_pairs.csv", row.names = FALSE)

## --- Plot 1: boxplot, same vs different condition ---
p1 <- ggplot(pair_df, aes(x = same_condition, y = r, fill = same_condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.6) +
  scale_x_discrete(labels = c("FALSE" = "Different\ncondition", "TRUE" = "Same\ncondition")) +
  scale_fill_manual(values = c("FALSE" = "#F48FB1", "TRUE" = "#4A148C"), guide = "none") +
  labs(x = NULL, y = "Pairwise Pearson r (log intensity, all regions)",
       title = "Do same-condition samples show more similar expression patterns?",
       subtitle = sprintf("Group-level, n = %d samples | Wilcoxon p = %.4f | mean r same = %.3f vs. different = %.3f",
                           length(samples), test_result$p.value,
                           summary_stats$mean_r[summary_stats$same_condition],
                           summary_stats$mean_r[!summary_stats$same_condition])) +
  theme_minimal(base_size = 10.5)

ggsave("~/Desktop/condition_pattern_similarity_boxplot.png", p1, width = 6, height = 5, units = "in", dpi = 300)
ggsave("~/Desktop/condition_pattern_similarity_boxplot.pdf", p1, width = 6, height = 5, units = "in")

## --- Plot 2: correlation heatmap, samples ordered/annotated by condition ---
library(ComplexHeatmap)
library(circlize)

ord <- order(sample_cond$condition[match(samples, sample_cond$sample_id)])
cond_ordered <- sample_cond$condition[match(samples, sample_cond$sample_id)][ord]

ht <- Heatmap(
  cor_mat[ord, ord],
  name = "Pearson r",
  col = colorRamp2(c(min(cor_mat), 1), c("#FDF0F5", "#4A148C")),
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_split = cond_ordered, column_split = cond_ordered,
  row_title_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 8),
  row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 7),
  column_names_rot = 45,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", cor_mat[ord, ord][i, j]), x, y, gp = gpar(fontsize = 6))
  }
)

pdf("~/Desktop/condition_pattern_similarity_heatmap.pdf", width = 8, height = 7.5)
draw(ht)
dev.off()
png("~/Desktop/condition_pattern_similarity_heatmap.png", width = 8, height = 7.5, units = "in", res = 300)
draw(ht)
dev.off()


# ---------------------------------------------------------------------------
# fs_acute_vs_chronic_pooled
# ---------------------------------------------------------------------------
# Forced swim Acute vs Chronic, ALL animals pooled (both sexes) -- now
# meaningful after the DT4 sex-label correction, since chronic forced swim
# is no longer all-female (was n=9 Female, 0 Male; now 6 Female, 3 Male) and
# acute forced swim remains mostly male (3 Female, 10 Male). Sex is still
# not perfectly balanced between the two conditions, but both now contain
# both sexes, so a direct duration comparison is possible for the first
# time. See fs_acute_vs_chronic_female_only.R for the sex-restricted
# control analysis (unaffected by this correction, still valid).
# Mirrors condition_matrix_stats.R's whole-brain/division/region hierarchy,
# Wilcoxon rank-sum (2 groups) at each level, BH-FDR within each resolution.

library(dplyr)
library(tidyr)

setwd("~/Desktop/Bachelor")
source("prepare_heatmap_data.R")   # mat, long, major_division_per_row, condition_order

keep <- rownames(mat) != "FS"
major_division_per_row <- major_division_per_row[keep]
mat <- mat[keep, , drop = FALSE]
regions <- rownames(mat)

long_fs <- long %>%
  filter(condition %in% c("forced swim · Acute", "forced swim · Chronic")) %>%
  filter(region_acronym != "FS") %>%
  mutate(condition = droplevels(condition))

cat("=== n animals (both sexes pooled) ===\n")
print(long_fs %>% distinct(animal_id_original, sex, condition) %>% count(condition, sex))

## ============================================================
## 1a. Whole-brain level (Wilcoxon)
## ============================================================
animal_overall <- long_fs %>%
  group_by(animal_id_original, channel, condition) %>%
  summarise(overall_median = median(intensity, na.rm = TRUE), .groups = "drop")

wt_wb <- wilcox.test(overall_median ~ condition, data = animal_overall)
cat("\n=== 1a. Whole-brain (pooled), Acute vs Chronic FS ===\n")
cat(sprintf("n animal x channel samples: %d | W = %s, p = %.4f\n",
            nrow(animal_overall), wt_wb$statistic, wt_wb$p.value))
print(animal_overall %>% group_by(condition) %>%
        summarise(n = n(), median_overall = median(overall_median), .groups = "drop"))

write.csv(
  animal_overall %>% group_by(condition) %>%
    summarise(n = n(), median_overall = median(overall_median), .groups = "drop") %>%
    mutate(W = unname(wt_wb$statistic), p_value = wt_wb$p.value),
  "fs_acute_vs_chronic_POOLED_wholebrain.csv", row.names = FALSE
)

## ============================================================
## 1b. Division-level (whole-brain-median normalized, log10, Wilcoxon)
## ============================================================
long_norm <- long_fs %>%
  left_join(animal_overall, by = c("animal_id_original", "channel", "condition")) %>%
  mutate(norm_value = intensity / overall_median)

division_scores <- long_norm %>%
  group_by(animal_id_original, channel, condition,
           major_division = major_division_per_row[match(region_acronym, rownames(mat))]) %>%
  summarise(div_value = median(norm_value, na.rm = TRUE), n_regions = n(), .groups = "drop") %>%
  filter(!is.na(major_division)) %>%
  mutate(log_div_value = log10(div_value + 0.001))

divisions <- unique(division_scores$major_division)
division_results <- vector("list", length(divisions))

for (i in seq_along(divisions)) {
  dv <- divisions[i]
  d <- division_scores %>% filter(major_division == dv)
  if (n_distinct(d$condition) < 2 || any(table(d$condition) < 2)) next
  wt <- wilcox.test(log_div_value ~ condition, data = d)
  medians <- d %>% group_by(condition) %>% summarise(m = median(log_div_value), .groups = "drop")
  division_results[[i]] <- data.frame(
    major_division = dv, W = unname(wt$statistic),
    n_regions_pooled = round(mean(d$n_regions)), p_value = wt$p.value,
    t(setNames(medians$m, paste0("median_log_", medians$condition)))
  )
}

division_df <- bind_rows(division_results) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  arrange(p_value)

write.csv(division_df, "fs_acute_vs_chronic_POOLED_division.csv", row.names = FALSE)

cat("\n=== 1b. Division-level (pooled), n =", nrow(division_df), "divisions ===\n")
print(as.data.frame(division_df[, c("major_division", "W", "p_value", "p_adj_BH")]))
cat(sprintf("\nRaw p < 0.05: %d / %d | BH-FDR p < 0.05: %d / %d\n",
            sum(division_df$p_value < 0.05), nrow(division_df),
            sum(division_df$p_adj_BH < 0.05), nrow(division_df)))

## ============================================================
## 1c. Region level (Wilcoxon)
## ============================================================
region_results <- vector("list", length(regions))

for (i in seq_along(regions)) {
  reg <- regions[i]
  d <- long_fs %>% filter(region_acronym == reg)
  if (n_distinct(d$condition) < 2 || any(table(d$condition) < 2) || nrow(d) < 5) next
  wt <- suppressWarnings(wilcox.test(intensity ~ condition, data = d))
  medians <- d %>% group_by(condition) %>% summarise(m = median(intensity, na.rm = TRUE), .groups = "drop")
  region_results[[i]] <- data.frame(
    region_acronym = reg, W = unname(wt$statistic),
    n_total = nrow(d), p_value = wt$p.value,
    t(setNames(medians$m, paste0("median_", medians$condition)))
  )
}

region_df <- bind_rows(region_results) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj_BH, p_value)

write.csv(region_df, "fs_acute_vs_chronic_POOLED_region.csv", row.names = FALSE)

cat("\n=== 1c. Region-level (pooled), per region ===\n")
cat("Total tests:", nrow(region_df), "\n")
cat("Raw p < 0.05:", sum(region_df$p_value < 0.05, na.rm = TRUE), "\n")
cat("BH-FDR p < 0.05:", sum(region_df$p_adj_BH < 0.05, na.rm = TRUE), "\n")

sig_raw <- region_df %>% filter(p_value < 0.05)
acute_col <- grep("^median_forced.swim...Acute$", names(sig_raw), value = TRUE)
chronic_col <- grep("^median_forced.swim...Chronic$", names(sig_raw), value = TRUE)
cat(sprintf("Direction among raw-significant regions: Acute > Chronic: %d | Chronic > Acute: %d\n",
            sum(sig_raw[[acute_col]] > sig_raw[[chronic_col]]),
            sum(sig_raw[[chronic_col]] > sig_raw[[acute_col]])))

cat("\nTop 20 by raw p-value:\n")
print(head(as.data.frame(region_df %>% arrange(p_value)), 20))

cat("\nSaved:\n")
cat("  fs_acute_vs_chronic_POOLED_wholebrain.csv\n")
cat("  fs_acute_vs_chronic_POOLED_division.csv\n")
cat("  fs_acute_vs_chronic_POOLED_region.csv\n")


# ---------------------------------------------------------------------------
# fs_acute_vs_chronic_female_only
# ---------------------------------------------------------------------------
# Forced swim Acute vs Chronic, restricted to FEMALE animals only -- removes
# the sex confound noted in the general (non-sex-stratified) acute-vs-chronic
# comparison, since chronic forced swim is all-female (n=9) while acute
# forced swim is mostly male (n=3 Female, 10 Male). Comparing females only
# isolates the duration (acute vs chronic) effect from the sex effect.
# Mirrors condition_matrix_stats.R's whole-brain/division/region hierarchy,
# Wilcoxon rank-sum (2 groups) at each level, BH-FDR within each resolution.

library(dplyr)
library(tidyr)

setwd("~/Desktop/Bachelor")
source("prepare_heatmap_data.R")   # mat, long, major_division_per_row, condition_order

keep <- rownames(mat) != "FS"
major_division_per_row <- major_division_per_row[keep]
mat <- mat[keep, , drop = FALSE]
regions <- rownames(mat)

long_f <- long %>%
  filter(sex == "Female", condition %in% c("forced swim · Acute", "forced swim · Chronic")) %>%
  filter(region_acronym != "FS") %>%
  mutate(condition = droplevels(condition))

cat("=== n animals (Female only) ===\n")
print(long_f %>% distinct(animal_id_original, condition) %>% count(condition))

## ============================================================
## 1a. Whole-brain level (Wilcoxon)
## ============================================================
animal_overall_f <- long_f %>%
  group_by(animal_id_original, channel, condition) %>%
  summarise(overall_median = median(intensity, na.rm = TRUE), .groups = "drop")

wt_wb <- wilcox.test(overall_median ~ condition, data = animal_overall_f)
cat("\n=== 1a. Whole-brain (Female only), Acute vs Chronic FS ===\n")
cat(sprintf("n animal x channel samples: %d | W = %s, p = %.4f\n",
            nrow(animal_overall_f), wt_wb$statistic, wt_wb$p.value))
print(animal_overall_f %>% group_by(condition) %>%
        summarise(n = n(), median_overall = median(overall_median), .groups = "drop"))

write.csv(
  animal_overall_f %>% group_by(condition) %>%
    summarise(n = n(), median_overall = median(overall_median), .groups = "drop") %>%
    mutate(W = unname(wt_wb$statistic), p_value = wt_wb$p.value),
  "fs_acute_vs_chronic_FEMALE_wholebrain.csv", row.names = FALSE
)

## ============================================================
## 1b. Division-level (whole-brain-median normalized, log10, Wilcoxon)
## ============================================================
long_norm_f <- long_f %>%
  left_join(animal_overall_f, by = c("animal_id_original", "channel", "condition")) %>%
  mutate(norm_value = intensity / overall_median)

division_scores_f <- long_norm_f %>%
  group_by(animal_id_original, channel, condition,
           major_division = major_division_per_row[match(region_acronym, rownames(mat))]) %>%
  summarise(div_value = median(norm_value, na.rm = TRUE), n_regions = n(), .groups = "drop") %>%
  filter(!is.na(major_division)) %>%
  mutate(log_div_value = log10(div_value + 0.001))

divisions <- unique(division_scores_f$major_division)
division_results_f <- vector("list", length(divisions))

for (i in seq_along(divisions)) {
  dv <- divisions[i]
  d <- division_scores_f %>% filter(major_division == dv)
  if (n_distinct(d$condition) < 2 || any(table(d$condition) < 2)) next
  wt <- wilcox.test(log_div_value ~ condition, data = d)
  medians <- d %>% group_by(condition) %>% summarise(m = median(log_div_value), .groups = "drop")
  division_results_f[[i]] <- data.frame(
    major_division = dv, W = unname(wt$statistic),
    n_regions_pooled = round(mean(d$n_regions)), p_value = wt$p.value,
    t(setNames(medians$m, paste0("median_log_", medians$condition)))
  )
}

division_df_f <- bind_rows(division_results_f) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  arrange(p_value)

write.csv(division_df_f, "fs_acute_vs_chronic_FEMALE_division.csv", row.names = FALSE)

cat("\n=== 1b. Division-level (Female only), n =", nrow(division_df_f), "divisions ===\n")
print(as.data.frame(division_df_f[, c("major_division", "W", "p_value", "p_adj_BH")]))
cat(sprintf("\nRaw p < 0.05: %d / %d | BH-FDR p < 0.05: %d / %d\n",
            sum(division_df_f$p_value < 0.05), nrow(division_df_f),
            sum(division_df_f$p_adj_BH < 0.05), nrow(division_df_f)))

## ============================================================
## 1c. Region level (Wilcoxon)
## ============================================================
region_results_f <- vector("list", length(regions))

for (i in seq_along(regions)) {
  reg <- regions[i]
  d <- long_f %>% filter(region_acronym == reg)
  if (n_distinct(d$condition) < 2 || any(table(d$condition) < 2) || nrow(d) < 5) next
  wt <- suppressWarnings(wilcox.test(intensity ~ condition, data = d))
  medians <- d %>% group_by(condition) %>% summarise(m = median(intensity, na.rm = TRUE), .groups = "drop")
  region_results_f[[i]] <- data.frame(
    region_acronym = reg, W = unname(wt$statistic),
    n_total = nrow(d), p_value = wt$p.value,
    t(setNames(medians$m, paste0("median_", medians$condition)))
  )
}

region_df_f <- bind_rows(region_results_f) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj_BH, p_value)

write.csv(region_df_f, "fs_acute_vs_chronic_FEMALE_region.csv", row.names = FALSE)

cat("\n=== 1c. Region-level (Female only), per region ===\n")
cat("Total tests:", nrow(region_df_f), "\n")
cat("Raw p < 0.05:", sum(region_df_f$p_value < 0.05, na.rm = TRUE), "\n")
cat("BH-FDR p < 0.05:", sum(region_df_f$p_adj_BH < 0.05, na.rm = TRUE), "\n")

sig_raw <- region_df_f %>% filter(p_value < 0.05)
acute_col <- grep("^median_forced.swim...Acute$", names(sig_raw), value = TRUE)
chronic_col <- grep("^median_forced.swim...Chronic$", names(sig_raw), value = TRUE)
cat(sprintf("Direction among raw-significant regions: Acute > Chronic: %d | Chronic > Acute: %d\n",
            sum(sig_raw[[acute_col]] > sig_raw[[chronic_col]]),
            sum(sig_raw[[chronic_col]] > sig_raw[[acute_col]])))

cat("\nTop 15 by raw p-value:\n")
print(head(as.data.frame(region_df_f %>% arrange(p_value)), 15))

cat("\nSaved:\n")
cat("  fs_acute_vs_chronic_FEMALE_wholebrain.csv\n")
cat("  fs_acute_vs_chronic_FEMALE_division.csv\n")
cat("  fs_acute_vs_chronic_FEMALE_region.csv\n")

