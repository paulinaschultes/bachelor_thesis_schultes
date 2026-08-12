# =============================================================================
# Combined ranking of high-intensity, low-variability regions
#
# Consolidated from the following original scripts:
#   - top_regions_low_deviation.R
#   - top_regions_per_condition.R
#   - sex_top_regions_low_deviation.R
#   - condition_low_deviation_and_histology.R
# =============================================================================


# ---------------------------------------------------------------------------
# top_regions_low_deviation
# ---------------------------------------------------------------------------
# Top regions across all 9 groups: highest median intensity combined with
# lowest cross-group deviation (CV = SD/median), separately for Cre and tTA.
# "Top" = low combined rank of (rank by high median) + (rank by low CV) --
# i.e. regions that are both strongly and *consistently* labeled regardless
# of stress condition.
#
# Note: region "FS" was excluded -- it is the injection region itself, and
# all other regions' intensities are already expressed as a ratio relative
# to FS. FS is therefore trivially 1.0 (SD = 0) in every animal by
# construction, not because it is a biologically stable region -- including
# it would trivially "win" any low-deviation ranking for the wrong reason.
# Only regions present in all 9 groups are included, so the "across all
# groups" comparison is well-defined.

library(dplyr)
library(ggplot2)
library(ggrepel)
library(grid)

df <- read.csv("~/Desktop/ALL_DT.csv", stringsAsFactors = FALSE)

## --- Median per region x group x channel ---
group_med <- df %>%
  group_by(region_id, region_acronym, division, animal_group) %>%
  summarise(Cre = median(trapcre_intensity, na.rm = TRUE),
            tTA = median(traptta_intensity, na.rm = TRUE), .groups = "drop") %>%
  filter(region_acronym != "FS")

complete_regions <- group_med %>% count(region_acronym) %>% filter(n == 9) %>% pull(region_acronym)
gm <- group_med %>% filter(region_acronym %in% complete_regions)

division_lookup <- gm %>% distinct(region_acronym, division)

## --- Across-group median & CV per region, per channel ---
summarize_channel <- function(ch) {
  gm %>%
    group_by(region_acronym) %>%
    summarise(median_across = median(.data[[ch]], na.rm = TRUE),
              sd_across = sd(.data[[ch]], na.rm = TRUE), .groups = "drop") %>%
    mutate(cv = sd_across / median_across) %>%
    filter(is.finite(cv), median_across > 0) %>%
    mutate(rank_median = rank(-median_across), rank_cv = rank(cv),
           combined = rank_median + rank_cv) %>%
    arrange(combined) %>%
    left_join(division_lookup, by = "region_acronym")
}

cre_summary <- summarize_channel("Cre") %>% mutate(channel = "Cre")
tta_summary <- summarize_channel("tTA") %>% mutate(channel = "tTA")

top_n <- 20
cre_top <- cre_summary %>% slice_head(n = top_n)
tta_top <- tta_summary %>% slice_head(n = top_n)

write.csv(bind_rows(cre_top, tta_top),
          "~/Desktop/top_regions_low_deviation.csv", row.names = FALSE)

## --- Scatter: median vs CV, top regions highlighted ---
plot_data <- bind_rows(cre_summary %>% mutate(channel = "Cre"),
                        tta_summary %>% mutate(channel = "tTA")) %>%
  mutate(is_top = combined <= sort(combined)[top_n], .by = channel)

p <- ggplot(plot_data, aes(x = median_across, y = cv)) +
  geom_point(aes(color = is_top, size = is_top), alpha = 0.6) +
  geom_text_repel(data = plot_data %>% filter(is_top),
                   aes(label = region_acronym), size = 2.6, max.overlaps = 20,
                   segment.size = 0.3, color = "#4A148C") +
  scale_x_log10(name = "Median intensity across 9 groups (log scale)") +
  scale_y_log10(name = "Coefficient of variation (SD / median) across groups, log scale") +
  scale_color_manual(values = c("FALSE" = "#F48FB1", "TRUE" = "#4A148C"), guide = "none") +
  scale_size_manual(values = c("FALSE" = 1, "TRUE" = 2), guide = "none") +
  facet_wrap(~channel, scales = "free") +
  labs(title = sprintf("Top %d regions per channel: high median & low cross-group deviation", top_n),
       subtitle = "Highlighted = top regions by combined rank ('FS' = injection region, excluded: trivially 1.0 by definition)") +
  theme_minimal(base_size = 10) +
  theme(strip.text = element_text(face = "bold"))

ggsave("~/Desktop/top_regions_low_deviation.pdf", p, width = 10, height = 5.5, units = "in")
ggsave("~/Desktop/top_regions_low_deviation.png", p, width = 10, height = 5.5, units = "in", dpi = 300)

cat("\n=== Top", top_n, "Cre regions (high median, low CV) ===\n")
print(as.data.frame(cre_top %>% select(region_acronym, division, median_across, cv, combined)))
cat("\n=== Top", top_n, "tTA regions (high median, low CV) ===\n")
print(as.data.frame(tta_top %>% select(region_acronym, division, median_across, cv, combined)))


# ---------------------------------------------------------------------------
# top_regions_per_condition
# ---------------------------------------------------------------------------
# Top regions per stressor condition (not per group): pools the specific
# (group, channel) combinations that share the same real stress condition
# (e.g. "restraint" pools DT1_Cre, DT2_Cre, DT3_tTA, TRDT1_Cre, ... -- mixed
# groups AND mixed TRAP systems, wherever that condition was applied), then
# ranks regions by high median + low CV (= greatest agreement) within that
# pooled set.
#
# Condition assignment comes from group_condition_corrected.csv (one Cre
# condition + one tTA condition per GROUP, confirmed constant across all
# animals within a group via animalid.xlsx). NOTE: ALL_DT.csv's own
# trapcre_condition/traptta_condition columns are NOT used here -- they
# contain a data error where some individual animals' Cre/tTA condition
# labels are swapped relative to their groupmates.
#
# 'FS' (injection region, trivially 1.0 by definition -- see
# top_regions_low_deviation.R) is excluded throughout.

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

df <- read.csv("~/Desktop/ALL_DT.csv", stringsAsFactors = FALSE)
cond_map <- read.csv("~/Desktop/group_condition_corrected.csv", stringsAsFactors = FALSE)

## --- Median per region x group x channel (group-level is valid: condition
## is constant within a group) ---
group_med <- df %>%
  group_by(region_id, region_acronym, division, animal_group) %>%
  summarise(Cre = median(trapcre_intensity, na.rm = TRUE),
            tTA = median(traptta_intensity, na.rm = TRUE), .groups = "drop") %>%
  filter(region_acronym != "FS")

gm_long <- group_med %>%
  pivot_longer(cols = c(Cre, tTA), names_to = "channel", values_to = "value")

division_lookup <- group_med %>% distinct(region_acronym, division)
conditions <- sort(unique(cond_map$Condition))

## --- Per-condition: pool the qualifying (group, channel) samples per region ---
top_n <- 15
all_results <- list()

for (cond in conditions) {
  pairs <- cond_map %>% filter(Condition == cond) %>% rename(animal_group = Group, channel = Channel)
  pooled <- gm_long %>% inner_join(pairs, by = c("animal_group", "channel"))
  n_pairs <- nrow(pairs)

  summary_cond <- pooled %>%
    group_by(region_acronym) %>%
    summarise(n_samples = n(),
              median_across = median(value, na.rm = TRUE),
              sd_across = sd(value, na.rm = TRUE), .groups = "drop") %>%
    filter(n_samples == n_pairs) %>%
    mutate(cv = sd_across / median_across) %>%
    filter(is.finite(cv), median_across > 0) %>%
    mutate(rank_median = rank(-median_across), rank_cv = rank(cv),
           combined = rank_median + rank_cv) %>%
    arrange(combined) %>%
    left_join(division_lookup, by = "region_acronym") %>%
    mutate(condition = cond, n_pairs = n_pairs)

  all_results[[cond]] <- summary_cond
}

full_results <- bind_rows(all_results)
top_per_condition <- full_results %>% group_by(condition) %>% slice_head(n = top_n) %>% ungroup()

write.csv(top_per_condition %>%
            select(condition, region_acronym, division, median_across, cv, n_pairs, combined),
          "~/Desktop/top_regions_per_condition.csv", row.names = FALSE)

## --- Scatter: median vs CV per condition, top regions highlighted ---
full_results <- full_results %>% mutate(is_top = combined <= sort(combined)[top_n], .by = condition)

p <- ggplot(full_results, aes(x = median_across, y = cv)) +
  geom_point(aes(color = is_top, size = is_top), alpha = 0.6) +
  geom_text_repel(data = full_results %>% filter(is_top),
                   aes(label = region_acronym), size = 2.4, max.overlaps = 25,
                   segment.size = 0.3, color = "#4A148C") +
  scale_x_log10(name = "Median intensity (log scale)") +
  scale_y_log10(name = "CV across pooled samples (log scale)") +
  scale_color_manual(values = c("FALSE" = "#F48FB1", "TRUE" = "#4A148C"), guide = "none") +
  scale_size_manual(values = c("FALSE" = 0.8, "TRUE" = 1.8), guide = "none") +
  facet_wrap(~condition, scales = "free", ncol = 2) +
  labs(title = sprintf("Top %d regions per stressor\ncondition (high median, low CV)", top_n),
       subtitle = "Pooled across groups/channels sharing that condition\n(corrected group-level mapping; 'FS' excluded)") +
  theme_minimal(base_size = 8) +
  theme(strip.text = element_text(face = "bold", size = 7.5),
        plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 7.5))

ggsave("~/Desktop/top_regions_per_condition.pdf", p, width = 12, height = 7, units = "in")
ggsave("~/Desktop/top_regions_per_condition.png", p, width = 12, height = 7, units = "in", dpi = 300)
ggsave("~/Desktop/top_regions_per_condition_A4.pdf", p, width = 8.27, height = 11.69, units = "in")
ggsave("~/Desktop/top_regions_per_condition_A4.png", p, width = 8.27, height = 11.69, units = "in", dpi = 300)

for (cond in conditions) {
  cat("\n===", cond, "( n pairs =", unique(top_per_condition$n_pairs[top_per_condition$condition == cond]), ") ===\n")
  print(as.data.frame(top_per_condition %>% filter(condition == cond) %>%
                         select(region_acronym, division, median_across, cv)))
}


# ---------------------------------------------------------------------------
# sex_top_regions_low_deviation
# ---------------------------------------------------------------------------
# Top 40 regions per SEX: high median intensity + LOW deviation across the
# 4 (non-chronic-FS) conditions -- sex-specific counterpart to
# condition_low_deviation_and_histology.R part 1 (which pooled both sexes).
# "forced swim · Chronic" excluded (all-Female, no Male cell), same 4-
# condition set used throughout the condition x sex analysis, so Female and
# Male are compared on an equal footing (same 4 conditions each).

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

setwd("~/Desktop/Bachelor")
source("prepare_heatmap_data.R")   # long, major_division_per_row, mat, condition_order
long <- long %>% filter(region_acronym != "FS")

sex_order <- c("Female", "Male")
condition_order_no_chronicFS <- setdiff(condition_order, "social defeat · Acute")  # social defeat is now single-sex (Male only); excluded here so only conditions comparable across both sexes remain. Chronic FS is included since its sex labels were corrected.
top_n <- 40

major_cols <- c(
  "Isocortex"             = "#5B7C99", "Olfactory areas"       = "#8C6D46",
  "Hippocampal formation" = "#6B8F71", "Cortical subplate"     = "#9C6B85",
  "Striatum"              = "#4F8F8B", "Pallidum"              = "#7A6B9C",
  "Thalamus"              = "#B08B5A", "Hypothalamus"          = "#A85C5C",
  "Midbrain"              = "#6B6B99", "Pons"                  = "#8FA05E",
  "Medulla"               = "#5E8FA0", "Cerebellum"            = "#9C7A5B",
  "fiber tracts"          = "#A9A9A9", "ventricular systems"   = "#C4B7A6"
)
pink_scale <- c("#FDF0F5", "#FBC9DE", "#F48FB1", "#E05780", "#C2185B", "#8E24AA", "#4A148C")

region_summaries <- list()

for (sx in sex_order) {
  ## region x condition matrix, median intensity, THIS SEX only.
  ## Conditions with zero animals of this sex (e.g. "social defeat · Acute"
  ## is now all-Male after the DT4 sex correction) are dropped per-sex,
  ## rather than assumed present for both sexes.
  sx_conditions <- long %>% filter(sex == sx, condition %in% condition_order_no_chronicFS) %>%
    distinct(condition) %>% pull(condition) %>% as.character()
  sx_conditions <- condition_order_no_chronicFS[condition_order_no_chronicFS %in% sx_conditions]
  if (length(setdiff(condition_order_no_chronicFS, sx_conditions)) > 0) {
    cat(sprintf("%s: no animals for %s -- excluded from this sex's condition set\n",
                sx, paste(setdiff(condition_order_no_chronicFS, sx_conditions), collapse = ", ")))
  }

  agg_sx <- long %>%
    filter(sex == sx, condition %in% sx_conditions) %>%
    group_by(region_acronym, condition) %>%
    summarise(intensity = median(intensity, na.rm = TRUE), .groups = "drop")

  wide_sx <- agg_sx %>% pivot_wider(names_from = condition, values_from = intensity)
  mat_sx <- as.matrix(wide_sx[, sx_conditions])
  rownames(mat_sx) <- wide_sx$region_acronym
  mat_sx <- mat_sx[rowSums(is.na(mat_sx)) == 0, , drop = FALSE]   # regions present in all 4 conditions for this sex

  rs <- data.frame(
    region_acronym = rownames(mat_sx),
    major_division = major_division_per_row[match(rownames(mat_sx), rownames(mat))],
    median_across = apply(mat_sx, 1, median, na.rm = TRUE),
    sd_across = apply(mat_sx, 1, sd, na.rm = TRUE)
  ) %>%
    mutate(cv = sd_across / median_across) %>%
    filter(is.finite(cv), median_across > 0) %>%
    mutate(rank_median = rank(-median_across), rank_cv = rank(cv),
           combined = rank_median + rank_cv, sex = sx) %>%
    arrange(combined)

  region_summaries[[sx]] <- rs
  cat(sprintf("=== %s: top %d regions (high intensity, low deviation across %d conditions) ===\n",
              sx, top_n, length(sx_conditions)))
  print(head(rs %>% select(region_acronym, major_division, median_across, cv, combined), top_n))
}

region_summary_all <- bind_rows(region_summaries)
write.csv(region_summary_all, "sex_top_regions_low_deviation.csv", row.names = FALSE)

## --- Scatter plot: median vs CV, top 40 highlighted, faceted by sex ---
FONT_PT <- 14
mm_size <- function(pt) pt / (72.27 / 25.4)  # ggplot geom text 'size' aesthetic is in mm

top40_all <- region_summary_all %>% group_by(sex) %>% slice_head(n = top_n) %>% ungroup()

p <- ggplot(region_summary_all, aes(x = cv, y = median_across, color = combined)) +
  geom_point(alpha = 0.75, size = 2.2) +
  ggrepel::geom_text_repel(data = top40_all, aes(x = cv, y = median_across, label = region_acronym),
                            color = "#2A2A2A", size = mm_size(FONT_PT), fontface = "bold", max.overlaps = Inf,
                            box.padding = 0.4, point.padding = 0.2, force = 3, force_pull = 0.5,
                            min.segment.length = 0, segment.color = "grey50", inherit.aes = FALSE) +
  scale_color_gradientn(colours = rev(pink_scale), name = "Combined rank\n(darker = better)") +
  scale_y_log10() +
  coord_cartesian(xlim = c(0, 0.8), ylim = c(0.08, 3)) +
  facet_wrap(~sex, ncol = 2) +
  labs(title = sprintf("Top %d regions per sex: high intensity & low deviation across conditions",
                        top_n),
       subtitle = "4 conditions each (forced swim Acute/Chronic, restraint, tail suspension); social defeat excluded (single-sex, Male only)",
       x = "Coefficient of variation across conditions", y = "Median intensity across conditions (log scale)") +
  theme_minimal(base_size = FONT_PT) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = FONT_PT),
    plot.title = element_text(face = "bold", size = FONT_PT + 2),
    plot.subtitle = element_text(size = FONT_PT - 1),
    legend.position = "bottom",
    legend.text = element_text(size = FONT_PT),
    legend.title = element_text(size = FONT_PT),
    legend.key.width = unit(3, "cm")
  ) +
  guides(color = guide_colorbar(barwidth = unit(10, "cm"), barheight = unit(0.8, "cm")))

ggsave("sex_top_regions_low_deviation.png", p, width = 16, height = 9, units = "in", dpi = 300, limitsize = FALSE)

cat("\nSaved: sex_top_regions_low_deviation.csv, sex_top_regions_low_deviation.png\n")

## --- Overlap between the two sex-specific top-40 sets ---
f_top <- top40_all %>% filter(sex == "Female") %>% pull(region_acronym)
m_top <- top40_all %>% filter(sex == "Male") %>% pull(region_acronym)
common <- intersect(f_top, m_top)
cat(sprintf("\nShared in both sexes' top %d: %d regions (Jaccard = %.3f)\n",
            top_n, length(common), length(common) / length(union(f_top, m_top))))
print(common)


# ---------------------------------------------------------------------------
# condition_low_deviation_and_histology
# ---------------------------------------------------------------------------
# Two analyses using the condition-level matrix (prepare_heatmap_data.R):
#   1. Top regions with the LOWEST deviation across the 5 conditions
#      (opposite of condition-discriminating: high median intensity
#      combined with low coefficient of variation across conditions --
#      regions that are strongly AND consistently labeled regardless of
#      stressor). Condition-based update of top_regions_low_deviation.R,
#      which used animal_group x channel instead of condition.
#   2. Histological (major-division) intensity comparison, independent of
#      condition: are there systematic intensity differences between the
#      14 anatomical divisions themselves, pooling across all conditions?

library(dplyr)
library(ggplot2)
library(ggrepel)

setwd("~/Desktop/Bachelor")
source("prepare_heatmap_data.R")   # mat, long, major_division_per_row, condition_order

keep <- rownames(mat) != "FS"
mat <- mat[keep, , drop = FALSE]
major_division_per_row <- major_division_per_row[keep]
long <- long %>% filter(region_acronym != "FS")

## ============================================================
## 1. Top regions: high median intensity + LOW deviation ACROSS CONDITIONS
## ============================================================
region_summary <- data.frame(
  region_acronym = rownames(mat),
  major_division = major_division_per_row,
  median_across = apply(mat, 1, median, na.rm = TRUE),
  sd_across = apply(mat, 1, sd, na.rm = TRUE)
) %>%
  mutate(cv = sd_across / median_across) %>%
  filter(is.finite(cv), median_across > 0) %>%
  mutate(rank_median = rank(-median_across), rank_cv = rank(cv),
         combined = rank_median + rank_cv) %>%
  arrange(combined)

write.csv(region_summary, "condition_top_regions_low_deviation.csv", row.names = FALSE)

cat("=== 1. Top 20 regions: high intensity + low deviation ACROSS CONDITIONS ===\n")
print(head(region_summary, 20))

major_cols <- c(
  "Isocortex"             = "#5B7C99", "Olfactory areas"       = "#8C6D46",
  "Hippocampal formation" = "#6B8F71", "Cortical subplate"     = "#9C6B85",
  "Striatum"              = "#4F8F8B", "Pallidum"              = "#7A6B9C",
  "Thalamus"              = "#B08B5A", "Hypothalamus"          = "#A85C5C",
  "Midbrain"              = "#6B6B99", "Pons"                  = "#8FA05E",
  "Medulla"               = "#5E8FA0", "Cerebellum"            = "#9C7A5B",
  "fiber tracts"          = "#A9A9A9", "ventricular systems"   = "#C4B7A6"
)

pink_scale <- c("#FDF0F5", "#FBC9DE", "#F48FB1", "#E05780", "#C2185B", "#8E24AA", "#4A148C")

FONT_PT <- 14
mm_size <- function(pt) pt / (72.27 / 25.4)  # ggplot geom text 'size' aesthetic is in mm

top20 <- head(region_summary, 20)
p1 <- ggplot(region_summary, aes(x = cv, y = median_across, color = combined)) +
  geom_point(alpha = 0.75, size = 2.2) +
  ggrepel::geom_text_repel(data = top20, mapping = aes(x = cv, y = median_across, label = region_acronym),
                            color = "#2A2A2A", size = mm_size(FONT_PT), fontface = "bold", max.overlaps = Inf,
                            box.padding = 0.6, point.padding = 0.3, force = 4, force_pull = 0.5,
                            min.segment.length = 0, segment.color = "grey50", inherit.aes = FALSE) +
  scale_color_gradientn(colours = rev(pink_scale), name = "Combined rank\n(darker = better)") +
  scale_y_log10() +
  coord_cartesian(xlim = c(0, 1.2), ylim = c(0.08, 3)) +   # centered on the top-20 cluster (max CV 0.33, median range 0.37-2.07)
  labs(x = "Coefficient of variation across the 5 conditions", y = "Median intensity across conditions (log scale)") +
  theme_minimal(base_size = FONT_PT) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey93", linewidth = 0.3),
    plot.title = element_text(face = "bold", size = FONT_PT + 2),
    plot.subtitle = element_text(size = FONT_PT, color = "grey30"),
    axis.title = element_text(size = FONT_PT),
    axis.text = element_text(size = FONT_PT),
    legend.position = "bottom",
    legend.key.width = unit(3, "cm"),
    legend.text = element_text(size = FONT_PT),
    legend.title = element_text(size = FONT_PT)
  ) +
  guides(color = guide_colorbar(barwidth = unit(12, "cm"), barheight = unit(0.8, "cm")))
ggsave("condition_top_regions_low_deviation.png", p1, width = 12, height = 9, units = "in", dpi = 300, limitsize = FALSE)

cat("\nSaved: condition_top_regions_low_deviation.csv/.png\n")

## ============================================================
## 2. Histological (major-division) intensity comparison, condition-independent
## ============================================================
## Pool across ALL conditions/animals/channels: does overall intensity
## differ systematically between the 14 major anatomical divisions?
div_data <- long %>%
  mutate(major_division = major_division_per_row[match(region_acronym, rownames(mat))]) %>%
  filter(!is.na(major_division))

kt_div <- kruskal.test(intensity ~ major_division, data = div_data)
cat(sprintf("\n=== 2. Histological (major division) intensity comparison, ALL conditions pooled ===\n"))
cat(sprintf("Kruskal-Wallis across 14 divisions: chi-sq = %.2f, df = %d, p = %.2e\n",
            kt_div$statistic, kt_div$parameter, kt_div$p.value))

div_medians <- div_data %>%
  group_by(major_division) %>%
  summarise(median_intensity = median(intensity, na.rm = TRUE), n = n(), .groups = "drop") %>%
  arrange(desc(median_intensity))
print(div_medians)
write.csv(div_medians, "histological_division_intensity_medians.csv", row.names = FALSE)

## Pairwise post-hoc: Dunn-like via pairwise Wilcoxon, BH-corrected globally
div_pairs <- combn(unique(div_data$major_division), 2, simplify = FALSE)
div_pair_results <- vector("list", length(div_pairs))
for (i in seq_along(div_pairs)) {
  d1 <- div_pairs[[i]][1]; d2 <- div_pairs[[i]][2]
  v1 <- div_data$intensity[div_data$major_division == d1]
  v2 <- div_data$intensity[div_data$major_division == d2]
  wt <- suppressWarnings(wilcox.test(v1, v2))
  div_pair_results[[i]] <- data.frame(division1 = d1, division2 = d2, p_value = wt$p.value)
}
div_pair_df <- bind_rows(div_pair_results) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj_BH)
write.csv(div_pair_df, "histological_division_pairwise_wilcoxon.csv", row.names = FALSE)

cat(sprintf("\nPairwise division comparisons: %d total, %d significant after BH-FDR\n",
            nrow(div_pair_df), sum(div_pair_df$p_adj_BH < 0.05, na.rm = TRUE)))

p2 <- ggplot(div_data, aes(x = reorder(major_division, intensity, median), y = intensity, fill = major_division)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3) +
  scale_fill_manual(values = major_cols, guide = "none") +
  scale_y_log10() +
  coord_flip() +
  labs(title = "Histological (major division) intensity, all conditions pooled",
       subtitle = sprintf("Kruskal-Wallis chi-sq = %.1f, p = %.2e", kt_div$statistic, kt_div$p.value),
       x = NULL, y = "Intensity (log scale)") +
  theme_minimal(base_size = 10)
ggsave("histological_division_intensity_boxplot.png", p2, width = 9, height = 6, units = "in", dpi = 300)

cat("\nSaved: histological_division_intensity_medians.csv, histological_division_pairwise_wilcoxon.csv, histological_division_intensity_boxplot.png\n")

