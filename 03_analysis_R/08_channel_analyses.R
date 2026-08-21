# =============================================================================
# THESIS FIGURES PRODUCED BY THIS SCRIPT
#   Figure 26      No significant differences in whole-brain median intensity between reporter channels
#
# Note: figure numbers appearing in the comments further down refer to earlier
# drafts and are NOT the final numbering. The list above is authoritative.
# =============================================================================

# =============================================================================
# Reporter channel comparison (TRAP-Cre vs TRAP-tTA)
#
# Consolidated from the following original scripts:
#   - channel_bias_check.R
#   - channel_diff_within_condition_full.R
#   - channel_validity_boxplots.R
#   - channel_validity_visualization.R
#   - plot_cre_tta_correlation.R
# =============================================================================


# ---------------------------------------------------------------------------
# channel_bias_check
# ---------------------------------------------------------------------------
# Technical validity check: is there a systematic Cre-vs-tTA REPORTER/
# CHANNEL bias, independent of which stressor condition each channel
# happens to represent? Cre and tTA are two different transgenic reporter
# systems (different antibody/amplification chemistry in principle), so a
# whole-brain-level or region-level systematic shift between channels --
# regardless of condition -- would indicate a technical artifact rather
# than a biological (stressor) effect. This uses ALL animals (not just
# those with a mapped Condition), since it deliberately ignores condition
# and asks only about the channel itself.
#
#   1. Whole-brain level: paired test (Wilcoxon signed-rank + paired t)
#      comparing each animal's OWN overall (whole-brain median) Cre value
#      vs. its own overall tTA value. n = 27 animals, paired by animal.
#   2. Region level: same paired comparison repeated per region (paired by
#      animal), BH-FDR corrected, to check whether specific regions show a
#      channel-driven signal independent of condition.

library(dplyr)
library(tidyr)
library(ggplot2)

df <- read.csv("~/Desktop/ALL_DT.csv", stringsAsFactors = FALSE)

## --- 1. Whole-brain-level paired channel comparison ---
animal_overall <- df %>%
  select(animal_id_original, trapcre_intensity, traptta_intensity) %>%
  pivot_longer(cols = c(trapcre_intensity, traptta_intensity), names_to = "channel", values_to = "value") %>%
  mutate(channel = ifelse(channel == "trapcre_intensity", "Cre", "tTA")) %>%
  group_by(animal_id_original, channel) %>%
  summarise(overall_median = median(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = channel, values_from = overall_median)

cat("=== Whole-brain-level Cre vs tTA (paired by animal, n =", nrow(animal_overall), ") ===\n")
wt <- wilcox.test(animal_overall$Cre, animal_overall$tTA, paired = TRUE)
tt <- t.test(animal_overall$Cre, animal_overall$tTA, paired = TRUE)
cat(sprintf("Wilcoxon signed-rank: V = %.1f, p = %.4f\n", wt$statistic, wt$p.value))
cat(sprintf("Paired t-test: t = %.3f, df = %d, p = %.4f, mean diff (Cre - tTA) = %.4f\n",
            tt$statistic, tt$parameter, tt$p.value, tt$estimate))
cat(sprintf("Median Cre: %.4f | Median tTA: %.4f | Median Cre/tTA ratio: %.3f\n",
            median(animal_overall$Cre), median(animal_overall$tTA),
            median(animal_overall$Cre / animal_overall$tTA)))

write.csv(animal_overall, "~/Desktop/channel_bias_wholebrain.csv", row.names = FALSE)

p1 <- ggplot(animal_overall, aes(x = tTA, y = Cre)) +
  geom_point(size = 2.5, color = "#C2185B", alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  labs(title = "Whole-brain median intensity: Cre vs. tTA channel (paired by animal)",
       subtitle = sprintf("Wilcoxon signed-rank p = %.3f | n = %d animals | dashed line = equality", wt$p.value, nrow(animal_overall)),
       x = "tTA whole-brain median", y = "Cre whole-brain median") +
  theme_minimal(base_size = 11)
ggsave("~/Desktop/channel_bias_wholebrain.png", p1, width = 6.5, height = 6, units = "in", dpi = 300)

## --- 2. Region-level paired channel comparison, BH-FDR ---
region_wide <- df %>%
  filter(region_acronym != "FS") %>%
  select(animal_id_original, region_acronym, trapcre_intensity, traptta_intensity)

regions <- unique(region_wide$region_acronym)
results <- vector("list", length(regions))

for (i in seq_along(regions)) {
  reg <- regions[i]
  d <- region_wide %>% filter(region_acronym == reg) %>%
    filter(!is.na(trapcre_intensity), !is.na(traptta_intensity))
  if (nrow(d) < 5) next
  wt_r <- suppressWarnings(wilcox.test(d$trapcre_intensity, d$traptta_intensity, paired = TRUE))
  results[[i]] <- data.frame(region_acronym = reg, n = nrow(d),
                              median_Cre = median(d$trapcre_intensity), median_tTA = median(d$traptta_intensity),
                              p_value = wt_r$p.value)
}

region_results <- bind_rows(results) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"),
         log2_ratio = log2(median_Cre / median_tTA)) %>%
  arrange(p_adj_BH)

write.csv(region_results, "~/Desktop/channel_bias_per_region.csv", row.names = FALSE)

cat("\n=== Region-level Cre vs tTA (paired, n regions =", nrow(region_results), ") ===\n")
cat(sprintf("Regions with BH-FDR p < 0.05: %d / %d\n", sum(region_results$p_adj_BH < 0.05, na.rm = TRUE), nrow(region_results)))
cat(sprintf("Regions with raw p < 0.05: %d / %d\n", sum(region_results$p_value < 0.05, na.rm = TRUE), nrow(region_results)))
cat("\nTop 15 by raw p-value:\n")
print(as.data.frame(head(region_results, 15)))

cat("\nSaved:\n  ~/Desktop/channel_bias_wholebrain.csv/.png\n  ~/Desktop/channel_bias_per_region.csv\n")


# ---------------------------------------------------------------------------
# channel_diff_within_condition_full
# ---------------------------------------------------------------------------
# Channel (Cre vs tTA) difference test, WITHIN each condition -- extends the
# existing region-level-only test (condition_matrix_stats.R section 4b,
# condition_channel_diff_by_region.csv) to the same whole-brain/division/
# region hierarchy used for the condition and sex tests elsewhere.
#
# Unlike channel_bias_check.R (condition-independent, paired by animal,
# pools ALL animals regardless of which condition each channel happens to
# be tagged with -- some animals' Cre channel and tTA channel reflect
# DIFFERENT stress sessions), this restricts the comparison to animals
# whose channel is tagged with the SAME condition, so "channel" is not
# confounded with "which session was captured". Cre and tTA within one
# condition come from different animals here (unpaired), since only 3 of
# 5 conditions have both channels represented at all (forced swim chronic
# is Cre-only, social defeat is tTA-only -- no channel contrast possible
# for those two).

library(dplyr)
library(tidyr)

setwd("~/Desktop/Bachelor")
source("prepare_heatmap_data.R")   # mat, long, major_division_per_row, condition_order

keep <- rownames(mat) != "FS"
major_division_per_row <- major_division_per_row[keep]
mat <- mat[keep, , drop = FALSE]
regions <- rownames(mat)

channel_counts <- long %>% distinct(condition, channel) %>% count(condition)
testable_conditions <- as.character(channel_counts$condition[channel_counts$n == 2])
cat("Conditions with both channels (testable):", paste(testable_conditions, collapse = ", "), "\n")
cat("Conditions with only one channel (skipped):",
    paste(setdiff(as.character(condition_order), testable_conditions), collapse = ", "), "\n\n")

wholebrain_list <- vector("list", length(testable_conditions))
division_list   <- vector("list", length(testable_conditions))

for (ci in seq_along(testable_conditions)) {
  cond <- testable_conditions[ci]
  long_c <- long %>% filter(condition == cond, region_acronym != "FS")

  ## --- whole-brain: each animal x channel's own whole-brain median ---
  animal_overall_c <- long_c %>%
    group_by(animal_id_original, channel) %>%
    summarise(overall_median = median(intensity, na.rm = TRUE), .groups = "drop")

  wt_wb <- wilcox.test(overall_median ~ channel, data = animal_overall_c)
  medians_wb <- animal_overall_c %>% group_by(channel) %>%
    summarise(n = n(), median_overall = median(overall_median), .groups = "drop")
  wholebrain_list[[ci]] <- medians_wb %>%
    mutate(condition = cond, W = unname(wt_wb$statistic), p_value = wt_wb$p.value)

  ## --- division: whole-brain-median normalized per animal x channel ---
  long_norm_c <- long_c %>%
    left_join(animal_overall_c, by = c("animal_id_original", "channel")) %>%
    mutate(norm_value = intensity / overall_median)

  division_scores_c <- long_norm_c %>%
    group_by(animal_id_original, channel,
             major_division = major_division_per_row[match(region_acronym, rownames(mat))]) %>%
    summarise(div_value = median(norm_value, na.rm = TRUE), n_regions = n(), .groups = "drop") %>%
    filter(!is.na(major_division)) %>%
    mutate(log_div_value = log10(div_value + 0.001))

  divisions <- unique(division_scores_c$major_division)
  div_results <- vector("list", length(divisions))
  for (i in seq_along(divisions)) {
    dv <- divisions[i]
    d <- division_scores_c %>% filter(major_division == dv)
    if (n_distinct(d$channel) < 2 || any(table(d$channel) < 2)) next
    wt <- wilcox.test(log_div_value ~ channel, data = d)
    medians <- d %>% group_by(channel) %>% summarise(m = median(log_div_value), .groups = "drop")
    div_results[[i]] <- data.frame(
      condition = cond, major_division = dv, W = unname(wt$statistic),
      n_regions_pooled = round(mean(d$n_regions)), p_value = wt$p.value,
      t(setNames(medians$m, paste0("median_log_", medians$channel)))
    )
  }
  division_list[[ci]] <- bind_rows(div_results) %>%
    mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

  cat(sprintf("=== %s ===\n", cond))
  cat(sprintf("  whole-brain: n=%d, W=%.1f, p=%.4f (median Cre=%.4f, tTA=%.4f)\n",
              nrow(animal_overall_c), wt_wb$statistic, wt_wb$p.value,
              medians_wb$median_overall[medians_wb$channel == "Cre"],
              medians_wb$median_overall[medians_wb$channel == "tTA"]))
  cat(sprintf("  division:    %d/%d raw p<0.05, %d/%d BH-FDR p<0.05\n\n",
              sum(division_list[[ci]]$p_value < 0.05, na.rm = TRUE), nrow(division_list[[ci]]),
              sum(division_list[[ci]]$p_adj_BH < 0.05, na.rm = TRUE), nrow(division_list[[ci]])))
}

wholebrain_df <- bind_rows(wholebrain_list) %>%
  select(condition, channel, n, median_overall, W, p_value) %>%
  arrange(condition)
division_df <- bind_rows(division_list) %>% arrange(condition, p_value)

write.csv(wholebrain_df, "channel_diff_within_condition_wholebrain.csv", row.names = FALSE)
write.csv(division_df,   "channel_diff_within_condition_division.csv", row.names = FALSE)

cat("=== Summary: whole-brain Cre vs tTA per condition ===\n")
print(as.data.frame(wholebrain_df %>% distinct(condition, W, p_value)))

## Region-level test already exists (condition_matrix_stats.R, section 4b) --
## just re-summarize it here for a complete 3-level report.
region_df <- read.csv("condition_channel_diff_by_region.csv")
cat("\n=== Summary: region-level Cre vs tTA per condition (from condition_channel_diff_by_region.csv) ===\n")
print(region_df %>% group_by(condition) %>%
        summarise(n_tests = n(), raw_p_lt_05 = sum(p_value < 0.05, na.rm = TRUE),
                  BH_FDR_lt_05 = sum(p_adj_BH < 0.05, na.rm = TRUE), .groups = "drop"))

cat("\nSaved:\n")
cat("  channel_diff_within_condition_wholebrain.csv\n")
cat("  channel_diff_within_condition_division.csv\n")
cat("  (region-level: condition_channel_diff_by_region.csv, already existed)\n")


# ---------------------------------------------------------------------------
# channel_validity_boxplots
# ---------------------------------------------------------------------------
# Boxplots only: whole-brain intensity by channel (Cre vs tTA), WITHIN each
# testable condition (unpaired Wilcoxon) -- extracted from channel_validity_
# visualization.R panel B.

library(dplyr)
library(ggplot2)

setwd("~/Desktop/Bachelor")
source("prepare_heatmap_data.R")
long <- long %>% filter(region_acronym != "FS")

testable_conditions <- c("forced swim · Acute", "restraint · Acute", "tail suspension · Acute")
wb_within <- read.csv("channel_diff_within_condition_wholebrain.csv", stringsAsFactors = FALSE)

animal_overall_within <- long %>%
  filter(condition %in% testable_conditions) %>%
  group_by(animal_id_original, channel, condition) %>%
  summarise(overall_median = median(intensity, na.rm = TRUE), .groups = "drop")

p_labels <- wb_within %>% distinct(condition, p_value) %>%
  mutate(label = sprintf("p = %.3f", p_value))

p <- ggplot(animal_overall_within, aes(x = channel, y = overall_median, color = channel)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3, width = 0.5) +
  geom_jitter(width = 0.12, size = 9.6, alpha = 0.8) +
  geom_text(data = p_labels, aes(x = 1.5, y = max(animal_overall_within$overall_median) * 1.05, label = label),
            inherit.aes = FALSE, size = 21) +
  scale_color_manual(values = c("Cre" = "#C2185B", "tTA" = "#4A148C"), guide = "none") +
  facet_wrap(~condition, ncol = 3) +
  labs(title = "Whole-brain intensity by channel, within each testable condition",
       subtitle = "Unpaired Wilcoxon rank-sum test per condition",
       x = NULL, y = "Whole-brain median intensity") +
  theme_minimal(base_size = 72) +
  theme(plot.title = element_text(face = "bold", size = 78),
        plot.subtitle = element_text(size = 60),
        axis.title = element_text(size = 66),
        axis.text = element_text(size = 57),
        strip.text = element_text(face = "bold", size = 60))

ggsave("channel_validity_boxplots.png", p, width = 54, height = 27, units = "in", dpi = 150, limitsize = FALSE)
cat("Saved: channel_validity_boxplots.png\n")


# ---------------------------------------------------------------------------
# channel_validity_visualization
# ---------------------------------------------------------------------------
# Visualizes the channel validity check (thesis section 3.4): is there a
# systematic Cre-vs-tTA reporter bias that could confound the condition/sex
# results? Combines:
#   A. Whole-brain scatter, Cre vs tTA, paired by animal, condition-
#      independent (channel_bias_check.R design) -- all 27 animals.
#   B. Whole-brain intensity by channel, WITHIN each of the 3 testable
#      conditions (forced swim Acute, restraint Acute, tail suspension
#      Acute -- the only conditions tagged by both channels), unpaired
#      Wilcoxon per condition (channel_diff_within_condition_full.R design).
#   C. Volcano plots (region-level within-condition channel test,
#      condition_channel_diff_by_region.csv), one panel per testable
#      condition.

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(patchwork)

setwd("~/Desktop/Bachelor")

## ============================================================
## A. Whole-brain scatter, Cre vs tTA, paired by animal (condition-independent)
## ============================================================
wb_paired <- read.csv("channel_bias_wholebrain.csv", stringsAsFactors = FALSE)
wt_a <- wilcox.test(wb_paired$Cre, wb_paired$tTA, paired = TRUE)

pA <- ggplot(wb_paired, aes(x = tTA, y = Cre)) +
  geom_point(size = 2.2, color = "#C2185B", alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  labs(title = "A. Whole-brain Cre vs tTA (paired by animal, n=27, condition-independent)",
       subtitle = sprintf("Wilcoxon signed-rank: V=%.0f, p=%.3f | dashed = equality", wt_a$statistic, wt_a$p.value),
       x = "tTA whole-brain median", y = "Cre whole-brain median") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(size = 10, face = "bold"), plot.subtitle = element_text(size = 8))

## ============================================================
## B. Whole-brain intensity by channel, WITHIN each testable condition
## ============================================================
setwd("~/Desktop/Bachelor")
source("prepare_heatmap_data.R")
long <- long %>% filter(region_acronym != "FS")

testable_conditions <- c("forced swim · Acute", "restraint · Acute", "tail suspension · Acute")
wb_within <- read.csv("channel_diff_within_condition_wholebrain.csv", stringsAsFactors = FALSE)

animal_overall_within <- long %>%
  filter(condition %in% testable_conditions) %>%
  group_by(animal_id_original, channel, condition) %>%
  summarise(overall_median = median(intensity, na.rm = TRUE), .groups = "drop")

p_labels <- wb_within %>% distinct(condition, p_value) %>%
  mutate(label = sprintf("p = %.3f", p_value))

pB <- ggplot(animal_overall_within, aes(x = channel, y = overall_median, color = channel)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3, width = 0.5) +
  geom_jitter(width = 0.12, size = 1.8, alpha = 0.8) +
  geom_text(data = p_labels, aes(x = 1.5, y = max(animal_overall_within$overall_median) * 1.05, label = label),
            inherit.aes = FALSE, size = 3) +
  scale_color_manual(values = c("Cre" = "#C2185B", "tTA" = "#4A148C"), guide = "none") +
  facet_wrap(~condition, ncol = 3) +
  labs(title = "B. Whole-brain intensity by channel, WITHIN each testable condition (unpaired Wilcoxon)",
       x = NULL, y = "Whole-brain median intensity") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(size = 10, face = "bold"), strip.text = element_text(face = "bold", size = 8.5))

## ============================================================
## C. Volcano plots (region-level, within-condition), one per condition
## ============================================================
region_within <- read.csv("condition_channel_diff_by_region.csv", stringsAsFactors = FALSE)
pseudo <- 0.001
region_within <- region_within %>%
  mutate(log2FC = log2((median_Cre + pseudo) / (median_tTA + pseudo)),
         neglog10p = -log10(p_value),
         sig_raw = p_value < 0.05,
         sig_BH = p_adj_BH < 0.05)

make_volcano <- function(cond) {
  dd <- region_within %>% filter(condition == cond)
  top_label <- dd %>% filter(sig_raw) %>% arrange(p_value) %>% slice_head(n = 6)
  ggplot(dd, aes(x = log2FC, y = neglog10p)) +
    geom_point(aes(color = sig_raw), alpha = 0.6, size = 1.3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.35) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey60", linewidth = 0.3) +
    geom_text_repel(data = top_label, aes(label = region_acronym), size = 2.2, max.overlaps = 15) +
    scale_color_manual(values = c("TRUE" = "#C2185B", "FALSE" = "grey75"), guide = "none") +
    labs(title = cond, x = "log2FC (Cre / tTA)", y = expression(-log[10](italic(p)))) +
    theme_minimal(base_size = 9) +
    theme(plot.title = element_text(size = 9, face = "bold"))
}

pC <- wrap_plots(lapply(testable_conditions, make_volcano), ncol = 3) +
  plot_annotation(title = "C. Region-level channel difference (Cre vs tTA), within each testable condition",
                   theme = theme(plot.title = element_text(size = 10, face = "bold")))

## ============================================================
## Combined figure
## ============================================================
combined <- pA / pB / pC + plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(title = "Channel validity check: Cre vs tTA (thesis section 3.4)",
                   theme = theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5)))

ggsave("channel_validity_combined.png", combined, width = 11, height = 13, units = "in", dpi = 300)

cat("Saved: channel_validity_combined.png\n")
cat("\nSummary:\n")
cat(sprintf("A. Whole-brain paired (n=27): V=%.0f, p=%.3f\n", wt_a$statistic, wt_a$p.value))
print(wb_within %>% distinct(condition, p_value))
cat(sprintf("C. Region-level: %d tests, %d raw p<0.05, %d BH-FDR p<0.05\n",
            nrow(region_within), sum(region_within$sig_raw), sum(region_within$sig_BH)))


# ---------------------------------------------------------------------------
# plot_cre_tta_correlation
# ---------------------------------------------------------------------------
# Cre vs. tTA correlation per group: one panel per group (TRDT1, DT1...DT8), matching
# the MCC_style_groups.pdf panel layout. Each panel is a 2D density ("heatmap")
# scatter of median Cre vs. median tTA intensity across regions, with a fitted
# linear regression line and Pearson correlation annotated.

library(dplyr)
library(ggplot2)

## --- Load data ---
df <- read.csv("~/Desktop/ALL_DT.csv", stringsAsFactors = FALSE)

group_order <- c("TRDT1", paste0("DT", 1:8))   # all 9 groups
df <- df %>% filter(animal_group %in% group_order)

## --- Per-region median per group, per channel (same as MCC_style_groups.R) ---
agg <- df %>%
  group_by(region_acronym, animal_group) %>%
  summarise(
    Cre = median(trapcre_intensity, na.rm = TRUE),
    tTA = median(traptta_intensity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(is.finite(Cre), is.finite(tTA)) %>%
  mutate(animal_group = factor(animal_group, levels = group_order))

## --- Log10 transform (data is heavily right-skewed) ---
offset <- 0.001
agg <- agg %>% mutate(log_Cre = log10(Cre + offset), log_tTA = log10(tTA + offset))

## --- Pearson correlation per group (on log-transformed values) ---
cor_labels <- agg %>%
  group_by(animal_group) %>%
  summarise(
    r = cor(log_Cre, log_tTA, method = "pearson"),
    p = cor.test(log_Cre, log_tTA, method = "pearson")$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    label = sprintf("r = %.2f%s", r, ifelse(p < 0.001, " ***", ifelse(p < 0.01, " **", ifelse(p < 0.05, " *", "")))),
    log_Cre = min(agg$log_Cre), log_tTA = max(agg$log_tTA)
  )

## --- Plot: 2D density ("heatmap") scatter + regression line, per group ---
p <- ggplot(agg, aes(x = log_Cre, y = log_tTA)) +
  geom_hex(bins = 30) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.6) +
  geom_text(data = cor_labels, aes(label = label), hjust = 0, vjust = 1, size = 3, fontface = "bold") +
  facet_wrap(~animal_group, ncol = 2) +
  scale_fill_viridis_c(name = "n regions", option = "viridis") +
  labs(
    x = expression(log[10]~"(median Cre intensity)"),
    y = expression(log[10]~"(median tTA intensity)"),
    title = "Cre vs. tTA correlation per group (each point = 1 brain region)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(size = 11, face = "bold")
  )

ggsave("~/Desktop/Bachelor/cre_tta_correlation.pdf", p, width = 8.27, height = 11.69, units = "in")

## --- Print correlation summary to console ---
print(as.data.frame(cor_labels %>% select(animal_group, r, p)))

