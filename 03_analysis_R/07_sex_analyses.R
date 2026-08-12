# =============================================================================
# Sex-stratified comparisons
#
# Consolidated from the following original scripts:
#   - sex_comparison_DT1_TRDT1.R
#   - sex_comparison_DT1_TRDT1_normalized.R
#   - sex_comparison_curated_regions.R
#   - sex_comparison_curated_regions_raw.R
#   - sex_comparison_per_region_fdr.R
#   - sex_comparison_per_region_fdr_raw.R
#   - sex_comparison_top_regions_only.R
#   - sex_kruskal_by_condition.R
#   - sex_kruskal_stats.R
# =============================================================================


# ---------------------------------------------------------------------------
# sex_comparison_DT1_TRDT1
# ---------------------------------------------------------------------------
# Sex comparison: DT1 (female) vs. TRDT1 (male) -- these two groups share the
# IDENTICAL stressor conditions (both: tTA = forced swim acute, Cre =
# restraint acute; see group_condition_corrected.csv), so any systematic
# difference between them isolates a SEX effect, not a stressor confound.
#
# Method: per channel (Cre, tTA separately), treat each brain region as a
# paired observation (DT1's group median vs. TRDT1's group median for that
# region). First test normality of the paired differences (Shapiro-Wilk) to
# justify using the non-parametric Wilcoxon signed-rank test instead of a
# paired t-test. 'FS' excluded (injection region, trivially 1.0).

library(dplyr)
library(tidyr)
library(ggplot2)

df <- read.csv("~/Desktop/ALL_DT.csv", stringsAsFactors = FALSE)

## --- Confirm: DT1 and TRDT1 differ only in sex, same conditions ---
sex_check <- df %>% filter(animal_group %in% c("DT1", "TRDT1")) %>% distinct(animal_group, sex)
cat("Sex per group:\n"); print(as.data.frame(sex_check))

cond_map <- read.csv("~/Desktop/group_condition_corrected.csv", stringsAsFactors = FALSE)
cat("\nCondition per group (should be identical for DT1 and TRDT1):\n")
print(cond_map %>% filter(Group %in% c("DT1", "TRDT1")) %>% arrange(Channel))

## --- Median per region x group x channel ---
group_med <- df %>%
  filter(animal_group %in% c("DT1", "TRDT1"), region_acronym != "FS") %>%
  group_by(region_acronym, division, animal_group) %>%
  summarise(Cre = median(trapcre_intensity, na.rm = TRUE),
            tTA = median(traptta_intensity, na.rm = TRUE), .groups = "drop")

run_channel_test <- function(channel) {
  wide <- group_med %>%
    select(region_acronym, animal_group, value = all_of(channel)) %>%
    pivot_wider(names_from = animal_group, values_from = value) %>%
    filter(is.finite(DT1), is.finite(TRDT1))

  log_dt1   <- log10(wide$DT1 + 0.001)
  log_trdt1 <- log10(wide$TRDT1 + 0.001)
  diff <- log_dt1 - log_trdt1

  cat(sprintf("\n========== Channel: %s (n regions = %d) ==========\n", channel, nrow(wide)))

  ## --- Step 1: normality test on paired differences ---
  shap <- shapiro.test(diff)
  cat(sprintf("Shapiro-Wilk on paired differences (log DT1 - log TRDT1): W = %.4f, p = %.3g\n",
              shap$statistic, shap$p.value))
  cat(ifelse(shap$p.value < 0.05,
             "  -> differences are NOT normally distributed -> Wilcoxon signed-rank is the appropriate test.\n",
             "  -> differences look normally distributed -> a paired t-test would also be valid, Wilcoxon still reported as requested.\n"))

  ## --- Step 2: Wilcoxon matched-pairs signed-rank test ---
  wtest <- wilcox.test(log_dt1, log_trdt1, paired = TRUE)
  cat(sprintf("Wilcoxon signed-rank test: V = %.1f, p = %.3g\n", wtest$statistic, wtest$p.value))
  cat(sprintf("Median difference (DT1 - TRDT1, log10 scale): %.4f | direction: %s\n",
              median(diff), ifelse(median(diff) > 0, "DT1 (female) higher", "TRDT1 (male) higher")))

  list(wide = wide, diff = diff, shapiro = shap, wilcox = wtest, channel = channel)
}

res_cre <- run_channel_test("Cre")
res_tta <- run_channel_test("tTA")

## --- Save combined results table ---
summary_tbl <- bind_rows(
  data.frame(channel = "Cre", shapiro_W = res_cre$shapiro$statistic, shapiro_p = res_cre$shapiro$p.value,
             wilcoxon_V = res_cre$wilcox$statistic, wilcoxon_p = res_cre$wilcox$p.value,
             median_diff_log = median(res_cre$diff), n_regions = nrow(res_cre$wide)),
  data.frame(channel = "tTA", shapiro_W = res_tta$shapiro$statistic, shapiro_p = res_tta$shapiro$p.value,
             wilcoxon_V = res_tta$wilcox$statistic, wilcoxon_p = res_tta$wilcox$p.value,
             median_diff_log = median(res_tta$diff), n_regions = nrow(res_tta$wide))
)
write.csv(summary_tbl, "~/Desktop/sex_comparison_DT1_TRDT1_summary.csv", row.names = FALSE)

## --- Per-region paired differences, for follow-up screening ---
per_region <- bind_rows(
  res_cre$wide %>% mutate(channel = "Cre", log_diff = log10(DT1 + 0.001) - log10(TRDT1 + 0.001)),
  res_tta$wide %>% mutate(channel = "tTA", log_diff = log10(DT1 + 0.001) - log10(TRDT1 + 0.001))
) %>% arrange(channel, desc(abs(log_diff)))
write.csv(per_region, "~/Desktop/sex_comparison_DT1_TRDT1_per_region.csv", row.names = FALSE)

## --- Plots: QQ-plot (normality) + paired scatter, per channel ---
make_plots <- function(res) {
  qq_data <- data.frame(diff = res$diff)
  p_qq <- ggplot(qq_data, aes(sample = diff)) +
    stat_qq(color = "#8E24AA") + stat_qq_line(color = "grey40") +
    labs(title = sprintf("%s: QQ-plot of paired differences", res$channel),
         subtitle = sprintf("Shapiro-Wilk p = %.3g", res$shapiro$p.value)) +
    theme_minimal(base_size = 10)

  p_scatter <- ggplot(res$wide, aes(x = log10(TRDT1 + 0.001), y = log10(DT1 + 0.001))) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(color = "#C2185B", alpha = 0.5, size = 1.5) +
    labs(title = sprintf("%s: DT1 (female) vs. TRDT1 (male), per region", res$channel),
         subtitle = sprintf("Wilcoxon signed-rank p = %.3g", res$wilcox$p.value),
         x = "log10 TRDT1 (male)", y = "log10 DT1 (female)") +
    theme_minimal(base_size = 10)

  list(qq = p_qq, scatter = p_scatter)
}

plots_cre <- make_plots(res_cre)
plots_tta <- make_plots(res_tta)

library(patchwork)
combined <- (plots_cre$qq + plots_cre$scatter) / (plots_tta$qq + plots_tta$scatter)
ggsave("~/Desktop/sex_comparison_DT1_TRDT1.png", combined, width = 9, height = 8, units = "in", dpi = 300)
ggsave("~/Desktop/sex_comparison_DT1_TRDT1.pdf", combined, width = 9, height = 8, units = "in")


# ---------------------------------------------------------------------------
# sex_comparison_DT1_TRDT1_normalized
# ---------------------------------------------------------------------------
# Sex comparison: DT1 (female) vs. TRDT1 (male), NORMALIZED version.
# Same design as sex_comparison_DT1_TRDT1.R, but each sample (group x
# channel) is first divided by its OWN global median (across all regions),
# removing the between-group scaling/batch effect confirmed in the raw
# version (nearly all regions shifted the same direction -- a global, not
# region-specific, pattern). What's left, if anything, is a genuine
# REGION-SPECIFIC sex difference.
#
# DT1 and TRDT1 share identical stressor conditions (Cre = restraint acute,
# tTA = forced swim acute), so this isolates sex, not stressor.
# 'FS' excluded (injection region, trivially 1.0).

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

df <- read.csv("~/Desktop/ALL_DT.csv", stringsAsFactors = FALSE)

## --- Median per region x group x channel ---
group_med <- df %>%
  filter(animal_group %in% c("DT1", "TRDT1"), region_acronym != "FS") %>%
  group_by(region_acronym, division, animal_group) %>%
  summarise(Cre = median(trapcre_intensity, na.rm = TRUE),
            tTA = median(traptta_intensity, na.rm = TRUE), .groups = "drop")

## --- Normalize: divide by each (group, channel) sample's own global median
## across all regions ---
long <- group_med %>%
  pivot_longer(cols = c(Cre, tTA), names_to = "channel", values_to = "value")

overall <- long %>% group_by(animal_group, channel) %>%
  summarise(overall_median = median(value, na.rm = TRUE), .groups = "drop")
cat("Global (overall) median per group x channel (the effect being removed):\n")
print(as.data.frame(overall))

long <- long %>% left_join(overall, by = c("animal_group", "channel")) %>%
  mutate(norm_value = value / overall_median)

run_channel_test <- function(ch) {
  wide <- long %>% filter(channel == ch) %>%
    select(region_acronym, animal_group, norm_value) %>%
    pivot_wider(names_from = animal_group, values_from = norm_value) %>%
    filter(is.finite(DT1), is.finite(TRDT1))

  log_dt1   <- log10(wide$DT1 + 0.001)
  log_trdt1 <- log10(wide$TRDT1 + 0.001)
  diff <- log_dt1 - log_trdt1

  cat(sprintf("\n========== Channel: %s (n regions = %d), NORMALIZED ==========\n", ch, nrow(wide)))

  shap <- shapiro.test(diff)
  cat(sprintf("Shapiro-Wilk on paired differences: W = %.4f, p = %.3g\n", shap$statistic, shap$p.value))
  cat(ifelse(shap$p.value < 0.05,
             "  -> not normal -> Wilcoxon signed-rank appropriate.\n",
             "  -> looks normal -> paired t-test would also be valid; Wilcoxon still reported.\n"))

  wtest <- wilcox.test(log_dt1, log_trdt1, paired = TRUE)
  cat(sprintf("Wilcoxon signed-rank test: V = %.1f, p = %.3g\n", wtest$statistic, wtest$p.value))
  cat(sprintf("Median difference (DT1 - TRDT1, normalized log10 scale): %.4f | direction: %s\n",
              median(diff), ifelse(median(diff) > 0, "DT1 (female) higher", "TRDT1 (male) higher")))

  list(wide = wide, diff = diff, shapiro = shap, wilcox = wtest, channel = ch)
}

res_cre <- run_channel_test("Cre")
res_tta <- run_channel_test("tTA")

summary_tbl <- bind_rows(
  data.frame(channel = "Cre", shapiro_W = res_cre$shapiro$statistic, shapiro_p = res_cre$shapiro$p.value,
             wilcoxon_V = res_cre$wilcox$statistic, wilcoxon_p = res_cre$wilcox$p.value,
             median_diff_log = median(res_cre$diff), n_regions = nrow(res_cre$wide)),
  data.frame(channel = "tTA", shapiro_W = res_tta$shapiro$statistic, shapiro_p = res_tta$shapiro$p.value,
             wilcoxon_V = res_tta$wilcox$statistic, wilcoxon_p = res_tta$wilcox$p.value,
             median_diff_log = median(res_tta$diff), n_regions = nrow(res_tta$wide))
)
write.csv(summary_tbl, "~/Desktop/sex_comparison_DT1_TRDT1_normalized_summary.csv", row.names = FALSE)

## --- Per-region normalized differences (for follow-up screening / FDR if wanted) ---
per_region <- bind_rows(
  res_cre$wide %>% mutate(channel = "Cre", log_diff = log10(DT1 + 0.001) - log10(TRDT1 + 0.001)),
  res_tta$wide %>% mutate(channel = "tTA", log_diff = log10(DT1 + 0.001) - log10(TRDT1 + 0.001))
) %>% arrange(channel, desc(abs(log_diff)))
write.csv(per_region, "~/Desktop/sex_comparison_DT1_TRDT1_normalized_per_region.csv", row.names = FALSE)

## --- Plots ---
make_plots <- function(res) {
  qq_data <- data.frame(diff = res$diff)
  p_qq <- ggplot(qq_data, aes(sample = diff)) +
    stat_qq(color = "#8E24AA") + stat_qq_line(color = "grey40") +
    labs(title = sprintf("%s: QQ-plot of paired differences (normalized)", res$channel),
         subtitle = sprintf("Shapiro-Wilk p = %.3g", res$shapiro$p.value)) +
    theme_minimal(base_size = 10)

  p_scatter <- ggplot(res$wide, aes(x = log10(TRDT1 + 0.001), y = log10(DT1 + 0.001))) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(color = "#C2185B", alpha = 0.5, size = 1.5) +
    labs(title = sprintf("%s: DT1 (female) vs. TRDT1 (male), normalized", res$channel),
         subtitle = sprintf("Wilcoxon signed-rank p = %.3g", res$wilcox$p.value),
         x = "log10 TRDT1 (male, normalized)", y = "log10 DT1 (female, normalized)") +
    theme_minimal(base_size = 10)

  list(qq = p_qq, scatter = p_scatter)
}

plots_cre <- make_plots(res_cre)
plots_tta <- make_plots(res_tta)

combined <- (plots_cre$qq + plots_cre$scatter) / (plots_tta$qq + plots_tta$scatter)
ggsave("~/Desktop/sex_comparison_DT1_TRDT1_normalized.png", combined, width = 9, height = 8, units = "in", dpi = 300)
ggsave("~/Desktop/sex_comparison_DT1_TRDT1_normalized.pdf", combined, width = 9, height = 8, units = "in")


# ---------------------------------------------------------------------------
# sex_comparison_curated_regions
# ---------------------------------------------------------------------------
# DT1 (female) vs. TRDT1 (male): restrict the comparison to two curated
# region sets instead of testing everything:
#   1. "Top regions" -- high median AND low CV across all 9 groups
#      (from top_regions_low_deviation.csv)
#   2. "Most variable regions" -- highest CV across all 9 groups, among
#      regions with a reasonable minimum signal (median_across >= 25th
#      percentile) to avoid near-zero-median regions where CV explodes from
#      noise alone, not real biology.
#
# Each animal's per-region value is normalized by that animal's own
# whole-brain median first (removes individual global scaling, consistent
# with sex_comparison_per_region_fdr.R). Unpaired Wilcoxon rank-sum test per
# region (DT1 n=3 vs. TRDT1 n=4), BH-FDR applied WITHIN this curated set only.
#
# IMPORTANT CAVEAT: with n=3 vs n=4, the smallest possible p-value for this
# test is 2/choose(7,3) = 0.0571 -- already above 0.05. No single region can
# ever reach classical significance at this sample size, regardless of how
# few regions are tested or how they are selected. This script is therefore
# a descriptive/exploratory comparison, not a search for significant hits.

library(dplyr)
library(tidyr)
library(ggplot2)

df <- read.csv("~/Desktop/ALL_DT.csv", stringsAsFactors = FALSE)

## --- Region set 1: top regions (high median, low CV) ---
top_regions_df <- read.csv("~/Desktop/top_regions_low_deviation.csv", stringsAsFactors = FALSE)
top_regions <- unique(top_regions_df$region_acronym)

## --- Region set 2: most variable regions (high CV, excluding near-zero-median noise) ---
group_med_all <- df %>% filter(region_acronym != "FS") %>%
  group_by(region_id, region_acronym, animal_group) %>%
  summarise(Cre = median(trapcre_intensity, na.rm = TRUE),
            tTA = median(traptta_intensity, na.rm = TRUE), .groups = "drop")
complete_regions <- group_med_all %>% count(region_acronym) %>% filter(n == 9) %>% pull(region_acronym)
gm <- group_med_all %>% filter(region_acronym %in% complete_regions)

get_most_variable <- function(ch, n_top = 20) {
  s <- gm %>% group_by(region_acronym) %>%
    summarise(median_across = median(.data[[ch]], na.rm = TRUE),
              sd_across = sd(.data[[ch]], na.rm = TRUE), .groups = "drop") %>%
    mutate(cv = sd_across / median_across) %>%
    filter(is.finite(cv), median_across >= quantile(median_across, 0.25, na.rm = TRUE))
  s %>% arrange(desc(cv)) %>% slice_head(n = n_top) %>% pull(region_acronym)
}
most_variable <- unique(c(get_most_variable("Cre"), get_most_variable("tTA")))

cat("n top regions:", length(top_regions), "\n")
cat("n most-variable regions (signal-filtered):", length(most_variable), "\n")

curated_regions <- unique(c(top_regions, most_variable))
region_set_label <- setNames(
  ifelse(curated_regions %in% top_regions & curated_regions %in% most_variable, "Top & most variable",
         ifelse(curated_regions %in% top_regions, "Top region", "Most variable")),
  curated_regions
)
cat("n unique curated regions total:", length(curated_regions), "\n")

## --- Per-animal data, normalized, restricted to curated regions ---
sub <- df %>% filter(animal_group %in% c("DT1", "TRDT1"), region_acronym %in% curated_regions) %>%
  select(animal_group, animal_id_original, region_acronym, division,
         trapcre_intensity, traptta_intensity) %>%
  pivot_longer(cols = c(trapcre_intensity, traptta_intensity),
               names_to = "channel", values_to = "value") %>%
  mutate(channel = ifelse(channel == "trapcre_intensity", "Cre", "tTA"))

animal_overall <- df %>% filter(animal_group %in% c("DT1", "TRDT1")) %>%
  select(animal_id_original, trapcre_intensity, traptta_intensity) %>%
  pivot_longer(cols = c(trapcre_intensity, traptta_intensity), names_to = "channel", values_to = "value") %>%
  mutate(channel = ifelse(channel == "trapcre_intensity", "Cre", "tTA")) %>%
  group_by(animal_id_original, channel) %>%
  summarise(overall_median = median(value, na.rm = TRUE), .groups = "drop")

sub <- sub %>% left_join(animal_overall, by = c("animal_id_original", "channel")) %>%
  mutate(norm_value = value / overall_median)

## --- Per-region Wilcoxon rank-sum test (unpaired), within curated set only ---
run_region_tests <- function(ch) {
  d <- sub %>% filter(channel == ch)
  regions <- unique(d$region_acronym)
  out <- vector("list", length(regions))
  for (i in seq_along(regions)) {
    reg <- regions[i]
    dr <- d %>% filter(region_acronym == reg)
    dt1_vals   <- dr %>% filter(animal_group == "DT1") %>% pull(norm_value)
    trdt1_vals <- dr %>% filter(animal_group == "TRDT1") %>% pull(norm_value)
    if (length(dt1_vals) < 3 || length(trdt1_vals) < 3) { out[[i]] <- NULL; next }
    wt <- suppressWarnings(wilcox.test(dt1_vals, trdt1_vals))
    out[[i]] <- data.frame(region_acronym = reg, channel = ch,
                            median_DT1 = median(dt1_vals, na.rm = TRUE),
                            median_TRDT1 = median(trdt1_vals, na.rm = TRUE),
                            W = wt$statistic, p_value = wt$p.value)
  }
  bind_rows(out)
}

results <- bind_rows(run_region_tests("Cre"), run_region_tests("tTA")) %>%
  group_by(channel) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(region_set = region_set_label[region_acronym]) %>%
  arrange(channel, p_adj_BH, p_value)

write.csv(results, "~/Desktop/sex_comparison_curated_regions.csv", row.names = FALSE)

for (ch in c("Cre", "tTA")) {
  r <- results %>% filter(channel == ch)
  cat(sprintf("\n========== %s (n curated regions = %d) ==========\n", ch, nrow(r)))
  cat(sprintf("n regions raw p < 0.05: %d | BH-adj p < 0.05: %d | BH-adj p < 0.20: %d\n",
              sum(r$p_value < 0.05, na.rm = TRUE), sum(r$p_adj_BH < 0.05, na.rm = TRUE),
              sum(r$p_adj_BH < 0.20, na.rm = TRUE)))
  print(as.data.frame(head(r %>% select(region_acronym, region_set, median_DT1, median_TRDT1, p_value, p_adj_BH), 15)))
}

## --- Plot: all curated regions, faceted, colored by group ---
cond_cols <- c("DT1" = "#C2185B", "TRDT1" = "#4A148C")

make_plot <- function(ch) {
  d <- sub %>% filter(channel == ch, region_acronym %in%
                         (results %>% filter(channel == ch) %>% pull(region_acronym)))
  ord <- results %>% filter(channel == ch) %>% pull(region_acronym)
  d <- d %>% mutate(region_acronym = factor(region_acronym, levels = ord))

  ggplot(d, aes(x = animal_group, y = norm_value, color = animal_group)) +
    geom_jitter(width = 0.1, size = 1.8, alpha = 0.85) +
    stat_summary(fun = median, geom = "crossbar", width = 0.4, color = "grey20", linewidth = 0.3) +
    facet_wrap(~region_acronym, scales = "free_y", ncol = 5) +
    scale_color_manual(values = cond_cols) +
    labs(title = sprintf("%s: curated regions (top + most variable), DT1 vs TRDT1", ch),
         subtitle = "Normalized per animal. Ordered by BH-adjusted p-value (best first)",
         x = NULL, y = "Normalized intensity", color = "Group") +
    theme_minimal(base_size = 8) +
    theme(strip.text = element_text(face = "bold", size = 7),
          axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

p_cre <- make_plot("Cre")
p_tta <- make_plot("tTA")

ggsave("~/Desktop/sex_comparison_curated_regions_Cre.png", p_cre, width = 13, height = 10, units = "in", dpi = 300)
ggsave("~/Desktop/sex_comparison_curated_regions_tTA.png", p_tta, width = 13, height = 10, units = "in", dpi = 300)


# ---------------------------------------------------------------------------
# sex_comparison_curated_regions_raw
# ---------------------------------------------------------------------------
# RAW version of sex_comparison_curated_regions.R -- no global/batch
# normalization. DT1 (female) vs. TRDT1 (male), restricted to curated
# regions (top regions + most variable regions), raw intensity values,
# unpaired Wilcoxon rank-sum test per region, BH-FDR within the curated set.

library(dplyr)
library(tidyr)
library(ggplot2)

df <- read.csv("~/Desktop/ALL_DT.csv", stringsAsFactors = FALSE)

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

sub <- df %>% filter(animal_group %in% c("DT1", "TRDT1"), region_acronym %in% curated_regions) %>%
  select(animal_group, animal_id_original, region_acronym, division,
         trapcre_intensity, traptta_intensity) %>%
  pivot_longer(cols = c(trapcre_intensity, traptta_intensity),
               names_to = "channel", values_to = "value") %>%
  mutate(channel = ifelse(channel == "trapcre_intensity", "Cre", "tTA"))

run_region_tests <- function(ch) {
  d <- sub %>% filter(channel == ch)
  regions <- unique(d$region_acronym)
  out <- vector("list", length(regions))
  for (i in seq_along(regions)) {
    reg <- regions[i]
    dr <- d %>% filter(region_acronym == reg)
    dt1_vals   <- dr %>% filter(animal_group == "DT1") %>% pull(value)
    trdt1_vals <- dr %>% filter(animal_group == "TRDT1") %>% pull(value)
    if (length(dt1_vals) < 3 || length(trdt1_vals) < 3) { out[[i]] <- NULL; next }
    wt <- suppressWarnings(wilcox.test(dt1_vals, trdt1_vals))
    out[[i]] <- data.frame(region_acronym = reg, channel = ch,
                            median_DT1 = median(dt1_vals, na.rm = TRUE),
                            median_TRDT1 = median(trdt1_vals, na.rm = TRUE),
                            W = wt$statistic, p_value = wt$p.value)
  }
  bind_rows(out)
}

results <- bind_rows(run_region_tests("Cre"), run_region_tests("tTA")) %>%
  group_by(channel) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(region_set = region_set_label[region_acronym]) %>%
  arrange(channel, p_adj_BH, p_value)

write.csv(results, "~/Desktop/sex_comparison_curated_regions_raw.csv", row.names = FALSE)

for (ch in c("Cre", "tTA")) {
  r <- results %>% filter(channel == ch)
  cat(sprintf("\n========== %s RAW (n curated regions = %d) ==========\n", ch, nrow(r)))
  cat(sprintf("n regions raw p < 0.05: %d | BH-adj p < 0.05: %d | BH-adj p < 0.20: %d\n",
              sum(r$p_value < 0.05, na.rm = TRUE), sum(r$p_adj_BH < 0.05, na.rm = TRUE),
              sum(r$p_adj_BH < 0.20, na.rm = TRUE)))
  print(as.data.frame(head(r %>% select(region_acronym, region_set, median_DT1, median_TRDT1, p_value, p_adj_BH), 15)))
}

## --- Plot ---
cond_cols <- c("DT1" = "#C2185B", "TRDT1" = "#4A148C")

make_plot <- function(ch) {
  d <- sub %>% filter(channel == ch, region_acronym %in%
                         (results %>% filter(channel == ch) %>% pull(region_acronym)))
  ord <- results %>% filter(channel == ch) %>% pull(region_acronym)
  d <- d %>% mutate(region_acronym = factor(region_acronym, levels = ord))

  ggplot(d, aes(x = animal_group, y = value, color = animal_group)) +
    geom_jitter(width = 0.1, size = 1.8, alpha = 0.85) +
    stat_summary(fun = median, geom = "crossbar", width = 0.4, color = "grey20", linewidth = 0.3) +
    facet_wrap(~region_acronym, scales = "free_y", ncol = 5) +
    scale_color_manual(values = cond_cols) +
    labs(title = sprintf("%s (RAW, no normalization): curated regions, DT1 vs TRDT1", ch),
         subtitle = "Ordered by BH-adjusted p-value (best first)",
         x = NULL, y = "Raw intensity", color = "Group") +
    theme_minimal(base_size = 8) +
    theme(strip.text = element_text(face = "bold", size = 7),
          axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

p_cre <- make_plot("Cre")
p_tta <- make_plot("tTA")

ggsave("~/Desktop/sex_comparison_curated_regions_raw_Cre.png", p_cre, width = 13, height = 10, units = "in", dpi = 300)
ggsave("~/Desktop/sex_comparison_curated_regions_raw_tTA.png", p_tta, width = 13, height = 10, units = "in", dpi = 300)


# ---------------------------------------------------------------------------
# sex_comparison_per_region_fdr
# ---------------------------------------------------------------------------
# Which specific regions carry the DT1 (female) vs. TRDT1 (male) difference?
# Per-region test: for each region, compare DT1's individual animals
# (n = 3) vs. TRDT1's individual animals (n = 4) -- unpaired Wilcoxon
# rank-sum test (Mann-Whitney), since these are different animals, not
# matched pairs. Each animal's per-region value is first normalized by that
# animal's OWN whole-brain median (removes each individual's global scaling,
# consistent with sex_comparison_DT1_TRDT1_normalized.R). Benjamini-Hochberg
# FDR correction applied across all regions per channel.
#
# 'FS' excluded (injection region, trivially 1.0).

library(dplyr)
library(tidyr)
library(ggplot2)

df <- read.csv("~/Desktop/ALL_DT.csv", stringsAsFactors = FALSE)

sub <- df %>% filter(animal_group %in% c("DT1", "TRDT1"), region_acronym != "FS") %>%
  select(animal_group, animal_id_original, region_acronym, division,
         trapcre_intensity, traptta_intensity) %>%
  pivot_longer(cols = c(trapcre_intensity, traptta_intensity),
               names_to = "channel", values_to = "value") %>%
  mutate(channel = ifelse(channel == "trapcre_intensity", "Cre", "tTA"))

## --- Normalize: each animal's value / that animal's own whole-brain median ---
animal_overall <- sub %>% group_by(animal_id_original, channel) %>%
  summarise(overall_median = median(value, na.rm = TRUE), .groups = "drop")

sub <- sub %>% left_join(animal_overall, by = c("animal_id_original", "channel")) %>%
  mutate(norm_value = value / overall_median)

division_lookup <- sub %>% distinct(region_acronym, division)

## --- Per-region Wilcoxon rank-sum test (unpaired), per channel ---
run_region_tests <- function(ch) {
  d <- sub %>% filter(channel == ch)
  regions <- unique(d$region_acronym)
  out <- vector("list", length(regions))

  for (i in seq_along(regions)) {
    reg <- regions[i]
    dr <- d %>% filter(region_acronym == reg)
    dt1_vals   <- dr %>% filter(animal_group == "DT1") %>% pull(norm_value)
    trdt1_vals <- dr %>% filter(animal_group == "TRDT1") %>% pull(norm_value)
    if (length(dt1_vals) < 3 || length(trdt1_vals) < 3 ||
        all(is.na(dt1_vals)) || all(is.na(trdt1_vals))) { out[[i]] <- NULL; next }

    wt <- suppressWarnings(wilcox.test(dt1_vals, trdt1_vals))
    out[[i]] <- data.frame(
      region_acronym = reg, channel = ch,
      median_DT1 = median(dt1_vals, na.rm = TRUE),
      median_TRDT1 = median(trdt1_vals, na.rm = TRUE),
      W = wt$statistic, p_value = wt$p.value
    )
  }
  bind_rows(out)
}

res_cre <- run_region_tests("Cre")
res_tta <- run_region_tests("tTA")

results <- bind_rows(res_cre, res_tta) %>%
  group_by(channel) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  left_join(division_lookup, by = "region_acronym") %>%
  arrange(channel, p_adj_BH, p_value)

write.csv(results, "~/Desktop/sex_comparison_per_region_fdr.csv", row.names = FALSE)

for (ch in c("Cre", "tTA")) {
  r <- results %>% filter(channel == ch)
  cat(sprintf("\n========== %s (n regions tested = %d) ==========\n", ch, nrow(r)))
  cat(sprintf("n regions raw p < 0.05: %d\n", sum(r$p_value < 0.05, na.rm = TRUE)))
  cat(sprintf("n regions BH-adjusted p < 0.05: %d\n", sum(r$p_adj_BH < 0.05, na.rm = TRUE)))
  cat(sprintf("n regions BH-adjusted p < 0.10: %d\n", sum(r$p_adj_BH < 0.10, na.rm = TRUE)))
  cat(sprintf("n regions BH-adjusted p < 0.20: %d\n", sum(r$p_adj_BH < 0.20, na.rm = TRUE)))
  cat("\nTop 15 by BH-adjusted p-value:\n")
  print(as.data.frame(head(r %>% select(region_acronym, division, median_DT1, median_TRDT1, p_value, p_adj_BH), 15)))
}

## --- Plot: top regions per channel (by raw p, since BH may leave nothing) ---
plot_top <- function(ch, n_top = 9) {
  top_regs <- results %>% filter(channel == ch) %>% slice_head(n = n_top) %>% pull(region_acronym)
  d <- sub %>% filter(channel == ch, region_acronym %in% top_regs) %>%
    mutate(region_acronym = factor(region_acronym, levels = top_regs))

  ggplot(d, aes(x = animal_group, y = norm_value, color = animal_group)) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.85) +
    stat_summary(fun = median, geom = "crossbar", width = 0.4, color = "grey20", linewidth = 0.3) +
    facet_wrap(~region_acronym, scales = "free_y", ncol = 3) +
    scale_color_manual(values = c("DT1" = "#C2185B", "TRDT1" = "#4A148C")) +
    labs(title = sprintf("%s: top regions by (uncorrected) rank-sum p-value", ch),
         subtitle = "DT1 = female (n=3), TRDT1 = male (n=4); normalized per animal",
         x = NULL, y = "Normalized intensity", color = "Group") +
    theme_minimal(base_size = 9) +
    theme(strip.text = element_text(face = "bold"))
}

p_cre <- plot_top("Cre")
p_tta <- plot_top("tTA")

ggsave("~/Desktop/sex_comparison_per_region_fdr_Cre_top.png", p_cre, width = 9, height = 8, units = "in", dpi = 300)
ggsave("~/Desktop/sex_comparison_per_region_fdr_tTA_top.png", p_tta, width = 9, height = 8, units = "in", dpi = 300)


# ---------------------------------------------------------------------------
# sex_comparison_per_region_fdr_raw
# ---------------------------------------------------------------------------
# RAW version of sex_comparison_per_region_fdr.R -- no global/batch
# normalization. Per region: DT1 (female, n=3) vs. TRDT1 (male, n=4),
# unpaired Wilcoxon rank-sum test on raw intensity values, BH-FDR across
# all regions per channel.
#
# 'FS' excluded (injection region, trivially 1.0).

library(dplyr)
library(tidyr)
library(ggplot2)

df <- read.csv("~/Desktop/ALL_DT.csv", stringsAsFactors = FALSE)

sub <- df %>% filter(animal_group %in% c("DT1", "TRDT1"), region_acronym != "FS") %>%
  select(animal_group, animal_id_original, region_acronym, division,
         trapcre_intensity, traptta_intensity) %>%
  pivot_longer(cols = c(trapcre_intensity, traptta_intensity),
               names_to = "channel", values_to = "value") %>%
  mutate(channel = ifelse(channel == "trapcre_intensity", "Cre", "tTA"))

division_lookup <- sub %>% distinct(region_acronym, division)

run_region_tests <- function(ch) {
  d <- sub %>% filter(channel == ch)
  regions <- unique(d$region_acronym)
  out <- vector("list", length(regions))
  for (i in seq_along(regions)) {
    reg <- regions[i]
    dr <- d %>% filter(region_acronym == reg)
    dt1_vals   <- dr %>% filter(animal_group == "DT1") %>% pull(value)
    trdt1_vals <- dr %>% filter(animal_group == "TRDT1") %>% pull(value)
    if (length(dt1_vals) < 3 || length(trdt1_vals) < 3 ||
        all(is.na(dt1_vals)) || all(is.na(trdt1_vals))) { out[[i]] <- NULL; next }
    wt <- suppressWarnings(wilcox.test(dt1_vals, trdt1_vals))
    out[[i]] <- data.frame(region_acronym = reg, channel = ch,
                            median_DT1 = median(dt1_vals, na.rm = TRUE),
                            median_TRDT1 = median(trdt1_vals, na.rm = TRUE),
                            W = wt$statistic, p_value = wt$p.value)
  }
  bind_rows(out)
}

results <- bind_rows(run_region_tests("Cre"), run_region_tests("tTA")) %>%
  group_by(channel) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  left_join(division_lookup, by = "region_acronym") %>%
  arrange(channel, p_adj_BH, p_value)

write.csv(results, "~/Desktop/sex_comparison_per_region_fdr_raw.csv", row.names = FALSE)

for (ch in c("Cre", "tTA")) {
  r <- results %>% filter(channel == ch)
  cat(sprintf("\n========== %s RAW (n regions tested = %d) ==========\n", ch, nrow(r)))
  cat(sprintf("n regions raw p < 0.05: %d\n", sum(r$p_value < 0.05, na.rm = TRUE)))
  cat(sprintf("n regions BH-adjusted p < 0.05: %d\n", sum(r$p_adj_BH < 0.05, na.rm = TRUE)))
  cat(sprintf("n regions BH-adjusted p < 0.20: %d\n", sum(r$p_adj_BH < 0.20, na.rm = TRUE)))
  cat("\nTop 15 by BH-adjusted p-value:\n")
  print(as.data.frame(head(r %>% select(region_acronym, division, median_DT1, median_TRDT1, p_value, p_adj_BH), 15)))
}


# ---------------------------------------------------------------------------
# sex_comparison_top_regions_only
# ---------------------------------------------------------------------------
# Are the TOP regions (high median, low CV across all 9 groups; from
# top_regions_low_deviation.csv) the same between DT1 (female) and TRDT1
# (male)? Per region: unpaired Wilcoxon rank-sum test (DT1 n=3 vs. TRDT1
# n=4), raw intensity (no batch/global normalization, per prior decision),
# BH-FDR applied across just this top-region set (not all 578 regions).
#
# IMPORTANT STATISTICAL CAVEAT: a non-significant Wilcoxon test does NOT
# prove the two groups are "the same" -- it only means no difference was
# detected (absence of evidence != evidence of absence). With n=3 vs n=4,
# the test has very little power to detect real differences even if they
# exist. A true equivalence claim would need a dedicated equivalence test
# (e.g. TOST) with a pre-specified margin, not a standard Wilcoxon test.
# Results below are reported as requested, with this caveat attached.

library(dplyr)
library(tidyr)
library(ggplot2)

df <- read.csv("~/Desktop/ALL_DT.csv", stringsAsFactors = FALSE)

top_regions_df <- read.csv("~/Desktop/top_regions_low_deviation.csv", stringsAsFactors = FALSE)
top_regions <- unique(top_regions_df$region_acronym)
cat("n top regions (union across Cre/tTA lists):", length(top_regions), "\n")

sub <- df %>% filter(animal_group %in% c("DT1", "TRDT1"), region_acronym %in% top_regions) %>%
  select(animal_group, animal_id_original, region_acronym, division,
         trapcre_intensity, traptta_intensity) %>%
  pivot_longer(cols = c(trapcre_intensity, traptta_intensity),
               names_to = "channel", values_to = "value") %>%
  mutate(channel = ifelse(channel == "trapcre_intensity", "Cre", "tTA"))

run_region_tests <- function(ch) {
  d <- sub %>% filter(channel == ch)
  regions <- unique(d$region_acronym)
  out <- vector("list", length(regions))
  for (i in seq_along(regions)) {
    reg <- regions[i]
    dr <- d %>% filter(region_acronym == reg)
    dt1_vals   <- dr %>% filter(animal_group == "DT1") %>% pull(value)
    trdt1_vals <- dr %>% filter(animal_group == "TRDT1") %>% pull(value)
    if (length(dt1_vals) < 3 || length(trdt1_vals) < 3) { out[[i]] <- NULL; next }
    wt <- suppressWarnings(wilcox.test(dt1_vals, trdt1_vals))
    out[[i]] <- data.frame(region_acronym = reg, channel = ch,
                            median_DT1 = median(dt1_vals, na.rm = TRUE),
                            median_TRDT1 = median(trdt1_vals, na.rm = TRUE),
                            fold_DT1_over_TRDT1 = median(dt1_vals, na.rm = TRUE) / median(trdt1_vals, na.rm = TRUE),
                            W = wt$statistic, p_value = wt$p.value)
  }
  bind_rows(out)
}

results <- bind_rows(run_region_tests("Cre"), run_region_tests("tTA")) %>%
  group_by(channel) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(verdict = ifelse(p_adj_BH < 0.05, "significantly DIFFERENT",
                           "no significant difference detected (NOT proof of equality)")) %>%
  arrange(channel, p_adj_BH, p_value)

write.csv(results, "~/Desktop/sex_comparison_top_regions_only.csv", row.names = FALSE)

for (ch in c("Cre", "tTA")) {
  r <- results %>% filter(channel == ch)
  cat(sprintf("\n========== %s (n top regions = %d) ==========\n", ch, nrow(r)))
  cat(sprintf("n regions raw p < 0.05: %d\n", sum(r$p_value < 0.05, na.rm = TRUE)))
  cat(sprintf("n regions BH-adjusted p < 0.05 (= significantly different): %d\n", sum(r$p_adj_BH < 0.05, na.rm = TRUE)))
  cat(sprintf("n regions BH-adjusted p >= 0.05 (= no difference detected, NOT proof of sameness): %d\n", sum(r$p_adj_BH >= 0.05, na.rm = TRUE)))
  cat("\nAll regions, sorted by BH-adjusted p-value:\n")
  print(as.data.frame(r %>% select(region_acronym, median_DT1, median_TRDT1, fold_DT1_over_TRDT1, p_value, p_adj_BH)))
}

## --- Plot ---
cond_cols <- c("DT1" = "#C2185B", "TRDT1" = "#4A148C")
make_plot <- function(ch) {
  ord <- results %>% filter(channel == ch) %>% pull(region_acronym)
  d <- sub %>% filter(channel == ch) %>% mutate(region_acronym = factor(region_acronym, levels = ord))
  ggplot(d, aes(x = animal_group, y = value, color = animal_group)) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.85) +
    stat_summary(fun = median, geom = "crossbar", width = 0.4, color = "grey20", linewidth = 0.3) +
    facet_wrap(~region_acronym, scales = "free_y", ncol = 5) +
    scale_color_manual(values = cond_cols) +
    labs(title = sprintf("%s: top regions, DT1 (female) vs TRDT1 (male)", ch),
         subtitle = "Raw values, ordered by BH-adjusted p-value (most different first)",
         x = NULL, y = "Raw intensity", color = "Group") +
    theme_minimal(base_size = 9) +
    theme(strip.text = element_text(face = "bold"), axis.text.x = element_blank(), axis.ticks.x = element_blank())
}
ggsave("~/Desktop/sex_comparison_top_regions_only_Cre.png", make_plot("Cre"), width = 11, height = 7, units = "in", dpi = 300)
ggsave("~/Desktop/sex_comparison_top_regions_only_tTA.png", make_plot("tTA"), width = 11, height = 7, units = "in", dpi = 300)


# ---------------------------------------------------------------------------
# sex_kruskal_by_condition
# ---------------------------------------------------------------------------
# Sex (Male vs Female) effect, tested SEPARATELY WITHIN each condition --
# complements sex_kruskal_stats.R (which pools all conditions together to
# ask "is sex a confound overall?"). Here we ask "within condition X, does
# sex still show up as a difference?", mirroring the same whole-brain /
# division / region hierarchy and BH-FDR scheme, but run once per condition.
#
# "forced swim · Chronic" is excluded: it has 9 Female and 0 Male animals
# (see plot_heatmap_ALL_DT_by_sex.R), so a sex comparison is not possible
# for that condition -- there is no Male group to compare against.
#   1a. whole-brain level (1 test per condition)
#   1b. division level (14 tests per condition)
#   1c. region level (669 tests per condition)
# BH-FDR applied within each condition x resolution combination separately
# (i.e. the division-level p-values for "restraint" are corrected only
# against the other 13 division tests for "restraint", not pooled with
# other conditions' division tests).

library(dplyr)
library(tidyr)

setwd("~/Desktop/Bachelor")
source("prepare_heatmap_data.R")   # mat, long, major_division_per_row, condition_order

keep <- rownames(mat) != "FS"
major_division_per_row <- major_division_per_row[keep]
mat <- mat[keep, , drop = FALSE]
long <- long %>% filter(region_acronym != "FS")
regions <- rownames(mat)

## Only conditions with both sexes represented can be tested.
sex_counts_per_condition <- long %>% distinct(animal_id_original, sex, condition) %>%
  count(condition, sex) %>% tidyr::pivot_wider(names_from = sex, values_from = n, values_fill = 0)
cat("=== Animals per condition x sex ===\n")
print(as.data.frame(sex_counts_per_condition))

testable_conditions <- sex_counts_per_condition %>%
  filter(Female > 0, Male > 0) %>% pull(condition) %>% as.character()
skipped_conditions <- setdiff(as.character(unique(long$condition)), testable_conditions)
cat("\nTestable conditions (both sexes present):", paste(testable_conditions, collapse = ", "), "\n")
cat("Skipped (single-sex, no comparison possible):", paste(skipped_conditions, collapse = ", "), "\n\n")

wholebrain_list <- vector("list", length(testable_conditions))
division_list   <- vector("list", length(testable_conditions))
region_list     <- vector("list", length(testable_conditions))

for (ci in seq_along(testable_conditions)) {
  cond <- testable_conditions[ci]
  long_c <- long %>% filter(condition == cond)

  ## --- 1a. whole-brain ---
  animal_overall_c <- long_c %>%
    group_by(animal_id_original, channel, sex) %>%
    summarise(overall_median = median(intensity, na.rm = TRUE), .groups = "drop")

  kt_wb <- kruskal.test(overall_median ~ sex, data = animal_overall_c)
  medians_wb <- animal_overall_c %>% group_by(sex) %>%
    summarise(n = n(), median_overall = median(overall_median), .groups = "drop")
  wholebrain_list[[ci]] <- medians_wb %>%
    mutate(condition = cond, chi_sq = unname(kt_wb$statistic), df = unname(kt_wb$parameter),
           p_value = kt_wb$p.value)

  ## --- 1b. division-level (whole-brain-median normalized per animal x channel) ---
  long_norm_c <- long_c %>%
    left_join(animal_overall_c, by = c("animal_id_original", "channel", "sex")) %>%
    mutate(norm_value = intensity / overall_median)

  division_scores_c <- long_norm_c %>%
    group_by(animal_id_original, channel, sex,
             major_division = major_division_per_row[match(region_acronym, rownames(mat))]) %>%
    summarise(div_value = median(norm_value, na.rm = TRUE), n_regions = n(), .groups = "drop") %>%
    filter(!is.na(major_division)) %>%
    mutate(log_div_value = log10(div_value + 0.001))

  divisions <- unique(division_scores_c$major_division)
  div_results <- vector("list", length(divisions))
  for (i in seq_along(divisions)) {
    dv <- divisions[i]
    d <- division_scores_c %>% filter(major_division == dv)
    if (n_distinct(d$sex) < 2 || any(table(d$sex) < 2)) next
    kt <- kruskal.test(log_div_value ~ sex, data = d)
    medians <- d %>% group_by(sex) %>% summarise(m = median(log_div_value), .groups = "drop")
    div_results[[i]] <- data.frame(
      condition = cond, major_division = dv, chi_sq = unname(kt$statistic), df = unname(kt$parameter),
      n_regions_pooled = round(mean(d$n_regions)), p_value = kt$p.value,
      t(setNames(medians$m, paste0("median_log_", medians$sex)))
    )
  }
  division_list[[ci]] <- bind_rows(div_results) %>%
    mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

  ## --- 1c. region-level ---
  reg_results <- vector("list", length(regions))
  for (i in seq_along(regions)) {
    reg <- regions[i]
    d <- long_c %>% filter(region_acronym == reg)
    if (n_distinct(d$sex) < 2 || any(table(d$sex) < 2) || nrow(d) < 5) next
    kt <- suppressWarnings(kruskal.test(intensity ~ sex, data = d))
    medians <- d %>% group_by(sex) %>% summarise(m = median(intensity, na.rm = TRUE), .groups = "drop")
    reg_results[[i]] <- data.frame(
      condition = cond, region_acronym = reg, chi_sq = unname(kt$statistic), df = unname(kt$parameter),
      n_total = nrow(d), p_value = kt$p.value,
      t(setNames(medians$m, paste0("median_", medians$sex)))
    )
  }
  region_list[[ci]] <- bind_rows(reg_results) %>%
    mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

  cat(sprintf("=== %s ===\n", cond))
  cat(sprintf("  whole-brain: chi-sq=%.3f, p=%.4f\n", kt_wb$statistic, kt_wb$p.value))
  cat(sprintf("  division:    %d/%d raw p<0.05, %d/%d BH-FDR p<0.05\n",
              sum(division_list[[ci]]$p_value < 0.05, na.rm = TRUE), nrow(division_list[[ci]]),
              sum(division_list[[ci]]$p_adj_BH < 0.05, na.rm = TRUE), nrow(division_list[[ci]])))
  cat(sprintf("  region:      %d/%d raw p<0.05, %d/%d BH-FDR p<0.05\n\n",
              sum(region_list[[ci]]$p_value < 0.05, na.rm = TRUE), nrow(region_list[[ci]]),
              sum(region_list[[ci]]$p_adj_BH < 0.05, na.rm = TRUE), nrow(region_list[[ci]])))
}

wholebrain_df <- bind_rows(wholebrain_list) %>%
  select(condition, sex, n, median_overall, chi_sq, df, p_value) %>%
  arrange(p_value)
division_df <- bind_rows(division_list) %>% arrange(condition, p_value)
region_df   <- bind_rows(region_list) %>% arrange(condition, p_adj_BH, p_value)

write.csv(wholebrain_df, "sex_kruskal_by_condition_wholebrain.csv", row.names = FALSE)
write.csv(division_df,   "sex_kruskal_by_condition_division.csv", row.names = FALSE)
write.csv(region_df,     "sex_kruskal_by_condition_region.csv", row.names = FALSE)

cat("=== Summary: whole-brain sex effect per condition ===\n")
print(as.data.frame(wholebrain_df %>% distinct(condition, chi_sq, df, p_value)))

cat("\nSaved:\n")
cat("  sex_kruskal_by_condition_wholebrain.csv\n")
cat("  sex_kruskal_by_condition_division.csv\n")
cat("  sex_kruskal_by_condition_region.csv\n")


# ---------------------------------------------------------------------------
# sex_kruskal_stats
# ---------------------------------------------------------------------------
# Sex (Male vs Female) effect, mirroring condition_matrix_stats.R's
# whole-brain/division/region Kruskal-Wallis hierarchy exactly, but grouping
# by Sex instead of Condition -- checks whether sex is a confound
# independent of stressor condition, the same way channel_bias_check.R
# checks reporter channel. Sex is binary (14 Male, 13 Female), so
# Kruskal-Wallis here is mathematically equivalent to a Wilcoxon rank-sum
# test -- kept as Kruskal-Wallis purely for consistency of naming/method
# with the condition analysis.
#   1a. whole-brain level (1 test)
#   1b. division level (14 tests)
#   1c. region level (669 tests)
# BH-FDR applied within each resolution separately.

library(dplyr)
library(tidyr)

setwd("~/Desktop/Bachelor")
source("prepare_heatmap_data.R")   # mat, long (now includes sex), major_division_per_row

keep <- rownames(mat) != "FS"
major_division_per_row <- major_division_per_row[keep]
mat <- mat[keep, , drop = FALSE]
long <- long %>% filter(region_acronym != "FS")
regions <- rownames(mat)

## ============================================================
## 1a. Whole-brain-level Kruskal-Wallis (Sex)
## ============================================================
animal_overall_sex <- long %>%
  group_by(animal_id_original, channel, sex) %>%
  summarise(overall_median = median(intensity, na.rm = TRUE), .groups = "drop")

kt_wb_sex <- kruskal.test(overall_median ~ sex, data = animal_overall_sex)
cat("=== 1a. Whole-brain-level Kruskal-Wallis (Sex) ===\n")
cat(sprintf("n animal x channel samples: %d | chi-sq = %.3f, df = %d, p = %.4f\n",
            nrow(animal_overall_sex), kt_wb_sex$statistic, kt_wb_sex$parameter, kt_wb_sex$p.value))
print(animal_overall_sex %>% group_by(sex) %>%
        summarise(n = n(), median_overall = median(overall_median), .groups = "drop"))

write.csv(
  animal_overall_sex %>% group_by(sex) %>%
    summarise(n = n(), median_overall = median(overall_median), .groups = "drop") %>%
    mutate(chi_sq = unname(kt_wb_sex$statistic), df = unname(kt_wb_sex$parameter), p_value = kt_wb_sex$p.value),
  "sex_kruskal_wholebrain.csv", row.names = FALSE
)

## ============================================================
## 1b. Division-level Kruskal-Wallis (Sex)
## ============================================================
long_norm_sex <- long %>%
  left_join(animal_overall_sex, by = c("animal_id_original", "channel", "sex")) %>%
  mutate(norm_value = intensity / overall_median)

division_scores_sex <- long_norm_sex %>%
  group_by(animal_id_original, channel, sex, major_division = major_division_per_row[match(region_acronym, rownames(mat))]) %>%
  summarise(div_value = median(norm_value, na.rm = TRUE), n_regions = n(), .groups = "drop") %>%
  filter(!is.na(major_division)) %>%
  mutate(log_div_value = log10(div_value + 0.001))

divisions <- unique(division_scores_sex$major_division)
division_results_sex <- vector("list", length(divisions))

for (i in seq_along(divisions)) {
  dv <- divisions[i]
  d <- division_scores_sex %>% filter(major_division == dv)
  if (n_distinct(d$sex) < 2) next
  kt <- kruskal.test(log_div_value ~ sex, data = d)
  medians <- d %>% group_by(sex) %>% summarise(m = median(log_div_value), .groups = "drop")
  division_results_sex[[i]] <- data.frame(
    major_division = dv, chi_sq = unname(kt$statistic), df = unname(kt$parameter),
    n_regions_pooled = round(mean(d$n_regions)), p_value = kt$p.value,
    t(setNames(medians$m, paste0("median_log_", medians$sex)))
  )
}

division_kruskal_sex_df <- bind_rows(division_results_sex) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  arrange(p_value)

write.csv(division_kruskal_sex_df, "sex_kruskal_by_division.csv", row.names = FALSE)

cat("\n=== 1b. Division-level Kruskal-Wallis (Sex), n =", nrow(division_kruskal_sex_df), "divisions ===\n")
print(as.data.frame(division_kruskal_sex_df[, c("major_division", "chi_sq", "df", "p_value", "p_adj_BH")]))
cat(sprintf("\nDivisions with raw p < 0.05: %d / %d\n", sum(division_kruskal_sex_df$p_value < 0.05), nrow(division_kruskal_sex_df)))
cat(sprintf("Divisions with BH-FDR p < 0.05: %d / %d\n", sum(division_kruskal_sex_df$p_adj_BH < 0.05), nrow(division_kruskal_sex_df)))

## ============================================================
## 1c. Region-level Kruskal-Wallis (Sex)
## ============================================================
kruskal_results_sex <- vector("list", length(regions))

for (i in seq_along(regions)) {
  reg <- regions[i]
  d <- long %>% filter(region_acronym == reg)
  if (n_distinct(d$sex) < 2 || nrow(d) < 5) next
  kt <- suppressWarnings(kruskal.test(intensity ~ sex, data = d))
  medians <- d %>% group_by(sex) %>% summarise(m = median(intensity, na.rm = TRUE), .groups = "drop")
  kruskal_results_sex[[i]] <- data.frame(
    region_acronym = reg, chi_sq = unname(kt$statistic), df = unname(kt$parameter),
    n_total = nrow(d), p_value = kt$p.value,
    t(setNames(medians$m, paste0("median_", medians$sex)))
  )
}

kruskal_sex_df <- bind_rows(kruskal_results_sex) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj_BH, p_value)

write.csv(kruskal_sex_df, "sex_kruskal_by_region.csv", row.names = FALSE)

cat("\n=== 1c. Region-level Kruskal-Wallis (Sex), per region ===\n")
cat("Total tests:", nrow(kruskal_sex_df), "\n")
cat("Raw p < 0.05:", sum(kruskal_sex_df$p_value < 0.05, na.rm = TRUE), "\n")
cat("BH-FDR p < 0.05:", sum(kruskal_sex_df$p_adj_BH < 0.05, na.rm = TRUE), "\n")
cat("\nTop 15 by BH-adjusted p-value:\n")
print(head(as.data.frame(kruskal_sex_df), 15))

cat("\nSaved:\n")
cat("  sex_kruskal_wholebrain.csv\n")
cat("  sex_kruskal_by_division.csv\n")
cat("  sex_kruskal_by_region.csv\n")

