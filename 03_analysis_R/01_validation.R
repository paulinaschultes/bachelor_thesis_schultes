# =============================================================================
# Validation of TractQuant against manual CTCF quantification
#
# Consolidated from the following original scripts:
#   - validation_analysis.R
#   - validation_analysis_Cre.R
#   - validation_multi_rater.R
#   - validation_figure5_large_font.R
#   - validation_figure6_large_font.R
# =============================================================================


# ---------------------------------------------------------------------------
# validation_analysis
# ---------------------------------------------------------------------------
# Program validation: does the automated pipeline (traptta_fraction) agree
# with Paulina's manual intensity counts? Method-comparison analysis:
# scatterplot + identity/regression line, paired t-test, Pearson r, R²,
# bias and SD (+ Bland-Altman limits of agreement).

library(readxl)
library(dplyr)
library(ggplot2)
library(patchwork)

## --- Load & clean ---
## Note: the "Intensity count Paulina" column mixes decimal-comma (European,
## e.g. "0,137") and decimal-point (e.g. "3.79E-3") formatting -- both are
## normalized here before converting to numeric.
path <- "~/Desktop/validation1.xlsx"
df <- read_excel(path, sheet = "Tabelle1")
names(df) <- c("animal", "region", "program", "manual")
df <- df %>%
  mutate(
    program = as.numeric(program),
    manual  = as.numeric(gsub(",", ".", manual))
  )

stopifnot(sum(is.na(df$program)) == 0, sum(is.na(df$manual)) == 0)

## --- Stats ---
n <- nrow(df)
pearson <- cor.test(df$program, df$manual, method = "pearson")
r  <- unname(pearson$estimate)
r2 <- r^2
lm_fit <- lm(manual ~ program, data = df)
r2_lm <- summary(lm_fit)$r.squared

paired_t <- t.test(df$program, df$manual, paired = TRUE)

diff <- df$program - df$manual        # bias convention: program - manual
bias <- mean(diff)
sd_diff <- sd(diff)
loa_lower <- bias - 1.96 * sd_diff
loa_upper <- bias + 1.96 * sd_diff

## --- Colors: pink-to-violet qualitative, matching the rest of the thesis ---
animal_cols <- setNames(
  c("#C2185B", "#E05780", "#F48FB1", "#8E24AA", "#6A4C93", "#4A148C"),
  unique(df$animal)
)

## --- Panel 1: scatterplot, program vs. manual ---
lims <- range(c(df$program, df$manual))
p1 <- ggplot(df, aes(x = program, y = manual, color = animal)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black", linewidth = 0.6) +
  geom_point(size = 1.8, alpha = 0.85) +
  scale_color_manual(values = animal_cols, name = "Animal") +
  coord_equal(xlim = lims, ylim = lims) +
  annotate("text", x = lims[1], y = lims[2],
           label = sprintf("r = %.3f\nR² = %.3f\nn = %d", r, r2, n),
           hjust = 0, vjust = 1, size = 3.3, fontface = "bold") +
  labs(x = "Program (traptta_fraction)", y = "Manual count (Paulina)",
       title = "Program vs. manual quantification") +
  theme_minimal(base_size = 10)

## --- Panel 2: Bland-Altman (mean vs. difference) ---
ba_df <- df %>% mutate(mean_val = (program + manual) / 2, diff_val = program - manual)

p2 <- ggplot(ba_df, aes(x = mean_val, y = diff_val, color = animal)) +
  geom_hline(yintercept = bias, color = "#C2185B", linewidth = 0.7) +
  geom_hline(yintercept = loa_upper, color = "#8E24AA", linetype = "dashed") +
  geom_hline(yintercept = loa_lower, color = "#8E24AA", linetype = "dashed") +
  geom_point(size = 1.8, alpha = 0.85) +
  scale_color_manual(values = animal_cols, name = "Animal") +
  annotate("text", x = -Inf, y = bias, label = sprintf("bias = %.3f", bias),
           hjust = -0.05, vjust = -0.5, size = 3, color = "#C2185B", fontface = "bold") +
  annotate("text", x = -Inf, y = loa_upper, label = sprintf("+1.96 SD = %.3f", loa_upper),
           hjust = -0.05, vjust = -0.5, size = 2.8, color = "#8E24AA") +
  annotate("text", x = -Inf, y = loa_lower, label = sprintf("-1.96 SD = %.3f", loa_lower),
           hjust = -0.05, vjust = 1.3, size = 2.8, color = "#8E24AA") +
  labs(x = "Mean of program & manual", y = "Difference (program - manual)",
       title = "Bland-Altman: bias & limits of agreement") +
  theme_minimal(base_size = 10)

combined <- (p1 + p2) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("~/Desktop/Bachelor/validation_scatter_blandaltman.pdf", combined,
       width = 11, height = 6, units = "in")

## --- Console summary ---
cat("\n================ Validation summary ================\n")
cat(sprintf("n paired observations: %d\n", n))
cat(sprintf("Pearson r: %.4f  (95%% CI %.4f - %.4f), p = %.3g\n",
            r, pearson$conf.int[1], pearson$conf.int[2], pearson$p.value))
cat(sprintf("R^2 (r^2): %.4f   R^2 (lm manual~program): %.4f\n", r2, r2_lm))
cat(sprintf("Paired t-test: t(%d) = %.3f, p = %.3g\n", paired_t$parameter, paired_t$statistic, paired_t$p.value))
cat(sprintf("  mean difference (program - manual) = %.4f, 95%% CI %.4f - %.4f\n",
            paired_t$estimate, paired_t$conf.int[1], paired_t$conf.int[2]))
cat(sprintf("Bias (mean program - manual): %.4f\n", bias))
cat(sprintf("SD of differences: %.4f\n", sd_diff))
cat(sprintf("95%% limits of agreement: %.4f to %.4f\n", loa_lower, loa_upper))
cat("======================================================\n")


# ---------------------------------------------------------------------------
# validation_analysis_Cre
# ---------------------------------------------------------------------------
# Program validation, Cre channel: does the automated pipeline (trapcre_fraction)
# agree with Paulina's manual intensity counts? Method-comparison analysis:
# scatterplot + identity/regression line, paired t-test, Pearson r, R²,
# bias and SD (+ Bland-Altman limits of agreement).
# Same analysis/style as validation_analysis.R (tTA channel), applied to
# Channel_Cre.xlsx: 6 animals x 15 regions = 90 paired observations.
# Sarah/Christina have no entries for this channel, so no multi-rater panel.

library(readxl)
library(dplyr)
library(ggplot2)
library(patchwork)

## --- Load & clean ---
## Data table starts at row 6 (header) / row 7 (first data row); only the
## animal/region/program/manual columns (D:G) are used.
path <- "~/Desktop/Channel_Cre.xlsx"
raw <- read_excel(path, sheet = "Tabelle1", skip = 5)

df <- raw %>%
  transmute(
    animal  = `Animal number`,
    region  = `Brain Region`,
    program = as.numeric(trapcre_fraction),
    manual  = as.numeric(`Intensity count Paulina`)
  ) %>%
  filter(!is.na(animal), !is.na(program), !is.na(manual))

stopifnot(nrow(df) == 90)

## --- Stats ---
n <- nrow(df)
pearson <- cor.test(df$program, df$manual, method = "pearson")
r  <- unname(pearson$estimate)
r2 <- r^2
lm_fit <- lm(manual ~ program, data = df)
r2_lm <- summary(lm_fit)$r.squared

paired_t <- t.test(df$program, df$manual, paired = TRUE)

diff <- df$program - df$manual        # bias convention: program - manual
bias <- mean(diff)
sd_diff <- sd(diff)
loa_lower <- bias - 1.96 * sd_diff
loa_upper <- bias + 1.96 * sd_diff

## --- Colors: pink-to-violet qualitative, matching the rest of the thesis ---
animal_cols <- setNames(
  c("#C2185B", "#E05780", "#F48FB1", "#8E24AA", "#6A4C93", "#4A148C"),
  unique(df$animal)
)

## --- Panel 1: scatterplot, program vs. manual ---
lims <- range(c(df$program, df$manual))
p1 <- ggplot(df, aes(x = program, y = manual, color = animal)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black", linewidth = 0.6) +
  geom_point(size = 1.8, alpha = 0.85) +
  scale_color_manual(values = animal_cols, name = "Animal") +
  coord_equal(xlim = lims, ylim = lims) +
  annotate("text", x = lims[1], y = lims[2],
           label = sprintf("r = %.3f\nR² = %.3f\nn = %d", r, r2, n),
           hjust = 0, vjust = 1, size = 3.3, fontface = "bold") +
  labs(x = "Program (trapcre_fraction)", y = "Manual count (Paulina)",
       title = "Program vs. manual quantification (Cre channel)") +
  theme_minimal(base_size = 10)

## --- Panel 2: Bland-Altman (mean vs. difference) ---
ba_df <- df %>% mutate(mean_val = (program + manual) / 2, diff_val = program - manual)

p2 <- ggplot(ba_df, aes(x = mean_val, y = diff_val, color = animal)) +
  geom_hline(yintercept = bias, color = "#C2185B", linewidth = 0.7) +
  geom_hline(yintercept = loa_upper, color = "#8E24AA", linetype = "dashed") +
  geom_hline(yintercept = loa_lower, color = "#8E24AA", linetype = "dashed") +
  geom_point(size = 1.8, alpha = 0.85) +
  scale_color_manual(values = animal_cols, name = "Animal") +
  annotate("text", x = -Inf, y = bias, label = sprintf("bias = %.3f", bias),
           hjust = -0.05, vjust = -0.5, size = 3, color = "#C2185B", fontface = "bold") +
  annotate("text", x = -Inf, y = loa_upper, label = sprintf("+1.96 SD = %.3f", loa_upper),
           hjust = -0.05, vjust = -0.5, size = 2.8, color = "#8E24AA") +
  annotate("text", x = -Inf, y = loa_lower, label = sprintf("-1.96 SD = %.3f", loa_lower),
           hjust = -0.05, vjust = 1.3, size = 2.8, color = "#8E24AA") +
  labs(x = "Mean of program & manual", y = "Difference (program - manual)",
       title = "Bland-Altman: bias & limits of agreement (Cre channel)") +
  theme_minimal(base_size = 10)

combined <- (p1 + p2) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("~/Desktop/Bachelor/validation_scatter_blandaltman_Cre.pdf", combined,
       width = 11, height = 6, units = "in")
ggsave("~/Desktop/Bachelor/validation_scatter_blandaltman_Cre.png", combined,
       width = 11, height = 6, units = "in", dpi = 300)

## --- Console summary ---
cat("\n================ Validation summary (Cre channel) ================\n")
cat(sprintf("n paired observations: %d\n", n))
cat(sprintf("Pearson r: %.4f  (95%% CI %.4f - %.4f), p = %.3g\n",
            r, pearson$conf.int[1], pearson$conf.int[2], pearson$p.value))
cat(sprintf("R^2 (r^2): %.4f   R^2 (lm manual~program): %.4f\n", r2, r2_lm))
cat(sprintf("Paired t-test: t(%d) = %.3f, p = %.3g\n", paired_t$parameter, paired_t$statistic, paired_t$p.value))
cat(sprintf("  mean difference (program - manual) = %.4f, 95%% CI %.4f - %.4f\n",
            paired_t$estimate, paired_t$conf.int[1], paired_t$conf.int[2]))
cat(sprintf("Bias (mean program - manual): %.4f\n", bias))
cat(sprintf("SD of differences: %.4f\n", sd_diff))
cat(sprintf("95%% limits of agreement: %.4f to %.4f\n", loa_lower, loa_upper))
cat("=====================================================================\n")


# ---------------------------------------------------------------------------
# validation_multi_rater
# ---------------------------------------------------------------------------
# Multi-rater validation: program vs. three independent manual counters
# (Paulina, Christina, Sarah), across 8 regions. Not every rater counted
# every region (missing values allowed). Follows the same visual style as
# validation_analysis.R (pink-violet palette, ggplot2).
#
#   Panel 1: one line for "program" (reference) across regions, plus
#            colored points for each rater -- log scale y-axis since
#            values span ~3 orders of magnitude (0.005 to 2.6).
#   Panel 2: deviation from the per-region mean (across whichever
#            sources have a value for that region), as a % deviation --
#            shows which rater/program tends to over- or under-count
#            relative to the group consensus, per region.

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

df <- read.csv("~/Desktop/validation_multi_rater.csv", stringsAsFactors = FALSE)
region_order <- df$region  # keep original row order

long <- df %>%
  pivot_longer(cols = -region, names_to = "source", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(region = factor(region, levels = region_order),
         source = factor(source, levels = c("program", "Paulina", "Christina", "Sarah")))

source_cols <- c("program" = "grey20", "Paulina" = "#C2185B", "Christina" = "#8E24AA", "Sarah" = "#4A148C")

## --- Panel 1: program line + colored rater points, log scale ---
program_line <- long %>% filter(source == "program")

p1 <- ggplot(long, aes(x = region, y = value, color = source, group = source)) +
  geom_line(data = program_line, aes(group = 1), color = "grey20", linewidth = 0.8) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = source_cols, name = NULL) +
  scale_y_log10() +
  labs(title = "Program vs. three manual raters, per region",
       subtitle = "Line = program (reference) | points = program & each rater | log scale",
       x = NULL, y = "Value (log scale)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 0), legend.position = "bottom")

## --- Panel 2: % deviation from the per-region mean (across all
## available sources for that region) ---
dev_df <- long %>%
  group_by(region) %>%
  mutate(region_mean = mean(value), pct_dev = (value - region_mean) / region_mean * 100) %>%
  ungroup()

p2 <- ggplot(dev_df, aes(x = region, y = pct_dev, fill = source)) +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.5) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  scale_fill_manual(values = source_cols, name = NULL) +
  labs(title = "Deviation from per-region mean",
       subtitle = "Mean = average of all available sources (program + raters) for that region",
       x = NULL, y = "Deviation from mean (%)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

combined <- (p1 / p2) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave("~/Desktop/validation_multi_rater.png", combined, width = 9, height = 9, units = "in", dpi = 300)
ggsave("~/Desktop/validation_multi_rater.pdf", combined, width = 9, height = 9, units = "in")

## --- Console summary: deviation table + simple agreement stats ---
write.csv(dev_df %>% select(region, source, value, region_mean, pct_dev) %>% arrange(region, source),
          "~/Desktop/validation_multi_rater_deviations.csv", row.names = FALSE)

cat("=== Deviation from per-region mean (%) ===\n")
print(as.data.frame(dev_df %>% select(region, source, value, region_mean, pct_dev) %>%
                       mutate(across(where(is.numeric), ~round(., 3)))))

cat("\n=== Mean absolute %% deviation per source (lower = closer to group consensus) ===\n")
summary_tbl <- dev_df %>% group_by(source) %>%
  summarise(n_regions = n(), mean_abs_pct_dev = mean(abs(pct_dev)), .groups = "drop") %>%
  arrange(mean_abs_pct_dev)
print(as.data.frame(summary_tbl))

cat("\nSaved:\n  ~/Desktop/validation_multi_rater.png/.pdf\n  ~/Desktop/validation_multi_rater_deviations.csv\n")


# ---------------------------------------------------------------------------
# validation_figure5_large_font
# ---------------------------------------------------------------------------
# Rebuilds Figure 5 (TractQuant vs. manual quantification validation
# scatterplots) with font size set to 12pt throughout -- addresses the
# reviewer comment "Increase the font, too small". Source data verified to
# reproduce the exact reported statistics (Cre: rho=0.938, p=2.65e-42; tTA:
# rho=0.888, p=1.66e-31; n=90 each) via Spearman rank correlation, matching
# the text.
#
# Data sources:
#   Cre: ~/Desktop/Channel_Cre.xlsx (skip 5 header rows)
#   tTA: ~/Desktop/validation1.csv (semicolon-delimited)

library(readxl)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)

setwd("~/Desktop")

## --- Load Cre channel ---
raw_cre <- read_excel("Validation/Channel_Cre.xlsx", sheet = "Tabelle1", skip = 5)
df_cre <- data.frame(
  animal  = raw_cre[["Animal number"]],
  region  = raw_cre[["Brain Region"]],
  program = as.numeric(raw_cre[["trapcre_fraction"]]),
  manual  = as.numeric(raw_cre[["Intensity count Paulina"]])
) %>% filter(!is.na(animal), !is.na(program), !is.na(manual))
stopifnot(nrow(df_cre) == 90)

## --- Load tTA channel ---
df_tta <- read.csv("validation1.csv", sep = ";", stringsAsFactors = FALSE, header = TRUE)
names(df_tta) <- c("animal", "region", "program", "manual")
df_tta <- df_tta %>%
  mutate(program = as.numeric(program), manual = as.numeric(gsub(",", ".", manual))) %>%
  filter(!is.na(program), !is.na(manual))
stopifnot(nrow(df_tta) == 90)

## --- Stats (Spearman, matching thesis text) ---
ct_cre <- cor.test(df_cre$program, df_cre$manual, method = "spearman")
ct_tta <- cor.test(df_tta$program, df_tta$manual, method = "spearman")
cat(sprintf("Cre: rho=%.3f, p=%.3g, n=%d\n", ct_cre$estimate, ct_cre$p.value, nrow(df_cre)))
cat(sprintf("tTA: rho=%.3f, p=%.3g, n=%d\n", ct_tta$estimate, ct_tta$p.value, nrow(df_tta)))

## --- Plot builder: LOESS fit + CI band, region labels, 12pt fonts ---
FONT_PT <- 12
mm_size <- function(pt) pt / (72.27 / 25.4)  # ggplot geom text 'size' aesthetic is in mm

make_panel <- function(df, rho, pval, n, point_color, ribbon_color, xlab, title) {
  ggplot(df, aes(x = program, y = manual)) +
    geom_smooth(method = "loess", se = TRUE, color = "black", fill = ribbon_color, linewidth = 0.9, alpha = 0.35) +
    geom_point(color = point_color, size = 2.2, alpha = 0.85) +
    ggrepel::geom_text_repel(aes(label = region), size = mm_size(FONT_PT), max.overlaps = 25,
                              segment.size = 0.4, segment.color = "grey60", color = "grey20") +
    annotate("text", x = -Inf, y = Inf, hjust = 0, vjust = 1.3,
             label = sprintf("rho = %.3f\np = %.2e\nn = %d\n6x15 regions", rho, pval, n),
             size = mm_size(FONT_PT), fontface = "bold") +
    labs(title = title, x = xlab, y = "Manual count") +
    theme_minimal(base_size = FONT_PT) +
    theme(
      plot.title = element_text(size = FONT_PT, face = "bold", hjust = 0.5),
      axis.title = element_text(size = FONT_PT),
      axis.text = element_text(size = FONT_PT),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

pA <- make_panel(df_cre, ct_cre$estimate, ct_cre$p.value, nrow(df_cre),
                  point_color = "#2E6F9E", ribbon_color = "#AFCBE3",
                  xlab = "Program (trapcre_fraction)", title = "TRAP-Cre")

pB <- make_panel(df_tta, ct_tta$estimate, ct_tta$p.value, nrow(df_tta),
                  point_color = "#C2185B", ribbon_color = "#F4B8CE",
                  xlab = "Program (traptta_fraction)", title = "TRAP-tTA")

combined <- (pA | pB) +
  plot_annotation(
    title = "Figure 5: TractQuant strongly correlates with manual quantification",
    tag_levels = "A",
    theme = theme(plot.title = element_text(size = FONT_PT + 2, face = "bold", hjust = 0.5))
  )

ggsave("~/Desktop/Bachelor/Figure5_validation_large_font.png", combined, width = 12, height = 6, units = "in", dpi = 300, limitsize = FALSE)
ggsave("~/Desktop/Bachelor/Figure5_validation_large_font.pdf", combined, width = 12, height = 6, units = "in", limitsize = FALSE)

cat("\nSaved: Figure5_validation_large_font.png / .pdf\n")


# ---------------------------------------------------------------------------
# validation_figure6_large_font
# ---------------------------------------------------------------------------
# Rebuilds Figure 6 (Bland-Altman analysis) with font size set to 12pt
# throughout -- same request/style as Figure 5 (validation_figure5_large_font.R).
# Heading intentionally omits "percentile-based" wording. Uses the same
# verified data sources (Cre: Channel_Cre.xlsx; tTA: validation1.csv).

library(readxl)
library(dplyr)
library(ggplot2)
library(patchwork)

setwd("~/Desktop")

## --- Load Cre channel ---
raw_cre <- read_excel("Validation/Channel_Cre.xlsx", sheet = "Tabelle1", skip = 5)
df_cre <- data.frame(
  animal  = raw_cre[["Animal number"]],
  region  = raw_cre[["Brain Region"]],
  program = as.numeric(raw_cre[["trapcre_fraction"]]),
  manual  = as.numeric(raw_cre[["Intensity count Paulina"]])
) %>% filter(!is.na(animal), !is.na(program), !is.na(manual))
stopifnot(nrow(df_cre) == 90)

## --- Load tTA channel ---
df_tta <- read.csv("validation1.csv", sep = ";", stringsAsFactors = FALSE, header = TRUE)
names(df_tta) <- c("animal", "region", "program", "manual")
df_tta <- df_tta %>%
  mutate(program = as.numeric(program), manual = as.numeric(gsub(",", ".", manual))) %>%
  filter(!is.na(program), !is.na(manual))
stopifnot(nrow(df_tta) == 90)

## --- Bland-Altman prep: diff = program - manual, mean = (program+manual)/2,
## percentile-based limits of agreement (2.5th / 97.5th pct of the diffs) ---
prep_ba <- function(df) {
  df %>% mutate(diff = program - manual, mean_val = (program + manual) / 2)
}
ba_cre <- prep_ba(df_cre)
ba_tta <- prep_ba(df_tta)

stats_ba <- function(d) {
  list(bias = median(d$diff), lo = quantile(d$diff, 0.025, names = FALSE),
       hi = quantile(d$diff, 0.975, names = FALSE))
}
s_cre <- stats_ba(ba_cre)
s_tta <- stats_ba(ba_tta)
cat(sprintf("Cre: bias=%.3f, limits=[%.3f, %.3f]\n", s_cre$bias, s_cre$lo, s_cre$hi))
cat(sprintf("tTA: bias=%.3f, limits=[%.3f, %.3f]\n", s_tta$bias, s_tta$lo, s_tta$hi))

## --- Plot builder: 12pt fonts, median-bias + limit-of-agreement lines ---
FONT_PT <- 12
mm_size <- function(pt) pt / (72.27 / 25.4)  # ggplot geom text 'size' aesthetic is in mm

make_ba_panel <- function(d, s, point_color, title) {
  ggplot(d, aes(x = mean_val, y = diff)) +
    geom_point(color = point_color, size = 2.2, alpha = 0.85) +
    geom_hline(yintercept = s$bias, color = "#B71C1C", linewidth = 0.9) +
    geom_hline(yintercept = s$hi, color = "#2A2A2A", linetype = "dashed", linewidth = 0.7) +
    geom_hline(yintercept = s$lo, color = "#2A2A2A", linetype = "dashed", linewidth = 0.7) +
    annotate("text", x = -Inf, y = Inf, hjust = 0, vjust = 1.3,
             label = "n = 90\n6x15 regions", size = mm_size(FONT_PT), fontface = "bold") +
    annotate("text", x = Inf, y = s$bias, hjust = 1.05, vjust = -0.6,
             label = sprintf("median bias = %.3f", s$bias), size = mm_size(FONT_PT), fontface = "bold", color = "#B71C1C") +
    annotate("text", x = Inf, y = s$hi, hjust = 1.05, vjust = -0.5,
             label = sprintf("97.5th pct = %.3f", s$hi), size = mm_size(FONT_PT)) +
    annotate("text", x = Inf, y = s$lo, hjust = 1.05, vjust = 1.4,
             label = sprintf("2.5th pct = %.3f", s$lo), size = mm_size(FONT_PT)) +
    labs(title = title, x = "Mean of program & manual", y = "Difference (program - manual)") +
    theme_minimal(base_size = FONT_PT) +
    theme(
      plot.title = element_text(size = FONT_PT, face = "bold", hjust = 0.5),
      axis.title = element_text(size = FONT_PT),
      axis.text = element_text(size = FONT_PT),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

pA <- make_ba_panel(ba_cre, s_cre, "#2E6F9E", "Bland-Altman: TRAP-Cre")
pB <- make_ba_panel(ba_tta, s_tta, "#C2185B", "Bland-Altman: TRAP-tTA")

combined <- (pA | pB) +
  plot_annotation(
    title = "Figure 6: Bland-Altman analysis of TractQuant indicates a channel-dependent bias",
    tag_levels = "A",
    theme = theme(plot.title = element_text(size = FONT_PT + 2, face = "bold", hjust = 0.5))
  )

ggsave("~/Desktop/Bachelor/Figure6_validation_large_font.png", combined, width = 12, height = 6, units = "in", dpi = 300, limitsize = FALSE)
ggsave("~/Desktop/Bachelor/Figure6_validation_large_font.pdf", combined, width = 12, height = 6, units = "in", limitsize = FALSE)

cat("\nSaved: Figure6_validation_large_font.png / .pdf\n")

