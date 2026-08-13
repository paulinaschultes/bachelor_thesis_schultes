# =============================================================================
# Validation of TractQuant against manual CTCF quantification
#
# Consolidated from the following original scripts:
#   - validation_figure5_large_font.R  (Figure 11 - Spearman correlation with LOESS fit, per reporter channel)
#   - validation_figure6_large_font.R  (Figure 12 - Bland-Altman analysis, median bias and 2.5/97.5 percentiles)
# =============================================================================


# ---------------------------------------------------------------------------
# validation_figure5_large_font
# Figure 11 - Spearman correlation with LOESS fit, per reporter channel
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
# Figure 12 - Bland-Altman analysis, median bias and 2.5/97.5 percentiles
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

