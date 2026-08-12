# =============================================================================
# Post-hoc power analysis and result tables
#
# Consolidated from the following original scripts:
#   - condition_posthoc_power.R
#   - condition_table_pretty.R
#   - plot_combined_groups.R
# =============================================================================


# ---------------------------------------------------------------------------
# condition_posthoc_power
# ---------------------------------------------------------------------------
# Post-hoc power analysis for the condition Kruskal-Wallis tests
# (condition_matrix_stats.R, section 3.2.2.1): quantifies WHY the
# whole-brain, division, and region-level tests found no BH-significant
# condition effect, by (a) computing the achieved power for the actually
# observed effect size at each resolution, and (b) the minimum effect size
# (Cohen's f) that COULD have been detected with 80% power given the actual
# sample sizes -- i.e. what the study was and wasn't sensitive to.
#
# Kruskal-Wallis effect size: eta-squared_H = (H - k + 1) / (n - k)
# (standard rank-based analogue of eta-squared for KW; H = chi-sq
# statistic, k = number of groups, n = total sample size), converted to
# Cohen's f = sqrt(eta2 / (1 - eta2)) for use with pwr::pwr.anova.test
# (KW and one-way ANOVA have ~95.5% asymptotic relative efficiency, so the
# ANOVA-based power approximation is standard practice for KW power).

library(dplyr)
library(pwr)

setwd("~/Desktop/Bachelor")

eta2_H <- function(H, k, n) (H - k + 1) / (n - k)
cohens_f <- function(eta2) sqrt(pmax(eta2, 0) / (1 - pmax(eta2, 0)))

k <- 5   # 5 conditions

## ============================================================
## 1. Whole-brain level
## ============================================================
## Observed: chi-sq = 7.01, n = 54 (group sizes 13, 9, 18, 6, 8)
group_n_wb <- c(13, 9, 18, 6, 8)
n_wb <- sum(group_n_wb)
n_harm_wb <- k / sum(1 / group_n_wb)   # harmonic mean group size, for unequal n

H_wb <- 7.012
eta2_wb <- eta2_H(H_wb, k, n_wb)
f_wb <- cohens_f(eta2_wb)

power_wb <- pwr.anova.test(k = k, n = n_harm_wb, f = f_wb, sig.level = 0.05)$power
mdes_wb  <- pwr.anova.test(k = k, n = n_harm_wb, power = 0.80, sig.level = 0.05)$f

cat("=== 1. Whole-brain level ===\n")
cat(sprintf("Observed: chi-sq = %.3f, n = %d (harmonic mean per group = %.2f)\n", H_wb, n_wb, n_harm_wb))
cat(sprintf("Observed effect size: eta2_H = %.4f -> Cohen's f = %.3f\n", eta2_wb, f_wb))
cat(sprintf("Achieved power (alpha=.05) for this observed effect: %.3f\n", power_wb))
cat(sprintf("Minimum effect size (Cohen's f) detectable at 80%% power with n=%.2f/group: f = %.3f\n\n", n_harm_wb, mdes_wb))

## ============================================================
## 2. Division level (14 tests, same n = 54 per test; report using the
##    strongest observed division effect, Olfactory areas, as the
##    best-case scenario, plus the MDES which applies to all 14)
## ============================================================
div <- read.csv("condition_kruskal_by_division.csv")
div_best <- div[which.min(div$p_value), ]

eta2_div <- eta2_H(div_best$chi_sq, k, n_wb)
f_div <- cohens_f(eta2_div)
power_div <- pwr.anova.test(k = k, n = n_harm_wb, f = f_div, sig.level = 0.05)$power
mdes_div <- mdes_wb   # same n as whole-brain -> same MDES

cat("=== 2. Division level ===\n")
cat(sprintf("Strongest observed division effect: %s, chi-sq = %.3f, p = %.3f\n",
            div_best$major_division, div_best$chi_sq, div_best$p_value))
cat(sprintf("Observed effect size: eta2_H = %.4f -> Cohen's f = %.3f\n", eta2_div, f_div))
cat(sprintf("Achieved power (alpha=.05) for this best-case observed effect: %.3f\n", power_div))
cat(sprintf("Minimum effect size (Cohen's f) detectable at 80%% power: f = %.3f (same n as whole-brain)\n\n", mdes_div))

## ============================================================
## 3. Region level (669 tests, n varies 42-54 due to missing values;
##    report using the strongest observed region, DMH, plus the MDES
##    range across the observed n)
## ============================================================
reg <- read.csv("condition_kruskal_by_region.csv")
reg_best <- reg[which.min(reg$p_value), ]

n_reg_range <- range(reg$n_total, na.rm = TRUE)
n_reg_median <- median(reg$n_total, na.rm = TRUE)

eta2_reg <- eta2_H(reg_best$chi_sq, k, reg_best$n_total)
f_reg <- cohens_f(eta2_reg)
power_reg <- pwr.anova.test(k = k, n = reg_best$n_total / k, f = f_reg, sig.level = 0.05)$power
mdes_reg_median <- pwr.anova.test(k = k, n = n_reg_median / k, power = 0.80, sig.level = 0.05)$f

cat("=== 3. Region level ===\n")
cat(sprintf("n per region test ranges %d-%d (median %.0f)\n", n_reg_range[1], n_reg_range[2], n_reg_median))
cat(sprintf("Strongest observed region effect: %s, chi-sq = %.3f, p = %.4f, n = %d\n",
            reg_best$region_acronym, reg_best$chi_sq, reg_best$p_value, reg_best$n_total))
cat(sprintf("Observed effect size: eta2_H = %.4f -> Cohen's f = %.3f\n", eta2_reg, f_reg))
cat(sprintf("Achieved power (alpha=.05) for this single strongest region: %.3f\n", power_reg))
cat(sprintf("Minimum effect size (Cohen's f) detectable at 80%% power (median n): f = %.3f\n\n", mdes_reg_median))

## ============================================================
## Summary table
## ============================================================
summary_df <- data.frame(
  level = c("Whole-brain", "Division (best case)", "Region (best case, DMH)"),
  n = c(n_wb, n_wb, reg_best$n_total),
  observed_f = round(c(f_wb, f_div, f_reg), 3),
  achieved_power = round(c(power_wb, power_div, power_reg), 3),
  mdes_f_80pct = round(c(mdes_wb, mdes_div, mdes_reg_median), 3)
)
write.csv(summary_df, "condition_posthoc_power_summary.csv", row.names = FALSE)
cat("=== Summary ===\n")
print(summary_df)
cat("\nSaved: condition_posthoc_power_summary.csv\n")


# ---------------------------------------------------------------------------
# condition_table_pretty
# ---------------------------------------------------------------------------
# Nicely formatted, color-coded version of channel_condition_table.csv for
# insertion into the Bachelor thesis (Word). One compact block per stress
# condition (same pink-to-violet palette as the thesis figures), packed into
# two columns so the whole figure stays small (each block sized to its own
# row count -- no padding/empty space).

library(dplyr)
library(tidyr)
library(grid)

df <- read.csv("~/Desktop/ALL_DT.csv", stringsAsFactors = FALSE)
distinct_cond <- df %>% distinct(animal_group, animal_id_original, trapcre_condition, traptta_condition)

long <- distinct_cond %>%
  pivot_longer(cols = c(trapcre_condition, traptta_condition), names_to = "channel_col", values_to = "condition") %>%
  mutate(channel = ifelse(channel_col == "trapcre_condition", "Cre", "tTA"),
         condition = gsub("TRAP-(Cre|tTA) . ", "", condition))

per_group <- long %>%
  group_by(Group = animal_group, Channel = channel, Condition = condition) %>%
  summarise(n_animals = n(), .groups = "drop") %>%
  arrange(Condition, Group, Channel)

## --- Condition colors: same pink-to-violet family as the thesis figures ---
cond_cols <- c(
  "forced swim · Acute"     = "#FBC9DE",
  "forced swim · Chronic"   = "#F48FB1",
  "restraint · Acute"       = "#E05780",
  "social defeat · Acute"   = "#C2185B",
  "tail suspension · Acute" = "#8E24AA"
)
text_cols <- c(
  "forced swim · Acute"     = "#3A3A3A",
  "forced swim · Chronic"   = "#3A3A3A",
  "restraint · Acute"       = "#FFFFFF",
  "social defeat · Acute"   = "#FFFFFF",
  "tail suspension · Acute" = "#FFFFFF"
)

conditions <- names(cond_cols)
blocks <- lapply(conditions, function(cond) per_group %>% filter(Condition == cond) %>% select(Group, Channel, n_animals))
names(blocks) <- conditions

## Row height is fixed in absolute units (npc per "text line"); each block's
## own height = (nrow + 2) rows (title + header + data), so blocks are never
## padded with empty space.
line_h  <- 0.042
gap_h   <- 0.02
col_widths <- c(0.40, 0.32, 0.28)   # Group | Channel | n (fractions of column width)

draw_block <- function(cond, x0, y_top, width) {
  rows <- blocks[[cond]]
  fill <- cond_cols[[cond]]
  txtcol <- text_cols[[cond]]
  col_x <- x0 + c(0, cumsum(col_widths))[1:3] * width

  y <- y_top
  grid.rect(x = x0, y = y, width = width, height = line_h, just = c("left", "top"), gp = gpar(fill = fill, col = "white"))
  grid.text(cond, x = x0 + 0.01, y = y - line_h/2, just = c("left", "center"), gp = gpar(fontsize = 8.5, fontface = "bold", col = txtcol))
  y <- y - line_h

  headers <- c("Group", "Channel", "n")
  for (j in seq_along(headers)) {
    grid.rect(x = col_x[j], y = y, width = col_widths[j] * width, height = line_h, just = c("left", "top"), gp = gpar(fill = "#EDEDED", col = "white"))
    grid.text(headers[j], x = col_x[j] + 0.01, y = y - line_h/2, just = c("left", "center"), gp = gpar(fontsize = 7, fontface = "bold", col = "#333333"))
  }
  y <- y - line_h

  for (i in seq_len(nrow(rows))) {
    r <- rows[i, ]
    shade <- if (i %% 2 == 0) "#FAFAFA" else "white"
    for (j in 1:3) {
      grid.rect(x = col_x[j], y = y, width = col_widths[j] * width, height = line_h, just = c("left", "top"), gp = gpar(fill = shade, col = "#E0E0E0"))
    }
    grid.text(r$Group,   x = col_x[1] + 0.01, y = y - line_h/2, just = c("left", "center"), gp = gpar(fontsize = 7, fontface = "bold"))
    grid.text(r$Channel, x = col_x[2] + 0.01, y = y - line_h/2, just = c("left", "center"), gp = gpar(fontsize = 7))
    grid.text(as.character(r$n_animals), x = col_x[3] + 0.01, y = y - line_h/2, just = c("left", "center"), gp = gpar(fontsize = 7))
    y <- y - line_h
  }

  block_height <- line_h * (nrow(rows) + 2)
  grid.rect(x = x0, y = y_top, width = width, height = block_height, just = c("left", "top"), gp = gpar(fill = NA, col = "#BBBBBB"))
  block_height
}

## --- Greedy column packing: always add the next block to the shorter column ---
n_cols <- 2
col_x0 <- c(0, 0.51)
col_width <- 0.48
col_cursor <- rep(1, n_cols)   # current top y (1 = page top, within the content viewport)

# order blocks largest-first for more balanced packing
block_sizes <- sapply(blocks, nrow) + 2
order_idx <- order(-block_sizes)

draw_all <- function() {
  grid.newpage()
  pushViewport(viewport(x = 0.5, y = 0.5, width = 0.95, height = 0.92))

  grid.text("Condition per group and channel (TRAP-Cre / TRAP-tTA)",
            x = 0, y = 1, just = c("left", "top"), gp = gpar(fontsize = 11, fontface = "bold"))

  pushViewport(viewport(y = 0.90, height = 0.90, just = "top"))
  cursor <- c(1, 1)
  for (idx in order_idx) {
    cond <- conditions[idx]
    cc <- which.max(cursor)   # column with more remaining space (i.e. lower content so far)
    h <- draw_block(cond, col_x0[cc], cursor[cc], col_width)
    cursor[cc] <- cursor[cc] - h - gap_h
  }
  upViewport()
  upViewport()
}

png("~/Desktop/condition_table_pretty.png", width = 6.5, height = 4.3, units = "in", res = 300)
draw_all()
dev.off()

pdf("~/Desktop/condition_table_pretty.pdf", width = 6.5, height = 4.3)
draw_all()
dev.off()


# ---------------------------------------------------------------------------
# plot_combined_groups
# ---------------------------------------------------------------------------
# Combined per-group overview: MCC-style median heatmap (Region x Cre/tTA)
# next to the Cre-vs-tTA correlation scatter, for each of the 9 groups
# (TRDT1, DT1...DT8). Both halves share the same pink-to-violet color scheme
# used in ALL_DT_heatmap.pdf (matches the viral construct diagram, Fig. 1).

library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)

setwd("~/Desktop/Bachelor")

## --- Shared pink-to-violet palette (same as ALL_DT_heatmap.pdf) ---
pink_violet <- c("#FDF0F5", "#FBC9DE", "#F48FB1", "#E05780", "#C2185B", "#8E24AA", "#4A148C")

## --- Load data ---
df <- read.csv("~/Desktop/ALL_DT.csv", stringsAsFactors = FALSE)
structures <- read.csv("~/Desktop/tractquant/allen_mouse_25um_v1.2/structures.csv", stringsAsFactors = FALSE)

group_order   <- c("TRDT1", paste0("DT", 1:8))   # all 9 groups
channel_order <- c("Cre", "tTA")
df <- df %>% filter(animal_group %in% group_order)

## --- Per-region median per group, per channel ---
agg <- df %>%
  group_by(region_id, region_acronym, animal_group) %>%
  summarise(
    Cre = median(trapcre_intensity, na.rm = TRUE),
    tTA = median(traptta_intensity, na.rm = TRUE),
    .groups = "drop"
  )

## --- Major anatomical division per region (Allen CCF ontology) ---
major_ids <- c(
  "315" = "Isocortex", "698" = "Olfactory areas", "1089" = "Hippocampal formation",
  "703" = "Cortical subplate", "477" = "Striatum", "803" = "Pallidum",
  "549" = "Thalamus", "1097" = "Hypothalamus", "313" = "Midbrain",
  "771" = "Pons", "354" = "Medulla", "512" = "Cerebellum",
  "1009" = "fiber tracts", "73" = "ventricular systems"
)
major_order <- unname(major_ids)

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
  ungroup() %>%
  mutate(major_division = factor(major_division, levels = major_order)) %>%
  arrange(major_division, region_acronym)

region_order_fixed <- region_lookup$region_acronym
division_per_row <- region_lookup$major_division

major_cols <- c(
  "Isocortex" = "#5B7C99", "Olfactory areas" = "#8C6D46", "Hippocampal formation" = "#6B8F71",
  "Cortical subplate" = "#9C6B85", "Striatum" = "#4F8F8B", "Pallidum" = "#7A6B9C",
  "Thalamus" = "#B08B5A", "Hypothalamus" = "#A85C5C", "Midbrain" = "#6B6B99",
  "Pons" = "#8FA05E", "Medulla" = "#5E8FA0", "Cerebellum" = "#9C7A5B",
  "fiber tracts" = "#A9A9A9", "ventricular systems" = "#C4B7A6"
)

## --- Heatmap matrices per group (rows = regions, fixed order; cols = Cre, tTA) ---
group_matrices <- lapply(group_order, function(g) {
  sub <- agg %>% filter(animal_group == g)
  m <- matrix(NA_real_, nrow = length(region_order_fixed), ncol = 2,
              dimnames = list(region_order_fixed, channel_order))
  idx <- match(sub$region_acronym, region_order_fixed)
  m[cbind(idx, 1)] <- sub$Cre
  m[cbind(idx, 2)] <- sub$tTA
  m
})
names(group_matrices) <- group_order

heat_vals <- unlist(lapply(group_matrices, as.vector))
heat_qs <- quantile(heat_vals, probs = c(0, 0.10, 0.25, 0.50, 0.75, 0.90, 0.99), na.rm = TRUE)
heat_col_fun <- colorRamp2(as.numeric(heat_qs), pink_violet)

## --- Scatter data: log10(Cre) vs log10(tTA), one point per region, per group ---
offset <- 0.001
scatter_data <- agg %>%
  filter(is.finite(Cre), is.finite(tTA)) %>%
  mutate(
    animal_group = factor(animal_group, levels = group_order),
    log_Cre = log10(Cre + offset),
    log_tTA = log10(tTA + offset)
  )

cor_stats <- scatter_data %>%
  group_by(animal_group) %>%
  summarise(
    r = cor(log_Cre, log_tTA, method = "pearson"),
    p = cor.test(log_Cre, log_tTA, method = "pearson")$p.value,
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("r = %.2f%s", r,
                          ifelse(p < 0.001, " ***", ifelse(p < 0.01, " **", ifelse(p < 0.05, " *", "")))))

make_scatter <- function(g) {
  d <- scatter_data %>% filter(animal_group == g)
  lab <- cor_stats %>% filter(animal_group == g) %>% pull(label)
  ggplot(d, aes(x = log_Cre, y = log_tTA)) +
    geom_hex(bins = 25) +
    geom_smooth(method = "lm", se = TRUE, color = "grey20", linewidth = 0.5) +
    annotate("text", x = -Inf, y = Inf, label = lab, hjust = -0.05, vjust = 1.3,
             size = 2.6, fontface = "bold") +
    scale_fill_gradientn(colours = pink_violet, name = "n") +
    labs(x = expression(log[10]~"Cre"), y = expression(log[10]~"tTA")) +
    theme_minimal(base_size = 7) +
    theme(legend.position = "none", plot.margin = margin(2, 2, 2, 2))
}

## --- Heatmap panel per group (MCC style: no clustering, division color blocks) ---
make_heatmap_panel <- function(g) {
  m <- group_matrices[[g]]
  row_anno <- rowAnnotation(
    div = division_per_row, col = list(div = major_cols),
    show_annotation_name = FALSE, show_legend = FALSE,
    simple_anno_size = unit(1.2, "mm")
  )
  Heatmap(
    m, col = heat_col_fun, name = "median",
    cluster_rows = FALSE, cluster_columns = FALSE,
    row_split = division_per_row, row_gap = unit(0, "mm"), row_title = NULL,
    border = TRUE, rect_gp = gpar(col = NA),
    show_row_names = FALSE,
    show_column_names = TRUE, column_names_gp = gpar(fontsize = 5.5), column_names_rot = 0,
    left_annotation = row_anno,
    show_heatmap_legend = FALSE
  )
}

ht_list <- lapply(group_order, make_heatmap_panel)
names(ht_list) <- group_order
gg_list <- lapply(group_order, make_scatter)
names(gg_list) <- group_order

## --- Draw: 5x2 grid of groups (9 groups, last cell empty), each split into
## [heatmap | scatter], + shared legend ---
draw_all <- function() {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 2, widths = unit(c(5, 1.1), "null"))))

  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  pushViewport(viewport(layout = grid.layout(nrow = ceiling(length(group_order) / 2), ncol = 2)))
  for (i in seq_along(group_order)) {
    g <- group_order[i]
    r <- ceiling(i / 2); cc <- ((i - 1) %% 2) + 1
    pushViewport(viewport(layout.pos.row = r, layout.pos.col = cc))

    grid.text(g, x = 0.5, y = 0.98, just = c("center", "top"), gp = gpar(fontsize = 9, fontface = "bold"))
    pushViewport(viewport(y = 0.46, height = 0.9, layout = grid.layout(nrow = 1, ncol = 2, widths = unit(c(1, 2), "null"))))

    pushViewport(viewport(layout.pos.col = 1))
    draw(ht_list[[g]], newpage = FALSE)
    upViewport()

    pushViewport(viewport(layout.pos.col = 2))
    print(gg_list[[g]], newpage = FALSE, vp = viewport())
    upViewport()

    upViewport()
    upViewport()
  }
  upViewport()
  upViewport()

  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
  color_legend <- Legend(col_fun = heat_col_fun, title = "Median intensity /\nhexbin density",
                          legend_height = unit(3, "cm"),
                          title_gp = gpar(fontsize = 7.5, fontface = "bold"),
                          labels_gp = gpar(fontsize = 6.5))
  division_legend <- Legend(labels = major_order, title = "Major division",
                             legend_gp = gpar(fill = major_cols[major_order]),
                             title_gp = gpar(fontsize = 7.5, fontface = "bold"),
                             labels_gp = gpar(fontsize = 6), ncol = 1)
  pd <- packLegend(color_legend, division_legend, direction = "vertical", gap = unit(4, "mm"))
  draw(pd, x = unit(0.05, "npc"), y = unit(0.5, "npc"), just = "left")
  upViewport()
}

a4_width <- 8.27; a4_height <- 11.69
pdf("~/Desktop/Bachelor/combined_groups_heatmap_scatter.pdf", width = a4_width, height = a4_height)
draw_all()
dev.off()

