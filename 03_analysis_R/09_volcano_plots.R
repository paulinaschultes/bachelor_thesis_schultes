# =============================================================================
# THESIS FIGURES PRODUCED BY THIS SCRIPT
#   Figure 18      Pairwise volcano plots show scattered nominal effects
#   Figure 21      Volcano plot of sex-specific candidate regions
#
# Note: figure numbers appearing in the comments further down refer to earlier
# drafts and are NOT the final numbering. The list above is authoritative.
# =============================================================================

# =============================================================================
# Volcano plots for condition and sex comparisons
#
# Consolidated from the following original scripts:
#   - volcano_condition_pairwise.R
#   - volcano_sex_by_condition.R
#   - volcano_sex_region.R
#   - fs_acute_vs_chronic_pooled_volcano.R
# =============================================================================


# ---------------------------------------------------------------------------
# volcano_condition_pairwise
# ---------------------------------------------------------------------------
# Volcano plots per CONDITION PAIR, with INDEPENDENT statistics per pair --
# fixes the previous version, which read from
# condition_pairwise_wilcoxon_by_region.csv (condition_matrix_stats.R
# section 1d), where BH-FDR was applied GLOBALLY across all 10 pairs x 669
# regions = 6690 tests at once. That shared correction makes it unfairly
# hard for any single, specific comparison (e.g. forced swim Acute vs
# Chronic) to reach significance, since it competes against unrelated pairs
# (e.g. social defeat vs tail suspension) for the same FDR budget.
#
# Here, each of the 10 condition pairs gets its own Wilcoxon rank-sum test
# per region AND its own BH-FDR correction across just that pair's regions
# (588, fiber tracts excluded) -- same design as sex_kruskal_by_region.csv /
# the sex volcano (volcano_sex_region.R), just repeated once per condition pair.

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(patchwork)

setwd("~/Desktop/Bachelor")
source("prepare_heatmap_data.R")   # long, condition_order
long <- long %>% filter(region_acronym != "FS")

pseudo <- 0.001
cond_pairs <- t(combn(condition_order, 2))
regions <- unique(long$region_acronym)

pair_results <- vector("list", nrow(cond_pairs))

for (p in seq_len(nrow(cond_pairs))) {
  c1 <- cond_pairs[p, 1]; c2 <- cond_pairs[p, 2]
  res <- vector("list", length(regions))
  for (i in seq_along(regions)) {
    reg <- regions[i]
    v1 <- long %>% filter(region_acronym == reg, condition == c1) %>% pull(intensity)
    v2 <- long %>% filter(region_acronym == reg, condition == c2) %>% pull(intensity)
    if (length(v1) < 2 || length(v2) < 2) next
    wt <- suppressWarnings(wilcox.test(v1, v2))
    res[[i]] <- data.frame(region_acronym = reg, condition1 = c1, condition2 = c2,
                            n1 = length(v1), n2 = length(v2),
                            median1 = median(v1), median2 = median(v2),
                            p_value = wt$p.value)
  }
  pair_df <- bind_rows(res) %>% mutate(p_adj_BH = p.adjust(p_value, method = "BH"))
  pair_results[[p]] <- pair_df
  cat(sprintf("[%s vs %s] n=%d regions | raw p<0.05: %d | BH-FDR p<0.05 (independent, n=%d tests): %d\n",
              c1, c2, nrow(pair_df), sum(pair_df$p_value < 0.05, na.rm = TRUE),
              nrow(pair_df), sum(pair_df$p_adj_BH < 0.05, na.rm = TRUE)))
}

all_pairs_df <- bind_rows(pair_results)
write.csv(all_pairs_df, "condition_pairwise_wilcoxon_independent_BH.csv", row.names = FALSE)

d <- all_pairs_df %>%
  mutate(pair = paste(condition1, "vs", condition2),
         log2FC = log2((median1 + pseudo) / (median2 + pseudo)),
         neglog10p = -log10(p_value),
         sig_raw = p_value < 0.05,
         sig_BH = p_adj_BH < 0.05)

pairs <- unique(d$pair)

panel_letters <- setNames(letters[seq_along(pairs)], pairs)

## Short labels for the compact grid's panel titles -- the full condition
## names ("forced swim · Acute vs restraint · Acute") are too long to fit
## a two-column panel at large font size and were overflowing into the
## neighbouring panel's title.
cond_abbrev <- c(
  "forced swim · Acute" = "FS Acute", "forced swim · Chronic" = "FS Chronic",
  "restraint · Acute" = "Restraint", "social defeat · Acute" = "Social defeat",
  "tail suspension · Acute" = "Tail suspension"
)
pair_short <- setNames(
  sprintf("%s vs %s", cond_abbrev[sapply(strsplit(pairs, " vs "), `[`, 1)],
                       cond_abbrev[sapply(strsplit(pairs, " vs "), `[`, 2)]),
  pairs
)

## Points in pink (single family, matching the pink-violet scheme used
## elsewhere in the thesis): darker/more saturated pink = higher -log10(p).
## Direction (c1 vs c2 higher) is still readable from the x-axis sign.
## Triangle = also BH-FDR significant, circle = raw only.
FONT_PT <- 14
mm_size <- function(pt) pt / (72.27 / 25.4)  # ggplot geom text 'size' aesthetic is in mm

make_volcano <- function(pair_label, compact = FALSE) {
  dd <- d %>% filter(pair == pair_label) %>%
    mutate(direction = ifelse(log2FC >= 0, "c1 higher", "c2 higher"))
  c1 <- unique(dd$condition1); c2 <- unique(dd$condition2)
  top_label <- dd %>% filter(sig_raw) %>% arrange(p_value) %>% slice_head(n = if (compact) 6 else 12)

  pt_sz <- if (compact) 1.8 else 2.4
  ggplot(dd, aes(x = log2FC, y = neglog10p)) +
    geom_point(data = dd %>% filter(!sig_raw), color = "grey75", alpha = 0.55,
               size = if (compact) 1.2 else 1.4) +
    geom_point(data = dd %>% filter(sig_raw),
               aes(alpha = neglog10p, shape = sig_BH), fill = "#C2185B", color = "white", stroke = 0.4, size = pt_sz) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey60", linewidth = 0.35) +
    geom_text_repel(data = top_label, aes(label = region_acronym), fontface = "bold", color = "#1A1A1A",
                     size = mm_size(FONT_PT), max.overlaps = Inf, min.segment.length = 0, segment.size = 0.4,
                     box.padding = 0.7, point.padding = 0.4, force = 10, force_pull = 0.3, seed = 42) +
    scale_alpha_continuous(range = c(0.55, 1), guide = "none") +
    scale_shape_manual(values = c("TRUE" = 24, "FALSE" = 21), guide = "none") +
    labs(title = if (compact) sprintf("(%s) %s", panel_letters[pair_label], pair_short[pair_label]) else NULL,
         x = if (compact) NULL else sprintf("log2FC (%s / %s)", c1, c2),
         y = if (compact) NULL else expression(-log[10](italic(p)))) +
    theme_minimal(base_size = FONT_PT) +
    theme(plot.title = element_text(size = FONT_PT + 2, face = "bold"),
          axis.title = element_text(size = FONT_PT),
          axis.text = element_text(size = FONT_PT),
          plot.margin = margin(15, 25, 15, 15))
}

## --- Individual, full-size volcano plots (one PNG per pair) ---
for (pl in pairs) {
  fname <- sprintf("volcano_condition_%s.png", gsub("[^A-Za-z0-9]+", "_", pl))
  ggsave(fname, make_volcano(pl, compact = FALSE), width = 10, height = 7, units = "in", dpi = 300)
}

## --- Compact combined grid: all 10 pairs on one page, no overall title ---
plots <- lapply(pairs, make_volcano, compact = TRUE)
n_regions_per_pair <- n_distinct(d$region_acronym)
combined <- wrap_plots(plots, ncol = 2)
ggsave("volcano_condition_pairwise_combined.png", combined, width = 14, height = 22, units = "in", dpi = 300, limitsize = FALSE)

cat("\n=== Summary: BH-FDR significant regions per pair (independent correction) ===\n")
print(d %>% group_by(pair) %>% summarise(n_BH_sig = sum(sig_BH), .groups = "drop") %>% arrange(desc(n_BH_sig)))

cat("\nSaved:\n")
cat("  condition_pairwise_wilcoxon_independent_BH.csv\n")
cat("  volcano_condition_<pair>.png (x10)\n")
cat("  volcano_condition_pairwise_combined.png\n")


# ---------------------------------------------------------------------------
# volcano_sex_by_condition
# ---------------------------------------------------------------------------
# Volcano plots: sex effect per region, WITHIN each condition -- the
# condition-specific counterpart to volcano_sex_region.R (which pooled all
# conditions). Uses the already-computed per-condition sex tests
# (sex_kruskal_by_condition.R -> sex_kruskal_by_condition_region.csv), each
# with its own independent BH-FDR correction (per condition, 636-666
# regions) -- no new testing, just visualization.

library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)

setwd("~/Desktop/Bachelor")
d <- read.csv("sex_kruskal_by_condition_region.csv", stringsAsFactors = FALSE)

pseudo <- 0.001
d <- d %>%
  mutate(log2FC = log2((median_Male + pseudo) / (median_Female + pseudo)),
         neglog10p = -log10(p_value),
         sig_raw = p_value < 0.05,
         sig_BH = p_adj_BH < 0.05)

conditions <- unique(d$condition)
cat("Conditions:", paste(conditions, collapse = ", "), "\n")

panel_letters <- setNames(letters[seq_along(conditions)], conditions)

## Direction-based colour scheme (pink family, matching the pink-violet
## style used elsewhere): crimson-pink = higher in Male, violet-pink =
## higher in Female, both scaled darker with stronger -log10(p). Triangle
## = also BH-FDR significant, circle = raw only.
FONT_PT <- 14
mm_size <- function(pt) pt / (72.27 / 25.4)  # ggplot geom text 'size' aesthetic is in mm

make_volcano <- function(cond, compact = FALSE) {
  dd <- d %>% filter(condition == cond) %>%
    mutate(direction = ifelse(log2FC >= 0, "Male higher", "Female higher"))
  top_label <- dd %>% filter(sig_raw) %>% arrange(p_value) %>% slice_head(n = if (compact) 8 else 12)

  pt_sz <- if (compact) 1.6 else 2.4
  ggplot(dd, aes(x = log2FC, y = neglog10p)) +
    geom_point(data = dd %>% filter(!sig_raw), color = "grey75", alpha = 0.55,
               size = if (compact) 1.0 else 1.4) +
    geom_point(data = dd %>% filter(sig_raw, direction == "Male higher"),
               aes(alpha = neglog10p, shape = sig_BH), fill = "#C2185B", color = "white", stroke = 0.4, size = pt_sz) +
    geom_point(data = dd %>% filter(sig_raw, direction == "Female higher"),
               aes(alpha = neglog10p, shape = sig_BH), fill = "#8E24AA", color = "white", stroke = 0.4, size = pt_sz) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey60", linewidth = 0.35) +
    geom_text_repel(data = top_label, aes(label = region_acronym), fontface = "bold", color = "#1A1A1A",
                     size = mm_size(FONT_PT), max.overlaps = Inf, min.segment.length = 0, segment.size = 0.3,
                     box.padding = 0.7, point.padding = 0.4, force = 8, force_pull = 0.4,
                     seed = 42) +
    scale_alpha_continuous(range = c(0.55, 1), guide = "none") +
    scale_shape_manual(values = c("TRUE" = 24, "FALSE" = 21), guide = "none") +
    labs(title = if (compact) sprintf("(%s) %s", panel_letters[cond], cond) else NULL,
         x = "log2FC (Male / Female)", y = expression(-log[10](italic(p)))) +
    theme_minimal(base_size = FONT_PT) +
    theme(plot.title = element_text(size = FONT_PT + 2, face = "bold"),
          axis.title = element_text(size = FONT_PT),
          axis.text = element_text(size = FONT_PT))
}

for (cond in conditions) {
  fname <- sprintf("volcano_sex_%s.png", gsub("[^A-Za-z0-9]+", "_", cond))
  ggsave(fname, make_volcano(cond, compact = FALSE), width = 10, height = 7, units = "in", dpi = 300)
}

plots <- lapply(conditions, make_volcano, compact = TRUE)
combined <- wrap_plots(plots, ncol = 2)
ggsave("volcano_sex_by_condition_combined.png", combined, width = 14, height = 10, units = "in", dpi = 300, limitsize = FALSE)

cat("\n=== Summary: sex-effect volcano per condition ===\n")
print(d %>% group_by(condition) %>%
        summarise(n_regions = n(), n_raw_sig = sum(sig_raw), n_BH_sig = sum(sig_BH), .groups = "drop"))

cat("\nSaved: volcano_sex_<condition>.png (x4), volcano_sex_by_condition_combined.png\n")


# ---------------------------------------------------------------------------
# volcano_sex_region
# ---------------------------------------------------------------------------
# Volcano plot: sex effect per region, from sex_kruskal_by_region.csv
# (sex_kruskal_stats.R). x-axis = signed log2 fold-change (Male / Female
# median intensity), y-axis = -log10(raw p-value). Complements the KW test
# table with a single visual summary of all 669 regions at once.

library(dplyr)
library(ggplot2)
library(ggrepel)

setwd("~/Desktop/Bachelor")
d <- read.csv("sex_kruskal_by_region.csv", stringsAsFactors = FALSE)

pseudo <- 0.001
d <- d %>%
  mutate(log2FC = log2((median_Male + pseudo) / (median_Female + pseudo)),
         neglog10p = -log10(p_value),
         sig_raw = p_value < 0.05,
         sig_BH = p_adj_BH < 0.05)

cat("n regions:", nrow(d), "| raw p<0.05:", sum(d$sig_raw), "| BH-FDR p<0.05:", sum(d$sig_BH), "\n")

top_label <- d %>% filter(sig_raw) %>% arrange(p_value) %>% slice_head(n = 16)

## Direction-based colour scheme (pink family): crimson-pink = higher in
## Male, violet-pink = higher in Female, both scaled darker with stronger
## -log10(p).
d <- d %>% mutate(direction = ifelse(log2FC >= 0, "Male higher", "Female higher"))

p <- ggplot(d, aes(x = log2FC, y = neglog10p)) +
  geom_point(data = d %>% filter(!sig_raw), color = "grey75", alpha = 0.55, size = 2.6) +
  geom_point(data = d %>% filter(sig_raw, direction == "Male higher"),
             aes(alpha = neglog10p), shape = 21, fill = "#C2185B", color = "white", stroke = 0.4, size = 4.4) +
  geom_point(data = d %>% filter(sig_raw, direction == "Female higher"),
             aes(alpha = neglog10p), shape = 21, fill = "#8E24AA", color = "white", stroke = 0.4, size = 4.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey60", linewidth = 0.4) +
  geom_text_repel(data = top_label, aes(label = region_acronym), size = 18, fontface = "bold",
                   color = "#1A1A1A", max.overlaps = Inf, min.segment.length = 0, segment.size = 0.4,
                   box.padding = 0.7, point.padding = 0.4, force = 8, force_pull = 0.4, seed = 42) +
  scale_alpha_continuous(range = c(0.55, 1), guide = "none") +
  labs(x = "log2 fold-change (Male / Female median intensity)",
       y = expression(-log[10](italic(p)))) +
  theme_minimal(base_size = 66) +
  theme(
    axis.title = element_text(size = 60),
    axis.text = element_text(size = 52)
  )

ggsave("volcano_sex_by_region.png", p, width = 42, height = 24, units = "in", dpi = 300)
cat("Saved: volcano_sex_by_region.png\n")


# ---------------------------------------------------------------------------
# fs_acute_vs_chronic_pooled_volcano
# ---------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(ggrepel)

setwd("~/Desktop/Bachelor")
d <- read.csv("fs_acute_vs_chronic_POOLED_region.csv", stringsAsFactors = FALSE)

pseudo <- 0.001
d <- d %>%
  mutate(log2FC = log2((median_forced.swim...Chronic + pseudo) / (median_forced.swim...Acute + pseudo)),
         neglog10p = -log10(p_value),
         sig_raw = p_value < 0.05)

top_label <- d %>% filter(sig_raw) %>% arrange(p_value) %>% slice_head(n = 12)

p <- ggplot(d, aes(x = log2FC, y = neglog10p)) +
  geom_point(data = d %>% filter(!sig_raw), color = "grey75", alpha = 0.55, size = 5) +
  geom_point(data = d %>% filter(sig_raw), aes(alpha = neglog10p), shape = 21,
             fill = "#8E24AA", color = "white", stroke = 0.6, size = 9) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey60", linewidth = 0.6) +
  geom_text_repel(data = top_label, aes(label = region_acronym), size = 18, fontface = "bold",
                   color = "#1A1A1A", max.overlaps = Inf, min.segment.length = 0, segment.size = 0.5) +
  scale_alpha_continuous(range = c(0.55, 1), guide = "none") +
  labs(x = "log2 fold-change (Chronic / Acute median intensity)",
       y = expression(-log[10](italic(p))),
       title = "Forced swim: Chronic vs Acute (both sexes pooled)",
       subtitle = "Region-level Wilcoxon rank-sum, raw p; none survive BH-FDR (588 tests)") +
  theme_minimal(base_size = 62) +
  theme(
    plot.title = element_text(size = 65, face = "bold"),
    plot.subtitle = element_text(size = 46, color = "grey30"),
    axis.title = element_text(size = 54),
    axis.text = element_text(size = 47)
  )

ggsave("fs_acute_vs_chronic_pooled_volcano.png", p, width = 34, height = 22, units = "in", dpi = 200)
cat("Saved: fs_acute_vs_chronic_pooled_volcano.png\n")

