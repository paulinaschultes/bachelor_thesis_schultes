# =============================================================================
# Retrograde tracing: cell counts, candidate input regions and circuit diagrams
#
# Consolidated from the following original scripts:
#   - retrograde_common.R
#   - combine_ret1_refatlasregions.R
#   - combine_ret1_refatlasregions_xlsx.R
#   - combine_ret3_refatlasregions.R
#   - ret1_input_regions_low_deviation.R
#   - input_regions_low_deviation_both_cohorts.R
#   - ret1_vs_ret3_wilcoxon_and_overlap.R
#   - ret1_vs_ret3_estimation_plots.R
#   - fs_circuit_diagram.R
#   - fs_circuit_diagram_by_condition.R
#   - fs_circuit_diagram_common.R
# =============================================================================


# ---------------------------------------------------------------------------
# retrograde_common
# ---------------------------------------------------------------------------
# Shared setup for all retrograde (RET1/RET3) analysis scripts:
#   - loads the combined per-animal RefAtlasRegions tables
#   - maps Region ID -> acronym + major division (same Allen CCF top-level
#     lookup as prepare_heatmap_data.R) and EXCLUDES fiber tracts and
#     ventricular systems, consistent with every anterograde analysis in
#     this thesis. Filtering by major-division membership (not a name
#     regex) matters here: many gray-matter nuclei have "tract" in their
#     name (e.g. "Nucleus of the lateral olfactory tract") and would be
#     wrongly dropped by a naive text match.

library(dplyr)

structures <- read.csv("~/Desktop/tractquant/allen_mouse_25um_v1.2/structures.csv", stringsAsFactors = FALSE)

major_ids <- c(
  "315"  = "Isocortex", "698"  = "Olfactory areas", "1089" = "Hippocampal formation",
  "703"  = "Cortical subplate", "477"  = "Striatum", "803"  = "Pallidum",
  "549"  = "Thalamus", "1097" = "Hypothalamus", "313"  = "Midbrain",
  "771"  = "Pons", "354"  = "Medulla", "512"  = "Cerebellum",
  "1009" = "fiber tracts", "73"   = "ventricular systems"
)

get_major_division <- function(path) {
  ids <- strsplit(gsub("^/|/$", "", path), "/")[[1]]
  hit <- ids[ids %in% names(major_ids)]
  if (length(hit) == 0) return(NA_character_)
  major_ids[[hit[length(hit)]]]
}

structures <- structures %>%
  rowwise() %>%
  mutate(major_division = get_major_division(structure_id_path)) %>%
  ungroup()

EXCLUDED_DIVISIONS <- c("fiber tracts", "ventricular systems")

load_retrograde <- function(path, exclude_fiber = TRUE) {
  df <- read.csv(path, stringsAsFactors = FALSE) %>%
    left_join(structures %>% select(id, acronym, major_division), by = c("Region.ID" = "id"))
  if (exclude_fiber) {
    n_before <- n_distinct(df$Region.ID)
    df <- df %>% filter(is.na(major_division) | !(major_division %in% EXCLUDED_DIVISIONS))
    n_after <- n_distinct(df$Region.ID)
    cat(sprintf("  %s: excluded %d fiber-tract/ventricular regions (%d -> %d)\n",
                basename(path), n_before - n_after, n_before, n_after))
  }
  df
}


# ---------------------------------------------------------------------------
# combine_ret1_refatlasregions
# ---------------------------------------------------------------------------
# Combines the per-animal whole-series RefAtlasRegions.csv Nutil output
# (one row per Allen atlas region, whole-animal totals across all sections)
# for RET1_1 - RET1_4 into a single long-format table with an animal_id
# column, analogous in spirit to how ALL_DT.csv combines the anterograde
# TractQuant output across animals.

library(dplyr)
library(readr)

files <- c(
  RET1_1 = "/Volumes/Extreme_SSD/RET/RET1/RET1_1/nutil/Reports/Ret1_RefAtlasRegions/Ret1_RefAtlasRegions.csv",
  RET1_2 = "/Volumes/Extreme_SSD/RET/RET1/RET1_2/nutil_output/Reports/Ret2_RefAtlasRegions/Ret2_RefAtlasRegions.csv",
  RET1_3 = "/Volumes/Extreme_SSD/RET/RET1/RET1_3/nutil/Reports/Ret3_RefAtlasRegions/Ret3_RefAtlasRegions.csv",
  RET1_4 = "/Volumes/Extreme_SSD/RET/RET1/RET1_4/nutil/Reports/Ret4_RefAtlasRegions/Ret4_RefAtlasRegions.csv"
)

read_one <- function(animal_id, path) {
  df <- read_delim(path, delim = ";", show_col_types = FALSE)
  df <- df[, !grepl("^\\.\\.\\.", names(df))]   # drop trailing empty columns
  df$animal_id <- animal_id
  df
}

combined <- bind_rows(lapply(names(files), function(a) read_one(a, files[[a]])))
combined <- combined %>% relocate(animal_id, .before = 1)

out_path <- "~/Desktop/Bachelor/retrograde/RET1_RefAtlasRegions_combined.csv"
write.csv(combined, out_path, row.names = FALSE)

cat("Animals combined:", paste(names(files), collapse = ", "), "\n")
cat("n rows total:", nrow(combined), "(", nrow(combined) / length(files), "regions x", length(files), "animals)\n")
cat("Saved:", out_path, "\n\n")

cat("=== Sanity check: total object count per animal ===\n")
print(combined %>% group_by(animal_id) %>% summarise(total_objects = sum(`Object count`, na.rm = TRUE)))


# ---------------------------------------------------------------------------
# combine_ret1_refatlasregions_xlsx
# ---------------------------------------------------------------------------
# Same combined RET1_1-4 RefAtlasRegions table as
# combine_ret1_refatlasregions.R, saved as .xlsx instead of/in addition to
# .csv (one sheet with all animals stacked long-format, plus one sheet per
# animal for easier manual inspection).

library(dplyr)
library(readr)
library(writexl)

files <- c(
  RET1_1 = "/Volumes/Extreme_SSD/RET/RET1/RET1_1/nutil/Reports/Ret1_RefAtlasRegions/Ret1_RefAtlasRegions.csv",
  RET1_2 = "/Volumes/Extreme_SSD/RET/RET1/RET1_2/nutil_output/Reports/Ret2_RefAtlasRegions/Ret2_RefAtlasRegions.csv",
  RET1_3 = "/Volumes/Extreme_SSD/RET/RET1/RET1_3/nutil/Reports/Ret3_RefAtlasRegions/Ret3_RefAtlasRegions.csv",
  RET1_4 = "/Volumes/Extreme_SSD/RET/RET1/RET1_4/nutil/Reports/Ret4_RefAtlasRegions/Ret4_RefAtlasRegions.csv"
)

read_one <- function(animal_id, path) {
  df <- read_delim(path, delim = ";", show_col_types = FALSE)
  df <- df[, !grepl("^\\.\\.\\.", names(df))]
  df$animal_id <- animal_id
  df %>% relocate(animal_id, .before = 1)
}

per_animal <- lapply(names(files), function(a) read_one(a, files[[a]]))
names(per_animal) <- names(files)
combined <- bind_rows(per_animal)

sheets <- c(list(All_animals = combined), per_animal)

out_path <- path.expand("~/Desktop/Bachelor/retrograde/RET1_RefAtlasRegions_combined.xlsx")
write_xlsx(sheets, out_path)

cat("Saved:", out_path, "\n")
cat("Sheets:", paste(names(sheets), collapse = ", "), "\n\n")
cat("=== Sanity check: total object count per animal ===\n")
print(combined %>% group_by(animal_id) %>% summarise(total_objects = sum(`Object count`, na.rm = TRUE)))


# ---------------------------------------------------------------------------
# combine_ret3_refatlasregions
# ---------------------------------------------------------------------------
# Combines the per-animal whole-series RefAtlasRegions.csv Nutil output
# for RET3_1 and RET3_3 into a single long-format table with an
# animal_id column -- same design as combine_ret1_refatlasregions.R.
# Paths updated after the Nutil object-colour fix (black/0,0,0 selected
# explicitly) was rerun for both animals.

library(dplyr)
library(readr)

files <- c(
  RET3_1 = "/Volumes/PortableSSD/RET3/RET_3_1/nutil/Reports/xx_RefAtlasRegions/xx_RefAtlasRegions.csv",
  RET3_3 = "/Volumes/PortableSSD/RET3/RET_3_3/nutil/Reports/xxx_RefAtlasRegions/xxx_RefAtlasRegions.csv"
)

read_one <- function(animal_id, path) {
  df <- read_delim(path, delim = ";", show_col_types = FALSE)
  df <- df[, !grepl("^\\.\\.\\.", names(df))]   # drop trailing empty columns
  df$animal_id <- animal_id
  df
}

combined <- bind_rows(lapply(names(files), function(a) read_one(a, files[[a]])))
combined <- combined %>% relocate(animal_id, .before = 1)

out_path <- "~/Desktop/Bachelor/retrograde/RET3_RefAtlasRegions_combined.csv"
write.csv(combined, out_path, row.names = FALSE)

cat("Animals combined:", paste(names(files), collapse = ", "), "\n")
cat("n rows total:", nrow(combined), "(", nrow(combined) / length(files), "regions x", length(files), "animals)\n")
cat("Saved:", out_path, "\n\n")

cat("=== Sanity check: total object count per animal ===\n")
print(combined %>% group_by(animal_id) %>% summarise(total_objects = sum(`Object count`, na.rm = TRUE)))


# ---------------------------------------------------------------------------
# ret1_input_regions_low_deviation
# ---------------------------------------------------------------------------
# Identifies candidate "input regions" from the RET1 retrograde tracing data:
# regions with consistently HIGH signal (Load) and LOW variability across
# the 4 animals -- same "combined rank" logic as sex_top_regions_low_deviation.R
# / Figure 12 for the anterograde data (high intensity rank + low CV rank).
#
# RET1_4 shows 0 objects in every single region (a known Nutil pipeline
# issue, diagnosed separately -- ilastik export format mismatch). Including
# it would artificially deflate every mean and inflate every CV, so the
# primary analysis uses n=3 (RET1_1-3); RET1_4 is kept in a parallel n=4
# run only for comparison/transparency.

library(dplyr)
library(tidyr)
library(ggplot2)

setwd("~/Desktop/Bachelor/retrograde")
df <- read.csv("RET1_RefAtlasRegions_combined.csv", stringsAsFactors = FALSE)

## All 4 animals share one condition -- no condition split needed, just
## summarise across animals per region.

summarise_regions <- function(data, label) {
  data %>%
    group_by(Region.ID, Region.Name) %>%
    summarise(
      n_animals   = n(),
      mean_load   = mean(Load, na.rm = TRUE),
      sd_load     = sd(Load, na.rm = TRUE),
      cv_load     = sd_load / mean_load,
      mean_count  = mean(Object.count, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(mean_load > 0, is.finite(cv_load)) %>%          # drop regions with no signal in any animal
    mutate(
      rank_intensity = rank(-mean_load),                    # higher load = better rank
      rank_cv         = rank(cv_load),                       # lower CV = better rank
      combined_rank   = rank_intensity + rank_cv,
      set = label
    ) %>%
    arrange(combined_rank)
}

## --- Primary: n=3 (RET1_1-3), RET1_4 excluded ---
df_n3 <- df %>% filter(animal_id != "RET1_4")
res_n3 <- summarise_regions(df_n3, "n=3 (RET1_4 excluded)")

## --- Comparison only: n=4, all animals ---
res_n4 <- summarise_regions(df, "n=4 (all animals)")

write.csv(res_n3, "RET1_input_regions_low_deviation_n3.csv", row.names = FALSE)
write.csv(res_n4, "RET1_input_regions_low_deviation_n4_comparison.csv", row.names = FALSE)

cat("=== Top 20 candidate input regions (n=3, RET1_1-3, ranked by combined rank) ===\n")
print(as.data.frame(head(res_n3 %>% select(Region.Name, mean_load, cv_load, mean_count, combined_rank), 20)))

cat("\n=== For comparison: same regions in the n=4 (RET1_4-included) version ===\n")
top20_names <- head(res_n3$Region.Name, 20)
print(as.data.frame(res_n4 %>% filter(Region.Name %in% top20_names) %>%
                       select(Region.Name, mean_load, cv_load, combined_rank) %>%
                       arrange(combined_rank)))

## --- Plot, same style as Figure 12 (Combined Rank) ---
p <- ggplot(res_n3, aes(x = cv_load, y = mean_load, color = combined_rank)) +
  geom_point(alpha = 0.7) +
  ggrepel::geom_text_repel(data = head(res_n3, 20), aes(label = Region.Name),
                            size = 2.8, max.overlaps = 30, color = "black") +
  scale_color_gradient(low = "#8E24AA", high = "#FDF0F5", name = "Combined\nrank\n(darker=better)") +
  scale_y_log10() +
  labs(title = "RET1 retrograde tracing: candidate input regions (n=3, RET1_4 excluded)",
       subtitle = "High mean Load + low coefficient of variation across animals",
       x = "Coefficient of variation across animals", y = "Mean Load (log scale)") +
  theme_minimal(base_size = 11)
ggsave("RET1_input_regions_low_deviation.png", p, width = 9, height = 7, dpi = 300)

cat("\nSaved:\n")
cat("  RET1_input_regions_low_deviation_n3.csv (primary, n=3)\n")
cat("  RET1_input_regions_low_deviation_n4_comparison.csv (comparison, n=4)\n")
cat("  RET1_input_regions_low_deviation.png\n")


# ---------------------------------------------------------------------------
# input_regions_low_deviation_both_cohorts
# ---------------------------------------------------------------------------
# Top candidate "input regions" (high mean Load, low coefficient of
# variation across animals) for both retrograde cohorts:
#   RET1 (n=3: RET1_1-3, RET1_4 excluded -- 0 objects, known Nutil issue)
#     = acute restraint, Female
#   RET3 (n=2: RET3_1, RET3_3, after the Nutil object-colour fix)
#     = acute forced swim, Male
# Same "combined rank" logic as sex_top_regions_low_deviation.R /
# Figure 12 for the anterograde data. Fiber tracts and ventricular systems
# excluded via retrograde_common.R (Allen CCF major-division lookup),
# consistent with every anterograde analysis in this thesis.

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

setwd("~/Desktop/Bachelor/retrograde")
source("retrograde_common.R")

summarise_regions <- function(data, label) {
  data %>%
    group_by(Region.ID, Region.Name, acronym) %>%
    summarise(
      n_animals  = n(),
      mean_load  = mean(Load, na.rm = TRUE),
      sd_load    = sd(Load, na.rm = TRUE),
      cv_load    = sd_load / mean_load,
      mean_count = mean(Object.count, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(mean_load > 0, is.finite(cv_load)) %>%
    mutate(
      rank_intensity = rank(-mean_load),
      rank_cv        = rank(cv_load),
      combined_rank  = rank_intensity + rank_cv,
      cohort = label
    ) %>%
    arrange(combined_rank)
}

cat("=== Loading (fiber tracts + ventricular systems excluded) ===\n")
ret1 <- load_retrograde("RET1_RefAtlasRegions_combined.csv") %>% filter(animal_id != "RET1_4")
ret3 <- load_retrograde("RET3_RefAtlasRegions_combined.csv")

res_ret1 <- summarise_regions(ret1, "RET1: acute restraint . Female (n=3)")
res_ret3 <- summarise_regions(ret3, "RET3: acute forced swim . Male (n=2)")

write.csv(res_ret1, "input_regions_RET1_restraint_female_n3.csv", row.names = FALSE)
write.csv(res_ret3, "input_regions_RET3_forcedswim_male_n2.csv", row.names = FALSE)

cat("\n=== Top 15: RET1, acute restraint . Female (n=3) ===\n")
print(as.data.frame(head(res_ret1 %>% select(acronym, Region.Name, mean_load, cv_load, mean_count, combined_rank), 15)))

cat("\n=== Top 15: RET3, acute forced swim . Male (n=2) ===\n")
print(as.data.frame(head(res_ret3 %>% select(acronym, Region.Name, mean_load, cv_load, mean_count, combined_rank), 15)))

top20_ret1 <- head(res_ret1$Region.Name, 20)
top20_ret3 <- head(res_ret3$Region.Name, 20)
overlap_names <- intersect(top20_ret1, top20_ret3)
cat("\n=== Regions in BOTH cohorts' top 20 ===\n")
print(if (length(overlap_names) == 0) "(none)" else overlap_names)

## --- Plots: ONE clean, wide, high-legibility panel per cohort ---
## (base_size, label size, canvas width all match the established
## sex_top_regions_low_deviation.R convention -- base_size=12-14, label
## size=3.5-4 bold, wide canvas rather than tall stacked facets)
make_plot <- function(res, top_n, title, fname) {
  top <- head(res, top_n)
  p <- ggplot(res, aes(x = cv_load, y = mean_load, color = combined_rank)) +
    geom_point(alpha = 0.65, size = 10.2) +
    geom_text_repel(data = top, aes(label = acronym), size = 18, fontface = "bold",
                     color = "#1A1A1A", max.overlaps = Inf, min.segment.length = 0,
                     segment.size = 0.7, box.padding = 0.4) +
    scale_color_gradient(low = "#8E24AA", high = "#FDF0F5", name = "Combined\nrank") +
    scale_y_log10() +
    labs(title = title,
         subtitle = "High mean Load + low coefficient of variation across animals (fiber tracts excluded)",
         x = "Coefficient of variation across animals", y = "Mean Load (log scale)") +
    theme_minimal(base_size = 90) +
    theme(
      plot.title = element_text(face = "bold", size = 90),
      plot.subtitle = element_text(size = 66, color = "grey30"),
      axis.title = element_text(size = 84),
      axis.text = element_text(size = 72),
      legend.title = element_text(size = 78),
      legend.text = element_text(size = 66)
    )
  ggsave(fname, p, width = 60, height = 36, dpi = 150, limitsize = FALSE)
  cat("Saved:", fname, "\n")
}

make_plot(res_ret1, 20, "RET1: acute restraint, Female (n=3) -- candidate input regions",
          "input_regions_low_deviation_RET1.png")
make_plot(res_ret3, 20, "RET3: acute forced swim, Male (n=2) -- candidate input regions",
          "input_regions_low_deviation_RET3.png")

cat("\nSaved CSVs:\n")
cat("  input_regions_RET1_restraint_female_n3.csv\n")
cat("  input_regions_RET3_forcedswim_male_n2.csv\n")
cat("Saved plots:\n")
cat("  input_regions_low_deviation_RET1.png\n")
cat("  input_regions_low_deviation_RET3.png\n")


# ---------------------------------------------------------------------------
# ret1_vs_ret3_wilcoxon_and_overlap
# ---------------------------------------------------------------------------
# Retrograde tracing (input regions): two analyses building on
# input_regions_low_deviation_both_cohorts.R
#
#   1. Region-level Wilcoxon rank-sum test, RET1 (acute restraint, Female,
#      n=3) vs RET3 (acute forced swim, Male, n=2), on Load -- same design
#      as the anterograde condition-pairwise volcano plots (independent
#      BH-FDR per test family). IMPORTANT: this comparison is doubly
#      confounded (sex AND stressor both differ between cohorts, cannot be
#      disentangled) and, with n=3 vs n=2, the exact Wilcoxon rank-sum test
#      has a hard floor of p=0.2 -- no region can reach even nominal
#      p<0.05 by construction. Reported as a ranked/exploratory comparison,
#      not a significance-hunting exercise.
#
#   2. Input/output overlap: regions with detected signal in BOTH the
#      anterograde (output, DT/TRAP tracing FROM the injection site) and
#      retrograde (input, RET, projections TO the injection site) datasets
#      -- candidate reciprocally-connected regions.
#
# Fiber tracts and ventricular systems excluded via retrograde_common.R,
# consistent with every anterograde analysis in this thesis.

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

setwd("~/Desktop/Bachelor/retrograde")
source("retrograde_common.R")

cat("=== Loading (fiber tracts + ventricular systems excluded) ===\n")
ret1 <- load_retrograde("RET1_RefAtlasRegions_combined.csv") %>% filter(animal_id != "RET1_4")
ret3 <- load_retrograde("RET3_RefAtlasRegions_combined.csv")

## ============================================================
## 1. Wilcoxon rank-sum test per region, RET1 vs RET3
## ============================================================
regions_common <- intersect(unique(ret1$Region.ID), unique(ret3$Region.ID))

wilcox_results <- vector("list", length(regions_common))
for (i in seq_along(regions_common)) {
  rid <- regions_common[i]
  v1 <- ret1$Load[ret1$Region.ID == rid]
  v3 <- ret3$Load[ret3$Region.ID == rid]
  if (length(v1) < 2 || length(v3) < 2) next
  if (sum(v1) == 0 && sum(v3) == 0) next   # no signal in either cohort -- not testable/meaningful

  wt <- suppressWarnings(wilcox.test(v1, v3))
  wilcox_results[[i]] <- data.frame(
    Region.ID = rid,
    Region.Name = ret1$Region.Name[ret1$Region.ID == rid][1],
    acronym = ret1$acronym[ret1$Region.ID == rid][1],
    median_RET1_restraint_F = median(v1), median_RET3_forcedswim_M = median(v3),
    p_value = wt$p.value
  )
}
wilcox_df <- bind_rows(wilcox_results) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    pseudo = min(median_RET1_restraint_F[median_RET1_restraint_F > 0], median_RET3_forcedswim_M[median_RET3_forcedswim_M > 0], na.rm = TRUE) / 2,
    log2FC = log2((median_RET3_forcedswim_M + pseudo) / (median_RET1_restraint_F + pseudo)),
    neglog10p = -log10(p_value),
    sig_raw = p_value < 0.05
  ) %>%
  arrange(p_value)

write.csv(wilcox_df, "RET1_vs_RET3_wilcoxon_by_region.csv", row.names = FALSE)

cat("=== RET1 (restraint, F, n=3) vs RET3 (forced swim, M, n=2): Wilcoxon per region ===\n")
cat("n regions tested:", nrow(wilcox_df), "\n")
cat("Minimum achievable exact p-value at n=3 vs n=2:", min(wilcox_df$p_value), "\n")
cat("Regions with raw p<0.05:", sum(wilcox_df$sig_raw), "(structurally impossible at this n -- see note)\n")
cat("Regions at the p-floor (p=0.2, maximally separated):", sum(wilcox_df$p_value <= 0.2 + 1e-9), "\n\n")

cat("=== Top 20 by p-value (i.e. most consistently separated between cohorts) ===\n")
print(as.data.frame(head(wilcox_df %>% select(acronym, Region.Name, median_RET1_restraint_F, median_RET3_forcedswim_M, log2FC, p_value), 20)))

## --- Volcano plot: direction-coloured (distinct from the pink-violet
## intensity gradient used elsewhere), no title/subtitle (see figure caption) ---
top_labels <- wilcox_df %>% filter(p_value <= 0.2 + 1e-9) %>% slice_max(order_by = abs(log2FC), n = 24)
wilcox_df <- wilcox_df %>% mutate(direction = ifelse(log2FC >= 0, "RET3 higher", "RET1 higher"))

p_volcano <- ggplot(wilcox_df, aes(x = log2FC, y = neglog10p)) +
  geom_point(data = wilcox_df %>% filter(p_value > 0.2 + 1e-9), color = "grey75", alpha = 0.55, size = 2.6) +
  geom_point(data = wilcox_df %>% filter(p_value <= 0.2 + 1e-9, direction == "RET3 higher"),
             aes(alpha = neglog10p), shape = 21, fill = "#F57C00", color = "white", stroke = 0.4, size = 4.4) +
  geom_point(data = wilcox_df %>% filter(p_value <= 0.2 + 1e-9, direction == "RET1 higher"),
             aes(alpha = neglog10p), shape = 21, fill = "#00897B", color = "white", stroke = 0.4, size = 4.4) +
  geom_text_repel(data = top_labels, aes(label = acronym), size = 5.2, fontface = "bold",
                   color = "#1A1A1A", max.overlaps = Inf, min.segment.length = 0, segment.size = 0.4) +
  scale_alpha_continuous(range = c(0.55, 1), guide = "none") +
  geom_hline(yintercept = -log10(0.2), linetype = "dashed", color = "grey40", linewidth = 0.6) +
  labs(x = "log2FC (RET3 forced swim, Male / RET1 restraint, Female)", y = expression(-log[10](p))) +
  theme_minimal(base_size = 16) +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13)
  )
ggsave("RET1_vs_RET3_wilcoxon_volcano.png", p_volcano, width = 13, height = 8.5, dpi = 300)


## ============================================================
## 2. Input/output overlap: regions detected in BOTH datasets
## ============================================================
anterograde <- read.csv("~/Desktop/Bachelor/condition_region_intensity_for_brainrender.csv", stringsAsFactors = FALSE) %>%
  mutate(mean_output_intensity = rowMeans(across(`X1`:`X5`), na.rm = TRUE)) %>%
  select(acronym = region_acronym, mean_output_intensity)

retrograde_combined <- bind_rows(
  ret1 %>% mutate(cohort = "RET1 (restraint, F)"),
  ret3 %>% mutate(cohort = "RET3 (forced swim, M)")
) %>%
  filter(!is.na(acronym)) %>%
  group_by(acronym) %>%
  summarise(mean_input_load = mean(Load, na.rm = TRUE), mean_input_count = mean(Object.count, na.rm = TRUE), .groups = "drop")

overlap <- inner_join(anterograde, retrograde_combined, by = "acronym") %>%
  filter(mean_output_intensity > 0, mean_input_load > 0) %>%
  mutate(
    output_rank = rank(-mean_output_intensity),
    input_rank  = rank(-mean_input_load),
    reciprocal_rank = output_rank + input_rank
  ) %>%
  arrange(reciprocal_rank)

write.csv(overlap, "input_output_overlap_regions.csv", row.names = FALSE)

cat("\n\n=== Input/output overlap (fiber tracts excluded) ===\n")
cat("n regions with detected signal in BOTH anterograde (output) and retrograde (input) data:", nrow(overlap), "\n")
cat("of", nrow(anterograde), "output regions and", nrow(retrograde_combined), "input regions with any signal\n\n")

cat("=== Top 20 candidate reciprocally-connected regions (strong in both) ===\n")
print(as.data.frame(head(overlap %>% select(acronym, mean_output_intensity, mean_input_load, reciprocal_rank), 20)))

## --- Scatter plot: larger fonts/canvas, matching thesis convention ---
top_recip <- head(overlap, 20)
p_overlap <- ggplot(overlap, aes(x = mean_output_intensity, y = mean_input_load)) +
  geom_point(aes(color = reciprocal_rank), alpha = 0.7, size = 9.9) +
  geom_text_repel(data = top_recip, aes(label = acronym), size = 18, fontface = "bold",
                   max.overlaps = Inf, color = "#1A1A1A", min.segment.length = 0, segment.size = 0.7) +
  scale_color_gradient(low = "#8E24AA", high = "#FDF0F5", name = "Reciprocal\nrank") +
  scale_x_log10() + scale_y_log10() +
  labs(
    title = "Candidate reciprocally-connected regions: anterograde output vs retrograde input",
    subtitle = sprintf("%d regions with detected signal in both datasets (fiber tracts excluded)", nrow(overlap)),
    x = "Mean anterograde output intensity (log scale, all conditions pooled)",
    y = "Mean retrograde input Load (log scale, RET1+RET3 pooled)"
  ) +
  theme_minimal(base_size = 90) +
  theme(
    plot.title = element_text(face = "bold", size = 90),
    plot.subtitle = element_text(size = 66, color = "grey30"),
    axis.title = element_text(size = 84),
    axis.text = element_text(size = 72)
  )
ggsave("input_output_overlap_scatter.png", p_overlap, width = 60, height = 38, dpi = 150, limitsize = FALSE)

cat("\nSaved:\n")
cat("  RET1_vs_RET3_wilcoxon_by_region.csv, RET1_vs_RET3_wilcoxon_volcano.png\n")
cat("  input_output_overlap_regions.csv, input_output_overlap_scatter.png\n")


# ---------------------------------------------------------------------------
# ret1_vs_ret3_estimation_plots
# ---------------------------------------------------------------------------
# Gardner-Altman estimation plots (dabestr): bootstrap mean difference +
# 95% CI, RET1 (acute restraint, Female, n=3) vs RET3 (acute forced swim,
# Male, n=2), Load. Same design and package as
# condition_estimation_plots_and_cluster_stability.R section 1 -- the
# effect-size complement to the Wilcoxon test, which is uninformative here
# since the exact test cannot go below p=0.2 at this sample size. An
# estimation plot still shows HOW BIG the difference is and how much the
# animal-level values overlap, even when a p-value can't be trusted.
#
# Regions shown: the 4 regions common to both cohorts' top-20 low-deviation
# list, plus the next 4 highest-combined-rank regions unique to each
# cohort (8 total) -- same regions named in the results text.

library(dplyr)
library(dabestr)
library(ggplot2)
library(patchwork)

setwd("~/Desktop/Bachelor/retrograde")
source("retrograde_common.R")

ret1 <- load_retrograde("RET1_RefAtlasRegions_combined.csv", exclude_fiber = TRUE) %>%
  filter(animal_id != "RET1_4")
ret3 <- load_retrograde("RET3_RefAtlasRegions_combined.csv", exclude_fiber = TRUE)

shared_regions <- c(
  "Basolateral amygdalar nucleus, anterior part",
  "Paraventricular nucleus of the thalamus",
  "Cortical amygdalar area, anterior part",
  "Parataenial nucleus"
)
ret1_unique <- c("Infralimbic area, layer 5", "Prelimbic area, layer 5")
ret3_unique <- c("Basolateral amygdalar nucleus, posterior part", "Central amygdalar nucleus, lateral part")

regions_to_plot <- c(shared_regions, ret1_unique, ret3_unique)

long_df <- bind_rows(
  ret1 %>% filter(Region.Name %in% regions_to_plot) %>%
    transmute(Region.Name, cohort = "RET1 restraint,F", Load),
  ret3 %>% filter(Region.Name %in% regions_to_plot) %>%
    transmute(Region.Name, cohort = "RET3 f.swim,M", Load)
)

group_levels <- c("RET1 restraint,F", "RET3 f.swim,M")
plots <- list()
effect_results <- vector("list", length(regions_to_plot))

for (i in seq_along(regions_to_plot)) {
  reg <- regions_to_plot[i]
  d <- long_df %>% filter(Region.Name == reg)
  if (n_distinct(d$cohort) < 2) next
  db <- d %>% dabestr::load(x = cohort, y = Load, idx = group_levels)
  db_diff <- mean_diff(db)

  boot <- db_diff$boot_result
  effect_results[[i]] <- data.frame(
    Region.Name = reg,
    mean_diff_RET3_minus_RET1 = boot$difference,
    bca_ci_low = boot$bca_ci_low,
    bca_ci_high = boot$bca_ci_high
  )

  p <- dabest_plot(db_diff, TRUE, raw_marker_size = 1.3, es_marker_size = 1.8) +
    ggtitle(reg) +
    theme(plot.title = element_text(size = 10, face = "bold"))
  plots[[reg]] <- p
}

effect_df <- bind_rows(effect_results)
write.csv(effect_df, "RET1_vs_RET3_estimation_effect_sizes.csv", row.names = FALSE)

cat("=== RET1 vs RET3: bootstrap mean difference in Load (RET3 - RET1) ===\n")
print(effect_df)

combined <- wrap_plots(plots, ncol = 2) +
  plot_annotation(
    title = "RET1 (restraint, Female, n=3) vs RET3 (forced swim, Male, n=2): estimation plots",
    subtitle = "Bootstrap mean difference in retrograde Load + 95% CI -- effect-size complement to the underpowered Wilcoxon test above",
    caption = "Left 4 panels: regions shared between both cohorts' top-20 input regions. Right 4 panels: top region unique to each cohort.",
    theme = theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      plot.caption = element_text(size = 9, hjust = 0.5)
    )
  )
ggsave("RET1_vs_RET3_estimation_plots.png", combined, width = 20, height = 11, dpi = 300, limitsize = FALSE)

cat("\nSaved: RET1_vs_RET3_estimation_effect_sizes.csv, RET1_vs_RET3_estimation_plots.png\n")


# ---------------------------------------------------------------------------
# fs_circuit_diagram
# ---------------------------------------------------------------------------
# Circuit-style diagram of the FS injection site's connectivity: FS at the
# centre, arrows pointing IN from the top retrograde-traced input regions
# (afferents) and arrows pointing OUT to the top anterograde output
# regions (efferents, from Figure 12's combined-rank list), with regions
# that are BOTH strong input AND strong output (from
# input_output_overlap_regions.csv) drawn as bidirectional/reciprocal
# edges in a distinct colour. This is the schematic-circuit complement to
# the brainrender 3D projection maps (Figure 24/25) and the scatter-style
# overlap plot (Figure 28) -- same underlying data, laid out the way a
# connectivity diagram is conventionally read in circuit neuroscience.

library(dplyr)
library(igraph)
library(ggraph)
library(tidygraph)
library(ggplot2)
library(ggrepel)

setwd("~/Desktop/Bachelor/retrograde")
source("retrograde_common.R")

N_INPUT_ONLY  <- 8
N_OUTPUT_ONLY <- 8

overlap <- read.csv("input_output_overlap_regions.csv", stringsAsFactors = FALSE)
reciprocal <- overlap %>% slice_min(reciprocal_rank, n = 6) %>% pull(acronym)

## --- input-only candidates: top retrograde regions not already reciprocal ---
ret1 <- load_retrograde("RET1_RefAtlasRegions_combined.csv", exclude_fiber = TRUE) %>% filter(animal_id != "RET1_4")
ret3 <- load_retrograde("RET3_RefAtlasRegions_combined.csv", exclude_fiber = TRUE)
input_ranked <- bind_rows(ret1, ret3) %>%
  filter(!is.na(acronym)) %>%
  group_by(acronym) %>%
  summarise(mean_input_load = mean(Load, na.rm = TRUE), .groups = "drop") %>%
  filter(mean_input_load > 0, !(acronym %in% reciprocal)) %>%
  arrange(desc(mean_input_load)) %>%
  slice_head(n = N_INPUT_ONLY) %>%
  pull(acronym)

## --- output-only candidates: top anterograde regions not already reciprocal ---
## (source table has one row per channel -- collapse to best/min rank per region first)
output_ranked <- read.csv("~/Desktop/Bachelor/top_regions_low_deviation.csv", stringsAsFactors = FALSE) %>%
  filter(!(region_acronym %in% reciprocal)) %>%
  group_by(region_acronym) %>%
  summarise(combined = min(combined), .groups = "drop") %>%
  arrange(combined) %>%
  slice_head(n = N_OUTPUT_ONLY) %>%
  pull(region_acronym)

cat("Reciprocal (input+output):", paste(reciprocal, collapse = ", "), "\n")
cat("Input-only:", paste(input_ranked, collapse = ", "), "\n")
cat("Output-only:", paste(output_ranked, collapse = ", "), "\n")

## --- build edge list: FS <-> each region ---
edges <- bind_rows(
  data.frame(from = input_ranked, to = "FS", type = "Input only (retrograde)"),
  data.frame(from = "FS", to = output_ranked, type = "Output only (anterograde)"),
  data.frame(from = reciprocal, to = "FS", type = "Reciprocal (both)"),
  data.frame(from = "FS", to = reciprocal, type = "Reciprocal (both)")
)

all_nodes <- unique(c("FS", input_ranked, output_ranked, reciprocal))
node_role <- case_when(
  all_nodes == "FS" ~ "Injection site (FS)",
  all_nodes %in% reciprocal ~ "Reciprocal (both)",
  all_nodes %in% input_ranked ~ "Input only (retrograde)",
  TRUE ~ "Output only (anterograde)"
)
nodes <- data.frame(name = all_nodes, role = node_role)

## --- manual circular layout: FS centred, inputs on the left arc, outputs
## on the right arc, reciprocal regions on a top/bottom arc ---
n_in <- length(input_ranked); n_out <- length(output_ranked); n_rec <- length(reciprocal)
layout_df <- data.frame(name = all_nodes, x = 0, y = 0)
layout_df$x[layout_df$name == "FS"] <- 0
layout_df$y[layout_df$name == "FS"] <- 0

angle_in  <- seq(100, 260, length.out = n_in)
angle_out <- seq(-80, 80, length.out = n_out)
angle_rec_top <- seq(10, 80, length.out = ceiling(n_rec/2))
angle_rec_bot <- seq(190, 260, length.out = floor(n_rec/2))
angle_rec <- c(angle_rec_top, angle_rec_bot)

r_in <- 3; r_out <- 3; r_rec <- 4.6
for (i in seq_along(input_ranked)) {
  layout_df$x[layout_df$name == input_ranked[i]] <- r_in * cos(angle_in[i] * pi/180)
  layout_df$y[layout_df$name == input_ranked[i]] <- r_in * sin(angle_in[i] * pi/180)
}
for (i in seq_along(output_ranked)) {
  layout_df$x[layout_df$name == output_ranked[i]] <- r_out * cos(angle_out[i] * pi/180)
  layout_df$y[layout_df$name == output_ranked[i]] <- r_out * sin(angle_out[i] * pi/180)
}
for (i in seq_along(reciprocal)) {
  layout_df$x[layout_df$name == reciprocal[i]] <- r_rec * cos(angle_rec[i] * pi/180)
  layout_df$y[layout_df$name == reciprocal[i]] <- r_rec * sin(angle_rec[i] * pi/180)
}

g <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
layout_manual <- create_layout(g, layout = "manual", x = layout_df$x[match(nodes$name, layout_df$name)],
                                 y = layout_df$y[match(nodes$name, layout_df$name)])

edge_cols <- c("Input only (retrograde)" = "#2E7D32", "Output only (anterograde)" = "#8E24AA",
               "Reciprocal (both)" = "#C2185B")
node_cols <- c("Injection site (FS)" = "black", "Input only (retrograde)" = "#2E7D32",
               "Output only (anterograde)" = "#8E24AA", "Reciprocal (both)" = "#C2185B")

p <- ggraph(layout_manual) +
  geom_edge_arc(aes(color = type), curvature = 0.15, strength = 0.15,
                arrow = arrow(length = unit(3.2, "mm"), type = "closed"),
                start_cap = circle(4.5, "mm"), end_cap = circle(4.5, "mm"), linewidth = 0.8, alpha = 0.85) +
  geom_node_point(aes(color = role, size = ifelse(role == "Injection site (FS)", 9, 5.5))) +
  geom_node_text(aes(label = name), repel = TRUE, size = 5, fontface = "bold", max.overlaps = Inf,
                  color = "black", bg.color = "white", bg.r = 0.12) +
  scale_edge_color_manual(values = edge_cols, name = "Connection type") +
  scale_color_manual(values = node_cols, name = "Region type") +
  scale_size_identity() +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  labs(
    title = "Circuit diagram: FS injection site -- input (retrograde) and output (anterograde) regions",
    subtitle = "Arrows into FS = retrograde-traced afferents; arrows out of FS = anterograde-traced efferents; reciprocal regions show both"
  ) +
  theme_void(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold", size = 17, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey30"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 30, 10, 30)
  ) +
  guides(color = guide_legend(override.aes = list(size = 5), title.position = "top"),
         edge_color = guide_legend(title.position = "top"))

ggsave("fs_circuit_diagram.png", p, width = 14, height = 12, dpi = 300)
cat("\nSaved: fs_circuit_diagram.png\n")


# ---------------------------------------------------------------------------
# fs_circuit_diagram_by_condition
# ---------------------------------------------------------------------------
# Condition-matched circuit diagrams: unlike fs_circuit_diagram.R (which
# pools RET1+RET3 retrograde input and all-condition anterograde output),
# this builds ONE diagram per retrograde cohort, pairing that cohort's
# retrograde input with the ANTEROGRADE OUTPUT DATA FROM THE SAME
# STRESSOR CONDITION ONLY:
#   RET1 (retrograde input, acute restraint, Female, n=3)
#     <-> anterograde output, acute restraint condition (column "3", sexes pooled)
#   RET3 (retrograde input, acute forced swim, Male, n=2)
#     <-> anterograde output, acute forced swim condition (column "1", sexes pooled)
# This is a stricter match than fs_circuit_diagram.R's pooled version:
# every edge in a given diagram now comes from the SAME stressor exposure,
# not just "signal detected somewhere across all conditions".

library(dplyr)
library(igraph)
library(ggraph)
library(tidygraph)
library(ggplot2)
library(ggrepel)

setwd("~/Desktop/Bachelor/retrograde")
source("retrograde_common.R")

N_INPUT_ONLY  <- 8
N_OUTPUT_ONLY <- 8
N_RECIPROCAL  <- 6

anterograde_by_cond <- read.csv("~/Desktop/Bachelor/condition_region_intensity_for_brainrender.csv", stringsAsFactors = FALSE)

build_cohort_diagram <- function(retro_df, cohort_label, anterograde_col, out_file, title) {
  retro_ranked <- retro_df %>%
    filter(!is.na(acronym), Load > 0) %>%
    group_by(acronym) %>%
    summarise(mean_input_load = mean(Load, na.rm = TRUE), .groups = "drop")

  antero_ranked <- anterograde_by_cond %>%
    transmute(acronym = region_acronym, output_intensity = .data[[anterograde_col]]) %>%
    filter(output_intensity > 0)

  overlap <- inner_join(antero_ranked, retro_ranked, by = "acronym") %>%
    mutate(output_rank = rank(-output_intensity), input_rank = rank(-mean_input_load),
           reciprocal_rank = output_rank + input_rank) %>%
    arrange(reciprocal_rank)
  reciprocal <- head(overlap$acronym, N_RECIPROCAL)

  input_only <- retro_ranked %>% filter(!(acronym %in% reciprocal)) %>%
    arrange(desc(mean_input_load)) %>% slice_head(n = N_INPUT_ONLY) %>% pull(acronym)
  output_only <- antero_ranked %>% filter(!(acronym %in% reciprocal)) %>%
    arrange(desc(output_intensity)) %>% slice_head(n = N_OUTPUT_ONLY) %>% pull(acronym)

  cat(sprintf("\n=== %s ===\n", cohort_label))
  cat("Reciprocal:", paste(reciprocal, collapse = ", "), "\n")
  cat("Input-only:", paste(input_only, collapse = ", "), "\n")
  cat("Output-only:", paste(output_only, collapse = ", "), "\n")

  edges <- bind_rows(
    data.frame(from = input_only, to = "FS", type = "Input only (retrograde)"),
    data.frame(from = "FS", to = output_only, type = "Output only (anterograde)"),
    data.frame(from = reciprocal, to = "FS", type = "Reciprocal (both)"),
    data.frame(from = "FS", to = reciprocal, type = "Reciprocal (both)")
  )
  all_nodes <- unique(c("FS", input_only, output_only, reciprocal))
  node_role <- case_when(
    all_nodes == "FS" ~ "Injection site (FS)",
    all_nodes %in% reciprocal ~ "Reciprocal (both)",
    all_nodes %in% input_only ~ "Input only (retrograde)",
    TRUE ~ "Output only (anterograde)"
  )
  nodes <- data.frame(name = all_nodes, role = node_role)

  n_in <- length(input_only); n_out <- length(output_only); n_rec <- length(reciprocal)
  layout_df <- data.frame(name = all_nodes, x = 0, y = 0)
  angle_in  <- seq(100, 260, length.out = n_in)
  angle_out <- seq(-80, 80, length.out = n_out)
  angle_rec_top <- seq(10, 80, length.out = ceiling(n_rec/2))
  angle_rec_bot <- seq(190, 260, length.out = floor(n_rec/2))
  angle_rec <- c(angle_rec_top, angle_rec_bot)
  r_in <- 3; r_out <- 3; r_rec <- 4.6
  for (i in seq_along(input_only)) {
    layout_df$x[layout_df$name == input_only[i]] <- r_in * cos(angle_in[i]*pi/180)
    layout_df$y[layout_df$name == input_only[i]] <- r_in * sin(angle_in[i]*pi/180)
  }
  for (i in seq_along(output_only)) {
    layout_df$x[layout_df$name == output_only[i]] <- r_out * cos(angle_out[i]*pi/180)
    layout_df$y[layout_df$name == output_only[i]] <- r_out * sin(angle_out[i]*pi/180)
  }
  for (i in seq_along(reciprocal)) {
    layout_df$x[layout_df$name == reciprocal[i]] <- r_rec * cos(angle_rec[i]*pi/180)
    layout_df$y[layout_df$name == reciprocal[i]] <- r_rec * sin(angle_rec[i]*pi/180)
  }

  g <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
  layout_manual <- create_layout(g, layout = "manual",
                                   x = layout_df$x[match(nodes$name, layout_df$name)],
                                   y = layout_df$y[match(nodes$name, layout_df$name)])

  edge_cols <- c("Input only (retrograde)" = "#2E7D32", "Output only (anterograde)" = "#8E24AA",
                 "Reciprocal (both)" = "#C2185B")
  node_cols <- c("Injection site (FS)" = "black", "Input only (retrograde)" = "#2E7D32",
                 "Output only (anterograde)" = "#8E24AA", "Reciprocal (both)" = "#C2185B")

  p <- ggraph(layout_manual) +
    geom_edge_arc(aes(color = type), strength = 0.15,
                  arrow = arrow(length = unit(3.2, "mm"), type = "closed"),
                  start_cap = circle(4.5, "mm"), end_cap = circle(4.5, "mm"), linewidth = 0.9, alpha = 0.85) +
    geom_node_point(aes(color = role, size = ifelse(role == "Injection site (FS)", 9, 5.5))) +
    geom_node_text(aes(label = name), repel = TRUE, size = 5.4, fontface = "bold", max.overlaps = Inf,
                    color = "black", bg.color = "white", bg.r = 0.12) +
    scale_edge_color_manual(values = edge_cols, name = "Connection type") +
    scale_color_manual(values = node_cols, name = "Region type") +
    scale_size_identity() +
    scale_x_continuous(expand = expansion(mult = 0.15)) +
    scale_y_continuous(expand = expansion(mult = 0.15)) +
    labs(title = title,
         subtitle = "Arrows into FS = retrograde input; arrows out of FS = anterograde output; both from the SAME stressor condition") +
    theme_void(base_size = 15) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey30"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(10, 30, 10, 30)
    ) +
    guides(color = guide_legend(override.aes = list(size = 5), title.position = "top"),
           edge_color = guide_legend(title.position = "top"))

  ggsave(out_file, p, width = 12, height = 10.5, dpi = 300)
  cat("Saved:", out_file, "\n")
  list(reciprocal = reciprocal, input_only = input_only, output_only = output_only)
}

ret1 <- load_retrograde("RET1_RefAtlasRegions_combined.csv", exclude_fiber = TRUE) %>% filter(animal_id != "RET1_4")
ret3 <- load_retrograde("RET3_RefAtlasRegions_combined.csv", exclude_fiber = TRUE)

res_ret1 <- build_cohort_diagram(ret1, "RET1: acute restraint, Female", "X3",
                                   "fs_circuit_diagram_RET1_restraint.png",
                                   "FS circuit: acute restraint, Female (RET1, n=3)")
res_ret3 <- build_cohort_diagram(ret3, "RET3: acute forced swim, Male", "X1",
                                   "fs_circuit_diagram_RET3_forcedswim.png",
                                   "FS circuit: acute forced swim, Male (RET3, n=2)")

cat("\n=== Reciprocal regions shared between the two condition-matched diagrams ===\n")
print(intersect(res_ret1$reciprocal, res_ret3$reciprocal))


# ---------------------------------------------------------------------------
# fs_circuit_diagram_common
# ---------------------------------------------------------------------------
# "Common network" (Section 3.6): ONE combined circuit diagram showing
# both retrograde cohorts (RET1 restraint-Female + RET3 forced-swim-Male,
# pooled) against anterograde OUTPUT data restricted to ONLY the two
# matching conditions (acute restraint + acute forced swim), NOT all 5
# pooled conditions -- so both sides of the comparison are drawn from the
# same two stressor exposures, a direct, condition-matched comparison
# rather than diluting the output side with unrelated conditions (social
# defeat, tail suspension, chronic forced swim) that have no retrograde
# counterpart at all.

library(dplyr)
library(igraph)
library(ggraph)
library(tidygraph)
library(ggplot2)
library(ggrepel)

setwd("~/Desktop/Bachelor/retrograde")
source("retrograde_common.R")

N_INPUT_ONLY  <- 8
N_OUTPUT_ONLY <- 8
N_RECIPROCAL  <- 6

ret1 <- load_retrograde("RET1_RefAtlasRegions_combined.csv", exclude_fiber = TRUE) %>% filter(animal_id != "RET1_4")
ret3 <- load_retrograde("RET3_RefAtlasRegions_combined.csv", exclude_fiber = TRUE)

## input: RET1 + RET3 pooled (mean Load per region across both cohorts)
input_ranked <- bind_rows(ret1, ret3) %>%
  filter(!is.na(acronym), Load > 0) %>%
  group_by(acronym) %>%
  summarise(mean_input_load = mean(Load, na.rm = TRUE), .groups = "drop")

## output: anterograde intensity, restricted to ONLY restraint (X3) and
## forced swim acute (X1) -- the two conditions with a retrograde match --
## averaged across just those two, not all five.
antero <- read.csv("~/Desktop/Bachelor/condition_region_intensity_for_brainrender.csv", stringsAsFactors = FALSE)
output_ranked <- antero %>%
  transmute(acronym = region_acronym, output_intensity = (X1 + X3) / 2) %>%
  filter(output_intensity > 0)

overlap <- inner_join(output_ranked, input_ranked, by = "acronym") %>%
  mutate(output_rank = rank(-output_intensity), input_rank = rank(-mean_input_load),
         reciprocal_rank = output_rank + input_rank) %>%
  arrange(reciprocal_rank)
reciprocal <- head(overlap$acronym, N_RECIPROCAL)

input_only <- input_ranked %>% filter(!(acronym %in% reciprocal)) %>%
  arrange(desc(mean_input_load)) %>% slice_head(n = N_INPUT_ONLY) %>% pull(acronym)
output_only <- output_ranked %>% filter(!(acronym %in% reciprocal)) %>%
  arrange(desc(output_intensity)) %>% slice_head(n = N_OUTPUT_ONLY) %>% pull(acronym)

cat("=== Common network (input: RET1+RET3 pooled; output: restraint+forced swim only) ===\n")
cat("Reciprocal:", paste(reciprocal, collapse = ", "), "\n")
cat("Input-only:", paste(input_only, collapse = ", "), "\n")
cat("Output-only:", paste(output_only, collapse = ", "), "\n")

write.csv(overlap, "fs_common_network_overlap.csv", row.names = FALSE)

edges <- bind_rows(
  data.frame(from = input_only, to = "FS", type = "Input only (retrograde)"),
  data.frame(from = "FS", to = output_only, type = "Output only (anterograde)"),
  data.frame(from = reciprocal, to = "FS", type = "Reciprocal (both)"),
  data.frame(from = "FS", to = reciprocal, type = "Reciprocal (both)")
)
all_nodes <- unique(c("FS", input_only, output_only, reciprocal))
node_role <- case_when(
  all_nodes == "FS" ~ "Injection site (FS)",
  all_nodes %in% reciprocal ~ "Reciprocal (both)",
  all_nodes %in% input_only ~ "Input only (retrograde)",
  TRUE ~ "Output only (anterograde)"
)
nodes <- data.frame(name = all_nodes, role = node_role)

n_in <- length(input_only); n_out <- length(output_only); n_rec <- length(reciprocal)
layout_df <- data.frame(name = all_nodes, x = 0, y = 0)
angle_in  <- seq(100, 260, length.out = n_in)
angle_out <- seq(-80, 80, length.out = n_out)
angle_rec_top <- seq(10, 80, length.out = ceiling(n_rec/2))
angle_rec_bot <- seq(190, 260, length.out = floor(n_rec/2))
angle_rec <- c(angle_rec_top, angle_rec_bot)
r_in <- 3; r_out <- 3; r_rec <- 4.6
for (i in seq_along(input_only)) {
  layout_df$x[layout_df$name == input_only[i]] <- r_in * cos(angle_in[i]*pi/180)
  layout_df$y[layout_df$name == input_only[i]] <- r_in * sin(angle_in[i]*pi/180)
}
for (i in seq_along(output_only)) {
  layout_df$x[layout_df$name == output_only[i]] <- r_out * cos(angle_out[i]*pi/180)
  layout_df$y[layout_df$name == output_only[i]] <- r_out * sin(angle_out[i]*pi/180)
}
for (i in seq_along(reciprocal)) {
  layout_df$x[layout_df$name == reciprocal[i]] <- r_rec * cos(angle_rec[i]*pi/180)
  layout_df$y[layout_df$name == reciprocal[i]] <- r_rec * sin(angle_rec[i]*pi/180)
}

g <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
layout_manual <- create_layout(g, layout = "manual",
                                 x = layout_df$x[match(nodes$name, layout_df$name)],
                                 y = layout_df$y[match(nodes$name, layout_df$name)])

edge_cols <- c("Input only (retrograde)" = "#2E7D32", "Output only (anterograde)" = "#8E24AA",
               "Reciprocal (both)" = "#C2185B")
node_cols <- c("Injection site (FS)" = "black", "Input only (retrograde)" = "#2E7D32",
               "Output only (anterograde)" = "#8E24AA", "Reciprocal (both)" = "#C2185B")

p <- ggraph(layout_manual) +
  geom_edge_arc(aes(color = type), strength = 0.15,
                arrow = arrow(length = unit(3.2, "mm"), type = "closed"),
                start_cap = circle(4.5, "mm"), end_cap = circle(4.5, "mm"), linewidth = 0.9, alpha = 0.85) +
  geom_node_point(aes(color = role, size = ifelse(role == "Injection site (FS)", 9, 5.5))) +
  geom_node_text(aes(label = name), repel = TRUE, size = 28.5, fontface = "bold", max.overlaps = Inf,
                  color = "black", bg.color = "white", bg.r = 0.12) +
  scale_edge_color_manual(values = edge_cols, name = "Connection type") +
  scale_color_manual(values = node_cols, name = "Region type") +
  scale_size_identity() +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  labs(subtitle = "Input: RET1 (restraint,F) + RET3 (forced swim,M) pooled | Output: anterograde restraint + forced swim only (not all 5 conditions)") +
  theme_void(base_size = 78) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.subtitle = element_text(size = 63, hjust = 0.5, color = "grey30"),
    legend.title = element_text(size = 66, face = "bold"),
    legend.text = element_text(size = 60),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 30, 10, 30)
  ) +
  guides(color = guide_legend(override.aes = list(size = 12), title.position = "top"),
         edge_color = guide_legend(title.position = "top"))

ggsave("fs_circuit_diagram_common.png", p, width = 50, height = 43, dpi = 150, limitsize = FALSE)
cat("\nSaved: fs_circuit_diagram_common.png, fs_common_network_overlap.csv\n")

