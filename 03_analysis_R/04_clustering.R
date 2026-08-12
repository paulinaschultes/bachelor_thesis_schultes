# =============================================================================
# Hierarchical clustering, dendrograms, cluster agreement and stability
#
# Consolidated from the following original scripts:
#   - clustering_ari.R
#   - condition_clustering_agreement.R
#   - condition_clustering_jaccard.R
#   - condition_estimation_plots_and_cluster_stability.R
#   - condition_per_condition_cluster_stability.R
#   - condition_region_dendrograms_A4.R
#   - condition_region_dendrograms_by_condition_v2.R
#   - condition_sex_cluster_stability.R
#   - condition_sex_clustering_agreement.R
#   - condition_sex_region_clustering_zscore.R
#   - condition_sex_region_dendrograms_A4.R
#   - condition_sex_region_dendrograms_detailed_v2.R
#   - sex_cluster_stability.R
#   - sex_region_clustering_zscore.R
#   - sex_region_dendrograms_A4.R
#   - sex_region_dendrograms_by_sex_v2.R
# =============================================================================


# ---------------------------------------------------------------------------
# clustering_ari
# ---------------------------------------------------------------------------
# Unsupervised clustering + Adjusted Rand Index (ARI): complements the PCA
# (pca_overview.R) and PERMANOVA (permanova_test.R) by asking, very
# concretely: if an algorithm clusters the animal x channel samples into 5
# groups WITHOUT ever seeing the condition labels, how well does that
# recovered grouping match the TRUE stressor conditions? ARI = 0 means "no
# better than random", ARI = 1 means "perfect recovery". A permutation test
# (shuffling condition labels) gives a p-value for whether the observed ARI
# exceeds what pure chance would produce for a random 5-way partition of
# this size.
#
# Same input matrix as pca_overview.R / permanova_test.R (whole-brain-
# median-normalized, log10, condition-group-mean-imputed). Hierarchical
# clustering (Ward's method D2) on Euclidean distance, cut into k = 5.

library(dplyr)
library(tidyr)
library(mclust)   # for adjustedRandIndex
library(ggplot2)

df <- read.csv("~/Desktop/ALL_DT.csv", stringsAsFactors = FALSE)
cond_map <- read.csv("~/Desktop/group_condition_corrected.csv", stringsAsFactors = FALSE)

long <- df %>%
  filter(region_acronym != "FS") %>%
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

long <- long %>%
  left_join(animal_overall, by = c("animal_id_original", "Channel" = "channel")) %>%
  mutate(log_value = log10(value / overall_median + 0.001))

wide <- long %>% select(sample_id, Condition, region_acronym, log_value) %>%
  pivot_wider(names_from = region_acronym, values_from = log_value)

sample_meta <- wide %>% select(sample_id, Condition)
mat <- as.matrix(wide %>% select(-sample_id, -Condition))
rownames(mat) <- wide$sample_id

for (j in seq_len(ncol(mat))) {
  col <- mat[, j]
  if (any(is.na(col))) {
    for (cond in unique(sample_meta$Condition)) {
      idx <- which(sample_meta$Condition == cond)
      cond_mean <- mean(col[idx], na.rm = TRUE)
      if (is.nan(cond_mean)) cond_mean <- mean(col, na.rm = TRUE)
      col[idx][is.na(col[idx])] <- cond_mean
    }
    mat[, j] <- col
  }
}
mat <- mat[, colSums(is.na(mat)) == 0, drop = FALSE]
meta <- sample_meta[match(rownames(mat), sample_meta$sample_id), ]
cat("Clustering input:", nrow(mat), "samples x", ncol(mat), "regions\n")

## --- Hierarchical clustering (Ward D2, Euclidean), cut into k = 5 ---
d <- dist(mat, method = "euclidean")
hc <- hclust(d, method = "ward.D2")
k <- length(unique(meta$Condition))
clusters <- cutree(hc, k = k)

observed_ari <- adjustedRandIndex(clusters, meta$Condition)
cat(sprintf("\n=== Observed Adjusted Rand Index (hierarchical clustering vs true Condition) ===\n"))
cat(sprintf("ARI = %.4f (0 = chance level, 1 = perfect recovery)\n", observed_ari))

## --- Permutation test: shuffle condition labels, recompute ARI, build null ---
set.seed(42)
n_perm <- 9999
null_ari <- numeric(n_perm)
for (i in seq_len(n_perm)) {
  shuffled <- sample(meta$Condition)
  null_ari[i] <- adjustedRandIndex(clusters, shuffled)
}
p_value <- (sum(null_ari >= observed_ari) + 1) / (n_perm + 1)

cat(sprintf("Permutation null: mean = %.4f, 95th percentile = %.4f\n", mean(null_ari), quantile(null_ari, 0.95)))
cat(sprintf("Permutation p-value: %.4f\n", p_value))

write.csv(data.frame(observed_ari = observed_ari, null_mean = mean(null_ari),
                      null_95pct = quantile(null_ari, 0.95), p_value = p_value),
          "~/Desktop/clustering_ari_result.csv", row.names = FALSE)

## --- Confusion-style table: cluster assignment vs true condition ---
tab <- table(Cluster = clusters, TrueCondition = meta$Condition)
cat("\n=== Cluster assignment vs. true condition ===\n")
print(tab)
write.csv(as.data.frame.matrix(tab), "~/Desktop/clustering_confusion_table.csv")

## --- Plot 1: null distribution with observed ARI marked ---
p1 <- ggplot(data.frame(ari = null_ari), aes(x = ari)) +
  geom_histogram(bins = 60, fill = "#F48FB1", alpha = 0.85) +
  geom_vline(xintercept = observed_ari, color = "#4A148C", linewidth = 1) +
  annotate("text", x = observed_ari, y = Inf, label = sprintf("observed ARI = %.3f\np = %.3f", observed_ari, p_value),
           hjust = ifelse(observed_ari > median(null_ari), 1.05, -0.05), vjust = 1.5, color = "#4A148C", size = 3.5) +
  labs(title = "Permutation null distribution of Adjusted Rand Index",
       subtitle = "Hierarchical clustering (k=5) vs. randomly shuffled condition labels (9999 permutations)",
       x = "Adjusted Rand Index", y = "Count") +
  theme_minimal(base_size = 11)
ggsave("~/Desktop/clustering_ari_permutation.png", p1, width = 8, height = 5.5, units = "in", dpi = 300)

## --- Plot 2: dendrogram colored by true condition ---
library(ggdendro)
cond_cols <- c(
  "forced swim · Acute"     = "#FBC9DE",
  "forced swim · Chronic"   = "#F48FB1",
  "restraint · Acute"       = "#E05780",
  "social defeat · Acute"   = "#C2185B",
  "tail suspension · Acute" = "#8E24AA"
)
dend_data <- dendro_data(hc)
label_meta <- meta[match(dend_data$labels$label, meta$sample_id), ]
dend_data$labels$Condition <- label_meta$Condition

p2 <- ggplot() +
  geom_segment(data = dend_data$segments, aes(x = x, y = y, xend = xend, yend = yend), color = "grey50", linewidth = 0.4) +
  geom_text(data = dend_data$labels, aes(x = x, y = -max(dend_data$segments$y) * 0.03, label = label, color = Condition),
            angle = 90, hjust = 1, size = 2.2) +
  scale_color_manual(values = cond_cols) +
  labs(title = "Hierarchical clustering (Ward D2) of whole-brain profiles",
       subtitle = sprintf("Leaf color = true stressor condition | ARI vs true labels = %.3f (p = %.3f)", observed_ari, p_value),
       x = NULL, y = "Height") +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "bottom", plot.margin = margin(10, 10, 60, 10))
ggsave("~/Desktop/clustering_dendrogram.png", p2, width = 12, height = 6.5, units = "in", dpi = 300)

cat("\nSaved:\n  ~/Desktop/clustering_ari_result.csv\n  ~/Desktop/clustering_confusion_table.csv\n")
cat("  ~/Desktop/clustering_ari_permutation.png\n  ~/Desktop/clustering_dendrogram.png\n")


# ---------------------------------------------------------------------------
# condition_clustering_agreement
# ---------------------------------------------------------------------------
# Agreement between the 5 per-condition hierarchical clusterings
# (condition_region_clusters_within_condition.csv, built by
# condition_matrix_stats.R section 3b -- top 40 highest-intensity regions
# per condition, clustered into k=5 groups). Each pair of conditions shares
# ~27-33 of their 40 top regions (largely the same high-baseline regions,
# e.g. striatum/cerebellum, regardless of condition). For each pair, the
# Adjusted Rand Index (ARI) is computed on the shared regions' cluster
# labels -- ARI = 0 means no better than chance agreement, ARI = 1 means
# identical clustering.

library(dplyr)
library(mclust)
library(ggplot2)
library(reshape2)

setwd("~/Desktop/Bachelor")
df <- read.csv("condition_region_clusters_within_condition.csv")
conditions <- unique(df$condition)
n <- length(conditions)

ari_mat <- matrix(NA, n, n, dimnames = list(conditions, conditions))
n_common_mat <- matrix(NA, n, n, dimnames = list(conditions, conditions))

for (i in seq_len(n)) {
  for (j in seq_len(n)) {
    if (i == j) { ari_mat[i, j] <- 1; next }
    d1 <- df %>% filter(condition == conditions[i])
    d2 <- df %>% filter(condition == conditions[j])
    common <- intersect(d1$region_acronym, d2$region_acronym)
    c1 <- d1$cluster[match(common, d1$region_acronym)]
    c2 <- d2$cluster[match(common, d2$region_acronym)]
    ari_mat[i, j] <- adjustedRandIndex(c1, c2)
    n_common_mat[i, j] <- length(common)
  }
}

cat("=== Adjusted Rand Index between per-condition clusterings (on shared top-40 regions) ===\n")
print(round(ari_mat, 3))
write.csv(as.data.frame(ari_mat), "condition_clustering_agreement_ari.csv")

cat("\nMean off-diagonal ARI:", round(mean(ari_mat[upper.tri(ari_mat)]), 3), "\n")

## --- Plot: ARI heatmap ---
ari_long <- melt(ari_mat, varnames = c("Condition1", "Condition2"), value.name = "ARI")

FONT_PT <- 12
mm_size <- function(pt) pt / (72.27 / 25.4)  # ggplot geom text 'size' aesthetic is in mm

p <- ggplot(ari_long, aes(x = Condition2, y = Condition1, fill = ARI)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", ARI), color = ARI > 0.55), size = mm_size(FONT_PT), show.legend = FALSE) +
  scale_fill_gradientn(colors = c("#F6EEF9", "#DCC6E8", "#B98CD1", "#9558B8", "#6E3A93", "#4A148C"),
                        limits = c(-0.05, 1)) +
  scale_color_manual(values = c("TRUE" = "white", "FALSE" = "#3A3A3A")) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = FONT_PT) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1, size = FONT_PT),
    axis.text.y = element_text(size = FONT_PT)
  )

ggsave("condition_clustering_agreement_ari.png", p, width = 14, height = 11, units = "in", dpi = 300)
cat("\nSaved: condition_clustering_agreement_ari.csv/.png\n")


# ---------------------------------------------------------------------------
# condition_clustering_jaccard
# ---------------------------------------------------------------------------
# Jaccard index between the 5 per-condition top-40-region sets
# (condition_region_clusters_within_condition.csv, same source used by
# condition_clustering_agreement.R for the ARI comparison). While the ARI
# measures whether the SHARED regions get grouped the same way (partition
# agreement), the Jaccard index instead measures how much the underlying
# region SETS themselves overlap, independent of any cluster labels:
# J(A,B) = |A intersect B| / |A union B|.

library(dplyr)
library(ggplot2)
library(reshape2)

setwd("~/Desktop/Bachelor")
df <- read.csv("condition_region_clusters_within_condition.csv")
conditions <- unique(df$condition)
n <- length(conditions)

region_sets <- lapply(conditions, function(cond) df$region_acronym[df$condition == cond])
names(region_sets) <- conditions

jaccard_mat <- matrix(NA, n, n, dimnames = list(conditions, conditions))
for (i in seq_len(n)) {
  for (j in seq_len(n)) {
    a <- region_sets[[i]]; b <- region_sets[[j]]
    jaccard_mat[i, j] <- length(intersect(a, b)) / length(union(a, b))
  }
}

cat("=== Jaccard index between per-condition top-40 region sets ===\n")
print(round(jaccard_mat, 3))
write.csv(as.data.frame(jaccard_mat), "condition_clustering_agreement_jaccard.csv")

cat("\nMean off-diagonal Jaccard:", round(mean(jaccard_mat[upper.tri(jaccard_mat)]), 3), "\n")

## --- Plot: Jaccard heatmap (same pink-violet style as the ARI heatmap) ---
jaccard_long <- melt(jaccard_mat, varnames = c("Condition1", "Condition2"), value.name = "Jaccard")

FONT_PT <- 12
mm_size <- function(pt) pt / (72.27 / 25.4)  # ggplot geom text 'size' aesthetic is in mm

p <- ggplot(jaccard_long, aes(x = Condition2, y = Condition1, fill = Jaccard)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", Jaccard), color = Jaccard > 0.55), size = mm_size(FONT_PT), show.legend = FALSE) +
  scale_fill_gradientn(colors = c("#F6EEF9", "#DCC6E8", "#B98CD1", "#9558B8", "#6E3A93", "#4A148C"),
                        limits = c(0, 1)) +
  scale_color_manual(values = c("TRUE" = "white", "FALSE" = "#3A3A3A")) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = FONT_PT) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1, size = FONT_PT),
    axis.text.y = element_text(size = FONT_PT)
  )

ggsave("condition_clustering_agreement_jaccard.png", p, width = 14, height = 11, units = "in", dpi = 300)
cat("\nSaved: condition_clustering_agreement_jaccard.csv/.png\n")


# ---------------------------------------------------------------------------
# condition_estimation_plots_and_cluster_stability
# ---------------------------------------------------------------------------
# Two analyses building on the condition Kruskal-Wallis (condition_matrix_stats.R)
# and the hierarchical clustering of conditions (section 3 there, which found
# two clusters: {forced swim Acute/Chronic, restraint} vs {social defeat,
# tail suspension}):
#
#   1. Estimation plots (Gardner-Altman, dabestr): for the 18 regions with
#      nominal region-level significance (raw p < 0.05 in the condition
#      Kruskal-Wallis), show the actual EFFECT SIZE (bootstrap mean
#      difference + 95% BCa CI) between the two hierarchically-derived
#      condition clusters, pooling all animals/channels within each
#      cluster. This answers "how big are these effects", complementing
#      the p-value-only view.
#   2. Cluster stability ("how much would the clustering deviate under
#      resampling"): regions are bootstrap-resampled (n=1000), the
#      condition clustering (Ward D2 on z-scored profiles, cut at k=2) is
#      recomputed each time, and for every pair of conditions we record how
#      often they land in the same cluster -- a standard bootstrap support
#      measure for hierarchical clustering (same logic as pvclust).

library(dplyr)
library(tidyr)
library(dabestr)
library(ggplot2)
library(patchwork)

setwd("~/Desktop/Bachelor")
source("prepare_heatmap_data.R")   # mat, long, condition_order, major_division_per_row

keep <- rownames(mat) != "FS"
mat <- mat[keep, , drop = FALSE]
long <- long %>% filter(region_acronym != "FS")

## ============================================================
## 1. Estimation plots: Cluster A (FS acute/chronic + restraint) vs
##    Cluster B (social defeat + tail suspension), for the 18 nominally
##    significant regions from the condition Kruskal-Wallis test.
## ============================================================
cluster_A <- c("forced swim · Acute", "forced swim · Chronic", "restraint · Acute")
cluster_B <- c("social defeat · Acute", "tail suspension · Acute")

kw <- read.csv("condition_kruskal_by_region.csv") %>% arrange(p_value)
top_regions <- head(kw$region_acronym[kw$p_value < 0.05], 18)

long_cl <- long %>%
  mutate(cluster_group = case_when(
    condition %in% cluster_A ~ "Cluster A",
    condition %in% cluster_B ~ "Cluster B"
  )) %>%
  filter(!is.na(cluster_group))

## Cluster A = forced swim (acute, chronic) + acute restraint
## Cluster B = acute social defeat + acute tail suspension
group_levels <- c("Cluster A", "Cluster B")

plots <- list()
effect_results <- vector("list", length(top_regions))

for (i in seq_along(top_regions)) {
  reg <- top_regions[i]
  d <- long_cl %>% filter(region_acronym == reg)
  db <- d %>% dabestr::load(x = cluster_group, y = intensity, idx = group_levels)
  db_diff <- mean_diff(db)

  boot <- db_diff$boot_result
  effect_results[[i]] <- data.frame(
    region_acronym = reg,
    mean_diff_B_minus_A = boot$difference,
    bca_ci_low = boot$bca_ci_low,
    bca_ci_high = boot$bca_ci_high,
    permtest_p = db_diff$permtest_pvals$pval_permtest
  )

  p <- dabest_plot(db_diff, TRUE, raw_marker_size = 0.9, es_marker_size = 1.2) +
    ggtitle(reg) + theme(plot.title = element_text(size = 9, face = "bold"))
  plots[[reg]] <- p
}

effect_df <- bind_rows(effect_results) %>% arrange(permtest_p)
write.csv(effect_df, "condition_cluster_effect_sizes.csv", row.names = FALSE)

cat("=== 1. Estimation plots: Cluster A vs Cluster B, 18 top regions ===\n")
print(effect_df)

combined <- wrap_plots(plots, ncol = 3) +
  plot_annotation(
    caption = "Cluster A = forced swim (acute, chronic) + acute restraint  |  Cluster B = acute social defeat + acute tail suspension",
    theme = theme(plot.caption = element_text(size = 10, hjust = 0.5))
  )
ggsave("condition_cluster_estimation_plots.png", combined,
       width = 15, height = ceiling(length(top_regions) / 3) * 4.5, units = "in", dpi = 300, limitsize = FALSE)
cat("\nSaved: condition_cluster_estimation_plots.png, condition_cluster_effect_sizes.csv\n")

## ============================================================
## 2. Bootstrap stability of the condition hierarchical clustering
## ============================================================
mat_z <- t(scale(t(mat)))
mat_z <- mat_z[rowSums(is.na(mat_z)) == 0, , drop = FALSE]

set.seed(42)
n_boot <- 1000
n_regions <- nrow(mat_z)
cocluster_count <- matrix(0, 5, 5, dimnames = list(condition_order, condition_order))

for (b in seq_len(n_boot)) {
  idx <- sample(seq_len(n_regions), n_regions, replace = TRUE)
  boot_mat <- mat_z[idx, , drop = FALSE]
  d <- dist(t(boot_mat), method = "euclidean")
  hc <- hclust(d, method = "ward.D2")
  cl <- cutree(hc, k = 2)
  same <- outer(cl, cl, `==`)
  cocluster_count <- cocluster_count + same
}

cocluster_prop <- cocluster_count / n_boot
write.csv(as.data.frame(cocluster_prop), "condition_cluster_bootstrap_stability.csv")

cat("\n=== 2. Bootstrap stability of condition clustering (k=2), n =", n_boot, "resamples ===\n")
cat("Proportion of resamples in which each condition PAIR fell in the same cluster:\n")
print(round(cocluster_prop, 3))

## Observed clustering on the full (non-resampled) data, for reference
hc_obs <- hclust(dist(t(mat_z), method = "euclidean"), method = "ward.D2")
cl_obs <- cutree(hc_obs, k = 2)
cat("\nObserved k=2 cluster assignment (full data):\n")
print(cl_obs)

cat("\nSaved: condition_cluster_bootstrap_stability.csv\n")


# ---------------------------------------------------------------------------
# condition_per_condition_cluster_stability
# ---------------------------------------------------------------------------
# Bootstrap stability of the per-condition region clusterings (Figure 10 /
# condition_matrix_stats.R section 3b): for each condition, animals are
# resampled with replacement (1000 iterations), the top-40-region clustering
# recomputed each time, and compared to the original (observed) clustering
# via Adjusted Rand Index (ARI). Unlike the condition-level 2-cluster
# stability test (which resampled REGIONS and was ~99.9-100% stable), this
# resamples ANIMALS -- the relevant unit of replication for a clustering
# that is itself built from animal-to-animal covariation -- and is
# expected to be noisier given the much smaller per-condition sample sizes
# (n = 6-18 animals).

library(dplyr)
library(tidyr)
library(mclust)

setwd("~/Desktop/Bachelor")
source("prepare_heatmap_data.R")
long <- long %>% filter(region_acronym != "FS")

top_n_regions <- 40
k_regions_per_cond <- 5
n_boot <- 1000
set.seed(42)

stability_results <- vector("list", length(condition_order))

for (cond in condition_order) {
  d <- long %>% filter(condition == cond) %>% select(region_acronym, animal_id_original, intensity)
  wide_cond <- d %>% pivot_wider(names_from = animal_id_original, values_from = intensity)
  mat_cond <- as.matrix(wide_cond[, -1])
  rownames(mat_cond) <- wide_cond$region_acronym
  mat_cond <- mat_cond[rowSums(is.na(mat_cond)) == 0, , drop = FALSE]

  med_per_region <- apply(mat_cond, 1, median, na.rm = TRUE)
  top_regions <- names(sort(med_per_region, decreasing = TRUE))[seq_len(min(top_n_regions, length(med_per_region)))]
  mat_cond_top <- mat_cond[top_regions, , drop = FALSE]

  ## Observed clustering (same as condition_matrix_stats.R section 3b)
  mat_z_obs <- t(scale(t(mat_cond_top)))
  mat_z_obs <- mat_z_obs[rowSums(is.na(mat_z_obs)) == 0, , drop = FALSE]
  hc_obs <- hclust(dist(mat_z_obs, method = "euclidean"), method = "ward.D2")
  k <- min(k_regions_per_cond, nrow(mat_z_obs) - 1)
  clusters_obs <- cutree(hc_obs, k = k)

  n_animals <- ncol(mat_cond_top)
  ari_vals <- numeric(n_boot)

  for (b in seq_len(n_boot)) {
    idx <- sample(seq_len(n_animals), n_animals, replace = TRUE)
    boot_mat <- mat_cond_top[, idx, drop = FALSE]
    ## regions with zero variance in this bootstrap draw can't be z-scored;
    ## match them to observed-set regions for a fair ARI comparison
    boot_z <- t(scale(t(boot_mat)))
    valid <- rowSums(is.na(boot_z)) == 0
    common_regions <- intersect(rownames(mat_z_obs), rownames(boot_z)[valid])
    if (length(common_regions) < k + 1) { ari_vals[b] <- NA; next }

    hc_boot <- hclust(dist(boot_z[common_regions, , drop = FALSE], method = "euclidean"), method = "ward.D2")
    clusters_boot <- cutree(hc_boot, k = k)

    ari_vals[b] <- adjustedRandIndex(clusters_obs[common_regions], clusters_boot[names(clusters_boot)])
  }

  stability_results[[cond]] <- data.frame(
    condition = cond, n_animals = n_animals, n_regions = nrow(mat_z_obs),
    mean_ari = mean(ari_vals, na.rm = TRUE), median_ari = median(ari_vals, na.rm = TRUE),
    sd_ari = sd(ari_vals, na.rm = TRUE),
    pct_ari_above_0.5 = mean(ari_vals > 0.5, na.rm = TRUE) * 100,
    n_valid_boots = sum(!is.na(ari_vals))
  )

  cat(sprintf("%s: n=%d animals, mean ARI = %.3f (sd %.3f), %.0f%% of resamples ARI>0.5\n",
              cond, n_animals, mean(ari_vals, na.rm = TRUE), sd(ari_vals, na.rm = TRUE),
              mean(ari_vals > 0.5, na.rm = TRUE) * 100))
}

stability_df <- bind_rows(stability_results) %>% arrange(desc(mean_ari))
write.csv(stability_df, "condition_per_condition_cluster_stability.csv", row.names = FALSE)

cat("\n=== Summary: per-condition clustering stability (bootstrap over animals, n=1000) ===\n")
print(stability_df)
cat("\nSaved: condition_per_condition_cluster_stability.csv\n")


# ---------------------------------------------------------------------------
# condition_region_dendrograms_A4
# ---------------------------------------------------------------------------
# Per-condition region dendrograms (Figure 9), laid out for a DIN A4 page
# (8.27 x 11.69 in), panels labeled (a)-(e) matching the thesis text order
# (a: forced swim Acute, b: forced swim Chronic, c: restraint Acute,
# d: social defeat Acute, e: tail suspension Acute). Same clustering method
# as condition_matrix_stats.R section 3b (top 40 regions, Ward D2,
# z-scored, k=5), with bootstrap ARI stability annotated per panel from
# condition_per_condition_cluster_stability.csv.

library(dplyr)
library(tidyr)
library(dendextend)

setwd("~/Desktop/Bachelor")
source("prepare_heatmap_data.R")   # long, major_division_per_row, mat, condition_order
long <- long %>% filter(region_acronym != "FS")

top_n_regions <- 40
k_regions_per_cond <- 5

major_cols <- c(
  "Isocortex"             = "#5B7C99", "Olfactory areas"       = "#8C6D46",
  "Hippocampal formation" = "#6B8F71", "Cortical subplate"     = "#9C6B85",
  "Striatum"              = "#4F8F8B", "Pallidum"              = "#7A6B9C",
  "Thalamus"              = "#B08B5A", "Hypothalamus"          = "#A85C5C",
  "Midbrain"              = "#6B6B99", "Pons"                  = "#8FA05E",
  "Medulla"               = "#5E8FA0", "Cerebellum"            = "#9C7A5B",
  "fiber tracts"          = "#A9A9A9", "ventricular systems"   = "#C4B7A6"
)

hc_by_condition <- vector("list", length(condition_order))
names(hc_by_condition) <- condition_order

for (cond in condition_order) {
  d <- long %>% filter(condition == cond) %>% select(region_acronym, animal_id_original, intensity)
  wide_cond <- d %>% pivot_wider(names_from = animal_id_original, values_from = intensity)
  mat_cond <- as.matrix(wide_cond[, -1])
  rownames(mat_cond) <- wide_cond$region_acronym
  mat_cond <- mat_cond[rowSums(is.na(mat_cond)) == 0, , drop = FALSE]

  med_per_region <- apply(mat_cond, 1, median, na.rm = TRUE)
  top_regions <- names(sort(med_per_region, decreasing = TRUE))[seq_len(min(top_n_regions, length(med_per_region)))]
  mat_cond_top <- mat_cond[top_regions, , drop = FALSE]

  mat_cond_z <- t(scale(t(mat_cond_top)))
  mat_cond_z <- mat_cond_z[rowSums(is.na(mat_cond_z)) == 0, , drop = FALSE]

  hc_cond <- hclust(dist(mat_cond_z, method = "euclidean"), method = "ward.D2")
  hc_cond$height <- hc_cond$height / max(hc_cond$height)  # unify height scale across conditions
  hc_by_condition[[cond]] <- hc_cond
}

## Bootstrap ARI stability (already computed, condition_per_condition_cluster_stability.R)
stability_path <- "condition_per_condition_cluster_stability.csv"
stability_df <- if (file.exists(stability_path)) read.csv(stability_path) else NULL

panel_letters <- setNames(letters[seq_along(condition_order)], condition_order)

## --- DIN A4 portrait: 8.27 x 11.69 in, 2 columns x 3 rows (5 conditions + legend) ---
draw_fig9 <- function() {
  par(mfrow = c(3, 2), oma = c(1, 0.5, 3, 0.5))
  for (cond in condition_order) {
    dend <- as.dendrogram(hc_by_condition[[cond]])
    leaf_divs <- major_division_per_row[match(labels(dend), rownames(mat))]
    labels_colors(dend) <- major_cols[as.character(leaf_divs)]

    stab_line <- ""
    if (!is.null(stability_df) && cond %in% stability_df$condition) {
      s <- stability_df[stability_df$condition == cond, ]
      stab_line <- sprintf("\nARI = %.2f +/- %.2f | ARI>0.5: %.0f%%", s$mean_ari, s$sd_ari, s$pct_ari_above_0.5)
    }

    par(mar = c(5, 4, 10, 9), cex = 1.6)
    plot(dend, horiz = TRUE,
         main = sprintf("(%s) %s\n\n(top %d regions)%s", panel_letters[cond], cond, top_n_regions, stab_line),
         xlab = "Height", cex.main = 2.2, cex.axis = 2.1, cex.lab = 2.3)
  }
  ## 6th grid cell: legend for the division color code
  par(mar = c(2, 0, 0, 0), cex = 1)
  plot.new()
  legend("center", legend = names(major_cols), fill = major_cols, cex = 2.4,
         bty = "n", title = "Major division", ncol = 1)
  mtext("Per-condition region dendrograms (top 40 regions, Ward D2, z-scored)",
        side = 3, outer = TRUE, cex = 2.5, font = 2, line = 1.2)
}

pdf("condition_region_dendrograms_A4.pdf", width = 24, height = 40)
draw_fig9()
dev.off()

png("condition_region_dendrograms_A4.png", width = 24, height = 40, units = "in", res = 150)
draw_fig9()
dev.off()
par(mfrow = c(1, 1), cex = 1)

cat("Saved: condition_region_dendrograms_A4.pdf\n")
cat("Panel labels: ", paste(sprintf("(%s)=%s", panel_letters, names(panel_letters)), collapse = "; "), "\n")


# ---------------------------------------------------------------------------
# condition_region_dendrograms_by_condition_v2
# ---------------------------------------------------------------------------
# Regenerates condition_region_dendrograms_by_condition.pdf (originally
# from condition_matrix_stats.R section 3b) with division-color leaf labels
# + a legend on EACH page, matching condition_region_dendrograms_A4.R's
# style, but one full page per condition instead of a combined grid.

library(dplyr)
library(tidyr)
library(dendextend)

setwd("~/Desktop/Bachelor")
source("prepare_heatmap_data.R")
long <- long %>% filter(region_acronym != "FS")

top_n_regions <- 40

major_cols <- c(
  "Isocortex"             = "#5B7C99", "Olfactory areas"       = "#8C6D46",
  "Hippocampal formation" = "#6B8F71", "Cortical subplate"     = "#9C6B85",
  "Striatum"              = "#4F8F8B", "Pallidum"              = "#7A6B9C",
  "Thalamus"              = "#B08B5A", "Hypothalamus"          = "#A85C5C",
  "Midbrain"              = "#6B6B99", "Pons"                  = "#8FA05E",
  "Medulla"               = "#5E8FA0", "Cerebellum"            = "#9C7A5B",
  "fiber tracts"          = "#A9A9A9", "ventricular systems"   = "#C4B7A6"
)

hc_by_condition <- vector("list", length(condition_order))
names(hc_by_condition) <- condition_order

for (cond in condition_order) {
  d <- long %>% filter(condition == cond) %>% select(region_acronym, animal_id_original, intensity)
  wide_cond <- d %>% pivot_wider(names_from = animal_id_original, values_from = intensity)
  mat_cond <- as.matrix(wide_cond[, -1])
  rownames(mat_cond) <- wide_cond$region_acronym
  mat_cond <- mat_cond[rowSums(is.na(mat_cond)) == 0, , drop = FALSE]

  med_per_region <- apply(mat_cond, 1, median, na.rm = TRUE)
  top_regions <- names(sort(med_per_region, decreasing = TRUE))[seq_len(min(top_n_regions, length(med_per_region)))]
  mat_cond_top <- mat_cond[top_regions, , drop = FALSE]

  mat_cond_z <- t(scale(t(mat_cond_top)))
  mat_cond_z <- mat_cond_z[rowSums(is.na(mat_cond_z)) == 0, , drop = FALSE]

  hc_cond <- hclust(dist(mat_cond_z, method = "euclidean"), method = "ward.D2")
  hc_cond$height <- hc_cond$height / max(hc_cond$height)  # unify height scale across conditions
  hc_by_condition[[cond]] <- hc_cond
}

pdf("condition_region_dendrograms_by_condition.pdf", width = 7, height = 9)
for (cond in condition_order) {
  dend <- as.dendrogram(hc_by_condition[[cond]])
  leaf_divs <- major_division_per_row[match(labels(dend), rownames(mat))]
  labels_colors(dend) <- major_cols[as.character(leaf_divs)]

  layout(matrix(c(1, 2), nrow = 2), heights = c(0.82, 0.18))
  par(mar = c(4, 2, 4.5, 10))
  plot(dend, horiz = TRUE, main = sprintf("%s\n\n(top %d regions by median intensity)", cond, top_n_regions),
       xlab = "Height", cex.main = 0.9)

  par(mar = c(0.5, 1, 1, 1))
  plot.new()
  legend("center", legend = names(major_cols), fill = major_cols, cex = 0.75,
         bty = "n", title = "Major division", ncol = 4)
}
dev.off()

cat("Saved (overwritten): condition_region_dendrograms_by_condition.pdf\n")


# ---------------------------------------------------------------------------
# condition_sex_cluster_stability
# ---------------------------------------------------------------------------
# Bootstrap stability of the condition x sex region clusterings
# (condition_sex_region_clustering_zscore.R), same method as
# condition_per_condition_cluster_stability.R / sex_cluster_stability.R:
# resample samples with replacement (1000 iterations) within each of the 8
# condition x sex cells, recompute the top-40(or fewer)-region clustering
# each time, compare to the observed clustering via Adjusted Rand Index.
# "forced swim · Chronic" excluded (all-Female, no Male cell), matching the
# 4x2 = 8-cell grid used throughout this analysis.

library(dplyr)
library(tidyr)
library(mclust)

setwd("~/Desktop/Bachelor")
source("prepare_heatmap_data.R")
long <- long %>% filter(region_acronym != "FS")

top_n_regions <- 40
k_regions <- 5
n_boot <- 1000
set.seed(42)
sex_order <- c("Female", "Male")
condition_order_no_chronicFS <- setdiff(condition_order, "social defeat · Acute")  # social defeat is now single-sex (Male only); excluded here so only conditions comparable across both sexes remain. Chronic FS is included since its sex labels were corrected.

stability_results <- list()

for (cond in condition_order_no_chronicFS) {
  for (sx in sex_order) {
    key <- paste(cond, sx, sep = " | ")
    d <- long %>% filter(condition == cond, sex == sx) %>%
      select(region_acronym, animal_id_original, channel, intensity) %>%
      mutate(sample_id = paste(animal_id_original, channel, sep = "_"))
    wide_cell <- d %>% select(region_acronym, sample_id, intensity) %>%
      pivot_wider(names_from = sample_id, values_from = intensity)
    mat_cell <- as.matrix(wide_cell[, -1]); rownames(mat_cell) <- wide_cell$region_acronym
    mat_cell <- mat_cell[rowSums(is.na(mat_cell)) == 0, , drop = FALSE]

    if (ncol(mat_cell) < 2) {
      cat(sprintf("Skipping %s: only %d sample(s), clustering needs >=2\n", key, ncol(mat_cell)))
      next
    }

    med_per_region <- apply(mat_cell, 1, median, na.rm = TRUE)
    n_top <- min(top_n_regions, length(med_per_region))
    top_regions <- names(sort(med_per_region, decreasing = TRUE))[seq_len(n_top)]
    mat_cell_top <- mat_cell[top_regions, , drop = FALSE]

    mat_z_obs <- t(scale(t(mat_cell_top)))
    mat_z_obs <- mat_z_obs[rowSums(is.na(mat_z_obs)) == 0, , drop = FALSE]
    hc_obs <- hclust(dist(mat_z_obs, method = "euclidean"), method = "ward.D2")
    k <- min(k_regions, nrow(mat_z_obs) - 1)
    clusters_obs <- cutree(hc_obs, k = k)

    n_samples <- ncol(mat_cell_top)
    ari_vals <- numeric(n_boot)

    for (b in seq_len(n_boot)) {
      idx <- sample(seq_len(n_samples), n_samples, replace = TRUE)
      boot_mat <- mat_cell_top[, idx, drop = FALSE]
      boot_z <- t(scale(t(boot_mat)))
      valid <- rowSums(is.na(boot_z)) == 0
      common_regions <- intersect(rownames(mat_z_obs), rownames(boot_z)[valid])
      if (length(common_regions) < k + 1) { ari_vals[b] <- NA; next }

      hc_boot <- hclust(dist(boot_z[common_regions, , drop = FALSE], method = "euclidean"), method = "ward.D2")
      clusters_boot <- cutree(hc_boot, k = k)

      ari_vals[b] <- adjustedRandIndex(clusters_obs[common_regions], clusters_boot[names(clusters_boot)])
    }

    stability_results[[key]] <- data.frame(
      condition = cond, sex = sx, n_samples = n_samples, n_regions = nrow(mat_z_obs),
      mean_ari = mean(ari_vals, na.rm = TRUE), median_ari = median(ari_vals, na.rm = TRUE),
      sd_ari = sd(ari_vals, na.rm = TRUE),
      pct_ari_above_0.5 = mean(ari_vals > 0.5, na.rm = TRUE) * 100,
      n_valid_boots = sum(!is.na(ari_vals))
    )

    cat(sprintf("%-25s %-6s: n=%d samples, mean ARI = %.3f (sd %.3f), %.0f%% of resamples ARI>0.5\n",
                cond, sx, n_samples, mean(ari_vals, na.rm = TRUE), sd(ari_vals, na.rm = TRUE),
                mean(ari_vals > 0.5, na.rm = TRUE) * 100))
  }
}

stability_df <- bind_rows(stability_results) %>% arrange(desc(mean_ari))
write.csv(stability_df, "condition_sex_cluster_stability.csv", row.names = FALSE)

cat("\n=== Summary: condition x sex clustering stability (bootstrap, n=1000) ===\n")
print(stability_df)
cat("\nSaved: condition_sex_cluster_stability.csv\n")


# ---------------------------------------------------------------------------
# condition_sex_clustering_agreement
# ---------------------------------------------------------------------------
# ARI and Jaccard agreement between the 8 condition x sex region clusterings
# (condition_sex_region_clusters.csv, from condition_sex_region_clustering_
# zscore.R -- top 40 regions per condition x sex cell, k=5 clusters).
# Same method as condition_clustering_agreement.R / condition_clustering_
# jaccard.R (which did this for the 5 conditions pooled across sex), just
# applied to the 8 condition x sex cells instead.

library(dplyr)
library(mclust)
library(ggplot2)
library(reshape2)

setwd("~/Desktop/Bachelor")
df <- read.csv("condition_sex_region_clusters.csv", stringsAsFactors = FALSE)
df$cell <- paste(df$condition, df$sex, sep = " | ")
cells <- unique(df$cell)
n <- length(cells)

## --- ARI: agreement of cluster labels on shared top-40 regions ---
ari_mat <- matrix(NA, n, n, dimnames = list(cells, cells))
n_common_mat <- matrix(NA, n, n, dimnames = list(cells, cells))

for (i in seq_len(n)) {
  for (j in seq_len(n)) {
    if (i == j) { ari_mat[i, j] <- 1; next }
    d1 <- df %>% filter(cell == cells[i])
    d2 <- df %>% filter(cell == cells[j])
    common <- intersect(d1$region_acronym, d2$region_acronym)
    c1 <- d1$cluster[match(common, d1$region_acronym)]
    c2 <- d2$cluster[match(common, d2$region_acronym)]
    ari_mat[i, j] <- adjustedRandIndex(c1, c2)
    n_common_mat[i, j] <- length(common)
  }
}

cat("=== ARI between the 8 condition x sex clusterings (on shared top-40 regions) ===\n")
print(round(ari_mat, 3))
write.csv(as.data.frame(ari_mat), "condition_sex_clustering_agreement_ari.csv")
cat("\nMean off-diagonal ARI:", round(mean(ari_mat[upper.tri(ari_mat)]), 3), "\n")

## Top 3 most similar and least similar pairs
ari_pairs <- which(upper.tri(ari_mat), arr.ind = TRUE)
ari_pairs_df <- data.frame(cell1 = cells[ari_pairs[, 1]], cell2 = cells[ari_pairs[, 2]],
                            ari = ari_mat[ari_pairs]) %>% arrange(desc(ari))
cat("\nTop 5 most similar cell pairs (ARI):\n"); print(head(ari_pairs_df, 5))
cat("\nTop 5 least similar cell pairs (ARI):\n"); print(tail(ari_pairs_df, 5))

## --- Jaccard: overlap of the top-40 region SETS ---
region_sets <- lapply(cells, function(cl) df$region_acronym[df$cell == cl])
names(region_sets) <- cells

jaccard_mat <- matrix(NA, n, n, dimnames = list(cells, cells))
for (i in seq_len(n)) {
  for (j in seq_len(n)) {
    a <- region_sets[[i]]; b <- region_sets[[j]]
    jaccard_mat[i, j] <- length(intersect(a, b)) / length(union(a, b))
  }
}

cat("\n=== Jaccard between the 8 condition x sex top-40 region sets ===\n")
print(round(jaccard_mat, 3))
write.csv(as.data.frame(jaccard_mat), "condition_sex_clustering_agreement_jaccard.csv")
cat("\nMean off-diagonal Jaccard:", round(mean(jaccard_mat[upper.tri(jaccard_mat)]), 3), "\n")

## --- Heatmap plots (pink scale, matching condition_lowdev_network.R etc.) ---
pink_scale <- c("#FDF0F5", "#FBC9DE", "#F48FB1", "#E05780", "#C2185B", "#8E24AA", "#4A148C")
FONT_PT <- 12
mm_size <- function(pt) pt / (72.27 / 25.4)  # ggplot geom text 'size' aesthetic is in mm

plot_heat <- function(mat, value_name, fname, limits) {
  long <- melt(mat, varnames = c("Cell1", "Cell2"), value.name = value_name)
  p <- ggplot(long, aes(x = Cell2, y = Cell1, fill = .data[[value_name]])) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", .data[[value_name]]), color = .data[[value_name]] > 0.55),
              size = mm_size(FONT_PT), show.legend = FALSE) +
    scale_fill_gradientn(colors = pink_scale, limits = limits) +
    scale_color_manual(values = c("TRUE" = "white", "FALSE" = "#3A3A3A")) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = FONT_PT) +
    theme(panel.background = element_rect(fill = "white", color = NA),
          plot.background  = element_rect(fill = "white", color = NA),
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1, size = FONT_PT),
          axis.text.y = element_text(size = FONT_PT),
          plot.margin = margin(10, 10, 40, 10))
  ggsave(fname, p, width = 14, height = 13, units = "in", dpi = 300)
}

plot_heat(ari_mat, "ARI", "condition_sex_clustering_agreement_ari.png", c(-0.05, 1))
plot_heat(jaccard_mat, "Jaccard", "condition_sex_clustering_agreement_jaccard.png", c(0, 1))

cat("\nSaved:\n")
cat("  condition_sex_clustering_agreement_ari.csv/.png\n")
cat("  condition_sex_clustering_agreement_jaccard.csv/.png\n")


# ---------------------------------------------------------------------------
# condition_sex_region_clustering_zscore
# ---------------------------------------------------------------------------
# Hierarchical clustering (Ward D2, Euclidean, z-scored) of regions WITHIN
# each CONDITION x SEX cell -- crosses condition_matrix_stats.R section 3b
# (per-condition dendrograms) with sex_region_clustering_zscore.R (per-sex
# dendrograms): 4 conditions x 2 sexes = 8 dendrograms. "forced swim ·
# Chronic" is excluded (all-Female, 0 Male -- no Male cell possible), same
# exclusion as plot_heatmap_ALL_DT_by_sex_no_chronicFS.R, giving a clean
# 4x2 grid instead of an unbalanced 5x2 with one empty cell.
#
# Method (identical to the per-condition-only / per-sex-only versions):
#   - top_n_regions (40, or fewer if a cell has < 40 regions with complete
#     data) highest-median-intensity regions for that condition x sex cell
#   - each region's row z-scored across that cell's own animal x channel
#     samples
#   - Ward D2 hierarchical clustering (Euclidean distance) on the z-scored
#     profiles, cut into k_regions clusters
#
# NOTE on sample size: several cells are very small (forced swim acute
# Female n=3, social defeat Female/Male n=3 each, tail suspension n=4 each)
# -- z-scoring and clustering 40 regions across only 3-4 samples is highly
# unstable and should be treated as exploratory/descriptive, not as a
# robust clustering result.

library(dplyr)
library(tidyr)
library(dendextend)

setwd("~/Desktop/Bachelor")
source("prepare_heatmap_data.R")   # long, major_division_per_row, mat, condition_order

long <- long %>% filter(region_acronym != "FS")

top_n_regions <- 40
k_regions <- 5
sex_order <- c("Female", "Male")
condition_order_no_chronicFS <- setdiff(condition_order, "social defeat · Acute")  # social defeat is now single-sex (Male only); excluded here so only conditions comparable across both sexes remain. Chronic FS is included since its sex labels were corrected.

major_cols <- c(
  "Isocortex"             = "#5B7C99", "Olfactory areas"       = "#8C6D46",
  "Hippocampal formation" = "#6B8F71", "Cortical subplate"     = "#9C6B85",
  "Striatum"              = "#4F8F8B", "Pallidum"              = "#7A6B9C",
  "Thalamus"              = "#B08B5A", "Hypothalamus"          = "#A85C5C",
  "Midbrain"              = "#6B6B99", "Pons"                  = "#8FA05E",
  "Medulla"               = "#5E8FA0", "Cerebellum"            = "#9C7A5B",
  "fiber tracts"          = "#A9A9A9", "ventricular systems"   = "#C4B7A6"
)

cluster_results <- list()
hc_by_cell <- list()

for (cond in condition_order_no_chronicFS) {
  for (sx in sex_order) {
    key <- paste(cond, sx, sep = " | ")
    d <- long %>% filter(condition == cond, sex == sx) %>%
      select(region_acronym, animal_id_original, channel, intensity) %>%
      mutate(sample_id = paste(animal_id_original, channel, sep = "_"))
    wide_cell <- d %>% select(region_acronym, sample_id, intensity) %>%
      pivot_wider(names_from = sample_id, values_from = intensity)
    mat_cell <- as.matrix(wide_cell[, -1])
    rownames(mat_cell) <- wide_cell$region_acronym
    mat_cell <- mat_cell[rowSums(is.na(mat_cell)) == 0, , drop = FALSE]

    if (ncol(mat_cell) < 2) {
      cat(sprintf("Skipping %s: only %d sample(s), clustering needs >=2\n", key, ncol(mat_cell)))
      next
    }

    med_per_region <- apply(mat_cell, 1, median, na.rm = TRUE)
    n_top <- min(top_n_regions, length(med_per_region))
    top_regions <- names(sort(med_per_region, decreasing = TRUE))[seq_len(n_top)]
    mat_cell_top <- mat_cell[top_regions, , drop = FALSE]

    mat_cell_z <- t(scale(t(mat_cell_top)))
    mat_cell_z <- mat_cell_z[rowSums(is.na(mat_cell_z)) == 0, , drop = FALSE]

    hc_cell <- hclust(dist(mat_cell_z, method = "euclidean"), method = "ward.D2")
    clusters_cell <- cutree(hc_cell, k = min(k_regions, nrow(mat_cell_z) - 1))

    hc_by_cell[[key]] <- hc_cell
    cluster_results[[key]] <- data.frame(
      condition = cond, sex = sx, region_acronym = names(clusters_cell), cluster = clusters_cell,
      median_intensity = med_per_region[names(clusters_cell)],
      major_division = major_division_per_row[match(names(clusters_cell), rownames(mat))]
    )
    cat(sprintf("%-25s %-6s: %d samples (animal x channel), %d regions clustered, cluster sizes: %s\n",
                cond, sx, ncol(mat_cell_top), nrow(mat_cell_z), paste(table(clusters_cell), collapse = ", ")))
  }
}

cluster_df <- bind_rows(cluster_results) %>% arrange(condition, sex, desc(median_intensity))
write.csv(cluster_df, "condition_sex_region_clusters.csv", row.names = FALSE)

## --- Detailed, readable multi-page PDF: one page per condition x sex cell ---
pdf("condition_sex_region_dendrograms_detailed.pdf", width = 7, height = 9)
for (cond in condition_order_no_chronicFS) {
  for (sx in sex_order) {
    key <- paste(cond, sx, sep = " | ")
    if (is.null(hc_by_cell[[key]])) { cat("Skipping plot for", key, "- no clustering result\n"); next }
    dend <- as.dendrogram(hc_by_cell[[key]])
    par(mar = c(4, 2, 3, 10))
    plot(dend, horiz = TRUE, main = sprintf("%s\n%s (top %d regions)", cond, sx, top_n_regions),
         xlab = "Height", cex.main = 0.9)
  }
}
dev.off()

## --- Compact single-page combined figure: 4 columns (conditions) x 2 rows
## (Female top, Male bottom) -- same column order as the by-sex heatmap ---
pdf("condition_sex_region_dendrograms_compact.pdf", width = 11, height = 7.5)
par(mfrow = c(2, 4), oma = c(2, 0.5, 2, 0.5))
for (sx in sex_order) {
  for (cond in condition_order_no_chronicFS) {
    key <- paste(cond, sx, sep = " | ")
    if (is.null(hc_by_cell[[key]])) { cat("Skipping plot for", key, "- no clustering result\n"); next }
    dend <- as.dendrogram(hc_by_cell[[key]])
    leaf_divs <- major_division_per_row[match(labels(dend), rownames(mat))]
    labels_colors(dend) <- major_cols[as.character(leaf_divs)]
    dend <- set(dend, "labels_cex", 0.38)
    par(mar = c(1, 0.3, 2.3, 3))
    plot(dend, horiz = TRUE,
         main = sprintf("%s\n%s", cond, sx), xlab = "",
         cex.main = 0.62, cex.axis = 0.5)
  }
}
mtext("Height", side = 1, outer = TRUE, cex = 0.7, line = 0.6)
mtext("Condition x Sex region dendrograms (top 40 regions, Ward D2, z-scored)",
      side = 3, outer = TRUE, cex = 0.85, font = 2, line = 0.3)
dev.off()

cat("\n=== Saved ===\n")
cat("  condition_sex_region_clusters.csv\n")
cat("  condition_sex_region_dendrograms_detailed.pdf   (8 pages, readable)\n")
cat("  condition_sex_region_dendrograms_compact.pdf     (1 page, 8 panels)\n")


# ---------------------------------------------------------------------------
# condition_sex_region_dendrograms_A4
# ---------------------------------------------------------------------------
# Per-condition-x-sex region dendrograms, laid out for a DIN A4 page
# (8.27 x 11.69 in), panels labeled (a)-(h) -- same style/method as
# condition_region_dendrograms_A4.R / sex_region_dendrograms_A4.R (top 40
# regions, Ward D2, z-scored), with bootstrap ARI stability annotated per
# panel from condition_sex_cluster_stability.csv. "forced swim · Chronic"
# excluded (all-Female, no Male cell -- see plot_heatmap_ALL_DT_by_sex_no_
# chronicFS.R), giving a clean 4 conditions x 2 sexes = 8-panel grid.
# Layout: rows = condition, columns = Female | Male (reading order a-h
# matches condition_order, Female before Male within each condition).

library(dplyr)
library(tidyr)
library(dendextend)

setwd("~/Desktop/Bachelor")
source("prepare_heatmap_data.R")
long <- long %>% filter(region_acronym != "FS")

top_n_regions <- 40
k_regions <- 5
sex_order <- c("Female", "Male")
condition_order_no_chronicFS <- setdiff(condition_order, "social defeat · Acute")  # social defeat is now single-sex (Male only); excluded here so only conditions comparable across both sexes remain. Chronic FS is included since its sex labels were corrected.

major_cols <- c(
  "Isocortex"             = "#5B7C99", "Olfactory areas"       = "#8C6D46",
  "Hippocampal formation" = "#6B8F71", "Cortical subplate"     = "#9C6B85",
  "Striatum"              = "#4F8F8B", "Pallidum"              = "#7A6B9C",
  "Thalamus"              = "#B08B5A", "Hypothalamus"          = "#A85C5C",
  "Midbrain"              = "#6B6B99", "Pons"                  = "#8FA05E",
  "Medulla"               = "#5E8FA0", "Cerebellum"            = "#9C7A5B",
  "fiber tracts"          = "#A9A9A9", "ventricular systems"   = "#C4B7A6"
)

cells <- expand.grid(sex = sex_order, condition = condition_order_no_chronicFS,
                      stringsAsFactors = FALSE)
## row-major reading order: condition changes slowest, sex fastest -> a-h
cells <- cells[order(match(cells$condition, condition_order_no_chronicFS),
                      match(cells$sex, sex_order)), ]
cells$label <- letters[seq_len(nrow(cells))]

hc_by_cell <- list()
regions_by_cell <- list()
for (i in seq_len(nrow(cells))) {
  cond <- cells$condition[i]; sx <- cells$sex[i]
  key <- paste(cond, sx, sep = " | ")

  d <- long %>% filter(condition == cond, sex == sx) %>%
    select(region_acronym, animal_id_original, channel, intensity) %>%
    mutate(sample_id = paste(animal_id_original, channel, sep = "_"))
  wide_cell <- d %>% select(region_acronym, sample_id, intensity) %>%
    pivot_wider(names_from = sample_id, values_from = intensity)
  mat_cell <- as.matrix(wide_cell[, -1]); rownames(mat_cell) <- wide_cell$region_acronym
  mat_cell <- mat_cell[rowSums(is.na(mat_cell)) == 0, , drop = FALSE]

  if (ncol(mat_cell) < 2) {
    cat(sprintf("Skipping %s: only %d sample(s), clustering needs >=2\n", key, ncol(mat_cell)))
    next
  }

  med_per_region <- apply(mat_cell, 1, median, na.rm = TRUE)
  n_top <- min(top_n_regions, length(med_per_region))
  top_regions <- names(sort(med_per_region, decreasing = TRUE))[seq_len(n_top)]
  mat_cell_top <- mat_cell[top_regions, , drop = FALSE]

  mat_z <- t(scale(t(mat_cell_top)))
  mat_z <- mat_z[rowSums(is.na(mat_z)) == 0, , drop = FALSE]

  hc_cell <- hclust(dist(mat_z, method = "euclidean"), method = "ward.D2")
  hc_cell$height <- hc_cell$height / max(hc_cell$height)  # unify height scale across cells
  hc_by_cell[[key]] <- hc_cell
  regions_by_cell[[key]] <- rownames(mat_z)
}

stability_path <- "condition_sex_cluster_stability.csv"
stability_df <- if (file.exists(stability_path)) read.csv(stability_path) else NULL

## --- One page, 3 rows x 3 cols (8 panels + legend in the last cell);
## drawn once per device so PDF and PNG match ---
draw_page <- function() {
  layout(matrix(c(1:8, 9), nrow = 3, ncol = 3, byrow = TRUE))
  par(oma = c(1, 0.5, 4, 0.5))

  for (i in seq_len(nrow(cells))) {
    cond <- cells$condition[i]; sx <- cells$sex[i]; lab <- cells$label[i]
    key <- paste(cond, sx, sep = " | ")
    if (is.null(hc_by_cell[[key]])) { cat("Skipping plot for", key, "- no clustering result\n"); next }

    dend <- as.dendrogram(hc_by_cell[[key]])
    leaf_divs <- major_division_per_row[match(labels(dend), rownames(mat))]
    labels_colors(dend) <- major_cols[as.character(leaf_divs)]
    dend <- set(dend, "labels_cex", 14 / 12)  # 14pt, matching the other sex-specific figures

    stab_line <- ""
    if (!is.null(stability_df)) {
      s <- stability_df[stability_df$condition == cond & stability_df$sex == sx, ]
      if (nrow(s) == 1) stab_line <- sprintf("\nARI=%.2f+/-%.2f (%.0f%%>0.5)", s$mean_ari, s$sd_ari, s$pct_ari_above_0.5)
    }

    par(mar = c(3, 3, 12, 5), cex = 14 / 12)
    plot(dend, horiz = TRUE,
         main = sprintf("(%s) %s\n\n%s%s", lab, cond, sx, stab_line),
         xlab = "", cex.main = 1.3, cex.axis = 1.3)
  }

  ## legend only lists divisions actually present among the plotted leaves
  ## (drops dead entries like "fiber tracts", which prepare_heatmap_data.R
  ## already excludes from the data entirely -- listing it would wrongly
  ## imply it was a plotted-but-absent category rather than never present)
  all_plotted_regions <- unique(unlist(regions_by_cell))
  divs_present <- as.character(unique(major_division_per_row[match(all_plotted_regions, rownames(mat))]))
  divs_present <- names(major_cols)[names(major_cols) %in% divs_present]  # keep major_cols' own order

  par(mar = c(0.5, 1, 2, 1), cex = 1)
  plot.new()
  legend("center", legend = divs_present, fill = major_cols[divs_present], cex = 1.6,
         bty = "n", title = "Major division", ncol = 2, y.intersp = 1.3)

  mtext("Condition x Sex region dendrograms (top 40 regions, Ward D2, z-scored)",
        side = 3, outer = TRUE, cex = 1.4, font = 2, line = 1.2)
}

## 3x3 grid (8 panels + legend) instead of 2 columns: fewer rows (3 instead
## of 4) means each row can be shorter for the same per-leaf spacing. 40
## leaves/panel at cex=1.3 needs ~40*1.3*12/72 = 8.7in of plot height per
## panel, plus ~2.7in of margin (mar top=12 + bottom=3 lines at this cex) =
## ~11.4in per row; 3 rows of equal height => ~34in total canvas height
## (down from 66in at 4 rows x 2 cols). Width upped slightly since panels
## are now narrower (3 across instead of 2).
pdf("condition_sex_region_dendrograms_A4.pdf", width = 27, height = 42)
draw_page()
dev.off()

png("condition_sex_region_dendrograms_A4.png", width = 27, height = 42, units = "in", res = 150)
draw_page()
dev.off()

cat("Saved: condition_sex_region_dendrograms_A4.pdf, condition_sex_region_dendrograms_A4.png\n")
cat("Panel labels:\n")
print(cells[, c("label", "condition", "sex")])


# ---------------------------------------------------------------------------
# condition_sex_region_dendrograms_detailed_v2
# ---------------------------------------------------------------------------
# Regenerates condition_sex_region_dendrograms_detailed.pdf (originally
# from condition_sex_region_clustering_zscore.R) with division-color leaf
# labels + a legend on EACH page, one full page per condition x sex cell
# (8 pages) instead of the black-and-white version.

library(dplyr)
library(tidyr)
library(dendextend)

setwd("~/Desktop/Bachelor")
source("prepare_heatmap_data.R")
long <- long %>% filter(region_acronym != "FS")

top_n_regions <- 40
sex_order <- c("Female", "Male")
condition_order_no_chronicFS <- setdiff(condition_order, "social defeat · Acute")  # social defeat is now single-sex (Male only); excluded here so only conditions comparable across both sexes remain. Chronic FS is included since its sex labels were corrected.

major_cols <- c(
  "Isocortex"             = "#5B7C99", "Olfactory areas"       = "#8C6D46",
  "Hippocampal formation" = "#6B8F71", "Cortical subplate"     = "#9C6B85",
  "Striatum"              = "#4F8F8B", "Pallidum"              = "#7A6B9C",
  "Thalamus"              = "#B08B5A", "Hypothalamus"          = "#A85C5C",
  "Midbrain"              = "#6B6B99", "Pons"                  = "#8FA05E",
  "Medulla"               = "#5E8FA0", "Cerebellum"            = "#9C7A5B",
  "fiber tracts"          = "#A9A9A9", "ventricular systems"   = "#C4B7A6"
)

hc_by_cell <- list()
for (cond in condition_order_no_chronicFS) {
  for (sx in sex_order) {
    key <- paste(cond, sx, sep = " | ")
    d <- long %>% filter(condition == cond, sex == sx) %>%
      select(region_acronym, animal_id_original, channel, intensity) %>%
      mutate(sample_id = paste(animal_id_original, channel, sep = "_"))
    wide_cell <- d %>% select(region_acronym, sample_id, intensity) %>%
      pivot_wider(names_from = sample_id, values_from = intensity)
    mat_cell <- as.matrix(wide_cell[, -1]); rownames(mat_cell) <- wide_cell$region_acronym
    mat_cell <- mat_cell[rowSums(is.na(mat_cell)) == 0, , drop = FALSE]

    if (ncol(mat_cell) < 2) {
      cat(sprintf("Skipping %s: only %d sample(s), clustering needs >=2\n", key, ncol(mat_cell)))
      next
    }

    med_per_region <- apply(mat_cell, 1, median, na.rm = TRUE)
    n_top <- min(top_n_regions, length(med_per_region))
    top_regions <- names(sort(med_per_region, decreasing = TRUE))[seq_len(n_top)]
    mat_cell_top <- mat_cell[top_regions, , drop = FALSE]

    mat_z <- t(scale(t(mat_cell_top)))
    mat_z <- mat_z[rowSums(is.na(mat_z)) == 0, , drop = FALSE]

    hc_cell <- hclust(dist(mat_z, method = "euclidean"), method = "ward.D2")
    hc_cell$height <- hc_cell$height / max(hc_cell$height)  # unify height scale across cells
    hc_by_cell[[key]] <- hc_cell
  }
}

pdf("condition_sex_region_dendrograms_detailed.pdf", width = 7, height = 9)
for (cond in condition_order_no_chronicFS) {
  for (sx in sex_order) {
    key <- paste(cond, sx, sep = " | ")
    if (is.null(hc_by_cell[[key]])) { cat("Skipping plot for", key, "- no clustering result\n"); next }
    dend <- as.dendrogram(hc_by_cell[[key]])
    leaf_divs <- major_division_per_row[match(labels(dend), rownames(mat))]
    labels_colors(dend) <- major_cols[as.character(leaf_divs)]

    layout(matrix(c(1, 2), nrow = 2), heights = c(0.82, 0.18))
    par(mar = c(4, 2, 4.5, 10))
    plot(dend, horiz = TRUE, main = sprintf("%s\n%s\n\n(top %d regions)", cond, sx, top_n_regions),
         xlab = "Height", cex.main = 0.9)

    par(mar = c(0.5, 1, 1, 1))
    plot.new()
    legend("center", legend = names(major_cols), fill = major_cols, cex = 0.75,
           bty = "n", title = "Major division", ncol = 4)
  }
}
dev.off()

cat("Saved (overwritten): condition_sex_region_dendrograms_detailed.pdf\n")


# ---------------------------------------------------------------------------
# sex_cluster_stability
# ---------------------------------------------------------------------------
# Bootstrap stability of the per-sex region clusterings
# (sex_region_clustering_zscore.R), same method as
# condition_per_condition_cluster_stability.R: resample samples with
# replacement (1000 iterations), recompute the top-40-region clustering
# each time, compare to the observed clustering via Adjusted Rand Index.
#
# Resampling unit here is animal x channel (matching sex_region_clustering_
# zscore.R's sample_id = animal_id_original + channel), not animal alone --
# consistent with how the observed dendrograms were actually built.

library(dplyr)
library(tidyr)
library(mclust)

setwd("~/Desktop/Bachelor")
source("prepare_heatmap_data.R")
long <- long %>% filter(region_acronym != "FS")

top_n_regions <- 40
k_regions <- 5
n_boot <- 1000
set.seed(42)
sex_order <- c("Female", "Male")

stability_results <- vector("list", length(sex_order))

for (sx in sex_order) {
  d <- long %>% filter(sex == sx) %>%
    select(region_acronym, animal_id_original, channel, intensity) %>%
    mutate(sample_id = paste(animal_id_original, channel, sep = "_"))
  wide_sx <- d %>% select(region_acronym, sample_id, intensity) %>%
    pivot_wider(names_from = sample_id, values_from = intensity)
  mat_sx <- as.matrix(wide_sx[, -1]); rownames(mat_sx) <- wide_sx$region_acronym
  mat_sx <- mat_sx[rowSums(is.na(mat_sx)) == 0, , drop = FALSE]

  med_per_region <- apply(mat_sx, 1, median, na.rm = TRUE)
  top_regions <- names(sort(med_per_region, decreasing = TRUE))[seq_len(min(top_n_regions, length(med_per_region)))]
  mat_sx_top <- mat_sx[top_regions, , drop = FALSE]

  mat_z_obs <- t(scale(t(mat_sx_top)))
  mat_z_obs <- mat_z_obs[rowSums(is.na(mat_z_obs)) == 0, , drop = FALSE]
  hc_obs <- hclust(dist(mat_z_obs, method = "euclidean"), method = "ward.D2")
  k <- min(k_regions, nrow(mat_z_obs) - 1)
  clusters_obs <- cutree(hc_obs, k = k)

  n_samples <- ncol(mat_sx_top)
  ari_vals <- numeric(n_boot)

  for (b in seq_len(n_boot)) {
    idx <- sample(seq_len(n_samples), n_samples, replace = TRUE)
    boot_mat <- mat_sx_top[, idx, drop = FALSE]
    boot_z <- t(scale(t(boot_mat)))
    valid <- rowSums(is.na(boot_z)) == 0
    common_regions <- intersect(rownames(mat_z_obs), rownames(boot_z)[valid])
    if (length(common_regions) < k + 1) { ari_vals[b] <- NA; next }

    hc_boot <- hclust(dist(boot_z[common_regions, , drop = FALSE], method = "euclidean"), method = "ward.D2")
    clusters_boot <- cutree(hc_boot, k = k)

    ari_vals[b] <- adjustedRandIndex(clusters_obs[common_regions], clusters_boot[names(clusters_boot)])
  }

  stability_results[[sx]] <- data.frame(
    sex = sx, n_samples = n_samples, n_regions = nrow(mat_z_obs),
    mean_ari = mean(ari_vals, na.rm = TRUE), median_ari = median(ari_vals, na.rm = TRUE),
    sd_ari = sd(ari_vals, na.rm = TRUE),
    pct_ari_above_0.5 = mean(ari_vals > 0.5, na.rm = TRUE) * 100,
    n_valid_boots = sum(!is.na(ari_vals))
  )

  cat(sprintf("%s: n=%d samples, mean ARI = %.3f (sd %.3f), %.0f%% of resamples ARI>0.5\n",
              sx, n_samples, mean(ari_vals, na.rm = TRUE), sd(ari_vals, na.rm = TRUE),
              mean(ari_vals > 0.5, na.rm = TRUE) * 100))
}

stability_df <- bind_rows(stability_results) %>% arrange(desc(mean_ari))
write.csv(stability_df, "sex_cluster_stability.csv", row.names = FALSE)

cat("\n=== Summary: per-sex clustering stability (bootstrap, n=1000) ===\n")
print(stability_df)
cat("\nSaved: sex_cluster_stability.csv\n")


# ---------------------------------------------------------------------------
# sex_region_clustering_zscore
# ---------------------------------------------------------------------------
# Hierarchical clustering (Ward D2, Euclidean, z-scored) of regions WITHIN
# each SEX -- the sex-specific counterpart to condition_matrix_stats.R
# section 3b (condition_region_dendrograms_by_condition.pdf /
# condition_region_clusters_within_condition.csv), which did the same thing
# grouped by condition. Here the grouping variable is Sex (Male, Female)
# instead of Condition, pooling across all 5 conditions (same pooling used
# by sex_kruskal_stats.R for the sex-effect testing).
#
# Method (identical to the per-condition version):
#   - top_n_regions (40) highest-median-intensity regions for that sex
#     (median across that sex's animals)
#   - each region's row z-scored across that sex's own animals
#   - Ward D2 hierarchical clustering (Euclidean distance) on the z-scored
#     profiles, cut into k_regions_per_sex clusters

library(dplyr)
library(tidyr)
library(dendextend)

setwd("~/Desktop/Bachelor")
source("prepare_heatmap_data.R")   # long, major_division_per_row, mat

long <- long %>% filter(region_acronym != "FS")

top_n_regions <- 40
k_regions_per_sex <- 5
sex_order <- c("Female", "Male")

major_cols <- c(
  "Isocortex"             = "#5B7C99", "Olfactory areas"       = "#8C6D46",
  "Hippocampal formation" = "#6B8F71", "Cortical subplate"     = "#9C6B85",
  "Striatum"              = "#4F8F8B", "Pallidum"              = "#7A6B9C",
  "Thalamus"              = "#B08B5A", "Hypothalamus"          = "#A85C5C",
  "Midbrain"              = "#6B6B99", "Pons"                  = "#8FA05E",
  "Medulla"               = "#5E8FA0", "Cerebellum"            = "#9C7A5B",
  "fiber tracts"          = "#A9A9A9", "ventricular systems"   = "#C4B7A6"
)

per_sex_clusters <- vector("list", length(sex_order))
hc_by_sex <- vector("list", length(sex_order))
names(per_sex_clusters) <- sex_order
names(hc_by_sex) <- sex_order

for (sx in sex_order) {
  d <- long %>% filter(sex == sx) %>%
    select(region_acronym, animal_id_original, channel, intensity) %>%
    mutate(sample_id = paste(animal_id_original, channel, sep = "_"))
  wide_sx <- d %>% select(region_acronym, sample_id, intensity) %>%
    pivot_wider(names_from = sample_id, values_from = intensity)
  mat_sx <- as.matrix(wide_sx[, -1])
  rownames(mat_sx) <- wide_sx$region_acronym
  mat_sx <- mat_sx[rowSums(is.na(mat_sx)) == 0, , drop = FALSE]

  med_per_region <- apply(mat_sx, 1, median, na.rm = TRUE)
  top_regions <- names(sort(med_per_region, decreasing = TRUE))[seq_len(min(top_n_regions, length(med_per_region)))]
  mat_sx_top <- mat_sx[top_regions, , drop = FALSE]

  mat_sx_z <- t(scale(t(mat_sx_top)))
  mat_sx_z <- mat_sx_z[rowSums(is.na(mat_sx_z)) == 0, , drop = FALSE]

  hc_sx <- hclust(dist(mat_sx_z, method = "euclidean"), method = "ward.D2")
  clusters_sx <- cutree(hc_sx, k = min(k_regions_per_sex, nrow(mat_sx_z) - 1))

  hc_by_sex[[sx]] <- hc_sx
  per_sex_clusters[[sx]] <- data.frame(
    sex = sx, region_acronym = names(clusters_sx), cluster = clusters_sx,
    median_intensity = med_per_region[names(clusters_sx)],
    major_division = major_division_per_row[match(names(clusters_sx), rownames(mat))]
  )
  cat(sprintf("\n%s: %d samples (animal x channel), top %d regions clustered, cluster sizes: %s\n",
              sx, ncol(mat_sx_top), nrow(mat_sx_z), paste(table(clusters_sx), collapse = ", ")))
}

per_sex_clusters_df <- bind_rows(per_sex_clusters) %>% arrange(sex, desc(median_intensity))
write.csv(per_sex_clusters_df, "sex_region_clusters_within_sex.csv", row.names = FALSE)

## --- Combined figure: Female + Male dendrograms side by side, leaves
## colored by major division (same convention as Figure 9/10) ---
png("sex_region_dendrograms_combined.png", width = 6, height = 10, units = "in", res = 300)
par(mfrow = c(2, 1))
for (sx in sex_order) {
  dend <- as.dendrogram(hc_by_sex[[sx]])
  leaf_divs <- major_division_per_row[match(labels(dend), rownames(mat))]
  labels_colors(dend) <- major_cols[as.character(leaf_divs)]

  par(mar = c(3, 0.5, 3, 8), cex = 0.65)
  plot(dend, horiz = TRUE,
       main = sprintf("%s (top %d regions)", sx, top_n_regions),
       xlab = "Height", cex.main = 1.1, cex.axis = 0.8, cex.lab = 0.9)
}
dev.off()
par(mfrow = c(1, 1), cex = 1)

pdf("sex_region_dendrograms_by_sex.pdf", width = 7, height = 9)
for (sx in sex_order) {
  dend <- as.dendrogram(hc_by_sex[[sx]])
  par(mar = c(4, 2, 3, 10))
  plot(dend, horiz = TRUE, main = sprintf("%s\n(top %d regions by median intensity)", sx, top_n_regions),
       xlab = "Height", cex.main = 0.9)
}
dev.off()

cat("\n=== Saved ===\n")
cat("  sex_region_clusters_within_sex.csv\n")
cat("  sex_region_dendrograms_combined.png\n")
cat("  sex_region_dendrograms_by_sex.pdf\n")


# ---------------------------------------------------------------------------
# sex_region_dendrograms_A4
# ---------------------------------------------------------------------------
# Per-sex region dendrograms, laid out for a DIN A4 page (8.27 x 11.69 in),
# panels labeled (a) Female, (b) Male -- same style/method as
# condition_region_dendrograms_A4.R (top 40 regions, Ward D2, z-scored),
# with bootstrap ARI stability annotated per panel from
# sex_cluster_stability.csv.

library(dplyr)
library(tidyr)
library(dendextend)

setwd("~/Desktop/Bachelor")
source("prepare_heatmap_data.R")
long <- long %>% filter(region_acronym != "FS")

top_n_regions <- 40
k_regions <- 5
sex_order <- c("Female", "Male")

major_cols <- c(
  "Isocortex"             = "#5B7C99", "Olfactory areas"       = "#8C6D46",
  "Hippocampal formation" = "#6B8F71", "Cortical subplate"     = "#9C6B85",
  "Striatum"              = "#4F8F8B", "Pallidum"              = "#7A6B9C",
  "Thalamus"              = "#B08B5A", "Hypothalamus"          = "#A85C5C",
  "Midbrain"              = "#6B6B99", "Pons"                  = "#8FA05E",
  "Medulla"               = "#5E8FA0", "Cerebellum"            = "#9C7A5B",
  "fiber tracts"          = "#A9A9A9", "ventricular systems"   = "#C4B7A6"
)

hc_by_sex <- vector("list", length(sex_order))
names(hc_by_sex) <- sex_order

for (sx in sex_order) {
  d <- long %>% filter(sex == sx) %>%
    select(region_acronym, animal_id_original, channel, intensity) %>%
    mutate(sample_id = paste(animal_id_original, channel, sep = "_"))
  wide_sx <- d %>% select(region_acronym, sample_id, intensity) %>%
    pivot_wider(names_from = sample_id, values_from = intensity)
  mat_sx <- as.matrix(wide_sx[, -1]); rownames(mat_sx) <- wide_sx$region_acronym
  mat_sx <- mat_sx[rowSums(is.na(mat_sx)) == 0, , drop = FALSE]

  med_per_region <- apply(mat_sx, 1, median, na.rm = TRUE)
  top_regions <- names(sort(med_per_region, decreasing = TRUE))[seq_len(min(top_n_regions, length(med_per_region)))]
  mat_sx_top <- mat_sx[top_regions, , drop = FALSE]

  mat_sx_z <- t(scale(t(mat_sx_top)))
  mat_sx_z <- mat_sx_z[rowSums(is.na(mat_sx_z)) == 0, , drop = FALSE]

  hc_sx <- hclust(dist(mat_sx_z, method = "euclidean"), method = "ward.D2")
  hc_sx$height <- hc_sx$height / max(hc_sx$height)  # unify height scale across sexes
  hc_by_sex[[sx]] <- hc_sx
}

stability_path <- "sex_cluster_stability.csv"
stability_df <- if (file.exists(stability_path)) read.csv(stability_path) else NULL

panel_letters <- setNames(letters[seq_along(sex_order)], sex_order)

## --- DIN A4 portrait: 8.27 x 11.69 in, 3 rows (Female, Male, legend) --
## more vertical room per panel than the 5-condition version since only
## 2 dendrograms are needed.
pdf("sex_region_dendrograms_A4.pdf", width = 8.27, height = 11.69)
layout(matrix(c(1, 2, 3), nrow = 3), heights = c(0.42, 0.42, 0.16))
par(oma = c(1, 0.5, 3, 0.5))

for (sx in sex_order) {
  dend <- as.dendrogram(hc_by_sex[[sx]])
  leaf_divs <- major_division_per_row[match(labels(dend), rownames(mat))]
  labels_colors(dend) <- major_cols[as.character(leaf_divs)]

  stab_line <- ""
  if (!is.null(stability_df) && sx %in% stability_df$sex) {
    s <- stability_df[stability_df$sex == sx, ]
    stab_line <- sprintf("\nARI = %.2f +/- %.2f | ARI>0.5: %.0f%%", s$mean_ari, s$sd_ari, s$pct_ari_above_0.5)
  }

  par(mar = c(3, 0.5, 4.2, 8), cex = 0.6)
  plot(dend, horiz = TRUE,
       main = sprintf("(%s) %s\n\n(top %d regions)%s", panel_letters[sx], sx, top_n_regions, stab_line),
       xlab = "Height", cex.main = 1, cex.axis = 0.85, cex.lab = 0.95)
}

par(mar = c(1, 1, 1, 1), cex = 1)
plot.new()
legend("center", legend = names(major_cols), fill = major_cols, cex = 0.9,
       bty = "n", title = "Major division", ncol = 3)

mtext("Per-sex region dendrograms (top 40 regions, Ward D2, z-scored)",
      side = 3, outer = TRUE, cex = 0.95, font = 2, line = 1.2)
dev.off()

cat("Saved: sex_region_dendrograms_A4.pdf\n")
cat("Panel labels: ", paste(sprintf("(%s)=%s", panel_letters, names(panel_letters)), collapse = "; "), "\n")


# ---------------------------------------------------------------------------
# sex_region_dendrograms_by_sex_v2
# ---------------------------------------------------------------------------
# Regenerates sex_region_dendrograms_by_sex.pdf (originally from
# sex_region_clustering_zscore.R) with division-color leaf labels + a
# legend on EACH page, matching sex_region_dendrograms_A4.R's style, but
# one full page per sex instead of a combined grid.

library(dplyr)
library(tidyr)
library(dendextend)

setwd("~/Desktop/Bachelor")
source("prepare_heatmap_data.R")
long <- long %>% filter(region_acronym != "FS")

top_n_regions <- 40
sex_order <- c("Female", "Male")

major_cols <- c(
  "Isocortex"             = "#5B7C99", "Olfactory areas"       = "#8C6D46",
  "Hippocampal formation" = "#6B8F71", "Cortical subplate"     = "#9C6B85",
  "Striatum"              = "#4F8F8B", "Pallidum"              = "#7A6B9C",
  "Thalamus"              = "#B08B5A", "Hypothalamus"          = "#A85C5C",
  "Midbrain"              = "#6B6B99", "Pons"                  = "#8FA05E",
  "Medulla"               = "#5E8FA0", "Cerebellum"            = "#9C7A5B",
  "fiber tracts"          = "#A9A9A9", "ventricular systems"   = "#C4B7A6"
)

hc_by_sex <- vector("list", length(sex_order))
names(hc_by_sex) <- sex_order

for (sx in sex_order) {
  d <- long %>% filter(sex == sx) %>%
    select(region_acronym, animal_id_original, channel, intensity) %>%
    mutate(sample_id = paste(animal_id_original, channel, sep = "_"))
  wide_sx <- d %>% select(region_acronym, sample_id, intensity) %>%
    pivot_wider(names_from = sample_id, values_from = intensity)
  mat_sx <- as.matrix(wide_sx[, -1]); rownames(mat_sx) <- wide_sx$region_acronym
  mat_sx <- mat_sx[rowSums(is.na(mat_sx)) == 0, , drop = FALSE]

  med_per_region <- apply(mat_sx, 1, median, na.rm = TRUE)
  top_regions <- names(sort(med_per_region, decreasing = TRUE))[seq_len(min(top_n_regions, length(med_per_region)))]
  mat_sx_top <- mat_sx[top_regions, , drop = FALSE]

  mat_sx_z <- t(scale(t(mat_sx_top)))
  mat_sx_z <- mat_sx_z[rowSums(is.na(mat_sx_z)) == 0, , drop = FALSE]

  hc_sx <- hclust(dist(mat_sx_z, method = "euclidean"), method = "ward.D2")
  hc_sx$height <- hc_sx$height / max(hc_sx$height)  # unify height scale across sexes
  hc_by_sex[[sx]] <- hc_sx
}

pdf("sex_region_dendrograms_by_sex.pdf", width = 7, height = 9)
for (sx in sex_order) {
  dend <- as.dendrogram(hc_by_sex[[sx]])
  leaf_divs <- major_division_per_row[match(labels(dend), rownames(mat))]
  labels_colors(dend) <- major_cols[as.character(leaf_divs)]

  layout(matrix(c(1, 2), nrow = 2), heights = c(0.82, 0.18))
  par(mar = c(4, 2, 4.5, 10))
  plot(dend, horiz = TRUE, main = sprintf("%s\n\n(top %d regions by median intensity)", sx, top_n_regions),
       xlab = "Height", cex.main = 0.9)

  par(mar = c(0.5, 1, 1, 1))
  plot.new()
  legend("center", legend = names(major_cols), fill = major_cols, cex = 0.75,
         bty = "n", title = "Major division", ncol = 4)
}
dev.off()

cat("Saved (overwritten): sex_region_dendrograms_by_sex.pdf\n")

