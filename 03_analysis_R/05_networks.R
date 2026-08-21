# =============================================================================
# THESIS FIGURES PRODUCED BY THIS SCRIPT
#   Figure 16B     Co-variation network of the consistent projection core
#   Figure 24B     Network structures of the sex-specific top regions
#
# Note: figure numbers appearing in the comments further down refer to earlier
# drafts and are NOT the final numbering. The list above is authoritative.
# =============================================================================

# =============================================================================
# Co-variation networks and hub metrics
#
# Consolidated from the following original scripts:
#   - condition_lowdev_network.R
#   - condition_cluster_network_hubs.R
#   - condition_per_condition_networks.R
#   - nullmodel_corrected_coactivation_networks.R
#   - nullmodel_corrected_coactivation_networks_normalized.R
#   - sex_coactivation_networks.R
#   - sex_lowdev_network.R
#   - coactive_projection_regions.R
#   - fs_acute_vs_chronic_female_networks.R
# =============================================================================


# ---------------------------------------------------------------------------
# condition_lowdev_network
# ---------------------------------------------------------------------------
# Coactivation network restricted to the top 20 "low deviation, high
# intensity" regions (condition_top_regions_low_deviation.csv / Figure 13) --
# same method as condition_cluster_network_hubs.R (Spearman |rho|>0.6 across
# all animal x channel samples, whole-brain-median normalized), but applied
# only to these 20 consistently-engaged regions instead of a z-score cluster.

library(dplyr)
library(tidyr)
library(igraph)
library(ggraph)
library(ggplot2)

setwd("~/Desktop/Bachelor")
source("prepare_heatmap_data.R")
long <- long %>% filter(region_acronym != "FS")

top20 <- read.csv("condition_top_regions_low_deviation.csv") %>% head(20)

major_cols <- c(
  "Isocortex"             = "#5B7C99", "Olfactory areas"       = "#8C6D46",
  "Hippocampal formation" = "#6B8F71", "Cortical subplate"     = "#9C6B85",
  "Striatum"              = "#4F8F8B", "Pallidum"              = "#7A6B9C",
  "Thalamus"              = "#B08B5A", "Hypothalamus"          = "#A85C5C",
  "Midbrain"              = "#6B6B99", "Pons"                  = "#8FA05E",
  "Medulla"               = "#5E8FA0", "Cerebellum"            = "#9C7A5B",
  "fiber tracts"          = "#A9A9A9", "ventricular systems"   = "#C4B7A6"
)

sample_overall <- long %>%
  group_by(animal_id_original, channel) %>%
  summarise(overall_median = median(intensity, na.rm = TRUE), .groups = "drop")

long_norm <- long %>%
  left_join(sample_overall, by = c("animal_id_original", "channel")) %>%
  mutate(norm_value = intensity / overall_median,
         sample_id = paste(animal_id_original, channel, sep = "_"))

wide <- long_norm %>%
  filter(region_acronym %in% top20$region_acronym) %>%
  select(sample_id, region_acronym, norm_value) %>%
  pivot_wider(names_from = region_acronym, values_from = norm_value)

mat_samples <- as.matrix(wide[, -1])
rownames(mat_samples) <- wide$sample_id
mat_samples <- mat_samples[, colSums(!is.na(mat_samples)) >= 10, drop = FALSE]

cor_mat <- cor(mat_samples, method = "spearman", use = "pairwise.complete.obs")
cor_thresh <- 0.6
adj <- (abs(cor_mat) > cor_thresh) & (abs(cor_mat) < 1)
diag(adj) <- FALSE

g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
deg <- as.numeric(igraph::degree(g)); names(deg) <- igraph::V(g)$name
eig <- tryCatch(igraph::eigen_centrality(g)$vector, error = function(e) rep(NA, igraph::vcount(g)))
divs <- top20$major_division[match(names(deg), top20$region_acronym)]

hub_df <- data.frame(region_acronym = names(deg), major_division = divs,
                      degree = deg, eigen_centrality = round(eig, 3)) %>%
  arrange(desc(degree), desc(eigen_centrality))
write.csv(hub_df, "condition_lowdev_network_hubs.csv", row.names = FALSE)

cat(sprintf("=== Low-deviation-region network (n=%d regions, %d edges at |rho|>%.1f) ===\n",
            length(deg), sum(adj) / 2, cor_thresh))
print(hub_df)

## Drop fully isolated regions (degree = 0) AND weakly-tethered ones
## (degree = 1, e.g. a single edge stretching far from the core) so the
## force-directed layout is not dominated by empty space; they are listed
## separately instead of being scattered across (or stretching) the canvas.
isolated <- names(deg)[deg <= 1]
g_connected <- igraph::delete_vertices(g, isolated)

igraph::V(g_connected)$degree <- deg[igraph::V(g_connected)$name]
igraph::V(g_connected)$division <- as.character(divs[match(igraph::V(g_connected)$name, names(deg))])

cat(sprintf("\nIsolated regions (degree = 0, excluded from plot): %s\n", paste(isolated, collapse = ", ")))

pink_scale <- c("#FDF0F5", "#FBC9DE", "#F48FB1", "#E05780", "#C2185B", "#8E24AA", "#4A148C")
FONT_PT <- 14
mm_size <- function(pt) pt / (72.27 / 25.4)  # ggplot geom text 'size' aesthetic is in mm

p <- ggraph(g_connected, layout = "fr") +
  geom_edge_link(alpha = 0.35, color = "grey55", width = 0.6) +
  geom_node_point(aes(size = degree, color = degree)) +
  geom_node_text(aes(label = name), nudge_y = 0.22, vjust = 0, size = mm_size(FONT_PT),
                 color = "#2A2A2A", fontface = "bold") +
  scale_color_gradientn(colours = pink_scale, name = "Degree") +
  scale_size_continuous(range = c(6, 16), guide = "none") +
  scale_x_continuous(expand = expansion(mult = 0.25)) +
  scale_y_continuous(expand = expansion(mult = 0.25)) +
  coord_fixed() +
  theme_void(base_size = FONT_PT) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.position = "bottom",
    legend.text = element_text(size = FONT_PT),
    legend.title = element_text(size = FONT_PT),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  guides(color = guide_colorbar(barwidth = unit(10, "cm"), barheight = unit(0.8, "cm")))

ggsave("condition_lowdev_network.png", p, width = 12, height = 12, units = "in", dpi = 300, limitsize = FALSE)
cat("\nSaved: condition_lowdev_network.png, condition_lowdev_network_hubs.csv\n")


# ---------------------------------------------------------------------------
# condition_cluster_network_hubs
# ---------------------------------------------------------------------------
# Network hub analysis WITHIN each of the 6 region clusters identified by
# hierarchical clustering of z-scored condition profiles
# (condition_region_clusters.csv, condition_matrix_stats.R section 3).
# Rationale: regions in the same cluster share a similar condition-response
# profile: it is a reasonable next question whether they also co-vary
# together across ANIMALS (a coactivation network, same logic as
# nullmodel_corrected_coactivation_networks.R), and if so, which regions act as
# network hubs -- plausible circuit-level integrators within that module,
# rather than isolated nodes.
#
# Network: nodes = regions within a cluster; edges = pairs with |Spearman
# rho| > 0.6 across all animal x channel samples (whole-brain-median
# normalized, to remove global scaling); hub metrics = degree and
# eigenvector centrality.

library(dplyr)
library(tidyr)
library(igraph)
library(ggraph)
library(ggplot2)

setwd("~/Desktop/Bachelor")
source("prepare_heatmap_data.R")   # long, mat

keep <- rownames(mat) != "FS"
long <- long %>% filter(region_acronym != "FS")

clusters <- read.csv("condition_region_clusters.csv")   # region_acronym, cluster, major_division

## --- Whole-brain-median normalization per animal x channel (removes global
## scaling before correlating regions with each other) ---
sample_overall <- long %>%
  group_by(animal_id_original, channel) %>%
  summarise(overall_median = median(intensity, na.rm = TRUE), .groups = "drop")

long_norm <- long %>%
  left_join(sample_overall, by = c("animal_id_original", "channel")) %>%
  mutate(norm_value = intensity / overall_median,
         sample_id = paste(animal_id_original, channel, sep = "_"))

wide <- long_norm %>%
  select(sample_id, region_acronym, norm_value) %>%
  pivot_wider(names_from = region_acronym, values_from = norm_value)

mat_samples <- as.matrix(wide[, -1])
rownames(mat_samples) <- wide$sample_id

cor_thresh <- 0.6
hub_results <- list()

for (cl in sort(unique(clusters$cluster))) {
  regs_cl <- clusters$region_acronym[clusters$cluster == cl]
  regs_cl <- regs_cl[regs_cl %in% colnames(mat_samples)]
  sub <- mat_samples[, regs_cl, drop = FALSE]
  sub <- sub[, colSums(!is.na(sub)) >= 10, drop = FALSE]   # need >=10 valid samples for a meaningful correlation
  if (ncol(sub) < 3) { cat(sprintf("\nCluster %d: skipped, only %d regions with >=10 valid samples\n", cl, ncol(sub))); next }

  cor_mat <- cor(sub, method = "spearman", use = "pairwise.complete.obs")
  adj <- (abs(cor_mat) > cor_thresh) & (abs(cor_mat) < 1)
  diag(adj) <- FALSE

  g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
  deg <- as.numeric(igraph::degree(g))
  names(deg) <- igraph::V(g)$name
  eig <- tryCatch(igraph::eigen_centrality(g)$vector, error = function(e) rep(NA, igraph::vcount(g)))

  hub_df <- data.frame(
    cluster = cl, region_acronym = names(deg), degree = deg, eigen_centrality = round(eig, 3),
    n_regions_in_cluster = ncol(sub)
  ) %>% arrange(desc(degree), desc(eigen_centrality))

  hub_results[[as.character(cl)]] <- hub_df

  cat(sprintf("\n=== Cluster %d (%d regions, %d edges at |rho|>%.1f) -- top 5 hubs ===\n",
              cl, ncol(sub), sum(adj) / 2, cor_thresh))
  print(head(hub_df, 5))

  ## Plot small clusters (<=40 regions) as an actual graph
  if (ncol(sub) <= 40 && sum(adj) > 0) {
    igraph::V(g)$degree <- deg
    p <- ggraph(g, layout = "fr") +
      geom_edge_link(alpha = 0.3) +
      geom_node_point(aes(size = degree), color = "#C2185B") +
      geom_node_text(aes(label = name), repel = TRUE, size = 3) +
      labs(title = sprintf("Cluster %d coactivation network (n=%d regions, |rho|>%.1f)", cl, ncol(sub), cor_thresh)) +
      theme_void(base_size = 10)
    ggsave(sprintf("condition_cluster_network_cluster%d.png", cl), p, width = 7, height = 6, units = "in", dpi = 300)
  }
}

hub_all <- bind_rows(hub_results)
write.csv(hub_all, "condition_cluster_network_hubs.csv", row.names = FALSE)

cat("\n\n=== Overall top 15 hub regions across all clusters (by degree) ===\n")
print(head(hub_all %>% arrange(desc(degree)), 15))

cat("\nSaved: condition_cluster_network_hubs.csv, condition_cluster_network_cluster<k>.png (small clusters)\n")


# ---------------------------------------------------------------------------
# condition_per_condition_networks
# ---------------------------------------------------------------------------
# Per-condition coactivation networks: same method as
# condition_cluster_network_hubs.R (Spearman |rho|>0.6 across animals,
# degree/eigenvector centrality), but applied separately to each
# condition's own top-40 regions and own animals -- the network-graph
# counterpart to the per-condition dendrograms (Figure 9/10), which can
# only show each region in a single branch. A network can show a region
# connected to multiple sub-groups at once.

library(dplyr)
library(tidyr)
library(igraph)
library(ggraph)
library(ggplot2)
library(patchwork)

setwd("~/Desktop/Bachelor")
source("prepare_heatmap_data.R")
long <- long %>% filter(region_acronym != "FS")

top_n_regions <- 40
cor_thresh <- 0.6

major_cols <- c(
  "Isocortex"             = "#5B7C99", "Olfactory areas"       = "#8C6D46",
  "Hippocampal formation" = "#6B8F71", "Cortical subplate"     = "#9C6B85",
  "Striatum"              = "#4F8F8B", "Pallidum"              = "#7A6B9C",
  "Thalamus"              = "#B08B5A", "Hypothalamus"          = "#A85C5C",
  "Midbrain"              = "#6B6B99", "Pons"                  = "#8FA05E",
  "Medulla"               = "#5E8FA0", "Cerebellum"            = "#9C7A5B",
  "fiber tracts"          = "#A9A9A9", "ventricular systems"   = "#C4B7A6"
)

hub_results <- list()
plots <- list()

for (cond in condition_order) {
  d <- long %>% filter(condition == cond) %>% select(region_acronym, animal_id_original, intensity)
  wide_cond <- d %>% pivot_wider(names_from = animal_id_original, values_from = intensity)
  mat_cond <- as.matrix(wide_cond[, -1])
  rownames(mat_cond) <- wide_cond$region_acronym
  mat_cond <- mat_cond[rowSums(is.na(mat_cond)) == 0, , drop = FALSE]

  med_per_region <- apply(mat_cond, 1, median, na.rm = TRUE)
  top_regions <- names(sort(med_per_region, decreasing = TRUE))[seq_len(min(top_n_regions, length(med_per_region)))]
  mat_cond_top <- mat_cond[top_regions, , drop = FALSE]

  cor_mat <- cor(t(mat_cond_top), method = "spearman", use = "pairwise.complete.obs")
  adj <- (abs(cor_mat) > cor_thresh) & (abs(cor_mat) < 1)
  diag(adj) <- FALSE

  g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
  deg <- as.numeric(igraph::degree(g)); names(deg) <- igraph::V(g)$name
  eig <- tryCatch(igraph::eigen_centrality(g)$vector, error = function(e) rep(NA, igraph::vcount(g)))

  divs <- major_division_per_row[match(names(deg), rownames(mat))]

  hub_df <- data.frame(
    condition = cond, region_acronym = names(deg), major_division = divs,
    degree = deg, eigen_centrality = round(eig, 3), n_animals = ncol(mat_cond_top)
  ) %>% arrange(desc(degree), desc(eigen_centrality))
  hub_results[[cond]] <- hub_df

  cat(sprintf("\n=== %s (%d animals, %d regions, %d edges at |rho|>%.1f) -- top 5 hubs ===\n",
              cond, ncol(mat_cond_top), length(deg), sum(adj) / 2, cor_thresh))
  print(head(hub_df, 5))

  igraph::V(g)$degree <- deg
  igraph::V(g)$division <- as.character(divs)
  p <- ggraph(g, layout = "fr") +
    geom_edge_link(alpha = 0.25, color = "grey50") +
    geom_node_point(aes(size = degree, color = division)) +
    geom_node_text(aes(label = name), repel = TRUE, size = 2.3, max.overlaps = 30) +
    scale_color_manual(values = major_cols, name = "Major division") +
    scale_size_continuous(range = c(1, 5), guide = "none") +
    labs(title = cond) +
    theme_void(base_size = 9) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10), legend.position = "none")
  plots[[cond]] <- p
}

hub_all <- bind_rows(hub_results)
write.csv(hub_all, "condition_per_condition_network_hubs.csv", row.names = FALSE)

cat("\n\n=== Top hub per condition (by degree) ===\n")
print(hub_all %>% group_by(condition) %>% slice_max(degree, n = 1, with_ties = TRUE) %>%
        select(condition, region_acronym, major_division, degree, eigen_centrality))

## Combined figure: 5 networks + 1 shared legend, same 2x3 layout as Figure 9/10
legend_plot <- ggplot(data.frame(x = 1, y = seq_along(major_cols), division = names(major_cols)),
                       aes(x, y, color = division)) +
  geom_point(size = 3) +
  scale_color_manual(values = major_cols, name = "Major division") +
  theme_void() +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(legend.position = "right", legend.title = element_text(face = "bold"))
legend_only <- cowplot::get_legend(legend_plot)

combined <- wrap_plots(c(plots, list(cowplot::ggdraw(legend_only))), ncol = 3) +
  plot_annotation(title = "Per-condition coactivation networks (top 40 regions, |rho| > 0.6)",
                   theme = theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5)))

ggsave("condition_per_condition_networks_combined.png", combined, width = 13, height = 9, units = "in", dpi = 300)

cat("\nSaved: condition_per_condition_network_hubs.csv, condition_per_condition_networks_combined.png\n")


# ---------------------------------------------------------------------------
# nullmodel_corrected_coactivation_networks
# ---------------------------------------------------------------------------
# Co-activation network analysis per stressor condition, using a null-model-corrected approach
# (2026):
#   1. Pairwise Spearman correlations between regions, across the animals
#      belonging to one condition.
#   2. Edges retained only if p < 0.01 (unadjusted) AND >= 3 paired obs.
#   3. Null-model correction: 1,000 degree-preserving network permutations
#      (igraph edge rewiring, preserving each node's degree); the mean edge
#      weight across these permutations is subtracted from every observed
#      edge weight. Edges with negligible corrected weight are dropped.
#   4. Louvain community detection on the corrected network.
#   5. Visualization: top 100 nodes by eigenvector centrality.
#
# One network per condition (5 total). Regions restricted to those with
# n >= 3 valid observations within that condition (same rule applied
# for region inclusion). 'FS' excluded (injection region, trivially 1.0).

library(dplyr)
library(tidyr)
library(igraph)
library(ggraph)
library(tidygraph)
library(ggplot2)

df <- read.csv("~/Desktop/ALL_DT.csv", stringsAsFactors = FALSE)
cond_map <- read.csv("~/Desktop/group_condition_corrected.csv", stringsAsFactors = FALSE)

## --- Node set: restricted to a meaningful region subset, analogous to
## the "FDR-significant regions from Level 1" -- since no region survived
## FDR here (see nullmodel_corrected_pairwise_stats.csv), we substitute the top
## candidates across all 10 pairwise comparisons (lowest p-value per region,
## union across comparisons) plus the previously curated top/most-variable
## regions. This keeps the network at a computationally tractable size
## (all-region networks produced 100k+ edges -- infeasible for 1000
## permutations) while still reflecting "regions of interest" rather than
## an arbitrary/random subset.
pairwise_stats <- read.csv("~/Desktop/nullmodel_corrected_pairwise_stats.csv", stringsAsFactors = FALSE)
top_by_pvalue <- pairwise_stats %>% group_by(region_acronym) %>%
  summarise(min_p = min(p_value), .groups = "drop") %>%
  arrange(min_p) %>% slice_head(n = 80) %>% pull(region_acronym)

top_regions_df <- read.csv("~/Desktop/top_regions_low_deviation.csv", stringsAsFactors = FALSE)
most_var_df <- read.csv("~/Desktop/most_variable_regions.csv", stringsAsFactors = FALSE)
node_regions <- unique(c(top_by_pvalue, top_regions_df$region_acronym, most_var_df$region_acronym))
cat("n node regions for network analysis:", length(node_regions), "\n")

long <- df %>%
  filter(region_acronym != "FS", region_acronym %in% node_regions) %>%
  select(animal_group, animal_id_original, region_acronym,
         trapcre_intensity, traptta_intensity) %>%
  pivot_longer(cols = c(trapcre_intensity, traptta_intensity),
               names_to = "channel", values_to = "value") %>%
  mutate(channel = ifelse(channel == "trapcre_intensity", "Cre", "tTA"),
         sample_id = paste(animal_id_original, channel, sep = "_")) %>%
  rename(Group = animal_group, Channel = channel) %>%
  left_join(cond_map, by = c("Group", "Channel")) %>%
  filter(!is.na(Condition))

conditions <- sort(unique(long$Condition))

## --- Fast Spearman correlation + p-value matrix (t-approximation) ---
spearman_with_p <- function(mat) {
  r <- cor(mat, method = "spearman", use = "pairwise.complete.obs")
  not_na <- !is.na(mat)
  n_pairs <- crossprod(not_na)   # n valid paired observations per region pair
  t_stat <- r * sqrt((n_pairs - 2) / (1 - r^2))
  p <- 2 * pt(-abs(t_stat), df = n_pairs - 2)
  diag(p) <- 1
  list(r = r, p = p, n = n_pairs)
}

build_network <- function(cond, n_perm = 1000, min_obs = 3) {
  d <- long %>% filter(Condition == cond)
  mat <- d %>% select(sample_id, region_acronym, value) %>%
    pivot_wider(names_from = region_acronym, values_from = value) %>%
    tibble::column_to_rownames("sample_id") %>% as.matrix()

  ## keep only regions with >= min_obs non-NA values
  keep <- colSums(!is.na(mat)) >= min_obs
  mat <- mat[, keep, drop = FALSE]
  cat(sprintf("[%s] n samples = %d, n regions (n>=%d) = %d\n", cond, nrow(mat), min_obs, ncol(mat)))

  sp <- spearman_with_p(mat)

  ## --- Threshold edges: p < 0.01 unadjusted, n_pairs >= 3 ---
  edge_mask <- (sp$p < 0.01) & (sp$n >= min_obs) & upper.tri(sp$r)
  edges <- which(edge_mask, arr.ind = TRUE)
  if (nrow(edges) == 0) {
    cat(sprintf("[%s] No edges survived thresholding.\n", cond))
    return(NULL)
  }
  edge_df <- data.frame(
    from = colnames(mat)[edges[, 1]], to = colnames(mat)[edges[, 2]],
    weight = sp$r[edges]
  )
  cat(sprintf("[%s] n edges after thresholding: %d\n", cond, nrow(edge_df)))

  g <- graph_from_data_frame(edge_df, directed = FALSE)

  ## --- Null model: 1000 degree-preserving rewirings; subtract mean null
  ## edge weight (using the same weight multiset, reshuffled onto rewired
  ## topology) from every observed edge weight ---
  null_means <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    g_null <- tryCatch(rewire(g, with = keeping_degseq(niter = min(ecount(g) * 10, 2000))),
                        error = function(e) NULL)
    if (is.null(g_null)) { null_means[i] <- NA; next }
    E(g_null)$weight <- sample(edge_df$weight)   # reshuffle the same weight values onto rewired topology
    null_means[i] <- mean(E(g_null)$weight)
  }
  null_mean <- mean(null_means, na.rm = TRUE)
  cat(sprintf("[%s] null-model mean edge weight (%d permutations): %.4f\n", cond, n_perm, null_mean))

  E(g)$weight_corrected <- E(g)$weight - null_mean
  g <- delete_edges(g, which(abs(E(g)$weight_corrected) < 0.01))
  cat(sprintf("[%s] n edges after null-model correction: %d\n", cond, ecount(g)))
  if (ecount(g) == 0) return(NULL)

  ## --- Louvain community detection ---
  comm <- cluster_louvain(g, weights = abs(E(g)$weight_corrected))
  V(g)$community <- membership(comm)

  ## --- Eigenvector centrality, keep top 100 nodes for visualization ---
  eig <- eigen_centrality(g, weights = abs(E(g)$weight_corrected))$vector
  V(g)$eigen_centrality <- eig

  list(graph = g, n_communities = length(unique(membership(comm))), null_mean = null_mean)
}

results <- list()
for (cond in conditions) {
  results[[cond]] <- build_network(cond)
}

## --- Summary ---
summary_tbl <- lapply(names(results), function(cond) {
  r <- results[[cond]]
  if (is.null(r)) return(data.frame(condition = cond, n_nodes = 0, n_edges = 0, n_communities = 0, null_mean = NA))
  data.frame(condition = cond, n_nodes = vcount(r$graph), n_edges = ecount(r$graph),
             n_communities = r$n_communities, null_mean = r$null_mean)
}) %>% bind_rows()
print(summary_tbl)
write.csv(summary_tbl, "~/Desktop/nullmodel_corrected_network_summary.csv", row.names = FALSE)

## --- Save hub tables (top nodes by eigenvector centrality, with community) ---
hub_tables <- lapply(names(results), function(cond) {
  r <- results[[cond]]
  if (is.null(r)) return(NULL)
  data.frame(condition = cond, region = V(r$graph)$name,
             community = V(r$graph)$community, eigen_centrality = V(r$graph)$eigen_centrality) %>%
    arrange(desc(eigen_centrality))
}) %>% bind_rows()
write.csv(hub_tables, "~/Desktop/nullmodel_corrected_network_hubs.csv", row.names = FALSE)

saveRDS(results, "~/Desktop/nullmodel_corrected_networks.rds")

## --- Plot each network (top 100 nodes by eigenvector centrality) ---
plot_network <- function(cond, r) {
  if (is.null(r)) return(NULL)
  g <- r$graph
  eig_named <- setNames(V(g)$eigen_centrality, V(g)$name)
  top_nodes <- names(sort(eig_named, decreasing = TRUE))[1:min(100, vcount(g))]
  g_sub <- induced_subgraph(g, top_nodes)
  # Fruchterman-Reingold requires non-negative weights; drop negative-weight
  # edges before layout so only positive co-activation is shown/laid out on.
  g_sub <- delete_edges(g_sub, E(g_sub)[weight_corrected <= 0])
  # ggraph's "fr" layout auto-detects an edge attribute literally named
  # "weight" and uses it for the layout unless told otherwise; that is
  # still the raw (possibly negative) Spearman r, not weight_corrected.
  # Drop it so the layout falls back to unweighted FR -- the actual
  # correlation strength/sign is encoded visually via aes() below, not
  # via the layout algorithm.
  g_sub <- delete_edge_attr(g_sub, "weight")

  tg <- as_tbl_graph(g_sub) %>%
    activate(nodes) %>% mutate(community = factor(community))

  ggraph(tg, layout = "fr") +
    geom_edge_link(aes(width = abs(weight_corrected), color = weight_corrected > 0), alpha = 0.4) +
    geom_node_point(aes(size = eigen_centrality, color = community)) +
    geom_node_text(aes(label = name), size = 2, repel = TRUE, max.overlaps = 30) +
    scale_edge_width(range = c(0.1, 1.2), guide = "none") +
    scale_edge_color_manual(values = c("TRUE" = "#C2185B", "FALSE" = "#8E24AA"), guide = "none") +
    scale_size(range = c(1, 6), guide = "none") +
    labs(title = sprintf("Co-activation network: %s", cond),
         subtitle = sprintf("Top %d nodes by eigenvector centrality, colored by Louvain community", min(100, vcount(g)))) +
    theme_graph(base_family = "sans") +
    theme(legend.position = "bottom")
}

for (cond in conditions) {
  p <- plot_network(cond, results[[cond]])
  if (is.null(p)) next
  fname <- sprintf("~/Desktop/nullmodel_corrected_network_%s.png", gsub("[^A-Za-z0-9]+", "_", cond))
  ggsave(fname, p, width = 10, height = 9, units = "in", dpi = 300)
  cat("Saved:", fname, "\n")
}


# ---------------------------------------------------------------------------
# nullmodel_corrected_coactivation_networks_normalized
# ---------------------------------------------------------------------------
# NORMALIZED version of nullmodel_corrected_coactivation_networks.R: each animal's
# per-region value is first divided by that animal's own whole-brain median
# (same per-animal normalization used earlier), removing each animal's
# global scaling before computing Spearman correlations. This tests how much
# of the raw network's structure (null-model mean edge weight 0.74-0.96) was
# driven by the global-scale confound rather than genuine co-activation.
#
# Co-activation network analysis per stressor condition, using a null-model-corrected approach
# (2026):
#   1. Pairwise Spearman correlations between regions, across the animals
#      belonging to one condition.
#   2. Edges retained only if p < 0.01 (unadjusted) AND >= 3 paired obs.
#   3. Null-model correction: 1,000 degree-preserving network permutations
#      (igraph edge rewiring, preserving each node's degree); the mean edge
#      weight across these permutations is subtracted from every observed
#      edge weight. Edges with negligible corrected weight are dropped.
#   4. Louvain community detection on the corrected network.
#   5. Visualization: top 100 nodes by eigenvector centrality.
#
# One network per condition (5 total). Regions restricted to those with
# n >= 3 valid observations within that condition (same rule applied
# for region inclusion). 'FS' excluded (injection region, trivially 1.0).

library(dplyr)
library(tidyr)
library(igraph)
library(ggraph)
library(tidygraph)
library(ggplot2)

df <- read.csv("~/Desktop/ALL_DT.csv", stringsAsFactors = FALSE)
cond_map <- read.csv("~/Desktop/group_condition_corrected.csv", stringsAsFactors = FALSE)

## --- Node set: restricted to a meaningful region subset, analogous to
## the "FDR-significant regions from Level 1" -- since no region survived
## FDR here (see nullmodel_corrected_pairwise_stats.csv), we substitute the top
## candidates across all 10 pairwise comparisons (lowest p-value per region,
## union across comparisons) plus the previously curated top/most-variable
## regions. This keeps the network at a computationally tractable size
## (all-region networks produced 100k+ edges -- infeasible for 1000
## permutations) while still reflecting "regions of interest" rather than
## an arbitrary/random subset.
pairwise_stats <- read.csv("~/Desktop/nullmodel_corrected_pairwise_stats.csv", stringsAsFactors = FALSE)
top_by_pvalue <- pairwise_stats %>% group_by(region_acronym) %>%
  summarise(min_p = min(p_value), .groups = "drop") %>%
  arrange(min_p) %>% slice_head(n = 80) %>% pull(region_acronym)

top_regions_df <- read.csv("~/Desktop/top_regions_low_deviation.csv", stringsAsFactors = FALSE)
most_var_df <- read.csv("~/Desktop/most_variable_regions.csv", stringsAsFactors = FALSE)
node_regions <- unique(c(top_by_pvalue, top_regions_df$region_acronym, most_var_df$region_acronym))
cat("n node regions for network analysis:", length(node_regions), "\n")

long <- df %>%
  filter(region_acronym != "FS", region_acronym %in% node_regions) %>%
  select(animal_group, animal_id_original, region_acronym,
         trapcre_intensity, traptta_intensity) %>%
  pivot_longer(cols = c(trapcre_intensity, traptta_intensity),
               names_to = "channel", values_to = "value") %>%
  mutate(channel = ifelse(channel == "trapcre_intensity", "Cre", "tTA"),
         sample_id = paste(animal_id_original, channel, sep = "_")) %>%
  rename(Group = animal_group, Channel = channel) %>%
  left_join(cond_map, by = c("Group", "Channel")) %>%
  filter(!is.na(Condition))

## --- Normalize: divide by each animal x channel's own whole-brain median
## (computed across ALL regions in ALL_DT.csv, not just the restricted node
## set, so the normalization factor reflects the animal's true global level) ---
animal_overall <- df %>%
  select(animal_id_original, trapcre_intensity, traptta_intensity) %>%
  pivot_longer(cols = c(trapcre_intensity, traptta_intensity), names_to = "channel", values_to = "value") %>%
  mutate(channel = ifelse(channel == "trapcre_intensity", "Cre", "tTA")) %>%
  group_by(animal_id_original, channel) %>%
  summarise(overall_median = median(value, na.rm = TRUE), .groups = "drop")

long <- long %>%
  left_join(animal_overall, by = c("animal_id_original", "Channel" = "channel")) %>%
  mutate(value = value / overall_median)

conditions <- sort(unique(long$Condition))

## --- Fast Spearman correlation + p-value matrix (t-approximation) ---
spearman_with_p <- function(mat) {
  r <- cor(mat, method = "spearman", use = "pairwise.complete.obs")
  not_na <- !is.na(mat)
  n_pairs <- crossprod(not_na)   # n valid paired observations per region pair
  t_stat <- r * sqrt((n_pairs - 2) / (1 - r^2))
  p <- 2 * pt(-abs(t_stat), df = n_pairs - 2)
  diag(p) <- 1
  list(r = r, p = p, n = n_pairs)
}

build_network <- function(cond, n_perm = 1000, min_obs = 3) {
  d <- long %>% filter(Condition == cond)
  mat <- d %>% select(sample_id, region_acronym, value) %>%
    pivot_wider(names_from = region_acronym, values_from = value) %>%
    tibble::column_to_rownames("sample_id") %>% as.matrix()

  ## keep only regions with >= min_obs non-NA values
  keep <- colSums(!is.na(mat)) >= min_obs
  mat <- mat[, keep, drop = FALSE]
  cat(sprintf("[%s] n samples = %d, n regions (n>=%d) = %d\n", cond, nrow(mat), min_obs, ncol(mat)))

  sp <- spearman_with_p(mat)

  ## --- Threshold edges: p < 0.01 unadjusted, n_pairs >= 3 ---
  edge_mask <- (sp$p < 0.01) & (sp$n >= min_obs) & upper.tri(sp$r)
  edges <- which(edge_mask, arr.ind = TRUE)
  if (nrow(edges) == 0) {
    cat(sprintf("[%s] No edges survived thresholding.\n", cond))
    return(NULL)
  }
  edge_df <- data.frame(
    from = colnames(mat)[edges[, 1]], to = colnames(mat)[edges[, 2]],
    weight = sp$r[edges]
  )
  cat(sprintf("[%s] n edges after thresholding: %d\n", cond, nrow(edge_df)))

  g <- graph_from_data_frame(edge_df, directed = FALSE)

  ## --- Null model: 1000 degree-preserving rewirings; subtract mean null
  ## edge weight (using the same weight multiset, reshuffled onto rewired
  ## topology) from every observed edge weight ---
  null_means <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    g_null <- tryCatch(rewire(g, with = keeping_degseq(niter = min(ecount(g) * 10, 2000))),
                        error = function(e) NULL)
    if (is.null(g_null)) { null_means[i] <- NA; next }
    E(g_null)$weight <- sample(edge_df$weight)   # reshuffle the same weight values onto rewired topology
    null_means[i] <- mean(E(g_null)$weight)
  }
  null_mean <- mean(null_means, na.rm = TRUE)
  cat(sprintf("[%s] null-model mean edge weight (%d permutations): %.4f\n", cond, n_perm, null_mean))

  E(g)$weight_corrected <- E(g)$weight - null_mean
  g <- delete_edges(g, which(abs(E(g)$weight_corrected) < 0.01))
  cat(sprintf("[%s] n edges after null-model correction: %d\n", cond, ecount(g)))
  if (ecount(g) == 0) return(NULL)

  ## --- Louvain community detection ---
  comm <- cluster_louvain(g, weights = abs(E(g)$weight_corrected))
  V(g)$community <- membership(comm)

  ## --- Eigenvector centrality, keep top 100 nodes for visualization ---
  eig <- eigen_centrality(g, weights = abs(E(g)$weight_corrected))$vector
  V(g)$eigen_centrality <- eig

  list(graph = g, n_communities = length(unique(membership(comm))), null_mean = null_mean)
}

results <- list()
for (cond in conditions) {
  results[[cond]] <- build_network(cond)
}

## --- Summary ---
summary_tbl <- lapply(names(results), function(cond) {
  r <- results[[cond]]
  if (is.null(r)) return(data.frame(condition = cond, n_nodes = 0, n_edges = 0, n_communities = 0, null_mean = NA))
  data.frame(condition = cond, n_nodes = vcount(r$graph), n_edges = ecount(r$graph),
             n_communities = r$n_communities, null_mean = r$null_mean)
}) %>% bind_rows()
print(summary_tbl)
write.csv(summary_tbl, "~/Desktop/nullmodel_corrected_network_summary_normalized.csv", row.names = FALSE)

## --- Save hub tables (top nodes by eigenvector centrality, with community) ---
hub_tables <- lapply(names(results), function(cond) {
  r <- results[[cond]]
  if (is.null(r)) return(NULL)
  data.frame(condition = cond, region = V(r$graph)$name,
             community = V(r$graph)$community, eigen_centrality = V(r$graph)$eigen_centrality) %>%
    arrange(desc(eigen_centrality))
}) %>% bind_rows()
write.csv(hub_tables, "~/Desktop/nullmodel_corrected_network_hubs_normalized.csv", row.names = FALSE)

saveRDS(results, "~/Desktop/nullmodel_corrected_networks_normalized.rds")

## --- Plot each network (top 100 nodes by eigenvector centrality) ---
plot_network <- function(cond, r) {
  if (is.null(r)) return(NULL)
  g <- r$graph
  eig_named <- setNames(V(g)$eigen_centrality, V(g)$name)
  top_nodes <- names(sort(eig_named, decreasing = TRUE))[1:min(100, vcount(g))]
  g_sub <- induced_subgraph(g, top_nodes)

  tg <- as_tbl_graph(g_sub) %>%
    activate(nodes) %>% mutate(community = factor(community))

  ggraph(tg, layout = "fr", weights = abs(weight_corrected)) +
    geom_edge_link(aes(width = abs(weight_corrected), color = weight_corrected > 0), alpha = 0.4) +
    geom_node_point(aes(size = eigen_centrality, color = community)) +
    geom_node_text(aes(label = name), size = 2, repel = TRUE, max.overlaps = 30) +
    scale_edge_width(range = c(0.1, 1.2), guide = "none") +
    scale_edge_color_manual(values = c("TRUE" = "#C2185B", "FALSE" = "#8E24AA"), guide = "none") +
    scale_size(range = c(1, 6), guide = "none") +
    labs(title = sprintf("Co-activation network (normalized): %s", cond),
         subtitle = sprintf("Top %d nodes by eigenvector centrality, colored by Louvain community", min(100, vcount(g)))) +
    theme_graph(base_family = "sans") +
    theme(legend.position = "bottom")
}

for (cond in conditions) {
  p <- plot_network(cond, results[[cond]])
  if (is.null(p)) next
  fname <- sprintf("~/Desktop/nullmodel_corrected_network_normalized_%s.png", gsub("[^A-Za-z0-9]+", "_", cond))
  ggsave(fname, p, width = 10, height = 9, units = "in", dpi = 300)
  cat("Saved:", fname, "\n")
}


# ---------------------------------------------------------------------------
# sex_coactivation_networks
# ---------------------------------------------------------------------------
# Per-sex coactivation networks: same method as condition_per_condition_networks.R
# (Spearman |rho|>0.6 across animal x channel samples, degree/eigenvector
# centrality), but applied separately to each SEX's own top-40 regions and
# own animals -- the sex-specific counterpart to the per-condition networks,
# to compare against sex_region_clustering_zscore.R's per-sex dendrograms.

library(dplyr)
library(tidyr)
library(igraph)
library(ggraph)
library(ggplot2)
library(patchwork)

setwd("~/Desktop/Bachelor")
source("prepare_heatmap_data.R")
long <- long %>% filter(region_acronym != "FS")

top_n_regions <- 40
cor_thresh <- 0.6
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

hub_results <- list()
plots <- list()

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

  cor_mat <- cor(t(mat_sx_top), method = "spearman", use = "pairwise.complete.obs")
  adj <- (abs(cor_mat) > cor_thresh) & (abs(cor_mat) < 1)
  diag(adj) <- FALSE

  g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
  deg <- as.numeric(igraph::degree(g)); names(deg) <- igraph::V(g)$name
  eig <- tryCatch(igraph::eigen_centrality(g)$vector, error = function(e) rep(NA, igraph::vcount(g)))

  divs <- major_division_per_row[match(names(deg), rownames(mat))]

  hub_df <- data.frame(
    sex = sx, region_acronym = names(deg), major_division = divs,
    degree = deg, eigen_centrality = round(eig, 3), n_samples = ncol(mat_sx_top)
  ) %>% arrange(desc(degree), desc(eigen_centrality))
  hub_results[[sx]] <- hub_df

  cat(sprintf("\n=== %s (%d samples, %d regions, %d edges at |rho|>%.1f) -- top 5 hubs ===\n",
              sx, ncol(mat_sx_top), length(deg), sum(adj) / 2, cor_thresh))
  print(head(hub_df, 5))

  igraph::V(g)$degree <- deg
  igraph::V(g)$division <- as.character(divs)
  p <- ggraph(g, layout = "fr") +
    geom_edge_link(alpha = 0.25, color = "grey50") +
    geom_node_point(aes(size = degree, color = division)) +
    geom_node_text(aes(label = name), repel = TRUE, size = 2.3, max.overlaps = 30) +
    scale_color_manual(values = major_cols, name = "Major division") +
    scale_size_continuous(range = c(1, 5), guide = "none") +
    labs(title = sx) +
    theme_void(base_size = 9) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11), legend.position = "none")
  plots[[sx]] <- p
}

hub_all <- bind_rows(hub_results)
write.csv(hub_all, "sex_coactivation_network_hubs.csv", row.names = FALSE)

cat("\n\n=== Top hub per sex (by degree) ===\n")
print(hub_all %>% group_by(sex) %>% slice_max(degree, n = 1, with_ties = TRUE) %>%
        select(sex, region_acronym, major_division, degree, eigen_centrality))

## Combined figure: Female + Male networks + shared legend
legend_plot <- ggplot(data.frame(x = 1, y = seq_along(major_cols), division = names(major_cols)),
                       aes(x, y, color = division)) +
  geom_point(size = 3) +
  scale_color_manual(values = major_cols, name = "Major division") +
  theme_void() +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(legend.position = "right", legend.title = element_text(face = "bold"))
legend_only <- cowplot::get_legend(legend_plot)

combined <- wrap_plots(c(plots, list(cowplot::ggdraw(legend_only))), ncol = 3) +
  plot_annotation(title = "Per-sex coactivation networks (top 40 regions, |rho| > 0.6)",
                   theme = theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5)))

ggsave("sex_coactivation_networks_combined.png", combined, width = 13, height = 5, units = "in", dpi = 300)

cat("\nSaved: sex_coactivation_network_hubs.csv, sex_coactivation_networks_combined.png\n")


# ---------------------------------------------------------------------------
# sex_lowdev_network
# ---------------------------------------------------------------------------
# Coactivation network per SEX, restricted to that sex's own top-40 "low
# deviation, high intensity" regions (sex_top_regions_low_deviation.csv),
# pooling animals across all 4 (non-chronic-FS) conditions -- sex-specific
# counterpart to condition_lowdev_network.R (which pooled both sexes and
# all 5 conditions). Same method: Spearman |rho|>0.6, whole-brain-median
# normalized, across all animal x channel samples of that sex.

library(dplyr)
library(tidyr)
library(igraph)
library(ggraph)
library(ggplot2)
library(patchwork)

setwd("~/Desktop/Bachelor")
source("prepare_heatmap_data.R")
long <- long %>% filter(region_acronym != "FS")

sex_order <- c("Female", "Male")
condition_order_no_chronicFS <- setdiff(condition_order, "social defeat · Acute")  # social defeat is now single-sex (Male only); excluded here so only conditions comparable across both sexes remain. Chronic FS is included since its sex labels were corrected.
top40_df <- read.csv("sex_top_regions_low_deviation.csv", stringsAsFactors = FALSE) %>%
  group_by(sex) %>% slice_head(n = 40) %>% ungroup()

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
cor_thresh <- 0.6

hub_results <- list()
plots <- list()
FONT_PT <- 14
mm_size <- function(pt) pt / (72.27 / 25.4)  # ggplot geom text 'size' aesthetic is in mm

for (sx in sex_order) {
  top40_sx <- top40_df %>% filter(sex == sx)

  sample_overall <- long %>%
    filter(sex == sx, condition %in% condition_order_no_chronicFS) %>%
    group_by(animal_id_original, channel) %>%
    summarise(overall_median = median(intensity, na.rm = TRUE), .groups = "drop")

  long_norm <- long %>%
    filter(sex == sx, condition %in% condition_order_no_chronicFS) %>%
    left_join(sample_overall, by = c("animal_id_original", "channel")) %>%
    mutate(norm_value = intensity / overall_median,
           sample_id = paste(animal_id_original, channel, sep = "_"))

  wide <- long_norm %>%
    filter(region_acronym %in% top40_sx$region_acronym) %>%
    select(sample_id, region_acronym, norm_value) %>%
    pivot_wider(names_from = region_acronym, values_from = norm_value)

  mat_samples <- as.matrix(wide[, -1])
  rownames(mat_samples) <- wide$sample_id
  mat_samples <- mat_samples[, colSums(!is.na(mat_samples)) >= 10, drop = FALSE]

  cor_mat <- cor(mat_samples, method = "spearman", use = "pairwise.complete.obs")
  adj <- (abs(cor_mat) > cor_thresh) & (abs(cor_mat) < 1)
  diag(adj) <- FALSE

  g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
  deg <- as.numeric(igraph::degree(g)); names(deg) <- igraph::V(g)$name
  eig <- tryCatch(igraph::eigen_centrality(g)$vector, error = function(e) rep(NA, igraph::vcount(g)))
  divs <- top40_sx$major_division[match(names(deg), top40_sx$region_acronym)]

  hub_df <- data.frame(sex = sx, region_acronym = names(deg), major_division = divs,
                        degree = deg, eigen_centrality = round(eig, 3),
                        n_samples = nrow(mat_samples)) %>%
    arrange(desc(degree), desc(eigen_centrality))
  hub_results[[sx]] <- hub_df

  cat(sprintf("=== %s low-deviation network (n=%d samples, %d regions, %d edges at |rho|>%.1f) ===\n",
              sx, nrow(mat_samples), length(deg), sum(adj) / 2, cor_thresh))
  print(head(hub_df, 8))

  isolated <- names(deg)[deg <= 1]
  g_connected <- igraph::delete_vertices(g, isolated)
  igraph::V(g_connected)$degree <- deg[igraph::V(g_connected)$name]
  igraph::V(g_connected)$division <- as.character(divs[match(igraph::V(g_connected)$name, names(deg))])
  cat(sprintf("Isolated (degree<=1, excluded from plot): %s\n\n", paste(isolated, collapse = ", ")))

  p <- ggraph(g_connected, layout = "fr") +
    geom_edge_link(alpha = 0.35, color = "grey55", width = 0.6) +
    geom_node_point(aes(size = degree, color = degree)) +
    geom_node_text(aes(label = name), nudge_y = 0.22, vjust = 0, size = mm_size(FONT_PT),
                   color = "#2A2A2A", fontface = "bold") +
    scale_color_gradientn(colours = pink_scale, name = "Degree") +
    scale_size_continuous(range = c(3, 10), guide = "none") +
    scale_x_continuous(expand = expansion(mult = 0.25)) +
    scale_y_continuous(expand = expansion(mult = 0.25)) +
    coord_fixed() +
    labs(title = sprintf("%s (n=%d)", sx, nrow(mat_samples))) +
    theme_void(base_size = FONT_PT) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = FONT_PT + 2),
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          legend.position = "bottom",
          legend.text = element_text(size = FONT_PT),
          legend.title = element_text(size = FONT_PT))

  ggsave(sprintf("sex_lowdev_network_%s.png", sx), p, width = 12, height = 12, units = "in", dpi = 300, limitsize = FALSE)
  plots[[sx]] <- p
}

hub_all <- bind_rows(hub_results)
write.csv(hub_all, "sex_lowdev_network_hubs.csv", row.names = FALSE)

combined <- wrap_plots(plots, ncol = 2) +
  plot_annotation(title = "Low-deviation-region coactivation networks by sex (top 40 regions, |rho| > 0.6)",
                   theme = theme(plot.title = element_text(face = "bold", size = FONT_PT + 2, hjust = 0.5)))
ggsave("sex_lowdev_network_combined.png", combined, width = 20, height = 12, units = "in", dpi = 300, limitsize = FALSE)

cat("\nSaved: sex_lowdev_network_hubs.csv, sex_lowdev_network_Female.png, sex_lowdev_network_Male.png, sex_lowdev_network_combined.png\n")


# ---------------------------------------------------------------------------
# coactive_projection_regions
# ---------------------------------------------------------------------------
# Top "coactive projection regions" per stressor condition: the hub nodes
# (highest eigenvector centrality) from the null-model-corrected co-activation
# networks (nullmodel_corrected_coactivation_networks(.R/_normalized.R)).
#
# Eigenvector centrality = how strongly & how centrally a region is
# co-activated with other well-connected regions in that condition's
# network, after null-model correction. High centrality = a "hub" region
# whose activity pattern is tightly coupled to the rest of the network for
# that particular stressor.
#
# Shown for BOTH raw and normalized networks (raw dominated by the global
# scale confound, normalized reflects genuine regional co-activation -- see
# nullmodel_corrected_coactivation_networks_normalized.R).

library(dplyr)
library(ggplot2)

df <- read.csv("~/Desktop/ALL_DT.csv", stringsAsFactors = FALSE)
structures <- read.csv("~/Desktop/tractquant/allen_mouse_25um_v1.2/structures.csv", stringsAsFactors = FALSE)

## --- Major anatomical division (Allen CCF), same derivation as
## prepare_heatmap_data.R -- NOT the coarse 3-value ALL_DT.csv$division col ---
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
  distinct(region_id, region_acronym, region_name) %>%
  left_join(structures %>% select(id, structure_id_path), by = c("region_id" = "id")) %>%
  rowwise() %>%
  mutate(division = get_major_division(structure_id_path)) %>%
  ungroup() %>%
  select(region_acronym, region_name, division)

hubs_raw  <- read.csv("~/Desktop/nullmodel_corrected_network_hubs.csv", stringsAsFactors = FALSE)
hubs_norm <- read.csv("~/Desktop/nullmodel_corrected_network_hubs_normalized.csv", stringsAsFactors = FALSE)

top_n <- 10

make_top_table <- function(hubs, version_label) {
  hubs %>%
    group_by(condition) %>%
    arrange(desc(eigen_centrality), .by_group = TRUE) %>%
    slice_head(n = top_n) %>%
    mutate(rank = row_number()) %>%
    ungroup() %>%
    left_join(region_lookup, by = c("region" = "region_acronym")) %>%
    transmute(condition, rank, region_acronym = region, region_name, division,
              community, eigen_centrality = round(eigen_centrality, 3),
              version = version_label)
}

top_raw  <- make_top_table(hubs_raw,  "raw")
top_norm <- make_top_table(hubs_norm, "normalized")

combined <- bind_rows(top_raw, top_norm)
write.csv(combined, "~/Desktop/coactive_projection_regions_top10.csv", row.names = FALSE)

## --- Colored table for the thesis (normalized version, one table per
## condition, pink-violet palette matching the project style), built with
## ggplot2 (no headless-Chrome dependency, unlike gt::gtsave/webshot2) ---
pal_fun <- colorRampPalette(c("#FDF0F5", "#F48FB1", "#C2185B", "#8E24AA"))

plot_table <- function(d, cond) {
  d <- d %>% arrange(rank) %>%
    mutate(region_name = ifelse(is.na(region_name), "", region_name),
           division = ifelse(is.na(division), "", division))

  d$y <- nrow(d):1
  header <- data.frame(x = c(0.5, 3.2, 9.4, 12.2, 13.6),
                        label = c("Region", "Full name", "Division", "Comm.", "Eig. cent."),
                        y = nrow(d) + 1)

  ggplot(d, aes(xmin = 0, xmax = 15, ymin = y - 0.5, ymax = y + 0.5)) +
    geom_rect(aes(fill = eigen_centrality), color = "white", linewidth = 0.6) +
    scale_fill_gradientn(colors = pal_fun(100), limits = range(top_norm$eigen_centrality), guide = "none") +
    geom_text(aes(x = 0.5, y = y, label = rank), hjust = 0, size = 3.2, fontface = "bold") +
    geom_text(aes(x = 1.4, y = y, label = region_acronym), hjust = 0, size = 3.2, fontface = "bold") +
    geom_text(aes(x = 3.2, y = y, label = region_name), hjust = 0, size = 2.6) +
    geom_text(aes(x = 9.4, y = y, label = division), hjust = 0, size = 2.6) +
    geom_text(aes(x = 12.2, y = y, label = community), hjust = 0.5, size = 3) +
    geom_text(aes(x = 13.6, y = y, label = sprintf("%.3f", eigen_centrality)), hjust = 0.5, size = 3) +
    geom_text(data = header, aes(x = x, y = y, label = label), inherit.aes = FALSE,
              hjust = 0, size = 3.2, fontface = "bold", color = "#4A148C") +
    coord_cartesian(xlim = c(0, 15), ylim = c(0.3, nrow(d) + 1.5), clip = "off") +
    labs(title = sprintf("Coactive projection regions: %s", cond),
         subtitle = "Top 10 hub regions by eigenvector centrality (normalized network)") +
    theme_void() +
    theme(plot.title = element_text(face = "bold", size = 13, hjust = 0, margin = margin(b = 4)),
          plot.subtitle = element_text(size = 9, color = "grey30", margin = margin(b = 8)),
          plot.margin = margin(10, 14, 10, 14))
}

for (cond in unique(top_norm$condition)) {
  d <- top_norm %>% filter(condition == cond) %>%
    select(rank, region_acronym, region_name, division, community, eigen_centrality)

  p <- plot_table(d, cond)
  fname_safe <- gsub("[^A-Za-z0-9]+", "_", cond)
  ggsave(sprintf("~/Desktop/coactive_projection_regions_%s.png", fname_safe), p,
         width = 9, height = 5, units = "in", dpi = 300)
  cat("Saved table for:", cond, "\n")
}

cat("\n\n=== TOP", top_n, "COACTIVE (NORMALIZED) REGIONS PER CONDITION ===\n")
for (cond in unique(top_norm$condition)) {
  cat("\n---", cond, "---\n")
  print(top_norm %>% filter(condition == cond) %>%
          select(rank, region_acronym, division, community, eigen_centrality) %>% as.data.frame())
}


# ---------------------------------------------------------------------------
# fs_acute_vs_chronic_female_networks
# ---------------------------------------------------------------------------
# Coactivation networks (Spearman |rho|>0.6, top 40 regions) for forced
# swim Acute vs Chronic, restricted to FEMALE animals only -- same method as
# condition_per_condition_networks.R / sex_coactivation_networks.R, but
# crossing condition x sex to isolate the acute-vs-chronic (duration)
# comparison from the sex confound (chronic forced swim is all-female,
# acute forced swim is mostly male -- see fs_acute_vs_chronic_female_only.R
# for the corresponding statistical tests).
#
# CAVEAT: forced swim Acute Female has only n=3 animal x channel samples.
# With n=3, Spearman rho can only take the values {-1, -0.5, 0.5, 1}, so
# |rho|>0.6 effectively means "rho = +-1" -- any two regions with a
# perfectly monotonic (or perfectly reversed) rank across just 3 points
# will pass threshold. This network is therefore expected to be extremely
# dense and should be treated as illustrative/exploratory only, not as
# evidence of genuine coactivation. Chronic Female (n=9) is more (though
# still sparsely) powered.

library(dplyr)
library(tidyr)
library(igraph)
library(ggraph)
library(ggplot2)
library(patchwork)

setwd("~/Desktop/Bachelor")
source("prepare_heatmap_data.R")
long <- long %>% filter(region_acronym != "FS")

top_n_regions <- 40
cor_thresh <- 0.6
cells <- list(
  c(condition = "forced swim · Acute",   sex = "Female"),
  c(condition = "forced swim · Chronic", sex = "Female")
)

major_cols <- c(
  "Isocortex"             = "#5B7C99", "Olfactory areas"       = "#8C6D46",
  "Hippocampal formation" = "#6B8F71", "Cortical subplate"     = "#9C6B85",
  "Striatum"              = "#4F8F8B", "Pallidum"              = "#7A6B9C",
  "Thalamus"              = "#B08B5A", "Hypothalamus"          = "#A85C5C",
  "Midbrain"              = "#6B6B99", "Pons"                  = "#8FA05E",
  "Medulla"               = "#5E8FA0", "Cerebellum"            = "#9C7A5B",
  "fiber tracts"          = "#A9A9A9", "ventricular systems"   = "#C4B7A6"
)

hub_results <- list()
plots <- list()

for (cell in cells) {
  cond <- cell["condition"]; sx <- cell["sex"]
  label <- paste(cond, sx, sep = " | ")

  d <- long %>% filter(condition == cond, sex == sx) %>%
    select(region_acronym, animal_id_original, channel, intensity) %>%
    mutate(sample_id = paste(animal_id_original, channel, sep = "_"))
  wide_cell <- d %>% select(region_acronym, sample_id, intensity) %>%
    pivot_wider(names_from = sample_id, values_from = intensity)
  mat_cell <- as.matrix(wide_cell[, -1])
  rownames(mat_cell) <- wide_cell$region_acronym
  mat_cell <- mat_cell[rowSums(is.na(mat_cell)) == 0, , drop = FALSE]

  med_per_region <- apply(mat_cell, 1, median, na.rm = TRUE)
  n_top <- min(top_n_regions, length(med_per_region))
  top_regions <- names(sort(med_per_region, decreasing = TRUE))[seq_len(n_top)]
  mat_cell_top <- mat_cell[top_regions, , drop = FALSE]

  cor_mat <- cor(t(mat_cell_top), method = "spearman", use = "pairwise.complete.obs")
  adj <- (abs(cor_mat) > cor_thresh) & (abs(cor_mat) < 1)
  diag(adj) <- FALSE

  g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
  deg <- as.numeric(igraph::degree(g)); names(deg) <- igraph::V(g)$name
  eig <- tryCatch(igraph::eigen_centrality(g)$vector, error = function(e) rep(NA, igraph::vcount(g)))

  divs <- major_division_per_row[match(names(deg), rownames(mat))]

  hub_df <- data.frame(
    condition = cond, sex = sx, region_acronym = names(deg), major_division = divs,
    degree = deg, eigen_centrality = round(eig, 3), n_samples = ncol(mat_cell_top)
  ) %>% arrange(desc(degree), desc(eigen_centrality))
  hub_results[[label]] <- hub_df

  cat(sprintf("\n=== %s (%d samples, %d regions, %d edges at |rho|>%.1f, of max %d) -- top 5 hubs ===\n",
              label, ncol(mat_cell_top), length(deg), sum(adj) / 2, cor_thresh, n_top * (n_top - 1) / 2))
  print(head(hub_df, 5))

  igraph::V(g)$degree <- deg
  igraph::V(g)$division <- as.character(divs)
  p <- ggraph(g, layout = "fr") +
    geom_edge_link(alpha = 0.25, color = "grey50") +
    geom_node_point(aes(size = degree, color = division)) +
    geom_node_text(aes(label = name), repel = TRUE, size = 2.3, max.overlaps = 30) +
    scale_color_manual(values = major_cols, name = "Major division") +
    scale_size_continuous(range = c(1, 5), guide = "none") +
    labs(title = sprintf("%s (%s, n=%d)", cond, sx, ncol(mat_cell_top))) +
    theme_void(base_size = 9) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10.5), legend.position = "none")
  plots[[label]] <- p
}

hub_all <- bind_rows(hub_results)
write.csv(hub_all, "fs_acute_vs_chronic_FEMALE_network_hubs.csv", row.names = FALSE)

cat("\n\n=== Top hub per condition (Female only, by degree) ===\n")
print(hub_all %>% group_by(condition) %>% slice_max(degree, n = 1, with_ties = TRUE) %>%
        select(condition, region_acronym, major_division, degree, eigen_centrality))

## Combined figure: Acute + Chronic (Female) networks + shared legend
legend_plot <- ggplot(data.frame(x = 1, y = seq_along(major_cols), division = names(major_cols)),
                       aes(x, y, color = division)) +
  geom_point(size = 3) +
  scale_color_manual(values = major_cols, name = "Major division") +
  theme_void() +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(legend.position = "right", legend.title = element_text(face = "bold"))
legend_only <- cowplot::get_legend(legend_plot)

combined <- wrap_plots(c(plots, list(cowplot::ggdraw(legend_only))), ncol = 3) +
  plot_annotation(title = "Forced swim Acute vs Chronic coactivation networks (Female only, top 40 regions, |rho| > 0.6)",
                   theme = theme(plot.title = element_text(face = "bold", size = 12.5, hjust = 0.5)))

ggsave("fs_acute_vs_chronic_FEMALE_networks_combined.png", combined, width = 13, height = 5, units = "in", dpi = 300)

cat("\nSaved: fs_acute_vs_chronic_FEMALE_network_hubs.csv, fs_acute_vs_chronic_FEMALE_networks_combined.png\n")

