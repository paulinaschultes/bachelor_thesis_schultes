# =============================================================================
# Venn diagrams of shared top regions
#
#   Figure 26 — top 40 high-intensity, low-variability regions per sex
#   Figure 29 — top 20 candidate input regions per retrograde cohort
#
# Both diagrams show the overlap of two ranked region sets. The ranking itself
# is computed elsewhere (combined rank of high median intensity and low
# coefficient of variation); this script only selects the top n per group,
# intersects them and draws the diagram.
#
# Requires: dplyr, ggplot2, ggvenn
#   install.packages(c("dplyr", "ggplot2", "ggvenn"))
# =============================================================================

library(dplyr)
library(ggplot2)
library(ggvenn)

base_dir <- "~/Desktop/Bachelor"
out_dir  <- base_dir

# Colours used throughout the thesis figures
col_pink   <- "#D6547F"
col_violet <- "#6B4C9A"


# -----------------------------------------------------------------------------
# Helper: take the n best regions of a group by combined rank
# -----------------------------------------------------------------------------
top_regions <- function(df, group_col, region_col, rank_col, n = 40) {
  df %>%
    group_by(.data[[group_col]]) %>%
    arrange(.data[[rank_col]], .by_group = TRUE) %>%
    slice_head(n = n) %>%
    ungroup() %>%
    select(group = all_of(group_col), region = all_of(region_col))
}


# -----------------------------------------------------------------------------
# Helper: draw a two-set Venn diagram and report the overlap
# -----------------------------------------------------------------------------
draw_venn <- function(sets, fill, title, file, width = 5.5, height = 4.5) {
  shared <- intersect(sets[[1]], sets[[2]])
  cat("\n", title, "\n", sep = "")
  cat("  ", names(sets)[1], ": ", length(sets[[1]]), " regions\n", sep = "")
  cat("  ", names(sets)[2], ": ", length(sets[[2]]), " regions\n", sep = "")
  cat("  shared: ", length(shared), "  -> ", paste(sort(shared), collapse = ", "),
      "\n", sep = "")

  p <- ggvenn(
    sets,
    fill_color      = fill,
    fill_alpha      = 0.45,
    stroke_color    = "grey25",
    stroke_size     = 0.5,
    set_name_size   = 5,
    text_size       = 5,
    show_percentage = FALSE
  ) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))

  ggsave(file.path(out_dir, file), p, width = width, height = height, dpi = 300)
  invisible(shared)
}


# =============================================================================
# Figure 26 — sex-specific top 40 regions
# =============================================================================
sex_df <- read.csv(file.path(base_dir, "sex_top_regions_low_deviation.csv"),
                   stringsAsFactors = FALSE)

sex_top <- top_regions(sex_df, "sex", "region_acronym", "combined", n = 40)

sex_sets <- list(
  Female = sex_top$region[sex_top$group == "Female"],
  Male   = sex_top$region[sex_top$group == "Male"]
)

shared_sex <- draw_venn(
  sex_sets,
  fill  = c(col_pink, col_violet),
  title = "Top 40 low-deviation regions per sex",
  file  = "venn_sex_top40.png"
)

write.csv(data.frame(region = sort(shared_sex)),
          file.path(out_dir, "venn_sex_top40_shared.csv"), row.names = FALSE)


# =============================================================================
# Figure 29 — retrograde cohorts, top 20 candidate input regions
# =============================================================================
ret1 <- read.csv(file.path(base_dir, "retrograde",
                           "input_regions_RET1_restraint_female_n3.csv"),
                 stringsAsFactors = FALSE)
ret3 <- read.csv(file.path(base_dir, "retrograde",
                           "input_regions_RET3_forcedswim_male_n2.csv"),
                 stringsAsFactors = FALSE)

ret_sets <- list(
  `Restraint (female)`   = ret1 %>% arrange(combined_rank) %>% slice_head(n = 20) %>% pull(acronym),
  `Forced swim (male)`   = ret3 %>% arrange(combined_rank) %>% slice_head(n = 20) %>% pull(acronym)
)

shared_ret <- draw_venn(
  ret_sets,
  fill  = c(col_pink, col_violet),
  title = "Top 20 candidate input regions per retrograde cohort",
  file  = "venn_retrograde_top20.png"
)

write.csv(data.frame(region = sort(shared_ret)),
          file.path(out_dir, "venn_retrograde_top20_shared.csv"), row.names = FALSE)


# =============================================================================
# Optional: proportional (area-weighted) version with eulerr
# Use this if the reviewer prefers circle areas to reflect set sizes.
# =============================================================================
if (requireNamespace("eulerr", quietly = TRUE)) {
  library(eulerr)

  euler_plot <- function(sets, fill, file) {
    a <- setdiff(sets[[1]], sets[[2]])
    b <- setdiff(sets[[2]], sets[[1]])
    ab <- intersect(sets[[1]], sets[[2]])
    counts <- c(length(a), length(b), length(ab))
    names(counts) <- c(names(sets)[1], names(sets)[2],
                       paste(names(sets)[1], names(sets)[2], sep = "&"))
    fit <- euler(counts)
    png(file.path(out_dir, file), width = 1600, height = 1300, res = 300)
    print(plot(fit, fills = list(fill = fill, alpha = 0.45),
               edges = list(col = "grey25"), quantities = TRUE))
    dev.off()
  }

  euler_plot(sex_sets, c(col_pink, col_violet), "venn_sex_top40_proportional.png")
  euler_plot(ret_sets, c(col_pink, col_violet), "venn_retrograde_top20_proportional.png")
}

cat("\nDone. Figures and shared-region tables written to ", out_dir, "\n", sep = "")
