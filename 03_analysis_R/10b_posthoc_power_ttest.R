# =============================================================================
# THESIS FIGURES PRODUCED BY THIS SCRIPT
#   (no main figure) post-hoc power analysis reported in section 2.10 (two-sample t-test, d = 1.2, Bonferroni over 588 regions)
#
# Note: figure numbers appearing in the comments further down refer to earlier
# drafts and are NOT the final numbering. The list above is authoritative.
# =============================================================================

# =============================================================================
# Post-hoc power analysis (two-sample t-test, Bonferroni-corrected)
#
# Reproduces the power analysis reported in section 2.10 of the thesis:
#   "a post-hoc power analysis (using a two-sample t-test, effect size d = 1.2,
#    alpha = 0.05, and Bonferroni correction for 588 regions, resulting in an
#    alpha of 8.5e-05) was conducted to determine the detection limit of this
#    study."
#
# Rationale: the region-wise comparisons are run across 588 brain regions
# (fibre tracts excluded). Bonferroni correction over that many tests gives
# alpha = 0.05 / 588 = 8.50e-05. This script quantifies, for the group sizes
# actually realised in the study, (a) the achieved power at the assumed large
# effect (d = 1.2), (b) the sample size that WOULD have been required for 80%
# power, and (c) the effect size that the largest cohort could have detected.
#
# Only base R (stats::power.t.test) is required.
#
# Outputs (written to OUT_DIR):
#   posthoc_power_ttest.csv     - power per group size
#   posthoc_power_ttest.png     - power curve with the realised group sizes
# =============================================================================

OUT_DIR <- path.expand("~/Desktop/Bachelor")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

N_REGIONS <- 588       # brain regions tested, fibre tracts excluded
ALPHA_RAW <- 0.05
DELTA     <- 1.2       # assumed effect size (Cohen's d), sd = 1
ALPHA     <- ALPHA_RAW / N_REGIONS

cat(sprintf("Bonferroni-corrected alpha = %.3g / %d = %.3e\n\n",
            ALPHA_RAW, N_REGIONS, ALPHA))

# ---------------------------------------------------------------------------
# (a) achieved power across the realised group sizes
#     n = 6-18  pooled anterograde conditions
#     n = 2-4   sex-specific and retrograde cohorts
# ---------------------------------------------------------------------------
grid <- data.frame(n_per_group = c(2, 3, 4, 6, 8, 10, 12, 14, 16, 18))
grid$cohort <- ifelse(grid$n_per_group <= 4,
                      "sex-specific / retrograde",
                      "pooled anterograde")

grid$power_bonferroni <- vapply(grid$n_per_group, function(n)
  power.t.test(n = n, delta = DELTA, sd = 1,
               sig.level = ALPHA, type = "two.sample")$power,
  numeric(1))

grid$power_uncorrected <- vapply(grid$n_per_group, function(n)
  power.t.test(n = n, delta = DELTA, sd = 1,
               sig.level = ALPHA_RAW, type = "two.sample")$power,
  numeric(1))

print(grid, row.names = FALSE, digits = 3)

# ---------------------------------------------------------------------------
# (b) what would have been required, (c) what the largest cohort could detect
# ---------------------------------------------------------------------------
n_needed <- power.t.test(delta = DELTA, sd = 1, sig.level = ALPHA,
                         power = 0.80, type = "two.sample")$n
d_detect <- power.t.test(n = 18, sd = 1, sig.level = ALPHA,
                         power = 0.80, type = "two.sample")$delta
p_max    <- grid$power_bonferroni[grid$n_per_group == 18]

cat(sprintf("\nLargest realised group (n = 18): power = %.1f%%\n", 100 * p_max))
cat(sprintf("n required for 80%% power at d = %.1f: %.1f per group\n", DELTA, n_needed))
cat(sprintf("Effect size detectable at 80%% power with n = 18: d = %.2f\n", d_detect))
cat(sprintf("For comparison, power at n = 18 without correction: %.1f%%\n",
            100 * grid$power_uncorrected[grid$n_per_group == 18]))

write.csv(grid, file.path(OUT_DIR, "posthoc_power_ttest.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# power curve
# ---------------------------------------------------------------------------
ns <- 2:40
pw <- vapply(ns, function(n)
  power.t.test(n = n, delta = DELTA, sd = 1,
               sig.level = ALPHA, type = "two.sample")$power,
  numeric(1))

png(file.path(OUT_DIR, "posthoc_power_ttest.png"),
    width = 1800, height = 1300, res = 300)
par(mar = c(4.5, 4.5, 3, 1))
plot(ns, pw, type = "l", lwd = 2, col = "#6B4C9A", ylim = c(0, 1),
     xlab = "n per group", ylab = "Power",
     main = sprintf("Post-hoc power, d = %.1f, alpha = %.2e", DELTA, ALPHA))
abline(h = 0.80, lty = 2, col = "grey50")
points(c(4, 18), pw[c(4, 18) - 1], pch = 19, col = "#D6547F")
text(18, pw[17], sprintf(" n = 18: %.1f%%", 100 * p_max), pos = 4, cex = 0.8)
text(4, pw[3], sprintf(" n = 4: %.1f%%", 100 * pw[3]), pos = 4, cex = 0.8)
abline(v = n_needed, lty = 3, col = "grey50")
text(n_needed, 0.05, sprintf(" n = %.0f for 80%%", n_needed), pos = 4, cex = 0.8)
dev.off()

cat(sprintf("\nWritten to %s\n", OUT_DIR))
