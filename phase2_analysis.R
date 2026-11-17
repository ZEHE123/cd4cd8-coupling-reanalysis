# =============================================================================
# Phase 2 re-analysis: Oldham's method and the variance-ratio
# (Pitman-Morgan) test applied to published CD4/CD8 ratio summary data
#
# Project : CD4/CD8 ratio and mathematical coupling (HDAT9901, UNSW)
# Student : z5519901
# Input   : data/phase2_data.csv  (one row per study group)
# Output  : results/table4.csv             Table 4 (data used), one row per group
#           results/phase2_results.csv      main table, one row per group
#           results/sensitivity_r.csv       every group across r = 0.3 to 0.9
#           results/gidari_rounding.csv     rounding analysis for Gidari 2023
#           results/taramasso_feasible.csv  feasible region for Taramasso 2022
#           results/monsalvo_rounding.csv   rounding analysis for Monsalvo 2019
#           results/simulation_check.csv    formula check on simulated data
#           results/session_info.txt        R and package versions
# Run     : Rscript phase2_analysis.R   (from the repository root), or open
#           cd4cd8-coupling-reanalysis.Rproj in RStudio and source this file.
# Needs   : base R only (no packages).
# =============================================================================


## ---- 0. Setup ---------------------------------------------------------------

data_file <- file.path("data", "phase2_data.csv")
if (!file.exists(data_file)) {
  stop("data/phase2_data.csv not found. Open cd4cd8-coupling-reanalysis.Rproj, ",
       "or setwd() to the repository root, then run again.")
}
dir.create("results", showWarnings = FALSE)

r_grid <- seq(0.3, 0.9, by = 0.1)  # sensitivity range for the unreported r
alpha  <- 0.05


## ---- 1. Functions -----------------------------------------------------------
# Notation: s1 = baseline SD, s2 = follow-up SD, sd_change = SD of the change,
# r = correlation between baseline and follow-up, n = number of pairs.

# SD from the interquartile range: Wan et al. (2014), median + IQR + n only
# (their scenario C3). Like IQR/1.35, this assumes approximate normality;
# it differs from IQR/1.35 only by a small-sample correction.
sd_from_iqr <- function(q1, q3, n) {
  (q3 - q1) / (2 * qnorm((0.75 * n - 0.125) / (n + 0.25)))
}

# Baseline-follow-up correlation recovered from the SD of the change:
# var(change) = s1^2 + s2^2 - 2 r s1 s2
r_from_sd_change <- function(s1, s2, sd_change) {
  (s1^2 + s2^2 - sd_change^2) / (2 * s1 * s2)
}

# Oldham's method: correlation between the change (y2 - y1) and the average
# ((y1 + y2) / 2), written in summary statistics. Under the null hypothesis of
# no true relation between change and initial level it is expected to be 0.
# Positive = people starting higher gained more; negative = people starting
# lower gained more.
oldham_r <- function(s1, s2, r) {
  (s2^2 - s1^2) / sqrt((s1^2 + s2^2)^2 - 4 * r^2 * s1^2 * s2^2)
}

# Variance-ratio test (Pitman 1939; Morgan 1939; Tu & Gilthorpe 2007).
# H0: var(baseline) = var(follow-up). t on n - 2 df.
# This t is identical to the t-test of Oldham's correlation (checked in
# section 2), so the two methods give one test, seen from two sides.
pitman_morgan <- function(s1, s2, r, n) {
  vr <- s2^2 / s1^2
  t  <- (vr - 1) * sqrt(n - 2) / (2 * sqrt(vr * (1 - r^2)))
  list(vr = vr, t = t, df = n - 2, p = 2 * pt(-abs(t), df = n - 2))
}

# Smallest r at which p reaches alpha (p falls as r rises when VR != 1).
# NA if p is below alpha for every r, or above it for every r.
r_at_alpha <- function(s1, s2, n, alpha = 0.05) {
  f <- function(r) pitman_morgan(s1, s2, r, n)$p - alpha
  lo <- 0; hi <- 0.9999
  if (f(lo) * f(hi) > 0) return(NA_real_)
  uniroot(f, c(lo, hi), tol = 1e-8)$root
}

# Follow-up SDs compatible with a known baseline SD and SD of change, for a
# given r (used when the follow-up SD is not reported). Real roots need
# r >= sqrt(1 - sd_change^2 / s1^2).
feasible_s2 <- function(s1, sd_change, r) {
  disc <- sd_change^2 - s1^2 * (1 - r^2)
  if (disc < 0) return(numeric(0))
  s2 <- unique(r * s1 + c(-1, 1) * sqrt(disc))
  s2[s2 > 0]
}


## ---- 2. Check the formulas on simulated individual-level data ---------------
# Summary-statistic formulas must reproduce what the raw-data calculation
# gives. If any line fails, the script stops here.

set.seed(20261102)
n_sim <- 500
z1 <- rnorm(n_sim); z2 <- rnorm(n_sim)
y1 <- 0.9 + 0.30 * z1
y2 <- 1.0 + 0.45 * (0.6 * z1 + sqrt(1 - 0.6^2) * z2)
s1_sim <- sd(y1); s2_sim <- sd(y2); r_sim <- cor(y1, y2)
avg <- (y1 + y2) / 2; chg <- y2 - y1

sim_check <- data.frame(
  quantity = c("Oldham r", "t statistic", "r between baseline and follow-up"),
  from_individual_data = c(cor(avg, chg),
                           unname(cor.test(avg, chg)$statistic),
                           r_sim),
  from_summary_formula = c(oldham_r(s1_sim, s2_sim, r_sim),
                           pitman_morgan(s1_sim, s2_sim, r_sim, n_sim)$t,
                           r_from_sd_change(s1_sim, s2_sim, sd(chg))))
stopifnot(isTRUE(all.equal(sim_check$from_individual_data,
                           sim_check$from_summary_formula)))
cat("Section 2: formula checks on simulated data passed.\n\n")
rm(z1, z2, y1, y2, avg, chg)


## ---- 3. Read the data ----------------------------------------------------------

# encoding = "UTF-8" (not fileEncoding) keeps the "ü" in Ergün readable
# whatever the computer's locale is.
d <- read.csv(data_file, encoding = "UTF-8", stringsAsFactors = FALSE)
stopifnot(nrow(d) == 7, !anyNA(d$n_baseline))

# Median (IQR) rows: convert all three IQRs with the same method.
is_iqr <- d$value_type == "median_iqr"
d$sd_method <- ifelse(is_iqr, "Wan et al. (2014) from IQR", "Reported")
d$bl_sd[is_iqr]  <- sd_from_iqr(d$bl_q1[is_iqr],  d$bl_q3[is_iqr],  d$n_baseline[is_iqr])
d$fu_sd[is_iqr]  <- sd_from_iqr(d$fu_q1[is_iqr],  d$fu_q3[is_iqr],  d$n_baseline[is_iqr])
d$chg_sd[is_iqr] <- sd_from_iqr(d$chg_q1[is_iqr], d$chg_q3[is_iqr], d$n_change[is_iqr])

# Recover r wherever baseline, follow-up and change SDs are all available.
d$r_recovered <- with(d, ifelse(!is.na(bl_sd) & !is.na(fu_sd) & !is.na(chg_sd),
                                r_from_sd_change(bl_sd, fu_sd, chg_sd), NA_real_))
stopifnot(all(is.na(d$r_recovered) | abs(d$r_recovered) <= 1))


## ---- 3b. Table 4: the summary data used ----------------------------------------
# Built from the data file as text, so values keep the decimals printed in the
# papers (e.g. "1.0", not 1).

txt <- read.csv(data_file, encoding = "UTF-8", colClasses = "character")
pm_sign <- " \u00b1 "
show_value <- function(type, m, s, med, q1, q3, n_txt = "", time = "") {
  v <- ifelse(type == "median_iqr",
              ifelse(med == "", "", sprintf("median %s (IQR %s to %s)", med, q1, q3)),
              ifelse(m == "", "", paste0(m, pm_sign, s)))
  v <- ifelse(v != "" & n_txt != "", paste0(v, " (n = ", n_txt, ")"), v)
  v <- ifelse(v != "" & time != "", paste0(v, ", ", time), v)
  ifelse(v == "", "-", v)
}
n_differs <- txt$n_change != "" & txt$n_change != txt$n_baseline
table4 <- data.frame(
  study = txt$study, group = txt$group, n = txt$n_baseline,
  baseline = show_value(txt$value_type, txt$bl_mean, txt$bl_sd, txt$bl_median,
                        txt$bl_q1, txt$bl_q3, ifelse(n_differs, txt$n_baseline, "")),
  follow_up = show_value(txt$value_type, txt$fu_mean, txt$fu_sd, txt$fu_median,
                         txt$fu_q1, txt$fu_q3, time = txt$followup),
  change = show_value(txt$value_type, txt$chg_mean, txt$chg_sd, txt$chg_median,
                      txt$chg_q1, txt$chg_q3, ifelse(n_differs, txt$n_change, "")),
  source = txt$source)
table4$follow_up[table4$follow_up == "-"] <- "Not reported"
cat("Table 4: summary data used (CD4/CD8 ratio; mean \u00b1 SD unless stated)\n")
print(table4, row.names = FALSE, right = FALSE)
cat("\n")
rm(txt)


## ---- 4. Analyse each study group ---------------------------------------------

# Sign of Oldham's r only; read it together with p (a sign with p >= 0.05
# is not evidence of anything).
oldham_sign <- function(o) {
  ifelse(abs(o) < 1e-12, "zero",
         ifelse(o > 0, "positive (higher start, larger gain)",
                       "negative (lower start, larger gain)"))
}

results <- list()
sens    <- list()

for (i in seq_len(nrow(d))) {
  g <- d[i, ]

  # (a) Follow-up SD not reported: the test is not identifiable.
  if (is.na(g$fu_sd)) {
    r_min <- sqrt(1 - g$chg_sd^2 / g$bl_sd^2)
    results[[i]] <- data.frame(
      id = g$id, study = g$study, group = g$group,
      coupled_claim = g$coupled_claim, n = g$n_baseline,
      bl_sd = g$bl_sd, fu_sd = NA, sd_method = g$sd_method,
      variance_ratio = NA, r_basis = "not identifiable",
      r_low = r_min, r_high = NA,
      oldham_at_r_low = NA, oldham_at_r_high = NA,
      p_at_r_low = NA, p_at_r_high = NA, r_at_p05 = NA,
      oldham_sign = NA,
      result = sprintf("Not testable: follow-up SD not reported (data only bound r >= %.2f)", r_min))
    next
  }

  n  <- g$n_baseline
  s1 <- g$bl_sd
  s2 <- g$fu_sd

  # Sensitivity grid for every group with both SDs.
  pm_grid <- lapply(r_grid, function(r) pitman_morgan(s1, s2, r, n))
  sens[[length(sens) + 1]] <- data.frame(
    id = g$id, study = g$study, group = g$group, r = r_grid, r_basis = "assumed",
    variance_ratio = sapply(pm_grid, `[[`, "vr"),
    oldham_r = oldham_r(s1, s2, r_grid),
    t = sapply(pm_grid, `[[`, "t"), df = n - 2,
    p = sapply(pm_grid, `[[`, "p"))

  # Main result: the recovered r if available, otherwise the sensitivity range.
  if (!is.na(g$r_recovered)) {
    r_use <- g$r_recovered
    r_basis <- "recovered from SD of change"
    pm <- pitman_morgan(s1, s2, r_use, n)
    sens[[length(sens) + 1]] <- data.frame(
      id = g$id, study = g$study, group = g$group, r = r_use, r_basis = "recovered",
      variance_ratio = pm$vr, oldham_r = oldham_r(s1, s2, r_use),
      t = pm$t, df = pm$df, p = pm$p)
    r_lo <- r_hi <- r_use
    p_lo <- p_hi <- pm$p
    result <- if (pm$p < alpha) "p < 0.05 at the recovered r" else "p >= 0.05 at the recovered r"
  } else {
    r_basis <- "sensitivity range"
    r_lo <- min(r_grid); r_hi <- max(r_grid)
    p_lo <- pitman_morgan(s1, s2, r_lo, n)$p
    p_hi <- pitman_morgan(s1, s2, r_hi, n)$p
    p_all <- sapply(pm_grid, `[[`, "p")
    result <- if (all(p_all < alpha)) {
      sprintf("p < 0.05 at every r from %.1f to %.1f", r_lo, r_hi)
    } else if (all(p_all >= alpha)) {
      sprintf("p >= 0.05 at every r from %.1f to %.1f", r_lo, r_hi)
    } else {
      sprintf("Depends on r: p < 0.05 only if r >= %.2f", r_at_alpha(s1, s2, n, alpha))
    }
  }

  results[[i]] <- data.frame(
    id = g$id, study = g$study, group = g$group,
    coupled_claim = g$coupled_claim, n = n,
    bl_sd = s1, fu_sd = s2, sd_method = g$sd_method,
    variance_ratio = s2^2 / s1^2, r_basis = r_basis,
    r_low = r_lo, r_high = r_hi,
    oldham_at_r_low = oldham_r(s1, s2, r_lo),
    oldham_at_r_high = oldham_r(s1, s2, r_hi),
    p_at_r_low = p_lo, p_at_r_high = p_hi,
    r_at_p05 = r_at_alpha(s1, s2, n, alpha),
    oldham_sign = oldham_sign(oldham_r(s1, s2, r_lo)),
    result = result)
}

results <- do.call(rbind, results)
sens    <- do.call(rbind, sens)


## ---- 5. Taramasso 2022: feasible region when the follow-up SD is missing -----
# Two unknowns (follow-up SD and r), one equation. For each r above the
# lower bound, list the follow-up SDs that fit, and Oldham's r for each.
# Note: baseline SD is from n = 120 and SD of change from n = 91.

tar <- d[d$id == "R0612", ]
r_min_tar <- sqrt(1 - tar$chg_sd^2 / tar$bl_sd^2)
r_seq <- seq(ceiling(r_min_tar * 1000) / 1000, 0.999, by = 0.001)
taramasso <- do.call(rbind, lapply(r_seq, function(r) {
  s2 <- feasible_s2(tar$bl_sd, tar$chg_sd, r)
  if (length(s2) == 0) return(NULL)
  data.frame(r = r, fu_sd = s2, variance_ratio = s2^2 / tar$bl_sd^2,
             oldham_r = oldham_r(tar$bl_sd, s2, r))
}))


## ---- 6. Gidari 2023: rounding analysis -----------------------------------------
# SDs are reported to one decimal place (0.8 and 0.8), so each true SD lies
# in [0.75, 0.85]. Test the corners and the reported values across r.

gid <- d[d$id == "R0646", ]
sd_vals <- c(0.75, 0.80, 0.85)
gidari <- expand.grid(bl_sd = sd_vals, fu_sd = sd_vals, r = r_grid)
gidari_pm <- mapply(function(a, b, r) unlist(pitman_morgan(a, b, r, gid$n_baseline)),
                    gidari$bl_sd, gidari$fu_sd, gidari$r)
gidari$variance_ratio <- gidari_pm["vr", ]
gidari$oldham_r <- oldham_r(gidari$bl_sd, gidari$fu_sd, gidari$r)
gidari$t <- gidari_pm["t", ]
gidari$p <- gidari_pm["p", ]
gidari$n <- gid$n_baseline


## ---- 7. Monsalvo 2019: does the IQR-to-SD method matter? ------------------------
# The same divisor is applied to the baseline, follow-up and change IQRs, so it
# cancels in the variance ratio, the recovered r, Oldham's r and the test.
# Re-run with IQR/1.35 to show this.

mon <- d[d$id == "R0406", ]
s1_alt <- (mon$bl_q3 - mon$bl_q1) / 1.35
s2_alt <- (mon$fu_q3 - mon$fu_q1) / 1.35
sd_alt <- (mon$chg_q3 - mon$chg_q1) / 1.35
r_alt  <- r_from_sd_change(s1_alt, s2_alt, sd_alt)
p_alt  <- pitman_morgan(s1_alt, s2_alt, r_alt, mon$n_baseline)$p
p_wan  <- results$p_at_r_low[results$id == "R0406"]
stopifnot(isTRUE(all.equal(r_alt, mon$r_recovered)), isTRUE(all.equal(p_alt, p_wan)))

# Correlation between baseline and change implied by the summary SDs:
# corr(y1, y2 - y1) = (r * s2 - s1) / sd_change. This is the coupled
# baseline-change correlation, comparable with the published Spearman rho.
mon_bc_corr <- with(mon, (r_recovered * fu_sd - bl_sd) / chg_sd)

# Rounding: quartiles are reported to two decimals, so each could be off by
# up to 0.005 and each IQR width by up to 0.01. Vary the three widths over
# five values each (step 0.005): 5 x 5 x 5 = 125 combinations.
w_bl  <- (mon$bl_q3  - mon$bl_q1)  + seq(-0.01, 0.01, by = 0.005)
w_fu  <- (mon$fu_q3  - mon$fu_q1)  + seq(-0.01, 0.01, by = 0.005)
w_chg <- (mon$chg_q3 - mon$chg_q1) + seq(-0.01, 0.01, by = 0.005)
eta   <- 2 * qnorm((0.75 * mon$n_baseline - 0.125) / (mon$n_baseline + 0.25))
monsalvo <- expand.grid(iqr_bl = w_bl, iqr_fu = w_fu, iqr_chg = w_chg)
monsalvo$bl_sd  <- monsalvo$iqr_bl  / eta
monsalvo$fu_sd  <- monsalvo$iqr_fu  / eta
monsalvo$chg_sd <- monsalvo$iqr_chg / eta
monsalvo$r <- r_from_sd_change(monsalvo$bl_sd, monsalvo$fu_sd, monsalvo$chg_sd)
monsalvo$variance_ratio <- monsalvo$fu_sd^2 / monsalvo$bl_sd^2
monsalvo$oldham_r <- oldham_r(monsalvo$bl_sd, monsalvo$fu_sd, monsalvo$r)
monsalvo$p <- mapply(function(a, b, r) pitman_morgan(a, b, r, mon$n_baseline)$p,
                     monsalvo$bl_sd, monsalvo$fu_sd, monsalvo$r)
stopifnot(nrow(monsalvo) == 125, all(abs(monsalvo$r) < 1))
mon_sig_neg <- sum(monsalvo$oldham_r < 0 & monsalvo$p < alpha)
mon_sig_pos <- sum(monsalvo$oldham_r > 0 & monsalvo$p < alpha)


## ---- 8. Save -------------------------------------------------------------------

round_df <- function(x, digits = 4) {
  num <- vapply(x, is.numeric, logical(1))
  x[num] <- lapply(x[num], function(v) signif(v, digits))
  x
}
write.csv(table4,              "results/table4.csv",            row.names = FALSE)
write.csv(round_df(results),   "results/phase2_results.csv",     row.names = FALSE)
write.csv(round_df(sens),      "results/sensitivity_r.csv",      row.names = FALSE)
write.csv(round_df(gidari),    "results/gidari_rounding.csv",    row.names = FALSE)
write.csv(round_df(taramasso), "results/taramasso_feasible.csv", row.names = FALSE)
write.csv(round_df(monsalvo),  "results/monsalvo_rounding.csv",  row.names = FALSE)
write.csv(round_df(sim_check, 8), "results/simulation_check.csv", row.names = FALSE)
writeLines(c(paste("Run on", format(Sys.time(), "%Y-%m-%d %H:%M %Z")), "",
             capture.output(sessionInfo())),
           "results/session_info.txt")


## ---- 9. Console summary ----------------------------------------------------------

f2 <- function(x) ifelse(is.na(x), "-", formatC(x, format = "f", digits = 2))
fp <- function(p) ifelse(is.na(p), "-", ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3)))

cat("Phase 2 results (variance ratio = follow-up variance / baseline variance)\n")
cat(strrep("-", 78), "\n")
for (k in seq_len(nrow(results))) {
  x <- results[k, ]
  cat(sprintf("%s %s (%s), n = %d\n", x$id, x$study, x$group, x$n))
  if (is.na(x$variance_ratio)) {
    cat("   ", x$result, "\n")
    next
  }
  r_txt <- if (x$r_low == x$r_high) sprintf("r = %.3f (%s)", x$r_low, x$r_basis)
           else sprintf("r = %.1f to %.1f", x$r_low, x$r_high)
  cat(sprintf("   VR = %s | %s | Oldham r = %s to %s | p = %s to %s\n",
              f2(x$variance_ratio), r_txt, f2(x$oldham_at_r_low), f2(x$oldham_at_r_high),
              fp(x$p_at_r_low), fp(x$p_at_r_high)))
  cat("   ", x$result, "| sign of Oldham r:", x$oldham_sign, "\n")
}
cat(strrep("-", 78), "\n")

cat(sprintf("Taramasso 2022: r >= %.3f; across the feasible region Oldham r runs from %.2f to %.2f\n",
            r_min_tar, min(taramasso$oldham_r), max(taramasso$oldham_r)))
corner <- gidari[gidari$variance_ratio == max(gidari$variance_ratio), ]
cat(sprintf("Gidari 2023 rounding: VR %.2f to %.2f; at the extremes p = %s (r = 0.3) and %s (r = 0.9)\n",
            min(gidari$variance_ratio), max(gidari$variance_ratio),
            fp(corner$p[corner$r == 0.3]), fp(corner$p[abs(corner$r - 0.9) < 1e-9])))
cat(sprintf("Monsalvo 2019: Wan (2014) and IQR/1.35 give the same r (%.3f) and p (%.3f)\n",
            mon$r_recovered, p_wan))
cat(sprintf("Monsalvo 2019: implied baseline-change correlation %.3f (published Spearman rho -0.25)\n",
            mon_bc_corr))
cat(sprintf("Monsalvo 2019 rounding (125 combinations): r %.3f to %.3f; Oldham r %.2f to %.2f; p %s to %s\n",
            min(monsalvo$r), max(monsalvo$r), min(monsalvo$oldham_r), max(monsalvo$oldham_r),
            fp(min(monsalvo$p)), fp(max(monsalvo$p))))
cat(sprintf("   significant negative (the published direction): %d; significant positive: %d\n",
            mon_sig_neg, mon_sig_pos))


## ---- 10. Check against the values reported in the dissertation -----------------
# Values from Table 5 and Section 6.3 (rounded as printed there).
# Prints OK or CHECK. If a value is CHECK, find out which side is wrong
# before changing either.

get <- function(id, grp, col) results[results$id == id & results$group == grp, col]
report <- data.frame(
  item = c("Monsalvo s1", "Monsalvo s2", "Monsalvo s_d", "Monsalvo r recovered",
           "Monsalvo VR", "Monsalvo Oldham r", "Monsalvo p",
           "Trizzino naive VR", "Trizzino naive Oldham (r=0.3)", "Trizzino naive Oldham (r=0.9)",
           "Trizzino naive p (r=0.3)",
           "Trizzino exp. VR", "Trizzino exp. Oldham (r=0.3)", "Trizzino exp. Oldham (r=0.9)",
           "Gidari VR",
           "Ergun non-INSTI VR", "Ergun non-INSTI Oldham (r=0.3)", "Ergun non-INSTI Oldham (r=0.9)",
           "Ergun non-INSTI p (r=0.3)", "Ergun non-INSTI p (r=0.9)",
           "Ergun INSTI VR", "Ergun INSTI Oldham (r=0.3)", "Ergun INSTI Oldham (r=0.9)",
           "Ergun non-INSTI r where p = 0.05",
           "Taramasso r lower bound",
           "Gidari rounding VR (lowest)", "Gidari rounding VR (highest)",
           "Gidari rounding p (r=0.3)", "Gidari rounding p (r=0.9)",
           "Monsalvo baseline-change corr.",
           "Monsalvo rounding: significant negative"),
  reported = c(0.380, 0.388, 0.119, 0.952,
               1.04, 0.06, 0.32,
               3.45, 0.57, 0.83, 0.002,
               1.80, 0.30, 0.57,
               1.00,
               1.35, 0.16, 0.33, 0.17, 0.003,
               2.45, 0.44, 0.73,
               0.75,
               0.91,
               0.78, 1.28,
               0.25, 0.01,
               -0.09,
               0),
  digits = c(3, 3, 3, 3,
             2, 2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 2,
             2, 2, 2, 2, 2, 2, 2, 0),
  script = c(mon$bl_sd, mon$fu_sd, mon$chg_sd, mon$r_recovered,
             get("R0406", "All dual regimens", "variance_ratio"),
             get("R0406", "All dual regimens", "oldham_at_r_low"),
             get("R0406", "All dual regimens", "p_at_r_low"),
             get("R0818", "ART-naive", "variance_ratio"),
             get("R0818", "ART-naive", "oldham_at_r_low"),
             get("R0818", "ART-naive", "oldham_at_r_high"),
             get("R0818", "ART-naive", "p_at_r_low"),
             get("R0818", "ART-experienced", "variance_ratio"),
             get("R0818", "ART-experienced", "oldham_at_r_low"),
             get("R0818", "ART-experienced", "oldham_at_r_high"),
             get("R0646", "ART-experienced aged >=60", "variance_ratio"),
             get("R0780", "Non-INSTI", "variance_ratio"),
             get("R0780", "Non-INSTI", "oldham_at_r_low"),
             get("R0780", "Non-INSTI", "oldham_at_r_high"),
             get("R0780", "Non-INSTI", "p_at_r_low"),
             get("R0780", "Non-INSTI", "p_at_r_high"),
             get("R0780", "INSTI", "variance_ratio"),
             get("R0780", "INSTI", "oldham_at_r_low"),
             get("R0780", "INSTI", "oldham_at_r_high"),
             get("R0780", "Non-INSTI", "r_at_p05"),
             r_min_tar,
             min(gidari$variance_ratio), max(gidari$variance_ratio),
             corner$p[corner$r == 0.3],
             corner$p[abs(corner$r - 0.9) < 1e-9],
             mon_bc_corr,
             mon_sig_neg))
stopifnot(length(report$script) == nrow(report), length(report$digits) == nrow(report))
report$match <- ifelse(round(report$script, report$digits) == round(report$reported, report$digits),
                       "OK", "CHECK")
cat("\nCheck against Table 5 and Section 6.3:\n")
print(report[, c("item", "reported", "script", "match")], row.names = FALSE, digits = 4)
cat(sprintf("\n%d of %d values match the dissertation.\n", sum(report$match == "OK"), nrow(report)))
