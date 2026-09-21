#!/usr/bin/env Rscript

# ---------------------------------------------------------------------------
# Genome-wide pi, dxy and Fst from windowed pixy output, with bootstrap CIs
#
# Takes the three per-window tables written by rule pixy_pi_dxy_fst
# (workflow/rules/3_diversity_stats.smk): pixy_pi.txt, pixy_dxy.txt and
# pixy_fst.txt, and collapses each to one genome-wide estimate per population
# (pi) or population pair (dxy, Fst), with a 95% confidence interval from
# bootstrap resampling of windows with replacement.
#
# Point estimates:
#   pi, dxy -> sum(count_diffs) / sum(count_comparisons) across windows, i.e.
#              the same ratio pixy itself uses per window (pixy/calc.py),
#              just pooled across windows before dividing. This is the
#              unbiased genome-wide estimator; averaging the per-window
#              avg_pi/avg_dxy values instead would let windows with a small
#              denominator (little usable data) swing the estimate.
#   Fst     -> a no_snps-weighted mean of the per-window avg_*_fst column.
#              pixy's Fst output only carries the raw variance components
#              (Weir-Cockerham a/b/c, or Hudson's numerator/denominator) when
#              --fst_components is passed, which rule pixy_pi_dxy_fst does
#              not do, so an exact pooled estimate can't be recomputed here.
#              A SNP-count-weighted mean of the per-window estimates is the
#              standard approximation for genome-wide Fst from pixy output.
#
# Bootstrap: for each population/pair, its windows are resampled with
# replacement (same number of windows drawn each time), n_boot times; the
# point-estimate formula above is recomputed on every resample, and the 95%
# CI is the 2.5th/97.5th percentile of the resulting distribution.
#
# Windows with an undefined statistic (pixy writes "NA": no comparisons for
# pi/dxy, or no SNPs for Fst) are dropped before both the point estimate and
# the bootstrap for that population/pair.
#
# Usage:
#   Rscript pixy_genomewide_bootstrap.R <pixy_pi.txt> <pixy_dxy.txt> \
#       <pixy_fst.txt> [output_dir] [n_boot] [seed]
#
#   output_dir defaults to the directory of <pixy_pi.txt>
#   n_boot     defaults to 10000
#   seed       defaults to 42 (for reproducible CIs)
#
# Writes, into output_dir:
#   pixy_pi_genomewide.tsv   pop       n_windows genome_pi  ci_lower ci_upper n_boot
#   pixy_dxy_genomewide.tsv  pop1 pop2 n_windows genome_dxy ci_lower ci_upper n_boot
#   pixy_fst_genomewide.tsv  pop1 pop2 n_windows genome_fst ci_lower ci_upper n_boot
# ---------------------------------------------------------------------------

# ---- 1. Arguments -----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(
    "Usage: Rscript pixy_genomewide_bootstrap.R <pixy_pi.txt> <pixy_dxy.txt> ",
    "<pixy_fst.txt> [output_dir] [n_boot] [seed]"
  )
}
pi_file <- args[1]
dxy_file <- args[2]
fst_file <- args[3]
out_dir <- if (length(args) >= 4) args[4] else dirname(pi_file)
n_boot <- if (length(args) >= 5) as.integer(args[5]) else 10000L
seed <- if (length(args) >= 6) as.integer(args[6]) else 42L

for (f in c(pi_file, dxy_file, fst_file)) {
  if (!file.exists(f)) stop(sprintf("Input file does not exist: %s", f))
}
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

set.seed(seed)

# ---- 2. Bootstrap engine -----------------------------------------------------
# stat_fn(idx) computes the point statistic from a set of row indices into the
# windows for one population/pair. Shared by pi, dxy and Fst below so the
# window-resampling logic itself only lives in one place.
bootstrap_ci <- function(n_windows, stat_fn, n_boot, ci_level = 0.95) {
  if (n_windows == 0) {
    return(c(estimate = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_))
  }
  point_est <- stat_fn(seq_len(n_windows))
  boot_est <- vapply(seq_len(n_boot), function(b) {
    stat_fn(sample.int(n_windows, n_windows, replace = TRUE))
  }, numeric(1))
  alpha <- (1 - ci_level) / 2
  bounds <- quantile(boot_est, probs = c(alpha, 1 - alpha), na.rm = TRUE, names = FALSE)
  c(estimate = point_est, ci_lower = bounds[1], ci_upper = bounds[2])
}

# Genome-wide pi/dxy statistic: ratio of summed counts (see header comment).
ratio_of_sums <- function(numerator, denominator) {
  function(idx) sum(numerator[idx]) / sum(denominator[idx])
}

# Genome-wide Fst statistic: SNP-count-weighted mean (see header comment).
weighted_mean_fn <- function(values, weights) {
  function(idx) weighted.mean(values[idx], weights[idx])
}

read_pixy <- function(path) {
  read.delim(path, header = TRUE, stringsAsFactors = FALSE, na.strings = "NA")
}

# do.call(rbind, list()) returns NULL, which write.table() silently turns into
# a blank line with no header. If every population/pair ends up with zero
# usable windows (e.g. an all-NA input), fall back to an empty frame with the
# right column names so the output file still has a proper header.
bind_rows_or_empty <- function(rows, col_names) {
  if (length(rows) == 0) {
    return(as.data.frame(setNames(
      replicate(length(col_names), character(0), simplify = FALSE), col_names
    )))
  }
  do.call(rbind, rows)
}

require_cols <- function(df, cols, path) {
  missing_cols <- setdiff(cols, names(df))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "%s is missing expected column(s): %s",
      path, paste(missing_cols, collapse = ", ")
    ))
  }
}

# ---- 3. Genome-wide pi (per population) --------------------------------------
summarise_pi <- function(path, n_boot) {
  df <- read_pixy(path)
  require_cols(df, c("pop", "count_diffs", "count_comparisons"), path)
  df <- df[!is.na(df$count_diffs) & !is.na(df$count_comparisons) & df$count_comparisons > 0, ]

  pops <- sort(unique(df$pop))
  rows <- lapply(pops, function(p) {
    sub <- df[df$pop == p, ]
    message(sprintf("  pi: %s (%d windows)", p, nrow(sub)))
    res <- bootstrap_ci(nrow(sub), ratio_of_sums(sub$count_diffs, sub$count_comparisons), n_boot)
    data.frame(
      pop = p,
      n_windows = nrow(sub),
      genome_pi = res[["estimate"]],
      ci_lower = res[["ci_lower"]],
      ci_upper = res[["ci_upper"]],
      n_boot = n_boot
    )
  })
  bind_rows_or_empty(rows, c("pop", "n_windows", "genome_pi", "ci_lower", "ci_upper", "n_boot"))
}

# ---- 4. Genome-wide dxy (per population pair) ---------------------------------
summarise_dxy <- function(path, n_boot) {
  df <- read_pixy(path)
  require_cols(df, c("pop1", "pop2", "count_diffs", "count_comparisons"), path)
  df <- df[!is.na(df$count_diffs) & !is.na(df$count_comparisons) & df$count_comparisons > 0, ]

  pairs <- unique(df[, c("pop1", "pop2")])
  rows <- lapply(seq_len(nrow(pairs)), function(i) {
    p1 <- pairs$pop1[i]
    p2 <- pairs$pop2[i]
    sub <- df[df$pop1 == p1 & df$pop2 == p2, ]
    message(sprintf("  dxy: %s x %s (%d windows)", p1, p2, nrow(sub)))
    res <- bootstrap_ci(nrow(sub), ratio_of_sums(sub$count_diffs, sub$count_comparisons), n_boot)
    data.frame(
      pop1 = p1,
      pop2 = p2,
      n_windows = nrow(sub),
      genome_dxy = res[["estimate"]],
      ci_lower = res[["ci_lower"]],
      ci_upper = res[["ci_upper"]],
      n_boot = n_boot
    )
  })
  bind_rows_or_empty(rows, c("pop1", "pop2", "n_windows", "genome_dxy", "ci_lower", "ci_upper", "n_boot"))
}

# ---- 5. Genome-wide Fst (per population pair) ----------------------------------
summarise_fst <- function(path, n_boot) {
  df <- read_pixy(path)
  require_cols(df, c("pop1", "pop2", "no_snps"), path)

  # The averaged column is named after --fst_type ("avg_wc_fst" or
  # "avg_hudson_fst"), so it is picked out by pattern instead of a fixed name.
  fst_col <- grep("^avg_.*_fst$", names(df), value = TRUE)
  if (length(fst_col) != 1) {
    stop(sprintf(
      "Expected exactly one 'avg_<type>_fst' column in %s, found: %s",
      path, paste(fst_col, collapse = ", ")
    ))
  }
  df$fst_value <- df[[fst_col]]
  df <- df[!is.na(df$fst_value) & !is.na(df$no_snps) & df$no_snps > 0, ]

  pairs <- unique(df[, c("pop1", "pop2")])
  rows <- lapply(seq_len(nrow(pairs)), function(i) {
    p1 <- pairs$pop1[i]
    p2 <- pairs$pop2[i]
    sub <- df[df$pop1 == p1 & df$pop2 == p2, ]
    message(sprintf("  fst: %s x %s (%d windows)", p1, p2, nrow(sub)))
    res <- bootstrap_ci(nrow(sub), weighted_mean_fn(sub$fst_value, sub$no_snps), n_boot)
    data.frame(
      pop1 = p1,
      pop2 = p2,
      n_windows = nrow(sub),
      genome_fst = res[["estimate"]],
      ci_lower = res[["ci_lower"]],
      ci_upper = res[["ci_upper"]],
      n_boot = n_boot
    )
  })
  bind_rows_or_empty(rows, c("pop1", "pop2", "n_windows", "genome_fst", "ci_lower", "ci_upper", "n_boot"))
}

# ---- 6. Run and write --------------------------------------------------------
message(sprintf("Bootstrapping with n_boot = %d, seed = %d", n_boot, seed))

message("Computing genome-wide pi...")
pi_out <- summarise_pi(pi_file, n_boot)
message("Computing genome-wide dxy...")
dxy_out <- summarise_dxy(dxy_file, n_boot)
message("Computing genome-wide Fst...")
fst_out <- summarise_fst(fst_file, n_boot)

write.table(
  pi_out, file.path(out_dir, "pixy_pi_genomewide.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  dxy_out, file.path(out_dir, "pixy_dxy_genomewide.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  fst_out, file.path(out_dir, "pixy_fst_genomewide.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

message(sprintf(
  "Wrote pixy_{pi,dxy,fst}_genomewide.tsv to: %s",
  normalizePath(out_dir, mustWork = FALSE)
))
