#!/usr/bin/env Rscript

# ---------------------------------------------------------------------------
# Windowed Fst, dxy and pi along the genome, modern populations only
#
# Takes the per-window tables written by rule pixy_pi_dxy_fst
# (workflow/rules/3_diversity_stats.smk), pixy_fst.txt, pixy_dxy.txt and
# pixy_pi.txt, and plots the statistic (y) against the genomic windows (x,
# chromosomes laid end to end). Each statistic goes to its own single-page PDF
# with one stacked panel per pair of modern populations (Fst, dxy) or per
# modern population (pi), on a shared x- and y-axis so they can be compared
# directly. Each statistic has its own y-axis range. For pi the axis is cut
# at the 99.9th percentile of the plotted windows so the bulk of the data uses
# the panel height; the few windows above the cut are drawn as triangles at
# the top edge of their panel instead of at their true value.
#
# Windows with an undefined statistic (pixy writes "NA") are not drawn. dxy
# and pi windows with fewer than min_sites genotyped sites (no_sites) are
# dropped too: the few windows with next to no data give extreme values that
# otherwise set the y-axis. Fst windows with fewer than min_snps SNPs
# (no_snps; pixy's Fst table has no no_sites column) are dropped for the same
# reason: Fst from a handful of SNPs is mostly noise.
#
# Usage:
#   Rscript plot_pixy_modern_windows.R <pixy_dir> [out_prefix] [modern_pops] \
#       [min_sites] [min_snps]
#
#   pixy_dir     folder holding pixy_fst.txt, pixy_dxy.txt and pixy_pi.txt
#   out_prefix   defaults to <pixy_dir>/pixy_modern_windows
#   modern_pops  comma-separated population codes; panels are stacked in this
#                order. Defaults to Spa,Kaz,Swe,Nth,Dev
#   min_sites    minimum no_sites for a dxy or pi window to be plotted;
#                defaults to 100, 0 keeps every window
#   min_snps     minimum no_snps for an Fst window to be plotted; defaults to
#                10, 0 keeps every window
#
# Writes <out_prefix>_fst.pdf, <out_prefix>_dxy.pdf and <out_prefix>_pi.pdf
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(ggplot2))

# ---- 1. Arguments -----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop(
    "Usage: Rscript plot_pixy_modern_windows.R <pixy_dir> [out_prefix] ",
    "[modern_pops] [min_sites] [min_snps]"
  )
}
pixy_dir <- args[1]
out_prefix <- if (length(args) >= 2) {
  args[2]
} else {
  file.path(pixy_dir, "pixy_modern_windows")
}
modern_pops <- if (length(args) >= 3) {
  strsplit(args[3], ",")[[1]]
} else {
  c("Spa", "Kaz", "Swe", "Nth", "Dev")
}
min_sites <- if (length(args) >= 4) as.numeric(args[4]) else 100
min_snps <- if (length(args) >= 5) as.numeric(args[5]) else 10

fst_file <- file.path(pixy_dir, "pixy_fst.txt")
dxy_file <- file.path(pixy_dir, "pixy_dxy.txt")
pi_file <- file.path(pixy_dir, "pixy_pi.txt")
for (f in c(fst_file, dxy_file, pi_file)) {
  if (!file.exists(f)) stop(sprintf("Input file does not exist: %s", f))
}
if (length(modern_pops) < 2) stop("Need at least two modern populations")

# ---- 2. Read windows and place them on a genome-wide x-axis -----------------
# value_pattern picks the statistic column: avg_pi, avg_dxy, or avg_hudson_fst
# / avg_wc_fst depending on the fst_type pixy was run with. Works for both the
# per-pair tables (pop1, pop2) and the per-population pi table (pop); the
# returned `panel` column is the pair or population each window belongs to.
# Windows with fewer than min_count in count_col (no_sites or no_snps) are
# dropped.
read_windows <- function(file, value_pattern, count_col, min_count = 0) {
  d <- read.delim(file, stringsAsFactors = FALSE)
  value_col <- grep(value_pattern, names(d), value = TRUE)
  if (length(value_col) != 1) {
    stop(sprintf("Could not identify the statistic column in %s", file))
  }
  by_pair <- "pop1" %in% names(d)
  file_pops <- if (by_pair) c(d$pop1, d$pop2) else d$pop
  missing_pops <- setdiff(modern_pops, file_pops)
  if (length(missing_pops) > 0) {
    stop(sprintf(
      "Population(s) not found in %s: %s",
      file, paste(missing_pops, collapse = ", ")
    ))
  }
  d$chromosome <- as.character(d$chromosome)

  # chromosomes in numeric order (names are plain integers), anything
  # non-numeric after them; each one offset by the lengths of those before it
  chr_end <- tapply(d$window_pos_2, d$chromosome, max)
  chr_num <- suppressWarnings(as.numeric(names(chr_end)))
  chr_end <- chr_end[order(chr_num, names(chr_end))]
  chroms <- data.frame(
    chromosome = names(chr_end),
    offset = c(0, head(cumsum(as.numeric(chr_end)), -1)),
    stringsAsFactors = FALSE
  )
  chroms$end <- chroms$offset + as.numeric(chr_end)
  chroms$mid <- chroms$offset + as.numeric(chr_end) / 2
  chroms$shade <- factor(seq_len(nrow(chroms)) %% 2)

  if (by_pair) {
    d <- d[d$pop1 %in% modern_pops & d$pop2 %in% modern_pops, ]
    # pixy lists each pair once, in the poplist's order; label every pair in
    # modern_pops order instead
    i1 <- match(d$pop1, modern_pops)
    i2 <- match(d$pop2, modern_pops)
    pairs <- combn(modern_pops, 2)
    panel_levels <- sprintf("%s vs %s", pairs[1, ], pairs[2, ])
    d$panel <- sprintf(
      "%s vs %s", modern_pops[pmin(i1, i2)], modern_pops[pmax(i1, i2)]
    )
  } else {
    d <- d[d$pop %in% modern_pops, ]
    panel_levels <- modern_pops
    d$panel <- d$pop
  }
  d$panel <- factor(d$panel, levels = panel_levels)

  d <- d[!is.na(d[[value_col]]), ]
  if (min_count > 0) {
    keep <- d[[count_col]] >= min_count
    message(sprintf(
      "%s: dropped %d of %d windows with %s < %g",
      basename(file), sum(!keep), length(keep), count_col, min_count
    ))
    d <- d[keep, ]
  }
  missing_panels <- setdiff(panel_levels, as.character(d$panel))
  if (length(missing_panels) > 0) {
    stop(sprintf(
      "No windows in %s for: %s",
      file, paste(missing_panels, collapse = ", ")
    ))
  }

  i <- match(d$chromosome, chroms$chromosome)
  d$value <- d[[value_col]]
  d$x <- chroms$offset[i] + (d$window_pos_1 + d$window_pos_2) / 2
  d$shade <- chroms$shade[i]
  list(windows = d, chroms = chroms)
}

# ---- 3. Plot ----------------------------------------------------------------
# one panel per pair / population, stacked in modern_pops order
# y_cap, if given, is the top of the y-axis: windows above it are drawn as
# triangles at y_cap rather than at their true value.
plot_windows <- function(w, y_lab, out_file, y_cap = NULL) {
  d <- w$windows
  d$capped <- FALSE
  if (!is.null(y_cap)) {
    d$capped <- d$value > y_cap
    message(sprintf(
      "%s: %d of %d windows above the y-axis cut (%.4g) drawn as triangles",
      basename(out_file), sum(d$capped), nrow(d), y_cap
    ))
    d$value <- pmin(d$value, y_cap)
  }
  chroms <- w$chroms
  n_panels <- nlevels(d$panel)
  # every other chromosome gets a light background band
  bands <- chroms[chroms$shade == "0", ]

  p <- ggplot(d, aes(x = x, y = value, color = shade)) +
    annotate(
      "rect",
      xmin = bands$offset, xmax = bands$end, ymin = -Inf, ymax = Inf,
      fill = "grey50", alpha = 0.15
    ) +
    geom_hline(yintercept = 0, color = "grey30", size = 0.3) +
    geom_point(data = d[!d$capped, ], size = 0.4, alpha = 0.7) +
    geom_point(data = d[d$capped, ], shape = 17, size = 1.3) +
    facet_grid(panel ~ .) +
    scale_color_manual(values = c("#444e86", "#9aa3cf"), guide = "none") +
    scale_x_continuous(
      breaks = chroms$mid,
      labels = chroms$chromosome,
      expand = expansion(mult = 0.005),
      guide = guide_axis(n.dodge = 2)
    ) +
    xlab("Chromosome") +
    ylab(y_lab) +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(size = 7),
      strip.text = element_text(size = 10, face = "bold"),
      strip.background = element_rect(fill = "grey92")
    )

  if (!is.null(y_cap)) {
    p <- p + scale_y_continuous(
      limits = c(NA, y_cap),
      expand = expansion(mult = c(0.05, 0.03))
    )
  }

  pdf(out_file, width = 14, height = 1.6 * n_panels + 1)
  print(p)
  invisible(dev.off())
  message(sprintf("Wrote %s (%d panels)", out_file, n_panels))
}

fst <- read_windows(fst_file, "^avg_.*fst$", "no_snps", min_snps)
dxy <- read_windows(dxy_file, "^avg_dxy$", "no_sites", min_sites)
pi <- read_windows(pi_file, "^avg_pi$", "no_sites", min_sites)

plot_windows(fst, expression(bold(F[ST])), paste0(out_prefix, "_fst.pdf"))
plot_windows(dxy, expression(bold(d[XY])), paste0(out_prefix, "_dxy.pdf"))
# pi axis cut at the 99.9th percentile across all populations
plot_windows(
  pi, expression(bold(pi)), paste0(out_prefix, "_pi.pdf"),
  y_cap = unname(quantile(pi$windows$value, 0.999))
)
