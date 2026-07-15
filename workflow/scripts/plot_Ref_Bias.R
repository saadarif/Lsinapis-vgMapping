#!/usr/bin/env Rscript

# ---------------------------------------------------------------------------
# Reference-bias boxplots across genotyping folders
#
# For each "genotyping*" folder found directly under a provided directory,
# reads TWO ref_bias files:
#   * individual calls: ...indCall.*.allsites.fmiss0.0.bcf.stats.ref_bias
#   * joint calls:      ...jointCall.*.allsites.fmiss0.0.bcf.stats.ref_bias
# (where * is variable text).
#
# Each sample is classified as Historical (prefix "OX" or "NHM") or Modern,
# and by country of origin from its prefix:
#   LsSpa -> Spain, LsKaz -> Kazakhstan, LsSwe -> Sweden, else Great Britain.
#
# Output is a single figure with one boxplot panel per
# (call type x genotyping method) combination. Boxplots are grouped by
# Historical vs Modern; every sample is overlaid as a point coloured by country.
#
# Usage:
#   Rscript ref_bias_boxplots.R <input_directory> [output_plot.png]
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(ggplot2))

# ---- 1. Arguments ---------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript ref_bias_boxplots.R <input_directory> [output_plot.png]")
}
base_dir <- args[1]
out_file <- if (length(args) >= 2) args[2] else "ref_bias_boxplots.png"

if (!dir.exists(base_dir)) {
  stop(sprintf("Input directory does not exist: %s", base_dir))
}

# ---- 2. Find the three genotyping folders ---------------------------------
subdirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
geno_dirs <- subdirs[grepl("^genotyping_", basename(subdirs))]

if (length(geno_dirs) == 0) {
  stop("No folders starting with 'genotyping' found in the provided directory.")
}
message(sprintf("Found %d genotyping folder(s):", length(geno_dirs)))
message(paste0("  - ", basename(geno_dirs), collapse = "\n"))

# ---- 3. Read both call-type files from each folder ------------------------
# The wildcard between the call tag and "allsites" may be empty (e.g.
# ".indCall.allsites...") or contain extra text (".indCall.notrans.allsites..."),
# so ".*" is used to match zero or more characters.
suffix <- "\\..*allsites\\.fmiss0\\.0\\.bcf\\.stats\\.ref_bias$"
call_defs <- list(
  list(label = "Individual call", pattern = paste0("indCall", suffix)),
  list(label = "Joint call", pattern = paste0("jointCall", suffix))
)

read_ref_bias <- function(folder, call_label, pattern) {
  f <- list.files(folder, pattern = pattern, full.names = TRUE)

  if (length(f) == 0) {
    warning(sprintf(
      "No '%s' file (%s) in %s -- skipping.",
      call_label,
      pattern,
      folder
    ))
    return(NULL)
  }
  if (length(f) > 1) {
    warning(sprintf(
      "Multiple '%s' files in %s; using the first: %s",
      call_label,
      folder,
      basename(f[1])
    ))
    f <- f[1]
  }

  df <- read.table(
    f,
    header = FALSE,
    stringsAsFactors = FALSE,
    col.names = c("sample", "ref_bias")
  )
  df$ref_bias <- as.numeric(df$ref_bias)

  # Historical if the sample starts with OX or NHM, otherwise Modern
  df$group <- ifelse(grepl("^(OX|NHM)", df$sample), "Historical", "Modern")

  # Country of origin from the sample prefix (default: Great Britain)
  df$country <- "Great Britain"
  df$country[grepl("^LsSpa", df$sample)] <- "Spain"
  df$country[grepl("^LsKaz", df$sample)] <- "Kazakhstan"
  df$country[grepl("^LsSwe", df$sample)] <- "Sweden"

  df$folder <- basename(folder)
  df$callType <- call_label
  df
}

records <- list()
for (d in geno_dirs) {
  for (cd in call_defs) {
    records[[length(records) + 1L]] <- read_ref_bias(d, cd$label, cd$pattern)
  }
}
all_data <- do.call(rbind, records)

if (is.null(all_data) || nrow(all_data) == 0) {
  stop("No data could be read from any genotyping folder.")
}

# Fix factor orderings for consistent panels, axes and colours
all_data$folder <- factor(all_data$folder, levels = basename(geno_dirs))
all_data$callType <- factor(
  all_data$callType,
  levels = c("Individual call", "Joint call")
)
all_data$group <- factor(all_data$group, levels = c("Historical", "Modern"))
all_data$country <- factor(
  all_data$country,
  levels = c("Great Britain", "Spain", "Sweden", "Kazakhstan")
)

# ---- 4. Plot: call type (rows) x genotyping method (columns) --------------
# Boxplots grouped by Historical vs Modern; every sample overlaid as a point
# coloured by country. outlier.shape = NA avoids drawing outliers twice, since
# all points are already shown.
country_cols <- c(
  "Great Britain" = "#444e86",
  "Spain" = "#ce4717",
  "Sweden" = "#ff7730",
  "Kazakhstan" = "#25eb81"
)

p <- ggplot(all_data, aes(x = group, y = ref_bias)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, fill = "grey95") +
  geom_point(
    aes(colour = country),
    position = position_jitter(width = 0.18, height = 0, seed = 42),
    size = 2.0,
    alpha = 0.8
  ) +
  facet_grid(callType ~ folder) +
  scale_colour_manual(values = country_cols, drop = FALSE) +
  labs(
    x = NULL,
    y = "Reference bias",
    colour = "Country of origin",
    title = "Reference bias by sample type, call method and genotyping run"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey90"),
    panel.grid.major.x = element_blank()
  )

# ---- 5. Save the plot -----------------------------------------------------
# Write with R's built-in grDevices devices (png/pdf) rather than letting
# ggsave pick an add-on device such as ragg or Cairo. Those add-on devices
# throw "Graphics API version mismatch" when they were compiled against a
# different R version than the one currently running; the built-in devices
# ship with R and cannot mismatch.
n_col <- nlevels(droplevels(all_data$folder)) # genotyping methods present
n_row <- nlevels(droplevels(all_data$callType)) # call types present
plot_w <- 3.5 * n_col # inches
plot_h <- 3.2 * n_row + 1 # inches (+1 for the bottom legend)
ext <- tolower(tools::file_ext(out_file))

if (ext == "pdf") {
  pdf(out_file, width = plot_w, height = plot_h)
} else {
  if (ext != "png") {
    out_file <- paste0(tools::file_path_sans_ext(out_file), ".png")
    message(sprintf(
      "Unrecognised extension; writing a PNG to %s instead.",
      out_file
    ))
  }
  png(out_file, width = plot_w, height = plot_h, units = "in", res = 300)
}
print(p)
invisible(dev.off())

message(sprintf("Saved plot to: %s", normalizePath(out_file, mustWork = FALSE)))
