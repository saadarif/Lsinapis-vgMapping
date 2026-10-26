sink(file(snakemake@log[[1]], open = "wt"), type = "message")

# GONE2 Ne backward in time for one population VCF: the run without -x and the
# run with -x (metapopulation) as two lines. Each Ne file is a tab-separated
# table with a header, generation in the first column and Ne in the second
# (Ne_diploids without -x, N_T with -x). Further columns (e.g. min/max with
# GONE2 -E) are ignored. A header-only table means GONE2 gave no estimate for
# that run (the reason is in its STATS file); the run is then named in the
# subtitle instead of drawn.

library(ggplot2)

runs <- c(nox = "without -x", x = "with -x (metapopulation)")

read_ne <- function(path, run) {
	d <- read.table(path, header = TRUE, sep = "\t", check.names = FALSE)
	if (ncol(d) < 2) {
		stop(path, ": expected a header and at least two columns (generation, Ne)")
	}
	message(run, ": ", basename(path), ", columns ", paste(colnames(d), collapse = ", "), ", ", nrow(d), " rows")
	data.frame(generation = d[[1]], ne = d[[2]], run = rep(run, nrow(d)))
}

ne <- do.call(rbind, lapply(names(runs), function(r) read_ne(snakemake@input[[r]], runs[[r]])))
ne$run <- factor(ne$run, levels = unname(runs))

no_estimate <- setdiff(unname(runs), unique(as.character(ne$run)))
subtitle <- snakemake@params[["subtitle"]]
if (length(no_estimate) > 0) {
	subtitle <- paste0(subtitle, "\nNo GONE2 estimate ", paste(no_estimate, collapse = " or "), " (see the STATS file)")
}

p <- ggplot(ne, aes(x = generation, y = ne, colour = run)) +
	scale_colour_manual(values = setNames(c("#1b6ca8", "#d1495b"), unname(runs)), drop = TRUE) +
	labs(
		x = "Generations ago",
		y = "Ne (log scale)",
		colour = NULL,
		title = snakemake@params[["title"]],
		subtitle = subtitle
	) +
	theme_bw() +
	theme(legend.position = "bottom", plot.subtitle = element_text(size = 7))
if (nrow(ne) > 0) {
	p <- p + geom_line(linewidth = 0.8) + scale_y_log10()
} else {
	p <- p + annotate("text", x = 0, y = 0, label = "No GONE2 estimate for either run")
}

pdf(snakemake@output[["pdf"]], width = 7, height = 5)
print(p)
invisible(dev.off())
