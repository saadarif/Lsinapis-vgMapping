sink(file(snakemake@log[[1]], open = "wt"), type = "message")

# GONE2 Ne trajectories of every subset of individuals (9b): the full set of
# n_max individuals in bold black, smaller subsets as thin lines in lighter
# colours by number of individuals. Page 1 has every subset, page 2 one panel
# per subset size (each with the full set for reference). Ne files are
# tab-separated with a header, generation in the first column and Ne in the
# second (Ne_diploids without -x, N_T with -x); a header-only file means GONE2
# gave no estimate for that subset (reason in its STATS file).

library(ggplot2)

files <- unlist(snakemake@input[["ne"]])
ids <- unlist(snakemake@params[["subsets"]])
n_max <- as.integer(snakemake@params[["n_max"]])
subsets <- read.table(snakemake@input[["subsets"]], header = TRUE, sep = "\t", colClasses = "character")

read_ne <- function(path, id) {
	d <- read.table(path, header = TRUE, sep = "\t", check.names = FALSE)
	if (ncol(d) < 2) {
		stop(path, ": expected a header and at least two columns (generation, Ne)")
	}
	data.frame(generation = d[[1]], ne = d[[2]], subset = rep(id, nrow(d)))
}

ne <- do.call(rbind, Map(read_ne, files, ids))
ne$n <- as.integer(subsets$n[match(ne$subset, subsets$subset)])
# Light colours per subset size, from orange (fewest individuals) to blue
sizes <- sort(unique(ne$n[ne$n < n_max]))
size_colours <- setNames(colorRampPalette(c("#fdae6b", "#9ecae1"))(max(length(sizes), 1))[seq_along(sizes)], sizes)
no_estimate <- setdiff(ids, unique(ne$subset))
message(length(ids), " subsets, ", length(no_estimate), " without a GONE2 estimate",
	if (length(no_estimate) > 0) paste0(": ", paste(no_estimate, collapse = ", ")) else "")

full <- ne[ne$n == n_max, ]
part <- ne[ne$n < n_max, ]

subtitle <- paste0(
	snakemake@params[["subtitle"]],
	"\nBlack: all ", n_max, " individuals. ",
	length(unique(part$subset)), " smaller subsets with an estimate",
	if (length(no_estimate) > 0) paste0(", ", length(no_estimate), " without (see their STATS files)") else ""
)

trajectories <- function(part, full) {
	p <- ggplot(mapping = aes(x = generation, y = ne))
	if (nrow(part) > 0) {
		p <- p +
			geom_line(data = part, aes(group = subset, colour = factor(n, levels = sizes)), linewidth = 0.4, alpha = 0.7) +
			scale_colour_manual(values = size_colours, name = "individuals") +
			guides(colour = guide_legend(override.aes = list(linewidth = 1.5, alpha = 1)))
	}
	if (nrow(full) > 0) {
		p <- p + geom_line(data = full, colour = "black", linewidth = 1.3)
	}
	if (nrow(part) + nrow(full) > 0) {
		p <- p + scale_y_log10()
	}
	p +
		labs(x = "Generations ago", y = "Ne (log scale)", title = snakemake@params[["title"]], subtitle = subtitle) +
		theme_bw() +
		theme(legend.position = "bottom", plot.subtitle = element_text(size = 7))
}

pdf(snakemake@output[["pdf"]], width = 8, height = 6)
print(trajectories(part, full))
if (nrow(part) > 0) {
	full_panels <- do.call(rbind, lapply(sizes, function(k) if (nrow(full) > 0) transform(full, n = k) else full))
	p2 <- trajectories(part, full_panels) +
		facet_wrap(~ n, labeller = labeller(n = function(k) paste(k, "individuals"))) +
		theme(legend.position = "none")
	print(p2)
}
invisible(dev.off())
