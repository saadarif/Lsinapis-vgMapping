sink(file(snakemake@log[[1]], open = "wt"), type = "message")

# GONE2 Ne trajectories for one subset size n (9b): the full set of n_max
# individuals as a dark line on top, every subset of n individuals as a light
# line. Ne files are tab-separated with a header, generation in the first
# column and Ne in the second (Ne_diploids without -x, N_T with -x); a
# header-only file means GONE2 gave no estimate for that subset (reason in its
# STATS file).

library(ggplot2)

n <- as.integer(snakemake@params[["n"]])
n_max <- as.integer(snakemake@params[["n_max"]])
ids <- unlist(snakemake@params[["subsets"]])
sub_files <- unlist(snakemake@input[["subsets"]])

read_ne <- function(path, id) {
	d <- read.table(path, header = TRUE, sep = "\t", check.names = FALSE)
	if (ncol(d) < 2) {
		stop(path, ": expected a header and at least two columns (generation, Ne)")
	}
	data.frame(generation = d[[1]], ne = d[[2]], subset = rep(id, nrow(d)))
}

full <- read_ne(unlist(snakemake@input[["full"]]), "full")
part <- if (length(sub_files) > 0) do.call(rbind, Map(read_ne, sub_files, ids)) else full[0, ]
no_estimate <- setdiff(ids, unique(part$subset))
message("n = ", n, ": ", length(ids), " subsets, ", length(no_estimate), " without a GONE2 estimate",
	if (length(no_estimate) > 0) paste0(": ", paste(no_estimate, collapse = ", ")) else "")

full_label <- paste0("all ", n_max, " individuals")
part_label <- paste0("subsets of ", n, " individuals (", length(ids) - length(no_estimate), ")")
colours <- setNames(c("#08306b", "#9ecae1"), c(full_label, part_label))

subtitle <- snakemake@params[["subtitle"]]
if (nrow(full) == 0) {
	subtitle <- paste0(subtitle, "\nNo GONE2 estimate for the full set (see its STATS file)")
}
if (length(no_estimate) > 0) {
	subtitle <- paste0(subtitle, "\n", length(no_estimate), " of ", length(ids), " subsets without a GONE2 estimate (see their STATS files)")
}

p <- ggplot(mapping = aes(x = generation, y = ne))
if (nrow(part) > 0) {
	p <- p + geom_line(data = part, aes(group = subset, colour = part_label), linewidth = 0.4, alpha = 0.8)
}
if (nrow(full) > 0) {
	p <- p + geom_line(data = full, aes(colour = full_label), linewidth = 1.3)
}
if (nrow(part) + nrow(full) > 0) {
	p <- p + scale_y_log10()
} else {
	p <- p + annotate("text", x = 0, y = 0, label = "No GONE2 estimate")
}
p <- p +
	scale_colour_manual(values = colours, breaks = names(colours), name = NULL) +
	labs(x = "Generations ago", y = "Ne (log scale)", title = snakemake@params[["title"]], subtitle = subtitle) +
	theme_bw() +
	theme(legend.position = "bottom", plot.subtitle = element_text(size = 7))

pdf(snakemake@output[["pdf"]], width = 7, height = 5)
print(p)
invisible(dev.off())
