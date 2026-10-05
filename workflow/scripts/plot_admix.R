sink(file(snakemake@log[[1]], open="wt"), type = "message")

# Adapted from PopGLen (workflow/scripts/plot_admix.R, Zachary J. Nolen):
# plots NGSadmix admixture proportions for every K and writes a convergence
# summary table, including which K is the best fit.

library(data.table)
library(reshape2)
library(ggplot2)
library(ggtext)
library(cowplot)
library(dplyr)

plot_admix <- function(kvalues, qoptlist, pop, optsumm) {

	admix <- list()

	# Generate a plot for each K
	for (i in c(1:length(kvalues))) {
		k <- kvalues[i]
		# One row per individual, in beagle column order, which is the order
		# the poplist was written in
		qopt <- read.table(qoptlist[i])
		if (nrow(qopt) != nrow(pop)) {
			stop(
				"K=", k, ": ", nrow(qopt), " individuals in the .qopt file but ",
				nrow(pop), " in the poplist"
			)
		}
		qopt$sample <- pop$sample
		qopt <- reshape2::melt(
			inner_join(qopt, pop, by = "sample"),
			id.vars = c("sample", "population", "time")
		)
		qopt$time <- factor(
			qopt$time,
			levels = c("historical", "modern"),
			labels = c("Historical", "Modern")
		)
		qopt$population <- as.factor(qopt$population)

		# Show strip only on first plot
		if (i == 1) {
			striptext <- element_text(size = 8)
		} else {
			striptext <- element_text(size = 0)
		}

		# Show sample names only on last
		if (i == length(kvalues)) {
			xtext <- element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6)
		} else {
			xtext <- element_blank()
		}

		# Bold K values of converged Ks, and mark the best-fit K
		ylabel <- paste0("K = ", k)
		if (optsumm$bestK[optsumm$k == k] == "Y") {
			ylabel <- paste0(ylabel, " (best)")
		}
		if (optsumm$converged[optsumm$k == k] == "Y") {
			ytext <- element_markdown(lineheight = 1.2, face = "bold")
		} else {
			ytext <- element_markdown(lineheight = 1.2)
		}

		# Make plot for given K
		admix[[i]] <- ggplot(qopt, aes(value, fill = variable, x = sample)) +
			geom_bar(
				position = "fill",
				stat = "identity",
				width = 1,
				color = "black",
				linewidth = 0.2
			) +
			facet_grid(
				~ population,
				scales = "free",
				space = "free"
			) +
			ylab(ylabel) +
			guides(fill = "none") +
			scale_y_continuous(expand = c(0,0), breaks = seq(0, 1, by = 0.25)) +
			theme_classic() +
			theme(
				axis.title.y = ytext,
				strip.text.x = striptext,
				strip.background = element_blank(),
				axis.title.x = element_blank(),
				axis.text.x = xtext,
				axis.line.x = element_blank(),
				axis.ticks.x = element_blank(),
				plot.margin = unit(c(0, 0, 0.1, 0), "cm")
			)
	}

	# Arrange plots in a grid with some rough evenness hopefully. Will probably be
	# a bit off due to variation in sample name length. First panel is taller for
	# the population strip, last one for the sample names.
	if (length(admix) == 1) {
		relheights <- 1
	} else {
		relheights <- c(1, rep(0.88, length(admix)-2), 1.25)
	}

	plot_grid(
		plotlist = admix,
		nrow = length(admix),
		ncol = 1,
		rel_heights = relheights
	)
}

# Read in necessary population values
pop <- as.data.frame(
	fread(
		snakemake@input[["poplist"]],
		header = TRUE,
		sep = "\t",
		select = c("sample", "population", "time")
	)
)

kvals <- as.integer(unlist(snakemake@params[["kvals"]]))
conv <- as.integer(snakemake@params[["conv"]])
thresh <- as.numeric(snakemake@params[["thresh"]])

# Summarize convergence output to table, use it for plotting. Each line of an
# optimization wrapper log is one replicate: replicate number, seed, likelihood
optwraps <- c()

for (i in c(1:length(kvals))) {
	k <- kvals[i]
	optwrap <- read.table(snakemake@input[["optwrap"]][i], header = FALSE)
	optwrap$k <- k
	optwraps <- rbind(optwraps, optwrap)
}

# deltaLike is the log-likelihood range of the top `conv` replicates; a K is
# converged when that is within `thresh`. NA (not converged) if fewer than
# `conv` replicates were run.
optsumm <- optwraps %>%
	group_by(k) %>%
	summarize(
		niter = max(V1),
		bestLike = max(V3),
		deltaLike = diff(
			range(
				sort(
					V3,
					decreasing = TRUE
				)[1:conv]
			)
		)
	) %>%
	as.data.frame()

optsumm$converged <- ifelse(
	!is.na(optsumm$deltaLike) & optsumm$deltaLike <= thresh,
	"Y",
	"N"
)

# Best-fit K: the highest K whose replicates converged
convks <- optsumm$k[optsumm$converged == "Y"]
optsumm$bestK <- "N"
if (length(convks) > 0) {
	optsumm$bestK[optsumm$k == max(convks)] <- "Y"
	message("Best-fit K (highest converged K): ", max(convks))
} else {
	message("No K converged, so no best-fit K is reported")
}

write.table(optsumm, snakemake@output[["convsumm"]], quote = FALSE, sep = "\t",
	row.names = FALSE)

# Run the plotting function
plot_admix(
	kvals,
	snakemake@input[["qopts"]],
	pop,
	optsumm
)

# Set reasonable dimensions based on number of samples and Ks
if (nrow(pop) <= 50) {
	plotwidth = 6
} else {
	plotwidth = (6/45)*nrow(pop)
}

plotheight = 1.65 + (1.44 * max(length(kvals)-2, 0)) + 2

# Save plot
ggsave(
	snakemake@output[["plot"]],
	width = plotwidth,
	height = plotheight,
	units = "in"
)
