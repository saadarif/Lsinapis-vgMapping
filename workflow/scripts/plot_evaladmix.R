sink(file(snakemake@log[[1]], open="wt"), type = "message")

# Adapted from PopGLen (workflow/scripts/plot_evaladmix.R, Zachary J. Nolen):
# heatmap of evalAdmix's pairwise correlation of residuals for one K.
# Upper triangle is the correlation for each pair of individuals, lower triangle
# the mean for each pair of populations.

# evalAdmix's plotting functions, shipped in the evaladmix container
source("/usr/local/bin/visFuns.R")

plot_evaladmix <- function(pop,qopt,corres,k) {

	# read population labels and estimated admixture proportions. Rows of the
	# poplist, the .qopt file and the correlation matrix are all individuals
	# in beagle column order
	pop<-read.table(pop, header = TRUE, sep = "\t")
	q<-read.table(qopt)

	# read in correlation matrix
	r<-as.matrix(read.table(corres))

	if (nrow(pop) != nrow(q) || nrow(pop) != nrow(r)) {
		stop(
			"K=", k, ": ", nrow(pop), " individuals in the poplist, ", nrow(q),
			" in the .qopt file and ", nrow(r), " in the .corres file"
		)
	}

	# order according to population
	ord<-orderInds(pop = as.vector(pop[,2]), q = q)

	# Plot correlation of residuals
	plotCorRes(cor_mat = r, pop = as.vector(pop[,2]), ord=ord,
		title=paste0("Evaluation of admixture proportions with K=",k),
		max_z=0.1, min_z=-0.1)

}

pdf(snakemake@output[["plot"]], width = 11.8, height = 9.8)
plot_evaladmix(snakemake@input[["pops"]],snakemake@input[["qopt"]],
	snakemake@input[["corres"]],snakemake@wildcards[["kvalue"]])
dev.off()
