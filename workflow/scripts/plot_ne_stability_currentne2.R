sink(file(snakemake@log[[1]], open = "wt"), type = "message")

# currentNe2 estimates of every subset of individuals (9b): one table of all
# estimates, and Ne against the number of individuals with a linear regression
# (Ne ~ n, least squares), one panel per estimate type.
#
# A currentNe2 output file has one block per estimate, each starting with a
# "# Ne estimation ..." line (e.g. "by integration over the whole genome",
# "based only on LD between chromosomes"), followed by "# Ne point estimate:"
# and the 50% / 90% CI limits, each value on the line after its label. When no
# estimate can be made there is no such block, only a line saying why.

library(ggplot2)

files <- unlist(snakemake@input[["txt"]])
ids <- unlist(snakemake@params[["subsets"]])
n_max <- as.integer(snakemake@params[["n_max"]])
subsets <- read.table(snakemake@input[["subsets"]], header = TRUE, sep = "\t", colClasses = "character")

parse_currentne2 <- function(path, id) {
	x <- readLines(path)
	heads <- grep("^# Ne estimation", x)
	if (length(heads) == 0) {
		why <- grep("cannot be estimated", x, value = TRUE)
		return(data.frame(
			subset = id, estimate = NA_character_, ne = NA_real_,
			ci50_lower = NA_real_, ci50_upper = NA_real_, ci90_lower = NA_real_, ci90_upper = NA_real_,
			note = if (length(why) > 0) sub("^#\\s*", "", why[1]) else "no Ne estimate in the output"
		))
	}
	ends <- c(heads[-1] - 1, length(x))
	do.call(rbind, lapply(seq_along(heads), function(k) {
		block <- x[heads[k]:ends[k]]
		value_after <- function(pattern) {
			i <- grep(pattern, block)
			if (length(i) == 0) NA_real_ else suppressWarnings(as.numeric(block[i[1] + 1]))
		}
		label <- sub("[.:]\\s*$", "", sub("^# Ne estimation\\s*", "", block[1]))
		data.frame(
			subset = id, estimate = label, ne = value_after("^# Ne point estimate"),
			ci50_lower = value_after("^# Lower (limit|bound) of 50% CI"),
			ci50_upper = value_after("^# Upper (limit|bound) of 50% CI"),
			ci90_lower = value_after("^# Lower (limit|bound) of 90% CI"),
			ci90_upper = value_after("^# Upper (limit|bound) of 90% CI"),
			note = ""
		)
	}))
}

est <- do.call(rbind, Map(parse_currentne2, files, ids))
rownames(est) <- NULL
est$n <- as.integer(subsets$n[match(est$subset, subsets$subset)])
est$dropped <- subsets$dropped[match(est$subset, subsets$subset)]
est <- est[, c("subset", "n", "dropped", "estimate", "ne", "ci50_lower", "ci50_upper", "ci90_lower", "ci90_upper", "note")]
write.table(est, snakemake@output[["tsv"]], sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

ok <- est[!is.na(est$ne), ]
n_none <- length(setdiff(ids, unique(ok$subset)))
message(length(ids), " subsets, ", n_none, " without a currentNe2 estimate")

# Linear regression per estimate type, when there are at least 3 estimates
# over at least 2 subset sizes
fits <- do.call(rbind, lapply(split(ok, ok$estimate), function(d) {
	if (nrow(d) < 3 || length(unique(d$n)) < 2) {
		return(data.frame(estimate = d$estimate[1], label = "too few estimates for a regression"))
	}
	f <- summary(lm(ne ~ n, data = d))
	cf <- f$coefficients
	message(d$estimate[1], ": Ne = ", signif(cf[1, 1], 4), " + ", signif(cf[2, 1], 4), " * n, R2 ", signif(f$r.squared, 3), ", p ", signif(cf[2, 4], 3))
	data.frame(
		estimate = d$estimate[1],
		label = sprintf("Ne = %.4g %s %.4g n\nR² = %.3f, p = %.3g", cf[1, 1], ifelse(cf[2, 1] < 0, "-", "+"), abs(cf[2, 1]), f$r.squared, cf[2, 4])
	)
}))

subtitle <- paste0(
	snakemake@params[["subtitle"]],
	"\nBlack: all ", n_max, " individuals; line: linear regression of Ne on n (95% CI band)",
	if (n_none > 0) paste0("\n", n_none, " of ", length(ids), " subsets without an estimate (see the summary table)") else ""
)

p <- ggplot(ok, aes(x = n, y = ne))
if (nrow(ok) > 0) {
	p <- p +
		geom_point(data = ok[ok$n < n_max, ], colour = "#6baed6", alpha = 0.5, size = 1.5,
			position = position_jitter(width = 0.12, height = 0, seed = 1)) +
		geom_smooth(method = "lm", formula = y ~ x, colour = "#d1495b", fill = "#d1495b", alpha = 0.15, linewidth = 0.8) +
		geom_point(data = ok[ok$n == n_max, ], colour = "black", size = 3) +
		geom_text(data = fits, aes(x = -Inf, y = Inf, label = label), hjust = -0.05, vjust = 1.2, size = 3, inherit.aes = FALSE) +
		facet_wrap(~ estimate, scales = "free_y", ncol = 1)
} else {
	p <- p + annotate("text", x = n_max, y = 0, label = "No currentNe2 estimate for any subset")
}
p <- p +
	scale_x_continuous(breaks = sort(unique(est$n))) +
	labs(x = "Number of individuals", y = "Ne (currentNe2 point estimate)", title = snakemake@params[["title"]], subtitle = subtitle) +
	theme_bw() +
	theme(plot.subtitle = element_text(size = 7))

pdf(snakemake@output[["pdf"]], width = 7, height = if (length(unique(ok$estimate)) > 1) 8 else 5)
print(p)
invisible(dev.off())
