# ==============================================================================
# POPULATION STRUCTURE: INDIVIDUAL ADMIXTURE PROPORTIONS (NGSadmix + evalAdmix)
# Adapted from PopGlen by Zachary Nolen
#   (i)   NGSadmix for every K from 1 to params: run_admixture: max_k, on the
#         exact beagle file PCAngsd runs on (pca_exclude_beagle in
#         6a_Strucure_angsdPCA.smk: LD-pruned, unlinked SNPs, with params:
#         run_pca: exclude_samples already dropped)
#   (ii)  Plot admixture proportions for all K and summarise convergence
#   (iii) evalAdmix on the best replicate of every K: pairwise correlation of
#         residuals as a check of model fit, plotted per K
#
# NGSadmix is an EM algorithm started from a random seed, so a single run can
# stop in a local optimum. workflow/scripts/ngsadmix.sh (from PopGLen) repeats
# each K: at least `minreps` replicates, then keeps going until the `conv`
# highest-likelihood replicates are within `thresh` log-likelihood units of each
# other, or `reps` replicates have been run. Only the highest-likelihood
# replicate is kept. The best-fit K is then taken as the highest K whose
# replicates met that convergence criterion (Pecnerova et al. 2021, Curr. Biol.),
# reported in the convergence summary table written by plot_admix.
#
# Like the PCA, this takes no called genotypes: the pruned beagle file's
# genotype likelihoods are used directly, and the same relatives are excluded
# for the same reason (close relatives otherwise get a cluster of their own).
# ==============================================================================
ADMIX_PARAMS = config.get("params", {}).get("run_admixture", {})
ADMIX_MAX_K = int(ADMIX_PARAMS.get("max_k", 5))
ADMIX_MINREPS = ADMIX_PARAMS.get("minreps", 20)
ADMIX_REPS = ADMIX_PARAMS.get("reps", 100)
ADMIX_THRESH = ADMIX_PARAMS.get("thresh", 2)
ADMIX_THREADS = ADMIX_PARAMS.get("threads", 4)
ADMIX_EXTRA = ADMIX_PARAMS.get("extra", "")
ADMIX_POPFILE = ADMIX_PARAMS.get("poplist", "Pops.txt")

# Number of top replicates compared for convergence. Not configurable:
# ngsadmix.sh hardcodes the comparison to the top 3 (head -n 3), its `conv`
# setting only changes the wording of its log messages.
ADMIX_CONV = 3

if ADMIX_MAX_K < 1:
    print("\nERROR: params: run_admixture: max_k must be an integer >= 1.\n")
    sys.exit(1)

ADMIX_KVALUES = list(range(1, ADMIX_MAX_K + 1))

# Same samples, in the same order, as the columns of the beagle file written by
# pca_exclude_beagle. NGSadmix writes one .qopt row per individual in beagle
# column order, so this is also the row order of every .qopt file.
ADMIX_SAMPLES = PCA_KEEP_SAMPLES

# Nothing of its own is added to the tag: K is in the file names, and the
# excluded samples are already carried by PCA_TAG.
ADMIX_TAG = f"{PRUNE_TAG}.{PCA_TAG}"
ADMIX_PREFIX = f"results/structure/ngsadmix/all.{REF_NAME}.{ADMIX_TAG}"

# Requested from the Snakefile target list (section 12). The K range is part of
# the summary file names so that raising max_k later produces new summary files
# rather than leaving stale ones in place (rerun-triggers is mtime only).
ADMIX_TARGETS = [
    f"{ADMIX_PREFIX}.K1-{ADMIX_MAX_K}.admix.pdf",
    f"{ADMIX_PREFIX}.K1-{ADMIX_MAX_K}.convergence_summary.tsv",
] + [f"{ADMIX_PREFIX}.K{k}.evaladmix.pdf" for k in ADMIX_KVALUES]


rule admix_poplist:
    """
    Sample / population / time table for the admixture plots, one row per
    individual in the same order as the beagle columns (and so as the .qopt
    rows). Population comes from params: run_admixture: poplist (two columns,
    sample_id and population, no header), time (modern/historical) from which
    sample sheet the sample was listed in.
    """
    input:
        pops=ADMIX_POPFILE,
    output:
        poplist="results/structure/ngsadmix/all.{ref_name}." + ADMIX_TAG + ".poplist.tsv",
    run:
        pop_of = {}
        with open(input.pops) as fh:
            for line in fh:
                fields = line.split()
                if len(fields) >= 2:
                    pop_of[fields[0]] = fields[1]
        missing = [sid for sid in ADMIX_SAMPLES if sid not in pop_of]
        if missing:
            raise ValueError(
                f"No population in {input.pops} for: {', '.join(missing)}"
            )
        time_of = dict(zip(samples_df["sample_id"], samples_df["source"]))
        with open(output.poplist, "w") as fh:
            fh.write("sample\tpopulation\ttime\n")
            for sid in ADMIX_SAMPLES:
                fh.write(f"{sid}\t{pop_of[sid]}\t{time_of[sid]}\n")


rule ngsadmix:
    """
    Individual admixture proportions for one K. ngsadmix.sh runs replicates
    until the top 3 are within `thresh` log-likelihood units of each other
    (checked once `minreps` have finished) or `reps` have been run, and keeps
    only the highest-likelihood replicate.

    .log is NGSadmix's own run log for the kept replicate (holds its 'best
    like='), tracked as output the same way PCAngsd's is in pca_pcangsd.
    optimization_wrapper.log has one line per replicate (replicate, seed,
    log-likelihood) and is what convergence is summarised from.
    """
    input:
        beagle="results/structure/pca/all.{ref_name}." + ADMIX_TAG + ".beagle.gz",
    output:
        qopt="results/structure/ngsadmix/all.{ref_name}." + ADMIX_TAG + ".K{kvalue}.qopt",
        fopt="results/structure/ngsadmix/all.{ref_name}." + ADMIX_TAG + ".K{kvalue}.fopt.gz",
        nglog="results/structure/ngsadmix/all.{ref_name}." + ADMIX_TAG + ".K{kvalue}.log",
        log="results/structure/ngsadmix/all.{ref_name}." + ADMIX_TAG + ".K{kvalue}.optimization_wrapper.log",
    params:
        prefix=lambda wildcards, output: output.qopt[: -len(".qopt")],
        extra=ADMIX_EXTRA,
        reps=ADMIX_REPS,
        minreps=ADMIX_MINREPS,
        thresh=ADMIX_THRESH,
        conv=ADMIX_CONV,
    wildcard_constraints:
        kvalue=r"\d+",
    log:
        "logs/structure/ngsadmix_{ref_name}.K{kvalue}." + ADMIX_TAG + ".log",
    benchmark:
        "benchmarks/structure/ngsadmix_{ref_name}.K{kvalue}." + ADMIX_TAG + ".benchmark"
    conda:
        "../envs/angsd.yaml"
    threads: ADMIX_THREADS
    script:
        "../scripts/ngsadmix.sh"


rule plot_admix:
    """
    Admixture proportions for every K in one figure, plus the convergence
    summary table: per K, how many replicates were run, the log-likelihood range
    of the top 3, whether that counts as converged, and which K is the best fit
    (the highest converged K).
    """
    input:
        qopts=expand(
            "results/structure/ngsadmix/all.{{ref_name}}." + ADMIX_TAG + ".K{kvalue}.qopt",
            kvalue=ADMIX_KVALUES,
        ),
        optwrap=expand(
            "results/structure/ngsadmix/all.{{ref_name}}." + ADMIX_TAG + ".K{kvalue}.optimization_wrapper.log",
            kvalue=ADMIX_KVALUES,
        ),
        poplist="results/structure/ngsadmix/all.{ref_name}." + ADMIX_TAG + ".poplist.tsv",
    output:
        plot="results/structure/ngsadmix/all.{ref_name}." + ADMIX_TAG + f".K1-{ADMIX_MAX_K}.admix.pdf",
        convsumm="results/structure/ngsadmix/all.{ref_name}." + ADMIX_TAG + f".K1-{ADMIX_MAX_K}.convergence_summary.tsv",
    params:
        kvals=ADMIX_KVALUES,
        thresh=ADMIX_THRESH,
        conv=ADMIX_CONV,
    log:
        "logs/structure/plot_admix_{ref_name}.log",
    conda:
        "../envs/r_admix.yaml"
    script:
        "../scripts/plot_admix.R"


rule evaladmix:
    """
    Model fit for one K with evalAdmix: the pairwise correlation of residuals
    between individuals, given the admixture proportions (.qopt) and ancestral
    allele frequencies (.fopt.gz) of the best NGSadmix replicate. With a good
    fit the residuals are uncorrelated (close to 0); individuals with positively
    correlated residuals share ancestry the K clusters do not account for.

    evalAdmix has no bioconda/conda-forge package, so this uses PopGLen's
    container (v0.961), which also ships the visFuns.R plotting functions
    plot_evaladmix needs.
    """
    input:
        beagle="results/structure/pca/all.{ref_name}." + ADMIX_TAG + ".beagle.gz",
        qopt="results/structure/ngsadmix/all.{ref_name}." + ADMIX_TAG + ".K{kvalue}.qopt",
        fopt="results/structure/ngsadmix/all.{ref_name}." + ADMIX_TAG + ".K{kvalue}.fopt.gz",
    output:
        corres="results/structure/ngsadmix/all.{ref_name}." + ADMIX_TAG + ".K{kvalue}.corres",
    wildcard_constraints:
        kvalue=r"\d+",
    log:
        "logs/structure/evaladmix_{ref_name}.K{kvalue}." + ADMIX_TAG + ".log",
    benchmark:
        "benchmarks/structure/evaladmix_{ref_name}.K{kvalue}." + ADMIX_TAG + ".benchmark"
    container:
        "docker://ghcr.io/zjnolen/evaladmix:0.961"
    threads: ADMIX_THREADS
    shell:
        """
        evalAdmix -beagle {input.beagle} -fname {input.fopt} \
            -qname {input.qopt} -o {output.corres} -P {threads} &> {log}
        """


rule plot_evaladmix:
    """
    Heatmap of evalAdmix's correlation of residuals for one K, individuals
    grouped by population.
    """
    input:
        corres="results/structure/ngsadmix/all.{ref_name}." + ADMIX_TAG + ".K{kvalue}.corres",
        qopt="results/structure/ngsadmix/all.{ref_name}." + ADMIX_TAG + ".K{kvalue}.qopt",
        pops="results/structure/ngsadmix/all.{ref_name}." + ADMIX_TAG + ".poplist.tsv",
    output:
        plot="results/structure/ngsadmix/all.{ref_name}." + ADMIX_TAG + ".K{kvalue}.evaladmix.pdf",
    wildcard_constraints:
        kvalue=r"\d+",
    log:
        "logs/structure/plot_evaladmix_{ref_name}.K{kvalue}.log",
    container:
        "docker://ghcr.io/zjnolen/evaladmix:0.961"
    script:
        "../scripts/plot_evaladmix.R"
