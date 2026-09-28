# ==============================================================================
# POPULATION STRUCTURE: PCA FROM GENOTYPE LIKELIHOODS (PCAngsd)
#   (i)  Drop any individuals listed under params: run_pca: exclude_samples from
#        the LD-pruned, unlinked-SNP beagle file built in 5_LD_estimation.smk
#   (ii) Run PCA directly on genotype likelihoods with PCAngsd
#
# Individuals are columns in an ANGSD beagle file, three per sample (one column
# per possible genotype: AA, Aa, aa), after the first three columns
# (marker, allele1, allele2). Column order follows the order ANGSD was given the
# BAMs in, which is KEEP_SAMPLES_REL (angsd_bamlist_relatedness in
# 4_relatedness.smk writes the bamlist/.ids files in that order, and ANGSD
# preserves it as the beagle column order) -- confirmed directly against a real
# beagle file's header (marker allele1 allele2 Ind0 Ind0 Ind0 Ind1 Ind1 Ind1 ...)
# and against a PCAngsd run reporting the expected sample count.
#
# Excluding related individuals before PCA is standard practice: close relatives
# otherwise pull principal components toward themselves (shared drift/IBD) rather
# than reflecting broader population structure. params: run_relatedness's
# ibsrelate.R0R1KING.tsv (4_relatedness.smk) is the usual source for who to list.
# ==============================================================================
PCA_PARAMS = config.get("params", {}).get("run_pca", {})
PCA_THREADS = PCA_PARAMS.get("threads", 8)
PCA_EXTRA = PCA_PARAMS.get("extra", "")

# Only samples both requested for exclusion AND actually present in the beagle
# file (i.e. in KEEP_SAMPLES_REL, set in 4_relatedness.smk) affect anything;
# anything else is silently ignored so the config entry doesn't need to be kept
# in sync with params: run_relatedness: drop_samples by hand.
PCA_EXCLUDE_SAMPLES = sorted(
    s for s in PCA_PARAMS.get("exclude_samples", []) if s in KEEP_SAMPLES_REL
)

PCA_TAG = "excl-" + ("-".join(PCA_EXCLUDE_SAMPLES) if PCA_EXCLUDE_SAMPLES else "none")
PCA_PREFIX = f"results/structure/pca/all.{REF_NAME}.{PRUNE_TAG}.{PCA_TAG}"

# 1-based column numbers to drop from the beagle file for every excluded sample:
# marker=1, allele1=2, allele2=3, then 3 genotype-likelihood columns per
# individual starting at column 4, in KEEP_SAMPLES_REL order.
PCA_EXCLUDE_COLS = []
for _sid in PCA_EXCLUDE_SAMPLES:
    _start = KEEP_SAMPLES_REL.index(_sid) * 3 + 4
    PCA_EXCLUDE_COLS.extend([_start, _start + 1, _start + 2])
PCA_EXCLUDE_COLS_STR = ",".join(str(c) for c in sorted(PCA_EXCLUDE_COLS))

# Requested from the Snakefile target list (section 11)
PCA_TARGETS = [
    f"{PCA_PREFIX}.cov",
]


rule pca_exclude_beagle:
    """
    Drops params: run_pca: exclude_samples from the LD-pruned beagle file before
    PCA. If nothing is excluded, this is just a pass-through symlink.
    """
    input:
        beagle="results/ld/all.{ref_name}." + PRUNE_TAG + ".pruned.beagle.gz",
    output:
        beagle="results/structure/pca/all.{ref_name}." + PRUNE_TAG + "." + PCA_TAG + ".beagle.gz",
    params:
        cols=PCA_EXCLUDE_COLS_STR,
    log:
        "logs/structure/pca_exclude_beagle_{ref_name}.log",
    # No conda: env: cut/zcat/gzip/ln are plain coreutils, same as e.g.
    # individual_heterozygosity in 3_diversity_stats.smk.
    shell:
        """
        if [ -n "{params.cols}" ]; then
            zcat {input.beagle} | cut -f{params.cols} --complement | gzip > {output.beagle} 2> {log}
        else
            ln -sf "$(realpath {input.beagle})" {output.beagle} 2> {log}
        fi
        """


rule pca_pcangsd:
    """
    PCA directly from genotype likelihoods with PCAngsd: no called genotypes are
    needed, uncertainty in low-coverage samples is propagated into the covariance
    matrix instead of being thrown away by a hard genotype call.

    pcangsd.log is a real output file PCAngsd itself writes next to the .cov
    matrix (its own run summary: arguments, sites loaded, convergence), tracked
    here as output rather than folded into the Snakemake log:, same as ANGSD's
    .arg files are tracked as output in 4_relatedness.smk.
    """
    input:
        beagle="results/structure/pca/all.{ref_name}." + PRUNE_TAG + "." + PCA_TAG + ".beagle.gz",
    output:
        cov="results/structure/pca/all.{ref_name}." + PRUNE_TAG + "." + PCA_TAG + ".cov",
        pcalog="results/structure/pca/all.{ref_name}." + PRUNE_TAG + "." + PCA_TAG + ".log",
    params:
        prefix=lambda wildcards, output: output.cov[: -len(".cov")],
        extra=PCA_EXTRA,
    log:
        "logs/structure/pca_pcangsd_{ref_name}.log",
    benchmark:
        "benchmarks/structure/pca_pcangsd_{ref_name}.benchmark"
    conda:
        "../envs/pcangsd.yaml"
    threads: PCA_THREADS
    shell:
        """
        pcangsd -b {input.beagle} -o {params.prefix} -t {threads} {params.extra} &> {log}
        """
