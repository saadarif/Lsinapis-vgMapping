# ==============================================================================
# POPULATION STRUCTURE: PCA FROM GENOTYPE LIKELIHOODS (PCAngsd)
#   (i)  Drop any individuals listed under params: run_pca: exclude_samples from
#        the LD-pruned, unlinked-SNP beagle file built in 5_LD_pruning.smk, then
#        every remaining individual of any population listed under params:
#        run_pca: exclude_populations
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
#
# Whole populations can be left out as well (params: run_pca: exclude_populations),
# e.g. to look at structure within a subset of populations without the most
# divergent ones dominating the first components. Population membership comes
# from params: run_pca: poplist.
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

# Populations to leave out entirely, applied after the individual exclusions
# above: every sample still left that belongs to one of them is dropped too.
# Empty (the default) keeps all populations.
PCA_EXCLUDE_POPS = PCA_PARAMS.get("exclude_populations", []) or []
if isinstance(PCA_EXCLUDE_POPS, str):
    PCA_EXCLUDE_POPS = [PCA_EXCLUDE_POPS]
PCA_EXCLUDE_POPS = sorted(set(PCA_EXCLUDE_POPS))
# Two columns, sample_id and population, no header
PCA_POPFILE = PCA_PARAMS.get("poplist", "Pops.txt")

PCA_EXCLUDE_POP_SAMPLES = []
if PCA_EXCLUDE_POPS:
    if not os.path.exists(PCA_POPFILE):
        print(
            f"\nERROR: params: run_pca: exclude_populations is set, but the poplist "
            f"'{PCA_POPFILE}' (params: run_pca: poplist) does not exist.\n"
        )
        sys.exit(1)
    _pca_pop_of = {}
    with open(PCA_POPFILE) as _fh:
        for _line in _fh:
            _fields = _line.split()
            if len(_fields) >= 2:
                _pca_pop_of[_fields[0]] = _fields[1]
    # A misspelt population would otherwise silently drop nobody
    _unknown_pops = [p for p in PCA_EXCLUDE_POPS if p not in set(_pca_pop_of.values())]
    if _unknown_pops:
        print(
            f"\nERROR: params: run_pca: exclude_populations lists "
            f"{', '.join(_unknown_pops)}, not found in {PCA_POPFILE}. Populations "
            f"there: {', '.join(sorted(set(_pca_pop_of.values())))}.\n"
        )
        sys.exit(1)
    _pca_remaining = [s for s in KEEP_SAMPLES_REL if s not in PCA_EXCLUDE_SAMPLES]
    _no_pop = [s for s in _pca_remaining if s not in _pca_pop_of]
    if _no_pop:
        print(
            f"\nERROR: no population in {PCA_POPFILE} for: {', '.join(_no_pop)}. "
            "Needed to apply params: run_pca: exclude_populations.\n"
        )
        sys.exit(1)
    PCA_EXCLUDE_POP_SAMPLES = [s for s in _pca_remaining if _pca_pop_of[s] in PCA_EXCLUDE_POPS]

# Everything dropped from the beagle file, and what is left in it, the latter in
# beagle column order (reused by 6b_Structure_NGadmix.smk, which runs on the
# same file).
PCA_DROP_SAMPLES = PCA_EXCLUDE_SAMPLES + PCA_EXCLUDE_POP_SAMPLES
PCA_KEEP_SAMPLES = [s for s in KEEP_SAMPLES_REL if s not in PCA_DROP_SAMPLES]
if not PCA_KEEP_SAMPLES:
    print(
        "\nERROR: params: run_pca: exclude_samples / exclude_populations leave no "
        "samples for the PCA.\n"
    )
    sys.exit(1)

# Excluded populations are named in the tag rather than their samples. The
# .exclpop- part is only added when a population is excluded, so file names are
# unchanged from before this option existed when it is left empty.
PCA_TAG = "excl-" + ("-".join(PCA_EXCLUDE_SAMPLES) if PCA_EXCLUDE_SAMPLES else "none")
if PCA_EXCLUDE_POPS:
    PCA_TAG += ".exclpop-" + "-".join(PCA_EXCLUDE_POPS)
PCA_PREFIX = f"results/structure/pca/all.{REF_NAME}.{PRUNE_TAG}.{PCA_TAG}"

# 1-based column numbers to drop from the beagle file for every dropped sample:
# marker=1, allele1=2, allele2=3, then 3 genotype-likelihood columns per
# individual starting at column 4, in KEEP_SAMPLES_REL order.
PCA_EXCLUDE_COLS = []
for _sid in PCA_DROP_SAMPLES:
    _start = KEEP_SAMPLES_REL.index(_sid) * 3 + 4
    PCA_EXCLUDE_COLS.extend([_start, _start + 1, _start + 2])
PCA_EXCLUDE_COLS_STR = ",".join(str(c) for c in sorted(PCA_EXCLUDE_COLS))

# Requested from the Snakefile target list (section 11)
PCA_TARGETS = [
    f"{PCA_PREFIX}.cov",
]


rule pca_exclude_beagle:
    """
    Drops params: run_pca: exclude_samples, and every individual of params:
    run_pca: exclude_populations, from the LD-pruned beagle file before PCA. If
    nothing is excluded, this is just a pass-through symlink.
    """
    input:
        beagle="results/ld/all.{ref_name}." + PRUNE_TAG + ".pruned.beagle.gz",
    output:
        beagle="results/structure/pca/all.{ref_name}." + PRUNE_TAG + "." + PCA_TAG + ".beagle.gz",
    params:
        cols=PCA_EXCLUDE_COLS_STR,
    log:
        "logs/structure/pca_exclude_beagle_{ref_name}." + PRUNE_TAG + "." + PCA_TAG + ".log",
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
        "logs/structure/pca_pcangsd_{ref_name}." + PRUNE_TAG + "." + PCA_TAG + ".log",
    benchmark:
        "benchmarks/structure/pca_pcangsd_{ref_name}." + PRUNE_TAG + "." + PCA_TAG + ".benchmark"
    conda:
        "../envs/pcangsd.yaml"
    threads: PCA_THREADS
    shell:
        """
        pcangsd -b {input.beagle} -o {params.prefix} -t {threads} {params.extra} &> {log}
        """
