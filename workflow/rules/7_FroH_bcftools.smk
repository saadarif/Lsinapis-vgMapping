# ==============================================================================
# RUNS OF HOMOZYGOSITY (bcftools roh)
#   Calls runs of homozygosity (RoH) per individual from the called genotypes of
#   the genotyping workflows (2a_call_genotypes_noTrans / 2b_call_genotypes_rescaled).
#
# Which BCFs RoH is called on is set under params: run_roh in the config:
#   datasets   -> "notrans", "rescaled", or both
#   call_types -> "indCall", "jointCall", or both
# and, for each of those, RoH is called at every missingness threshold under
# params: run_genotyping: maxMissing.
#
# Follows the bcftools roh settings of Nolen et al.
# (zjnolen/polyommatini-temporal-genomics, bcftools_roh.smk): the multi-sample
# BCF is the input, but each individual is analysed on its own -- a fixed
# default allele frequency (--AF-dflt) is used instead of one estimated from the
# samples in the file, which would otherwise depend on which populations and how
# many individuals happen to be in it, and hom-ref genotypes are ignored
# (--ignore-homref) so only sites where the individual carries an alternate
# allele inform the HMM. That is also why the biallelic (variant sites only)
# BCFs are used rather than the allsites ones.
# ==============================================================================
ROH_PARAMS = config.get("params", {}).get("run_roh", {})

ROH_DATASETS = ROH_PARAMS.get("datasets", ["notrans"])
ROH_CALL_TYPES = ROH_PARAMS.get("call_types", ["jointCall"])
ROH_REC_RATE = ROH_PARAMS.get("rec_rate", "1e-8") # bcftools roh -M (per bp)
ROH_EXTRA = ROH_PARAMS.get("extra", "")
ROH_THREADS = ROH_PARAMS.get("threads", 4)

ROH_SITE_TYPE = "biallelic"

ROH_TAG = f"roh.M{ROH_REC_RATE}"

# Every missingness threshold produced by the genotyping workflows
ROH_MISSING_VALS = config["params"]["run_genotyping"].get("maxMissing", [0.0])
if not isinstance(ROH_MISSING_VALS, list):
    ROH_MISSING_VALS = [ROH_MISSING_VALS]

# Per genotyping workflow: its config toggle and the BCF name stem up to the
# call type, built from that workflow's own constants (2a / 2b rule files).
ROH_DATASET_INFO = {
    "notrans": {
        "toggle": "run_genotyping_notrans",
        "stem": f"merged.all.{REF_NAME}.sitefilt.bQ{BASEQ_NT}.mq{MAPQ_NT}.snps5.noIndel.Q30.dp{MIN_DP_NT}-{MAX_DP_NT}.AB",
    },
    "rescaled": {
        "toggle": "run_genotyping_rescaled",
        "stem": f"merged.all.{REF_NAME}.sitefilt.bQ{BASEQ_RS}.mq{MAPQ_RS}.snps5.noIndel.Q30.dp{MIN_DP_RS}-{MAX_DP_RS}.AB",
    },
}

# Requested from the Snakefile target list (section 13). A dataset is only
# included if its genotyping workflow actually ran, otherwise the input BCF
# would not exist.
ROH_TARGETS = []
if config.get("run_roh", False):
    for _dataset in ROH_DATASETS:
        if _dataset not in ROH_DATASET_INFO:
            print(f"\nERROR: unknown 'datasets' entry '{_dataset}' under params: run_roh. Use 'notrans' and/or 'rescaled'.\n")
            sys.exit(1)
        if not config.get(ROH_DATASET_INFO[_dataset]["toggle"], False):
            continue
        for _c_type in ROH_CALL_TYPES:
            if _c_type not in ("indCall", "jointCall"):
                print(f"\nERROR: unknown 'call_types' entry '{_c_type}' under params: run_roh. Use 'indCall' and/or 'jointCall'.\n")
                sys.exit(1)
            for _m_val in ROH_MISSING_VALS:
                ROH_TARGETS.append(
                    f"results/roh/bcftools/{_dataset}/{ROH_DATASET_INFO[_dataset]['stem']}.{_c_type}.{_dataset}.{ROH_SITE_TYPE}.fmiss{_m_val}.{ROH_TAG}.regs.txt"
                )


rule bcftools_roh:
    """
    Runs of homozygosity per individual with bcftools roh.

      -G             use called genotypes (GT) rather than genotype likelihoods,
                     with this phred-scaled quality assumed for every genotype
      --ignore-homref  skip hom-ref genotypes
      --AF-dflt      allele frequency assumed at every site
      -M             constant recombination rate per bp
      -O r           write only the RoH regions (RG lines: sample, chromosome,
                     start, end, length in bp, number of markers, quality), not
                     the per-site HMM state of every individual

    The rule is generic: {dataset} is "notrans" or "rescaled" and {prefix}
    already encodes the call type / site type / missingness of the input BCF.
    """
    input:
        bcf = "results/genotyping_{dataset}/{prefix}.bcf",
        csi = "results/genotyping_{dataset}/{prefix}.bcf.csi",
    output:
        regs = "results/roh/bcftools/{dataset}/{prefix}." + ROH_TAG + ".regs.txt",
    wildcard_constraints:
        dataset = "notrans|rescaled",
        prefix = "[^/]+",
    params:
        rec_rate = ROH_REC_RATE,
        extra = ROH_EXTRA,
    log:
        "logs/roh/bcftools/{dataset}/bcftools_roh_{prefix}.log",
    benchmark:
        "benchmarks/roh/bcftools/{dataset}/bcftools_roh_{prefix}.benchmark"
    conda:
        "../envs/bcftools121.yaml"
    threads: ROH_THREADS
    shell:
        """
        bcftools roh -G30 --ignore-homref -M {params.rec_rate} \
            --AF-dflt 0.4 --threads {threads} {params.extra} \
            -O r -o {output.regs} {input.bcf} 2> {log}
        """
