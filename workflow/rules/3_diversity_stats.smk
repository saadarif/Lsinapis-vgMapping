# ==============================================================================
# DIVERSITY STATISTICS
#   (i)  Per-individual heterozygosity from the bcftools stats files produced by
#        the genotyping workflows (2a_call_genotypes_noTrans / 2b_call_genotypes_rescaled)
#   (ii) Windowed pi, dxy and Fst with pixy
#
# The rules are generic: {dataset} is either "notrans" or "rescaled" (which is
# also the tag carried in each workflow's file names), and {prefix} already
# encodes the call type / site type / missingness of the input file. Which of
# those combinations actually get built is decided by the target list in the
# Snakefile (section 8), driven by params: run_diversity_stats in the config.
# ==============================================================================
DIV_PARAMS = config.get("params", {}).get("run_diversity_stats", {})

PIXY_FST_TYPE = DIV_PARAMS.get("fst_type", "hudson")
PIXY_CORES = DIV_PARAMS.get("n_cores", 4)
PIXY_EXTRA = DIV_PARAMS.get("extra", "")


def get_poplist_for_dataset(wildcards):
    """The per-workflow poplist, already filtered down to the samples that were
    actually genotyped (i.e. with `drop_samples` removed)."""
    return FILTERED_POPLIST_NT if wildcards.dataset == "notrans" else FILTERED_POPLIST_RS


# ==============================================================================
# (i) INDIVIDUAL HETEROZYGOSITY
# ==============================================================================
rule individual_heterozygosity:
    """
    Per-sample heterozygosity from the PSC block of a bcftools stats file.
      hetperbp  = nHets / (nRefHom + nNonRefHom + nHets)  -- het sites per called site
      hetperalt = nHets / (2*nNonRefHom + nHets)          -- het sites per alternate allele
    """
    input:
        stats = "results/genotyping_{dataset}/{prefix}.bcf.stats",
    output:
        het = "results/diversity_stats/{dataset}/heterozygosity/{prefix}.bcf.stats.het",
    wildcard_constraints:
        dataset = "notrans|rescaled",
        prefix = "[^/]+",
    shell:
        """
        echo "sample\thetperbp\thetperalt" > {output.het}
        grep PSC {input.stats} | grep -v "#" | \
            awk '{{print $3"\t"$6/($6+$5+$4)"\t"$6/(2*$5+$6)}}' >> {output.het}
        """


# ==============================================================================
# (ii) PI, DXY AND FST WITH PIXY
# ==============================================================================
rule pixy_pi_dxy_fst:
    """
    Windowed pi (within population), dxy (between populations) and Fst with pixy.
    Requires a bgzipped + tabix indexed VCF that retains invariant sites, hence
    the allsites VCFs rather than the biallelic ones.
    """
    input:
        vcf = "results/genotyping_{dataset}/{prefix}.vcf.gz",
        tbi = "results/genotyping_{dataset}/{prefix}.vcf.gz.tbi",
        pops = get_poplist_for_dataset,
    output:
        pi = "results/diversity_stats/{dataset}/pixy/{prefix}.w{winsize}/pixy_pi.txt",
        dxy = "results/diversity_stats/{dataset}/pixy/{prefix}.w{winsize}/pixy_dxy.txt",
        fst = "results/diversity_stats/{dataset}/pixy/{prefix}.w{winsize}/pixy_fst.txt",
    wildcard_constraints:
        dataset = "notrans|rescaled",
        prefix = "[^/]+",
        winsize = r"\d+",
    params:
        fold = lambda wildcards, output: os.path.dirname(output.pi),
        fst_type = PIXY_FST_TYPE,
        extra = PIXY_EXTRA,
    log:
        "logs/diversity_stats/{dataset}/pixy_{prefix}.w{winsize}.log",
    benchmark:
        "benchmarks/diversity_stats/{dataset}/pixy_{prefix}.w{winsize}.benchmark"
    conda:
        "../envs/pixy.yaml"
    threads: PIXY_CORES
    shell:
        """
        pixy --stats pi fst dxy --vcf {input.vcf} --populations {input.pops} \
            --window_size {wildcards.winsize} --fst_type {params.fst_type} \
            --n_cores {threads} {params.extra} \
            --output_folder {params.fold} > {log} 2>&1
        """
