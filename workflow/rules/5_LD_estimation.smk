# ==============================================================================
# LINKAGE DISEQUILIBRIUM ESTIMATION AND PRUNING
#   (i)   Pairwise LD between SNPs with ngsLD, from the same all-samples beagle
#         genotype likelihoods built in 4_relatedness.smk
#   (ii)  Prune down to a set of SNPs in approximate linkage equilibrium with
#         prune_graph
#   (iii) Subset the beagle file to just those unlinked SNPs
#
# ngsLD (fgvieira/ngsLD) and prune_graph (fgvieira/prune_graph) have no
# bioconda/conda-forge package, so the two rules that need them use `container:`
# with the image PopGLen builds for exactly this gap, rather than `conda:`.
# `apptainer` needs to be added to profiles/default/config.yaml:
# software-deployment-method for this to run (kept conda for every other rule).
#
# The pruned beagle file feeds population structure analyses downstream
# (6a_Strucure_angsdPCA.smk, and NGSadmix later).
# ==============================================================================
LD_PARAMS = config.get("params", {}).get("run_ld_estimation", {})

LD_MAX_KB = LD_PARAMS.get("max_kb_dist", 4000)  # ngsLD --max_kb_dist (kb)
# prune_graph's own distance cutoff, in BASES rather than kb, and independent of
# LD_MAX_KB above: ngsLD computes LD out to LD_MAX_KB, but prune_graph can be
# told to only treat the closer subset of those pairs as "linked" edges, e.g. to
# compute LD over a wide window for decay analysis while pruning more tightly.
# Only has an effect while <= LD_MAX_KB * 1000; prune_graph never sees a pair
# further apart than that, since ngsLD didn't compute LD for it in the first place.
LD_PRUNE_MAX_BP = LD_PARAMS.get("prune_max_dist_bp", 50000)
if LD_PRUNE_MAX_BP > LD_MAX_KB * 1000:
    print(
        f"\nWARNING: params: run_ld_estimation: prune_max_dist_bp ({LD_PRUNE_MAX_BP} bp) "
        f"is larger than max_kb_dist ({LD_MAX_KB} kb = {LD_MAX_KB * 1000} bp); ngsLD never "
        "computes LD beyond max_kb_dist, so prune_max_dist_bp has no effect past that point.\n"
    )
LD_MIN_R2 = LD_PARAMS.get("min_r2", 0.1)  # prune_graph edge threshold
LD_THREADS = LD_PARAMS.get("threads", 10)

NGSLD_CONTAINER = "docker://ghcr.io/zjnolen/ngsld:1.2.0"

LD_TAG = f"{REL_TAG}.maxkb{LD_MAX_KB}"
LD_PREFIX = f"results/ld/all.{REF_NAME}.{LD_TAG}"
PRUNE_TAG = f"{LD_TAG}.prunebp{LD_PRUNE_MAX_BP}.minr2{LD_MIN_R2}"
PRUNE_PREFIX = f"results/ld/all.{REF_NAME}.{PRUNE_TAG}"

# Requested from the Snakefile target list (section 10)
LD_TARGETS = [
    f"{LD_PREFIX}.ld.gz",
    f"{PRUNE_PREFIX}.unlinked.pos",
    f"{PRUNE_PREFIX}.pruned.beagle.gz",
]


# ==============================================================================
# (i) PAIRWISE LD WITH ngsLD
# ==============================================================================
rule ngsld_estimate:
    """
    Pairwise LD between all SNP pairs within max_kb_dist of each other, computed
    straight from the all-samples beagle file the relatedness workflow builds
    (4_relatedness.smk, angsd_beagle_relatedness).

    The .pos file ngsLD needs (chr, 1-based pos, no header, same row order as the
    beagle) is derived from the beagle's own marker column: ANGSD names each
    marker "{chr}_{pos}", so splitting on the *last* underscore recovers both,
    which is safe even if a chromosome name itself contains an underscore since
    ANGSD never puts one in the position half.

    ngsLD's own --out only writes plain text; a 20 kb test window alone produced
    an 11 million line .ld file, so the output is piped through gzip here rather
    than using --out directly. prune_graph (ngsld_prune, below) reads it back via
    zcat, so nothing downstream needs it uncompressed.
    """
    input:
        beagle="results/relatedness/all.{ref_name}." + REL_TAG + ".beagle.gz",
    output:
        pos="results/ld/all.{ref_name}." + LD_TAG + ".pos",
        ld="results/ld/all.{ref_name}." + LD_TAG + ".ld.gz",
    params:
        nind=N_IND_REL,
        max_kb=LD_MAX_KB,
    log:
        "logs/ld/ngsld_estimate_{ref_name}.log",
    benchmark:
        "benchmarks/ld/ngsld_estimate_{ref_name}.benchmark"
    container:
        NGSLD_CONTAINER
    threads: LD_THREADS
    shell:
        """
        (zcat {input.beagle} | awk '{{print $1}}' | sed 's/\\(.*\\)_/\\1\\t/' \
            | tail -n +2 > {output.pos}

        nsites=$(wc -l < {output.pos})
        echo "Running ngsLD on $nsites sites for {params.nind} individuals"

        ngsLD --geno {input.beagle} --n_ind {params.nind} --n_sites $nsites \
            --pos {output.pos} --probs --n_threads {threads} \
            --max_kb_dist {params.max_kb} \
            | gzip > {output.ld}) &> {log}
        """


# ==============================================================================
# (ii) PRUNE TO UNLINKED SNPs WITH prune_graph
# ==============================================================================
rule ngsld_prune:
    """
    Prunes SNPs down to a set in approximate linkage equilibrium: SNPs within
    prune_max_dist_bp (bases, see LD_PRUNE_MAX_BP above -- independent of
    ngsld_estimate's max_kb_dist) with r2 >= min_r2 are treated as an edge
    (column_3 = Dist in bp, column_7 = r2 in the .ld file), and prune_graph
    iteratively drops the most connected node until no edges remain, leaving one
    representative SNP per linked cluster. SNPs that never appear in an edge
    (i.e. never within prune_max_dist_bp of another SNP that also passed
    max_kb_dist filtering when the .ld file was built) are, as a side effect,
    not written to the output either; this matches how prune_graph is used
    upstream in PopGLen.

    -out is one surviving node per line, formatted "chr:pos" (the same format
    ngsLD used for Pos1/Pos2), which is exactly the "positions to keep" format
    workflow/scripts/prune_beagle.py expects.
    """
    input:
        ld="results/ld/all.{ref_name}." + LD_TAG + ".ld.gz",
    output:
        pos="results/ld/all.{ref_name}." + PRUNE_TAG + ".unlinked.pos",
    params:
        max_bp=LD_PRUNE_MAX_BP,
        min_r2=LD_MIN_R2,
    log:
        "logs/ld/ngsld_prune_{ref_name}.log",
    benchmark:
        "benchmarks/ld/ngsld_prune_{ref_name}.benchmark"
    container:
        NGSLD_CONTAINER
    threads: LD_THREADS
    shell:
        """
        (zcat {input.ld} | prune_graph --n-threads {threads} \
            --weight-field column_7 \
            --weight-filter "column_3 <= {params.max_bp} && column_7 >= {params.min_r2}" \
            --out {output.pos}) &> {log}
        """


# ==============================================================================
# (iii) SUBSET THE BEAGLE FILE TO THE UNLINKED SNPs
# ==============================================================================
rule prune_beagle:
    """
    Subsets the all-samples beagle file down to the unlinked SNPs from
    ngsld_prune. See workflow/scripts/prune_beagle.py.
    """
    input:
        beagle="results/relatedness/all.{ref_name}." + REL_TAG + ".beagle.gz",
        pos="results/ld/all.{ref_name}." + PRUNE_TAG + ".unlinked.pos",
    output:
        beagle="results/ld/all.{ref_name}." + PRUNE_TAG + ".pruned.beagle.gz",
    log:
        "logs/ld/prune_beagle_{ref_name}.log",
    benchmark:
        "benchmarks/ld/prune_beagle_{ref_name}.benchmark"
    conda:
        "../envs/python.yaml"
    threads: 2
    script:
        "../scripts/prune_beagle.py"
