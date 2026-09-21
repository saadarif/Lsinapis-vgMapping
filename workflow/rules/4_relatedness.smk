# ==============================================================================
# RELATEDNESS FROM GENOTYPE LIKELIHOODS
#   (i)  A single beagle genotype likelihood file holding all samples, built
#        with ANGSD in one pass (GATK model, -doGlf 2)
#   (ii) R0, R1 and KING-robust kinship for every pair of samples with
#        ngsRelate, using the SFS based IBSrelate estimator of
#        Waples, Albrechtsen & Nielsen (2019) Mol Ecol 28:35-48
#
# ANGSD is restricted to the same sites-filter the genotyping workflows call over,
# i.e. config: site_filter_bed, which 2a/2b hand to `bcftools mpileup -R`.
# The BED is converted to an ANGSD sites file first: ANGSD reads a three column
# sites file as 1-based inclusive "chr start end", so the BED start needs a +1
# to cover exactly the same bases (this also keeps BED starts of 0 legal, ANGSD
# rejects a position of 0).
#
# Everything is produced for one sample set in one go, so there are no
# {call_type}/{site_type} wildcards here: the sample set and the filters come
# from params: run_relatedness in the config and are baked into the file names.
# ==============================================================================
import math

REL_PARAMS = config.get("params", {}).get("run_relatedness", {})

BASEQ_REL = config["baseQ"]
MAPQ_REL = config["mapQ"]
SITEFILTER_BED_REL = config["site_filter_bed"]

# Which BAM set the likelihoods are built from:
#   "dedup"    -> (default) historical deduplicated (unrescaled) + modern
#                 clipOverlap, the BAM set the no-transitions genotyping uses.
#                 Pair this with rm_trans: 1 so damaged sites are dropped
#                 instead of down-weighted.
#   "rescaled" -> historical mapDamage rescaled + modern clipOverlap, i.e. the
#                 final BAMs, with damage down-weighted in the base qualities
#                 instead of filtered out. Opt in explicitly with bam_stage:
#                 "rescaled" when that's what's wanted.
REL_BAM_STAGE = REL_PARAMS.get("bam_stage", "dedup")
if REL_BAM_STAGE not in ("rescaled", "dedup"):
    print(
        f"\nERROR: unknown 'bam_stage' entry '{REL_BAM_STAGE}' under "
        "params: run_relatedness. Use 'rescaled' or 'dedup'.\n"
    )
    sys.exit(1)

REL_GL_MODEL = REL_PARAMS.get("gl_model", 2)          # 2 = GATK model
REL_MIN_DP_IND = REL_PARAMS.get("minDP_ind", 3)       # -setMinDepthInd
REL_MAX_DP_IND = REL_PARAMS.get("maxDP_ind", 30)      # -setMaxDepthInd
REL_MIN_IND_FRAC = REL_PARAMS.get("min_ind_frac", 0.9)
REL_MIN_MAF = REL_PARAMS.get("min_maf", 0.05)
REL_SNP_PVAL = str(REL_PARAMS.get("snp_pval", "1e-6"))
REL_RM_TRANS = REL_PARAMS.get("rm_trans", 1)
REL_EXTRA = REL_PARAMS.get("extra", "")
REL_THREADS = REL_PARAMS.get("threads", 8)

# Its own drop_samples, independent of params: run_genotyping. Relatedness is
# often the analysis used to decide who to drop in the first place, so nobody
# is dropped by default here; set params: run_relatedness: drop_samples to
# exclude specific samples (e.g. known duplicates/outgroups) if needed.
DROP_SAMPLES_REL = REL_PARAMS.get("drop_samples", [])
KEEP_SAMPLES_REL = (
    samples_df[~samples_df["sample_id"].isin(DROP_SAMPLES_REL)]["sample_id"]
    .unique()
    .tolist()
)
N_IND_REL = len(KEEP_SAMPLES_REL)
# "present in 90% of the individuals" -> -minInd, counted after the per
# individual depth filter has blanked the individuals that fall outside
# setMinDepthInd/setMaxDepthInd at that site
MIN_IND_REL = math.ceil(REL_MIN_IND_FRAC * N_IND_REL)

REL_TAG = (
    f"{REL_BAM_STAGE}.sitefilt.bQ{BASEQ_REL}.mq{MAPQ_REL}"
    f".dpInd{REL_MIN_DP_IND}-{REL_MAX_DP_IND}.minInd{REL_MIN_IND_FRAC}"
    f".maf{REL_MIN_MAF}.pv{REL_SNP_PVAL}"
)
REL_PREFIX = f"results/relatedness/all.{REF_NAME}.{REL_TAG}"

# Requested from the Snakefile target list (section 9)
RELATEDNESS_TARGETS = [
    f"{REL_PREFIX}.beagle.gz",
    f"{REL_PREFIX}.ngsrelate.tsv",
    f"{REL_PREFIX}.ibsrelate.R0R1KING.tsv",
]


# ==============================================================================
# INPUT HELPER FUNCTIONS
# ==============================================================================
def get_bam_for_relatedness(sid):
    """BAM handed to ANGSD for one sample, following params: run_relatedness:
    bam_stage. Samples listed in subsample_samples use their subsampled BAM, as
    they do in the genotyping workflows."""
    src = samples_df[samples_df["sample_id"] == sid]["source"].iloc[0]

    if src == "modern":
        stage_str = ".clipped"
    else:
        stage_str = ".rescaled" if REL_BAM_STAGE == "rescaled" else ""

    subsample_samples = config.get("subsample_samples", [])
    target_dp = config.get("subsample_depth", 15)

    if sid in subsample_samples:
        return (
            f"results/mapping/{src}/{sid}.{REF_NAME}.merged.dedup.merged"
            f"{stage_str}.subs{target_dp}.q{MAPQ_REL}.bam"
        )
    return f"results/mapping/{src}/{sid}.{REF_NAME}.merged.dedup.merged{stage_str}.bam"


def get_all_bams_for_relatedness(wildcards):
    return [get_bam_for_relatedness(sid) for sid in KEEP_SAMPLES_REL]


def get_all_bais_for_relatedness(wildcards):
    return [get_bam_for_relatedness(sid) + ".bai" for sid in KEEP_SAMPLES_REL]


# ==============================================================================
# (i) GENOTYPE LIKELIHOODS FOR ALL SAMPLES (ANGSD BEAGLE FILE)
# ==============================================================================
rule angsd_reference:
    """
    Stage the reference for ANGSD. ANGSD refuses a reference whose .fa.fai is
    not strictly newer than the .fa and, unhelpfully, quits with exit status 0
    and no error before writing anything (aio::isNewer, called from
    abcGetFasta.cpp). The mtimes of the project reference currently trip that
    check, so ANGSD gets a symlink to the fasta next to a freshly written copy
    of the index: the shared reference directory is left untouched and only the
    few kB of the .fai are duplicated.

    Not relatedness-specific: any future ANGSD-based rule (PCA, admixture, ...)
    hits the same isNewer check, so this lives directly under results/ rather
    than under results/relatedness/ to be reused rather than duplicated.
    """
    input:
        fa=config["reference"] + ".fa",
        fai=config["reference"] + ".fa.fai",
    output:
        fa="results/angsd_ref/{ref_name}.fa",
        fai="results/angsd_ref/{ref_name}.fa.fai",
    log:
        "logs/angsd_ref/{ref_name}.log",
    shell:
        """
        (ln -sf "$(realpath {input.fa})" {output.fa}

        cp -f {input.fai} {output.fai}
        chmod u+w {output.fai}
        touch {output.fai}) &> {log}
        """


rule angsd_sites_from_bed:
    """
    Turn the genotyping site filter BED into an indexed ANGSD sites file, plus
    the matching -rf chromosome list so ANGSD does not walk scaffolds that the
    filter excludes entirely. `$2+1` converts the 0-based half-open BED interval
    into the 1-based inclusive interval ANGSD expects, keeping the sites
    identical to those passed to `bcftools mpileup -R` in 2a/2b.

    Not relatedness-specific: any future ANGSD-based rule (PCA, admixture, ...)
    that should run over the same site_filter_bed can depend on this directly,
    so it lives under results/ rather than under results/relatedness/.
    """
    input:
        bed=SITEFILTER_BED_REL,
    output:
        sites="results/angsd_sites/{ref_name}.sitefilt.sites",
        bin="results/angsd_sites/{ref_name}.sitefilt.sites.bin",
        idx="results/angsd_sites/{ref_name}.sitefilt.sites.idx",
        regions="results/angsd_sites/{ref_name}.sitefilt.rf",
    log:
        "logs/angsd_sites/{ref_name}.log",
    conda:
        "../envs/angsd.yaml"
    shell:
        """
        (awk -F '\\t' -v OFS='\\t' '{{print $1, $2+1, $3}}' {input.bed} | \
            sort -k1,1V -k2,2n > {output.sites}

        cut -f1 {output.sites} | uniq > {output.regions}

        angsd sites index {output.sites}

        # ANGSD warns and asks for a reindex when the binary index looks older
        # than the sites file, and both are written within the same second here
        touch {output.bin} {output.idx}) &> {log}
        """


rule angsd_bamlist_relatedness:
    """
    BAM list for ANGSD and the matching sample IDs for ngsRelate. Both are
    written in the same order so that individual i in the beagle file is line i
    of the ID file.
    """
    input:
        bams=get_all_bams_for_relatedness,
        bais=get_all_bais_for_relatedness,
    output:
        bamlist="results/relatedness/all.{ref_name}.bamlist",
        ids="results/relatedness/all.{ref_name}.ids",
    run:
        with open(output.bamlist, "w") as fh:
            for sid in KEEP_SAMPLES_REL:
                fh.write(get_bam_for_relatedness(sid) + "\n")
        with open(output.ids, "w") as fh:
            for sid in KEEP_SAMPLES_REL:
                fh.write(sid + "\n")


rule angsd_beagle_relatedness:
    """
    Genotype likelihoods for all samples in a single ANGSD run, written as a
    beagle file (-doGlf 2). Major/minor and allele frequencies are inferred from
    the data (-doMajorMinor 1 -doMaf 1) so that SNPs can be called with
    -SNP_pval / -minMaf.

    -doCounts 1 is what makes -setMinDepthInd / -setMaxDepthInd work: ANGSD
    zeroes the reads of any individual outside that depth range at that site
    before the likelihoods are computed, and -minInd is then applied to the
    individuals that still carry data.

    The reference comes from angsd_reference rather than straight from
    config: reference, see that rule for why.
    """
    input:
        bamlist="results/relatedness/all.{ref_name}.bamlist",
        bams=get_all_bams_for_relatedness,
        bais=get_all_bais_for_relatedness,
        ref="results/angsd_ref/{ref_name}.fa",
        fai="results/angsd_ref/{ref_name}.fa.fai",
        sites="results/angsd_sites/{ref_name}.sitefilt.sites",
        sites_bin="results/angsd_sites/{ref_name}.sitefilt.sites.bin",
        sites_idx="results/angsd_sites/{ref_name}.sitefilt.sites.idx",
        regions="results/angsd_sites/{ref_name}.sitefilt.rf",
    output:
        beagle="results/relatedness/all.{ref_name}." + REL_TAG + ".beagle.gz",
        mafs="results/relatedness/all.{ref_name}." + REL_TAG + ".mafs.gz",
        arg="results/relatedness/all.{ref_name}." + REL_TAG + ".arg",
    params:
        out=lambda wildcards, output: output.arg[: -len(".arg")],
        gl_model=REL_GL_MODEL,
        mapq=MAPQ_REL,
        baseq=BASEQ_REL,
        min_dp_ind=REL_MIN_DP_IND,
        max_dp_ind=REL_MAX_DP_IND,
        min_ind=MIN_IND_REL,
        min_maf=REL_MIN_MAF,
        snp_pval=REL_SNP_PVAL,
        rm_trans=REL_RM_TRANS,
        extra=REL_EXTRA,
    log:
        "logs/relatedness/angsd_beagle_{ref_name}.log",
    benchmark:
        "benchmarks/relatedness/angsd_beagle_{ref_name}.benchmark"
    conda:
        "../envs/angsd.yaml"
    threads: REL_THREADS
    shell:
        """
        angsd -bam {input.bamlist} -ref {input.ref} \
            -rf {input.regions} -sites {input.sites} \
            -GL {params.gl_model} -doGlf 2 -doMajorMinor 1 -doMaf 1 -doCounts 1 \
            -minMapQ {params.mapq} -minQ {params.baseq} \
            -setMinDepthInd {params.min_dp_ind} -setMaxDepthInd {params.max_dp_ind} \
            -minInd {params.min_ind} \
            -SNP_pval {params.snp_pval} -minMaf {params.min_maf} \
            -rmTrans {params.rm_trans} {params.extra} \
            -nThreads {threads} -out {params.out} &> {log}
        """


# ==============================================================================
# (ii) PAIRWISE RELATEDNESS WITH ngsRelate (IBSrelate, SFS BASED)
# ==============================================================================
rule ngsrelate_ibsrelate:
    """
    Pairwise relatedness for every sample pair from the beagle file. No allele
    frequency file is supplied (-L instead of -f), so R0, R1 and KING-robust
    kinship come from the two-dimensional SFS of each pair, which is the
    allele-frequency-free IBSrelate estimator of Waples et al. (2019).

    Two files are produced: the full ngsRelate table (the nine Jacquard
    coefficients, theta, the 2dsfs, ...) and the plain R0/R1/KING table. The
    columns are picked out by header name rather than by position, since the
    header gains the ida/idb columns only when -z is given.
    """
    input:
        beagle="results/relatedness/all.{ref_name}." + REL_TAG + ".beagle.gz",
        ids="results/relatedness/all.{ref_name}.ids",
    output:
        res="results/relatedness/all.{ref_name}." + REL_TAG + ".ngsrelate.tsv",
        tbl="results/relatedness/all.{ref_name}." + REL_TAG + ".ibsrelate.R0R1KING.tsv",
    params:
        nind=N_IND_REL,
    log:
        "logs/relatedness/ngsrelate_{ref_name}.log",
    benchmark:
        "benchmarks/relatedness/ngsrelate_{ref_name}.benchmark"
    conda:
        "../envs/ngsrelate.yaml"
    threads: REL_THREADS
    shell:
        """
        (nsites=$(zcat {input.beagle} | tail -n +2 | wc -l)
        echo "Running ngsRelate on $nsites sites for {params.nind} individuals"

        ngsRelate -G {input.beagle} -n {params.nind} -L $nsites \
            -z {input.ids} -p {threads} -O {output.res}

        awk -F '\\t' -v OFS='\\t' '
            NR == 1 {{
                for (i = 1; i <= NF; i++) col[$i] = i
                print "ind1", "ind2", "nSites", "R0", "R1", "KING"
                next
            }}
            {{
                print $col["ida"], $col["idb"], $col["nSites"], \
                      $col["R0"], $col["R1"], $col["KING"]
            }}
        ' {output.res} > {output.tbl}) &> {log}
        """
