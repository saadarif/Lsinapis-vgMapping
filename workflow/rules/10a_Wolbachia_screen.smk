# ==============================================================================
# WOLBACHIA SCREEN (competitive mapping with bwa mem + mosdepth)
#   (i)   Concatenate the nuclear reference (config: reference), the mtDNA and
#         the Wolbachia assembly (params: run_wolbachia_screen: mtdna_fasta /
#         wolbachia_fasta) into one FASTA and index it for bwa. Every header of
#         the mtDNA FASTA must contain "mitochondri" and every header of the
#         Wolbachia FASTA "wolbachia" (case ignored); see
#         workflow/scripts/concat_screen_reference.py for the checks. It also
#         writes the contig -> genome file (.stb) and one BED per genome.
#   (ii)  Map the trimmed reads of every modern sample, and of the blanks in
#         params: run_wolbachia_screen: blanks_csv, to it: bwa mem per run ->
#         merge runs of the same index -> Picard MarkDuplicates -> merge the
#         indices of a sample -> bam clipOverlap. The same steps as the modern
#         vg mapping (1.1_mapping.smk), with bwa mem in place of vg map.
#   (iii) Split each sample's BAM into nuclear, mitochondrial and Wolbachia BAMs
#         at MAPQ >= params: run_wolbachia_screen: mapq. Mapping to all three
#         genomes at once means a read that fits the mtDNA or Wolbachia better
#         than a nuclear insertion (NUMT) or vice versa goes to the better one,
#         and a read that fits both equally gets MAPQ 0 and is dropped here.
#   (iv)  mosdepth on each of those BAMs, then one table of reads, mean depth
#         and breadth per sample and genome, with the depth ratios between the
#         genomes (Wolbachia / nuclear, mtDNA / nuclear, Wolbachia / mtDNA).
#         No infection call is made; the table is for choosing thresholds.
#
# Only modern samples (and blanks) for now. Historical samples would need their
# own mapping rule (short single-end reads, e.g. bwa aln, and dedup instead of
# Picard), but every path already carries {source}, so the split / mosdepth /
# summary rules would take them as they are.
#
# The bwa mapping rules (bwa_*) are not tied to this reference: {cref} is looked
# up in BWA_REFS, so another competitive reference added to BWA_REFS is mapped
# by the same rules. Its FASTA just has to be written to BWA_REF_DIR/{cref}.fa.
#
# The inStrain screen (10b_InStrain_Wolbachia_screen.smk) maps separately, with
# bowtie2, but to the FASTA built here.
#
# Library (LB) in the read groups is the sample's index, so Picard marks
# duplicates across the lanes of an index (the lanes of an index are one PCR
# library). 1.1_mapping.smk sets LB to the lane at mapping and corrects it
# afterwards in modern_dedup_libraries.
# ==============================================================================
WOLB_PARAMS = config.get("params", {}).get("run_wolbachia_screen", {})
WOLB_MTDNA_FA = WOLB_PARAMS.get("mtdna_fasta", "")
WOLB_WOLBACHIA_FA = WOLB_PARAMS.get("wolbachia_fasta", "")
WOLB_MAPQ = WOLB_PARAMS.get("mapq", 30)
WOLB_THRESHOLDS = [int(t) for t in WOLB_PARAMS.get("breadth_thresholds", [1, 5, 10])]
WOLB_THREADS_BWA = WOLB_PARAMS.get("threads_bwa", config.get("threads_mapping", 16))
WOLB_THREADS_MOSDEPTH = WOLB_PARAMS.get("threads_mosdepth", 4)

WOLB_GENOMES = ["nuclear", "mitochondrion", "wolbachia"]

# Name of the competitive reference, carried in every file name downstream
WOLB_REF = f"{REF_NAME}_mtWolb"
WOLB_TAG = f"q{WOLB_MAPQ}"
WOLB_COV_TAG = f"{WOLB_TAG}.bt{'-'.join(str(t) for t in WOLB_THRESHOLDS)}"

# Competitive references for the bwa mapping rules: name -> FASTA
BWA_REF_DIR = "results/bwa_mapping/reference"
BWA_REFS = {WOLB_REF: f"{BWA_REF_DIR}/{WOLB_REF}.fa"}

# Samples mapped: every modern sample, plus the blanks. The blanks sheet has the
# same columns as modern_samples.tsv (fq2 may be left out or empty for
# single-end reads) and goes through the Snakefile's load_metadata.
WOLB_BLANKS_CSV = config.get("blanks_csv", "")
WOLB_SOURCES = "modern|blank"

_wolb_on = config.get("run_wolbachia_screen", False)
# 10b maps the same samples to the same reference, so the inputs are checked
# when either screen is on
_wolb_inputs_on = _wolb_on or config.get("run_instrain_wolbachia", False)
if _wolb_inputs_on:
    for _key, _path in [("mtdna_fasta", WOLB_MTDNA_FA), ("wolbachia_fasta", WOLB_WOLBACHIA_FA)]:
        if not _path or not os.path.exists(_path):
            print(f"\nERROR: params: run_wolbachia_screen: {_key} '{_path}' does not exist.\n")
            sys.exit(1)
    if WOLB_BLANKS_CSV and not os.path.exists(WOLB_BLANKS_CSV):
        print(f"\nERROR: blanks_csv '{WOLB_BLANKS_CSV}' does not exist.\n")
        sys.exit(1)
    if not WOLB_THRESHOLDS:
        print("\nERROR: params: run_wolbachia_screen: breadth_thresholds needs at least one value.\n")
        sys.exit(1)

wolb_blanks_df = load_metadata("blanks_csv", "blank")
_clash = set(wolb_blanks_df["sample_id"]) & set(samples_df["sample_id"])
if _wolb_inputs_on and _clash:
    print(f"\nERROR: blank sample_id(s) {sorted(_clash)} are also sample_ids in the sample sheets; give the blanks their own IDs.\n")
    sys.exit(1)

# fillna again after the concat: columns only one sheet has (e.g. rgGroup_ID,
# platform) are NaN for the other's rows
WOLB_SAMPLES_DF = pd.concat(
    [samples_df[samples_df["source"] == "modern"], wolb_blanks_df], ignore_index=True
).fillna("unknown")

# (source, sample_id) of every sample in the screen
WOLB_UNITS = list(
    WOLB_SAMPLES_DF[["source", "sample_id"]].drop_duplicates().itertuples(index=False, name=None)
)


def wolb_genome_bam(source, sample_id, genome):
    return (f"results/wolbachia_screen/bams/{source}/{sample_id}.{WOLB_REF}"
            f".merged.dedup.merged.clipped.{WOLB_TAG}.{genome}.bam")


def wolb_mosdepth_prefix(source, sample_id, genome):
    return (f"results/wolbachia_screen/mosdepth/{source}/{sample_id}.{WOLB_REF}"
            f".merged.dedup.merged.clipped.{WOLB_COV_TAG}.{genome}")


WOLB_SUMMARY_PREFIX = f"results/wolbachia_screen/wolbachia_screen.{WOLB_REF}.{WOLB_COV_TAG}"

# Requested from the Snakefile target list (section 14)
WOLB_SCREEN_TARGETS = []
if _wolb_on:
    WOLB_SCREEN_TARGETS = [
        f"{WOLB_SUMMARY_PREFIX}.coverage_long.tsv",
        f"{WOLB_SUMMARY_PREFIX}.coverage_wide.tsv",
    ] + [
        wolb_genome_bam(src, sid, g) for src, sid in WOLB_UNITS for g in WOLB_GENOMES
    ]


# ------------------------------------------------------------------------------
# Helper functions for the bwa mapping rules
# ------------------------------------------------------------------------------
def bwa_ref_fasta(wildcards):
    if wildcards.cref not in BWA_REFS:
        raise ValueError(f"No competitive reference '{wildcards.cref}' in BWA_REFS ({list(BWA_REFS)})")
    return BWA_REFS[wildcards.cref]


def bwa_run_row(wildcards):
    match = WOLB_SAMPLES_DF[(WOLB_SAMPLES_DF.sample_id == wildcards.sample_id) &
                           (WOLB_SAMPLES_DF.run_id == wildcards.run_id) &
                           (WOLB_SAMPLES_DF.source == wildcards.source)]
    if match.empty:
        raise ValueError(f"No run '{wildcards.run_id}' of {wildcards.source} sample "
                         f"'{wildcards.sample_id}' in the sample sheets")
    return match.iloc[0].to_dict()


def bwa_reads(wildcards):
    """fq1, plus fq2 when the run is paired-end."""
    row = bwa_run_row(wildcards)
    return [row["fq1"]] + ([row["fq2"]] if row["fq2"] != "unknown" else [])


def bwa_rg_string(wildcards):
    """Read group as in 1.1_mapping.smk, except LB is the index (see top)."""
    row = bwa_run_row(wildcards)
    rg_id = row.get("rgGroup_ID", "unknown")
    if rg_id == "unknown":
        rg_id = wildcards.run_id
    platform = row.get("platform", "unknown")
    if platform == "unknown":
        platform = "ILLUMINA"
    index = wildcards.run_id.split("_")[0]
    return (f"@RG\\tID:{wildcards.sample_id}_{rg_id}"
            f"\\tSM:{wildcards.sample_id}"
            f"\\tLB:{wildcards.sample_id}_{index}"
            f"\\tPL:{platform}")


def bwa_runs_for_index(wildcards):
    """Every run (lane) of a sample with this index, as in get_runs_for_index (1.1)."""
    subset = WOLB_SAMPLES_DF[(WOLB_SAMPLES_DF.sample_id == wildcards.sample_id) &
                            (WOLB_SAMPLES_DF.source == wildcards.source)]
    runs = [r for r in subset["run_id"] if r.split("_")[0] == wildcards.index]
    return expand("results/bwa_mapping/{source}/{sample_id}.{run_id}.{cref}.bam",
                  source=wildcards.source, sample_id=wildcards.sample_id,
                  run_id=runs, cref=wildcards.cref)


def bwa_indices_for_sample(wildcards):
    """Every deduplicated index BAM of a sample, as in get_indices_for_sample (1.1)."""
    subset = WOLB_SAMPLES_DF[(WOLB_SAMPLES_DF.sample_id == wildcards.sample_id) &
                            (WOLB_SAMPLES_DF.source == wildcards.source)]
    indices = sorted({r.split("_")[0] for r in subset["run_id"]})
    return expand("results/bwa_mapping/{source}/{sample_id}.{index}.{cref}.merged.dedup.bam",
                  source=wildcards.source, sample_id=wildcards.sample_id,
                  index=indices, cref=wildcards.cref)


# ------------------------------------------------------------------------------
# (i) Competitive reference
# ------------------------------------------------------------------------------
rule wolb_concat_reference:
    """
    Nuclear + mtDNA + Wolbachia in one FASTA, with the contig -> genome file
    and one BED per genome. Fails if a header does not name the genome its
    FASTA is given as, or a contig name occurs twice.
    """
    input:
        nuclear = f"{config['reference']}.fa",
        mtdna = WOLB_MTDNA_FA,
        wolbachia = WOLB_WOLBACHIA_FA,
    output:
        fasta = f"{BWA_REF_DIR}/{WOLB_REF}.fa",
        stb = f"{BWA_REF_DIR}/{WOLB_REF}.stb",
        bed_nuclear = f"{BWA_REF_DIR}/{WOLB_REF}.nuclear.bed",
        bed_mitochondrion = f"{BWA_REF_DIR}/{WOLB_REF}.mitochondrion.bed",
        bed_wolbachia = f"{BWA_REF_DIR}/{WOLB_REF}.wolbachia.bed",
    log:
        f"logs/wolbachia_screen/concat_reference/{WOLB_REF}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/concat_screen_reference.py"


rule bwa_index_reference:
    """bwa index + samtools faidx of a competitive reference."""
    input:
        fa = BWA_REF_DIR + "/{cref}.fa",
    output:
        multiext(BWA_REF_DIR + "/{cref}.fa", ".amb", ".ann", ".bwt", ".pac", ".sa", ".fai"),
    wildcard_constraints:
        cref = "[^./]+",
    log:
        "logs/bwa_mapping/index/{cref}.log",
    benchmark:
        "benchmarks/bwa_mapping/index/{cref}.json"
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        bwa index {input.fa} 2> {log}
        samtools faidx {input.fa} 2>> {log}
        """


# ------------------------------------------------------------------------------
# (ii) bwa mapping, same chain as the modern vg mapping in 1.1_mapping.smk.
#      Everything up to the clipped BAM is temp(), like the vg mapping chain
#      (see the note in 1.1_mapping.smk). The clipped BAM is kept.
# ------------------------------------------------------------------------------
rule bwa_map:
    """bwa mem per run (paired-end, or single-end when fq2 is missing)."""
    input:
        reads = bwa_reads,
        ref = bwa_ref_fasta,
        idx = lambda w: [bwa_ref_fasta(w) + ext for ext in (".amb", ".ann", ".bwt", ".pac", ".sa")],
    output:
        bam = temp("results/bwa_mapping/{source}/{sample_id}.{run_id}.{cref}.bam"),
        flagstat = "results/bwa_mapping/{source}/stats/bwa/{sample_id}.{run_id}.{cref}.flagstat.txt",
    wildcard_constraints:
        source = WOLB_SOURCES,
        cref = "[^./]+",
    params:
        rg = bwa_rg_string,
    threads: WOLB_THREADS_BWA
    log:
        "logs/bwa_mapping/bwa_mem/{source}/{sample_id}.{run_id}.{cref}.log",
    benchmark:
        "benchmarks/bwa_mapping/bwa_mem/{source}/{sample_id}.{run_id}.{cref}.json"
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        bwa mem -t {threads} -R "{params.rg}" {input.ref} {input.reads} 2> {log} | \
            samtools sort -@ 4 -o {output.bam} - 2>> {log}
        samtools flagstat {output.bam} > {output.flagstat}
        """


rule bwa_merge_index_lanes:
    """Merges runs with the same index. If only one run exists, it renames for consistency."""
    input:
        bwa_runs_for_index,
    output:
        temp("results/bwa_mapping/{source}/{sample_id}.{index}.{cref}.merged.bam"),
    wildcard_constraints:
        source = WOLB_SOURCES,
        index = "[^.]+",
        cref = "[^./]+",
    threads: 4
    conda:
        "../envs/bwa.yaml"
    shell:
        "samtools merge -@ {threads} {output} {input}"


rule bwa_dedup:
    """Picard MarkDuplicates, duplicates removed (as deduplicate_modern in 1.1)."""
    input:
        bam = "results/bwa_mapping/{source}/{sample_id}.{index}.{cref}.merged.bam",
    output:
        bam = temp("results/bwa_mapping/{source}/{sample_id}.{index}.{cref}.merged.dedup.bam"),
        metrics = "results/bwa_mapping/{source}/stats/merged_dedup/{sample_id}.{index}.{cref}.merged_picMetrics.txt",
    wildcard_constraints:
        source = WOLB_SOURCES,
        index = "[^.]+",
        cref = "[^./]+",
    threads: config.get("threads_dedup", 4)
    log:
        "logs/bwa_mapping/dedup/{source}/{sample_id}.{index}.{cref}.dedup.log",
    benchmark:
        "benchmarks/bwa_mapping/dedup/{source}/{sample_id}.{index}.{cref}.json"
    conda:
        "../envs/dedup.yaml"
    shell:
        """
        picard MarkDuplicates -Xmx10000m \
            I={input.bam} O={output.bam} M={output.metrics} \
            REMOVE_DUPLICATES=true &> {log}
        """


rule bwa_merge_sample_indices:
    """Merges different indices for the same sample. Renames if only one index exists."""
    input:
        bwa_indices_for_sample,
    output:
        bam = temp("results/bwa_mapping/{source}/{sample_id}.{cref}.merged.dedup.merged.bam"),
    wildcard_constraints:
        source = WOLB_SOURCES,
        cref = "[^./]+",
    threads: 4
    conda:
        "../envs/bwa.yaml"
    shell:
        "samtools merge -@ {threads} {output.bam} {input}"


rule bwa_clip_overlap:
    """bam clipOverlap, so overlapping mates aren't counted twice (as modern_clip_overlap in 1.1)."""
    input:
        "results/bwa_mapping/{source}/{sample_id}.{cref}.merged.dedup.merged.bam",
    output:
        bam = "results/bwa_mapping/{source}/{sample_id}.{cref}.merged.dedup.merged.clipped.bam",
        bai = "results/bwa_mapping/{source}/{sample_id}.{cref}.merged.dedup.merged.clipped.bam.bai",
        log = "results/bwa_mapping/{source}/stats/merged_dedup_merged_clipped/{sample_id}.{cref}.merged.dedup.merged.clipped.log",
    wildcard_constraints:
        source = WOLB_SOURCES,
        cref = "[^./]+",
    benchmark:
        "benchmarks/bwa_mapping/clip_overlap/{source}/{sample_id}.{cref}.json"
    conda:
        "../envs/bamutil.yaml"
    shell:
        """
        bam clipOverlap --in {input} --out {output.bam} --stats 2> {output.log}
        samtools index {output.bam}
        """


# ------------------------------------------------------------------------------
# (iii) + (iv) Split by genome, depth and breadth
# ------------------------------------------------------------------------------
rule wolb_split_genome:
    """
    Reads mapped to one genome's contigs at MAPQ >= mapq. -M with -L reads the
    BAM through its index, region by region, without writing a read twice.
    The nreads file is the number of primary mapped reads in the output BAM.
    """
    input:
        bam = f"results/bwa_mapping/{{source}}/{{sample_id}}.{WOLB_REF}.merged.dedup.merged.clipped.bam",
        bai = f"results/bwa_mapping/{{source}}/{{sample_id}}.{WOLB_REF}.merged.dedup.merged.clipped.bam.bai",
        bed = f"{BWA_REF_DIR}/{WOLB_REF}.{{genome}}.bed",
    output:
        bam = wolb_genome_bam("{source}", "{sample_id}", "{genome}"),
        bai = wolb_genome_bam("{source}", "{sample_id}", "{genome}") + ".bai",
        nreads = wolb_genome_bam("{source}", "{sample_id}", "{genome}") + ".nreads.txt",
    wildcard_constraints:
        source = WOLB_SOURCES,
        genome = "|".join(WOLB_GENOMES),
    params:
        mapq = WOLB_MAPQ,
    threads: 4
    log:
        f"logs/wolbachia_screen/split_genome/{{source}}/{{sample_id}}.{WOLB_REF}.{WOLB_TAG}.{{genome}}.log",
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        samtools view -@ {threads} -b -q {params.mapq} -M -L {input.bed} \
            -o {output.bam} {input.bam} 2> {log}
        samtools index {output.bam} 2>> {log}
        samtools view -@ {threads} -c -F 0x904 {output.bam} > {output.nreads} 2>> {log}
        """


rule wolb_mosdepth:
    """
    Depth per contig (summary "bases") and bases at >= each breadth threshold
    (thresholds.bed.gz) over one genome's BED. -n: no per-base output.
    """
    input:
        bam = wolb_genome_bam("{source}", "{sample_id}", "{genome}"),
        bai = wolb_genome_bam("{source}", "{sample_id}", "{genome}") + ".bai",
        bed = f"{BWA_REF_DIR}/{WOLB_REF}.{{genome}}.bed",
    output:
        summary = wolb_mosdepth_prefix("{source}", "{sample_id}", "{genome}") + ".mosdepth.summary.txt",
        thresholds = wolb_mosdepth_prefix("{source}", "{sample_id}", "{genome}") + ".thresholds.bed.gz",
        regions = wolb_mosdepth_prefix("{source}", "{sample_id}", "{genome}") + ".regions.bed.gz",
    wildcard_constraints:
        source = WOLB_SOURCES,
        genome = "|".join(WOLB_GENOMES),
    params:
        prefix = lambda w: wolb_mosdepth_prefix(w.source, w.sample_id, w.genome),
        mapq = WOLB_MAPQ,
        thresholds = ",".join(str(t) for t in WOLB_THRESHOLDS),
    threads: WOLB_THREADS_MOSDEPTH
    log:
        f"logs/wolbachia_screen/mosdepth/{{source}}/{{sample_id}}.{WOLB_REF}.{WOLB_COV_TAG}.{{genome}}.log",
    benchmark:
        f"benchmarks/wolbachia_screen/mosdepth/{{source}}/{{sample_id}}.{WOLB_REF}.{WOLB_COV_TAG}.{{genome}}.json"
    conda:
        "../envs/mosdepth.yaml"
    shell:
        """
        mosdepth -t {threads} -n -Q {params.mapq} --by {input.bed} \
            --thresholds {params.thresholds} {params.prefix} {input.bam} 2> {log}
        """


rule wolb_coverage_summary:
    """Reads, mean depth, breadth and depth ratios for every sample and genome."""
    input:
        summary = [wolb_mosdepth_prefix(s, i, g) + ".mosdepth.summary.txt"
                   for s, i in WOLB_UNITS for g in WOLB_GENOMES],
        thresholds = [wolb_mosdepth_prefix(s, i, g) + ".thresholds.bed.gz"
                      for s, i in WOLB_UNITS for g in WOLB_GENOMES],
        nreads = [wolb_genome_bam(s, i, g) + ".nreads.txt"
                 for s, i in WOLB_UNITS for g in WOLB_GENOMES],
    output:
        long = f"{WOLB_SUMMARY_PREFIX}.coverage_long.tsv",
        wide = f"{WOLB_SUMMARY_PREFIX}.coverage_wide.tsv",
    params:
        units = [(s, i, g) for s, i in WOLB_UNITS for g in WOLB_GENOMES],
        thresholds = WOLB_THRESHOLDS,
        mapq = WOLB_MAPQ,
    log:
        f"logs/wolbachia_screen/summary/wolbachia_screen.{WOLB_REF}.{WOLB_COV_TAG}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/wolbachia_screen_summary.py"
