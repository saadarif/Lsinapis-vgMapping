# ==============================================================================
# WOLBACHIA SCREEN WITH inSTRAIN (bowtie2 + inStrain profile)
#   (i)   the mapping FASTA: by default the nuclear reference (config:
#         reference) + mtDNA + Wolbachia (config: mtdna_fasta /
#         wolbachia_fasta) concatenated here, with the same header checks as
#         10a (workflow/scripts/concat_screen_reference.py), so reads still map
#         competitively; then its bowtie2 index
#   (ii)  per sample, every run (all lanes and indices) mapped in a single
#         bowtie2 call -> sorted BAM, one read group per sample
#   (iii) Picard MarkDuplicates, duplicates removed
#   (iv)  prodigal genes of the profile FASTA (unless genes_fna is given), for
#         inStrain -g
#   (v)   inStrain profile of each BAM against the profile FASTA only (by
#         default the Wolbachia FASTA): coverage, breadth, expected breadth,
#         nucleotide diversity, ANI to the reference, ...
#   (vi)  the genome_info table of every sample in one file
#
# Independent of 10a_Wolbachia_screen.smk: nothing defined there is used here
# (this file builds its own copy of the competitive FASTA, under
# results/wolbachia_screening/bowtie2_mapping/), so either file can be run or included alone.
#
# Samples: every modern sample, paired-end reads only. Historical samples are
# not included yet; every path carries {source} so they can be added.
#
# All runs of a sample share one read group, so Picard treats them as one
# library: duplicates are removed across indices as well as lanes.
#
# The mapping and profile FASTAs can be replaced (params: run_instrain_wolbachia),
# e.g. nuclear + mtDNA + several microbial genomes to map to, and just the
# microbial genomes, with an .stb saying which contig is which genome, to
# profile. The profile FASTA's contigs must be in the mapping FASTA under the
# same names, or inStrain finds no reads for them.
# ==============================================================================
IS_PARAMS = config.get("params", {}).get("run_instrain_wolbachia", {})
IS_MTDNA_FA = config.get("mtdna_fasta", "")
IS_WOLBACHIA_FA = config.get("wolbachia_fasta", "")
IS_REF_DIR = "results/wolbachia_screening/bowtie2_mapping/reference"
# Default mapping FASTA, built by is_concat_reference
IS_CONCAT_NAME = f"{REF_NAME}_mtWolb"
IS_CONCAT_FA = f"{IS_REF_DIR}/{IS_CONCAT_NAME}.fa"
IS_MAP_FA = IS_PARAMS.get("mapping_fasta", "") or IS_CONCAT_FA
IS_MAP_NAME = IS_PARAMS.get("mapping_name", "") or IS_CONCAT_NAME
IS_PROFILE_FA = IS_PARAMS.get("profile_fasta", "") or IS_WOLBACHIA_FA
IS_PROFILE_NAME = IS_PARAMS.get("profile_name", "") or "wolbachia"
IS_STB = IS_PARAMS.get("stb", "")
IS_GENES_FNA = IS_PARAMS.get("genes_fna", "")
IS_PRODIGAL_MODE = IS_PARAMS.get("prodigal_mode", "single")
IS_MIN_MAPQ = IS_PARAMS.get("min_mapq", 2)
IS_THREADS = IS_PARAMS.get("threads", 8)
IS_THREADS_BT2 = IS_PARAMS.get("threads_bowtie2", 16)
IS_BT2_EXTRA = IS_PARAMS.get("bowtie2_extra", "")
IS_EXTRA = IS_PARAMS.get("extra", "")

IS_SOURCES = "modern"

_is_on = config.get("run_instrain_wolbachia", False)
if _is_on:
    # mtdna_fasta is only needed for the default mapping FASTA; wolbachia_fasta
    # for that and as the default profile FASTA
    _needed = []
    if not IS_PARAMS.get("mapping_fasta", ""):
        _needed = [("mtdna_fasta", IS_MTDNA_FA), ("wolbachia_fasta", IS_WOLBACHIA_FA)]
    elif not IS_PARAMS.get("profile_fasta", ""):
        _needed = [("wolbachia_fasta", IS_WOLBACHIA_FA)]
    for _key, _path in _needed:
        if not _path or not os.path.exists(_path):
            print(f"\nERROR: {_key} '{_path}' (config.yaml, Reference Files) does not exist.\n")
            sys.exit(1)
    if IS_PARAMS.get("mapping_fasta", "") and not IS_PARAMS.get("mapping_name", ""):
        print("\nERROR: params: run_instrain_wolbachia: mapping_name is needed when mapping_fasta is set.\n")
        sys.exit(1)
    for _name in (IS_MAP_NAME, IS_PROFILE_NAME):
        if "." in _name or "/" in _name:
            print(f"\nERROR: params: run_instrain_wolbachia: name '{_name}' must not contain '.' or '/'.\n")
            sys.exit(1)
    for _key, _path in [("mapping_fasta", IS_PARAMS.get("mapping_fasta", "")),
                        ("profile_fasta", IS_PROFILE_FA), ("stb", IS_STB),
                        ("genes_fna", IS_GENES_FNA)]:
        if _path and not os.path.exists(_path):
            print(f"\nERROR: params: run_instrain_wolbachia: {_key} '{_path}' does not exist.\n")
            sys.exit(1)
    if IS_PRODIGAL_MODE not in ("single", "meta"):
        print("\nERROR: params: run_instrain_wolbachia: prodigal_mode must be 'single' or 'meta'.\n")
        sys.exit(1)

IS_BT2_INDEX = f"{IS_REF_DIR}/{IS_MAP_NAME}"
IS_BT2_INDEX_FILES = [f"{IS_BT2_INDEX}{ext}" for ext in
                      (".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")]
IS_BAM_STEM = f"results/wolbachia_screening/bowtie2_mapping/{{source}}/{{sample_id}}.{IS_MAP_NAME}.bt2"
IS_GENES = IS_GENES_FNA or (f"results/wolbachia_screening/instrain_wolbachia/reference/"
                            f"{IS_PROFILE_NAME}.prodigal_{IS_PRODIGAL_MODE}.genes.fna")
IS_TAG = f"{IS_MAP_NAME}.bt2.dedup.{IS_PROFILE_NAME}.mq{IS_MIN_MAPQ}"


def is_outdir(source, sample_id):
    return f"results/wolbachia_screening/instrain_wolbachia/{source}/{sample_id}.{IS_TAG}.IS"


def is_genome_info(source, sample_id):
    """inStrain names its tables after the output directory."""
    outdir = is_outdir(source, sample_id)
    return f"{outdir}/output/{os.path.basename(outdir)}_genome_info.tsv"


IS_SUMMARY = f"results/wolbachia_screening/instrain_wolbachia/instrain.{IS_TAG}.genome_info.tsv"

# Requested from the Snakefile target list (section 15)
IS_TARGETS = [IS_SUMMARY] if _is_on else []

# Every modern sample
IS_SAMPLES_DF = samples_df[samples_df["source"] == "modern"]

# (source, sample_id) of every sample profiled
IS_UNITS = list(
    IS_SAMPLES_DF[["source", "sample_id"]].drop_duplicates().itertuples(index=False, name=None)
)


def is_reads(wildcards):
    """fq1 and fq2 of every run (all lanes and indices) of a sample."""
    runs = IS_SAMPLES_DF[(IS_SAMPLES_DF.sample_id == wildcards.sample_id) &
                        (IS_SAMPLES_DF.source == wildcards.source)]
    if runs.empty:
        raise ValueError(f"No runs of {wildcards.source} sample '{wildcards.sample_id}' in the sample sheets")
    if (runs.fq2 == "unknown").any():
        raise ValueError(f"'{wildcards.sample_id}' has runs without fq2; the inStrain screen "
                         "maps paired-end reads only")
    return {"r1": list(runs.fq1), "r2": list(runs.fq2)}


def is_bowtie2_read_args(wildcards):
    """-1 / -2 as comma-separated lists, so every run goes into one bowtie2 call."""
    reads = is_reads(wildcards)
    return f"-1 {','.join(reads['r1'])} -2 {','.join(reads['r2'])}"


rule is_concat_reference:
    """
    Nuclear + mtDNA + Wolbachia in one FASTA (the default mapping FASTA), with
    the same header and contig-name checks as 10a. Only the FASTA is used here;
    the .stb and BEDs are written by the script anyway.
    """
    input:
        nuclear = f"{config['reference']}.fa",
        mtdna = IS_MTDNA_FA,
        wolbachia = IS_WOLBACHIA_FA,
    output:
        fasta = IS_CONCAT_FA,
        stb = f"{IS_REF_DIR}/{IS_CONCAT_NAME}.stb",
        bed_nuclear = f"{IS_REF_DIR}/{IS_CONCAT_NAME}.nuclear.bed",
        bed_mitochondrion = f"{IS_REF_DIR}/{IS_CONCAT_NAME}.mitochondrion.bed",
        bed_wolbachia = f"{IS_REF_DIR}/{IS_CONCAT_NAME}.wolbachia.bed",
    log:
        f"logs/wolbachia_screening/instrain_wolbachia/concat_reference/{IS_CONCAT_NAME}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/concat_screen_reference.py"


rule is_bowtie2_index:
    """bowtie2 index of the mapping FASTA."""
    input:
        fa = IS_MAP_FA,
    output:
        IS_BT2_INDEX_FILES,
    params:
        prefix = IS_BT2_INDEX,
    threads: IS_THREADS_BT2
    log:
        f"logs/wolbachia_screening/bowtie2_mapping/index/{IS_MAP_NAME}.log",
    benchmark:
        f"benchmarks/wolbachia_screening/bowtie2_mapping/index/{IS_MAP_NAME}.json"
    conda:
        "../envs/bowtie2.yaml"
    shell:
        "bowtie2-build --threads {threads} {input.fa} {params.prefix} > {log} 2>&1"


rule is_bowtie2_map:
    """
    Every run of a sample in one bowtie2 call (comma-separated -1/-2 lists),
    end-to-end mode with bowtie2's defaults plus bowtie2_extra. bowtie2's
    alignment summary is kept as a stats file.
    """
    input:
        reads = lambda w: is_reads(w)["r1"] + is_reads(w)["r2"],
        idx = IS_BT2_INDEX_FILES,
    output:
        bam = temp(IS_BAM_STEM + ".bam"),
        summary = f"results/wolbachia_screening/bowtie2_mapping/{{source}}/stats/{{sample_id}}.{IS_MAP_NAME}.bowtie2_summary.txt",
    wildcard_constraints:
        source = IS_SOURCES,
    params:
        index = IS_BT2_INDEX,
        reads = is_bowtie2_read_args,
        extra = IS_BT2_EXTRA,
    threads: IS_THREADS_BT2
    log:
        f"logs/wolbachia_screening/bowtie2_mapping/bowtie2/{{source}}/{{sample_id}}.{IS_MAP_NAME}.log",
    benchmark:
        f"benchmarks/wolbachia_screening/bowtie2_mapping/bowtie2/{{source}}/{{sample_id}}.{IS_MAP_NAME}.json"
    conda:
        "../envs/bowtie2.yaml"
    shell:
        """
        bowtie2 -p {threads} -x {params.index} {params.reads} \
            --rg-id {wildcards.sample_id} --rg SM:{wildcards.sample_id} \
            --rg LB:{wildcards.sample_id} --rg PL:ILLUMINA {params.extra} \
            2> {output.summary} | \
            samtools sort -@ 4 -o {output.bam} - 2> {log}
        """


rule is_dedup:
    """Picard MarkDuplicates, duplicates removed, then indexed for inStrain."""
    input:
        bam = IS_BAM_STEM + ".bam",
    output:
        bam = IS_BAM_STEM + ".dedup.bam",
        bai = IS_BAM_STEM + ".dedup.bam.bai",
        metrics = f"results/wolbachia_screening/bowtie2_mapping/{{source}}/stats/{{sample_id}}.{IS_MAP_NAME}.bt2.picMetrics.txt",
    wildcard_constraints:
        source = IS_SOURCES,
    threads: config.get("threads_dedup", 4)
    log:
        f"logs/wolbachia_screening/bowtie2_mapping/dedup/{{source}}/{{sample_id}}.{IS_MAP_NAME}.log",
    benchmark:
        f"benchmarks/wolbachia_screening/bowtie2_mapping/dedup/{{source}}/{{sample_id}}.{IS_MAP_NAME}.json"
    conda:
        "../envs/dedup.yaml"
    shell:
        """
        picard MarkDuplicates -Xmx10000m \
            I={input.bam} O={output.bam} M={output.metrics} \
            REMOVE_DUPLICATES=true &> {log}
        samtools index {output.bam} 2>> {log}
        """


rule is_prodigal:
    """Genes of the profile FASTA for inStrain -g (skipped when genes_fna is given)."""
    input:
        fa = IS_PROFILE_FA,
    output:
        fna = f"results/wolbachia_screening/instrain_wolbachia/reference/{IS_PROFILE_NAME}.prodigal_{IS_PRODIGAL_MODE}.genes.fna",
        faa = f"results/wolbachia_screening/instrain_wolbachia/reference/{IS_PROFILE_NAME}.prodigal_{IS_PRODIGAL_MODE}.genes.faa",
        gff = f"results/wolbachia_screening/instrain_wolbachia/reference/{IS_PROFILE_NAME}.prodigal_{IS_PRODIGAL_MODE}.genes.gff",
    params:
        mode = IS_PRODIGAL_MODE,
    log:
        f"logs/wolbachia_screening/instrain_wolbachia/prodigal/{IS_PROFILE_NAME}.prodigal_{IS_PRODIGAL_MODE}.log",
    conda:
        "../envs/prodigal.yaml"
    shell:
        """
        prodigal -i {input.fa} -d {output.fna} -a {output.faa} \
            -o {output.gff} -f gff -p {params.mode} 2> {log}
        """


rule is_profile:
    """
    inStrain profile against the profile FASTA only. --min_mapq: a pair is
    used if either read has MAPQ above it (strictly greater). -s is added when
    an .stb is given; without one all contigs are one genome ("all_scaffolds").
    """
    input:
        bam = IS_BAM_STEM + ".dedup.bam",
        bai = IS_BAM_STEM + ".dedup.bam.bai",
        fasta = IS_PROFILE_FA,
        genes = IS_GENES,
    output:
        genome_info = is_genome_info("{source}", "{sample_id}"),
    wildcard_constraints:
        source = IS_SOURCES,
    params:
        outdir = lambda w: is_outdir(w.source, w.sample_id),
        min_mapq = IS_MIN_MAPQ,
        stb = f"-s {IS_STB}" if IS_STB else "",
        extra = IS_EXTRA,
    threads: IS_THREADS
    log:
        f"logs/wolbachia_screening/instrain_wolbachia/profile/{{source}}/{{sample_id}}.{IS_TAG}.log",
    benchmark:
        f"benchmarks/wolbachia_screening/instrain_wolbachia/profile/{{source}}/{{sample_id}}.{IS_TAG}.json"
    conda:
        "../envs/instrain.yaml"
    shell:
        """
        inStrain profile {input.bam} {input.fasta} -o {params.outdir} \
            -p {threads} -g {input.genes} --min_mapq {params.min_mapq} \
            {params.stb} {params.extra} > {log} 2>&1
        """


rule is_genome_info_summary:
    """inStrain's genome_info table of every sample, with sample_id and source columns."""
    input:
        [is_genome_info(s, i) for s, i in IS_UNITS],
    output:
        IS_SUMMARY,
    params:
        units = IS_UNITS,
    log:
        f"logs/wolbachia_screening/instrain_wolbachia/summary/instrain.{IS_TAG}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/instrain_genome_info_summary.py"
