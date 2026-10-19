# ==============================================================================
# 1. HELPER FUNCTIONS FOR RULES
# ==============================================================================
#Move these to a common.smk file if these are resued elsewhere in the ...
#For now, they are specific to the mapping rules.
from types import SimpleNamespace

def get_mapping_input(wildcards):
    # Filter for the specific combination
    match = samples_df[(samples_df.sample_id == wildcards.sample_id) & 
                       (samples_df.run_id == wildcards.run_id)]
    
    if match.empty:
        # This error message will save you hours of debugging
        available_runs = samples_df[samples_df.sample_id == wildcards.sample_id]['run_id'].tolist()
        raise ValueError(
            f"\n[METADATA ERROR]\n"
            f"Snakemake is looking for Sample: '{wildcards.sample_id}' with Run: '{wildcards.run_id}'\n"
            f"In your file, Sample '{wildcards.sample_id}' only has these runs: {available_runs}\n"
            f"Check if an underscore in your ID is causing a wildcard split error!"
        )
    
    row = match.iloc[0]
    if row['source'] == 'modern' and row['fq2'] != "unknown":
        return {"r1": row['fq1'], "r2": row['fq2']}
    return {"r1": row['fq1']}

def get_run_metadata(wildcards):
    """Returns the full row of metadata for Read Group (@RG) construction."""
    return samples_df[(samples_df.sample_id == wildcards.sample_id) & 
                     (samples_df.run_id == wildcards.run_id)].iloc[0].to_dict()


def get_rg_string(wildcards):
    """Constructs the Read Group string by calling the metadata lookup."""
    meta = get_run_metadata(wildcards)

    #Split the run_id (e.g., '01_L2') and get the second part ('L2')
    # If run_id is '01_L2', run_suffix becomes 'L2'
    run_suffix = wildcards.run_id.split('_')[1] if '_' in wildcards.run_id else wildcards.run_id

    return (
        f"@RG\\tID:{wildcards.sample_id}_{meta.get('rgGroup_ID', 'lib1')}"
        f"\\tSM:{wildcards.sample_id}"
        f"\\tLB:{run_suffix}"
        f"\\tPL:{meta.get('platform', 'ILLUMINA')}")

# def get_runs_to_merge(wildcards):
#     """
#     Finds all runs for a sample and matches the correct 'source' tag.
#     """
#     subset = samples_df[samples_df.sample_id == wildcards.sample_id]
#     source = subset.iloc[0]['source']  # Get source (modern/historical)
    
#     return expand(
#         "mapped/{{sample_id}}.{run_id}.{{ref_name}}.{src}.reordered.bam",
#         run_id=subset['run_id'],
#         src=source)

def get_runs_for_index(wildcards):
    """Gathers all runs (lanes) that share the same sample_id and index.
     This is used for merging runs before deduplication, ensuring that only runs with the same index are merged together.
    """
    # 1. Get all runs for this sample and source
    subset = samples_df[(samples_df.sample_id == wildcards.sample_id) & 
                        (samples_df.source == wildcards.source)]
    
    # 2. Filter runs where the first part of run_id matches the wildcards.index
    # e.g., if run_id is '01_L2', r.split('_')[0] returns '01'
    matching_runs = [r for r in subset['run_id'] if r.split('_')[0] == wildcards.index]
    
    return expand("results/mapping/{source}/{sample_id}.{run_id}.{ref_name}.{source}.reordered.bam",
                  sample_id=wildcards.sample_id,
                  run_id=matching_runs,
                  ref_name=wildcards.ref_name,
                  source=wildcards.source)

def get_library_lb_map(wildcards):
    """
    "RG ID=LB;RG ID=LB;..." for every run of a modern sample, with LB set to
    {sample_id}_{index}. The RG IDs are rebuilt with get_rg_string, so they are
    exactly the ones written at mapping. Used by modern_dedup_libraries.
    """
    runs = samples_df[(samples_df.sample_id == wildcards.sample_id) &
                      (samples_df.source == "modern")]
    pairs = []
    for run_id in runs["run_id"]:
        rg = get_rg_string(SimpleNamespace(sample_id=wildcards.sample_id, run_id=run_id))
        rg_id = rg.split("\\t")[1][len("ID:"):]
        pairs.append(f"{rg_id}={wildcards.sample_id}_{run_id.split('_')[0]}")
    return ";".join(pairs)

def get_indices_for_sample(wildcards):
    """Gathers all deduplicated index files for a final sample merge."""
    # 1. Get all runs for this sample
    subset = samples_df[(samples_df.sample_id == wildcards.sample_id) & 
                        (samples_df.source == wildcards.source)]
    
    # 2. Extract unique indices from the run_ids
    unique_indices = list(set([r.split('_')[0] for r in subset['run_id']]))
    
    return expand("results/mapping/{source}/{sample_id}.{index}.{ref_name}.merged.dedup.bam",
                  source=wildcards.source,
                  sample_id=wildcards.sample_id,
                  index=unique_indices,
                  ref_name=REF_NAME)

# ==============================================================================
# 2. PROCESSING RULES (ALIGNMENT, DEDUP, CLIPPING)
# ==============================================================================
# NOTE ON temp() AND THE PER-RUN BAM CHAIN
# Every BAM between the fastqs and the per-sample merged BAM is wrapped in
# temp(), i.e. the outputs of map_modern / map_historical -> picard_reorder ->
# merge_index_lanes -> deduplicate_modern / deduplicate_historical. There is one
# of these per sequencing run per sample, so keeping them would multiply the disk
# footprint of the run for files that are never used again once
# merge_sample_indices has produced {sample_id}.{ref_name}.merged.dedup.merged.bam.
# Snakemake deletes each one as soon as the last job consuming it finishes.
#
# Their absence from results/mapping/ is therefore expected and does NOT by
# itself schedule any remapping: if the downstream merged BAM exists and is
# newer, the chain is up to date and Snakemake stops recursing there. What does
# schedule a remap is a provenance rerun trigger (a changed env yaml or rule
# body) firing on one of these rules, because that marks its output as needing
# regeneration and every dependent job then reruns as well. See
# profiles/default/config.yaml, which pins rerun-triggers to mtime to avoid that.
# The flip side of temp() is that once a remap IS triggered it starts from the
# fastqs, since none of the intermediates are on disk to resume from.
#
# Modern samples then go merged.dedup.merged -> modern_dedup_libraries (temp,
# duplicates between lanes of an index removed) -> modern_clip_overlap.

rule map_modern:
    """Mapping for paired-end modern reads."""
    input:
        r1 = lambda wildcards: samples_df[(samples_df.sample_id == wildcards.sample_id) & (samples_df.run_id == wildcards.run_id)].iloc[0]['fq1'],
        r2 = lambda wildcards: samples_df[(samples_df.sample_id == wildcards.sample_id) & (samples_df.run_id == wildcards.run_id)].iloc[0]['fq2'],
        ref = f"{config['reference']}.fa",
        xg  = f"{config['vg_prefix']}.xg",
        gcsa = f"{config['vg_prefix']}.gcsa"
    output:
        bam = temp("results/mapping/modern/{sample_id}.{run_id}.{ref_name}.modern.bam"),
        #TODO: add temp bai file
        flagstat = "results/mapping/modern/stats/vgmap/{sample_id}.{run_id}.{ref_name}.flagstat.txt"
    conda: "../envs/vg.yaml"
    threads: config.get("threads_mapping", 16)
    params:
        rg = get_rg_string
    log:
        "logs/mapping/vgmap/{sample_id}.{run_id}.{ref_name}.modern.log"
    benchmark:
        "benchmarks/mapping/vgmap/{sample_id}.{run_id}.{ref_name}.modern.json"
    shell:
        """
        
        RG_STR="{params.rg}"

        vg map -t {threads} --log-time \
          -x {input.xg} -g {input.gcsa} \
          -f {input.r1} -f {input.r2} \
          --surject-to bam | \
          samtools addreplacerg -u -r "$RG_STR" - | \
          samtools calmd -u - {input.ref} | \
          samtools sort -o {output.bam} &> {log}
        
        samtools index {output.bam}
        samtools flagstat {output.bam} > {output.flagstat}
        """

rule map_historical:
    """Mapping for single-end historical reads."""
    input:
        r1 = lambda wildcards: samples_df[(samples_df.sample_id == wildcards.sample_id) & (samples_df.run_id == wildcards.run_id)].iloc[0]['fq1'],
        ref = f"{config['reference']}.fa",
        xg  = f"{config['vg_prefix']}.xg",
        gcsa = f"{config['vg_prefix']}.gcsa",
        vg = f"{config['vg_prefix']}.vg"
    output:
        bam = temp("results/mapping/historical/{sample_id}.{run_id}.{ref_name}.historical.bam"),
        bai = temp("results/mapping/historical/{sample_id}.{run_id}.{ref_name}.historical.bam.bai"),
        flagstat = "results/mapping/historical/stats/vgmap/{sample_id}.{run_id}.{ref_name}.flagstat.txt"
    conda: "../envs/vg.yaml"
    threads: config.get("threads_mapping", 16)
    params:
        rg = get_rg_string
    log:
        "logs/mapping/vgmap/{sample_id}.{run_id}.{ref_name}.historical.log"
    benchmark:
        "benchmarks/mapping/vgmap/{sample_id}.{run_id}.{ref_name}.historical.json"
    shell:
        """
        RG_STR="{params.rg}"

        vg map -t {threads} -w 300 -k 15 --log-time \
          -x {input.xg} -g {input.gcsa} \
          -f {input.r1} \
          --surject-to bam | \
          samtools addreplacerg -u -r "$RG_STR" - | \
          samtools calmd -u - {input.ref} | \
          samtools sort -o {output.bam} &> {log}

        samtools index {output.bam}
        samtools flagstat {output.bam} > {output.flagstat}
        """

rule picard_reorder:
    """
    Uses picard to sort everything in the mapped BAMs in the same order as the
    reference fasta, which vg doesn't do automatically.
    """
    input:
        bam="results/mapping/{source}/{sample_id}.{run_id}.{ref_name}.{source}.bam",
        ref = f"{config['reference']}.fa",
        fai = f"{config['reference']}.fa.fai",
    output:
        bam=temp("results/mapping/{source}/{sample_id}.{run_id}.{ref_name}.{source}.reordered.bam")
    log:
        "logs/mapping/picard_reorder/{sample_id}.{run_id}.{ref_name}.{source}.log"
    conda: "../envs/vg.yaml"
    benchmark:
        "benchmarks/mapping/picard_reorder/{sample_id}.{run_id}.{ref_name}.{source}.json"
    shell:
        """
        picard ReorderSam -Xmx50000m \
            --INPUT {input.bam} \
            --SEQUENCE_DICTIONARY {input.ref} \
            --OUTPUT {output.bam} &> {log}
        """


rule merge_index_lanes:
    """Merges runs with the same index. If only one run exists, it renames for consistency."""
    input: get_runs_for_index
    output: temp("results/mapping/{source}/{sample_id}.{index,[^.]+}.{ref_name,[^.]+}.merged.bam")
    conda: "../envs/vg.yaml"
    threads: 4
    shell:
        "samtools merge -@ {threads} {output} {input}"



rule deduplicate_modern:
    """
    Deduplicates modern paired-end samples using Picard MarkDuplicates.
    """
    input: 
        bam = "results/mapping/modern/{sample_id}.{index}.{ref_name}.merged.bam"
    output:
        bam = temp("results/mapping/modern/{sample_id}.{index}.{ref_name}.merged.dedup.bam"),
        metrics = "results/mapping/modern/stats/merged_dedup/{sample_id}.{index}.{ref_name}.merged_picMetrics.txt"
    conda: "../envs/dedup.yaml"
    threads: config.get("threads_dedup", 4)
    log: "logs/mapping/dedup_modern/{sample_id}.{index}.{ref_name}.dedup.log"
    benchmark: "benchmarks/mapping/dedup_modern/{sample_id}.{index}.{ref_name}.json"
    shell:
        """
        picard MarkDuplicates -Xmx10000m \
            I={input.bam} O={output.bam} M={output.metrics} \
            REMOVE_DUPLICATES=true &> {log}

        """

rule deduplicate_historical:
    """
    Deduplicates historical single-end samples using the 'dedup' tool.
    """
    input: 
        bam = "results/mapping/historical/{sample_id}.{index}.{ref_name}.merged.bam"
    output:
        bam=temp("results/mapping/historical/dedup/{sample_id}.{index}.{ref_name}.merged_rmdup.bam"),  # Temporary output from dedup
        bamfin = temp("results/mapping/historical/{sample_id}.{index}.{ref_name}.merged.dedup.bam"),
        baifin = temp("results/mapping/historical/{sample_id}.{index}.{ref_name}.merged.dedup.bam.bai"),
        json = "results/mapping/historical/stats/merged_dedup/{sample_id}.{index}.{ref_name}.merged.dedup.json",
        hist = "results/mapping/historical/stats/merged_dedup/{sample_id}.{index}.{ref_name}.merged.dedup.hist",
        log = "results/mapping/historical/stats/merged_dedup/{sample_id}.{index}.{ref_name}.merged.dedup.log"
    conda: "../envs/dedup.yaml"
    threads: config.get("threads_dedup", 4)
    log: "logs/mapping/dedup_historical/{sample_id}.{index}.{ref_name}.dedup.log"
    benchmark: "benchmarks/mapping/dedup_historical/{sample_id}.{index}.{ref_name}.json"
    params:
        # dedup often outputs to a directory based on input name
        outdir=lambda w, output: os.path.dirname(output.bam),
    shell:
        """
        dedup -i {input.bam} -m -u -o {params.outdir} &> {log}
        samtools sort -o {output.bamfin} {output.bam}
        samtools index {output.bamfin}
        
        # Move metadata files to the final QC directory
        mv {params.outdir}/{wildcards.sample_id}.{wildcards.index}.{wildcards.ref_name}.merged.dedup.json {output.json}        
        mv {params.outdir}/{wildcards.sample_id}.{wildcards.index}.{wildcards.ref_name}.merged.hist {output.hist}
        mv {params.outdir}/{wildcards.sample_id}.{wildcards.index}.{wildcards.ref_name}.merged.log {output.log}
        """

rule merge_sample_indices:
    """Merges different indices for the same sample. Renames if only one index exists."""
    input: get_indices_for_sample
    output: 
        bam= "results/mapping/{source}/{sample_id}.{ref_name}.merged.dedup.merged.bam",
        bai= "results/mapping/{source}/{sample_id}.{ref_name}.merged.dedup.merged.bam.bai",
    conda: "../envs/vg.yaml"
    threads: 4
    shell:
        "samtools merge -@ {threads} {output.bam} {input}; samtools index {output.bam}"


rule modern_dedup_libraries:
    """
    Removes duplicates between the lanes of each library of a modern sample.

    The lanes of one index are one PCR library, so duplicates should be removed
    across them; different indices are different libraries and are not
    deduplicated against each other. deduplicate_modern runs on the merged
    lanes of an index, but Picard marks duplicates only within one read-group
    library (LB), and get_rg_string sets LB to the lane, so it removed
    duplicates within each lane only (the Picard metrics in stats/merged_dedup/
    list L1 and L2 as separate libraries).

    Here the LB of every read group is set to {sample_id}_{index} with
    samtools reheader (reads are not changed, only the header), then
    MarkDuplicates is run on the merged BAM. Duplicates within each lane are
    already gone, so this removes the ones between lanes of the same index.

    Runs on the existing merged.dedup.merged BAM, so adding this step reruns
    clipping and everything downstream of the modern BAMs without remapping.
    Historical samples do not go through it. The modern merged.dedup.merged
    BAM itself (and its stats/depth/qualimap) is still the within-lane one;
    the clipped BAM, which every downstream stage uses, is built from this
    rule's output.
    """
    input:
        bam = "results/mapping/modern/{sample_id}.{ref_name}.merged.dedup.merged.bam",
    output:
        header = temp("results/mapping/modern/{sample_id}.{ref_name}.merged.dedup.merged.relabelLB.header.sam"),
        relabelled = temp("results/mapping/modern/{sample_id}.{ref_name}.merged.dedup.merged.relabelLB.bam"),
        bam = temp("results/mapping/modern/{sample_id}.{ref_name}.merged.dedup.merged.libdedup.bam"),
        metrics = "results/mapping/modern/stats/merged_dedup_merged_libdedup/{sample_id}.{ref_name}.merged.dedup.merged.libdedup_picMetrics.txt",
    params:
        lb_map = get_library_lb_map,
    conda: "../envs/dedup.yaml"
    threads: config.get("threads_dedup", 4)
    log: "logs/mapping/dedup_modern_libraries/{sample_id}.{ref_name}.libdedup.log"
    benchmark: "benchmarks/mapping/dedup_modern_libraries/{sample_id}.{ref_name}.json"
    shell:
        """
        samtools view -H {input.bam} | \
            awk -F'\\t' -v OFS='\\t' -v map="{params.lb_map}" '
                BEGIN {{
                    n = split(map, pairs, ";")
                    for (i = 1; i <= n; i++) {{ split(pairs[i], kv, "="); lb[kv[1]] = kv[2] }}
                }}
                /^@RG/ {{
                    id = ""
                    for (i = 2; i <= NF; i++) if ($i ~ /^ID:/) id = substr($i, 4)
                    if (!(id in lb)) {{ print "read group ID " id " is not a run of this sample in the sample sheet" > "/dev/stderr"; exit 1 }}
                    found = 0
                    for (i = 2; i <= NF; i++) if ($i ~ /^LB:/) {{ $i = "LB:" lb[id]; found = 1 }}
                    if (!found) $0 = $0 OFS "LB:" lb[id]
                }}
                {{ print }}' > {output.header} 2> {log}
        samtools reheader {output.header} {input.bam} > {output.relabelled} 2>> {log}
        picard MarkDuplicates -Xmx10000m \
            I={output.relabelled} O={output.bam} M={output.metrics} \
            REMOVE_DUPLICATES=true &>> {log}
        """


rule modern_clip_overlap:
    input: "results/mapping/modern/{sample_id}.{ref_name}.merged.dedup.merged.libdedup.bam"
    output: 
        bam="results/mapping/modern/{sample_id}.{ref_name}.merged.dedup.merged.clipped.bam",
        bai="results/mapping/modern/{sample_id}.{ref_name}.merged.dedup.merged.clipped.bam.bai",
        log="results/mapping/modern/stats/dedup_merged_merged_clipped/{sample_id}.{ref_name}.merged.dedup.merged.clipped.log"
    conda: "../envs/bamutil.yaml"
    benchmark: "benchmarks/mapping/clip_overlap/{sample_id}.{ref_name}.json"
    shell: 
      """
        bam clipOverlap --in {input} --out {output.bam} --stats 2> {output.log}
        samtools index {output.bam}
      """

# ==============================================================================
# 3. QUALITY CONTROL, Rescaling & REPORTING
# ==============================================================================
rule bam_stats_dedup_merged:
    """Generates samtools stats in the source-specific QC folder for deduplicated and merged BAMs."""
    input:  "results/mapping/{source}/{sample_id}.{ref_name}.merged.dedup.merged.bam"
    output: "results/mapping/{source}/stats/merged_dedup_merged/{sample_id}.{ref_name}.merged.dedup.merged.stats.txt"
    log: "logs/mapping/bam_stats/{source}/{sample_id}.{ref_name}.merged.dedup.merged.stats.log"
    conda:  "../envs/vg.yaml"
    shell:  "samtools stats {input} 1> {output} 2> {log}"

rule mapdamage_dedup:
    """Runs mapDamage on historical deduplicated BAMs and stores results in the source-specific QC folder.
    Also creats a rescaled BAM this moved to main results folder for downstream processing."""
    input: 
        bam = "results/mapping/historical/{sample_id}.{ref_name}.merged.dedup.merged.bam",
        bai = "results/mapping/historical/{sample_id}.{ref_name}.merged.dedup.merged.bam.bai",
        ref = f"{config['reference']}.fa"
    output: 
        dir = directory("results/mapping/historical/stats/merged_dedup_merged/mapdamage/{sample_id}.{ref_name}"),
        rescaled_bam = "results/mapping/historical/{sample_id}.{ref_name}.merged.dedup.merged.rescaled.bam",
        rescaled_bai = "results/mapping/historical/{sample_id}.{ref_name}.merged.dedup.merged.rescaled.bam.bai",
    log: "logs/mapping/mapdamage/historical/{sample_id}.{ref_name}.merged.dedup.merged.mapdamage.log"
    conda: "../envs/mapdamage.yaml"
    shell: 
        """
        mapDamage -i {input.bam} -r {input.ref} -d {output.dir} \
        --merge-reference-sequences --rescale 2> {log} &&
        mv {output.dir}/{wildcards.sample_id}.{wildcards.ref_name}.merged.dedup.merged.rescaled.bam {output.rescaled_bam}
        samtools index {output.rescaled_bam} {output.rescaled_bai}
        """

rule calculate_depth_dedup:
    """Calculates mean depth for the deduplicated BAM file over specified BED regions."""
    input:
        bam = "results/mapping/{source}/{sample_id}.{ref_name}.merged.dedup.merged.bam"
    output:
        depth = "results/mapping/{source}/stats/merged_dedup_merged/{sample_id}.{ref_name}.merged.dedup.merged.regfilt.Q20.q30.depth.txt"
    params:
        bed = config.get("site_filter_bed", None),  # BED file for filtering sites, if provided
        mapQ = config.get("mapQ", 30),
        baseQ = config.get("baseQ", 20)
    conda: "../envs/vg.yaml"  # Reusing your environment that contains samtools
    log: "logs/mapping/depth/{source}/{sample_id}.{ref_name}.merged.dedup.merged.regfilt.Q20.q30.depth.log"
    threads: 2
    shell:
        """
        samtools depth -a -b {params.bed} -Q {params.mapQ} -q {params.baseQ} {input.bam} | \
        awk '{{sum+=$3; cnt++}} END {{if(cnt>0) print sum/cnt; else print 0}}' 1> {output.depth} 2> {log}
        """

rule bam_stats_final:
    """Generates samtools stats in the source-specific QC folder for clipped/rescaled BAMs."""
    input:  "results/mapping/{source}/{sample_id}.{ref_name}.merged.dedup.merged.{stage}.bam"
    output: "results/mapping/{source}/stats/merged_dedup_merged_{stage}/{sample_id}.{ref_name}.merged.dedup.merged.{stage}.stats.txt"
    wildcard_constraints:
        stage="clipped|rescaled"
    log: "logs/mapping/bam_stats/{source}/{sample_id}.{ref_name}.merged.dedup.merged.{stage}.stats.log"
    conda: "../envs/vg.yaml"
    shell: "samtools stats {input} 1> {output} 2> {log}"

rule calculate_depth_final:
    """Calculates mean depth for the final clipped/rescaled BAM files."""
    input:
        bam = "results/mapping/{source}/{sample_id}.{ref_name}.merged.dedup.merged.{stage}.bam"
    output:
        depth = "results/mapping/{source}/stats/merged_dedup_merged_{stage}/{sample_id}.{ref_name}.merged.dedup.merged.{stage}.regfilt.Q20.q30.depth.txt"
    params:
        bed = config.get("site_filter_bed", None),  # BED file for filtering sites, if provided
        mapQ = config.get("mapQ", 30),
        baseQ = config.get("baseQ", 20),
    wildcard_constraints:
        stage="clipped|rescaled"
    conda: "../envs/vg.yaml"
    log: "logs/mapping/depth/{source}/{sample_id}.{ref_name}.merged.dedup.merged.{stage}.regfilt.Q20.q30.depth.log"
    threads: 2
    shell:
        """
        samtools depth -a -b {params.bed} -Q {params.mapQ} -q {params.baseQ} {input.bam} | \
        awk '{{sum+=$3; cnt++}} END {{if(cnt>0) print sum/cnt; else print 0}}' 1> {output.depth} 2> {log}
        """

rule qualimap_dedup_merged:
    """Runs Qualimap and stores results in the source-specific QC folder for deduplicated and merged BAMs."""
    input: 
        bam = "results/mapping/{source}/{sample_id}.{ref_name}.merged.dedup.merged.bam",
        bai = "results/mapping/{source}/{sample_id}.{ref_name}.merged.dedup.merged.bam.bai"
    output: 
        txt = "results/mapping/{source}/stats/merged_dedup_merged/qualimap/{sample_id}.{ref_name}.merged.dedup.merged/genome_results.txt"
    params:
        outdir = "results/mapping/{source}/stats/merged_dedup_merged/qualimap/{sample_id}.{ref_name}.merged.dedup.merged"
    log: "logs/mapping/qualimap/{source}/{sample_id}.{ref_name}.merged.dedup.merged.qualimap.log"
    conda: "../envs/qualimap.yaml"
    shell: "qualimap bamqc -bam {input.bam} -outdir {params.outdir} --java-mem-size=8G --outformat HTML &> {log}"

rule qualimap_final:
    """Runs Qualimap and stores results in the source-specific QC folder. Handles both clipped and rescaled BAMs."""
    input: 
        bam = "results/mapping/{source}/{sample_id}.{ref_name}.merged.dedup.merged.{stage}.bam",
        bai = "results/mapping/{source}/{sample_id}.{ref_name}.merged.dedup.merged.{stage}.bam.bai"
    output: 
        txt = "results/mapping/{source}/stats/merged_dedup_merged_{stage}/qualimap/{sample_id}.{ref_name}.merged.dedup.merged.{stage}/genome_results.txt"
    params:
        outdir = "results/mapping/{source}/stats/merged_dedup_merged_{stage}/qualimap/{sample_id}.{ref_name}.merged.dedup.merged.{stage}"
    wildcard_constraints:
        stage="clipped|rescaled"
    log: "logs/mapping/qualimap/{source}/{sample_id}.{ref_name}.merged.dedup.merged.{stage}.qualimap.log"
    conda: "../envs/qualimap.yaml"
    shell: "qualimap bamqc -bam {input.bam} -outdir {params.outdir} --java-mem-size=8G --outformat HTML &> {log}"

rule multiqc_dedup:
    """Runs MultiQC for all deduplicated samples within a source category (modern or historical)."""
    input:
        stats = lambda w: expand("results/mapping/{source}/stats/merged_dedup_merged/{s}.{ref}.merged.dedup.merged.stats.txt", 
                                 source=w.source, 
                                 s=samples_df[samples_df.source == w.source]['sample_id'].unique(), 
                                 ref=REF_NAME),
        qmap  = lambda w: expand("results/mapping/{source}/stats/merged_dedup_merged/qualimap/{s}.{ref}.merged.dedup.merged/genome_results.txt", 
                                 source=w.source, 
                                 s=samples_df[samples_df.source == w.source]['sample_id'].unique(),
                                 ref=REF_NAME)
    output: "results/mapping/{source}/stats/merged_dedup_merged/multiqc_{source}_dedup_merged_report.html"
    log: "logs/mapping/multiqc/{source}/multiqc_{source}_dedup.log"
    conda:  "../envs/multiqc.yaml"
    shell:  
        """
        multiqc results/mapping/{wildcards.source}/stats/merged_dedup_merged/ --force -o results/mapping/{wildcards.source}/stats/merged_dedup_merged/ -n multiqc_{wildcards.source}_dedup_merged_report.html &> {log}
        """

rule multiqc_final:
    """Runs MultiQC for all final processed samples (clipped for modern, rescaled for historical) within a source category."""
    input:
        stats = lambda w: expand("results/mapping/{source}/stats/merged_dedup_merged_{stage}/{s}.{ref}.merged.dedup.merged.{stage}.stats.txt", 
                                 source=w.source, 
                                 stage=w.stage, 
                                 s=samples_df[samples_df.source == w.source]['sample_id'].unique(), 
                                 ref=REF_NAME),
        qmap  = lambda w: expand("results/mapping/{source}/stats/merged_dedup_merged_{stage}/qualimap/{s}.{ref}.merged.dedup.merged.{stage}/genome_results.txt", 
                                 source=w.source, 
                                 stage=w.stage, 
                                 s=samples_df[samples_df.source == w.source]['sample_id'].unique(), 
                                 ref=REF_NAME)
    output: "results/mapping/{source}/stats/merged_dedup_merged_{stage}/multiqc_{source}_{stage}_report.html"
    wildcard_constraints:
        stage="clipped|rescaled" 
    log: "logs/mapping/multiqc/{source}/multiqc_{source}_{stage}.log"
    conda:  "../envs/multiqc.yaml"
    shell:  
        """
        multiqc results/mapping/{wildcards.source}/stats/merged_dedup_merged_{wildcards.stage}/ --force -o results/mapping/{wildcards.source}/stats/merged_dedup_merged_{wildcards.stage}/ -n multiqc_{wildcards.source}_{wildcards.stage}_report.html &> {log}
        """