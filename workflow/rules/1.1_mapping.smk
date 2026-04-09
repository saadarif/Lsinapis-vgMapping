# ==============================================================================
# 1. HELPER FUNCTIONS FOR RULES
# ==============================================================================
#Move these to a common.smk file if these are resued elsewhere in the ... 
#For now, they are specific to the mapping rules.

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
# 2. PROCESSING RULES (ALIGNMENT, DEDUP, CLIPPING, MASKING)
# ==============================================================================

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


rule modern_clip_overlap:
    input: "results/mapping/modern/{sample_id}.{ref_name}.merged.dedup.merged.bam"
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

rule historical_mask:
    input: 
        bam="results/mapping/historical/{sample_id}.{ref_name}.merged.dedup.merged.bam", 
        bed=config["snp_mask_bed"],
    output: 
        bam="results/mapping/historical/{sample_id}.{ref_name}.merged.dedup.merged.masked.bam",
        bai="results/mapping/historical/{sample_id}.{ref_name}.merged.dedup.merged.masked.bam.bai",
        stats="results/mapping/historical/stats/merged_dedup_merged_masked/{sample_id}.{ref_name}.merged.dedup.merged.masked.bam_bamrefine_stats.tx"
    log: "logs/mapping/historical_mask/{sample_id}.{ref_name}.merged.dedup.merged.masked.log"
    benchmark: "benchmarks/mapping/historical_mask/{sample_id}.{ref_name}.json"
    params: 
        extra=config["params"]["bamrefine"]["pmd_length_thresholds"]  # length thresholds for bamrefine
    conda: "../envs/bamrefine.yaml"
    threads: config.get("threads_mask", 8)
    shell: 
        """ 
        bamrefine --snps {input.bed} --threads {threads} {params.extra} --add-tags {input.bam} {output.bam} &> {log}
        samtools index {output.bam}

        # Move the auto-generated stats file to the new target directory)
        mv results/mapping/historical/{wildcards.sample_id}.{wildcards.ref_name}.merged.dedup.merged.masked.bam_bamrefine_stats.txt {output.stats}
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
    """Generates samtools stats in the source-specific QC folder for clipped/masked BAMs."""
    input:  "results/mapping/{source}/{sample_id}.{ref_name}.merged.dedup.merged.{stage}.bam"
    output: "results/mapping/{source}/stats/merged_dedup_merged_{stage}/{sample_id}.{ref_name}.merged.dedup.merged.{stage}.stats.txt"
    wildcard_constraints:
        stage="clipped|masked|rescaled"
    log: "logs/mapping/bam_stats/{source}/{sample_id}.{ref_name}.merged.dedup.merged.{stage}.stats.log"
    conda: "../envs/vg.yaml"
    shell: "samtools stats {input} 1> {output} 2> {log}"

# rule mapdamage_historical_masked:
#     """Runs mapDamage on final masked BAMs for historical samples and stores results in the source-specific QC folder.
#     No rescaling is performed at this point. Only 50% of random reads are sampled to speed up the process."""
#     input: 
#         bam = "results/mapping/historical/{sample_id}.{ref_name}.merged.dedup.merged.masked.bam",
#         bai = "results/mapping/historical/{sample_id}.{ref_name}.merged.dedup.merged.masked.bam.bai",
#         ref = f"{config['reference']}.fa"
#     output: dir = directory("results/mapping/historical/stats/merged_dedup_merged_masked/mapdamage/{sample_id}.{ref_name}")
#     log: "logs/mapping/mapdamage/{sample_id}.{ref_name}.merged.dedup.merged.masked.mapdamage.log"
#     conda: "../envs/mapdamage.yaml"
#     shell: "mapDamage -i {input.bam} -r {input.ref} -d {output.dir} --downsample=0.5 --merge-reference-sequences &> {log}"

rule calculate_depth_final:
    """Calculates mean depth for the final clipped/masked BAM files."""
    input:
        bam = "results/mapping/{source}/{sample_id}.{ref_name}.merged.dedup.merged.{stage}.bam"
    output:
        depth = "results/mapping/{source}/stats/merged_dedup_merged_{stage}/{sample_id}.{ref_name}.merged.dedup.merged.{stage}.regfilt.Q20.q30.depth.txt"
    params:
        bed = config.get("site_filter_bed", None),  # BED file for filtering sites, if provided
        mapQ = config.get("mapQ", 30),
        baseQ = config.get("baseQ", 20),
    wildcard_constraints:
        stage="clipped|masked|rescaled"
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
    """Runs Qualimap and stores results in the source-specific QC folder. Handles both clipped and masked BAMs."""
    input: 
        bam = "results/mapping/{source}/{sample_id}.{ref_name}.merged.dedup.merged.{stage}.bam",
        bai = "results/mapping/{source}/{sample_id}.{ref_name}.merged.dedup.merged.{stage}.bam.bai"
    output: 
        txt = "results/mapping/{source}/stats/merged_dedup_merged_{stage}/qualimap/{sample_id}.{ref_name}.merged.dedup.merged.{stage}/genome_results.txt"
    params:
        outdir = "results/mapping/{source}/stats/merged_dedup_merged_{stage}/qualimap/{sample_id}.{ref_name}.merged.dedup.merged.{stage}"
    wildcard_constraints:
        stage="clipped|masked|rescaled"
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
    """Runs MultiQC for all final processed samples (clipped for modern, masked for historical) within a source category."""
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
        stage="clipped|masked|rescaled" 
    log: "logs/mapping/multiqc/{source}/multiqc_{source}_{stage}.log"
    conda:  "../envs/multiqc.yaml"
    shell:  
        """
        multiqc results/mapping/{wildcards.source}/stats/merged_dedup_merged_{wildcards.stage}/ --force -o results/mapping/{wildcards.source}/stats/merged_dedup_merged_{wildcards.stage}/ -n multiqc_{wildcards.source}_{wildcards.stage}_report.html &> {log}
        """