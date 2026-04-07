#==============================================================================
# VARIABLES & SETUP FOR NO-TRANSITIONS WORKFLOW
# ==============================================================================
MIN_DP_NT = config["params"]["run_genotyping"]["minDP"]
MAX_DP_NT = config["params"]["run_genotyping"]["maxDP"]
BASEQ_NT = config["baseQ"]
MAPQ_NT = config["mapQ"]
POPLIST_NT = config["params"]["run_genotyping"]["poplist"]
SITEFILTER_BED = config["site_filter_bed"]  # Using the more conservative site filter bed for no-transitions genotyping

DROP_SAMPLES_NT = config["params"]["run_genotyping"].get("drop_samples", [])
VALID_SAMPLES_DF_NT = samples_df[~samples_df['sample_id'].isin(DROP_SAMPLES_NT)] 
KEEP_SAMPLES_NT = VALID_SAMPLES_DF_NT['sample_id'].tolist()

# DYNAMICALLY FILTER THE POPLIST to ONLY INCLUDE KEPT SAMPLES
# ==============================================================================
ORIGINAL_POPLIST_NT = config["params"]["run_genotyping"]["poplist"]
FILTERED_POPLIST_NT = "results/genotyping/filtered_poplist.txt"

# Only attempt to filter if the original file actually exists
if os.path.exists(ORIGINAL_POPLIST_NT):
    # Read the poplist (using '\s+' handles both spaces and tabs)
    pop_df = pd.read_csv(ORIGINAL_POPLIST_NT, sep=r'\s+', header=None, names=['sample_id', 'pop_id'])
    
    # Filter the dataframe to only include samples in KEEP_SAMPLES
    filtered_pop_df = pop_df[pop_df['sample_id'].isin(KEEP_SAMPLES_NT)]
    
    # Create the output directory if it doesn't exist and save the filtered file
    os.makedirs(os.path.dirname(FILTERED_POPLIST_NT), exist_ok=True)
    filtered_pop_df.to_csv(FILTERED_POPLIST_NT, sep='\t', header=False, index=False)


# ==============================================================================
# INPUT HELPER FUNCTIONS
# ==============================================================================
def get_bam_for_genotyping_notrans(wildcards=None, sid=None):
    """Historical uses base deduplicated BAM; Modern uses clipped."""
    if sid is None:
        sid = wildcards.sample_id
    src = VALID_SAMPLES_DF_NT[VALID_SAMPLES_DF_NT['sample_id'] == sid]['source'].iloc[0]
    
    subsample_samples = config.get("subsample_samples", [])
    target_dp = config.get("subsample_depth", 15)
    
    # Pre-masking for historical means no stage suffix. Modern stays clipped.
    stage_str = ".clipped" if src == "modern" else ""
    mapq = config["mapQ"]
    if sid in subsample_samples:
        return f"results/mapping/{src}/{sid}.{REF_NAME}.merged.dedup.merged{stage_str}.subs{target_dp}.q{mapq}.bam"
    else:
        return f"results/mapping/{src}/{sid}.{REF_NAME}.merged.dedup.merged{stage_str}.bam"

def get_bai_for_genotyping_notrans(wildcards=None, sid=None):
    return get_bam_for_genotyping_notrans(wildcards, sid) + ".bai"

def get_all_bams_for_joint_calling_notrans(wildcards):
    return [get_bam_for_genotyping_notrans(sid=sid) for sid in KEEP_SAMPLES_NT]

def get_all_bais_for_joint_calling_notrans(wildcards):
    return [get_bai_for_genotyping_notrans(sid=sid) for sid in KEEP_SAMPLES_NT]

# ==============================================================================
# CORE CALLING RULES (With Awk transition filter injected)
# ==============================================================================
rule call_individual_genotypes_notrans:
    input:
        alignments = get_bam_for_genotyping_notrans,
        index = get_bai_for_genotyping_notrans,
        ref = config["reference"] + ".fa",
        targets = SITEFILTER_BED
    output:
        bcf = "results/genotyping_notrans/individual/{sample_id}.{ref_name}.sitefilt.bQ" + str(BASEQ_NT) + ".mq" + str(MAPQ_NT) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_NT) + "-" + str(MAX_DP_NT) + ".AB.notrans.bcf",
        idx = "results/genotyping_notrans/individual/{sample_id}.{ref_name}.sitefilt.bQ" + str(BASEQ_NT) + ".mq" + str(MAPQ_NT) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_NT) + "-" + str(MAX_DP_NT) + ".AB.notrans.bcf.csi"
    params:
        min_dp = MIN_DP_NT,
        max_dp = MAX_DP_NT,
        mapq = MAPQ_NT,
        baseq = BASEQ_NT
    log: "logs/genotyping_notrans/individual/{sample_id}_{ref_name}.log"
    benchmark: "benchmarks/genotyping_notrans/individual/{sample_id}_{ref_name}.benchmark"
    conda: "../envs/bcftools121.yaml"
    threads: 8
    shell:
        """
        bcftools mpileup --threads {threads} -f {input.ref} -R {input.targets} \
            -Ou -B --min-MQ {params.mapq} --min-BQ {params.baseq} -a "FORMAT/AD,FORMAT/DP,INFO/AD" {input.alignments} | \
            bcftools call -m -a GQ,GP -Ou | \
            bcftools filter -g 5 -i'QUAL >= 30' -Ou | \
            bcftools view -V indels -M2 -Ov | \
            awk -F '\\t' '/^#/ || !(($4 == "A" && $5 == "G") || ($4 == "G" && $5 == "A") || ($4 == "C" && $5 == "T") || ($4 == "T" && $5 == "C"))' | \
            bcftools view -Ou | \
            bcftools +setGT -Ou -- -t q -n . -i"FMT/DP<{params.min_dp} | FMT/DP>{params.max_dp}" | \
            bcftools +setGT -Ou -- -t q -n . -i'GT="het" & (FMT/AD[:0]/FMT/DP < 0.21 | FMT/AD[:0]/FMT/DP > 0.79 | FMT/AD[:1]/FMT/DP < 0.21 | FMT/AD[:1]/FMT/DP > 0.79)' | \
            bcftools +fill-tags -Ob -- -t all > {output.bcf}
        
        bcftools index -o {output.idx} {output.bcf}
        """

rule merge_genotypes_notrans:
    input:
        bcfs = lambda wildcards: expand("results/genotyping_notrans/individual/{sample_id}.{ref_name}.sitefilt.bQ" + str(BASEQ_NT) + ".mq" + str(MAPQ_NT) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_NT) + "-" + str(MAX_DP_NT) + ".AB.notrans.bcf", sample_id=KEEP_SAMPLES_NT, ref_name=REF_NAME),
        csis = lambda wildcards: expand("results/genotyping_notrans/individual/{sample_id}.{ref_name}.sitefilt.bQ" + str(BASEQ_NT) + ".mq" + str(MAPQ_NT) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_NT) + "-" + str(MAX_DP_NT) + ".AB.notrans.bcf.csi", sample_id=KEEP_SAMPLES_NT, ref_name=REF_NAME)
    output:
        merged_bcf = "results/genotyping_notrans/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ_NT) + ".mq" + str(MAPQ_NT) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_NT) + "-" + str(MAX_DP_NT) + ".AB.indCall.notrans.allsites.bcf",
        csi = "results/genotyping_notrans/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ_NT) + ".mq" + str(MAPQ_NT) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_NT) + "-" + str(MAX_DP_NT) + ".AB.indCall.notrans.allsites.bcf.csi",
        stats = "results/genotyping_notrans/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ_NT) + ".mq" + str(MAPQ_NT) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_NT) + "-" + str(MAX_DP_NT) + ".AB.indCall.notrans.allsites.bcf.stats",
    log: "logs/genotyping_notrans/merged_genotypes_{ref_name}.log"
    benchmark: "benchmarks/genotyping_notrans/merged_genotypes_{ref_name}.benchmark"
    conda: "../envs/bcftools121.yaml"
    threads: 4
    shell:
        """
        bcftools merge --force-samples -Ou {input.bcfs} | \
            bcftools +fill-tags -Ou -- -t all | \
            bcftools view -Ob -o {output.merged_bcf} 2> {log}
        bcftools index -o {output.csi} {output.merged_bcf} 2>> {log}
        bcftools stats -s - {output.merged_bcf} > {output.stats}
        """

rule joint_call_genotypes_notrans:
    input:
        bams = get_all_bams_for_joint_calling_notrans,
        bais = get_all_bais_for_joint_calling_notrans,
        ref = config["reference"] + ".fa",
        fai = f"{config['reference']}.fa.fai",
        targets = SITEFILTER_BED,
        poplist = FILTERED_POPLIST_NT #for HWE filtering
    output:
        bcf = "results/genotyping_notrans/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ_NT) + ".mq" + str(MAPQ_NT) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_NT) + "-" + str(MAX_DP_NT) + ".AB.jointCall.notrans.allsites.bcf",
        csi = "results/genotyping_notrans/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ_NT) + ".mq" + str(MAPQ_NT) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_NT) + "-" + str(MAX_DP_NT) + ".AB.jointCall.notrans.allsites.bcf.csi",
        stats = "results/genotyping_notrans/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ_NT) + ".mq" + str(MAPQ_NT) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_NT) + "-" + str(MAX_DP_NT) + ".AB.jointCall.notrans.allsites.bcf.stats",
    params:
        min_dp = MIN_DP_NT,
        max_dp = MAX_DP_NT,
        mapq = MAPQ_NT,
        baseq = BASEQ_NT
    log: "logs/genotyping_notrans/joint_call_genotypes_{ref_name}.log"
    benchmark: "benchmarks/genotyping_notrans/joint_call_genotypes_{ref_name}.benchmark"
    conda: "../envs/bcftools121.yaml"
    threads: 8
    shell:
        """
        bcftools mpileup --threads {threads} -f {input.ref} -R {input.targets} \
            -a "FORMAT/QS,FORMAT/AD,FORMAT/DP,INFO/AD" -B \
            --min-MQ {params.mapq} --min-BQ {params.baseq}   -Ou {input.bams} | \
            bcftools call -m -a GQ,GP -G {input.poplist} -Ou | \
            bcftools filter -g 5 -i'QUAL >= 30' -Ou | \
            bcftools view -V indels -M2 -Ov | \
            awk -F '\\t' '/^#/ || !(($4 == "A" && $5 == "G") || ($4 == "G" && $5 == "A") || ($4 == "C" && $5 == "T") || ($4 == "T" && $5 == "C"))' | \
            bcftools view -Ou | \
            bcftools +setGT -Ou -- -t q -n . -i"FMT/DP<{params.min_dp} | FMT/DP>{params.max_dp}" | \
            bcftools +setGT -Ou -- -t q -n . -i'GT="het" & (FMT/VAF < 0.21 | FMT/VAF > 0.79') | \
            bcftools +fill-tags -Ob -- -t all > {output.bcf} 2> {log}
            
        bcftools index -o {output.csi} {output.bcf} 2>> {log}
        bcftools stats -s - {output.bcf} > {output.stats}  
        """

# ==============================================================================
# POST-PROCESSING RULES (Tracks the .notrans. tag)
# ==============================================================================
rule filter_variants_notrans:
    input:
        bcf = "results/genotyping_notrans/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ_NT) + ".mq" + str(MAPQ_NT) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_NT) + "-" + str(MAX_DP_NT) + ".AB.{call_type}.notrans.allsites.bcf",
    output:
        bcf = "results/genotyping_notrans/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ_NT) + ".mq" + str(MAPQ_NT) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_NT) + "-" + str(MAX_DP_NT) + ".AB.{call_type}.notrans.biallelic.bcf",
        csi = "results/genotyping_notrans/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ_NT) + ".mq" + str(MAPQ_NT) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_NT) + "-" + str(MAX_DP_NT) + ".AB.{call_type}.notrans.biallelic.bcf.csi",
        stats = "results/genotyping_notrans/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ_NT) + ".mq" + str(MAPQ_NT) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_NT) + "-" + str(MAX_DP_NT) + ".AB.{call_type}.notrans.biallelic.bcf.stats",
    log: "logs/genotyping_notrans/filter_variants_{ref_name}_{call_type}.log"
    conda: "../envs/bcftools121.yaml"
    shell:
        """
        bcftools view -v snps -m 2 -M 2 -i 'MAF>0' -Ob -o {output.bcf} {input.bcf} 2> {log}
        bcftools index -o {output.csi} {output.bcf} 2>> {log}
        bcftools stats -s - {output.bcf} > {output.stats}
        """

rule filter_missingness_notrans:
    input:
        bcf = "results/genotyping_notrans/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ_NT) + ".mq" + str(MAPQ_NT) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_NT) + "-" + str(MAX_DP_NT) + ".AB.{call_type}.notrans.{site_type}.bcf",
        csi = "results/genotyping_notrans/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ_NT) + ".mq" + str(MAPQ_NT) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_NT) + "-" + str(MAX_DP_NT) + ".AB.{call_type}.notrans.{site_type}.bcf.csi"
    output:
        bcf = "results/genotyping_notrans/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ_NT) + ".mq" + str(MAPQ_NT) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_NT) + "-" + str(MAX_DP_NT) + ".AB.{call_type}.notrans.{site_type}.fmiss{max_miss}.bcf",
        csi = "results/genotyping_notrans/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ_NT) + ".mq" + str(MAPQ_NT) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_NT) + "-" + str(MAX_DP_NT) + ".AB.{call_type}.notrans.{site_type}.fmiss{max_miss}.bcf.csi",
        stats = "results/genotyping_notrans/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ_NT) + ".mq" + str(MAPQ_NT) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_NT) + "-" + str(MAX_DP_NT) + ".AB.{call_type}.notrans.{site_type}.fmiss{max_miss}.bcf.stats",
        modsites = temp("results/genotyping_notrans/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ_NT) + ".mq" + str(MAPQ_NT) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_NT) + "-" + str(MAX_DP_NT) + ".AB.{call_type}.notrans.{site_type}.fmiss{max_miss}.modsites.tmp"),
        histsites = temp("results/genotyping_notrans/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ_NT) + ".mq" + str(MAPQ_NT) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_NT) + "-" + str(MAX_DP_NT) + ".AB.{call_type}.notrans.{site_type}.fmiss{max_miss}.histsites.tmp"),
    params:
        modsamps = ",".join(VALID_SAMPLES_DF_NT[VALID_SAMPLES_DF_NT['source'] == 'modern']['sample_id'].tolist()),
        histsamps = ",".join(VALID_SAMPLES_DF_NT[VALID_SAMPLES_DF_NT['source'] == 'historical']['sample_id'].tolist())
    log: "logs/genotyping_notrans/filter_missingness_{ref_name}_{call_type}_{site_type}_fmiss{max_miss}.log"
    conda: "../envs/bcftools121.yaml"
    threads: 4
    shell:
        """
        bcftools view -s {params.histsamps} --force-samples -Ou {input.bcf} | \
            bcftools +fill-tags -Ou -- -t F_MISSING | \
            bcftools query -i 'F_MISSING <= {wildcards.max_miss}' \
                -f '%CHROM\\t%POS\\n' > {output.histsites} 2> {log}
        
        bcftools view -s {params.modsamps} --force-samples -Ou {input.bcf} | \
            bcftools +fill-tags -Ou -- -t F_MISSING | \
            bcftools query -i 'F_MISSING <= {wildcards.max_miss}' \
                -f '%CHROM\\t%POS\\n' > {output.modsites} 2>> {log}
        
        bcftools view -T {output.histsites} -Ou {input.bcf} | \
            bcftools view -T {output.modsites} -Ob -o {output.bcf} 2>> {log}
        
        bcftools index -o {output.csi} {output.bcf} 2>> {log}
        bcftools stats -s - {output.bcf} > {output.stats}
        """
