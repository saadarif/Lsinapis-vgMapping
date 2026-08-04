# Genotyping rules for rescaled BAMs from mapDamage
#Extract parameters for filtering
MIN_DP_RS = config["params"]["run_genotyping"]["minDP"]
MAX_DP_RS = config["params"]["run_genotyping"]["maxDP"]
BASEQ_RS = config["baseQ"]
MAPQ_RS = config["mapQ"]

# Handle samples to drop and keep, needed for merging samples
DROP_SAMPLES_RS = config["params"]["run_genotyping"].get("drop_samples", [])
VALID_SAMPLES_DF_RS = samples_df[~samples_df['sample_id'].isin(DROP_SAMPLES_RS)] #required for the rule filter_missingness
KEEP_SAMPLES_RS = VALID_SAMPLES_DF_RS['sample_id'].unique().tolist()

# DYNAMICALLY FILTER THE POPLIST_RS to ONLY INCLUDE KEPT SAMPLES
# ==============================================================================
ORIGINAL_POPLIST_RS = config["params"]["run_genotyping"]["poplist"] #original unfiltered poplist, which may contain samples we want to drop
FILTERED_POPLIST_RS = "results/genotyping_rescaled/filtered_poplist_rescaled.txt"

# Only attempt to filter if the original file actually exists
if os.path.exists(ORIGINAL_POPLIST_RS):
    # Read the POPLIST_RS (using '\s+' handles both spaces and tabs)
    pop_df = pd.read_csv(ORIGINAL_POPLIST_RS, sep=r'\s+', header=None, names=['sample_id', 'pop_id'])
    
    # Filter the dataframe to only include samples in KEEP_SAMPLES_RS
    filtered_pop_df = pop_df[pop_df['sample_id'].isin(KEEP_SAMPLES_RS)]
    
    # Create the output directory if it doesn't exist and save the filtered file
    os.makedirs(os.path.dirname(FILTERED_POPLIST_RS), exist_ok=True)
    filtered_pop_df.to_csv(FILTERED_POPLIST_RS, sep='\t', header=False, index=False)


# ==============================================================================
# HELPER FUNCTION FOR GENOTYPING INPUT
# ==============================================================================
def get_bam_for_genotyping_rescaled(wildcards):
    """
    Determines the correct BAM file to use for genotyping a given sample.
    Uses the subsampled BAM if the sample is in 'subsample_samples', 
    otherwise uses the full-depth final BAM.
    """
    sid = wildcards.sample_id
    # Look up the source (modern vs historical) for this sample
    src = samples_df[samples_df['sample_id'] == sid]['source'].iloc[0]
    stage = "rescaled" if src == "historical" else "clipped"
    
    subsample_samples = config.get("subsample_samples", [])
    target_dp = config.get("subsample_depth", 15)
    mapq = config["mapQ"]
    
    if sid in subsample_samples:
        # Return the subsampled BAM
        return f"results/mapping/{src}/{sid}.{REF_NAME}.merged.dedup.merged.{stage}.subs{target_dp}.q{mapq}.bam"
    else:
        # Return the standard final BAM
        return f"results/mapping/{src}/{sid}.{REF_NAME}.merged.dedup.merged.{stage}.bam"


def get_bai_for_genotyping_rescaled(wildcards):
    """Simply appends .bai to the dynamically selected BAM file."""
    return get_bam_for_genotyping_rescaled(wildcards) + ".bai"


# ==============================================================================
# HELPER FUNCTIONS FOR JOINT CALLING
# ==============================================================================
def get_all_bams_for_joint_calling_rescaled(wildcards):
    """Gathers the correct BAM files for all kept samples for joint calling."""
    bams = []
    subsample_samples = config.get("subsample_samples", [])
    target_dp = config.get("subsample_depth", 15)
    
    for sid in KEEP_SAMPLES_RS:
        src = VALID_SAMPLES_DF_RS[VALID_SAMPLES_DF_RS['sample_id'] == sid]['source'].iloc[0]
        stage = "rescaled" if src == "historical" else "clipped"
        mapq = config["mapQ"]
        if sid in subsample_samples:
            bams.append(f"results/mapping/{src}/{sid}.{REF_NAME}.merged.dedup.merged.{stage}.subs{target_dp}.q{mapq}.bam")
        else:
            bams.append(f"results/mapping/{src}/{sid}.{REF_NAME}.merged.dedup.merged.{stage}.bam")
    return bams

def get_all_bais_for_joint_calling_rescaled(wildcards):
    """Requires the indices for all joint calling BAMs."""
    return [bam + ".bai" for bam in get_all_bams_for_joint_calling_rescaled(wildcards)]

# ==============================================================================
# GENOTYPING RULES
# ==============================================================================
rule call_individual_genotypes_rescaled:
    """Calls genotypes for a single sample using bcftools.
    Applies filters based on mapping quality, base quality, and depth.
    """
    input:
        bam = get_bam_for_genotyping_rescaled,
        bai = get_bai_for_genotyping_rescaled,
        ref = config["reference"] + ".fa",  # Assumes reference path needs .fa appended
        fai = f"{config['reference']}.fa.fai",
        targets = config["site_filter_bed"]   # List of target regions for genotyping
    output:
        bcf = "results/genotyping_rescaled/individual/{sample_id}.{ref_name}.sitefilt.bQ" + str(BASEQ_RS) + ".mq" + str(MAPQ_RS) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_RS) +"-" +str(MAX_DP_RS)  + ".AB.rescaled.bcf",
        csi = "results/genotyping_rescaled/individual/{sample_id}.{ref_name}.sitefilt.bQ" + str(BASEQ_RS) + ".mq" + str(MAPQ_RS) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_RS) +"-" +str(MAX_DP_RS)  + ".AB.rescaled.bcf.csi"
    params:
        mapq = MAPQ_RS,
        baseq = BASEQ_RS,
        mindp = MIN_DP_RS,
        maxdp = MAX_DP_RS,
    log:
        "logs/genotyping_rescaled/individual/{sample_id}_{ref_name}.log"
    benchmark: "benchmarks/genotyping_rescaled/individual/{sample_id}_{ref_name}.benchmark"
    conda: "../envs/bcftools121.yaml" # Assuming bcftools is in this env, or point to a genotyping.yaml
    threads: 16
    shell:
        """
        bcftools mpileup --threads {threads} -f {input.ref} -R {input.targets} \
            -Ou -B --min-MQ {params.mapq} \
            --min-BQ {params.baseq} -a "FORMAT/AD,FORMAT/DP,INFO/AD" \
            {input.bam} | \
            bcftools call -m -f GQ,GP -Ou | \
            bcftools filter -g 5 -i'QUAL >= 30' -Ou | \
            bcftools view -V indels -Ou | \
            bcftools +setGT -Ou -- -t q -n . -i"FMT/DP<{params.mindp} | FMT/DP>{params.maxdp}" | \
            bcftools +setGT -Ou -- -t q -n . \
                -i'GT="het" & (FMT/AD[:0]/FMT/DP < 0.21 | FMT/AD[:0]/FMT/DP > .79 | FMT/AD[:1]/FMT/DP < 0.21 | FMT/AD[:1]/FMT/DP > .79)' | \
            bcftools view -i'GT!="./."' -Ou | \
            bcftools +fill-tags -Ob -- -t all > {output.bcf} 2> {log}
        bcftools index -o {output.csi} {output.bcf} 2>> {log}
        """

rule merge_genotypes_rescaled:
    """1. Merges valid individual BCF files, excluding dropped samples as specified in config. 
      2. Generates a merged BCF file for all samples."""
    input:
        bcfs = lambda wildcards: expand(
            "results/genotyping_rescaled/individual/{sample_id}.{ref_name}.sitefilt.bQ" + str(BASEQ_RS) + ".mq" + str(MAPQ_RS) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_RS) + "-" + str(MAX_DP_RS) + ".AB.rescaled.bcf", 
            sample_id=KEEP_SAMPLES_RS,
            ref_name=REF_NAME 
        ),
        csis = lambda wildcards: expand(
            "results/genotyping_rescaled/individual/{sample_id}.{ref_name}.sitefilt.bQ" + str(BASEQ_RS) + ".mq" + str(MAPQ_RS) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_RS) + "-" + str(MAX_DP_RS) + ".AB.rescaled.bcf.csi", 
            sample_id=KEEP_SAMPLES_RS,
            ref_name=REF_NAME  
        )
    output:
        # You can keep the merged output name shorter, or make it match the long string. 
        # I recommend keeping it long so you know exactly what filters were applied to the merged file!
        merged_bcf = "results/genotyping_rescaled/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ_RS) + ".mq" + str(MAPQ_RS) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_RS) + "-" + str(MAX_DP_RS) + ".AB.indCall.rescaled.allsites.bcf",
        csi = "results/genotyping_rescaled/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ_RS) + ".mq" + str(MAPQ_RS) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_RS) + "-" + str(MAX_DP_RS) + ".AB.indCall.rescaled.allsites.bcf.csi",
        stats = "results/genotyping_rescaled/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ_RS) + ".mq" + str(MAPQ_RS) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_RS) + "-" + str(MAX_DP_RS) + ".AB.indCall.rescaled.allsites.bcf.stats",
    log:
        "logs/genotyping_rescaled/merge_genotypes_{ref_name}.log"
    benchmark: "benchmarks/genotyping_rescaled/merge_genotypes_{ref_name}.benchmark"
    conda: "../envs/bcftools121.yaml"
    threads: 4
    shell:
        """
        bcftools merge --force-samples -Ou {input.bcfs} | \
            bcftools +fill-tags -Ou -- -t all | \
            bcftools filter -g 5 -Ou | \
            bcftools view -M2 -V indels -Ob -o {output.merged_bcf} 2> {log}
        bcftools index -o {output.csi} {output.merged_bcf} 2>> {log}
        bcftools stats -s - {output.merged_bcf} > {output.stats}
        """

rule joint_call_genotypes_rescaled:
    """Calls genotypes jointly across all valid samples."""
    input:
        bams = get_all_bams_for_joint_calling_rescaled,
        bais = get_all_bais_for_joint_calling_rescaled,
        ref = config["reference"] + ".fa",
        fai = f"{config['reference']}.fa.fai",
        targets=config["site_filter_bed"] ,  # List of target regions for genotyping
        poplist = FILTERED_POPLIST_RS #for HWE filtering
    output:
        bcf = "results/genotyping_rescaled/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ_RS) + ".mq" + str(MAPQ_RS) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_RS) + "-" + str(MAX_DP_RS) + ".AB.jointCall.rescaled.allsites.bcf",
        csi = "results/genotyping_rescaled/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ_RS) + ".mq" + str(MAPQ_RS) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_RS) + "-" + str(MAX_DP_RS) + ".AB.jointCall.rescaled.allsites.bcf.csi",
        stats = "results/genotyping_rescaled/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ_RS) + ".mq" + str(MAPQ_RS) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_RS) + "-" + str(MAX_DP_RS) + ".AB.jointCall.rescaled.allsites.bcf.stats",
    params:
        MIN_DP_RS = MIN_DP_RS,
        MAX_DP_RS = MAX_DP_RS,
        mapq = MAPQ_RS,
        baseq = BASEQ_RS,
    log:
        "logs/genotyping_rescaled/joint_call_genotypes_{ref_name}.log"
    benchmark: "benchmarks/genotyping_rescaled/joint_call_genotypes_{ref_name}.benchmark"
    conda: "../envs/bcftools121.yaml"
    threads: 8
    shell:
        """
        bcftools mpileup --threads {threads} -f {input.ref} -T {input.targets} \
            -a "FORMAT/AD,FORMAT/DP,INFO/AD" -B \
            --min-MQ {params.mapq} --min-BQ {params.baseq} -Ou {input.bams} | \
            bcftools call -m -a GQ,GP -G {input.poplist} -Ou | \
            bcftools filter -g 5 -i'QUAL >= 30' -Ou | \
            bcftools view -V indels -M2 -Ou | \
            bcftools +fill-tags -Ou -- -t all | \
            bcftools +setGT -Ou -- -t q -n . -i"FMT/DP<{params.MIN_DP_RS} | FMT/DP>{params.MAX_DP_RS}" | \
            bcftools +setGT -Ou -- -t q -n . -i'GT="het" & (FMT/VAF < 0.21 | FMT/VAF > 0.79)'  | \
            bcftools +fill-tags -Ob -- -t all > {output.bcf} 2> {log}
            
        bcftools index -o {output.csi} {output.bcf} 2>> {log}
        bcftools stats -s - {output.bcf} > {output.bcf}.stats
        """


rule filter_variants_rescaled:
    """Filters variant + invariant site bcf for only variant and biallelic SNP sites."""
    input:
        bcf = "results/genotyping_rescaled/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ_RS) + ".mq" + str(MAPQ_RS) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_RS) + "-" + str(MAX_DP_RS) + ".AB.{call_type}.rescaled.allsites.bcf",
        csi = "results/genotyping_rescaled/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ_RS) + ".mq" + str(MAPQ_RS) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_RS) + "-" + str(MAX_DP_RS) + ".AB.{call_type}.rescaled.allsites.bcf.csi"
    output:
        bcf = "results/genotyping_rescaled/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ_RS) + ".mq" + str(MAPQ_RS) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_RS) + "-" + str(MAX_DP_RS) + ".AB.{call_type}.rescaled.biallelic.bcf",
        csi = "results/genotyping_rescaled/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ_RS) + ".mq" + str(MAPQ_RS) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_RS) + "-" + str(MAX_DP_RS) + ".AB.{call_type}.rescaled.biallelic.bcf.csi",
        stats = "results/genotyping_rescaled/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ_RS) + ".mq" + str(MAPQ_RS) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_RS) + "-" + str(MAX_DP_RS) + ".AB.{call_type}.rescaled.biallelic.bcf.stats"
    log:
        "logs/genotyping_rescaled/filter_variants_{ref_name}_{call_type}.log"
    conda: "../envs/bcftools121.yaml"
    shell:
        """
        bcftools view -v snps -m 2 -M 2 -i 'MAF>0' -Ob -o {output.bcf} {input.bcf} 2> {log}
        bcftools index -o {output.csi} {output.bcf} 2>> {log}
        bcftools stats -s - {output.bcf} > {output.stats} 
        """

rule filter_missingness_rescaled:
    """Filters sites based on missingness threshold {max_miss} for either allsites or variants."""
    input:
        bcf = "results/genotyping_rescaled/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ_RS) + ".mq" + str(MAPQ_RS) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_RS) + "-" + str(MAX_DP_RS) + ".AB.{call_type}.rescaled.{site_type}.bcf",
        csi = "results/genotyping_rescaled/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ_RS) + ".mq" + str(MAPQ_RS) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_RS) + "-" + str(MAX_DP_RS) + ".AB.{call_type}.rescaled.{site_type}.bcf.csi"
    output:
        bcf = "results/genotyping_rescaled/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ_RS) + ".mq" + str(MAPQ_RS) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_RS) + "-" + str(MAX_DP_RS) + ".AB.{call_type}.rescaled.{site_type}.fmiss{max_miss}.bcf",
        csi = "results/genotyping_rescaled/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ_RS) + ".mq" + str(MAPQ_RS) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_RS) + "-" + str(MAX_DP_RS) + ".AB.{call_type}.rescaled.{site_type}.fmiss{max_miss}.bcf.csi",
        stats = "results/genotyping_rescaled/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ_RS) + ".mq" + str(MAPQ_RS) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_RS) + "-" + str(MAX_DP_RS) + ".AB.{call_type}.rescaled.{site_type}.fmiss{max_miss}.bcf.stats",
        modsites = temp("results/genotyping_rescaled/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ_RS) + ".mq" + str(MAPQ_RS) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_RS) + "-" + str(MAX_DP_RS) + ".AB.{call_type}.rescaled.{site_type}.fmiss{max_miss}.modsites.tmp"),
        histsites = temp("results/genotyping_rescaled/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ_RS) + ".mq" + str(MAPQ_RS) + ".snps5.noIndel.Q30.dp" + str(MIN_DP_RS) + "-" + str(MAX_DP_RS) + ".AB.{call_type}.rescaled.{site_type}.fmiss{max_miss}.histsites.tmp")
    params:
        modsamps = ",".join(VALID_SAMPLES_DF_RS[VALID_SAMPLES_DF_RS['source'] == 'modern']['sample_id'].unique().tolist()),
        histsamps = ",".join(VALID_SAMPLES_DF_RS[VALID_SAMPLES_DF_RS['source'] == 'historical']['sample_id'].unique().tolist())
    log:
        "logs/genotyping_rescaled/filter_missingness_{ref_name}_{call_type}_{site_type}_fmiss{max_miss}.log"
    conda: "../envs/bcftools121.yaml"
    threads: 4
    shell:
        """
        bcftools view -s {params.histsamps} --force-samples -Ou {input.bcf} | \\
            bcftools +fill-tags -Ou -- -t F_MISSING | \\
            bcftools query -i 'F_MISSING <= {wildcards.max_miss}' \\
                -f '%CHROM\\t%POS\\n' > {output.histsites} 2> {log}
        
        bcftools view -s {params.modsamps} --force-samples -Ou {input.bcf} | \\
            bcftools +fill-tags -Ou -- -t F_MISSING | \\
            bcftools query -i 'F_MISSING <= {wildcards.max_miss}' \\
                -f '%CHROM\\t%POS\\n' > {output.modsites} 2>> {log}
        
        bcftools view -T {output.histsites} -Ou {input.bcf} | \\
            bcftools view -T {output.modsites} -Ob -o {output.bcf} 2>> {log}
        
        bcftools index -o {output.csi} {output.bcf}
        bcftools stats -s - {output.bcf} > {output.stats}
        """

rule bcf2vcf_rescaled:
    """
    Makes a bcf a vcf when needed.
    """
    input:
        bcf="results/genotyping_rescaled/{prefix}.bcf",
        idx="results/genotyping_rescaled/{prefix}.bcf.csi",
    output:
        vcf="results/genotyping_rescaled/{prefix}.vcf.gz",
        tbi="results/genotyping_rescaled/{prefix}.vcf.gz.tbi",
    conda:
        "../envs/bcftools121.yaml"
    threads: 6
    shell:
        """
        bcftools view -Oz {input.bcf} > {output.vcf}
        tabix {output.vcf}
        """

rule bcf_ref_bias_rescaled:
    """
    Calculate reference bias (ref alleles / total alleles) per sample from calls
    """
    input:
        stats="results/genotyping_rescaled/{prefix}.bcf.stats",
    output:
        bias="results/genotyping_rescaled/{prefix}.bcf.stats.ref_bias",
    shell:
        """
        grep PSC {input.stats} | \
            grep -v "#" | \
            awk '{{print $3"\t"(2*$4+$6)/(2*($4+$5+$6))}}' > {output.bias}
        """