#Extract parameters for filtering
MIN_DP = config["params"]["call_genotypes"]["minDP"]
MAX_DP = config["params"]["call_genotypes"]["maxDP"]
BASEQ = config["baseQ"]
MAPQ = config["mapQ"]
POPLIST = config["parameters"]["call_genotypes"]["poplist"] #required for joint calling, can be a simple text file with sample_id and population columns 

# Handle samples to drop and keep, needed for merging samples
DROP_SAMPLES = config["parameters"]["call_genotypes"].get("drop_samples", [])
VALID_SAMPLES_DF = samples_df[~samples_df['sample_id'].isin(DROP_SAMPLES)] #required for the rule filter_missingness
KEEP_SAMPLES = VALID_SAMPLES_DF['sample_id'].tolist()

# ==============================================================================
# HELPER FUNCTION FOR GENOTYPING INPUT
# ==============================================================================
def get_bam_for_genotyping(wildcards):
    """
    Determines the correct BAM file to use for genotyping a given sample.
    Uses the subsampled BAM if the sample is in 'subsample_samples', 
    otherwise uses the full-depth final BAM.
    """
    sid = wildcards.sample_id
    # Look up the source (modern vs historical) for this sample
    src = samples_df[samples_df['sample_id'] == sid]['source'].iloc[0]
    stage = "masked" if src == "historical" else "clipped"
    
    subsample_samples = config.get("subsample_samples", [])
    target_dp = config.get("subsample_depth", 15)
    
    if sid in subsample_samples:
        # Return the subsampled BAM
        return f"results/mapping/{src}/{sid}.{REF_NAME}.merged.dedup.merged.{stage}.subs{target_dp}.sitefilt.Q20.q30.bam"
    else:
        # Return the standard final BAM
        return f"results/mapping/{src}/{sid}.{REF_NAME}.merged.dedup.merged.{stage}.bam"


def get_bai_for_genotyping(wildcards):
    """Simply appends .bai to the dynamically selected BAM file."""
    return get_bam_for_genotyping(wildcards) + ".bai"


# ==============================================================================
# HELPER FUNCTIONS FOR JOINT CALLING
# ==============================================================================
def get_all_bams_for_joint_calling(wildcards):
    """Gathers the correct BAM files for all kept samples for joint calling."""
    bams = []
    subsample_samples = config.get("subsample_samples", [])
    target_dp = config.get("subsample_depth", 15)
    
    for sid in KEEP_SAMPLES:
        src = VALID_SAMPLES_DF[VALID_SAMPLES_DF['sample_id'] == sid]['source'].iloc[0]
        stage = "masked" if src == "historical" else "clipped"
        
        if sid in subsample_samples:
            bams.append(f"results/mapping/{src}/{sid}.{REF_NAME}.merged.dedup.merged.{stage}.subs{target_dp}.sitefilt.Q20.q30.bam")
        else:
            bams.append(f"results/mapping/{src}/{sid}.{REF_NAME}.merged.dedup.merged.{stage}.bam")
    return bams

def get_all_bais_for_joint_calling(wildcards):
    """Requires the indices for all joint calling BAMs."""
    return [bam + ".bai" for bam in get_all_bams_for_joint_calling(wildcards)]

# ==============================================================================
# GENOTYPING RULES
# ==============================================================================
rule call_individual_genotypes:
    """Calls genotypes for a single sample using bcftools.
    Applies filters based on mapping quality, base quality, and depth.
    """
    input:
        bam = get_bam_for_genotyping,
        bai = get_bai_for_genotyping,
        ref = config["reference"] + ".fa",  # Assumes reference path needs .fa appended
        fai = f"{config['reference']}.fa.fai",
        sites=config["snp_mask_bed"]  # List of target regions for genotyping
    output:
        bcf = "results/genotyping/individual/{sample_id}.{ref_name}.sitefilt.bQ" + str(BASEQ) + ".mq" + str(MAPQ) + ".snps5.noIndel.Q30.dp" + str(MIN_DP) +"-" +str(MAX_DP)  + ".AB.bcf",
        csi = "results/genotyping/individual/{sample_id}.{ref_name}.sitefilt.bQ" + str(BASEQ) + ".mq" + str(MAPQ) + ".snps5.noIndel.Q30.dp" + str(MIN_DP) +"-" +str(MAX_DP)  + ".AB.bcf.csi"
    params:
        mapq = MAPQ,
        baseq = BASEQ,
        mindp = MIN_DP,
        maxdp = MAX_DP,
    conda: "../envs/bcftools121.yaml" # Assuming bcftools is in this env, or point to a genotyping.yaml
    threads: 16
    shell:
        """
        bcftools mpileup --threads {threads} -f {input.ref} -T {input.sites} \
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
            bcftools +fill-tags -Ob -- -t all > {output.bcf}
        bcftools index -o {output.csi} {output.bcf}
        """

rule merge_genotypes:
    """1. Merges valid individual BCF files, excluding dropped samples as specified in config. 
      2. Generates a merged BCF file for all samples."""
    input:
        bcfs = lambda wildcards: expand(
            "results/genotyping/individual/{sample_id}.{ref_name}.sitefilt.bQ" + str(BASEQ) + ".mq" + str(MAPQ) + ".snps5.noIndel.Q30.dp" + str(MIN_DP) + "-" + str(MAX_DP) + ".AB.bcf", 
            sample_id=KEEP_SAMPLES,
            ref_name=REF_NAME 
        ),
        csis = lambda wildcards: expand(
            "results/genotyping/individual/{sample_id}.{ref_name}.sitefilt.bQ" + str(BASEQ) + ".mq" + str(MAPQ) + ".snps5.noIndel.Q30.dp" + str(MIN_DP) + "-" + str(MAX_DP) + ".AB.bcf.csi", 
            sample_id=KEEP_SAMPLES,
            ref_name=REF_NAME  
        )
    output:
        # You can keep the merged output name shorter, or make it match the long string. 
        # I recommend keeping it long so you know exactly what filters were applied to the merged file!
        merged_bcf = "results/genotyping/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ) + ".mq" + str(MAPQ) + ".snps5.noIndel.Q30.dp" + str(MIN_DP) + "-" + str(MAX_DP) + ".AB.indCall.allsites.bcf",
        csi = "results/genotyping/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ) + ".mq" + str(MAPQ) + ".snps5.noIndel.Q30.dp" + str(MIN_DP) + "-" + str(MAX_DP) + ".AB.indCall.allsites.bcf.csi"
        stats = "results/genotyping/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ) + ".mq" + str(MAPQ) + ".snps5.noIndel.Q30.dp" + str(MIN_DP) + "-" + str(MAX_DP) + ".AB.indCall.allsites.bcf.stats"
    conda: "../envs/bcftools121.yaml"
    threads: 4
    shell:
        """
        bcftools merge --force-samples -Ou {input.bcfs} | \
            bcftools +fill-tags -Ou -- -t all | \
            bcftools filter -g 5 -Ou | \
            bcftools view -M2 -V indels -Ob -o {output.merged_bcf}
        bcftools index -o {output.csi} {output.merged_bcf}
        bcftools stats -s - {output.merged_bcf} > {output.stats}
        """

rule joint_call_genotypes:
    """Calls genotypes jointly across all valid samples."""
    input:
        bams = get_all_bams_for_joint_calling,
        bais = get_all_bais_for_joint_calling,
        ref = config["reference"] + ".fa",
        fai = f"{config['reference']}.fa.fai",
        sites=config["snp_mask_bed"],  # List of target regions for genotyping
        poplist = POPLIST #for HWE filtering
    output:
        bcf = "results/genotyping/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ) + ".mq" + str(MAPQ) + ".snps5.noIndel.Q30.dp" + str(MIN_DP) + "-" + str(MAX_DP) + ".AB.jointCall.allsites.bcf",
        csi = "results/genotyping/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ) + ".mq" + str(MAPQ) + ".snps5.noIndel.Q30.dp" + str(MIN_DP) + "-" + str(MAX_DP) + ".AB.jointCall.allsites.bcf.csi"
        stats = "results/genotyping/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ) + ".mq" + str(MAPQ) + ".snps5.noIndel.Q30.dp" + str(MIN_DP) + "-" + str(MAX_DP) + ".AB.jointCall.allsites.bcf.stats"
    params:
        min_dp = MIN_DP,
        max_dp = MAX_DP,
        mapq = MAPQ,
        baseq = BASEQ
    conda: "../envs/bcftools121.yaml"
    threads: 8
    resources:
        runtime="8h"
    shell:
        """
        bcftools mpileup --threads {threads} -f {input.ref} -T {input.sites} \\
            -a "FORMAT/AD,FORMAT/DP,INFO/AD" -B \\
            --min-MQ {params.mapq} --min-BQ {params.baseq} -Ou {input.bams} | \\
            bcftools call -m -a GQ,GP -G {input.poplist} -Ou | \\
            bcftools filter -g 5 -i'QUAL >= 30' -Ou | \\
            bcftools view -V indels -M2 -Ou | \\
            bcftools +fill-tags -Ou -- -t all | \\
            bcftools +setGT -Ou -- -t q -n . -i"FMT/DP<{params.min_dp} | FMT/DP>{params.max_dp}" | \\
            bcftools +setGT -Ou -- -t q -n . -i'GT="het" & (FMT/VAF < 0.21 | FMT/VAF > 0.79)'  | \\
            bcftools +fill-tags -Ob -- -t all > {output.bcf}
            
        bcftools index -o {output.csi} {output.bcf}
        bcftools stats -s - {output.bcf} > {output.bcf}.stats
        """


rule filter_variants:
    """Filters variant + invariant site bcf for only variant and biallelic SNP sites."""
    input:
        bcf = "results/genotyping/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ) + ".mq" + str(MAPQ) + ".snps5.noIndel.Q30.dp" + str(MIN_DP) + "-" + str(MAX_DP) + ".AB.{call_type}.allsites.bcf",
        csi = "results/genotyping/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ) + ".mq" + str(MAPQ) + ".snps5.noIndel.Q30.dp" + str(MIN_DP) + "-" + str(MAX_DP) + ".AB.{call_type}.allsites.bcf.csi"
    output:
        bcf = "results/genotyping/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ) + ".mq" + str(MAPQ) + ".snps5.noIndel.Q30.dp" + str(MIN_DP) + "-" + str(MAX_DP) + ".AB.{call_type}.biallelic.bcf",
        csi = "results/genotyping/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ) + ".mq" + str(MAPQ) + ".snps5.noIndel.Q30.dp" + str(MIN_DP) + "-" + str(MAX_DP) + ".AB.{call_type}.biallelic.bcf.csi",
        stats = "results/genotyping/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ) + ".mq" + str(MAPQ) + ".snps5.noIndel.Q30.dp" + str(MIN_DP) + "-" + str(MAX_DP) + ".AB.{call_type}.biallelic.bcf.stats"
    conda: "../envs/bcftools121.yaml"
    shell:
        """
        bcftools view -v snps -m 2 -M 2 -i 'MAF>0' -Ob -o {output.bcf} {input.bcf}
        bcftools index -o {output.csi} {output.bcf}
        bcftools stats -s - {output.bcf} > {output.stats}
        """

rule filter_missingness:
    """Filters sites based on missingness threshold {max_miss} for either allsites or variants."""
    input:
        bcf = "results/genotyping/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ) + ".mq" + str(MAPQ) + ".snps5.noIndel.Q30.dp" + str(MIN_DP) + "-" + str(MAX_DP) + ".AB.{call_type}.{site_type}.bcf",
        csi = "results/genotyping/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ) + ".mq" + str(MAPQ) + ".snps5.noIndel.Q30.dp" + str(MIN_DP) + "-" + str(MAX_DP) + ".AB.{call_type}.{site_type}.bcf.csi"
    output:
        bcf = "results/genotyping/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ) + ".mq" + str(MAPQ) + ".snps5.noIndel.Q30.dp" + str(MIN_DP) + "-" + str(MAX_DP) + ".AB.{call_type}.{site_type}.fmiss{max_miss}.bcf",
        csi = "results/genotyping/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ) + ".mq" + str(MAPQ) + ".snps5.noIndel.Q30.dp" + str(MIN_DP) + "-" + str(MAX_DP) + ".AB.{call_type}.{site_type}.fmiss{max_miss}.bcf.csi",
        stats = "results/genotyping/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ) + ".mq" + str(MAPQ) + ".snps5.noIndel.Q30.dp" + str(MIN_DP) + "-" + str(MAX_DP) + ".AB.{call_type}.{site_type}.fmiss{max_miss}.bcf.stats",
        modsites = temp("results/genotyping/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ) + ".mq" + str(MAPQ) + ".snps5.noIndel.Q30.dp" + str(MIN_DP) + "-" + str(MAX_DP) + ".AB.{call_type}.{site_type}.fmiss{max_miss}.modsites.tmp"),
        histsites = temp("results/genotyping/merged.all.{ref_name}.sitefilt.bQ" + str(BASEQ) + ".mq" + str(MAPQ) + ".snps5.noIndel.Q30.dp" + str(MIN_DP) + "-" + str(MAX_DP) + ".AB.{call_type}.{site_type}.fmiss{max_miss}.histsites.tmp")
    params:
        modsamps = ",".join(VALID_SAMPLES_DF[VALID_SAMPLES_DF['source'] == 'modern']['sample_id'].tolist()),
        histsamps = ",".join(VALID_SAMPLES_DF[VALID_SAMPLES_DF['source'] == 'historical']['sample_id'].tolist())
    conda: "../envs/bcftools121.yaml"
    threads: 4
    shell:
        """
        bcftools view -s {params.histsamps} --force-samples -Ou {input.bcf} | \\
            bcftools +fill-tags -Ou -- -t F_MISSING | \\
            bcftools query -i 'F_MISSING <= {wildcards.max_miss}' \\
                -f '%CHROM\\t%POS\\n' > {output.histsites}
        
        bcftools view -s {params.modsamps} --force-samples -Ou {input.bcf} | \\
            bcftools +fill-tags -Ou -- -t F_MISSING | \\
            bcftools query -i 'F_MISSING <= {wildcards.max_miss}' \\
                -f '%CHROM\\t%POS\\n' > {output.modsites}
        
        bcftools view -T {output.histsites} -Ou {input.bcf} | \\
            bcftools view -T {output.modsites} -Ob -o {output.bcf}
        
        bcftools index -o {output.csi} {output.bcf}
        bcftools stats -s - {output.bcf} > {output.stats}
        """

rule bcf2vcf:
    """
    Makes a bcf a vcf when needed.
    """
    input:
        bcf="results/{geno_dir}/{prefix}.bcf",
        idx="results/{geno_dir}/{prefix}.bcf.csi",
    output:
        vcf="results/{geno_dir}/{prefix}.vcf.gz",
        tbi="results/{geno_dir}/{prefix}.vcf.gz.tbi",
    wildcard_constraints:
        geno_dir="genotyping|genotyping_notrans"
    conda:
        "../envs/bcftools121.yaml"
    threads: 6
    shell:
        """
        bcftools view -Oz {input.bcf} > {output.vcf}
        tabix {output.vcf}
        """

rule bcf_ref_bias:
    """
    Calculate reference bias (ref alleles / total alleles) per sample from calls
    """
    input:
        stats="results/{geno_dir}/{prefix}.bcf.stats",
    output:
        bias="results/{geno_dir}/{prefix}.bcf.stats.ref_bias",
    wildcard_constraints:
        geno_dir="genotyping|genotyping_notrans"
    shell:
        """
        grep PSC {input.stats} | \
            grep -v "#" | \
            awk '{{print $3"\t"(2*$4+$6)/(2*($4+$5+$6))}}' > {output.bias}
        """