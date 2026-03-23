
# Extract the target depth from the config file once at the top for naming outputs
TARGET_DP = config.get("subsample_depth", 15)

rule subsample_dedup_bam:
    """Subsamples deduplicated BAM files based on target depth defined in config."""
    input:
        bam = "results/mapping/{source}/{sample_id}.{ref_name}.merged.dedup.merged.bam",
        depth = "results/mapping/{source}/stats/merged_dedup_merged/{sample_id}.{ref_name}.merged.dedup.merged.regfilt.Q20.q30.depth.txt"
    output:
        bam = "results/mapping/{source}/{sample_id}.{ref_name}.merged.dedup.merged.subs" + str(TARGET_DP) +".regfilt.Q20.q30.bam",
        bai = "results/mapping/{source}/{sample_id}.{ref_name}.merged.dedup.merged.subs" + str(TARGET_DP) + ".regfilt.Q20.q30.bam.bai"
    params:
        target_depth = config.get("subsample_depth", 15),
        mapq = config.get("mapQ", 30)
    log: "logs/mapping/subsample_dedup/{source}/{sample_id}.{ref_name}.subsample_dedup.log"
    conda: "../envs/vg.yaml"
    threads: 4
    shell:
        """
        MEAN_DEPTH=$(cat {input.depth})
        TARGET_DEPTH={params.target_depth}

        # Calculate fraction using awk (returns 1.0 if mean depth is lower than target)
        FRACTION=$(awk -v mean=$MEAN_DEPTH -v target=$TARGET_DEPTH 'BEGIN {{if(mean>target) print target/mean; else print 1.0}}')

        # Check if we need to downsample (fraction < 1.0)
        if awk -v frac=$FRACTION 'BEGIN {{if(frac < 1.0) exit 0; else exit 1}}'; then
            # samtools view -s uses the integer part as the random seed (42 here) 
            # and the decimal part as the percentage to keep.
            SEED_FRACTION=$(awk -v frac=$FRACTION 'BEGIN {{print 42 + frac}}')
            samtools view -h -F 4 --min-MQ {params.mapq} -@ {threads} -u {input.bam} |
              samtools view -@ {threads} -h -s $SEED_FRACTION -b > {output.bam} 2> {log}
        else
            ln -sf {input.bam} {output.bam} 2>> {log}
        fi
        
        samtools index -@ {threads} {output.bam} 2>> {log}
        """

rule subsample_final_bam:
    """Subsamples final masked/clipped BAM files."""
    input:
        bam = "results/mapping/{source}/{sample_id}.{ref_name}.merged.dedup.merged.{stage}.bam",
        depth = "results/mapping/{source}/stats/merged_dedup_merged_{stage}/{sample_id}.{ref_name}.merged.dedup.merged.{stage}.regfilt.Q20.q30.depth.txt"
    output:
        bam = "results/mapping/{source}/{sample_id}.{ref_name}.merged.dedup.merged.{stage}.subs" + str(TARGET_DP) + ".regfilt.Q20.q30.bam",
        bai = "results/mapping/{source}/{sample_id}.{ref_name}.merged.dedup.merged.{stage}.subs" + str(TARGET_DP) + ".regfilt.Q20.q30.bam.bai"
    params:
        target_depth = config.get("subsample_depth", 15),
        mapq = config.get("mapQ", 30),
    log: "logs/mapping/subsample_final/{source}/{sample_id}.{ref_name}.{stage}.subsample_final.log"
    conda: "../envs/vg.yaml"
    threads: 4
    shell:
        """
        MEAN_DEPTH=$(cat {input.depth})
        TARGET_DEPTH={params.target_depth}

        FRACTION=$(awk -v mean=$MEAN_DEPTH -v target=$TARGET_DEPTH 'BEGIN {{if(mean>target) print target/mean; else print 1.0}}')

        if awk -v frac=$FRACTION 'BEGIN {{if(frac < 1.0) exit 0; else exit 1}}'; then
            SEED_FRACTION=$(awk -v frac=$FRACTION 'BEGIN {{print 42 + frac}}')
             samtools view -h -F 4 --min-MQ {params.mapq} -@ {threads} -u {input.bam} |
              samtools view -@ {threads} -h -s $SEED_FRACTION -b > {output.bam} 2> {log}
        else
            ln -sf {input.bam} {output.bam} 2>> {log}
        fi
        
        samtools index -@ {threads} {output.bam} 2>> {log}
        """