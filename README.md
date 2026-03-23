## A snakemake pipeline for calling genotypes from modern and historical specimens using a sequence variation graph

This pipeline borrows heavily from https://github.com/zjnolen/polyommatini-temporal-genomics/tree/main  

The following steps were performed prior to running this pipeline

1) Used the [generode pipeline](https://github.com/NBISweden/GenErode) to generate bam files for modern specimens using `BWA mem` and historical ones using `BWA aln`
2) Called high quality genotypes from resulting modern BAM files above and used `vg` to create a sequence variation graph. This was conducted using `utlitlites/vg_construct.sh`. Additionally a catalog of the snps identified was also constructed for fitlering in this pipeline
3) Site fitlers were constucted using repeat masker, genmap (for mappability) and additional regionss that might be problematic for mapping.

