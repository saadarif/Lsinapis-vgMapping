# Lsinapis-vgMapping

Snakemake pipeline for population genomics of modern + historical *Leptidea
sinapis* specimens, mapped to a sequence variation graph. Built incrementally,
one new workflow stage per work session — see "Current status" below before
starting new work, and update it as stages are added.

## Running it

```bash
cd /mnt/sda/saad/WW_Trimmed_Reads_WWPaper/Lsinapis-vgMapping
snakemake -n -p        # dry run, full rendered shell commands
snakemake              # actual run
```

- No bare `snakemake` on `PATH` by default on this machine — use the
  `vg_mapping_smk_env` conda env (Snakemake 9.16.2) if it isn't already active.
- `profiles/default/config.yaml` is auto-detected from the repo root: sets
  `cores: 16`, `rerun-triggers: mtime` (edits to a rule alone do NOT trigger a
  rerun — force with `-R <rulename>`), and
  `software-deployment-method: [conda, apptainer]`.
- Background long runs: `nohup snakemake > snakemake_run.log 2>&1 &`.

## Layout

- `Snakefile` — loads `config.yaml`, builds `samples_df` from
  `modern_csv`/`historical_csv`, `include:`s every `workflow/rules/*.smk`, and
  assembles `get_final_targets()` as a numbered list of sections, one per
  analysis toggle (`run_<name>: TRUE/FALSE` in `config.yaml`). Each rule file
  builds its own `<NAME>_TARGETS` list; the Snakefile just extends into it
  when the toggle is on. Follow this pattern for new stages rather than
  inlining target logic in the Snakefile.
- `workflow/rules/*.smk` — one file per pipeline stage, numbered in
  dependency order: `1.1_mapping`, `1.2_subsampling`, `2a_call_genotypes_noTrans`,
  `2b_call_genotypes_rescaled`, `3_diversity_stats`, `4_relatedness`,
  `5_LD_pruning`, `6a_Strucure_angsdPCA` (typo in filename is intentional/
  established, don't "fix" it without renaming deliberately),
  `6b_Structure_NGadmix`, `7_FroH_bcftools`, `10a_Wolbachia_screen`,
  `10b_InStrain_Wolbachia_screen`. All `include:`d
  into one global namespace — later files freely reuse Python-level constants
  and functions from earlier ones (e.g. `REF_NAME`, `samples_df`,
  `KEEP_SAMPLES_REL`, `REL_TAG`, `PRUNE_TAG`).
- `workflow/envs/*.yaml` — one conda env per rule (or small group), pinned
  versions. Even coreutils-only rules (`grep`/`awk`/`cut`) skip `conda:`
  entirely rather than reusing an unrelated env — see
  `individual_heterozygosity` in `3_diversity_stats.smk`.
- `workflow/scripts/` — Python/R/bash scripts invoked via `script:`.
- `utilities/` — standalone helper scripts run by hand from the CLI, not
  wired into Snakemake (no rule, toggle or env). Includes
  `plot_pixy_modern_windows.R` (see session 2026-10-12 below).
- `workflow/schemas/config.schema.yaml` — JSON schema for `config.yaml`
  (every toggle + every `params: run_<name>` block). Not wired in: nothing
  calls `validate()` on it. Add the new toggle/params block here whenever a
  stage is added.
- `config.yaml` — one `run_<name>: TRUE/FALSE` toggle per stage under
  `# --- Analysis Toggles ---`, plus a `params: run_<name>: {...}` block per
  stage with defaults read via `.get()` (never hard `config["params"][...]`
  for new stage-specific params — keep old configs from breaking when a
  stage is added).

## Conventions worth knowing before adding a stage

- **Tag-chaining file names.** Each stage's output tag folds in every
  upstream tag plus its own params, so file names alone are provenance:
  `REL_TAG` (4_relatedness: bam_stage/site filter/bQ/mq/depth/minInd/maf/pval)
  → `LD_TAG` (adds `.maxkb{N}`) → `PRUNE_TAG` (adds `.prunebp{N}.minr2{N}`) →
  `PCA_TAG` (`excl-none` or `excl-<sorted-sample-ids>`, plus
  `.exclpop-<sorted-populations>` only when populations are excluded) → `ADMIX_TAG`
  (= `PRUNE_TAG.PCA_TAG`, K goes in the file name as `.K{kvalue}`). A new stage
  consuming an existing output should depend on the upstream tag/prefix
  variable directly (e.g. `PRUNE_TAG`), not re-derive or hardcode it. Log and benchmark files should also contain tags to they can be differentiated from different runs of the rule
- **conda: vs container:** default to `conda:` with a pinned bioconda/
  conda-forge env. Only reach for `container:` when the tool genuinely has no
  bioconda/conda-forge package (checked via `bioconda-recipes` GitHub
  listing, not just `mamba search` — search past cases: `ngsLD`,
  `prune_graph` and `evalAdmix` have neither). When a `container:` rule is added,
  `apptainer` must already be enabled in `profiles/default/config.yaml`
  (it is, as of `5_LD_pruning.smk`) — no further config needed as long as
  the container only touches paths under the repo working directory
  (apptainer auto-binds `$PWD`).
- **Validate real commands against real data before wiring into a rule.**
  The working pattern so far: build a small real slice of the actual
  upstream output (not synthetic data), run the tool's exact CLI by hand,
  confirm output shape/format, *then* write the rule and confirm with
  `snakemake -n -p` that the rendered command matches what was hand-tested.
  Don't re-litigate this for a tool already proven this way (e.g. don't
  re-pull/re-test the `ngsld` container for a new rule that reuses it).
- **Script rules** (`script:`) get their own conda env even if only stdlib
  is needed (see `workflow/envs/python.yaml`), for reproducibility
  consistency with every other rule — but avoid adding a dependency (e.g.
  pandas) an existing script doesn't need.

## Current status (update this section as stages land)

Built and wired into `Snakefile`/`config.yaml`, in dependency order:

1. `1.1_mapping.smk`, `1.2_subsampling.smk` — mapping to the variation graph
   + optional depth subsampling
2. `2a_call_genotypes_noTrans.smk`, `2b_call_genotypes_rescaled.smk` —
   bcftools-based genotyping (two parallel workflows: no-transitions vs.
   mapDamage-rescaled)
3. `3_diversity_stats.smk` — per-sample heterozygosity + pixy pi/dxy/Fst
4. `4_relatedness.smk` — one ANGSD beagle file for all samples
   (`gl_model: GATK`, dedup or rescaled BAM stage), then ngsRelate
   (IBSrelate/SFS) for pairwise R0/R1/KING
5. `5_LD_pruning.smk` — ngsLD (container) → prune_graph (local binary, path
   in `params.run_ld_estimation.prune_graph`) → `prune_beagle.py` →
   LD-pruned, unlinked-SNP beagle file
6. `6a_Strucure_angsdPCA.smk` — drops `params.run_pca.exclude_samples`, then
   every remaining individual of `params.run_pca.exclude_populations`
   (membership from `params.run_pca.poplist`, default `Pops.txt`), from
   the pruned beagle, then PCAngsd
7. `6b_Structure_NGadmix.smk` — NGSadmix for K = 1..`params.run_admixture.max_k`
   on `pca_exclude_beagle`'s output (same samples as the PCA), replicated to
   convergence by `ngsadmix.sh` (from PopGLen; NGSadmix is in
   `workflow/envs/angsd.yaml`), then `plot_admix.R` (plot + convergence summary
   with best-fit K = highest converged K), evalAdmix 0.961 (container) and
   `plot_evaladmix.R` per K. Tested end-to-end on a real slice (every 50th
   site of the relatedness beagle, K = 1-3) via a scratch Snakefile that
   `include:`s the rule file. The real pruned beagle now exists (see below)
   and `run_admixture` is TRUE in `config.yaml`.
8. `7_FroH_bcftools.smk` — `bcftools roh` (`-G30 --ignore-homref --AF-dflt 0.4
   -M <rec_rate> -O r`, settings from zjnolen/polyommatini-temporal-genomics
   `bcftools_roh.smk`) on the biallelic BCFs, for every `params.run_roh`
   dataset (notrans/rescaled) × call type (indCall/jointCall) × every
   `run_genotyping.maxMissing` value. Output is regions only (`RG` lines) in
   `results/roh/{dataset}/`. Hand-tested on the four real notrans jointCall
   BCFs (seconds each; 60 samples); not yet run through Snakemake because the
   repo was locked by the LD run. `rec_rate: "1e-8"` is bcftools' default, not
   a *L. sinapis* estimate. bcftools' "lines total/processed" log line looks
   low (about 11%) because with `--ignore-homref` it counts one sample's
   non-hom-ref sites; that is expected, not a filter problem.

9. `10a_Wolbachia_screen.smk` (`run_wolbachia_screen`, FALSE in
   `config.yaml` until the full run is started) — the nuclear reference +
   mtDNA + Wolbachia FASTAs (top-level `mtdna_fasta` / `wolbachia_fasta`) concatenated into
   `results/bwa_mapping/reference/{REF_NAME}_mtWolb.fa` (+ `.stb` and one
   BED per genome), then modern samples (paired-end) mapped with
   bwa mem → merge lanes per index → Picard MarkDuplicates → merge indices →
   clipOverlap (`results/bwa_mapping/{source}/`), split by genome at
   MAPQ ≥ `mapq` (`results/wolbachia_screen/bams/`), mosdepth per genome,
   and a long + wide table of reads / mean depth / breadth / depth ratios
   (`results/wolbachia_screen/wolbachia_screen.*.coverage_{long,wide}.tsv`).
   No infection call. Tested end-to-end on a slice (see session 2026-10-19).
10. `10b_InStrain_Wolbachia_screen.smk` (`run_instrain_wolbachia`, FALSE) —
   the user's spec: per sample, all paired-end runs (every lane and index)
   in one bowtie2 call to nuclear + mtDNA + Wolbachia (its own copy,
   `results/bowtie2_mapping/reference/{REF_NAME}_mtWolb.fa`),
   duplicates removed (Picard), then `inStrain profile <bam> <wolbachia.fa>
   -g <prodigal genes.fna> --min_mapq 2` against the Wolbachia FASTA only;
   all samples' `genome_info` tables in
   `results/instrain_wolbachia/instrain.*.genome_info.tsv`. Mapping FASTA,
   profile FASTA, `.stb` and genes file can be swapped in
   `params.run_instrain_wolbachia` (e.g. other microbial genomes). Tested
   end-to-end on the same slice (2 modern samples).
   - 10a and 10b are independent (user's requirement, 2026-10-19): neither
     uses a rule, constant or function of the other, and each runs or can be
     `include:`d alone. Shared inputs live outside both: `mtdna_fasta` /
     `wolbachia_fasta` (top-level config, Reference Files) and `samples_df`.
     Each builds its own concatenated FASTA with the same script
     (`concat_screen_reference.py`), so with both on it is built twice.
     Checked with scratch Snakefiles including 10a only, 10b only (also run
     end-to-end) and both.

**Session 2026-10-19 (stages 10a, 10b), things not obvious from the code:**
- The user wants three separate workflows: the bwa/mosdepth screen (10a),
  the inStrain screen (10b), both modern only for now, and later
  the consensus/tree workflow, which must work for modern AND historical
  samples. The user named the files 10a/10b (first written as 8a); they
  wrote "scrreen" for 10b, Claude used "screen".
- The bwa rules (`bwa_*`) take any reference by name: `{cref}` is looked up
  in `BWA_REFS`, and the FASTA must sit at `BWA_REF_DIR/{cref}.fa`. Only the
  screen's reference uses them now (10b maps with bowtie2).
- Header checks are in `concat_screen_reference.py`: every mtDNA header must
  contain "mitochondri", every Wolbachia header "wolbachia" (case ignored),
  no nuclear header either word, contig names unique. The mtDNA contig name
  is `ENA|FR990201|FR990201.1` (pipes are fine for bwa/samtools/mosdepth).
- Library design (user, 2026-10-19): the lanes of one index are the same
  PCR library, so duplicates are removed after merging lanes; different
  indices are different PCR libraries, so there is no dedup after merging
  indices.
- 10a read group `LB` is `{sample}_{index}`, matching that design. 1.1 sets
  `LB` to the lane at mapping (`get_rg_string`), and Picard groups by `LB`,
  not by file, so `deduplicate_modern` removed duplicates within each lane
  only (the 1.1 Picard metrics list `L1`/`L2` as separate libraries for the
  29 two-lane modern samples; the other 18 have one lane).
- Fixed in 1.1 (2026-10-19, user's request) by `modern_dedup_libraries`,
  between `merge_sample_indices` and `modern_clip_overlap`, modern only:
  `samtools reheader` sets every read group's `LB` to `{sample}_{index}`
  (RG IDs rebuilt from the sample sheet with `get_rg_string`; all 76 existing
  RG IDs matched), then Picard MarkDuplicates again. Its output is temp();
  the clipped BAM keeps its name, so no downstream path changed (the clipped
  name does not show the extra step). `get_rg_string` and the mapping rules
  were not changed, so the mapping params don't change.
  - Because the new intermediate is temp(), a plain `snakemake` sees the
    existing clipped BAMs as up to date ("Nothing to be done"). The rerun is
    opt-in: `snakemake -R modern_dedup_libraries`, which dry-runs to 494 jobs:
    47 × (dedup, clip, clipped stats/depth/qualimap, individual calls), M1-M4
    subsampling, joint calls/merges, and everything downstream including
    ngsLD (43 h, 716 GB). No historical per-sample job is scheduled.
  - Real effect is tiny: on chromosome 47 (pairs with both mates on it),
    LsKaz10 (HiSeq, 2 lanes, 1.7% within-lane duplication) lost 71 of
    439,294 pairs (0.016%), M11 (NovaSeq, 2 lanes, 27-29% within-lane) 277
    of 445,390 (0.06%). Control: MarkDuplicates on the same slice without
    the relabel removed 0. The earlier 207,446 vs 148,590 comparison came
    from a slice of the first reads of each lane and was not representative.
    Whether to rerun is the user's call (not run as of 2026-10-19).
- 10b maps all runs of a sample in one bowtie2 call (user's spec), so it has
  one read group per sample and dedups across indices too, which goes
  against the design above. No effect now (every modern sample has a single
  index, `01`); for multi-index samples (e.g. historical) it would need one
  bowtie2 call per index.
- 34 rows of `modern_samples.tsv` have a space instead of a tab between
  `rgGroup_ID` and `platform` (e.g. `C80HPANXX_L1_01 hiseq2000`), so their
  RG ID contains the platform and PL is missing (`unknown` in the vg BAMs,
  `ILLUMINA` in 10a). Harmless (IDs stay unique); the sheet was not edited.
- mosdepth's summary `total` row only covers contigs that had reads, and
  `regions.bed.gz` rounds the mean to 2 decimals, so the summary script
  computes depth as Σ per-contig `bases` / genome length (from the
  thresholds file, which lists every contig) instead.
- Real-repo dry runs: 10a alone 550 jobs, 10b alone 146 (no 10a job),
  both 695 (47 modern samples, 76 runs). The 10a clipped BAM and the
  per-genome BAMs are kept (the tree workflow will read the mtDNA/Wolbachia
  ones); the 10b deduplicated BAM is kept (inStrain input).
- Envs: `bwa.yaml` (bwa 0.7.19, samtools 1.21), `mosdepth.yaml` (0.3.14),
  `bowtie2.yaml` (2.5.5, samtools 1.21), `prodigal.yaml` (2.6.3),
  `instrain.yaml` (inStrain 1.10.0). inStrain needs the `defaults`
  channel: with conda-forge + bioconda only, the solve ran about 40 min
  without finishing; with `-c conda-forge -c bioconda -c defaults` (the
  user's suggestion) it solved in 2.5 min (python 3.8). The env has no
  samtools; inStrain reads the BAMs through pysam, and `is_dedup` indexes
  them first.
- inStrain `--min_mapq` is strict (MAPQ > value) and applies to either read
  of a pair. Without `-s` every contig is one genome, `all_scaffolds`.
- inStrain writes the 30 `genome_info` columns in a different order from
  run to run (seen between samples of one run), so
  `instrain_genome_info_summary.py` matches columns by name; matching by
  position failed on the second test run.
- No blanks (negative controls) exist for the modern samples; blanks only
  apply to historical samples (user, 2026-10-19). Blank / single-end support
  was written for 10a/10b and then removed at the user's request; both
  screens map paired-end modern reads only and stop with an error on a run
  without fq2. Blanks belong with the historical samples, if those are added.
- Hand test, LsKaz10 slice (3M pairs, slice reference): inStrain on the
  Wolbachia FASTA gave coverage 2.60x, breadth 0.821, breadth_expected
  0.899; 10a's mosdepth on the same reads (bwa) gave 2.69x, 0.824 at 1x.
  bowtie2 end-to-end aligned 26% of the slice reads, bwa mem 53-58% (bwa
  soft-clips reads from chromosomes missing from the slice into it).

**Session 2026-10-05 (stage 6b), things not obvious from the code:**
- `ngsadmix.sh` differs from upstream PopGLen by one line: its last line
  redirected to `"$snakemake_log[0]"`, which writes the log to a file
  literally named `<log>[0]`; now `"${snakemake_log[0]}"`. It also hardcodes
  "top 3" for convergence, so `conv` is fixed at 3 in the rule file
  (`ADMIX_CONV`), not a config option.
- `plot_admix.R` differs from upstream: `(6/45)*length()` crash for >50
  samples fixed, works for a single K, adds `bestLike`/`bestK` columns to the
  convergence summary. `plot_evaladmix.R` is upstream plus row-count checks
  and `dev.off()`; it needs `/usr/local/bin/visFuns.R` from the evalAdmix
  container, so its rule uses `container:` too.
- Admixture plots take populations from `Pops.txt` (all 61 samples), not
  the `Pops_GTS.txt` pixy uses (lacks `NHM091`).
- NGSadmix applies its own `-minMaf 0.05` after samples are excluded
  (45,894 → 45,188 sites on the test slice); evalAdmix applies the same
  filter, so the `.fopt.gz` and beagle stay consistent without extra flags.
- `workflow/envs/r_admix.yaml` and the evalAdmix image are already built
  under `.snakemake/` (the test used the repo's conda/apptainer prefixes).
- Testing pattern while the repo is locked by a running job: a scratch
  directory with a small Snakefile that stubs the upstream constants
  (`REF_NAME`, `PRUNE_TAG`, `PCA_TAG`, `KEEP_SAMPLES_REL`,
  `PCA_KEEP_SAMPLES`, `samples_df`) and `include:`s the one rule file,
  run with `--conda-prefix`/`--apptainer-prefix` pointing at the repo's
  `.snakemake/`.
- Test timing, 45,894 sites × 59 individuals, 4 threads, busy machine: 20
  replicates took under 1 min at K=1, about 7 min at K=2, about 10 min at
  K=3; all three converged at 20 replicates.

**Session 2026-10-12 (pixy window plots, no new pipeline stage):**
- `utilities/plot_pixy_modern_windows.R <pixy_dir> [out_prefix] [modern_pops]
  [min_sites] [min_snps]` reads `pixy_fst.txt`, `pixy_dxy.txt` and
  `pixy_pi.txt` from one pixy output folder and writes
  `<out_prefix>_{fst,dxy,pi}.pdf` (default prefix
  `<pixy_dir>/pixy_modern_windows`). Each PDF is one page: the statistic
  against genomic windows (chromosomes end to end, every other one on a grey
  band, line at y = 0), with one stacked panel per modern pair (Fst, dxy) or
  modern population (pi). Arguments are positional.
- "Modern" is a hardcoded default, `Spa,Kaz,Swe,Nth,Dev` (everything except
  `NthH`/`DevH`), which is also the panel order; the y-axis range is computed
  from the populations passed, so plots of different population sets are not
  on the same scale. Each statistic has its own y-range. pi's axis is cut at
  the 99.9th percentile of the plotted pi windows (0.0112 on `fmiss0.2`), and
  the windows above it (53 of 52,550) are drawn as triangles at the top edge;
  the user asked for the pi panels to use more of their height and accepted
  pi and dxy not sharing an axis. dxy is not cut.
- Default filters: dxy/pi windows need `no_sites >= 100` (`min_sites`), Fst
  windows need `no_snps >= 10` (`min_snps`; pixy's Fst table has no
  `no_sites`). Without them a few near-empty windows dominate: pi/dxy reach
  about 0.11 against a median of 0.004-0.005, and Fst from under 10 SNPs is
  mostly noise. On `fmiss0.2` modern windows this drops 460 of 105,560 (dxy),
  230 of 52,780 (pi) and 848 of 105,312 (Fst). The cutoffs are judgement
  calls, not derived from anything.
- The plot has no jitter. Negative Fst values are real (Hudson estimator):
  about 3% of plotted modern-pair windows, 24% for Kaz vs Swe.
- The user wanted per-window genome scans from the windowed tables, not a
  summary of the `pixy_*_genomewide.tsv` files (those come from
  `workflow/scripts/pixy_genomewide_bootstrap.R`, run by hand, and exist only
  for `fmiss0.1` and `fmiss0.2`).
- The script needs only `ggplot2` and is written to also run on the system R
  (`/usr/bin/Rscript`, R 4.5.2 with ggplot2 3.3.6): it uses `size` rather
  than `linewidth` and base `pdf()` rather than `ggsave()`, whose PNG device
  fails there with "Graphics API version mismatch". The pixy data has 45
  chromosomes: 1 and 4-47 (no 2 or 3).
- Only tested with output to a scratch directory; no plots have been written
  under `results/` by Claude.

**Real data facts** (2026-09-21 relatedness run, `bam_stage: dedup`): 61
samples, 2,294,746 genome-wide SNP sites. Chromosome names are plain
integers.

**Full-genome LD run:** `ngsld_estimate` with `max_kb_dist: 4000` took 43 h
and wrote a 716 GB `.ld.gz`. `ngsld_prune` with the containerised
`prune_graph` (50 kb / r2 0.1) ran 5+ days without finishing. Resolved (user,
reported 2026-10-12) by switching to a local build of a recent `prune_graph`
(`/mnt/sda/saad/tools/prune_graph/target/release/prune_graph`, set in
`config.yaml`); the pruned beagle in `results/ld/` is dated 2026-10-08 and
`results/structure/{pca,ngsadmix}` have output.

**Next planned stages** (plan as of 2026-10-19):

1. **Run 10a and 10b on the real data** (both toggles FALSE until then) and
   look at the coverage tables. Nothing for an infection call is built yet.
2. **Wolbachia/mtDNA consensus + trees**, its own workflow, for modern AND
   historical samples (historical mapping for the screen is not built yet).
   Earlier proposal (2026-10-12, not confirmed): infected / uninfected sets
   from the screen table (thresholds in config, a checkpoint); haploid calls
   (bcftools, ploidy 1) and a depth-masked consensus per sample for
   Wolbachia and mtDNA from the per-genome BAMs of 10a; IQ-TREE (ModelFinder
   + ultrafast bootstrap); an R script comparing the two trees (tanglegram,
   Robinson-Foulds, optionally AU tests). Open: infection thresholds; mtDNA
   tree for infected samples only or all; outgroups; whether there is enough
   within-species variation; multiple Wolbachia strains per sample;
   packaging of iqtree / R ape / phytools not checked.
3. RoH summary + plots, added to `7_FroH_bcftools.smk`: per-individual FRoH
   and plots from the `bcftools_roh` region files. Not designed yet: needs
   the autosome length for the denominator, a minimum RoH length and length
   bins; the reference repo's `plot_Froh.R`
4. One the RoH run is complete calculate individuals heterozygosity again as in 3_diversity_stats.smk but this time inside and outside RoHs of a specified size
   (zjnolen/polyommatini-temporal-genomics) is the likely starting point.
5. Effective population size with GONE2, modern samples only, as its own
   new workflow file. Nothing designed or checked yet: input genotype set,
   per-population grouping, and whether GONE2 has a bioconda/conda-forge
   package (check before choosing `conda:` vs `container:`) are all open.

See `HANDOFF_2026-09-28.md` for the detailed session log of stages 5–6a
(what was tested, how, and why each design choice was made).
