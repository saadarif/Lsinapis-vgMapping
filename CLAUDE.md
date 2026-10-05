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
  `5_LD_estimation`, `6a_Strucure_angsdPCA` (typo in filename is intentional/
  established, don't "fix" it without renaming deliberately),
  `6b_Structure_NGadmix`. All `include:`d
  into one global namespace — later files freely reuse Python-level constants
  and functions from earlier ones (e.g. `REF_NAME`, `samples_df`,
  `KEEP_SAMPLES_REL`, `REL_TAG`, `PRUNE_TAG`).
- `workflow/envs/*.yaml` — one conda env per rule (or small group), pinned
  versions. Even coreutils-only rules (`grep`/`awk`/`cut`) skip `conda:`
  entirely rather than reusing an unrelated env — see
  `individual_heterozygosity` in `3_diversity_stats.smk`.
- `workflow/scripts/` — Python/R/bash scripts invoked via `script:`.
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
  `PCA_TAG` (`excl-none` or `excl-<sorted-sample-ids>`) → `ADMIX_TAG`
  (= `PRUNE_TAG.PCA_TAG`, K goes in the file name as `.K{kvalue}`). A new stage
  consuming an existing output should depend on the upstream tag/prefix
  variable directly (e.g. `PRUNE_TAG`), not re-derive or hardcode it.
- **conda: vs container:** default to `conda:` with a pinned bioconda/
  conda-forge env. Only reach for `container:` when the tool genuinely has no
  bioconda/conda-forge package (checked via `bioconda-recipes` GitHub
  listing, not just `mamba search` — search past cases: `ngsLD`,
  `prune_graph` and `evalAdmix` have neither). When a `container:` rule is added,
  `apptainer` must already be enabled in `profiles/default/config.yaml`
  (it is, as of `5_LD_estimation.smk`) — no further config needed as long as
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
5. `5_LD_estimation.smk` — ngsLD (container) → prune_graph (container) →
   `prune_beagle.py` → LD-pruned, unlinked-SNP beagle file
6. `6a_Strucure_angsdPCA.smk` — drops `params.run_pca.exclude_samples` from
   the pruned beagle, then PCAngsd
7. `6b_Structure_NGadmix.smk` — NGSadmix for K = 1..`params.run_admixture.max_k`
   on `pca_exclude_beagle`'s output (same samples as the PCA), replicated to
   convergence by `ngsadmix.sh` (from PopGLen; NGSadmix is in
   `workflow/envs/angsd.yaml`), then `plot_admix.R` (plot + convergence summary
   with best-fit K = highest converged K), evalAdmix 0.961 (container) and
   `plot_evaladmix.R` per K. Tested end-to-end on a real slice (every 50th
   site of the relatedness beagle, K = 1-3) via a scratch Snakefile that
   `include:`s the rule file; not yet run on the real pruned beagle, which
   does not exist yet (see below). `run_admixture` is FALSE in `config.yaml`
   until `prune_graph` is sorted out.

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
  `PCA_EXCLUDE_SAMPLES`, `samples_df`) and `include:`s the one rule file,
  run with `--conda-prefix`/`--apptainer-prefix` pointing at the repo's
  `.snakemake/`.
- Test timing, 45,894 sites × 59 individuals, 4 threads, busy machine: 20
  replicates took under 1 min at K=1, about 7 min at K=2, about 10 min at
  K=3; all three converged at 20 replicates.

**Real data facts** (2026-09-21 relatedness run, `bam_stage: dedup`): 61
samples, 2,294,746 genome-wide SNP sites. Chromosome names are plain
integers.

**Full-genome LD run (as of 2026-10-05):** `ngsld_estimate` with
`max_kb_dist: 4000` took 43 h and wrote a 716 GB `.ld.gz`. `ngsld_prune`
(`prune_graph`, 50 kb / r2 0.1) started 2026-09-30 and was still running after
5+ days, with no progress output. Everything downstream (`prune_beagle`, PCA,
admixture) is waiting on it, and the repo's Snakemake lock is held while it
runs (`snakemake -n` still works). Open option, not decided: lower
`max_kb_dist` (pruning only uses pairs within 50 kb) and rerun.

**Next planned stage:** runs of homozygosity (RoH) with `bcftools roh`,
once the `prune_graph` runtime question is settled one way or the other. No
rule design decided yet (which genotyping workflow's calls to use, allele
frequency source, per-population vs. all samples are all open).
`workflow/envs/bcftools121.yaml` already exists; check `bcftools roh` in
that env before adding another.

See `HANDOFF_2026-09-28.md` for the detailed session log of stages 5–6a
(what was tested, how, and why each design choice was made).
