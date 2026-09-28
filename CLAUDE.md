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
  established, don't "fix" it without renaming deliberately). All `include:`d
  into one global namespace — later files freely reuse Python-level constants
  and functions from earlier ones (e.g. `REF_NAME`, `samples_df`,
  `KEEP_SAMPLES_REL`, `REL_TAG`, `PRUNE_TAG`).
- `workflow/envs/*.yaml` — one conda env per rule (or small group), pinned
  versions. Even coreutils-only rules (`grep`/`awk`/`cut`) skip `conda:`
  entirely rather than reusing an unrelated env — see
  `individual_heterozygosity` in `3_diversity_stats.smk`.
- `workflow/scripts/` — Python/R scripts invoked via `script:`.
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
  `PCA_TAG` (`excl-none` or `excl-<sorted-sample-ids>`). A new stage
  consuming an existing output should depend on the upstream tag/prefix
  variable directly (e.g. `PRUNE_TAG`), not re-derive or hardcode it.
- **conda: vs container:** default to `conda:` with a pinned bioconda/
  conda-forge env. Only reach for `container:` when the tool genuinely has no
  bioconda/conda-forge package (checked via `bioconda-recipes` GitHub
  listing, not just `mamba search` — search past cases: `ngsLD` and
  `prune_graph` have neither). When a `container:` rule is added,
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

**Real data facts** (2026-09-21 relatedness run, `bam_stage: dedup`): 61
samples, 2,294,746 genome-wide SNP sites. Chromosome names are plain
integers.

**Known open risk:** `ngsLD`/`prune_graph` at genome scale with a several-Mb
`max_kb_dist` may be slow — a 20kb/43k-site stress test put `prune_graph`
alone at 16 minutes for ~11M edges. Not yet run for real at full genome
scale as of the last update to this file.

**Next planned stage:** NGSadmix, on the exact same beagle file
`6a_Strucure_angsdPCA.smk`'s `pca_exclude_beagle` rule produces (i.e. the
input to `pca_pcangsd`, not its `.cov` output) — reuse `PRUNE_TAG`/`PCA_TAG`/
`PCA_EXCLUDE_SAMPLES` from that file rather than recomputing. NGSadmix ships
with ANGSD (`workflow/envs/angsd.yaml`) — verify it's actually on `PATH` in
that env before assuming it "just works" like PCAngsd did. Needs multiple
replicates per K for convergence checking; no rule design decided yet.

See `HANDOFF_2026-09-28.md` for the detailed session log of stages 5–6a
(what was tested, how, and why each design choice was made).
