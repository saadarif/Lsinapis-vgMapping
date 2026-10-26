# ==============================================================================
# STABILITY OF GONE2 / currentNe2 ESTIMATES TO SAMPLE SIZE AND COMPOSITION
#   One population's VCF (e.g. a population VCF from 9a_ContemporaryNE_est.smk,
#   results/contemporary_ne/<vcf>/population_vcfs/). GONE2 and currentNe2 are
#   run on the full set of n_max individuals, then on every subset with one
#   individual dropped, every subset with two dropped, and so on down to
#   subsets of n_min individuals (inclusive). Both programs get the same -r
#   and are run either all with -x or all without (params: run_ne_stability: x).
#
#   Each subset: its samples, then sites with a minor allele count >= 2 in the
#   subset, so sites that become singletons or monomorphic when individuals
#   are dropped are removed again (as in 9a). Missingness is not filtered
#   again.
#
#   Plots: all GONE2 Ne trajectories (full set in bold, smaller subsets in
#   lighter colours), and currentNe2 Ne against the number of individuals with
#   a linear regression; plus a table of all currentNe2 estimates.
#
# There are C(n_max, n) subsets of n individuals, e.g. for n_max = 11: 1 of 11,
# 11 of 10, 55 of 9, 165 of 8, 330 of 7, 462 of 6, so 1024 subsets in total
# with n_min = 6. Each is one GONE2 and one currentNe2 run.
#
# The VCF's sample names are read when the workflow is parsed (to list the
# subsets), so the VCF must exist before this workflow is switched on: run 9a
# first when using one of its population VCFs.
#
# Subsets are named n<size>.c<number>, numbered in itertools.combinations
# order over the sorted sample names; <vcf name>.n<min>-<max>.subsets.tsv
# lists the samples dropped and kept in each. The numbering depends only on
# the VCF's samples, so changing n_min reuses existing runs.
#
# Results go in results/contemporary_ne_stability/<vcf name>/ (logs/ and
# benchmarks/ likewise), the vcf name being its file name without .vcf.gz,
# .vcf or .bcf:
#   gone2/       GONE2 Ne, d2 and STATS files per subset
#   currentne2/  currentNe2 output per subset, and the summary table
#   plots/       GONE2 trajectories, currentNe2 Ne vs number of individuals
# ==============================================================================
import gzip
import itertools
import math

NES_PARAMS = config.get("params", {}).get("run_ne_stability", {})

NES_VCF = NES_PARAMS.get("vcf", "")
NES_N_MAX = NES_PARAMS.get("n_max")
NES_N_MIN = NES_PARAMS.get("n_min")
NES_X = bool(NES_PARAMS.get("x", False))
NES_REC_RATE = NES_PARAMS.get("rec_rate")  # cM/Mb, -r of both programs
NES_GONE2 = NES_PARAMS.get("gone2", "")
NES_CURRENTNE2 = NES_PARAMS.get("currentne2", "")
# Same seed for every subset: GONE2 is a genetic algorithm, and the same seed
# and thread count give the same result
NES_GONE2_SEED = NES_PARAMS.get("gone2_seed", 1)
NES_GONE2_EXTRA = NES_PARAMS.get("gone2_extra", "")
NES_CURRENTNE2_EXTRA = NES_PARAMS.get("currentne2_extra", "")
NES_THREADS = NES_PARAMS.get("threads", 8)

NES_MODE = "x" if NES_X else "nox"
NES_GONE2_TAG = f"gone2.r{NES_REC_RATE}.S{NES_GONE2_SEED}.{NES_MODE}"
NES_CURRENTNE2_TAG = f"currentne2.r{NES_REC_RATE}.{NES_MODE}"


def nes_error(message):
    print(f"\nERROR: params: run_ne_stability: {message}\n")
    sys.exit(1)


def nes_vcf_samples(path):
    """Sample names from the header of a .vcf, .vcf.gz or .bcf file."""
    if path.endswith(".bcf"):
        # BCF: BGZF-compressed "BCF\2\2", header length (uint32), header text
        with gzip.open(path, "rb") as fh:
            if fh.read(5) != b"BCF\x02\x02":
                nes_error(f"{path} is not a BCF file.")
            header = fh.read(int.from_bytes(fh.read(4), "little")).decode()
        lines = header.split("\n")
    else:
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt") as fh:
            lines = []
            for line in fh:
                lines.append(line)
                if line.startswith("#CHROM"):
                    break
    for line in lines:
        if line.startswith("#CHROM"):
            return line.rstrip("\n").split("\t")[9:]
    nes_error(f"no #CHROM line in the header of {path}.")


# Subset name -> its samples, n_max down to n_min. Only worked out when the
# workflow is on, so a VCF that doesn't exist yet doesn't stop other runs.
# The name comes from the file name alone, so the rules' paths are well formed
# (no empty folder name) even when the workflow is off or the VCF is not there.
NES_NAME = "none"
for _ext in (".vcf.gz", ".vcf", ".bcf"):
    if os.path.basename(NES_VCF).endswith(_ext) and len(os.path.basename(NES_VCF)) > len(_ext):
        NES_NAME = os.path.basename(NES_VCF)[: -len(_ext)]
        break
NES_SAMPLES = []
NES_SUBSETS = {}
if config.get("run_ne_stability", False):
    if NES_NAME == "none":
        nes_error(f"vcf '{NES_VCF}' must end in .vcf.gz, .vcf or .bcf.")
    if not os.path.exists(NES_VCF):
        nes_error(
            f"vcf '{NES_VCF}' does not exist. Its sample names are needed to list the "
            "subsets, so it has to exist first (e.g. run run_contemporary_ne first)."
        )
    NES_SAMPLES = sorted(nes_vcf_samples(NES_VCF))
    if len(set(NES_SAMPLES)) != len(NES_SAMPLES):
        nes_error(f"{NES_VCF} has duplicate sample names.")
    try:
        NES_N_MAX, NES_N_MIN = int(str(NES_N_MAX)), int(str(NES_N_MIN))
    except ValueError:
        nes_error(f"n_max and n_min must be whole numbers, not '{NES_N_MAX}' and '{NES_N_MIN}'.")
    if NES_N_MAX != len(NES_SAMPLES):
        nes_error(f"n_max is {NES_N_MAX}, but {NES_VCF} has {len(NES_SAMPLES)} samples: {', '.join(NES_SAMPLES)}.")
    # GONE2 runs with as few as 3 individuals
    if not 3 <= NES_N_MIN <= NES_N_MAX:
        nes_error(f"n_min must be from 3 to n_max ({NES_N_MAX}), not {NES_N_MIN}.")
    try:
        if float(NES_REC_RATE) <= 0:
            raise ValueError
    except (TypeError, ValueError):
        nes_error(f"rec_rate must be a positive number (cM/Mb), not '{NES_REC_RATE}'.")
    for _key, _exe in (("gone2", NES_GONE2), ("currentne2", NES_CURRENTNE2)):
        if not (os.path.isfile(_exe) and os.access(_exe, os.X_OK)):
            nes_error(f"{_key} '{_exe}' is not an executable file.")
    # Zero-padded to widths that depend only on n_max, so names don't change with n_min
    _w_n = len(str(NES_N_MAX))
    _w_c = len(str(math.comb(NES_N_MAX, NES_N_MAX // 2)))
    for _n in range(NES_N_MAX, NES_N_MIN - 1, -1):
        for _i, _keep in enumerate(itertools.combinations(NES_SAMPLES, _n), start=1):
            NES_SUBSETS[f"n{_n:0{_w_n}d}.c{_i:0{_w_c}d}"] = list(_keep)

NES_DIR = f"results/contemporary_ne_stability/{NES_NAME}"
NES_LOG_DIR = f"logs/contemporary_ne_stability/{NES_NAME}"
NES_BENCH_DIR = f"benchmarks/contemporary_ne_stability/{NES_NAME}"
NES_RANGE = f"n{NES_N_MIN}-{NES_N_MAX}"

NES_SUBSETS_TABLE = f"{NES_DIR}/{NES_NAME}.{NES_RANGE}.subsets.tsv"
NES_GONE2_NE = f"{NES_DIR}/gone2/{NES_NAME}.{{nes_subset}}.{NES_GONE2_TAG}.Ne.tsv"
NES_CURRENTNE2_TXT = f"{NES_DIR}/currentne2/{NES_NAME}.{{nes_subset}}.{NES_CURRENTNE2_TAG}.txt"
NES_GONE2_PLOT = f"{NES_DIR}/plots/{NES_NAME}.{NES_RANGE}.{NES_GONE2_TAG}.Ne_trajectories.pdf"
NES_CURRENTNE2_SUMMARY = f"{NES_DIR}/currentne2/{NES_NAME}.{NES_RANGE}.{NES_CURRENTNE2_TAG}.summary.tsv"
NES_CURRENTNE2_PLOT = f"{NES_DIR}/plots/{NES_NAME}.{NES_RANGE}.{NES_CURRENTNE2_TAG}.Ne_vs_n.pdf"

# Requested from the Snakefile target list (section 17)
NES_TARGETS = []
if config.get("run_ne_stability", False):
    NES_TARGETS = [NES_GONE2_PLOT, NES_CURRENTNE2_SUMMARY, NES_CURRENTNE2_PLOT]

NES_CONSTRAINTS = {"nes_subset": r"n[0-9]+\.c[0-9]+"}


rule nes_subsets_table:
    """
    One row per subset: name, number of individuals, samples dropped and kept.
    """
    input:
        vcf = NES_VCF,
    output:
        tsv = NES_SUBSETS_TABLE,
    run:
        with open(output.tsv, "w") as fh:
            fh.write("subset\tn\tdropped\tkept\n")
            for sid, keep in NES_SUBSETS.items():
                dropped = [s for s in NES_SAMPLES if s not in keep]
                fh.write(f"{sid}\t{len(keep)}\t{','.join(dropped) or 'none'}\t{','.join(keep)}\n")


rule nes_gone2:
    """
    GONE2 on one subset. The subset VCF (its samples, then minor allele count
    >= 2) is written uncompressed into the job's own temporary folder, as GONE2
    needs a .vcf and writes temporary files next to its output. As in 9a, a run
    that ends with only a STATS file (no Ne estimate, e.g. -x not converging)
    keeps STATS and writes header-only Ne and d2 tables; no STATS fails the job.
    """
    input:
        vcf = NES_VCF,
    output:
        ne = NES_GONE2_NE,
        d2 = f"{NES_DIR}/gone2/{NES_NAME}.{{nes_subset}}.{NES_GONE2_TAG}.d2.tsv",
        stats = f"{NES_DIR}/gone2/{NES_NAME}.{{nes_subset}}.{NES_GONE2_TAG}.STATS.txt",
    wildcard_constraints:
        **NES_CONSTRAINTS,
    params:
        keep = lambda wildcards: ",".join(NES_SUBSETS[wildcards.nes_subset]),
        exe = NES_GONE2,
        rec_rate = NES_REC_RATE,
        seed = NES_GONE2_SEED,
        extra = NES_GONE2_EXTRA,
        x = "-x" if NES_X else "",
        out = "_GONE2_mix_" if NES_X else "_GONE2_",
    log: f"{NES_LOG_DIR}/gone2/gone2_{NES_NAME}.{{nes_subset}}.{NES_GONE2_TAG}.log"
    benchmark: f"{NES_BENCH_DIR}/gone2/gone2_{NES_NAME}.{{nes_subset}}.{NES_GONE2_TAG}.benchmark"
    conda: "../envs/bcftools121.yaml"
    threads: NES_THREADS
    shell:
        """
        tmp=$(mktemp -d "$(dirname {output.ne})/tmp.gone2.XXXXXX")
        trap 'rm -rf "$tmp"' EXIT
        (bcftools view -s {params.keep} -Ou {input.vcf} | \
            bcftools view -i 'MAC>=2' -Ov -o "$tmp/subset.vcf") 2> {log}
        status=0
        {params.exe} {params.x} -r {params.rec_rate} -S {params.seed} -t {threads} {params.extra} \
            -o "$tmp/out" "$tmp/subset.vcf" >> {log} 2>&1 || status=$?
        if [ ! -f "$tmp/out{params.out}STATS" ]; then
            echo "ERROR: GONE2 exited with status $status and wrote no STATS file" >> {log}
            exit 1
        fi
        mv "$tmp/out{params.out}STATS" {output.stats}
        if [ -f "$tmp/out{params.out}Ne" ]; then
            mv "$tmp/out{params.out}Ne" {output.ne}
            mv "$tmp/out{params.out}d2" {output.d2}
        else
            echo "NO Ne ESTIMATE: GONE2 (exit status $status) wrote only the STATS file, which gives the reason. Ne and d2 tables written with a header only." >> {log}
            printf 'Generation\\tNe\\n' > {output.ne}
            printf 'c_bin\\tnumber_of_SNP_pairs\\tObserved_d2\\tPredicted_d2\\n' > {output.d2}
        fi
        """


rule nes_currentne2:
    """
    currentNe2 on one subset, the subset VCF made as for GONE2 in the job's
    own temporary folder (currentNe2 writes a progress file next to its input).
    """
    input:
        vcf = NES_VCF,
    output:
        txt = NES_CURRENTNE2_TXT,
    wildcard_constraints:
        **NES_CONSTRAINTS,
    params:
        keep = lambda wildcards: ",".join(NES_SUBSETS[wildcards.nes_subset]),
        exe = NES_CURRENTNE2,
        rec_rate = NES_REC_RATE,
        extra = NES_CURRENTNE2_EXTRA,
        x = "-x" if NES_X else "",
    log: f"{NES_LOG_DIR}/currentne2/currentne2_{NES_NAME}.{{nes_subset}}.{NES_CURRENTNE2_TAG}.log"
    benchmark: f"{NES_BENCH_DIR}/currentne2/currentne2_{NES_NAME}.{{nes_subset}}.{NES_CURRENTNE2_TAG}.benchmark"
    conda: "../envs/bcftools121.yaml"
    threads: NES_THREADS
    shell:
        """
        tmp=$(mktemp -d "$(dirname {output.txt})/tmp.currentne2.XXXXXX")
        trap 'rm -rf "$tmp"' EXIT
        (bcftools view -s {params.keep} -Ou {input.vcf} | \
            bcftools view -i 'MAC>=2' -Ov -o "$tmp/subset.vcf") 2> {log}
        {params.exe} {params.x} -r {params.rec_rate} -t {threads} {params.extra} \
            -o {output.txt} "$tmp/subset.vcf" >> {log} 2>&1
        """


rule nes_plot_gone2:
    """
    GONE2 Ne trajectories of every subset: the full set in bold, smaller
    subsets in lighter colours by number of individuals; second page one panel
    per subset size.
    """
    input:
        ne = expand(NES_GONE2_NE, nes_subset=list(NES_SUBSETS)),
        subsets = NES_SUBSETS_TABLE,
    output:
        pdf = NES_GONE2_PLOT,
    params:
        subsets = list(NES_SUBSETS),
        n_max = NES_N_MAX,
        title = f"GONE2 Ne, subsets of {NES_N_MIN}-{NES_N_MAX} individuals",
        subtitle = f"{NES_NAME}\nGONE2 {'-x ' if NES_X else ''}-r {NES_REC_RATE} -S {NES_GONE2_SEED}",
    log: f"{NES_LOG_DIR}/plots/plot_gone2_{NES_NAME}.{NES_RANGE}.{NES_GONE2_TAG}.log"
    conda: "../envs/r_ne_plot.yaml"
    script: "../scripts/plot_ne_stability_gone2.R"


rule nes_currentne2_summary:
    """
    Every currentNe2 estimate (point estimate and 50% / 90% CI, per estimate
    type in the output: whole genome, between chromosomes only) in one table,
    and Ne against the number of individuals with a linear regression.
    """
    input:
        txt = expand(NES_CURRENTNE2_TXT, nes_subset=list(NES_SUBSETS)),
        subsets = NES_SUBSETS_TABLE,
    output:
        tsv = NES_CURRENTNE2_SUMMARY,
        pdf = NES_CURRENTNE2_PLOT,
    params:
        subsets = list(NES_SUBSETS),
        n_max = NES_N_MAX,
        title = f"currentNe2 Ne, subsets of {NES_N_MIN}-{NES_N_MAX} individuals",
        subtitle = f"{NES_NAME}\ncurrentNe2 {'-x ' if NES_X else ''}-r {NES_REC_RATE}",
    log: f"{NES_LOG_DIR}/plots/currentne2_summary_{NES_NAME}.{NES_RANGE}.{NES_CURRENTNE2_TAG}.log"
    conda: "../envs/r_ne_plot.yaml"
    script: "../scripts/plot_ne_stability_currentne2.R"
