# ==============================================================================
# RECENT AND CONTEMPORARY EFFECTIVE POPULATION SIZE (GONE2, currentNe2)
#   From a multi-sample VCF given in the config (modern and/or historical
#   samples), for every population and missingness value listed:
#
#   1. Population VCF: the population's samples, biallelic SNPs with at most
#      maxMissing missing genotypes in that population and a minor allele count
#      of at least 2 (removes singletons, and sites that are monomorphic in the
#      population), then bcftools stats
#   2. GONE2 (Ne over the last 150 generations), run with and without -x
#      (-x: sample from a metapopulation of equal-sized subpopulations), and a
#      plot of both Ne trajectories
#   3. currentNe2 (contemporary Ne), run with and without -x
#
# Both programs are local installations (paths in params: run_contemporary_ne)
# and get the same constant recombination rate (-r, cM/Mb), so genetic
# distances come from the physical positions.
#
# Results for one VCF go in results/contemporary_ne/<vcf name>/ (logs/ and
# benchmarks/ likewise), the vcf name being its file name without .vcf.gz,
# .vcf or .bcf; every file name starts with it too. Pointing the config at
# another VCF gives a separate folder.
#   population_vcfs/  per-population sample lists, VCFs and bcftools stats
#   gone2/            GONE2 Ne, d2 and STATS files per mode, plus the Ne plot
#   currentne2/       currentNe2 output per mode
#
# Limits of the installed binaries (compile-time): 2,000,000 SNPs, 2000
# individuals, 1000 chromosomes. A larger population VCF makes them stop with
# "Reached max number of loci" in the log. GONE2 also leaves out SNPs with
# 20% or more missing genotypes by itself, whatever maxMissing is.
# ==============================================================================
import re

CNE_PARAMS = config.get("params", {}).get("run_contemporary_ne", {})

CNE_VCFS = CNE_PARAMS.get("vcf", []) or []
if isinstance(CNE_VCFS, str):
    CNE_VCFS = [CNE_VCFS]
CNE_POPS = CNE_PARAMS.get("populations", []) or []
if isinstance(CNE_POPS, str):
    CNE_POPS = [CNE_POPS]
CNE_MISSING_VALS = CNE_PARAMS.get("maxMissing", [0.0])
if not isinstance(CNE_MISSING_VALS, list):
    CNE_MISSING_VALS = [CNE_MISSING_VALS]
# Two columns, sample_id and population, no header
CNE_POPLIST = CNE_PARAMS.get("poplist", POPLIST_NT)
CNE_REC_RATE = CNE_PARAMS.get("rec_rate")  # cM/Mb, -r of both programs
CNE_GONE2 = CNE_PARAMS.get("gone2", "")
CNE_CURRENTNE2 = CNE_PARAMS.get("currentne2", "")
# GONE2 is a genetic algorithm: the same seed and thread count give the same result
CNE_GONE2_SEED = CNE_PARAMS.get("gone2_seed", 1)
CNE_GONE2_EXTRA = CNE_PARAMS.get("gone2_extra", "")
CNE_CURRENTNE2_EXTRA = CNE_PARAMS.get("currentne2_extra", "")
CNE_THREADS = CNE_PARAMS.get("threads", 8)

CNE_DIR = "results/contemporary_ne"
CNE_LOG_DIR = "logs/contemporary_ne"
CNE_BENCH_DIR = "benchmarks/contemporary_ne"

# "x": run with -x, "nox": without
CNE_MODES = ["nox", "x"]
CNE_GONE2_TAG = f"gone2.r{CNE_REC_RATE}.S{CNE_GONE2_SEED}"
CNE_CURRENTNE2_TAG = f"currentne2.r{CNE_REC_RATE}"


def cne_vcf_name(path):
    base = os.path.basename(path)
    for ext in (".vcf.gz", ".vcf", ".bcf"):
        if base.endswith(ext):
            return base[: -len(ext)]
    print(f"\nERROR: params: run_contemporary_ne: vcf '{path}' must end in .vcf.gz, .vcf or .bcf.\n")
    sys.exit(1)


# VCF name -> path, and population -> its samples in the poplist. Only worked
# out when the workflow is on, so a missing poplist doesn't stop other runs.
CNE_VCF_PATHS = {}
CNE_POP_SAMPLES = {}
if config.get("run_contemporary_ne", False):
    for _path in CNE_VCFS:
        _name = cne_vcf_name(_path)
        if _name in CNE_VCF_PATHS:
            print(f"\nERROR: params: run_contemporary_ne: two VCFs are both named '{_name}'.\n")
            sys.exit(1)
        CNE_VCF_PATHS[_name] = _path
    if not CNE_VCF_PATHS:
        print("\nERROR: params: run_contemporary_ne: vcf is empty.\n")
        sys.exit(1)
    if not CNE_POPS:
        print("\nERROR: params: run_contemporary_ne: populations is empty.\n")
        sys.exit(1)
    try:
        if float(CNE_REC_RATE) <= 0:
            raise ValueError
    except (TypeError, ValueError):
        print(f"\nERROR: params: run_contemporary_ne: rec_rate must be a positive number (cM/Mb), not '{CNE_REC_RATE}'.\n")
        sys.exit(1)
    for _key, _exe in (("gone2", CNE_GONE2), ("currentne2", CNE_CURRENTNE2)):
        if not (os.path.isfile(_exe) and os.access(_exe, os.X_OK)):
            print(f"\nERROR: params: run_contemporary_ne: {_key} '{_exe}' is not an executable file.\n")
            sys.exit(1)
    if not os.path.exists(CNE_POPLIST):
        print(f"\nERROR: poplist '{CNE_POPLIST}' (params: run_contemporary_ne: poplist) does not exist.\n")
        sys.exit(1)
    with open(CNE_POPLIST) as _fh:
        for _line in _fh:
            _fields = _line.split()
            if len(_fields) >= 2 and _fields[1] in CNE_POPS:
                CNE_POP_SAMPLES.setdefault(_fields[1], []).append(_fields[0])
    _missing_pops = [p for p in CNE_POPS if p not in CNE_POP_SAMPLES]
    if _missing_pops:
        print(f"\nERROR: params: run_contemporary_ne: populations {', '.join(_missing_pops)} not found in {CNE_POPLIST}.\n")
        sys.exit(1)
    # Population names go into file names
    _bad_names = [p for p in CNE_POPS if not re.fullmatch(r"[A-Za-z0-9_-]+", p)]
    if _bad_names:
        print(
            f"\nERROR: population names may only contain letters, digits, '_' and '-' "
            f"(they go into file names): {', '.join(_bad_names)}.\n"
        )
        sys.exit(1)

CNE_CONSTRAINTS = {
    "cne_vcf": "|".join(re.escape(n) for n in CNE_VCF_PATHS) if CNE_VCF_PATHS else r"[^/]+",
    "cne_pop": "|".join(CNE_POP_SAMPLES) if CNE_POP_SAMPLES else r"[A-Za-z0-9_-]+",
    "cne_miss": r"[0-9]+(\.[0-9]+)?",
    "cne_mode": "|".join(CNE_MODES),
}

# Folder and file name stem of one population VCF; everything downstream adds to it
CNE_POP_STEM = "{cne_vcf}/population_vcfs/{cne_vcf}.{cne_pop}.fmiss{cne_miss}.mac2"
CNE_RUN_STEM = "{cne_vcf}.{cne_pop}.fmiss{cne_miss}.mac2"

# Requested from the Snakefile target list (section 16)
CNE_TARGETS = []
if config.get("run_contemporary_ne", False):
    for _name in CNE_VCF_PATHS:
        for _pop in CNE_POPS:
            for _m_val in CNE_MISSING_VALS:
                _wc = dict(cne_vcf=_name, cne_pop=_pop, cne_miss=_m_val)
                _stem = CNE_RUN_STEM.format(**_wc)
                CNE_TARGETS.append(f"{CNE_DIR}/{CNE_POP_STEM.format(**_wc)}.vcf.gz.stats")
                CNE_TARGETS.append(f"{CNE_DIR}/{_name}/gone2/{_stem}.{CNE_GONE2_TAG}.Ne.pdf")
                for _mode in CNE_MODES:
                    CNE_TARGETS.append(f"{CNE_DIR}/{_name}/currentne2/{_stem}.{CNE_CURRENTNE2_TAG}.{_mode}.txt")


rule cne_population_samples:
    """
    The population's samples (from the poplist) that are in the VCF. Samples
    of the population missing from the VCF are listed in the log; none at all
    is an error.
    """
    input:
        vcf = lambda wildcards: CNE_VCF_PATHS[wildcards.cne_vcf],
    output:
        samples = f"{CNE_DIR}/{{cne_vcf}}/population_vcfs/{{cne_vcf}}.{{cne_pop}}.samples.txt",
    wildcard_constraints:
        **CNE_CONSTRAINTS,
    params:
        samples = lambda wildcards: " ".join(CNE_POP_SAMPLES[wildcards.cne_pop]),
    log: f"{CNE_LOG_DIR}/{{cne_vcf}}/population_vcfs/population_samples_{{cne_vcf}}.{{cne_pop}}.log"
    conda: "../envs/bcftools121.yaml"
    shell:
        """
        bcftools query -l {input.vcf} 2> {log} | \
            awk -v want="{params.samples}" 'BEGIN {{ n = split(want, w, " "); for (i = 1; i <= n; i++) keep[w[i]] = 1 }} ($1 in keep)' > {output.samples}
        for s in {params.samples}; do
            grep -qxF "$s" {output.samples} || echo "{wildcards.cne_pop} sample not in the VCF: $s" >> {log}
        done
        n=$(wc -l < {output.samples})
        echo "{wildcards.cne_pop}: $n of $(echo {params.samples} | wc -w) poplist samples are in the VCF" >> {log}
        if [ "$n" -eq 0 ]; then
            echo "ERROR: no {wildcards.cne_pop} sample is in the VCF" >> {log}
            exit 1
        fi
        """


rule cne_population_vcf:
    """
    One population's samples, INFO tags (AC, AN, MAF, F_MISSING, ...)
    recomputed for them, then biallelic SNPs with F_MISSING <= {cne_miss}
    and a minor allele count >= 2 in the population.
    """
    input:
        vcf = lambda wildcards: CNE_VCF_PATHS[wildcards.cne_vcf],
        samples = rules.cne_population_samples.output.samples,
    output:
        vcf = f"{CNE_DIR}/{CNE_POP_STEM}.vcf.gz",
        tbi = f"{CNE_DIR}/{CNE_POP_STEM}.vcf.gz.tbi",
        stats = f"{CNE_DIR}/{CNE_POP_STEM}.vcf.gz.stats",
    wildcard_constraints:
        **CNE_CONSTRAINTS,
    log: f"{CNE_LOG_DIR}/{{cne_vcf}}/population_vcfs/population_vcf_{CNE_RUN_STEM}.log"
    benchmark: f"{CNE_BENCH_DIR}/{{cne_vcf}}/population_vcfs/population_vcf_{CNE_RUN_STEM}.benchmark"
    conda: "../envs/bcftools121.yaml"
    threads: 2
    shell:
        """
        (bcftools view -S {input.samples} -Ou {input.vcf} | \
            bcftools +fill-tags -Ou -- -t all | \
            bcftools view --threads {threads} -v snps -m 2 -M 2 \
                -i 'F_MISSING<={wildcards.cne_miss} && MAC>=2' -Oz -o {output.vcf}) 2> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        bcftools stats -s - {output.vcf} > {output.stats} 2>> {log}
        """


rule cne_gone2:
    """
    GONE2 on one population VCF, with -x ({cne_mode} = x) or without (nox),
    unphased diploid genotypes (GONE2 -g 0, its default). Both programs need an
    uncompressed .vcf and write temporary files next to their input or output,
    so each job decompresses the VCF into its own temporary folder. GONE2 names
    its outputs <-o>_GONE2_Ne, _d2 and _STATS (_GONE2_mix_ with -x); they are
    renamed to the same names for both modes.

    GONE2 can finish without an Ne estimate and write only the STATS file,
    which gives the reason: e.g. with -x when the metapopulation search does
    not converge (exit status 0; seen on 6 Dev samples, 3 chromosomes), or
    with too few recombination bins (exit status 1). That is a result, not an
    error: STATS is kept, the Ne and d2 tables are written with a header only,
    and the log says so. No STATS file at all (input errors) fails the job.
    """
    input:
        vcf = f"{CNE_DIR}/{CNE_POP_STEM}.vcf.gz",
    output:
        ne = f"{CNE_DIR}/{{cne_vcf}}/gone2/{CNE_RUN_STEM}.{CNE_GONE2_TAG}.{{cne_mode}}.Ne.tsv",
        d2 = f"{CNE_DIR}/{{cne_vcf}}/gone2/{CNE_RUN_STEM}.{CNE_GONE2_TAG}.{{cne_mode}}.d2.tsv",
        stats = f"{CNE_DIR}/{{cne_vcf}}/gone2/{CNE_RUN_STEM}.{CNE_GONE2_TAG}.{{cne_mode}}.STATS.txt",
    wildcard_constraints:
        **CNE_CONSTRAINTS,
    params:
        exe = CNE_GONE2,
        rec_rate = CNE_REC_RATE,
        seed = CNE_GONE2_SEED,
        extra = CNE_GONE2_EXTRA,
        x = lambda wildcards: "-x" if wildcards.cne_mode == "x" else "",
        out = lambda wildcards: "_GONE2_mix_" if wildcards.cne_mode == "x" else "_GONE2_",
    log: f"{CNE_LOG_DIR}/{{cne_vcf}}/gone2/gone2_{CNE_RUN_STEM}.{CNE_GONE2_TAG}.{{cne_mode}}.log"
    benchmark: f"{CNE_BENCH_DIR}/{{cne_vcf}}/gone2/gone2_{CNE_RUN_STEM}.{CNE_GONE2_TAG}.{{cne_mode}}.benchmark"
    threads: CNE_THREADS
    shell:
        """
        tmp=$(mktemp -d "$(dirname {output.ne})/tmp.gone2.XXXXXX")
        trap 'rm -rf "$tmp"' EXIT
        gzip -dc {input.vcf} > "$tmp/population.vcf"
        status=0
        {params.exe} {params.x} -r {params.rec_rate} -S {params.seed} -t {threads} {params.extra} \
            -o "$tmp/out" "$tmp/population.vcf" > {log} 2>&1 || status=$?
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


rule cne_plot_gone2:
    """
    GONE2 Ne backward in time for one population VCF, with and without -x.
    """
    input:
        nox = f"{CNE_DIR}/{{cne_vcf}}/gone2/{CNE_RUN_STEM}.{CNE_GONE2_TAG}.nox.Ne.tsv",
        x = f"{CNE_DIR}/{{cne_vcf}}/gone2/{CNE_RUN_STEM}.{CNE_GONE2_TAG}.x.Ne.tsv",
    output:
        pdf = f"{CNE_DIR}/{{cne_vcf}}/gone2/{CNE_RUN_STEM}.{CNE_GONE2_TAG}.Ne.pdf",
    wildcard_constraints:
        **CNE_CONSTRAINTS,
    params:
        title = lambda wildcards: f"{wildcards.cne_pop}, maxMissing {wildcards.cne_miss}",
        subtitle = lambda wildcards: f"{wildcards.cne_vcf}\nGONE2 -r {CNE_REC_RATE} -S {CNE_GONE2_SEED}",
    log: f"{CNE_LOG_DIR}/{{cne_vcf}}/gone2/plot_gone2_{CNE_RUN_STEM}.{CNE_GONE2_TAG}.log"
    conda: "../envs/r_ne_plot.yaml"
    script: "../scripts/plot_gone2_ne.R"


rule cne_currentne2:
    """
    currentNe2 on one population VCF, with -x ({cne_mode} = x) or without
    (nox). Without -x, the number of full siblings (k) is estimated from the
    data (currentNe2 default). Decompressed into a temporary folder as for
    GONE2: currentNe2 writes a progress file next to its input, which the -x
    and nox jobs on the same VCF would otherwise share.
    """
    input:
        vcf = f"{CNE_DIR}/{CNE_POP_STEM}.vcf.gz",
    output:
        txt = f"{CNE_DIR}/{{cne_vcf}}/currentne2/{CNE_RUN_STEM}.{CNE_CURRENTNE2_TAG}.{{cne_mode}}.txt",
    wildcard_constraints:
        **CNE_CONSTRAINTS,
    params:
        exe = CNE_CURRENTNE2,
        rec_rate = CNE_REC_RATE,
        extra = CNE_CURRENTNE2_EXTRA,
        x = lambda wildcards: "-x" if wildcards.cne_mode == "x" else "",
    log: f"{CNE_LOG_DIR}/{{cne_vcf}}/currentne2/currentne2_{CNE_RUN_STEM}.{CNE_CURRENTNE2_TAG}.{{cne_mode}}.log"
    benchmark: f"{CNE_BENCH_DIR}/{{cne_vcf}}/currentne2/currentne2_{CNE_RUN_STEM}.{CNE_CURRENTNE2_TAG}.{{cne_mode}}.benchmark"
    threads: CNE_THREADS
    shell:
        """
        tmp=$(mktemp -d "$(dirname {output.txt})/tmp.currentne2.XXXXXX")
        trap 'rm -rf "$tmp"' EXIT
        gzip -dc {input.vcf} > "$tmp/population.vcf"
        {params.exe} {params.x} -r {params.rec_rate} -t {threads} {params.extra} \
            -o {output.txt} "$tmp/population.vcf" > {log} 2>&1
        """
