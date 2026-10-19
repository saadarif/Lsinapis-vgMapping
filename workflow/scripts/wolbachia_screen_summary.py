"""
Summarises the Wolbachia screen (workflow/rules/10a_Wolbachia_screen.smk):
per sample, the reads, mean depth and breadth of coverage on the nuclear,
mitochondrial and Wolbachia genomes, plus the depth ratios between them.

Called via Snakemake's `script:` directive. The input lists (mosdepth summary,
mosdepth thresholds, read count) and params.units are in the same order, one
entry per (source, sample_id, genome).

Depth and breadth are computed per genome from the per-contig mosdepth output
rather than taken from its "total" row, which only covers the contigs that
had reads and would overstate the depth of a genome with uncovered contigs:
  mean_depth  = sum of mosdepth "bases" over the genome's contigs / genome length
  breadth_Nx  = bases at depth >= N (mosdepth --thresholds) / genome length
The genome length is the sum of the contig lengths in the thresholds file,
which lists every contig of the genome's BED whether it had reads or not.

Outputs:
  long  one row per sample and genome
  wide  one row per sample, with the depth ratios
        mitochondrion/nuclear, wolbachia/nuclear and wolbachia/mitochondrion
        (NA when the denominator has no coverage)
"""

import gzip
import sys

sys.stderr = open(snakemake.log[0], "w")

GENOMES = ["nuclear", "mitochondrion", "wolbachia"]
RATIOS = [("mitochondrion", "nuclear"), ("wolbachia", "nuclear"), ("wolbachia", "mitochondrion")]

units = snakemake.params.units  # [(source, sample_id, genome), ...]
thresholds = [int(t) for t in snakemake.params.thresholds]
mapq = snakemake.params.mapq


def read_thresholds(path):
    """Returns (n contigs, genome length, bases at >= each threshold)."""
    n_contigs, length = 0, 0
    covered = [0] * len(thresholds)
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        expected = [f"{t}X" for t in thresholds]
        if header[4:] != expected:
            raise ValueError(f"{path}: threshold columns {header[4:]}, expected {expected}")
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            n_contigs += 1
            length += int(fields[2]) - int(fields[1])
            for i, value in enumerate(fields[4:]):
                covered[i] += int(value)
    return n_contigs, length, covered


def read_bases(path):
    """Sum of mosdepth 'bases' (total depth) over the per-contig rows of a summary."""
    total = 0
    with open(path) as fh:
        next(fh)  # header
        for line in fh:
            chrom, _, bases = line.split("\t")[:3]
            if chrom in ("total", "total_region") or chrom.endswith("_region"):
                continue
            total += int(bases)
    return total


def read_count(path):
    with open(path) as fh:
        return int(fh.read().strip())


def fmt(x, digits=6):
    return "NA" if x is None else f"{x:.{digits}g}"


rows = {}  # (source, sample_id) -> genome -> stats
for (source, sample, genome), summary, thr, nreads in zip(
    units, snakemake.input.summary, snakemake.input.thresholds, snakemake.input.nreads
):
    n_contigs, length, covered = read_thresholds(thr)
    if length == 0:
        raise ValueError(f"{thr}: genome {genome} has length 0")
    rows.setdefault((source, sample), {})[genome] = {
        "n_contigs": n_contigs,
        "length": length,
        "reads": read_count(nreads),
        "depth": read_bases(summary) / length,
        "breadth": [c / length for c in covered],
    }

breadth_cols = [f"breadth_{t}x" for t in thresholds]
reads_col = f"reads_mapq{mapq}"

with open(snakemake.output.long, "w") as fh:
    fh.write("\t".join(["sample_id", "source", "genome", "n_contigs", "length_bp", reads_col,
                        "mean_depth", *breadth_cols, "depth_ratio_to_nuclear"]) + "\n")
    for (source, sample), by_genome in sorted(rows.items()):
        nuc_depth = by_genome["nuclear"]["depth"]
        for genome in GENOMES:
            s = by_genome[genome]
            ratio = s["depth"] / nuc_depth if nuc_depth > 0 else None
            fh.write("\t".join([sample, source, genome, str(s["n_contigs"]), str(s["length"]),
                                str(s["reads"]), fmt(s["depth"]),
                                *[fmt(b) for b in s["breadth"]], fmt(ratio)]) + "\n")

with open(snakemake.output.wide, "w") as fh:
    cols = ["sample_id", "source"]
    for genome in GENOMES:
        cols += [f"{genome}_{reads_col}", f"{genome}_mean_depth"] + [f"{genome}_{c}" for c in breadth_cols]
    cols += [f"depth_ratio_{a}_to_{b}" for a, b in RATIOS]
    fh.write("\t".join(cols) + "\n")
    for (source, sample), by_genome in sorted(rows.items()):
        out = [sample, source]
        for genome in GENOMES:
            s = by_genome[genome]
            out += [str(s["reads"]), fmt(s["depth"])] + [fmt(b) for b in s["breadth"]]
        for a, b in RATIOS:
            den = by_genome[b]["depth"]
            out.append(fmt(by_genome[a]["depth"] / den if den > 0 else None))
        fh.write("\t".join(out) + "\n")

print(f"Summarised {len(rows)} samples x {len(GENOMES)} genomes", file=sys.stderr)
