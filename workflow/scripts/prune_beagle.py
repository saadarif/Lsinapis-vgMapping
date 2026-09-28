"""
Subsets a beagle genotype likelihood file (gzipped) down to the unlinked SNPs
picked by prune_graph in the ngsld_prune rule (workflow/rules/5_LD_estimation.smk).

Called via Snakemake's `script:` directive, so `snakemake.input`, `snakemake.output`
and `snakemake.log` are injected automatically; see that rule for how the two
inputs are produced.

This script is adapted fron Zach Nolen's PopGLen pipeline,
Original: https://github.com/zjnolen/PopGLen/blob/master/workflow/scripts/prune_beagle.py
"""

import gzip
import sys

sys.stderr = open(snakemake.log[0], "w")


def load_unlinked_markers(pos_path):
    """
    Reads the unlinked-sites file (prune_graph's -out), one site per line as
    "chr:pos", and returns the matching ANGSD beagle marker IDs ("chr_pos").
    rsplit on the last ':' keeps this correct even if a chromosome name itself
    contains a colon.
    """
    markers = set()
    with open(pos_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            chrom, pos = line.rsplit(":", 1)
            markers.add(f"{chrom}_{pos}")
    return markers


def prune_beagle(beagle_in_path, beagle_out_path, markers):
    """
    Streams the beagle file rather than loading it whole, keeping the header and
    only the rows whose marker column (column 1) is in `markers`.
    Returns the number of sites kept.
    """
    kept = 0
    with gzip.open(beagle_in_path, "rt") as fin, gzip.open(
        beagle_out_path, "wt"
    ) as fout:
        header = fin.readline()
        if not header:
            raise ValueError(f"{beagle_in_path} is empty, no header line found")
        fout.write(header)
        for line in fin:
            marker = line.split("\t", 1)[0]
            if marker in markers:
                fout.write(line)
                kept += 1
    return kept


markers = load_unlinked_markers(snakemake.input.pos)
n_kept = prune_beagle(snakemake.input.beagle, snakemake.output.beagle, markers)

print(
    f"Requested {len(markers)} unlinked sites, kept {n_kept} in the pruned beagle",
    file=sys.stderr,
)
if n_kept != len(markers):
    raise ValueError(
        f"Pruned beagle has {n_kept} sites but the unlinked-sites file listed "
        f"{len(markers)}; some requested markers were not found in the beagle "
        "file (mismatched beagle/pos inputs?)."
    )
