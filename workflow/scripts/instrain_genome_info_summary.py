"""
Concatenates the inStrain profile genome_info tables of every sample
(workflow/rules/10b_InStrain_Wolbachia_screen.smk) into one table, with
sample_id and source as the first two columns. inStrain's own columns are
kept as they are (coverage, breadth, breadth_expected, nucl_diversity, ...).

Called via Snakemake's `script:` directive; snakemake.input and
snakemake.params.units ([(source, sample_id), ...]) are in the same order.
"""

import sys

sys.stderr = open(snakemake.log[0], "w")

header = None
n_rows = 0
with open(snakemake.output[0], "w") as fout:
    for (source, sample), path in zip(snakemake.params.units, snakemake.input):
        with open(path) as fh:
            this_header = fh.readline().rstrip("\n")
            if header is None:
                header = this_header
                fout.write(f"sample_id\tsource\t{header}\n")
            elif this_header != header:
                raise ValueError(f"{path} has different columns from the first genome_info table")
            for line in fh:
                line = line.rstrip("\n")
                if line:
                    fout.write(f"{sample}\t{source}\t{line}\n")
                    n_rows += 1

print(f"{n_rows} genome rows from {len(snakemake.input)} samples", file=sys.stderr)
