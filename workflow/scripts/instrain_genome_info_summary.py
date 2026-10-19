"""
Concatenates the inStrain profile genome_info tables of every sample
(workflow/rules/10b_InStrain_Wolbachia_screen.smk) into one table, with
sample_id and source as the first two columns. inStrain's own columns are
kept as they are (coverage, breadth, breadth_expected, nucl_diversity, ...).

inStrain does not write the columns in the same order every run (same 30
columns, shuffled between samples), so rows are matched to the first table's
columns by name; a table with a different set of columns is an error.

Called via Snakemake's `script:` directive; snakemake.input and
snakemake.params.units ([(source, sample_id), ...]) are in the same order.
"""

import csv
import sys

sys.stderr = open(snakemake.log[0], "w")

columns = None
n_rows = 0
with open(snakemake.output[0], "w", newline="") as fout:
    writer = csv.writer(fout, delimiter="\t", lineterminator="\n")
    for (source, sample), path in zip(snakemake.params.units, snakemake.input):
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if columns is None:
                columns = reader.fieldnames
                writer.writerow(["sample_id", "source"] + columns)
            elif set(reader.fieldnames) != set(columns):
                raise ValueError(
                    f"{path} has different columns from the first genome_info table: "
                    f"extra {sorted(set(reader.fieldnames) - set(columns))}, "
                    f"missing {sorted(set(columns) - set(reader.fieldnames))}"
                )
            for row in reader:
                writer.writerow([sample, source] + [row[c] for c in columns])
                n_rows += 1

print(f"{n_rows} genome rows from {len(snakemake.input)} samples", file=sys.stderr)
