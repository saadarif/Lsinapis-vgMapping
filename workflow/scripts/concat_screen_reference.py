"""
Builds the competitive mapping reference of the Wolbachia screen
(workflow/rules/10a_Wolbachia_screen.smk): the nuclear reference, the mtDNA and
the Wolbachia assembly concatenated into one FASTA, plus the files that say
which genome each contig belongs to.

Called via Snakemake's `script:` directive, so `snakemake.input`,
`snakemake.output` and `snakemake.log` are injected automatically.

Outputs:
  fasta           the three FASTAs one after the other, headers unchanged
  stb             contig -> genome, two tab-separated columns, no header (the
                  scaffold-to-bin file inStrain's -s takes)
  bed_<genome>    one BED per genome covering every contig of it end to end,
                  4th column = genome (used to split the BAMs and by mosdepth)

The genome of each contig is checked against its FASTA header: every header in
the mtDNA FASTA has to contain "mitochondri" (mitochondrion / mitochondrial)
and every header in the Wolbachia FASTA "wolbachia" (case ignored), and no
nuclear header may contain either word, since a second copy of the mtDNA or
Wolbachia among the nuclear contigs would give their reads MAPQ 0. Contig
names (the header up to the first space, as bwa and samtools use them) must be
unique across the three files.
"""

import gzip
import sys

sys.stderr = open(snakemake.log[0], "w")

KEYWORDS = {"mitochondrion": "mitochondri", "wolbachia": "wolbachia"}

# (genome label, FASTA path), in the order they are written
SOURCES = [
    ("nuclear", snakemake.input.nuclear),
    ("mitochondrion", snakemake.input.mtdna),
    ("wolbachia", snakemake.input.wolbachia),
]


def check_header(header, genome, path):
    """Raises if the header does not name the genome its FASTA is supposed to hold."""
    text = header.lower()
    found = {g for g, word in KEYWORDS.items() if word in text}
    if genome == "nuclear":
        if found:
            raise ValueError(
                f"Nuclear reference {path} has a contig whose header mentions "
                f"{'/'.join(sorted(found))}: '{header}'. Remove it from the nuclear "
                "FASTA (the screen adds its own mtDNA and Wolbachia sequences)."
            )
    elif found != {genome}:
        raise ValueError(
            f"Header '{header}' in {path} should contain '{KEYWORDS[genome]}' "
            f"(case ignored) and not the keyword of the other genome, since the "
            f"file is given as the {genome} FASTA."
        )


def open_fasta(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


lengths = {}  # contig -> length, in file order
genome_of = {}  # contig -> genome

with open(snakemake.output.fasta, "w") as fout:
    for genome, path in SOURCES:
        name = None
        with open_fasta(path) as fin:
            for line in fin:
                line = line.rstrip("\r\n")
                if line.startswith(">"):
                    header = line[1:].strip()
                    check_header(header, genome, path)
                    name = header.split()[0] if header else ""
                    if not name:
                        raise ValueError(f"Empty FASTA header in {path}")
                    if name in genome_of:
                        raise ValueError(
                            f"Contig name '{name}' occurs more than once "
                            f"({genome_of[name]} and {genome}); names must be unique "
                            "across the nuclear, mtDNA and Wolbachia FASTAs."
                        )
                    genome_of[name] = genome
                    lengths[name] = 0
                elif line:
                    if name is None:
                        raise ValueError(f"{path} does not start with a FASTA header")
                    lengths[name] += len(line.strip())
                else:
                    continue
                # Written line by line with an explicit newline, so a file
                # without a trailing newline can't glue its last sequence line
                # to the next file's first header.
                fout.write(line + "\n")
        if name is None:
            raise ValueError(f"No sequences found in the {genome} FASTA {path}")

empty = [c for c, n in lengths.items() if n == 0]
if empty:
    raise ValueError(f"Contigs without sequence: {', '.join(empty)}")

with open(snakemake.output.stb, "w") as fh:
    for contig, genome in genome_of.items():
        fh.write(f"{contig}\t{genome}\n")

for genome, _ in SOURCES:
    with open(snakemake.output[f"bed_{genome}"], "w") as fh:
        for contig, n in lengths.items():
            if genome_of[contig] == genome:
                fh.write(f"{contig}\t0\t{n}\t{genome}\n")

for genome, path in SOURCES:
    contigs = [c for c in lengths if genome_of[c] == genome]
    total = sum(lengths[c] for c in contigs)
    print(f"{genome}: {len(contigs)} contigs, {total} bp, from {path}", file=sys.stderr)
