from Bio import SeqIO
from collections import defaultdict, Counter
import csv

all_sequences = defaultdict(list)
all_strains = defaultdict(list)

for year in [2023, 2024, 2025]:
    fasta = f"data/HA_{year}.fasta"
    metadata = f"data/HA_{year}_metadata.csv"

    records = list(SeqIO.parse(fasta, "fasta"))

    with open(metadata, newline="") as f:
        rows = list(csv.DictReader(f))

    print("\nYEAR:", year)
    print("Sequence count:", len(records))
    print("Metadata count:", len(rows))

    for r in records:
        all_sequences[str(r.seq).upper()].append(
            (year, r.id)
        )

    for row in rows:
        strain = row["Strain"].strip()
        date = row["Collection_Date"].strip()

        all_strains[strain].append(
            (year, row["Accession"])
        )

        print(row["Accession"], "|", date, "|", strain)

print("\n=== FINAL QUALITY CONTROL ===")

print("Total records:", sum(
    len(v) for v in all_sequences.values()
))

print("Unique protein sequences:", len(all_sequences))

print("Unique strain names:", len(all_strains))

print("\nIDENTICAL SEQUENCES ACROSS RECORDS:")

for seq, samples in all_sequences.items():
    if len(samples) > 1:
        print(samples)

print("\nREPEATED STRAIN NAMES:")

for strain, samples in all_strains.items():
    if len(samples) > 1:
        print(strain, samples)

print("\nQC completed.")
