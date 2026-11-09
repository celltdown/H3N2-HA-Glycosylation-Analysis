from Bio import SeqIO
from collections import Counter, defaultdict
import pandas as pd
import re

years = [2023, 2024, 2025]
records = {}
sequence_groups = defaultdict(list)
site_counts = defaultdict(Counter)

for year in years:
    path = f"data/HA_{year}.fasta"
    records[year] = list(SeqIO.parse(path, "fasta"))

    for record in records[year]:
        sequence_groups[str(record.seq).upper()].append(
            (year, record.id)
        )

print("SEQUENCE QUALITY CONTROL")

for year in years:
    print(year, ":", len(records[year]), "sequences")

print("Unique sequences across all years:",
      len(sequence_groups))

duplicates = {
    seq: samples
    for seq, samples in sequence_groups.items()
    if len(samples) > 1
}

print("Repeated amino acid sequences:",
      len(duplicates))

aligned = {
    record.id: str(record.seq).upper()
    for record in SeqIO.parse(
        "alignment/H3N2_HA_aligned.fasta", "fasta"
    )
}

for year in years:
    for record in records[year]:
        if record.id not in aligned:
            raise ValueError(
                "Missing alignment: " + record.id
            )

        aligned_seq = aligned[record.id]
        raw_seq = str(record.seq).upper()

        if aligned_seq.replace("-", "") != raw_seq:
            raise ValueError(
                "Sequence mismatch: " + record.id
            )

        position_map = {}
        ungapped_pos = 0

        for alignment_pos, aa in enumerate(
            aligned_seq, start=1
        ):
            if aa != "-":
                ungapped_pos += 1
                position_map[ungapped_pos] = alignment_pos

        positions = [
            position_map[m.start() + 1]
            for m in re.finditer(
                r"(?=N[^P][ST])", raw_seq
            )
        ]

        for pos in positions:
            site_counts[pos][year] += 1

rows = []

for pos in sorted(site_counts):
    row = {"Alignment_Position": pos}

    for year in years:
        row[str(year)] = site_counts[pos][year]

    row["Difference_2025_vs_2023"] = (
        row["2025"] - row["2023"]
    )

    rows.append(row)

df = pd.DataFrame(rows)

df.to_csv(
    "results/HA_motif_position_comparison.csv",
    index=False
)

print("\nMOTIF POSITION COMPARISON")
print(df.to_string(index=False))

changed = df[
    (df["2023"] != df["2024"]) |
    (df["2024"] != df["2025"])
]

print("\nPositions with different frequencies:",
      len(changed))

print("\nAnalysis completed successfully!")
