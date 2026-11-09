import csv
import re
from Bio import SeqIO
from collections import defaultdict

results = []
positions = [138, 149, 160]

for year in [2023, 2024, 2025]:
    fasta = {
        r.id: str(r.seq).upper()
        for r in SeqIO.parse(
            f"data/HA_{year}.fasta", "fasta"
        )
    }

    with open(f"data/HA_{year}_metadata.csv") as f:
        metadata = list(csv.DictReader(f))

    for row in metadata:
        accession = row["Accession"]
        strain = row["Strain"]
        sequence = fasta[accession]

        region = strain.split("/")[1]

        result = {
            "Year": year,
            "Region": region,
            "Accession": accession,
            "Strain": strain
        }

        for pos in positions:
            motif = sequence[pos-1:pos+2]

            valid = bool(
                re.fullmatch(r"N[^P][ST]", motif)
            )

            result[f"Position_{pos}"] = motif
            result[f"Motif_{pos}"] = int(valid)

        results.append(result)

output = "results/HA_geographic_motif_analysis.csv"

with open(output, "w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=results[0].keys()
    )
    writer.writeheader()
    writer.writerows(results)

print("\nGEOGRAPHIC MOTIF ANALYSIS")

groups = defaultdict(list)

for row in results:
    groups[(row["Year"], row["Region"])].append(row)

for (year, region), samples in sorted(groups.items()):
    print("\n", year, region, "n =", len(samples))

    for pos in positions:
        count = sum(
            r[f"Motif_{pos}"] for r in samples
        )
        print(f"Position {pos}: {count}/{len(samples)}")

print("\nSaved:", output)
