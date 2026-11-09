from Bio import SeqIO
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import re

rows = []

for year in [2023, 2024, 2025]:
    file = f"data/HA_{year}.fasta"

    for record in SeqIO.parse(file, "fasta"):
        seq = str(record.seq).upper().replace("-", "")

        sites = [
            m.start() + 1
            for m in re.finditer(r"(?=N[^P][ST])", seq)
        ]

        rows.append({
            "Year": year,
            "Accession": record.id,
            "Length": len(seq),
            "Number_of_Motifs": len(sites),
            "Motif_Positions": ";".join(map(str, sites))
        })

df = pd.DataFrame(rows)

df.to_csv(
    "results/HA_glycosylation_sites.csv",
    index=False
)

summary = df.groupby("Year")["Number_of_Motifs"].agg(
    ["count", "mean", "std", "min", "max"]
)

summary.to_csv(
    "results/HA_glycosylation_summary.csv"
)

print("\nGLYCOSYLATION SUMMARY")
print(summary.round(2))

plt.figure(figsize=(7, 5))

years = summary.index.astype(str)
means = summary["mean"]

plt.bar(years, means)

plt.xlabel("Collection Year")
plt.ylabel("Mean Number of N-X-S/T Motifs")
plt.title("Potential N-linked Glycosylation Motifs in H3N2 HA")

plt.tight_layout()
plt.savefig(
    "plots/HA_glycosylation_comparison.png",
    dpi=300
)

print("\nTotal sequences:", len(df))
print("Analysis completed successfully!")
