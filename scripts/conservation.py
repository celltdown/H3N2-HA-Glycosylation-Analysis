from Bio import AlignIO
from collections import Counter
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

alignment = AlignIO.read(
    "alignment/H3N2_HA_aligned.fasta", "fasta"
)

results = []

for position in range(alignment.get_alignment_length()):
    column = alignment[:, position]
    amino_acids = [aa for aa in column if aa != "-"]

    counts = Counter(amino_acids)
    major_aa, major_count = counts.most_common(1)[0]

    conservation = major_count / len(amino_acids) * 100

    results.append({
        "Position": position + 1,
        "Major_AA": major_aa,
        "Conservation_Percent": conservation,
        "Variable": conservation < 100
    })

df = pd.DataFrame(results)
df.to_csv("results/HA_conservation.csv", index=False)

plt.figure(figsize=(14, 4))
plt.plot(df["Position"], df["Conservation_Percent"],
         linewidth=0.8)
plt.xlabel("Amino Acid Position")
plt.ylabel("Conservation (%)")
plt.title("H3N2 HA Sequence Conservation (2023–2025)")
plt.ylim(0, 105)
plt.tight_layout()
plt.savefig("plots/HA_conservation.png", dpi=300)

print("Total positions:", len(df))
print("Fully conserved positions:",
      (df["Conservation_Percent"] == 100).sum())
print("Variable positions:",
      (df["Conservation_Percent"] < 100).sum())
print("Mean conservation:",
      round(df["Conservation_Percent"].mean(), 2), "%")
print("Analysis completed successfully!")
