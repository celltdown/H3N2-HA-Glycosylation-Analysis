import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv(
    "results/HA_motif_position_comparison.csv"
)

years = ["2023", "2024", "2025"]

matrix = df[years].to_numpy(dtype=float).T
positions = df["Alignment_Position"].astype(str).tolist()

fig, ax = plt.subplots(figsize=(13, 4.5))

image = ax.imshow(
    matrix,
    cmap="YlOrRd",
    vmin=0,
    vmax=10,
    aspect="auto"
)

ax.set_xticks(np.arange(len(positions)))
ax.set_xticklabels(positions)

ax.set_yticks(np.arange(3))
ax.set_yticklabels(years)

ax.set_xlabel("HA alignment position")
ax.set_ylabel("Collection year")

ax.set_title(
    "Predicted N-linked Glycosylation Motif Frequency "
    "in Influenza A(H3N2) HA"
)

for i in range(matrix.shape[0]):
    for j in range(matrix.shape[1]):
        value = int(matrix[i, j])
        color = "white" if value >= 7 else "black"
        ax.text(
            j, i, str(value),
            ha="center",
            va="center",
            color=color,
            fontsize=10
        )

cbar = fig.colorbar(image, ax=ax)
cbar.set_label("Sequences with motif (n/10)")

plt.tight_layout()

plt.savefig(
    "plots/HA_glycosylation_heatmap.png",
    dpi=300,
    bbox_inches="tight"
)

plt.savefig(
    "plots/HA_glycosylation_heatmap.pdf",
    bbox_inches="tight"
)

print("Heatmap generated successfully!")
print("Total motif positions:", len(positions))
print("Files saved in plots/")
