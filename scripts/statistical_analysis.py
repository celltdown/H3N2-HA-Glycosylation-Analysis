import pandas as pd
from scipy.stats import fisher_exact

df = pd.read_csv(
    "results/HA_motif_position_comparison.csv"
)

results = []

for _, row in df.iterrows():
    a = int(row["2023"])
    b = int(row["2025"])

    table = [
        [a, 10 - a],
        [b, 10 - b]
    ]

    odds_ratio, p = fisher_exact(
        table, alternative="two-sided"
    )

    results.append({
        "Alignment_Position": int(row["Alignment_Position"]),
        "2023_Positive": a,
        "2025_Positive": b,
        "Frequency_Difference": b/10 - a/10,
        "Odds_Ratio": odds_ratio,
        "P_Value": p
    })

stats = pd.DataFrame(results)

# Benjamini-Hochberg correction
m = len(stats)
order = stats["P_Value"].sort_values().index
sorted_p = stats.loc[order, "P_Value"].values

adjusted = [0.0] * m
previous = 1.0

for i in range(m - 1, -1, -1):
    value = sorted_p[i] * m / (i + 1)
    previous = min(previous, value)
    adjusted[i] = previous

stats.loc[order, "Adjusted_P"] = adjusted

stats = stats.sort_values("Adjusted_P")

stats.to_csv(
    "results/HA_statistical_analysis.csv",
    index=False
)

print("\nSTATISTICAL RESULTS")
print(stats.to_string(index=False))

print("\nPositions with adjusted p < 0.05:",
      (stats["Adjusted_P"] < 0.05).sum())

print("\nStatistical analysis completed!")
