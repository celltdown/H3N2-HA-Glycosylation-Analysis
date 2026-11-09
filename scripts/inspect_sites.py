from Bio import SeqIO
from collections import Counter

alignment = {
    r.id: str(r.seq).upper()
    for r in SeqIO.parse(
        "alignment/H3N2_HA_aligned.fasta", "fasta"
    )
}

positions = [79, 138, 149, 160]

for pos in positions:
    print("\n" + "=" * 45)
    print("ALIGNMENT POSITION:", pos)

    for year in [2023, 2024, 2025]:
        records = list(SeqIO.parse(
            f"data/HA_{year}.fasta", "fasta"
        ))

        aa_counts = Counter()
        motif_counts = Counter()

        for record in records:
            seq = alignment[record.id]

            aa = seq[pos - 1]
            region = seq[pos - 1:pos + 2]

            aa_counts[aa] += 1
            motif_counts[region] += 1

        print("\nYear:", year)
        print("Amino acids:", dict(aa_counts))
        print("Three-residue patterns:", dict(motif_counts))

print("\nAnalysis completed!")
