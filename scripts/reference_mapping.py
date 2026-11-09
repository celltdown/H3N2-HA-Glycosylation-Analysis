from Bio import AlignIO
import csv

original = AlignIO.read(
    "alignment/H3N2_HA_aligned.fasta", "fasta"
)

alignment = AlignIO.read(
    "alignment/H3N2_HA_with_reference.fasta", "fasta"
)

reference = [
    r for r in alignment
    if "P03437" in r.id
]

if len(reference) != 1:
    raise ValueError("Reference not found uniquely")

ref = reference[0]
original_sequences = {
    r.id: str(r.seq) for r in original
}
aligned_sequences = {
    r.id: str(r.seq) for r in alignment
}

for accession, sequence in original_sequences.items():
    if aligned_sequences.get(accession) != sequence:
        raise ValueError(
            "Original alignment changed: " + accession
        )

mapping = []
reference_position = 0

for alignment_pos, aa in enumerate(ref.seq, start=1):
    if aa != "-":
        reference_position += 1
        mature_position = reference_position - 16
    else:
        mature_position = None

    mapping.append({
        "Alignment_Position": alignment_pos,
        "Reference_AA": str(aa),
        "Reference_Position": (
            reference_position if aa != "-" else ""
        ),
        "H3_Mature_Position": (
            mature_position
            if mature_position is not None
            and mature_position > 0
            else ""
        )
    })

with open(
    "results/HA_reference_position_mapping.csv",
    "w",
    newline=""
) as f:
    writer = csv.DictWriter(
        f,
        fieldnames=mapping[0].keys()
    )
    writer.writeheader()
    writer.writerows(mapping)

print("\nREFERENCE MAPPING RESULTS")

for pos in [79, 138, 149, 160]:
    row = mapping[pos - 1]
    print(
        "Alignment:", pos,
        "| Reference:", row["Reference_Position"],
        "| H3 mature:", row["H3_Mature_Position"],
        "| AA:", row["Reference_AA"]
    )

print("\nReference length:", reference_position)
print("Mapping completed successfully!")
