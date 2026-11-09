import requests
import time
import re
import csv
from io import StringIO
from Bio import SeqIO

YEAR = "2025"
BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

query = (
    'H3N2[All Fields] AND hemagglutinin[Title] '
    'AND "Homo sapiens"[All Fields] '
    'AND 2025[All Fields] NOT partial[Title]'
)

session = requests.Session()

response = session.get(
    BASE + "esearch.fcgi",
    params={
        "db": "protein",
        "term": query,
        "retmax": 500,
        "retmode": "json"
    },
    timeout=60
)
response.raise_for_status()

ids = response.json()["esearchresult"]["idlist"]
print("Candidate IDs:", len(ids), flush=True)

selected = []
metadata = []
seen = set()

for i in range(0, len(ids), 20):
    batch = ids[i:i+20]

    response = session.get(
        BASE + "efetch.fcgi",
        params={
            "db": "protein",
            "id": ",".join(batch),
            "rettype": "gb",
            "retmode": "text"
        },
        timeout=90
    )
    response.raise_for_status()

    for record in SeqIO.parse(
        StringIO(response.text), "genbank"
    ):
        source = next(
            (f for f in record.features
             if f.type == "source"), None
        )

        if source is None:
            continue

        q = source.qualifiers

        host = q.get("host", [""])[0].lower()
        subtype = q.get("serotype", [""])[0].upper()
        date = q.get("collection_date", [""])[0]
        strain = q.get("strain", [""])[0]
        seq = str(record.seq).upper()

        if host not in ["human", "homo sapiens"]:
            continue
        if subtype != "H3N2":
            continue
        if not re.search(r"\b2025\b", date):
            continue
        if "partial" in record.description.lower():
            continue
        if not 550 <= len(seq) <= 575:
            continue
        if set(seq) - set("ACDEFGHIKLMNPQRSTVWY"):
            continue
        if seq in seen:
            continue

        seen.add(seq)
        selected.append(record)

        metadata.append([
            record.id,
            date,
            host,
            strain,
            len(seq)
        ])

        print(
            "SELECTED:",
            record.id,
            date,
            len(seq),
            flush=True
        )

        if len(selected) == 10:
            break

    if len(selected) == 10:
        break

    time.sleep(0.5)

SeqIO.write(
    selected,
    "data/HA_2025.fasta",
    "fasta"
)

with open("data/HA_2025_metadata.csv", "w") as f:
    writer = csv.writer(f)
    writer.writerow([
        "Accession",
        "Collection_Date",
        "Host",
        "Strain",
        "Length"
    ])
    writer.writerows(metadata)

print("\nFINAL RESULT:", len(selected), "sequences")
