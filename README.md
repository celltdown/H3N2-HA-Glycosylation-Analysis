# Influenza A(H3N2) Hemagglutinin Analysis

## Project Overview

Comparative bioinformatics analysis of hemagglutinin (HA) sequence conservation and putative N-linked glycosylation motifs in selected human influenza A(H3N2) sequences collected between 2023 and 2025.

This exploratory study investigates amino acid variability and differences in potential glycosylation motifs relevant to influenza antigenic research.

## Dataset

- Source: NCBI Protein Database
- Virus: Human influenza A(H3N2)
- Protein: Hemagglutinin (HA)
- Collection period: 2023-2025
- Total sequences: 30
- Sequences per year: 10
- Sequence length: 566 amino acids

Protein sequences were screened using available host, subtype, collection date, sequence length, and completeness information.

Identical amino acid sequences and repeated strain names were checked.

## Bioinformatics Workflow

1. Sequence retrieval from NCBI
2. Sequence filtering and metadata verification
3. Multiple sequence alignment using MAFFT
4. Sequence conservation analysis using Biopython
5. Detection of putative N-X-S/T glycosylation motifs
6. Comparison of motif frequencies across collection years
7. Inspection of amino acid substitutions
8. Geographic metadata analysis
9. Reference sequence alignment
10. Exploratory statistical analysis using Fisher's exact test and Benjamini-Hochberg correction

## Main Results

### Sequence Conservation

- Total alignment positions: 566
- Fully conserved positions: 508
- Variable positions: 58
- Mean positional conservation: 98.86%

### Putative Glycosylation Motifs

Mean motif counts per sequence:

| Year | Sequences | Mean | SD |
|------|-----------|------|----|
| 2023 | 10 | 12.8 | 0.42 |
| 2024 | 10 | 11.7 | 0.82 |
| 2025 | 10 | 12.1 | 0.32 |

Sixteen alignment positions were identified as potential motif start sites across the dataset.

### Selected Motif Differences

| Alignment position | 2023 | 2024 | 2025 |
|--------------------|------|------|------|
| 138 | 9/10 | 2/10 | 0/10 |
| 149 | 10/10 | 10/10 | 0/10 |
| 160 | 0/10 | 0/10 | 10/10 |

Exploratory Fisher's exact tests comparing 2023 and 2025 identified three positions with adjusted p-values below 0.05:

- Position 138: adjusted p = 0.000635
- Position 149: adjusted p = 0.000087
- Position 160: adjusted p = 0.000087

## Repository Structure

- `data/` - Protein sequences and sample metadata
- `alignment/` - MAFFT alignment files
- `scripts/` - Python analysis scripts
- `results/` - Statistical and sequence analysis tables
- `plots/` - Research figures

## Software

- Ubuntu Linux
- Python 3
- Biopython
- MAFFT
- SeqKit
- Pandas
- Matplotlib
- SciPy
- Requests

## Interpretation

The selected H3N2 sequences show high overall HA conservation alongside differences in putative N-linked glycosylation motifs.

These findings describe sequence-level patterns rather than experimentally confirmed glycosylation or immune escape.

The dataset is exploratory and was not geographically balanced or randomly sampled. Statistical comparisons should not be interpreted as global temporal trends.

## Data Availability

Sequence accession identifiers and collection metadata are provided in the `data/` directory.

Original sequence records can be retrieved from NCBI.

## Project

Celltdown Bioinformatics
