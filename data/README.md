# Data

No data is stored in this repository; `bash scripts/00_download_data.sh` downloads everything from the [QIIME 2 data resources](https://docs.qiime2.org/2024.10/tutorials/pd-mice/).

| Input | Description |
|---|---|
| `metadata.tsv` | 48 fecal samples: 12 mice (6 cages × 2), 4 time points (days 7, 14, 21, 49 after transplant); genotype (wild type / ASO susceptible) × donor (healthy / Parkinson's disease) |
| `manifest.tsv`, `demultiplexed_seqs/` | Demultiplexed single-end Illumina reads, 16S rRNA V4 (515F/806R), 150 bp |
| `gg-13-8-99-515-806-nb-classifier.qza` | Naive Bayes classifier, Greengenes 13_8, 99% OTUs, V4 region (scikit-learn 1.4.2) |
| `sepp-refs-gg-13-8.qza` | SEPP reference tree and alignment (Greengenes 13_8) for fragment insertion |

Original study: Sampson TR *et al.* (2016). Gut microbiota regulate motor deficits and neuroinflammation in a model of Parkinson's disease. *Cell* 167:1469–1480.
