# Pipeline configuration (paths relative to the repository root).

QIIME_RELEASE="2024.10"
TUTORIAL_URL="https://data.qiime2.org/${QIIME_RELEASE}/tutorials/pd-mice"
CLASSIFIER_URL="https://data.qiime2.org/classifiers/sklearn-1.4.2/greengenes/gg-13-8-99-515-806-nb-classifier.qza"
SEPP_REFS_URL="https://data.qiime2.org/classifiers/sepp-ref-dbs/sepp-refs-gg-13-8.qza"

DATA_DIR="data"
OUT_DIR="results"
THREADS=2

# DADA2 (single-end 150 bp V4 reads)
TRUNC_LEN=150

# Diversity
SAMPLING_DEPTH=2000          # rarefaction depth for core metrics
RAREFACTION_MAX_DEPTH=4250
