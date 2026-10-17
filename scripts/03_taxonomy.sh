#!/usr/bin/env bash
# Step 3 — Taxonomic classification of the ASVs with a naive Bayes classifier
# trained on Greengenes 13_8 (99% OTUs, 515F/806R V4 region), and taxa barplots.
source "$(dirname "$0")/utils/common.sh"
require qiime

out="${OUT_DIR}/03_taxonomy"
mkdir -p "$out"

qiime feature-classifier classify-sklearn \
    --i-reads "${OUT_DIR}/02_dada2/rep_seqs.qza" \
    --i-classifier "${DATA_DIR}/gg-13-8-99-515-806-nb-classifier.qza" \
    --p-n-jobs "$THREADS" \
    --o-classification "${out}/taxonomy.qza"
qiime metadata tabulate --m-input-file "${out}/taxonomy.qza" --o-visualization "${out}/taxonomy.qzv"

# Samples below the rarefaction depth are dropped for the barplot
qiime feature-table filter-samples --i-table "${OUT_DIR}/02_dada2/table.qza" \
    --p-min-frequency "$SAMPLING_DEPTH" --o-filtered-table "${out}/table_min_depth.qza"
qiime taxa barplot --i-table "${out}/table_min_depth.qza" --i-taxonomy "${out}/taxonomy.qza" \
    --m-metadata-file "$META" --o-visualization "${out}/taxa_barplot.qzv"
