#!/usr/bin/env bash
# Step 5 — Alpha rarefaction curves and core phylogenetic diversity metrics
# (Shannon, Faith PD, observed features, evenness; UniFrac and Bray-Curtis).
source "$(dirname "$0")/utils/common.sh"
require qiime

table="${OUT_DIR}/02_dada2/table.qza"
tree="${OUT_DIR}/04_phylogeny/tree.qza"
out="${OUT_DIR}/05_diversity"
mkdir -p "$out"

qiime diversity alpha-rarefaction \
    --i-table "$table" --i-phylogeny "$tree" --m-metadata-file "$META" \
    --p-min-depth 10 --p-max-depth "$RAREFACTION_MAX_DEPTH" \
    --o-visualization "${out}/alpha_rarefaction.qzv"
qiime tools export --input-path "${out}/alpha_rarefaction.qzv" --output-path "${out}/alpha_rarefaction_export"

rm -rf "${out}/core_metrics"
qiime diversity core-metrics-phylogenetic \
    --i-table "$table" --i-phylogeny "$tree" --m-metadata-file "$META" \
    --p-sampling-depth "$SAMPLING_DEPTH" --p-n-jobs-or-threads "$THREADS" \
    --output-dir "${out}/core_metrics"
