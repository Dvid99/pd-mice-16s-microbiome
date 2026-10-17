#!/usr/bin/env bash
# Step 2 — Quality filtering, denoising into amplicon sequence variants (ASVs)
# and chimera removal with DADA2 (reads truncated at TRUNC_LEN).
source "$(dirname "$0")/utils/common.sh"
require qiime

out="${OUT_DIR}/02_dada2"
mkdir -p "$out"

qiime dada2 denoise-single \
    --i-demultiplexed-seqs "${OUT_DIR}/01_import/demux_seqs.qza" \
    --p-trunc-len "$TRUNC_LEN" --p-n-threads "$THREADS" \
    --o-table "${out}/table.qza" \
    --o-representative-sequences "${out}/rep_seqs.qza" \
    --o-denoising-stats "${out}/stats.qza"

qiime metadata tabulate --m-input-file "${out}/stats.qza" --o-visualization "${out}/stats.qzv"
qiime feature-table summarize --i-table "${out}/table.qza" --m-sample-metadata-file "$META" \
    --o-visualization "${out}/table.qzv"
qiime feature-table tabulate-seqs --i-data "${out}/rep_seqs.qza" --o-visualization "${out}/rep_seqs.qzv"

# Per-sample stats as TSV (input, filtered, denoised, non-chimeric reads)
qiime tools export --input-path "${out}/stats.qza" --output-path "${out}/stats_export"
head -n 5 "${out}/stats_export/stats.tsv" | column -t -s $'\t'
