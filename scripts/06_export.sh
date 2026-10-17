#!/usr/bin/env bash
# Step 6 — Export QIIME 2 artifacts to plain files for R/phyloseq:
# feature table (TSV), taxonomy (TSV), tree (Newick), UniFrac distances.
source "$(dirname "$0")/utils/common.sh"
require qiime biom

out="${OUT_DIR}/06_export"
mkdir -p "$out"
tmp="$(mktemp -d)"

qiime tools export --input-path "${OUT_DIR}/02_dada2/table.qza" --output-path "$tmp/table"
biom convert -i "$tmp/table/feature-table.biom" -o "${out}/feature_table.tsv" --to-tsv
qiime tools export --input-path "${OUT_DIR}/03_taxonomy/taxonomy.qza" --output-path "$tmp/tax"
cp "$tmp/tax/taxonomy.tsv" "${out}/taxonomy.tsv"
qiime tools export --input-path "${OUT_DIR}/04_phylogeny/tree.qza" --output-path "$tmp/tree"
cp "$tmp/tree/tree.nwk" "${out}/tree.nwk"
qiime tools export --input-path "${OUT_DIR}/05_diversity/core_metrics/unweighted_unifrac_distance_matrix.qza" \
    --output-path "$tmp/uu"
cp "$tmp/uu/distance-matrix.tsv" "${out}/unweighted_unifrac.tsv"
cp "$META" "${out}/metadata.tsv"
rm -rf "$tmp"
ls -l "$out"
