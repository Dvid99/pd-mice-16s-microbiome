#!/usr/bin/env bash
# Step 4 — Phylogeny by fragment insertion (SEPP): ASVs are placed into the
# Greengenes 13_8 reference tree, which is more robust for short reads than a
# de novo tree.
source "$(dirname "$0")/utils/common.sh"
require qiime

out="${OUT_DIR}/04_phylogeny"
mkdir -p "$out"

qiime fragment-insertion sepp \
    --i-representative-sequences "${OUT_DIR}/02_dada2/rep_seqs.qza" \
    --i-reference-database "${DATA_DIR}/sepp-refs-gg-13-8.qza" \
    --p-threads "$THREADS" \
    --o-tree "${out}/tree.qza" --o-placements "${out}/placements.qza"
