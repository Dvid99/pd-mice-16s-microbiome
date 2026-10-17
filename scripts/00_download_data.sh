#!/usr/bin/env bash
# Step 0 — Download the public data of the QIIME 2 "Parkinson's Mouse" tutorial
# (48 fecal samples, 16S rRNA V4, single-end 150 bp) and the reference databases.
source "$(dirname "$0")/utils/common.sh"
require wget unzip

mkdir -p "$DATA_DIR"
get() { [[ -s "$2" ]] || { log "Downloading $(basename "$2")"; wget -q -O "$2" "$1"; }; }

get "${TUTORIAL_URL}/sample_metadata.tsv"     "$META"
get "${TUTORIAL_URL}/manifest"                "${DATA_DIR}/manifest.tsv"
get "${TUTORIAL_URL}/demultiplexed_seqs.zip"  "${DATA_DIR}/demultiplexed_seqs.zip"
get "$CLASSIFIER_URL"                         "${DATA_DIR}/gg-13-8-99-515-806-nb-classifier.qza"
get "$SEPP_REFS_URL"                          "${DATA_DIR}/sepp-refs-gg-13-8.qza"

# The manifest points to ./demultiplexed_seqs/<file>.fastq.gz (relative to data/)
[[ -d "${DATA_DIR}/demultiplexed_seqs" ]] || unzip -q -o "${DATA_DIR}/demultiplexed_seqs.zip" -d "$DATA_DIR"
log "Samples in manifest: $(( $(grep -vc '^#' "${DATA_DIR}/manifest.tsv") - 1 ))"
