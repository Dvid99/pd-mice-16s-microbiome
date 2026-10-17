#!/usr/bin/env bash
# Step 1 — Import the demultiplexed single-end reads (Phred33 manifest) and
# summarise per-sample read counts and quality.
source "$(dirname "$0")/utils/common.sh"
require qiime

out="${OUT_DIR}/01_import"
mkdir -p "$out"

# Manifest paths are relative to the manifest location -> import from data/
(cd "$DATA_DIR" && qiime tools import \
    --type "SampleData[SequencesWithQuality]" \
    --input-format SingleEndFastqManifestPhred33V2 \
    --input-path manifest.tsv \
    --output-path "../${out}/demux_seqs.qza")

qiime demux summarize --i-data "${out}/demux_seqs.qza" --o-visualization "${out}/demux_seqs.qzv"
qiime metadata tabulate --m-input-file "$META" --o-visualization "${out}/metadata.qzv"
