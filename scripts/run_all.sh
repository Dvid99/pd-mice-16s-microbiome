#!/usr/bin/env bash
# Run every step. Steps 0-6 and 8 need the QIIME 2 environment, step 7 the R one:
#   conda activate qiime2-amplicon-2024.10 && bash scripts/run_all.sh qiime
#   conda activate pd-mice-r               && bash scripts/run_all.sh r
source "$(dirname "$0")/utils/common.sh"

stage="${1:-all}"
if [[ "$stage" == "all" || "$stage" == "qiime" ]]; then
    for step in 00_download_data 01_import_demux 02_dada2 03_taxonomy 04_phylogeny 05_diversity 06_export; do
        log "===== ${step} ====="
        CONFIG="$CONFIG" bash "scripts/${step}.sh"
    done
    log "===== 08_plot_rarefaction ====="
    python scripts/08_plot_rarefaction.py "${OUT_DIR}/05_diversity/alpha_rarefaction_export/shannon.csv" \
        --group cage_id --out "${OUT_DIR}/07_phyloseq/figures/rarefaction_shannon_cage.png"
fi
if [[ "$stage" == "all" || "$stage" == "r" ]]; then
    log "===== 07_phyloseq_analysis ====="
    Rscript scripts/07_phyloseq_analysis.R
fi
