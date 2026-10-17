#!/usr/bin/env bash
# Shared helpers sourced by every step.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"
CONFIG="${CONFIG:-config/config.sh}"
# shellcheck source=../../config/config.sh
source "$CONFIG"

log() { printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*" >&2; }

require() {
    for tool in "$@"; do
        command -v "$tool" >/dev/null 2>&1 || { log "ERROR: '$tool' not found — activate the right conda environment (see README)"; exit 1; }
    done
}

# shellcheck disable=SC2034  # used by the scripts that source this file
META="${DATA_DIR}/metadata.tsv"
