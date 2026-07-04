#!/usr/bin/env bash
# Idempotent provisioning of the ChromBERT mirror for a given genome.
# Re-run safely after any mirror restart (installs/downloads are ephemeral unless
# the user snapshots the mirror). Does nothing if already provisioned.
#
#   ./scripts/setup_mirror.sh mm10     # or hg38
#
# Multi-species: genome is a parameter. ChromBERT natively supports BOTH hg38
# (6k regulators) and mm10 (5k) at 1 kb — so mouse needs NO liftOver; we just
# fetch the mm10 checkpoint + grid. liftOver stays a fallback only for models
# that genuinely lack a species (see docs/server_mirrors.md).
set -euo pipefail
GENOME="${1:?usage: setup_mirror.sh <hg38|mm10>}"
source "$(dirname "$0")/mirror_env.sh"

echo ">> [1/3] ensuring /root/.bashrc has PATH + ChromBERT env (for interactive sessions)"
# The repo scripts set these explicitly via mirror_py, so this is only for
# convenience in manual/interactive shells. Idempotent (guarded by a marker).
$SSH 'grep -q "ecr_navigator provisioning" ~/.bashrc 2>/dev/null || cat >> ~/.bashrc <<"RC"
# --- ecr_navigator provisioning ---
export PATH=/root/bin:$PATH
export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python
export HF_ENDPOINT=https://hf-mirror.com
# --- end ecr_navigator provisioning ---
RC'

echo ">> [2/3] ensuring bedtools (needed by chrombert_make_dataset)"
if ! $SSH 'source /opt/conda/etc/profile.d/conda.sh; conda activate base; which bedtools' >/dev/null 2>&1; then
  mirror_py "conda install -y -c bioconda -c conda-forge bedtools"
else
  echo "   bedtools present"
fi

echo ">> [3/3] ensuring ChromBERT $GENOME data (checkpoint + 1kb grid) in cache"
GRID="mm10_5k_1kb_region.bed"; [ "$GENOME" = "hg38" ] && GRID="hg38_6k_1kb_region.bed"
if ! $SSH "test -f ~/.cache/chrombert/data/config/$GRID"; then
  mirror_py "python -m chrombert.scripts.chrombert_prepare_env --genome $GENOME --resolution 1kb"
else
  echo "   $GENOME grid present"
fi
echo ">> setup complete for $GENOME"
