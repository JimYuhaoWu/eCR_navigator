#!/usr/bin/env bash
# Produce a ChromBERT region-embedding artifact for one cell state, driven from
# the repo (scripts live here; compute runs on the mirror; results come back).
#
#   ./scripts/run_chrombert_region_emb.sh \
#       --peaks data/MEF.e7_peaks.bed --genome mm10 --cell-state MEF \
#       --out artifacts/chrombert.MEF.mm10.npz
#
# Steps (all on the mirror): make_dataset (overlap peaks onto the 1 kb grid) ->
# get_region_emb (768-d TRN embedding per region) -> fetch hdf5 + tsv back ->
# hdf5_to_artifact.py converts to the .npz embedding-artifact contract.
set -euo pipefail
PEAKS="" GENOME="" STATE="" OUT=""
while [ $# -gt 0 ]; do case "$1" in
  --peaks) PEAKS="$2"; shift 2;;
  --genome) GENOME="$2"; shift 2;;
  --cell-state) STATE="$2"; shift 2;;
  --out) OUT="$2"; shift 2;;
  *) echo "unknown arg $1" >&2; exit 2;;
esac; done
: "${PEAKS:?--peaks required}"; : "${GENOME:?--genome required}"
: "${STATE:?--cell-state required}"; : "${OUT:?--out required}"
source "$(dirname "$0")/mirror_env.sh"

REMOTE="/tmp/ecr_nav_${STATE}_${GENOME}"
base="$(basename "$PEAKS")"
echo ">> uploading $PEAKS -> mirror:$REMOTE/"
$SSH "mkdir -p $REMOTE"
$SCP_BASE "$PEAKS" "${MIRROR_USER}@${MIRROR_IP}:$REMOTE/$base"

echo ">> make_dataset ($GENOME) + get_region_emb on mirror (GPU)"
mirror_py "cd $REMOTE && \
  python -m chrombert.scripts.chrombert_make_dataset $base -g $GENOME -o dataset.tsv && \
  python -m chrombert.scripts.chrombert_get_region_emb dataset.tsv -g $GENOME -o emb.hdf5"

echo ">> fetching results"
mkdir -p "$(dirname "$OUT")"
tmp="$(mktemp -d)"
$SCP_BASE "${MIRROR_USER}@${MIRROR_IP}:$REMOTE/dataset.tsv" "$tmp/dataset.tsv"
$SCP_BASE "${MIRROR_USER}@${MIRROR_IP}:$REMOTE/emb.hdf5"   "$tmp/emb.hdf5"

echo ">> converting to embedding artifact -> $OUT"
python "$(dirname "$0")/hdf5_to_artifact.py" \
  --hdf5 "$tmp/emb.hdf5" --dataset "$tmp/dataset.tsv" \
  --genome "$GENOME" --cell-state "$STATE" --out "$OUT"
rm -rf "$tmp"
echo ">> done: $OUT"
