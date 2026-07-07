#!/usr/bin/env bash
# Compute the aTPM (normalized per-peak ATAC) channel for GET's region_motif input.
# Runs on the ATAC DATA server (peilab2, wuyuhao@172.16.78.234:2020), where the
# bigWig signal tracks live. Output is a small per-peak table transferred to the
# GPU model instance and combined with the motif matrix (get_regionmotif_matrix.py).
#
# GET input per region = [282 motif scores | 1 aTPM] = 283. aTPM is the CELL-STATE
# channel: it makes MEF vs mES embeddings differ (peaks open in one state / closed
# in the other), which is the driver signal. GET trains aTPM on [0,1] (demo zarr:
# median ~0.02, max 1.0), so we normalize per state to that range.
#
# Inputs (symlinked under .../MEF_mESC_3/preprocess/):
#   Bw/MEF.bw, Bw/mES.bw          # ATAC signal tracks (mm10)
#   Peak/MEF.e7_peaks.bed, Peak/mES.e7_peaks.bed
# Tools on the data server: bedtools, bigWigAverageOverBed, python+numpy.
set -euo pipefail
PRE=/mnt3/wuyuhao/MEF_mESC_3/preprocess
OUT=${1:-/tmp/ecr_atpm}
mkdir -p "$OUT"; cd "$OUT"

# 1. union of both states' e7 peaks (the common region set both states are scored on)
cat "$PRE"/Peak/MEF.e7_peaks.bed "$PRE"/Peak/mES.e7_peaks.bed | cut -f1-3 \
  | sort -k1,1 -k2,2n | bedtools merge -i - > union.bed
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$1":"$2"-"$3}' union.bed > union_named.bed
echo "union peaks: $(wc -l < union_named.bed)"

# 2. per-peak mean signal (mean0 = sum/size, col 5) for each state
for s in MEF mES; do
  bigWigAverageOverBed "$PRE/Bw/$s.bw" union_named.bed "$s.tab"
done

# 3. normalize each state to [0,1] by its 99th percentile (robust to outlier peaks),
#    join by peak name, emit chrom,start,end,atpm_MEF,atpm_mES
python - <<'PY'
import numpy as np
def load(t):
    d={}
    for l in open(t):
        f=l.split("\t"); d[f[0]]=float(f[4])
    return d
mef,mes=load("MEF.tab"),load("mES.tab")
names=[l.split("\t")[3].strip() for l in open("union_named.bed")]
mv=np.array([mef[n] for n in names]); sv=np.array([mes[n] for n in names])
def norm(x):
    p=np.percentile(x[x>0],99) if (x>0).any() else 1.0
    return np.clip(x/p,0,1)
mn,sn=norm(mv),norm(sv)
with open("atpm_union.tsv","w") as o:
    o.write("chrom\tstart\tend\tatpm_MEF\tatpm_mES\n")
    for n,a,b in zip(names,mn,sn):
        c,se=n.split(":"); s,e=se.split("-")
        o.write(f"{c}\t{s}\t{e}\t{a:.5f}\t{b:.5f}\n")
print("wrote atpm_union.tsv:", len(names), "regions;",
      int((abs(mn-sn)>0.3).sum()), "differ >0.3 (candidate drivers)")
PY
echo "done -> $OUT/atpm_union.tsv (+ union_named.bed)"
