#!/usr/bin/env bash
# Compute the aTPM (normalized per-peak ATAC) channel for GET's region_motif input.
# Species-AGNOSTIC: aTPM depends only on bigWig + peaks, not the genome assembly —
# works for human or mouse states alike. Runs on the ATAC data server (where the
# bigWigs live); the small output table is transferred to the GET model instance.
#
# GET input per region = [282 motif scores | 1 aTPM] = 283. aTPM is the CELL-STATE
# channel: it makes state-A vs state-B embeddings differ (peaks open in one / closed
# in the other) = the driver signal. GET trains aTPM on [0,1]; we normalize per state
# to that range (99th-pct, robust to outlier peaks).
#
# Usage (two states A vs B, e.g. MEF vs mES, or two human cell types):
#   compute_atpm.sh <A_name> <A.bw> <A_peaks.bed> <B_name> <B.bw> <B_peaks.bed> <outdir>
# Output: <outdir>/atpm_union.tsv  (chrom,start,end,atpm_<A>,atpm_<B>) + union_named.bed
# Tools: bedtools, bigWigAverageOverBed, python+numpy.
set -euo pipefail
AN=${1:?A name}; ABW=${2:?A.bw}; APK=${3:?A peaks}
BN=${4:?B name}; BBW=${5:?B.bw}; BPK=${6:?B peaks}
OUT=${7:-/tmp/ecr_atpm}
mkdir -p "$OUT"; cd "$OUT"

# 1. union of both states' peaks (the common region set both states are scored on)
cat "$APK" "$BPK" | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -i - > union.bed
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$1":"$2"-"$3}' union.bed > union_named.bed
echo "union peaks: $(wc -l < union_named.bed)"

# 2. per-peak mean signal (mean0 = sum/size, col 5) for each state
bigWigAverageOverBed "$ABW" union_named.bed A.tab
bigWigAverageOverBed "$BBW" union_named.bed B.tab

# 3. normalize each state to [0,1] (99th pct), join by peak name -> table
AN="$AN" BN="$BN" python - <<'PY'
import os, numpy as np
AN, BN = os.environ["AN"], os.environ["BN"]
def load(t):
    d={}
    for l in open(t):
        f=l.split("\t"); d[f[0]]=float(f[4])
    return d
a,b=load("A.tab"),load("B.tab")
names=[l.split("\t")[3].strip() for l in open("union_named.bed")]
av=np.array([a[n] for n in names]); bv=np.array([b[n] for n in names])
def norm(x):
    p=np.percentile(x[x>0],99) if (x>0).any() else 1.0
    return np.clip(x/p,0,1)
an,bn=norm(av),norm(bv)
with open("atpm_union.tsv","w") as o:
    o.write(f"chrom\tstart\tend\tatpm_{AN}\tatpm_{BN}\n")
    for n,x,y in zip(names,an,bn):
        c,se=n.split(":"); s,e=se.split("-")
        o.write(f"{c}\t{s}\t{e}\t{x:.5f}\t{y:.5f}\n")
print("wrote atpm_union.tsv:", len(names), "regions;",
      int((abs(an-bn)>0.3).sum()), "differ >0.3 (candidate drivers)")
PY
echo "done -> $OUT/atpm_union.tsv (+ union_named.bed)"
