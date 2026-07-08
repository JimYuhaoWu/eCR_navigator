#!/usr/bin/env python
"""
Dump an AlphaGenome output-track ontology (per organism, per modality) to TSV — used
to pick which predicted tracks stand in for a cell state (e.g. mouse DNase NIH3T3 for
MEF, ES-E14 for mES).

AlphaGenome serves track metadata over its gRPC API (gdmscience.googleapis.com), so
this runs where that endpoint is reachable — NOT the model-zoo mirror (which can't
reach Google). Needs the `alphagenome` client (Python >= 3.10) and an API key.

Proxy gotcha: the client's gRPC channel honors the LOWERCASE `https_proxy`/`grpc_proxy`
env vars only. If your machine proxies HTTPS via uppercase `HTTPS_PROXY` (e.g. a local
Clash on 127.0.0.1:7890), `create()` hangs forever on channel-ready. Export the
lowercase vars too:

  ALPHAGENOME_API_KEY=... https_proxy=http://127.0.0.1:7890 grpc_proxy=http://127.0.0.1:7890 \
      python alphagenome_dump_tracks.py --organism mouse --out-dir .

Each `output_metadata(organism).<modality>` is a pandas DataFrame with columns incl.
`biosample_name`, `ontology_curie`, `biosample_life_stage`, `Assay title`.
"""
from __future__ import annotations

import argparse
import os
import sys

from alphagenome.models import dna_client

ORGANISMS = {
    "human": dna_client.Organism.HOMO_SAPIENS,
    "mouse": dna_client.Organism.MUS_MUSCULUS,
}
# keywords flagged as likely MEF / mES proxies when scanning a modality
FLAG = ["fibroblast", "3t3", "embryonic stem", "es-e14", "es-", "stem", "pluripotent",
        "ipsc", "epiblast"]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--api-key", default=os.environ.get("ALPHAGENOME_API_KEY"))
    ap.add_argument("--organism", choices=list(ORGANISMS), default="mouse")
    ap.add_argument("--modalities", nargs="+", default=["atac", "dnase"])
    ap.add_argument("--out-dir", default=".")
    a = ap.parse_args()
    if not a.api_key:
        sys.exit("need --api-key or ALPHAGENOME_API_KEY env var")

    client = dna_client.create(a.api_key)
    md = client.output_metadata(organism=ORGANISMS[a.organism])

    for name in a.modalities:
        df = getattr(md, name, None)
        if df is None:
            print(f"{name}: not present for {a.organism}")
            continue
        outp = os.path.join(a.out_dir, f"ag_{a.organism}_{name}.tsv")
        df.to_csv(outp, sep="\t", index=True)
        print(f"\n{name.upper()}: {len(df)} {a.organism} tracks -> {outp}")
        low = df.astype(str).apply(lambda r: " ".join(r).lower(), axis=1)
        for kw in FLAG:
            hits = df[low.str.contains(kw, regex=False)]
            if len(hits):
                bios = sorted(hits["biosample_name"].astype(str).unique())[:6]
                print(f"  [{kw}] {len(hits)} track(s): idx {list(hits.index)[:6]} e.g. {bios}")


if __name__ == "__main__":
    main()
