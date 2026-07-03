"""
eCR_navigator — driver-vs-passenger genomic region weighting for MEF→iPSC
reprogramming, and target nomination for engineered chromatin regulators.

Front-end of the eCR platform. Produces per-region driver scores that weight
off-target severity in eCR_predictor (Tier 2) and nominate target loci. The
stable output is the region-weight contract (see ecr_navigator.weights and
docs/region_weight_contract.md); the model that produces it (ecr_navigator.model)
is an open design decision — see CLAUDE.md.
"""
__version__ = "0.0.0"
