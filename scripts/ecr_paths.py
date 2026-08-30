"""Where the model data lives — resolved once, overridable by environment.

The mirrors bake `/yutiancheng/yuhao/...` into conda `conda-meta` records and
script shebangs, so those absolute paths cannot simply be rewritten: on the SH
(DCU) cluster the shared tree is bind-mounted to `/yutiancheng` so the original
paths still resolve. What CAN be fixed is our own hardcoded defaults, so the
pipeline can run wherever the tree happens to be mounted.

Set ECR_DATA_ROOT to relocate. The default is the historical path, so existing
invocations behave exactly as before.

    export ECR_DATA_ROOT=/data/ecr_staging/shared/yuhao

Deliberately imports nothing beyond the stdlib: it is imported by scripts that
run in the light navigator env (no torch), and by the tests.
"""
from __future__ import annotations

import os
import posixpath

DEFAULT_DATA_ROOT = "/yutiancheng/yuhao"
ENV_VAR = "ECR_DATA_ROOT"


def data_root() -> str:
    """Root of the shared model/data tree. Trailing slashes are stripped so
    join() results are stable."""
    return (os.environ.get(ENV_VAR) or DEFAULT_DATA_ROOT).rstrip("/") or "/"


def under(*parts: str) -> str:
    """Path beneath the data root, e.g. under('models', 'get') ->
    '/yutiancheng/yuhao/models/get'.

    Uses posixpath, not os.path: these are always POSIX paths on the cluster,
    and os.path.join would emit backslashes when this module is imported on
    Windows (where the repo is edited).
    """
    return posixpath.join(data_root(), *parts)


def model(*parts: str) -> str:
    """Path beneath <data_root>/models."""
    return under("models", *parts)
