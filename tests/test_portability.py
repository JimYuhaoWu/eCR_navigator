"""
Tests for the portability work: ECR_DATA_ROOT path resolution, backend-aware
flash-attention detection, and the provenance block added to the embedding
artifact contract.

All three are written so they run in the LIGHT navigator env, which has no
torch: the path helper is stdlib-only, the flash-attn decision is a pure
function fed the runtime facts, and runtime_provenance() takes an injectable
torch module. Real-hardware values are verified separately on the cluster.
"""
from __future__ import annotations

import json
import os
import tempfile
import types

import numpy as np

TAB = chr(9)
NL = chr(10)

from _runner import add_repo_paths, run

add_repo_paths()

import ecr_paths                                             # noqa: E402  (scripts/)
from ecr_runtime import (dtype_tag,                           # noqa: E402
                          flash_attn_supported, runtime_provenance)
from embedding_artifact import write_embedding_artifact      # noqa: E402


# --------------------------------------------------------------- ecr_paths

def _with_root(value):
    """Set/clear ECR_DATA_ROOT for one call; returns the previous value."""
    prev = os.environ.get(ecr_paths.ENV_VAR)
    if value is None:
        os.environ.pop(ecr_paths.ENV_VAR, None)
    else:
        os.environ[ecr_paths.ENV_VAR] = value
    return prev


def test_data_root_defaults_to_historical_path():
    prev = _with_root(None)
    try:
        assert ecr_paths.data_root() == "/yutiancheng/yuhao"
    finally:
        _with_root(prev)


def test_data_root_honours_env_var():
    prev = _with_root("/data/ecr_staging/shared/yuhao")
    try:
        assert ecr_paths.data_root() == "/data/ecr_staging/shared/yuhao"
    finally:
        _with_root(prev)


def test_data_root_strips_trailing_slash():
    prev = _with_root("/mnt/somewhere/")
    try:
        assert ecr_paths.data_root() == "/mnt/somewhere"
    finally:
        _with_root(prev)


def test_paths_are_posix_even_on_windows():
    """os.path.join emits backslashes when the repo is edited on Windows;
    these are always cluster paths, so they must stay POSIX."""
    prev = _with_root("/yutiancheng/yuhao")
    try:
        p = ecr_paths.model("get", "pretrain_fetal_adult")
        assert p == "/yutiancheng/yuhao/models/get/pretrain_fetal_adult"
        assert "\\" not in p
    finally:
        _with_root(prev)


def test_empty_env_var_falls_back_to_default():
    """An exported-but-empty ECR_DATA_ROOT must not yield a relative path."""
    prev = _with_root("")
    try:
        assert ecr_paths.data_root() == "/yutiancheng/yuhao"
    finally:
        _with_root(prev)


# ------------------------------------------------- flash-attention detection

def test_flash_needs_the_package():
    ok, why = flash_attn_supported(False, None, (8, 0))
    assert ok is False and "not importable" in why


def test_flash_on_hip_when_package_present():
    """gfx928 (DCU): Hygon ships a HIP flash_attn build."""
    ok, why = flash_attn_supported(True, "6.3.25211", (9, 2))
    assert ok is True and "HIP" in why


def test_hip_without_flash_package_is_false():
    """REGRESSION for the old capability[0] >= 8 check: on ROCm,
    get_device_capability() returns the gfx version, so a DCU with no usable
    flash_attn build reported (9, x) and the old logic wrongly said True."""
    ok, _ = flash_attn_supported(False, "6.3.25211", (9, 0))
    assert ok is False


def test_flash_on_ampere():
    ok, why = flash_attn_supported(True, None, (8, 0))
    assert ok is True and "sm_80" in why


def test_no_flash_on_volta():
    ok, why = flash_attn_supported(True, None, (7, 0))
    assert ok is False and "sm_70" in why


def test_no_cuda_device_is_false():
    ok, why = flash_attn_supported(True, None, None)
    assert ok is False and "no CUDA device" in why


# ------------------------------------------------------------- provenance

def _fake_torch(*, hip=None, cuda=None, available=False, name=None,
                gcn=None, major=8, minor=0):
    version = types.SimpleNamespace(hip=hip, cuda=cuda)
    props = types.SimpleNamespace(major=major, minor=minor)
    if gcn is not None:
        props.gcnArchName = gcn
    cuda_ns = types.SimpleNamespace(
        is_available=lambda: available,
        get_device_name=lambda i: name,
        get_device_properties=lambda i: props,
    )
    return types.SimpleNamespace(__version__="2.4.1+test", version=version,
                                 cuda=cuda_ns)


def test_provenance_without_torch_does_not_raise():
    """The light env has no torch; importing and calling must still work."""
    prov = runtime_provenance(torch_mod=None)
    for key in ("host", "torch", "torch_build", "device", "accelerator",
                "arch", "hip", "cuda"):
        assert key in prov


def test_provenance_records_dcu_arch():
    prov = runtime_provenance(_fake_torch(hip="6.3.25211", available=True,
                                          name="K100_AI",
                                          gcn="gfx928:sramecc+:xnack-"))
    assert prov["device"] == "cuda"
    assert prov["accelerator"] == "K100_AI"
    assert prov["arch"] == "gfx928:sramecc+:xnack-"
    assert prov["hip"] == "6.3.25211"


def test_provenance_records_nvidia_sm():
    prov = runtime_provenance(_fake_torch(cuda="12.4", available=True,
                                          name="NVIDIA A800-SXM4-80GB",
                                          major=8, minor=0))
    assert prov["arch"] == "sm_80"
    assert prov["cuda"] == "12.4"
    assert prov["hip"] is None


def test_provenance_cpu_when_no_accelerator():
    prov = runtime_provenance(_fake_torch(available=False))
    assert prov["device"] == "cpu" and prov["accelerator"] is None


# --------------------------------------------- artifact contract regression

def _artifact(**kw):
    n, d = 5, 4
    chrom = np.array(["chr%d" % (i % 2 + 1) for i in range(n)])
    start = np.arange(n, dtype=np.int64) * 1000
    end = start + 500
    emb = np.arange(n * d, dtype=np.float32).reshape(n, d)
    out = os.path.join(tempfile.mkdtemp(), "a.npz")
    write_embedding_artifact(out, chrom, start, end, emb, model="m",
                             cell_state="s", assembly="hg38", source="test", **kw)
    return out, chrom, start, end, emb


def test_arrays_are_unchanged_by_the_provenance_change():
    """The meta block gained a field; the DATA must be byte-identical."""
    out, chrom, start, end, emb = _artifact()
    z = np.load(out)
    assert np.array_equal(z["chrom"], chrom)
    assert np.array_equal(z["start"], start)
    assert np.array_equal(z["end"], end)
    assert np.array_equal(z["embedding"], emb)
    assert z["embedding"].dtype == np.float32


def test_existing_meta_fields_are_preserved():
    """Readers of the old contract must keep working."""
    out = _artifact()[0]
    meta = json.loads(str(np.load(out)["meta"]))
    for key in ("model", "cell_state", "assembly", "dim", "source", "has_signal"):
        assert key in meta
    assert meta["dim"] == 4 and meta["has_signal"] is False


def test_provenance_is_written_into_meta():
    out = _artifact()[0]
    meta = json.loads(str(np.load(out)["meta"]))
    assert "provenance" in meta and "host" in meta["provenance"]


def test_dtype_and_attn_are_recorded_when_supplied():
    out = _artifact(dtype="fp16", attn="flash")[0]
    prov = json.loads(str(np.load(out)["meta"]))["provenance"]
    assert prov["dtype"] == "fp16" and prov["attn"] == "flash"


def test_dtype_and_attn_absent_when_not_supplied():
    """Do not fabricate fields the caller did not assert."""
    out = _artifact()[0]
    prov = json.loads(str(np.load(out)["meta"]))["provenance"]
    assert "dtype" not in prov and "attn" not in prov


def test_caller_supplied_provenance_is_not_mutated():
    given = {"host": "elsewhere"}
    out = _artifact(provenance=given, dtype="fp32")[0]
    assert given == {"host": "elsewhere"}          # no in-place surprise
    prov = json.loads(str(np.load(out)["meta"]))["provenance"]
    assert prov["host"] == "elsewhere" and prov["dtype"] == "fp32"


# --------------------------------------------------------------- dtype tags

class _FakeTensor:
    def __init__(self, dtype):
        self.dtype = dtype


def test_dtype_tag_maps_the_common_dtypes():
    assert dtype_tag(_FakeTensor("torch.float32")) == "fp32"
    assert dtype_tag(_FakeTensor("torch.float16")) == "fp16"
    assert dtype_tag(_FakeTensor("torch.bfloat16")) == "bf16"
    assert dtype_tag(_FakeTensor("torch.float64")) == "fp64"


def test_dtype_tag_passes_through_the_unexpected():
    """An unmapped dtype must still produce something legible, not crash."""
    assert dtype_tag(_FakeTensor("torch.int8")) == "int8"


# ------------------------------------- ChromBERT: provenance through the hdf5

def _fake_chrombert_hdf5(path, provenance=None):
    """Minimal stand-in for chrombert_get_region_emb output."""
    import h5py
    with h5py.File(path, "w") as f:
        f.create_dataset("emb", data=np.zeros((3, 4), dtype=np.float16))
        f.create_dataset("region", data=np.array(
            [[0, 0, 1000, 0], [0, 1000, 2000, 1], [0, 2000, 3000, 2]], dtype=np.int64))
        if provenance is not None:
            f.attrs["ecr_provenance"] = json.dumps(provenance)
    return path


def _run_hdf5_to_artifact(tmp, provenance):
    """Drive hdf5_to_artifact.main() the way run_chrombert_region_emb.sh does."""
    import sys
    import hdf5_to_artifact
    h5 = _fake_chrombert_hdf5(os.path.join(tmp, "emb.hdf5"), provenance)
    tsv = os.path.join(tmp, "dataset.tsv")
    header = ["chrom", "start", "end", "build_region_index", "label"]
    rows = [["chr1", str(i * 1000), str((i + 1) * 1000), str(i), "1"]
            for i in range(3)]
    with open(tsv, "w") as fh:
        fh.write(TAB.join(header) + NL)
        for r in rows:
            fh.write(TAB.join(r) + NL)
    out = os.path.join(tmp, "a.npz")
    argv = sys.argv
    sys.argv = ["hdf5_to_artifact.py", "--hdf5", h5, "--dataset", tsv,
                "--genome", "hg38", "--cell-state", "MEF", "--out", out]
    try:
        hdf5_to_artifact.main()
    finally:
        sys.argv = argv
    return json.loads(str(np.load(out)["meta"]))


def test_chrombert_provenance_survives_the_hdf5_hop():
    """The GPU box stamps the hdf5; the converting host must pass it through
    unchanged rather than describing itself."""
    gpu = {"host": "gpu-box-01", "accelerator": "NVIDIA A800-SXM4-80GB",
           "arch": "sm_80", "cuda": "12.4", "hip": None}
    meta = _run_hdf5_to_artifact(tempfile.mkdtemp(), gpu)
    prov = meta["provenance"]
    assert prov["host"] == "gpu-box-01"          # NOT this machine
    assert prov["accelerator"] == "NVIDIA A800-SXM4-80GB"
    assert prov["arch"] == "sm_80"


def test_chrombert_records_fp16_not_the_cast():
    """emb is stored fp16 and cast to fp32 on read; the tag must say fp16."""
    meta = _run_hdf5_to_artifact(tempfile.mkdtemp(), {"host": "gpu-box-01"})
    assert meta["provenance"]["dtype"] == "fp16"
    assert meta["provenance"]["attn"] == "eager"


def test_chrombert_hdf5_without_provenance_does_not_fabricate_one():
    """Files written before this existed must not be stamped with the
    converting host, which would be a lie about where they came from."""
    meta = _run_hdf5_to_artifact(tempfile.mkdtemp(), None)
    prov = meta["provenance"]
    assert prov.get("accelerator") is None
    assert prov.get("dtype") == "fp16"


if __name__ == "__main__":
    run(globals())
