"""Runtime introspection: which accelerator are we on, and what can it do.

Deliberately stdlib-only. These helpers are imported by scripts that run in the
light navigator env (which has no torch) and by the tests, so torch is passed in
or imported lazily — never at module scope.

Background: the models were developed on NVIDIA A800s and now also run on Hygon
DCUs (gfx928) under DTK. The two report their capabilities differently, and
results differ slightly, so both facts need handling explicitly rather than
inferred from CUDA-shaped assumptions.
"""
from __future__ import annotations

import platform


def flash_attn_supported(has_flash_attn: bool, hip_version, capability):
    """Decide whether to take the flash-attention path.

    Pure function: every runtime fact is passed in, so the logic is testable
    without a GPU. Returns (bool, reason).

    has_flash_attn  whether `import flash_attn` succeeded
    hip_version     torch.version.hip  (None on NVIDIA / CPU)
    capability      torch.cuda.get_device_capability() or None

    On ROCm/HIP, `get_device_capability()` reports the *gfx* version, NOT an
    NVIDIA SM number: gfx928 comes back as (9, 2). So the older check
    `capability[0] >= 8` passed on DCU only by coincidence, and would also pass
    on, say, a gfx906 card with no usable flash-attention build. On that
    backend the honest signal is whether a working `flash_attn` is importable,
    so branch on the backend first.
    """
    if not has_flash_attn:
        return False, "flash_attn not importable"
    if hip_version:
        return True, f"HIP {hip_version}: flash_attn present"
    if capability is None:
        return False, "no CUDA device"
    ok = capability[0] >= 8                       # Ampere (sm_80) or newer
    return ok, f"CUDA sm_{capability[0]}{capability[1]}"


def detect_flash_attn(torch_mod=None):
    """flash_attn_supported() wired to the live torch runtime."""
    if torch_mod is None:
        import torch as torch_mod                 # noqa: PLC0415
    try:
        import flash_attn                         # noqa: F401,PLC0415
        has = True
    except Exception:
        has = False
    cap = (torch_mod.cuda.get_device_capability()
           if torch_mod.cuda.is_available() else None)
    return flash_attn_supported(has, torch_mod.version.hip, cap)


def runtime_provenance(torch_mod=None):
    """Describe the runtime that produced an embedding.

    DCU and A800 results differ slightly and systematically, so every artifact
    records how it was made. The *dtype* and *attention path* matter as much as
    the accelerator: measured on EpiAgent, like-for-like fp32 agreement between
    an A800 and a K100_AI is ~1.7e-06, while switching to fp16 autocast or a
    different flash-attention kernel each shift results by ~1e-01 at the tail.

    torch_mod is injectable for testing. When torch is absent (the light
    navigator env has none) this degrades to a host-only record rather than
    raising, so importing this module never requires torch.
    """
    prov = {"host": platform.node(), "torch": None, "torch_build": None,
            "device": "cpu", "accelerator": None, "arch": None,
            "hip": None, "cuda": None}
    if torch_mod is None:
        try:
            import torch as torch_mod             # noqa: PLC0415
        except Exception:
            return prov
    prov["torch"] = getattr(torch_mod, "__version__", None)
    # torch.__version__ drops the local build tag: the DTK build reports plain
    # "2.4.1", indistinguishable from stock torch. The distribution metadata
    # keeps it ("2.4.1+das.opt1.dtk25041"), which is the part that identifies
    # the build actually used.
    try:
        from importlib.metadata import version           # noqa: PLC0415
        prov["torch_build"] = version("torch")
    except Exception:
        prov["torch_build"] = None
    prov["hip"] = getattr(torch_mod.version, "hip", None)
    prov["cuda"] = getattr(torch_mod.version, "cuda", None)
    try:
        if torch_mod.cuda.is_available():
            props = torch_mod.cuda.get_device_properties(0)
            prov["device"] = "cuda"
            prov["accelerator"] = torch_mod.cuda.get_device_name(0)
            prov["arch"] = (getattr(props, "gcnArchName", None)
                            or f"sm_{props.major}{props.minor}")
    except Exception:
        pass                                      # keep the host-only record
    return prov
