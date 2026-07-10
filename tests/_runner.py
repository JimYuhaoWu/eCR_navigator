"""
Minimal stdlib test runner — pytest is not in the light navigator env, but these
test modules are written as plain `test_*` functions with bare asserts so they are
*also* pytest-discoverable if pytest is later added. Run a module directly:

    python tests/test_direction.py
    python tests/test_measured_signal.py

`@xfail(reason)` marks a test that encodes DESIRED (not-yet-implemented) behaviour —
it is expected to fail today (a living TODO for a known issue). If such a test starts
passing, the runner shouts XPASS so you remember to un-mark it.
"""
from __future__ import annotations

import os
import sys
import traceback


def xfail(reason):
    def deco(fn):
        fn._xfail_reason = reason
        return fn
    return deco


def run(namespace):
    tests = {n: f for n, f in namespace.items()
             if n.startswith("test_") and callable(f)}
    passed = failed = xfailed = xpassed = 0
    failures = []
    for name in sorted(tests):
        fn = tests[name]
        reason = getattr(fn, "_xfail_reason", None)
        try:
            fn()
        except AssertionError as e:
            if reason:
                xfailed += 1
                print(f"xfail  {name}  ({reason})")
            else:
                failed += 1
                failures.append(name)
                print(f"FAIL   {name}: {e}")
        except Exception as e:  # noqa: BLE001 — surface any unexpected error
            failed += 1
            failures.append(name)
            print(f"ERROR  {name}: {type(e).__name__}: {e}")
            traceback.print_exc()
        else:
            if reason:
                xpassed += 1
                print(f"XPASS  {name}  (known issue now passes — remove @xfail: {reason})")
            else:
                passed += 1
                print(f"ok     {name}")
    print(f"\n{passed} passed, {failed} failed, {xfailed} xfail, {xpassed} xpass")
    sys.exit(1 if failed else 0)


def add_repo_paths():
    """Put the repo root (for `ecr_navigator`) and scripts/ (for the mirror
    scripts' sibling imports) on sys.path, regardless of CWD."""
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)
    for p in (root, os.path.join(root, "scripts")):
        if p not in sys.path:
            sys.path.insert(0, p)
