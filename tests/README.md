# tests

Pure-stdlib + numpy tests for the navigator core (no torch, runs in the light env).
pytest is not required; each module self-runs:

```bash
python tests/test_direction.py        # direction channel: signed_delta, attach_direction, alignment, contract
python tests/test_measured_signal.py  # attach_measured_signal.map_signal / read_intensity
```

They are also plain `test_*` functions with bare asserts, so `pytest tests/` works
if pytest is added. A test decorated `@xfail(...)` encodes desired-but-unimplemented
behaviour (a known issue); it is expected to fail today and will shout `XPASS` when
the fix lands so it gets un-marked. See `test_uncovered_region_should_not_be_a_measured_zero`
(review Finding 1).
