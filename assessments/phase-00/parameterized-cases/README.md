# Phase 00 parameterized validation cases

Unseen inputs for the computational problems, so answers cannot be copied. Each case is a JSON file `{"lesson": "P00.NN", "seed": int, "inputs": {...}, "expected": {...}, "tolerance": {...}}`. Generate new ones with `generate.py` (which uses the reference solutions) and check a learner's implementation with `run.py`.

Cases shipped:

- `p00-03-isa-*.json`: altitude arrays → temperature, pressure, density, speed of sound (rtol 1e-4).
- `p00-14-propagation-*.json`: function spec and input uncertainties → combined standard uncertainty (rtol 1e-3).
- `p00-15-freefall-*.json`: heights, synthetic timing sets → mean times, uncertainties, normalized error (rtol 1e-3).

These are written once the reference solutions for lessons 03, 14 and 15 are merged; `generate.py` then produces them deterministically from a seed.
