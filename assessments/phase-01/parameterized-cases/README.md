# Phase 01 parameterized validation cases

Unseen inputs for the computational problems, so answers cannot be copied. Each case is a JSON file `{"lesson": "P01.NN", "seed": int, "inputs": {...}, "expected": {...}, "tolerance": {...}}`. Generate new ones with `generate.py` (which uses the reference solutions) and check a learner's implementation with `run.py`. Expected values for the linear-algebra lessons come from `numpy.linalg` as the oracle, and the tolerance is set so that a correct from-scratch implementation passes but a wrong pivot order, a sign error, or a missing normalisation does not.

Cases shipped:

- `p01-01-derivative-*.json`: a function spec, evaluation points and step sizes → forward, central and symbolic derivatives, and the step that minimises the central-difference error (rtol 1e-6 on derivatives; the optimal step must be within one decade).
- `p01-02-glide-*.json`: C_D0, K, mass, wing area, altitude → C_L*, (L/D)_max, best-glide speed, still-air range (rtol 1e-4).
- `p01-05-taylor-*.json`: angle range and polynomial degree → remainder bound and actual error of the small-angle approximations (rtol 1e-6; the actual error must not exceed the bound).
- `p01-13-elimination-*.json`: random well-conditioned n × n systems (n = 3–8), with and without a required row exchange → solution, pivots, determinant (rtol 1e-8; the determinant sign must be right).
- `p01-16-eigen-*.json`: random symmetric 3 × 3 inertia tensors → principal moments in ascending order and the principal-axis matrix (rtol 1e-8 on eigenvalues; eigenvectors compared up to sign through |vᵢ · v̂ᵢ| = 1 to 1e-8; the matrix must be orthogonal with determinant +1).
- `p01-17-polar-fit-*.json`: synthetic wind-tunnel C_L, C_D arrays with noise → C_D0, K, residual norm, (L/D)_max (rtol 1e-6 against `lstsq`; the residual must be orthogonal to the design columns to 1e-10).
- `p01-18-conditioning-*.json`: a stiffness matrix and a perturbed load → singular values, condition number, and the relative displacement error compared with the κ‖δb‖/‖b‖ bound (rtol 1e-8; the observed error must not exceed the bound).

These are written once the reference solutions for lessons 01, 02, 05, 13, 16, 17 and 18 are merged; `generate.py` then produces them deterministically from a seed.
