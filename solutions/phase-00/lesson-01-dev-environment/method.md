# P00.01  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. For each string, take the leading digits-and-dots run, split on ".", convert each piece to an integer. 2. For each comparison, parse both, pad to equal length with zeros, walk left to right to the first differing position; the larger integer there wins. Mark that position.

## 2.
1. State the claim. 2. Scenario A: two projects pinning incompatible versions of the same package; what breaks when both are global. 3. Scenario B: reproducing a result after upgrades; what `pip freeze` in a per-project environment gives you that a global install cannot. 4. Concede the strongest counterpoint (single-project machine, throwaway scripts) and bound when it is acceptable.

## 3.
1. Compare the Python shown with the one in `.venv`; note the editor ran a different interpreter. 2. Name the fix (select interpreter / activate). 3. Observe where `import numpy` sits in the file and conclude the table could not have printed before the traceback. 4. Describe the only code change that would make the transcript possible (moving the NumPy import inside the smoke-test function). 5. Say why the lesson's layout (import at top) is the safer design anyway.

## 4.
1. Copy `REQUIREMENTS`, add the two entries using the distribution names. 2. Run `check_packages(reqs)` and `format_report`. 3. Try the import name `yaml` as the key and observe `missing`. 4. Explain the two namespaces (distribution name on the index vs module name in `import`).

## 5.
1. Loop: `total = 0.0; for _ in range(10_000_000): total += 0.1`. 2. `np.full(10_000_000, 0.1).sum()`. 3. Relative error = |result − 1e6| / 1e6 for each. 4. Explain pairwise summation: error grows like log₂ n instead of n.

## 6.
1. Candidate A: different BLAS/LAPACK or CPU instruction set (last-bit differences in matrix products and `sum`). 2. Candidate B: different random seeds, operating-system line endings in data files, or locale decimal separators when reading CSV. 3. For each: order of magnitude of the discrepancy; how it grows (or does not) through an atmosphere table, a one-hour integration, a 1.5 safety margin. 4. Propose: print `sys.executable` and `numpy.show_config()`, report results with uncertainty and tolerances, record the exact environment with `pip freeze`.
