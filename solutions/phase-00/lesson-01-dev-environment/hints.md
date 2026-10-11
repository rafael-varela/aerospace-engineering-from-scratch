# P00.01  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Pad the shorter tuple with zeros, then find the first position where the two tuples differ. If they never differ, they are equal.
2. Think about what "it works" means a year later when the global NumPy has moved two major versions.
3. In `main.py`, `import numpy as np` is at the top of the file. Could the table have printed if that line had failed?
4. `importlib.metadata.version` looks up the name under which the package was *installed*, which is what `pip list` shows. Compare `pip list` with the name you type after `import`.
5. NumPy's `sum` adds the numbers in pairs (pairwise summation); a loop adds them to one growing total. Think about how big the rounding error of each addition is relative to the running total.
6. The version number describes the Python code; it says nothing about the compiled libraries underneath (BLAS, LAPACK, the C math library), the CPU, or the operating system.
