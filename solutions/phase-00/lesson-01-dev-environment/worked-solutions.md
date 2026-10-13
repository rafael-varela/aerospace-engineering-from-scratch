# P00.01  Worked solutions

Each solution states: assumptions, the reasoning, the result, at least one independent check, and the common mistakes.

## 1. Versions by hand

**Result.** `"3.11.4"` → (3, 11, 4); `"1.26"` → (1, 26); `"2.0.0rc1"` → (2, 0, 0); `"9.1.1"` → (9, 1, 1); `"0.23.post1"` → (0, 23).
Comparisons (padded): `"1.26.4" >= "1.26"`: (1, 26, 4) vs (1, 26, 0), decided at position 3, 4 ≥ 0, **true**. `"1.9.0" >= "1.26"`: (1, 9, 0) vs (1, 26, 0), position 2, 9 < 26, **false**. `"3.11" >= "3.11.0"`: (3, 11, 0) vs (3, 11, 0), no differing position, equal, **true**. `"10.0" >= "9.9"`: (10, 0) vs (9, 9), position 1, 10 > 9, **true**.
**Check.** Run `parse_version` and `meets_minimum` on each; the tests in `tests/test_main.py` cover exactly these cases.
**Common mistakes.** Treating the strings as decimals (1.26 < 1.9 as numbers, but version 1.26 is newer); keeping `rc1` as a component; forgetting to pad, which makes (1, 26) < (1, 26, 0) under Python's tuple ordering because a shorter prefix sorts first.

## 2. The case for isolation

**Argument (one acceptable version).** A global install works until two things need different versions of the same package. Project A, written in 2023, pins NumPy 1.26 because it relies on behaviour that NumPy 2 changed (for example, the printing of scalars and some default integer types); the course needs NumPy 2.5 for a newer SciPy. Installed globally, the second install silently overwrites the first and project A starts producing different output with no change to its code. A per-project environment keeps each on its own copy. The second scenario is time, not conflict: a drag number reported today must be reproducible when a reviewer asks in a year. `pip freeze > requirements-lock.txt` inside the project's environment records the exact versions that produced the number; a global environment that has been upgraded twenty times since has no such record. Isolation costs one command to create and one to activate.
**The colleague's strongest point.** On a machine used for exactly one project, or for throwaway one-off scripts whose output will never be cited, isolation adds a step for no benefit, and that is fine *as long as* nothing from that machine is ever reported as a result. The moment a number leaves the machine, it needs its environment recorded.
**Common mistakes.** Arguing only from authority ("best practice says so") rather than from a concrete failure; failing to concede the counterpoint, which makes the argument weaker, not stronger.

## 3. Diagnose from a transcript

**Diagnosis.** The editor's Run button used the system interpreter, Python 3.12.3, in which none of the course packages are installed; the terminal was using the activated `.venv`. The fix is one setting: "Python: Select Interpreter" in VS Code and choose the `.venv` interpreter (equivalently, activate the environment in the editor's terminal). The report itself, if it printed `sys.executable`, would have shown the two paths differ.
**Why the transcript is impossible as shown.** `import numpy as np` is at the top of `main.py`. If NumPy were missing, that import would raise `ModuleNotFoundError` immediately and *nothing* would print: no table, no header. A table followed by a traceback can only happen if the NumPy import is executed *after* the table is printed, which means someone moved `import numpy` inside `speed_of_sound_smoke()` (a "lazy import"). That design is actually reasonable for a check script, because it lets the report print before the smoke test fails, but then the smoke test should catch `ImportError` and print a FAIL line instead of a traceback.
**Check.** Reproduce it: deactivate the environment and run the script; observe that with the import at the top you get only the traceback.
**Common mistakes.** Blaming the student for "not installing NumPy" when `pip list` in the terminal shows it present; not noticing that the printed traceback contradicts the file layout.

## 4. Extend the check

**Result.** With `reqs = {**REQUIREMENTS, "PyYAML": "6.0", "requests": "2.31"}` the report gains two lines, both OK on a machine where `requirements.txt` was installed (observed on the author's machine: PyYAML 6.0.1, requests 2.34.2; your numbers will differ). With the key `"yaml"` instead of `"PyYAML"`, the line reads `yaml  missing  6.0  FAIL`.
**Explanation.** `importlib.metadata.version` searches the installed *distributions*, the names under which packages are published and which `pip list` prints (`PyYAML`, case-insensitive). `import yaml` uses the *module* name inside the distribution. Usually they coincide (`numpy`/`numpy`); when they do not, the check must use the distribution name, and the lesson's `REQUIREMENTS` is written that way.
**Check.** `python -m pip show PyYAML` prints the distribution name and version that `importlib.metadata` reports.
**Common mistakes.** Concluding that PyYAML is missing because `yaml` reported `missing`.

## 5. Floating point on your machine

**Result (author's machine; yours will differ in the trailing digits only).** `0.1 + 0.2` → `0.30000000000000004`, so `== 0.3` is `False`. Plain loop of 10⁷ additions of 0.1: 999 999.9998389754 (relative error 1.6 × 10⁻¹⁰). `np.float64` sum: 1 000 000.0 exactly on this machine (relative error below 10⁻¹⁶; on another machine expect something like 10⁻¹⁵). `np.float32` sum: 1 000 000.125 (relative error 1.3 × 10⁻⁷). Exact: 1 000 000.
**Explanation.** Each addition rounds to the nearest representable value, an error of about half a unit in the last place of the *result*. In the loop, the running total grows to 10⁶ while each addend is 0.1, so each of the 10⁷ roundings is about 10⁶ × 2⁻⁵³ ≈ 10⁻¹⁰ and they accumulate roughly linearly: 10⁷ × 10⁻¹⁷ relative to the total gives about 10⁻¹⁰. NumPy's `sum` uses pairwise summation: it adds neighbours, then pairs of sums, and so on; the partial sums stay small for most of the additions and the error grows like log₂(10⁷) ≈ 23 roundings, not 10⁷.
**Check.** `math.fsum([0.1] * 10_000_000)` gives the correctly rounded exact sum, 1 000 000.0, and is the reference.
**Common mistakes.** Calling the loop result "wrong"; it is the correct floating-point result of that sequence of operations. The lesson is to choose the sequence (and the tolerance) deliberately.

## 6. False confidence

**Situation A: different compiled backends.** Same NumPy version, one linked against OpenBLAS, the other against Apple Accelerate or Intel MKL. Matrix products, `sum`, `exp` and `sqrt` can differ in the last one or two bits: relative 10⁻¹⁶ to 10⁻¹⁵. Atmosphere table: irrelevant (eight figures are more than any instrument). One-hour trajectory integration at 10 Hz (36 000 steps): rounding differences accumulate and can be amplified by a sensitive (chaotic) dynamics; expect 10⁻¹² to 10⁻⁹ relative, still irrelevant for engineering unless the trajectory is near a separatrix (a knife-edge between two outcomes), where the lesson is that the *physics* is sensitive, not that a machine is wrong. Safety margin 1.5: irrelevant.
**Situation B: data-reading differences.** A CSV written on a machine whose locale uses a decimal comma, or with Windows line endings read by a parser that does not expect them; a random number generator without a fixed seed. These produce *large* differences (a parsed `1,225` becoming 1225, or a Monte Carlo spread of several percent) and the environment check is blind to all of them. For an atmosphere table a misparsed constant is a factor-of-1000 error; for a trajectory, fatal; for a 1.5 safety margin, decisive.
**What to add.** Print `sys.executable`, `platform.platform()` and `numpy.show_config()` in the report; pin exact versions in a lock file; set and record random seeds; read data with explicit `encoding`, `newline` and decimal settings (`P00.08`); and report every number with its uncertainty and the tolerance it was tested against, so a last-bit difference is seen for what it is.
**Common mistakes.** Treating all differences as equally alarming, or dismissing all of them because "floating point is always a bit off".

## Check-yourself answers

1. `import sys; print(sys.executable)`.
2. True, because strings compare character by character and `"9" > "2"`. Versions must be compared as tuples of integers: (1, 9) < (1, 26).
3. So the decision logic can be tested with a fake lookup covering present, outdated and missing packages, independently of what is actually installed (dependency injection).
4. `False`; machine epsilon for `float64` is about 2 × 10⁻¹⁶.
5. Any of: that two machines agree to the last digit; that the editor uses the same interpreter as the terminal; that compiled backends are identical; that data files will be read the same way.
6. `python -m pip` runs pip belonging to the interpreter named `python` in the active environment; a bare `pip` can belong to a different interpreter on the path.
7. a ∝ √R, so a 0.16 % error in R (287.5 vs 287.05) gives about 0.08 % in a, roughly 0.27 m/s. Larger than 0.05 m/s, so yes, caught.
