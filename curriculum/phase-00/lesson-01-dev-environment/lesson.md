# P00.01  Dev Environment for Engineers

> If it does not run on a clean machine, you do not have a result; you have an anecdote.

| | |
|---|---|
| **Learning outcomes** | 1. Install Python 3.11+, create an isolated virtual environment, and install the course's packages into it from `requirements.txt`. 2. Implement and test a script that reports the interpreter and package versions, compares them numerically against minimums, and runs a numerical smoke test. 3. Diagnose the three common environment failures (wrong interpreter, missing package, version too old) from a transcript. 4. Decide whether a machine is fit to reproduce a numerical result, and say what the check does *not* prove. |
| **Prerequisites** | None. This is the first lesson. Skill check: you can open a terminal (Command Prompt or PowerShell on Windows; Terminal on macOS or Linux), type a command and read its output. |
| **Engineering question** | Can this computer reproduce the course's numerical results, and how would we know before trusting a single number it prints? |
| **Assumptions and limits** | The check assumes the version of a package tells you how it behaves, which is true for the public interface and false for compiled numerical libraries underneath it (two installs of the same NumPy version can use different linear-algebra backends and differ in the last digit). It assumes a package that imports also works. It cannot see hardware, operating-system, or locale differences. It proves readiness for *this* course's tests, not for arbitrary software. |
| **Analytical work** | Version comparison as an ordering on tuples of integers (and why string comparison gets it wrong); the floating-point representation of 0.1 and why engineering code compares with tolerances; the ISA sea-level speed of sound as a known answer. |
| **Computational deliverable** | `starter/main.py`: `parse_version`, `meets_minimum`, `check_python`, `check_packages` (with an injectable version lookup so the logic is testable), `speed_of_sound_smoke`, `format_report`, and a `main()` that exits non-zero on failure. Tests prove numeric (not lexicographic) ordering, correct OK/FAIL flags for present, missing and outdated packages, and that *your* environment passes. |
| **Practical activity** | **computational experiment**: break your environment on purpose three ways and record what the check reports each time; measure how 0.1 + 0.2 and a long float32 sum behave on your machine. |
| **Equipment** | A computer running Windows 10+, macOS 12+, or any current Linux, with about 2 GB of free disk and an internet connection for the install. Tier: none. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | The report printed by `python main.py` on your machine; the pytest output; your three failure transcripts with the one-line diagnosis of each; the float experiment table. |
| **Assessment** | `assessments/phase-00/grading-rubric.md`; competencies 2, 3 and 4 (no physics is assessed here). Pass: all public tests pass on your machine; each failure transcript correctly diagnosed; judgment question in the problem set answered with a concrete example. |
| **References** | `@numpy-userguide` "Installation" (required-free); `@pytest-getstarted` "Get Started" > "Install pytest" (required-free). Everything about editors, virtual environments and Git in this lesson is original instruction; the installed tool's own `--help` is the authority on its flags. |
| **Extensions** | ★ Make `main.py` also print the BLAS/LAPACK backend NumPy is using (`numpy.show_config()`) and explain why two machines with identical versions can still differ in the 16th digit. |
| **Time** | unmeasured |

## 1. The engineering problem

Picture a wind-tunnel group that re-runs an analysis script from two years earlier to check a published drag number, and gets a different number. Nobody has changed the script. What has changed is everything around it: a newer version of a numerical library that fixed a rounding bug, a different default integer size on a new machine, a colleague's helper module that was quietly upgraded. The old number is not reproducible, so it is not defensible, and it has to be withdrawn. This scenario is invented, but every element of it happens routinely in engineering groups, which is why reproducibility is now treated as part of the result, not an afterthought.

Every later lesson in this course asks you to compute something and compare it with a prediction or a measurement. That comparison only means something if the computation is **reproducible**: the same inputs, on your machine today and on a reviewer's machine next year, give the same outputs. The first step is to know exactly what your machine is running: which Python, which NumPy, which SciPy, and whether they are new enough. The second step is to isolate the course's packages from everything else on the computer, so that installing some unrelated program next month cannot change your atmosphere table.

The engineering question is: can this computer reproduce the course's numerical results, and how would we know before trusting a single number it prints? Success criteria, fixed before we start: `python --version` reports 3.11 or newer; the environment check reports OK for every package in `requirements.txt`; the smoke test computes the ISA sea-level speed of sound as 340.29 m/s; and `pytest` on this lesson passes. If any of those fail, stop and fix the environment, because every result produced afterwards would be in doubt.

## 2. Understand and derive

### 2.1 Three layers: interpreter, packages, environment

**Python** is an interpreter: a program that reads your `.py` file and executes it line by line. It ships with a standard library (`math`, `csv`, `datetime`, and a few hundred more modules). Everything else the course uses is a **package** installed separately: NumPy (arrays), SciPy (numerical algorithms), Matplotlib (plots), SymPy (algebra), pint (units), pytest (tests). Packages are installed by a tool called `pip` (or the faster `uv`), which downloads them from the Python Package Index and puts them in a folder the interpreter searches when it sees `import numpy`.

That search folder is the problem. If every project on your computer shares one folder, project A needing NumPy 1.26 and project B needing NumPy 2.5 cannot coexist, and upgrading for B silently changes A. A **virtual environment** solves this: it is a private folder containing its own copy of the interpreter's entry point and its own package folder. Activate it and `python` means *that* Python with *those* packages; deactivate it and nothing has changed elsewhere. One environment per project is the rule.

```
  your computer
  ├── system Python 3.12  (leave it alone; the OS may depend on it)
  ├── project-A/.venv/    → Python 3.12 + NumPy 1.26 + ...
  └── aerospace-engineering-from-scratch/.venv/
                          → Python 3.13 + NumPy 2.5 + SciPy 1.18 + ...   ← this course
```

### 2.2 What a version number means, and how to compare two of them

A version like `1.26.4` is three integers: major.minor.patch. Convention (not law) says a patch release fixes bugs without changing behaviour you rely on, a minor release adds features without removing any, and a major release may break things. `requirements.txt` says `numpy>=1.26`: anything from 1.26.0 upward is acceptable, because the course only uses features that have been stable since then.

Here is the trap that this lesson's code must avoid. Versions are written as strings, and Python compares strings character by character: `"1.9" > "1.26"` is `True`, because `"9" > "2"`. Compared that way, NumPy 1.9 (from 2014) would pass a `>= 1.26` check and NumPy 1.26 would fail a `>= 1.9` one. The correct comparison converts each component to an integer and compares the tuples: `(1, 9) < (1, 26)` because 9 < 26 at the second position. Tuples of unequal length are padded with zeros, so `(1, 26)` equals `(1, 26, 0)`. Pre-release tags such as `2.0.0rc1` ("release candidate 1") are dropped for our purposes: we treat `2.0.0rc1` as `(2, 0, 0)`. This is a simplification (a release candidate is strictly *older* than the final release), and it is a safe one here because the course never pins to an exact release.

**Worked example.** Is installed `"1.26.4"` at least minimum `"1.9"`? Parse: (1, 26, 4) and (1, 9). Pad: (1, 26, 4) and (1, 9, 0). Compare position by position: 1 = 1, then 26 > 9, so yes. Is `"3.10.14"` at least `"3.11"`? (3, 10, 14) against (3, 11, 0): 3 = 3, then 10 < 11, so no. Python 3.10 fails the course's minimum because the code uses `X | Y` type hints and `match` statements that 3.10 lacks or handles differently.

### 2.3 Why engineering code never tests `a == b` on floats

Type `0.1 + 0.2` into Python. You get `0.30000000000000004`. Nothing is broken: a computer stores numbers in binary, and 0.1 in binary is a repeating fraction (like 1/3 in decimal), so it is rounded to the nearest of the 2⁵³ values a 64-bit float can hold. Adding two rounded numbers gives a result rounded again. The relative error is about 10⁻¹⁶, called **machine epsilon** for `float64` (2.2 × 10⁻¹⁶). Harmless for an atmosphere table; fatal for an `if p == p_expected:` test, which will fail for no physical reason.

The rule for the whole course: compare floating-point results with a **tolerance**, `abs(a - b) <= tol`, where `tol` is chosen from the physics and the arithmetic, not from habit. `P00.07` makes this precise with `pytest.approx`. In this lesson the smoke test checks the ISA sea-level speed of sound, a = √(γ R T₀) = √(1.4 × 287.05287 × 288.15) = 340.294 m/s, against 340.29 with a tolerance of 0.05 m/s: loose enough that no correct machine fails, tight enough that a wrong constant or a wrong `sqrt` fails.

### 2.4 What a test runner does

`pytest` is a program that finds files named `test_*.py`, runs every function in them whose name starts with `test_`, and reports which ones raised an exception. A test is nothing more than Python code with `assert` statements: `assert meets_minimum("1.26.4", "1.26")` passes silently or fails loudly. There is no magic, which is the point. `@pytest-getstarted` "Get Started" shows the two-line first test and the `pytest --version` check we use as the environment's own smoke test. In this repository, the tests for each lesson import your code through `conftest.py` at the repo root, so you always run `pytest` from there:

```text
python -m pytest curriculum/phase-00/lesson-01-dev-environment
```

A green run means every assertion held on your machine. A red run prints the assertion that failed and the values on each side, which is the most useful error message you will see all course. `P00.07` is entirely about writing good tests; for now you only need to run them.

### 2.5 What the check proves and what it does not

The check proves that a particular interpreter can import a particular set of packages of acceptable versions, and that one computation gives the known answer. It does *not* prove that two machines will agree to the last digit (compiled libraries differ), that the packages are uncorrupted, or that your editor is running the same interpreter as your terminal (the most common failure of all, see §4). A passing check is necessary, not sufficient. Treat it like a pre-flight walk-around: it catches the obvious, and you still fly carefully.

## 3. Build the computational model

There is no physics model in this lesson; the "model" is the environment itself and the script that inspects it. Each step is a command to run or a code block to complete in `starter/main.py`.

### Step 1: Install Python 3.11 or newer

Open a terminal and type `python --version` (on some systems `python3 --version`). If it prints 3.11 or higher, skip to Step 2. Otherwise install the current release from python.org (Windows: tick "Add python.exe to PATH" in the installer; macOS: the python.org installer or `brew install python`; Linux: your package manager, for example `sudo apt install python3 python3-venv`). Close and reopen the terminal and check again. The `@numpy-userguide` "Installation" page describes the same prerequisites from NumPy's point of view.

### Step 2: Create and activate a virtual environment

From the folder where you cloned (or will clone, in `P00.02`) this repository:

```text
python -m venv .venv                 # creates the folder .venv with a private Python
source .venv/bin/activate            # macOS / Linux
.venv\Scripts\activate               # Windows (PowerShell or Command Prompt)
```

Your prompt now starts with `(.venv)`. Type `python -c "import sys; print(sys.executable)"` and confirm the path points inside `.venv`. That one line is the definitive test of "which Python am I running?"; use it whenever something imports on the command line but not in your editor.

If you prefer `uv`, a faster drop-in tool: `uv venv` creates the environment and `uv pip install -r requirements.txt` fills it. Either tool is fine; the result is the same folder structure.

### Step 3: Install the course packages

```text
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

`requirements.txt` at the repo root lists each package with its minimum version. Writing `python -m pip` rather than bare `pip` guarantees you install into the Python you just activated, not some other one on the path. When it finishes, `python -m pip list` shows what is installed. Check pytest separately as `@pytest-getstarted` "Install pytest" suggests: `python -m pytest --version`.

### Step 4: Editor

Install Visual Studio Code (free) and its Python extension. Open the repository folder (File → Open Folder). Press Ctrl+Shift+P (⌘⇧P on macOS), type "Python: Select Interpreter" and pick the one inside `.venv`. The status bar now shows that interpreter, and the editor's "Run" button, its integrated terminal and its test explorer all use it. Any editor works, but it must be pointed at the `.venv` interpreter or you will spend an afternoon on "it works in the terminal but not in the editor".

Jupyter notebooks are optional in this course and never the primary artifact (`LESSON_TEMPLATE.md`). VS Code opens `.ipynb` files directly once the Python extension is installed; you do not need a separate install to read a companion notebook.

### Step 5: Parse and compare versions

Open `starter/main.py`. Implement `parse_version`:

```python
import re

def parse_version(text: str) -> tuple[int, ...]:
    m = re.match(r"\s*(\d+(?:\.\d+)*)", text)     # leading digits separated by dots
    if not m:
        raise ValueError(f"not a version string: {text!r}")
    return tuple(int(part) for part in m.group(1).split("."))
```

The regular expression grabs the longest run of `digits.digits.digits` at the start and ignores whatever follows (`rc1`, `.post1`, `+local`). Raising `ValueError` on garbage like `"latest"` is deliberate: a check script that treats a nonsense version as acceptable is worse than no check.

```python
def meets_minimum(installed: str, minimum: str) -> bool:
    a, b = parse_version(installed), parse_version(minimum)
    n = max(len(a), len(b))
    a += (0,) * (n - len(a))
    b += (0,) * (n - len(b))
    return a >= b
```

Python compares tuples element by element, left to right, exactly as §2.2 describes. The zero-padding makes `(1, 26)` and `(1, 26, 0)` equal.

### Step 6: Ask the environment what is installed

```python
import importlib.metadata

def installed_version(name: str) -> str | None:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return None
```

`importlib.metadata` reads the version recorded when the package was installed, without importing the package (importing Matplotlib takes a second and may need a display; reading its metadata takes a millisecond). A missing package gives `None`, not a crash: the report should list what is missing, not stop at the first problem.

### Step 7: The check itself, written to be testable

```python
def check_packages(requirements=REQUIREMENTS, lookup=installed_version):
    results = []
    for name, minimum in requirements.items():
        version = lookup(name)
        ok = version is not None and meets_minimum(version, minimum)
        results.append(CheckResult(name, version, minimum, ok))
    return results
```

Notice the `lookup` argument. By default it asks the real environment; a test can pass a dictionary's `.get` instead and pretend that SciPy is 1.9.3 or that pint is missing. This is **dependency injection**: the logic that decides OK or FAIL is separated from the code that touches the outside world, so the logic can be tested on every path, including failures you cannot easily produce on a healthy machine. You will use the same pattern in `P00.08` to test a flight-log parser on a string instead of a file. `check_python` works the same way, with the running interpreter's `sys.version_info` as the default.

### Step 8: Smoke test and report

```python
import numpy as np

def speed_of_sound_smoke() -> float:
    return float(np.sqrt(GAMMA * R_AIR * T0))     # 340.29 m/s at ISA sea level
```

One computation with a known answer, through the library you most depend on. `format_report` turns the list of results into a fixed-width table with an OK/FAIL column and the word `missing` where there is no version; `main()` prints it, prints the smoke test, and returns exit code 0 only when everything passes, so a script or a continuous-integration job can rely on it.

### Verify

```text
python curriculum/phase-00/lesson-01-dev-environment/starter/main.py
python -m pytest curriculum/phase-00/lesson-01-dev-environment
```

The tests check: `parse_version` on six shapes of version string plus one rejection; `meets_minimum` on the lexicographic trap (`"1.26.4" >= "1.9"`, `"10.0" >= "9.9"`), equality under padding, and transitivity on a sorted list; `check_packages` with an injected lookup covering present-and-new-enough, present-but-old, and missing; `check_python` with injected version tuples; and then the real thing: `test_real_environment_passes`, which fails if *your* machine is not ready. Finally the smoke test must give 340.29 ± 0.05 m/s and the report must mention every component with its status.

### Compare

Compare your report against two independent sources. `python -m pip list` is pip's own view of the environment: every version in your report must appear there. `python -c "import numpy; print(numpy.__version__)"` is the package's own view; it should agree with `importlib.metadata`. If `pip list` shows a package your report calls missing, you have two Pythons and are running the wrong one; `sys.executable` will show which.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Break the environment on purpose, three ways, and record the check's output each time. Write the prediction *before* running each one.

1. **Wrong interpreter.** Deactivate the environment (`deactivate`) and run `python main.py` with your system Python. Prediction: Python passes; every package is missing or old; exit code 1. (If your system Python happens to have all the packages, the check passes, which is itself a lesson: you cannot tell from a pass which Python you used. `sys.executable` can.)
2. **Missing package.** In a fresh throwaway environment (`python -m venv /tmp/scratch-env`, activate it), install only NumPy and run the check. Prediction: `numpy` OK, the rest `missing`, smoke test OK because it needs only NumPy.
3. **Too old.** In the same scratch environment, `python -m pip install "scipy<1.11"` if your platform offers such a wheel (if not, lower the minimum in a copy of `REQUIREMENTS` to something above what you have, which produces the same report path). Prediction: a FAIL line with a version shown, distinguishable from `missing`.

Then two measurements of your machine's arithmetic, which you will need in `P00.07`:

```python
import numpy as np
print(0.1 + 0.2 == 0.3, 0.1 + 0.2)
x32 = np.full(10_000_000, 0.1, dtype=np.float32).sum()
x64 = np.full(10_000_000, 0.1, dtype=np.float64).sum()
print(x32, x64, "exact: 1000000.0")
```

Prediction: the first line prints `False 0.30000000000000004`. The `float32` sum will be off from 1 000 000 by something visible in the seventh figure (NumPy's pairwise summation keeps it that small; a naive `float32` loop would be off in the second figure); the `float64` sum will be correct to at least fifteen figures, and on many machines exactly 1 000 000.0. Record your actual numbers. Raw data format: a plain-text table of (experiment, prediction, observed, diagnosis), kept with your lesson notes.

Delete the scratch environment afterwards; the point of isolation is that you can.

## 5. Validate and improve

Lay the expected and observed side by side:

| component | minimum | observed (author's machine) | status | explanation of any gap |
|---|---|---|---|---|
| python | 3.11 | 3.13.16 | OK | newer is fine: minor releases add, do not remove |
| numpy | 1.26 | 2.5.3 | OK | major version 2 changed some defaults (`P00.03` notes one); the course uses none that moved |
| scipy | 1.11 | 1.18.1 | OK | |
| matplotlib | 3.8 | 3.11.2 | OK | |
| sympy | 1.12 | 1.14.0 | OK | |
| pint | 0.23 | 0.26.1 | OK | |
| pytest | 8.0 | 9.1.1 | OK | |
| smoke test | 340.29 ± 0.05 m/s | 340.29 | OK | |

Every discrepancy you meet has one of a few causes, and "it just didn't work" is not one of them:

- **The report says `missing` but `pip list` shows the package.** Two interpreters. Diagnose with `sys.executable` in both places; fix by activating the environment or re-selecting the interpreter in the editor.
- **`FAIL` with a version shown.** Too old. `python -m pip install --upgrade <package>` inside the environment.
- **Python itself fails.** You installed 3.11+ but the terminal still finds the old one first. On Windows use `py -3.13 -m venv .venv` to pick the interpreter explicitly; on macOS/Linux use `python3.13 -m venv .venv`.
- **The smoke test fails but everything else passes.** Almost always a typo in a constant in your own code; genuinely broken NumPy builds are rare and fail at import, not at `sqrt`.

Improvement worth making: have `main()` also print `sys.executable`, so the report itself says which Python produced it. Improvement worth not making: trying to compare versions with a full specifier grammar (`~=`, `!=`, wildcards). The course does not need it, and the simple tuple ordering is correct for everything in `requirements.txt`.

## 6. Defend your work

Submit: the report your `main.py` prints; the pytest output; the three failure transcripts, each with a one-line diagnosis naming the cause from §5; the arithmetic table from §4; and three sentences stating what a passing check does and does not prove about your machine.

### Adversarial self-checks

- **What if a version string is `"2.0.0rc1"`? `"1.26.4+mkl"`? `""`?** The first two parse to (2, 0, 0) and (1, 26, 4); the third must raise, not pass. Try them.
- **Which assumption fails first in the real world?** "Same version means same behaviour." Two NumPy 2.5.3 installs with different BLAS backends can differ in the last few bits of a matrix product. That is why the course's tests use tolerances (`P00.07`) and why a reviewer re-running your code expects agreement to 10⁻¹², not bit-for-bit.
- **What measurement could falsify the check?** Run the same `main.py` on a second machine that reports OK everywhere and compare the printed speed of sound to all 17 digits. If they differ, the check was right (both environments are fit) and the lesson is that fitness is not identity.
- **Does the result survive a change in working directory? In shell?** Run from the repo root, from inside the lesson folder, from PowerShell and from Command Prompt. Same report, or you have a path dependence to remove.
- **Which decision changes if the tolerance on the smoke test doubled?** None here; 0.1 m/s still cannot be produced by a healthy machine. But if it were 10 m/s, a `GAMMA = 1.3` typo would pass. Tolerances must be tight enough to catch the mistakes you can plausibly make.

## Check yourself

1. What single Python expression tells you which interpreter is actually running your code?
2. Is `"1.9" >= "1.26"` true or false in Python, and why is that the wrong answer for versions?
3. Why does `check_packages` take a `lookup` argument instead of calling `importlib.metadata.version` directly?
4. What does `0.1 + 0.2 == 0.3` evaluate to, and what is the order of magnitude of machine epsilon for `float64`?
5. Name one thing a passing environment check does not prove.
6. Why does the lesson say `python -m pip install` rather than `pip install`?
7. The smoke test tolerance is 0.05 m/s on 340.29 m/s. If `R_AIR` were mistyped as 287.5, would the test catch it? (Estimate, do not compute exactly.)
