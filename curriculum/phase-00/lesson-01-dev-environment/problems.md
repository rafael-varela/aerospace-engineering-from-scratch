# P00.01  Problem set

Attempt every problem before opening `solutions/`. Record what you tried, including the exact commands and their output.

## 1. (Warm-up) Versions by hand

Without running code, write the tuple `parse_version` should return for each of: `"3.11.4"`, `"1.26"`, `"2.0.0rc1"`, `"9.1.1"`, `"0.23.post1"`. Then decide, by the padded tuple rule, whether each of the following is true: `"1.26.4" >= "1.26"`, `"1.9.0" >= "1.26"`, `"3.11" >= "3.11.0"`, `"10.0" >= "9.9"`. Answer form: five tuples and four true/false answers with the deciding position marked.

## 2. (Argument) The case for isolation

A colleague says: "Virtual environments are bureaucracy. I install everything globally and it works." Write a short argument (under 200 words) for one environment per project. Use at least one concrete scenario involving two projects with incompatible requirements, and one involving a result that must be reproduced a year later. Then state the strongest point *for* the colleague's position and say when it would actually be fine.

## 3. (Unfamiliar scenario) Diagnose from a transcript

A lab PC produces this when a student runs the environment check from inside VS Code's "Run" button, while the same script passes in the terminal:

```text
component      installed    minimum  status
python           3.12.3       3.11  OK
numpy            missing       1.26  FAIL
scipy            missing       1.11  FAIL
matplotlib       missing        3.8  FAIL
sympy            missing       1.12  FAIL
pint             missing       0.23  FAIL
pytest           missing        8.0  FAIL
Traceback (most recent call last):
  File "main.py", line 19, in <module>
    import numpy as np
ModuleNotFoundError: No module named 'numpy'
```

Explain what has happened, name the single command or setting that fixes it, and explain why the traceback appears *after* the table rather than before (hint: look at where `numpy` is first used in `main.py`). Then explain why the script as written could never have printed this exact output, and what change to `main.py` would make such an output possible. Answer form: a diagnosis paragraph and the fix.

## 4. (Numerical, uses your code) Extend the check

Add `pyyaml` (minimum `6.0`; its distribution name is `PyYAML`, its import name is `yaml`) and `requests` (minimum `2.31`) to a copy of `REQUIREMENTS` and run `check_packages` on it. Report the output. Then explain, using what `importlib.metadata.version` actually looks up, which of the two names (distribution name or import name) you must give it, and show what happens with the other one. Answer form: the report lines and a two-sentence explanation.

## 5. (Numerical) Floating point on your machine

Run the two-line arithmetic experiment from part 4 of the lesson. Then sum `0.1` ten million times in a plain Python `for` loop (not NumPy) and compare with `np.float64` sum and with the exact answer 1 000 000. Report all three results to full precision and the relative error of each. Explain why the NumPy sum is more accurate than the loop even though both use 64-bit floats. Answer form: a three-row table and a short explanation.

## 6. (Judgment) False confidence

Describe two situations in which the environment check prints OK on every line and the course's numerical results would still differ between your machine and a reviewer's. For each, say how large the difference could be (order of magnitude, with reasoning) and whether it matters for an atmosphere table, for a trajectory integrated over an hour of flight, and for a structural safety margin of 1.5. What would you add to the check, or to the way results are reported, to make the difference visible?
