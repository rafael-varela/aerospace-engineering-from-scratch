# P00.07  Testing Engineering Code

> A calculation you have not tested is a guess with a decimal point.

| | |
|---|---|
| **Learning outcomes** | 1. Explain, with a worked binary example, why two floating-point results that "should" be equal usually are not, and choose a relative or absolute tolerance for a comparison. 2. Implement a small lift-and-drag coefficient calculator that validates its inputs and raises `ValueError` on unphysical values. 3. Write a pytest suite that uses `pytest.approx`, `numpy.testing.assert_allclose`, `pytest.raises`, `pytest.mark.parametrize`, limiting cases and a regression test against a stored reference table. 4. Decide, for a given failing test, whether the code, the test or the stored reference is wrong. |
| **Prerequisites** | `P00.01`, `P00.03`; skill check: can write a Python function with a docstring and run it from the command line; can rearrange $L = \tfrac{1}{2}\rho V^2 S C_L$ for any one variable. |
| **Engineering question** | Can we trust the lift-to-drag ratio our drag-polar code reports for a light aircraft, to within 0.01 %, after every future edit to that code? |
| **Assumptions and limits** | Parabolic drag polar $C_D = C_{D0} + k C_L^2$ with constant $C_{D0}$ and $k$; incompressible flow (Mach below about 0.3); no stall, so $C_L$ is unbounded in the model. The floating-point discussion assumes IEEE 754 double precision (NumPy `float64`). |
| **Analytical work** | Derive the lift coefficient from the lift equation; derive $C_L^*$ and $(L/D)_{\max}$ from the parabolic polar; bound the rounding error of a repeated subtraction. |
| **Computational deliverable** | `starter/main.py`: `dynamic_pressure`, `lift_coefficient`, `drag_coefficient`, `drag_polar`, `lift_to_drag`, `max_lift_to_drag`, `polar_table`, `is_close`, `rocket_mass_closed_form`, `rocket_mass_stepwise`. Tests prove: analytic values to 1e-12 relative, limiting cases ($V \to 0$, $k \to 0$, $C_L \to 0$), input validation, the optimum beats every sampled point, and agreement with a stored table to 2e-7. |
| **Practical activity** | **computational experiment**: measure how rounding error in a stepwise rocket-mass loop grows with the number of steps and with the step size, and compare with the bound from part 2. |
| **Equipment** | None |
| **Safety** | None beyond normal computer use |
| **Evidence** | Your completed `main.py`, the passing test run (`pytest -v` output), the drift-versus-steps table and plot from part 4, and a half-page note on which tolerance you chose for each test and why. |
| **Assessment** | `assessments/phase-00/grading-rubric.md`; competencies 2 (analytical), 3 (computational) and 4 (interpretation); pass when all public tests pass *and* the tolerance note defends every tolerance with a number. |
| **References** | `@pytest-getstarted` "Get Started": "Create your first test", "Assert that a certain exception is raised", "Compare floating-point values with pytest.approx", "Group multiple tests in a class" (required-free); `@numpy-userguide` NumPy reference, `numpy.testing` (`assert_allclose`) (supplementary-free). |
| **Extensions** | ★ Add a `hypothesis`-style property test by hand: generate 1000 random $(C_{D0}, k)$ pairs with a seeded `numpy.random.default_rng` and check the optimum property for each. ★ Replace the parabolic polar with a tabulated polar read from a file and make the regression test cover interpolation. |
| **Time** | unmeasured |

## 1. The engineering problem

You have written a drag-polar calculator. It takes a lift coefficient and returns a drag coefficient, and from those two numbers a lift-to-drag ratio $L/D$. Someone on your team will use $L/D$ to size the fuel tanks of a light aircraft: a 1 % error in $L/D$ is a 1 % error in range. Next month a teammate will "tidy up" the code, and the month after that you will add compressibility corrections. How do you know the number is still right?

The answer in professional engineering software is the same as the answer in a wind tunnel: you do not trust a reading, you check it against a reference, and you repeat the check every time something changes. In code the reference is a *test*: a short function that runs your code on a case whose answer you know by other means and fails loudly if the two disagree. A **test suite** is a set of such tests run together by a tool called **pytest**. The success criterion for this lesson is a suite that passes against a correct calculator, fails against a broken one, and tells you *which* property broke.

What goes wrong without it is not hypothetical. Untested unit conversions lost the Mars Climate Orbiter in 1999 (you met that story in P00.05); untested sign conventions regularly invert control surfaces in simulators; and a test that compares two floating-point numbers with `==` will pass on one machine and fail on another for reasons that have nothing to do with physics. The centre of this lesson is that last trap: **why exact float equality fails and what to assert instead.**

## 2. Understand and derive

### 2.1 The calculator we are testing

Define a term each time it appears. The **dynamic pressure** $q$ is the kinetic energy per unit volume of the oncoming air:

$$q = \tfrac{1}{2}\rho V^2 \qquad [\mathrm{Pa}],$$

with air density $\rho$ in kg/m³ and true airspeed $V$ in m/s. Units check: $\mathrm{kg\,m^{-3}} \cdot \mathrm{m^2\,s^{-2}} = \mathrm{kg\,m^{-1}\,s^{-2}} = \mathrm{N\,m^{-2}} = \mathrm{Pa}$. Lift $L$ and drag $D$ are forces in newtons; dividing them by $q$ and the reference wing area $S$ (m²) gives the dimensionless **lift coefficient** and **drag coefficient**:

$$C_L = \frac{L}{qS}, \qquad C_D = \frac{D}{qS}.$$

For a whole aircraft at low speed, $C_D$ depends on $C_L$ through the **parabolic drag polar**

$$C_D = C_{D0} + k\,C_L^2, \qquad k = \frac{1}{\pi e\,AR},$$

where $C_{D0}$ is the drag at zero lift (skin friction and form drag), $AR$ is the wing aspect ratio and $e$ the Oswald efficiency factor. The model assumes the induced drag is exactly quadratic in $C_L$ and that $C_{D0}$ does not change with speed; it fails near stall (where $C_L$ stops growing) and above Mach 0.3 (where compressibility raises $C_{D0}$). We use $C_{D0} = 0.020$, $k = 0.050$ (roughly $e = 0.8$, $AR = 8$) as the lesson's default.

**Best lift-to-drag ratio.** Write $f(C_L) = L/D = C_L / (C_{D0} + k C_L^2)$. Differentiate with the quotient rule and set the result to zero:

$$f'(C_L) = \frac{(C_{D0} + k C_L^2) - C_L (2 k C_L)}{(C_{D0} + k C_L^2)^2} = \frac{C_{D0} - k C_L^2}{(C_{D0} + k C_L^2)^2} = 0 \;\Rightarrow\; C_L^* = \sqrt{C_{D0}/k}.$$

At the optimum the induced drag $kC_L^{*2}$ equals $C_{D0}$, so $C_D^* = 2C_{D0}$ and

$$(L/D)_{\max} = \frac{\sqrt{C_{D0}/k}}{2 C_{D0}} = \frac{1}{2\sqrt{C_{D0}\,k}}.$$

For the default polar: $C_L^* = \sqrt{0.4} = 0.6325$ and $(L/D)_{\max} = 1/(2\sqrt{0.001}) = 15.81$. Those two numbers are **analytic references**: a test can compare the code against them. Better still, a test can check the *property* "no sampled point on the polar beats the analytic optimum", which stays true for every $(C_{D0}, k)$ and does not depend on remembering 15.81.

### 2.2 Why `0.1 + 0.2 != 0.3`

A `float64` stores a number as $\pm\,1.f \times 2^{e}$ with 52 bits for the fraction $f$. Only fractions whose denominators are powers of two are exact. One tenth is $0.0001100110011\ldots_2$, a repeating binary fraction, so the computer stores the nearest representable value, about $0.1000000000000000055511$. The same happens to 0.2 and 0.3. Add the stored 0.1 and 0.2, round the sum to 52 bits, and you get $0.3000000000000000444$, while the stored 0.3 is $0.2999999999999999889$. They differ by one **ulp** (unit in the last place, the gap between neighbouring floats), which is $2^{-54} \approx 5.6 \times 10^{-17}$ near 0.3. So `0.1 + 0.2 == 0.3` is `False`, and no amount of care in your physics will change that.

The same thing happens in engineering code constantly. $\tfrac12 \cdot 1.225 \cdot 50^2$ is exactly 1531.25 in decimal, but 1.225 is not exact in binary (it is stored as $1.2250000000000000888\ldots$), so the product *may* land one ulp away from the literal `1531.25`. In this particular case the two roundings happen to cancel and Python returns exactly `1531.25`; change the density to 1.2 or the speed to 51 and they may not. You cannot tell in advance without doing the binary arithmetic, which is exactly why `==` is the wrong tool. A loop that subtracts a propellant mass of 0.1 kg per step a hundred times does not land on exactly $m_0 - 10$ kg: it lands on `39.99999999999986` for $m_0 = 50$.

**How large is the drift?** Each floating-point operation rounds its exact result to the nearest float, which introduces a relative error of at most $u = 2^{-53} \approx 1.1 \times 10^{-16}$ (half an ulp; `np.finfo(float).eps` is $2u$). After $n$ subtractions on a quantity of size about $m_0$, the accumulated error is bounded by roughly $n\,u\,m_0$, and in practice is smaller because errors partly cancel. For $m_0 = 50$ kg and $n = 100$: bound $\approx 100 \cdot 1.1 \times 10^{-16} \cdot 50 = 5.6 \times 10^{-13}$ kg. We will see about $1.4 \times 10^{-13}$ kg. Physically meaningless, numerically fatal to `==`.

### 2.3 What to assert instead: tolerances

Two numbers $a$ and $b$ are *close* if

$$|a - b| \le \max\big(\text{rel\_tol} \cdot \max(|a|, |b|),\ \text{abs\_tol}\big).$$

- The **relative tolerance** says "agree to this many significant figures". `rel_tol = 1e-12` means agreement to 12 figures, appropriate when both numbers come from a few arithmetic operations on exact inputs.
- The **absolute tolerance** is needed when the true answer is zero or nearly zero, because a relative tolerance on zero is zero: `is_close(0.0, 1e-12)` with any relative tolerance alone is `False`.

`pytest.approx(x, rel=..., abs=...)` implements this for scalars and small sequences (default `rel=1e-6`, `abs=1e-12`). `numpy.testing.assert_allclose(actual, desired, rtol=..., atol=...)` does it element-wise for arrays and prints the worst mismatch when it fails (its rule is `|a - d| <= atol + rtol * |d|`, a slightly different but equivalent-in-spirit form; note it is relative to `desired`, so put the reference second). Choosing the tolerance is an engineering decision, and the rule is: **as tight as the arithmetic allows, and justified by a number.** For an exact formula evaluated once, `1e-12` relative. For a value stored to 8 significant figures, `2e-7` relative. For a value that depends on a numerical method, the method's error estimate.

### 2.4 The five kinds of test in this lesson

1. **Analytic value.** Compare with a hand calculation: $q(1.225, 50) = 1531.25$ Pa.
2. **Limiting case.** $q \to 0$ as $V \to 0$; $C_D \to C_{D0}$ as $C_L \to 0$; $C_D$ independent of $C_L$ when $k = 0$. Limiting cases catch sign and factor errors that an analytic value at one point can miss.
3. **Invalid input raises.** A negative density or a zero wing area is a bug upstream; the function must refuse with `ValueError`, not return `inf` or `nan`. Silence is how bad numbers travel.
4. **Parametrized test.** The same property checked for several polars, so the test does not accidentally pass only for the default numbers.
5. **Regression test.** A table of outputs generated once by a reviewed version of the code and stored in the test file. It does not prove the code is right; it proves the code has not *changed*. When it fails, you must decide whether the change was intended.

## 3. Build the computational model

Open `starter/main.py`. Every function has a docstring that says what to return and when to raise. Work top to bottom; run `pytest curriculum/phase-00/lesson-07-testing-engineering-code -q` after each function and watch failures turn into passes.

### Step 1: input checks and dynamic pressure

Write two small helpers and use them everywhere. They accept scalars or arrays.

```python
import numpy as np

def _check_positive(name: str, value) -> None:
    arr = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must be finite, got {value!r}")
    if np.any(arr <= 0.0):
        raise ValueError(f"{name} must be strictly positive, got {value!r}")

def dynamic_pressure(rho, v):
    _check_positive("rho", rho)
    _check_nonnegative("v", v)      # same pattern, allows zero
    return 0.5 * np.asarray(rho, dtype=float) * np.asarray(v, dtype=float) ** 2
```

Two design choices to notice. `np.asarray(..., dtype=float)` makes the function work for a single float and for a whole array of speeds (vectorization, from P00.03). And the check for *finite* comes first: `nan <= 0.0` is `False`, so a NaN would otherwise slip through a positivity test.

### Step 2: coefficients and the polar

`lift_coefficient` divides by $qS$, so it must insist on strictly positive $V$ and $S$ (zero would divide by zero). `drag_coefficient` additionally rejects negative drag. `drag_polar` rejects negative $C_{D0}$ or $k$ and a non-finite $C_L$. `max_lift_to_drag` is the two formulas from part 2; because it divides by $\sqrt{C_{D0} k}$ it requires both strictly positive. `polar_table` stacks three columns with `np.column_stack`.

### Step 3: your own `is_close`

Implement the rule in 2.3 yourself before using the library versions; the point is to see that there is nothing magic inside `pytest.approx`.

```python
import math

def is_close(a, b, rel_tol=1e-9, abs_tol=0.0):
    if rel_tol < 0.0 or abs_tol < 0.0:
        raise ValueError("tolerances must be non-negative")
    if math.isnan(a) or math.isnan(b):
        return False
    if a == b:                       # catches equal infinities
        return True
    if math.isinf(a) or math.isinf(b):
        return False
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)
```

The test compares your function against `math.isclose` on a grid of cases, including the near-zero case where only an absolute tolerance can pass.

### Step 4: the rocket-mass drift demonstration

`rocket_mass_stepwise` subtracts `mdot * dt` in a loop and records every intermediate mass; `rocket_mass_closed_form` returns $m_0 - \dot m t$ directly. Both raise when the propellant would run out. The loop version is how a time-marching simulation depletes mass in Phase 05, and the test shows that it drifts from the closed form by a few $10^{-13}$ kg while remaining within tolerance.

### Step 5: read the test file as a worked example

`tests/test_main.py` is part of the lesson. Read it top to bottom and find each of the five kinds of test from 2.4. Note the comment above `REFERENCE_POLAR`: the table was produced once and its tolerance (`rtol=2e-7`) is justified by the 8 significant figures it was stored with.

### Verify

- **Analytic case:** `test_dynamic_pressure_analytic_value` and the optimum property in `test_max_lift_to_drag_is_the_maximum` (`k C_L^{*2} = C_{D0}` to 1e-12).
- **Limiting cases:** `test_dynamic_pressure_limits_and_scaling` ($V\to0$, $q \propto V^2$, $q \propto \rho$), `test_drag_polar_limiting_cases` ($C_L \to 0$, $k \to 0$, $C_{D0} \to 0$, evenness in $C_L$).
- **Conservation-like property:** `test_coefficients_round_trip` recovers the input force from the coefficient, and `test_rocket_mass_stepwise_drifts_but_stays_within_tolerance` checks that the drift is bounded by $n\,\varepsilon\, m_0$ and that mass never increases.
- **Error plot:** part 4 below.

### Compare

Run the same comparison with the library tools: `math.isclose`, `pytest.approx` and `numpy.testing.assert_allclose` all implement the tolerance rule of 2.3, with different defaults. `math.isclose` defaults to `rel_tol=1e-9, abs_tol=0`; `pytest.approx` to `rel=1e-6, abs=1e-12`; `assert_allclose` to `rtol=1e-7, atol=0`. The defaults differ by three orders of magnitude, which is one more reason to state your tolerance explicitly in every assertion.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Measure how the rounding drift of the stepwise mass depletion grows, and compare with the bound from 2.2.

1. Keep $m_0 = 50$ kg and the total burn $\dot m\,t_{\text{end}} = 10$ kg fixed.
2. For $n \in \{10, 100, 1000, 10^4, 10^5, 10^6\}$ steps, set $\mathrm{d}t = t_{\text{end}}/n$ and compute `drift = abs(rocket_mass_stepwise(...)[-1] - rocket_mass_closed_form(m0, mdot, t_end))`.
3. Record `n`, `drift`, and the bound `n * np.finfo(float).eps / 2 * m0` in a table. Write down your prediction first: does the drift grow like $n$, like $\sqrt{n}$, or not at all?
4. Plot drift and bound against $n$ on log-log axes (`ax.loglog`), save to `drift_vs_steps.png`.
5. Repeat with `mdot * dt` chosen so that the step is exactly representable in binary (for example 0.125 kg). What happens to the drift and why?

Raw data format: a CSV with columns `n, dt_s, drift_kg, bound_kg`, one row per $n$.

## 5. Validate and improve

Put the three things side by side: the analytic bound $n\,u\,m_0$, the measured drift, and the tolerance you chose in the test (`rtol=1e-12`). Before you ran it you may have predicted $\sqrt{n}$ growth from "random" rounding errors cancelling. The measurement says otherwise: with $m_0 = 50$ kg the drift is $1.4 \times 10^{-13}$ kg at $n = 100$, $2.3 \times 10^{-11}$ at $10^4$ and $3.2 \times 10^{-9}$ at $10^6$, a log-log slope close to 1, and only a factor of two to five below the bound. The reason is that the subtracted step `mdot * dt` carries the *same* representation error every time it is used, so the errors have the same sign and add; they are systematic, not random. "Rounding errors cancel" is an assumption, and here it fails. Write that discrepancy down and explain it; "numerical noise" is not an explanation.

The consequence for the test: at $n = 10^6$ the drift of $3.2 \times 10^{-9}$ kg on 40 kg is $8 \times 10^{-11}$ relative, so a `rtol=1e-12` assertion would fail even though the code is correct. That is the lesson's central judgment: **a tolerance is a claim about the arithmetic, and it must be re-derived when the arithmetic changes.** The right fix is to loosen the tolerance *with a reason written in the test* (for example `rtol = 10 * n * eps`), not to switch to `==`, and not to loosen it to 1e-3 "to be safe", which would hide a real factor-of-two bug in the next edit.

Explain the step-0.125 result in your note: a step that is a dyadic rational is represented exactly, and so is every partial sum while the mass stays within the same binade, so the drift is exactly zero. This is why simulation codes often choose time steps like $2^{-7}$ s.

## 6. Defend your work

Submit: `main.py`, the `pytest -v` log, `drift_vs_steps.png` with its CSV, and your tolerance note. Use the report template in `projects/_template/report-template.md` for the note.

### Adversarial self-checks

- What if the lift force doubles? $C_L$ doubles, $C_D$ does not: the round-trip test must still pass, and `lift_to_drag` must change. Does your suite notice if someone accidentally squares the lift?
- What if $V \to 0$? `dynamic_pressure` returns exactly 0.0 (no rounding: $0^2 = 0$ exactly) but `lift_coefficient` must raise. Why is it right for the two functions to behave differently?
- Which assumption fails first in the real world? The parabolic polar near stall. Which test would still pass with a stalled wing, and what new test would catch it?
- Does the result survive a change in solver or input? Replace `0.5 * rho * v ** 2` by `rho * v * v / 2`: which tests pass, which fail, and is a failure here a bug?
- Which decision changes if the uncertainty doubles? If your stored reference table had been rounded to 4 figures instead of 8, what tolerance would the regression test need, and could it still detect a 0.01 % change in $L/D$?

## Check yourself

1. In your own words, why is `0.1 + 0.2 == 0.3` false, and what is the size of the difference compared with 0.3?
2. Which tolerance would you use to compare $q = \tfrac12 \rho V^2$ with a hand value, and why not `1e-3`?
3. `is_close(0.0, 1e-15)` is `False` with the default tolerances. What single argument makes it `True`, and when is that the *right* thing to do?
4. A regression test fails after you corrected a genuine bug in `drag_polar`. What do you do with `REFERENCE_POLAR`, and what do you write in the commit message?
5. State the limiting case of $C_D$ as $k \to 0$ and the physical meaning of that limit.
6. Why does `lift_coefficient` reject $V = 0$ while `dynamic_pressure` accepts it?
7. For $C_{D0} = 0.03$ and $k = 0.06$, compute $C_L^*$ and $(L/D)_{\max}$ by hand, then check with `max_lift_to_drag`.
8. A colleague says "the parametrized test passes for four polars, so the formula is proven". What is wrong with that statement?
