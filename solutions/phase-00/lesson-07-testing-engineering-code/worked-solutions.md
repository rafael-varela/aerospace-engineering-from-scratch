# P00.07  Worked solutions

Each solution states: assumptions, units, the derivation, the result (with the tolerance that would be appropriate in a test), at least one independent check, and the common mistakes. Problems with a seed $s$ are worked for $s = 15$ and the general formula is given.

## 1. Coefficients in level flight

**Assumptions.** Level, unaccelerated flight so $L = W$; sea-level ISA density; parabolic polar with $C_{D0} = 0.020$, $k = 0.050$; incompressible flow at 40 m/s (Mach 0.12).

**Units.** $m$ in kg, $g$ in m/s², $W$ in N, $\rho$ in kg/m³, $V$ in m/s, $S$ in m², $q$ in Pa; coefficients dimensionless.

**Derivation.** For $s = 15$: $m = 750$ kg, $W = 750 \times 9.80665 = 7354.99$ N.
(a) $q = \tfrac12 \times 1.225 \times 40^2 = 980.0$ Pa.
(b) $C_L = W/(qS) = 7354.99/(980.0 \times 12.0) = 0.6254$. In general $C_L = (600+10s)\times 9.80665/11760 = (600 + 10 s) \times 8.339 \times 10^{-4}$.
(c) $C_D = 0.020 + 0.050 \times 0.6254^2 = 0.03956$; $L/D = 0.6254/0.03956 = 15.81$.
(d) $D = C_D\,q\,S = 0.03956 \times 980 \times 12 = 465.2$ N.

**Result.** $q = 980.0$ Pa, $C_L = 0.6254$, $C_D = 0.03956$, $L/D = 15.81$, $D = 465.2$ N. In a test: `lift_coefficient(7354.9875, 1.225, 40.0, 12.0) == pytest.approx(0.625424107, rel=1e-9)`.

**Check.** $D = W/(L/D) = 7354.99/15.81 = 465.2$ N, agreeing with (d). Also, $s = 15$ happens to put $C_L$ within 1 % of $C_L^* = 0.6325$, so $L/D$ must be within a hair of $(L/D)_{\max} = 15.81$: it is (15.8104 versus 15.8114).

**Common mistakes.** Using mass instead of weight for lift (a factor 9.8); forgetting the $\tfrac12$ in $q$ (a factor 2 in $C_L$ and about a factor 1.6 in $L/D$); reporting $C_L$ with units.

## 2. Two optima on the polar

**Assumptions.** $C_D = C_{D0} + kC_L^2$ for all $C_L > 0$; $C_{D0}, k > 0$ constants.

**Derivation.** (a) $f(C_L) = C_L/(C_{D0} + kC_L^2)$. Quotient rule: $f' = [(C_{D0} + kC_L^2) - C_L \cdot 2kC_L]/(C_{D0} + kC_L^2)^2 = (C_{D0} - kC_L^2)/(\cdot)^2$. The denominator is positive, so $f' = 0 \iff kC_L^2 = C_{D0}$, i.e. $C_L^* = \sqrt{C_{D0}/k}$, and the induced drag $kC_L^{*2}$ equals $C_{D0}$, giving $C_D^* = 2C_{D0}$. Then $(L/D)_{\max} = C_L^*/C_D^* = \sqrt{C_{D0}/k}/(2C_{D0}) = 1/(2\sqrt{C_{D0}k})$. Second derivative, or the sign change of $f'$ from positive to negative, confirms a maximum.

(b) At fixed weight, $V = \sqrt{2W/(\rho S C_L)} \propto C_L^{-1/2}$, and power $P = DV = (C_D/C_L)\,W\,V \propto C_D\,C_L^{-3/2}$. Let $h(C_L) = C_{D0}C_L^{-3/2} + kC_L^{1/2}$. Then $h' = -\tfrac32 C_{D0}C_L^{-5/2} + \tfrac12 kC_L^{-1/2} = 0 \Rightarrow kC_L^2 = 3C_{D0} \Rightarrow C_{L,\text{min power}} = \sqrt{3C_{D0}/k} = \sqrt3\,C_L^*$.

(c) The assumption is that the polar holds for every $C_L$, so $C_L$ can grow without limit. Near stall $C_L$ saturates at $C_{L,\max}$ and $C_D$ rises faster than quadratically. If $C_L^* > C_{L,\max}$ (a very low-drag, low-aspect-ratio case), the analytic optimum lies at an unreachable lift coefficient; a test comparing it with the sampled maximum on the *real* polar would fail, and the failure would be a true statement about the model's limits, not a bug.

**Result.** $C_L^* = \sqrt{C_{D0}/k}$, $(L/D)_{\max} = 1/(2\sqrt{C_{D0}k})$, $C_{L,\text{min power}} = \sqrt{3}\,C_L^*$. For the default polar: 0.6325, 15.81, 1.095.

**Check.** Dimensional: all quantities dimensionless. Limiting: $k \to 0$ gives $C_L^* \to \infty$ and $(L/D)_{\max} \to \infty$, as expected for a wing with no induced drag. Numerical: `lift_to_drag(np.linspace(0.05, 2, 400))` peaks at 15.811, within $10^{-4}$ relative of the formula.

**Common mistakes.** Differentiating $C_D$ instead of $C_L/C_D$; dropping the 2 in $(L/D)_{\max}$ (this is mutation 3 in problem 7); in (b) forgetting that $V$ depends on $C_L$.

## 3. Floating point by hand

**Derivation.** (a) $0.1 \times 2 = 0.2 \to 0$; $0.4 \to 0$; $0.8 \to 0$; $1.6 \to 1$, keep $0.6$; $1.2 \to 1$, keep $0.2$; and $0.2$ already appeared, so the pattern repeats: $0.1 = 0.0\overline{0011}_2$.

(b) $2^{10} = 1024 \le 1531.25 < 2048$, so ulp $= 2^{10-52} = 2^{-42} = 2.27 \times 10^{-13}$. $2^{-2} = 0.25 \le 0.3 < 0.5$, so ulp $= 2^{-54} = 5.55 \times 10^{-17}$. `math.ulp(1531.25)` and `math.ulp(0.3)` return exactly these.

(c) The stored 1.225 differs from the true value by at most half an ulp of 1.225, i.e. $2^{-53} \approx 1.1 \times 10^{-16}$ (the actual error is $8.9 \times 10^{-17}$). Multiplying by $\tfrac12 \times 2500 = 1250$ scales it to at most $1.4 \times 10^{-13}$ Pa, about 0.6 ulp of the result; the final rounding adds up to 0.5 ulp. Bound: about 1.1 ulp, so the result is within 2 ulps $= 4.5 \times 10^{-13}$ Pa of 1531.25. Running `0.5 * 1.225 * 50.0 ** 2 - 1531.25` gives exactly `0.0`: the two roundings happened to cancel. The bound was honest; the actual error was smaller.

(d) `rel=1e-12` allows $1.5 \times 10^{-9}$ Pa, three thousand times the bound: safe. `rel=1e-15` allows $1.5 \times 10^{-12}$ Pa, only about three times the bound: it would pass here but with no margin, and a slightly different input could legitimately fail it. A tolerance should sit at least an order of magnitude above the arithmetic bound; `1e-14` is the tightest defensible choice.

**Common mistakes.** Taking the ulp relative to 1 instead of relative to the number; believing the actual error is always the bound; concluding from the `0.0` that `==` would be safe.

## 4. Wind-tunnel balance in imperial units

**Assumptions.** The function is meant to return $C = F/(qS)$ with consistent imperial units: lbf, lbf/ft², ft².

**The three tests.**
1. *Known case.* `coefficient_from_balance(2.5, 1.0, 1.0) == pytest.approx(2.5, rel=1e-12)`: with unit $q$ and unit area the coefficient equals the force numerically, and the arithmetic is one division of exactly representable numbers, so 1e-12 is generous.
2. *Limiting / scaling case.* `coefficient_from_balance(0.0, 5.0, 2.0) == 0.0` exactly ($0/x$ is exact), and `coefficient_from_balance(3.0, 2.0, 1.5) == pytest.approx(coefficient_from_balance(3.0, 4.0, 1.5) * 2, rel=1e-12)`: doubling $q$ halves $C$.
3. *Invalid input.* `with pytest.raises(ValueError): coefficient_from_balance(3.0, 0.0, 1.5)` and the same for negative area and NaN force.

**Which test catches the stray 32.174.** Only test 1. Test 2 checks the *shape* of the dependence (zero stays zero; ratios are preserved), and a constant factor survives both. Test 3 is about refusing inputs. This is why every suite needs at least one absolute, hand-computed value: property tests alone are blind to an overall scale factor. (Problem 7, mutation 2, shows the same thing in the lesson's own suite.)

**Common mistakes.** Choosing a "known case" whose answer you got from the function itself; writing only ratio tests and declaring the function verified.

## 5. Measuring the drift

**Assumptions.** $m_0 = 50 + s = 65$ kg for $s = 15$; total burn 10 kg; `mdot * dt = 10 / n`.

**Result** (one run; your values differ only through $m_0$):

| $n$ | step [kg] | drift [kg] | bound $n\,u\,m_0$ [kg] |
|---|---|---|---|
| 10 | 1.0 | 0 | $7.2 \times 10^{-14}$ |
| $10^2$ | 0.1 | $7.1 \times 10^{-14}$ | $7.2 \times 10^{-13}$ |
| $10^3$ | 0.01 | $1.3 \times 10^{-12}$ | $7.2 \times 10^{-12}$ |
| $10^4$ | 0.001 | $1.6 \times 10^{-11}$ | $7.2 \times 10^{-11}$ |
| $10^5$ | $10^{-4}$ | $3.3 \times 10^{-10}$ | $7.2 \times 10^{-10}$ |
| $10^6$ | $10^{-5}$ | $3.2 \times 10^{-9}$ | $7.2 \times 10^{-9}$ |

The step 1.0 kg is exactly representable, so ten subtractions are exact and the drift is zero. Over the five non-zero points the log-log slope is $1.17 \pm 0.03$ (leave-one-out refits range from 1.15 to 1.21). The drift grows *linearly* with $n$, not like $\sqrt n$, and sits a factor 2 to 10 below the bound.

**Explanation.** The step $10/n$ is stored with a fixed representation error $\delta$; every subtraction applies the same $\delta$ with the same sign, so after $n$ steps the systematic part is $n\delta$. Random-walk cancellation would require independent errors with random signs, which these are not. With a step of exactly 0.125 kg and $n = 80$, every intermediate mass is a multiple of 0.125 inside $[40, 50]$, all exactly representable, and the drift is exactly 0.

**Check.** The bound is never violated. The $n = 100$ drift of $7.1 \times 10^{-14}$ kg is exactly 10 ulps of the final mass (ulp(55) $= 2^{-47} = 7.1 \times 10^{-15}$): a hundred subtractions each contributing a fraction of an ulp with the same sign, consistent with systematic accumulation.

**Common mistakes.** Including the zero-drift point in the log fit; comparing against `eps` instead of `eps/2` without saying so (both are acceptable if stated); concluding from one run that rounding errors "always" grow linearly: they do here because the error is systematic.

## 6. When the regression test fails

(a) **Neither is wrong.** You changed the model deliberately; the stored table describes the old model. Regenerate `REFERENCE_POLAR` from the new code *after* confirming, with `max_lift_to_drag` and the optimum property test, that the new polar is self-consistent, then commit the table with a message such as "Update reference polar for new wing: K_DEFAULT 0.050 -> 0.045 (AR 8.9); (L/D)max moves from 15.81 to 16.67; regenerated REFERENCE_POLAR". Keep the old table in the history, not in the file.

(b) A relative difference of $3 \times 10^{-7}$ is a thousand times larger than rounding (a handful of ulps, $\sim 10^{-16}$) and ten times larger than the table's 8-figure rounding ($\sim 10^{-8}$), yet far smaller than any formula error (which would be percent-level). The two likely causes are: a changed constant written with fewer digits (for example `0.3048` typed as `0.305` somewhere), or a changed order of operations that moved the result across the table's rounding boundary combined with an unnoticed change to the stored precision. The distinguishing measurement is to print the new value to 17 significant figures and compare with the *unrounded* value from the previous commit: a few-ulp difference points at operation order and means the table tolerance was too tight; a $10^{-7}$ difference in the full-precision value points at a changed constant.

(c) The reference's own error is up to about $10^{-8}$ relative (and `rtol` is a *per-step* request, so the global error can be somewhat larger). The tightest honest tolerance is about $10^{-7}$. Tighter would make the test fail on a correct result whenever the integrator's rounding changed (new SciPy version, new CPU), and you would then be forced to loosen it under pressure, which is how tolerances become arbitrary. A tolerance must be derived from the uncertainty of the reference, not from what happens to pass today.

## 7. Does the suite catch planted bugs?

| Mutation | Failing tests |
|---|---|
| 1. `cd0 + k * cl` | `test_drag_polar_limiting_cases`, all four `test_max_lift_to_drag_is_the_maximum[...]`, `test_polar_table_matches_stored_reference` (6) |
| 2. `rho * v ** 2` | `test_dynamic_pressure_analytic_value` only (1) |
| 3. missing 2 in $(L/D)_{\max}$ | all four `test_max_lift_to_drag_is_the_maximum[...]` (4) |
| 4. `< 0.0` for `<= 0.0` | three `test_lift_coefficient_rejects_invalid_inputs[...]` cases (zero rho, zero v, zero area) and `test_other_invalid_inputs_raise` (4) |
| 5. `min` for `max` in `is_close` | none (0) |

Mutation 2 is instructive: the scaling tests (`q \propto V^2`, `q \propto \rho`), the round-trip test and the vectorized test all pass with the wrong constant because they only check ratios. One hand-computed value is what saves the suite, exactly as in problem 4.

Mutation 5 escapes because every case in `test_is_close_matches_math_isclose` has $|a| \approx |b|$. A test that catches it:

```python
def test_is_close_uses_the_larger_magnitude():
    # |1 - 2| = 1 <= 0.5 * max(1, 2) = 1 -> close; with min(1, 2) the bound is 0.5 -> not close.
    assert impl.is_close(1.0, 2.0, rel_tol=0.5)
    assert impl.is_close(1.0, 2.0, rel_tol=0.5) == math.isclose(1.0, 2.0, rel_tol=0.5)
```

No tolerance is needed beyond the one under test, because the inputs are exact binary fractions and the comparison is a boolean.

**Common mistakes.** Running the mutated code against the solution tests in the wrong directory and seeing nothing fail; "fixing" a mutation by weakening the test that caught it.

## Check-yourself answers

1. 0.1, 0.2 and 0.3 are repeating binary fractions, so each is stored rounded to 53 significant bits; the rounded sum $0.1 + 0.2$ lands on the float just above the stored 0.3. The difference is one ulp near 0.3, $5.6 \times 10^{-17}$, about $2 \times 10^{-16}$ relative.
2. About `rel=1e-12`: the formula is three exact-input operations, so the arithmetic error is a few ulps ($\sim 10^{-16}$ relative), and 1e-12 leaves four orders of magnitude of margin. `1e-3` would let a 0.1 % error through, which is a physical error, not a rounding error, and would hide a mistyped density.
3. `abs_tol` (for example `abs_tol=1e-12`). It is right when the *expected* value is zero or near zero and you know the scale of the quantities that were subtracted to produce it; a relative tolerance cannot pass near zero because `rel * 0 = 0`.
4. Check the new output against an independent reference (the analytic optimum, a hand value), then regenerate `REFERENCE_POLAR` and commit it with a message that names the bug, the old and new values, and why the new values are right. Never edit individual table entries to make the test pass.
5. $C_D \to C_{D0}$, independent of $C_L$: a wing with no induced drag, which is the infinite-aspect-ratio (2-D aerofoil) limit.
6. `dynamic_pressure` at $V = 0$ is a legitimate value (zero, a parked aircraft), while `lift_coefficient` divides by $q$, so $V = 0$ makes the coefficient undefined; the function must refuse rather than return `inf` or `nan`.
7. $C_L^* = \sqrt{0.03/0.06} = 0.7071$; $(L/D)_{\max} = 1/(2\sqrt{0.0018}) = 11.79$. `max_lift_to_drag(0.03, 0.06)` returns `(0.7071067811865476, 11.78511301977579)`.
8. Passing tests show the absence of the specific failures those tests look for, for those inputs; they do not prove the formula. The parametrized test checks a property (the optimum bounds the sampled curve) for four polars; a bug that only appears for, say, $k > 1$ or that scales the output by a constant (problem 7, mutation 2) can survive. Proof comes from the derivation; tests guard it.
