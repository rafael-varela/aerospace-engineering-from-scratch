# P00.07  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. $W = m g$ with $m = 600 + 10s$. 2. $q = \tfrac12 \rho V^2$; check units give Pa. 3. $C_L = W/(qS)$. 4. $C_D = C_{D0} + k C_L^2$, then $L/D = C_L/C_D$. 5. $D = C_D\, q\, S$, or equivalently $D = W/(L/D)$: compute both and confirm they agree to rounding. 6. Call `lift_coefficient(W, 1.225, 40.0, 12.0)` and compare with `pytest.approx` at `rel=1e-12`.

## 2.
(a) 1. Write $f = C_L/(C_{D0}+kC_L^2)$. 2. Quotient rule; the denominator squared is positive, so only the numerator $C_{D0} - kC_L^2$ can vanish. 3. Solve for $C_L^*$; substitute back to get $C_D^*$ and $(L/D)_{\max}$. (b) 1. Power $\propto D V \propto C_D/C_L^{3/2}$ (state why: $V \propto C_L^{-1/2}$ at fixed weight). 2. Differentiate $C_{D0}C_L^{-3/2} + kC_L^{1/2}$ and set to zero. 3. Compare with (a) by taking the ratio. (c) Name the assumption (polar valid for all $C_L$; no stall) and reason about what a test comparing analytic $C_L^*$ with a sampled maximum would do if the sampled polar were capped at $C_{L,\max}$.

## 3.
(a) Doubling algorithm on 0.1; stop when a fractional part repeats. (b) $1531.25$ lies in $[2^{10}, 2^{11})$, $0.3$ in $[2^{-2}, 2^{-1})$; ulp $= 2^{e-52}$. (c) Representation error of 1.225 is at most half an ulp of 1.225; multiplied by $1250$ it becomes an absolute error in the product; add half an ulp of the result for the final rounding; convert to ulps of 1531.25. Then run the expression and subtract. (d) Compare the bound with $10^{-12} \times 1531.25$ and with $10^{-15} \times 1531.25$.

## 4.
1. Known case: choose inputs that make the answer an integer (unit $q$, unit area). 2. Limiting case: zero force gives zero coefficient exactly, or doubling $q$ halves the coefficient. 3. Invalid input: $q \le 0$, area $\le 0$, NaN must raise. 4. For the stray factor: ask which tests are sensitive to the absolute scale of the output and which only to ratios.

## 5.
1. Loop over $n$; `mdot = burn / n`, `dt = 1.0` (or any split with `mdot * dt = burn / n`). 2. `drift = abs(stepwise[-1] - closed_form)`. 3. `bound = n * eps / 2 * m0`. 4. Keep only non-zero drifts; `np.polyfit(np.log10(n), np.log10(drift), 1)` for the slope; estimate its uncertainty by refitting with one point removed at a time. 5. Repeat with `mdot * dt = 0.125` and $n = 80$.

## 6.
(a) Decide whether the change in output was intended; if so regenerate the table from the new code *after* checking it against the analytic optimum, and record both the reason and the old value in the commit. (b) Compute the relative difference; compare with the table's 8-figure rounding ($\sim 10^{-8}$) and with the tolerance ($2 \times 10^{-7}$). (c) The reference carries the integrator's error; the regression tolerance cannot be below it. Think about what happens when a tolerance is tighter than the reference's own uncertainty.

## 7.
1. Copy `main.py` to a scratch directory with the tests and `conftest.py`; apply one mutation; run `pytest -q`; record failing test names. 2. Repeat for each mutation. 3. For the uncaught one, construct two numbers whose closeness verdict depends on whether the larger or smaller magnitude multiplies the relative tolerance.
