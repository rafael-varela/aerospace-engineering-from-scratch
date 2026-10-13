# P00.06  SciPy Tour

> Never trust a library call you have not first checked against an answer you can derive.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the total impulse, peak thrust and burnout time of a simple thrust model by hand. 2. Implement and test the same quantities with `scipy.integrate`, `scipy.optimize`, `scipy.interpolate` and `scipy.linalg`, showing each result agrees with the derivation to a stated tolerance. 3. Measure the convergence rate of the trapezoidal rule and the error of three interpolation schemes against a known function. 4. Decide which SciPy submodule a given engineering question belongs to, and name the Phase 03 lesson where you will build that method yourself. |
| **Prerequisites** | `P00.03`, `P00.04`, `P00.05`; skill check: you can differentiate `t e^{−t/τ}` and integrate it by parts; you can solve a 3×3 linear system by hand. |
| **Engineering question** | A model-rocket motor's thrust was sampled 81 times over two seconds: what is its total impulse and letter class, when does thrust peak, when has it effectively burned out, and how sure are you of each number? |
| **Assumptions and limits** | The thrust samples are *synthetic*: generated from `F = a t e^{−t/τ}` plus Gaussian noise, so every answer has an analytic reference. Real motor curves have ignition delays and tail-off that this model does not represent. The ISA table is exact at its nodes; interpolation error between nodes is the quantity under study. The wind-tunnel balance is linear and its geometry is known exactly. |
| **Analytical work** | Closed-form impulse by parts; peak by differentiation; limits as `t_end → ∞`; the balance equations from force and moment equilibrium; why log-linear interpolation is exact for an isothermal atmosphere. |
| **Computational deliverable** | `main.py` with `total_impulse_samples` (trapezoid), `total_impulse_model` (quad), `trapezoid_convergence`, `fit_thrust_model` (curve_fit), `peak_thrust` (minimize_scalar), `burnout_time` (brentq), `interpolate_pressure` (three schemes), `resolve_balance` (solve). Tests prove: agreement with analytic values to 10⁻⁸, second-order trapezoid convergence, exact parameter recovery on noise-free data, monotone interpolation with log-linear beating linear, and an exact balance round trip. |
| **Practical activity** | **computational experiment**: a convergence study of the trapezoidal rule and an error comparison of three interpolation schemes against the closed-form atmosphere. |
| **Equipment** | None. |
| **Safety** | None beyond normal computer use. No motor is fired in this lesson; the data are simulated. |
| **Evidence** | Your `main.py`, passing tests, the convergence table, the interpolation error figure, and the report printed by `main()`. |
| **Assessment** | `assessments/phase-00/grading-rubric.md`; competencies 1, 3 and 4; pass: all public tests green, Problem 2 derivation complete, and a correct mapping of each question in Problem 6 to a submodule with a reason. |
| **References** | `@scipy-userguide` "Integration (scipy.integrate)" (required-free); "Optimization (scipy.optimize)" (required-free); "Interpolation (scipy.interpolate)" (required-free); "Linear Algebra (scipy.linalg)" (required-free). Supplementary: "Statistics (scipy.stats)", "Special Functions (scipy.special)", "Sparse Arrays (scipy.sparse)" (supplementary-free), for the questions in Problem 6 that fall outside the four submodules used here. |
| **Extensions** | ★ Replace the trapezoid with `scipy.integrate.simpson` and measure its convergence order. ★ Fit a two-phase model (linear rise, exponential tail) and compare the information criteria. |
| **Time** | unmeasured |

## 1. The engineering problem

A hobby-rocketry club has a test stand with a load cell. They sample thrust at 40 Hz and send you a file of 81 numbers covering two seconds. They want the motor's total impulse (which fixes its letter class and therefore which rockets it may legally fly in, under the safety code you meet in P00.13), its peak thrust (which fixes the launch-rod and fin loads), and its burnout time (which fixes when the altimeter should arm). They also want to know how much to believe each number given the noise on the load cell.

Each of those is a numerical task that SciPy does in one call: integrate the samples, fit a curve and find its maximum, find where it crosses a threshold. Later in the same week the club asks you to look up the air density at the launch altitude from a coarse table and to recover lift, drag and moment from the three load cells of a small wind-tunnel balance: an interpolation and a linear solve. Four questions, four SciPy submodules.

The engineering risk is not that SciPy gets the arithmetic wrong; it is that you call the wrong function, pass the wrong argument, or trust a result whose accuracy you never checked. This lesson protects against that in the only reliable way: every call is made first on a case whose answer you have derived by hand, and the agreement is asserted in a test. Because real load-cell data has no exact answer, the lesson *synthesizes* the data from a known model with known noise, and says so throughout. The success criterion is that every number `main()` prints matches a value you can write on paper, to a tolerance you chose.

## 2. Understand and derive

### 2.1 A thrust model with closed-form answers

Take `F(t) = a t e^{−t/τ}` for `t ≥ 0`, with `a` in N/s and `τ` in s. It rises linearly from zero, peaks, and decays: a crude but serviceable shape for a small composite motor. Its advantages are entirely mathematical.

*Peak.* `dF/dt = a e^{−t/τ} (1 − t/τ) = 0` at `t = τ`, where `F_max = a τ / e`. With `a = 60 N/s`, `τ = 0.3 s`: peak **6.62 N at 0.30 s**.

*Total impulse to time T.* Integrate by parts with `u = t`, `dv = e^{−t/τ} dt`:

```
∫₀ᵀ a t e^{−t/τ} dt = a [ −τ t e^{−t/τ} ]₀ᵀ + a τ ∫₀ᵀ e^{−t/τ} dt
                   = a τ² [ 1 − e^{−T/τ} (1 + T/τ) ].
```

Units: `(N/s)(s²) = N·s`. As `T → ∞` the bracket → 1 and `I → a τ² = 5.40 N·s`; at `T = 2 s`, `T/τ = 6.67`, `e^{−6.67}(1 + 6.67) = 0.0098`, so **I = 5.347 N·s**. A hobby-motor class table doubles each letter: 2.5–5 N·s is "B", 5–10 N·s is "C". This is a C motor, and a 1 % error in impulse would not change that; a 7 % error, near the boundary, might.

*Burnout.* There is no instant where thrust is zero; define burnout as the time after the peak when `F` falls below a threshold. For `F = 1 N`: solve `60 t e^{−t/0.3} = 1`, which has no closed form but is a monotone function on `[τ, T]`, so a bracketing root finder converges without fail: **t = 1.309 s**.

These three numbers are the reference the tests compare against.

### 2.2 How wrong is the trapezoidal rule?

With samples at spacing `h`, the trapezoidal rule approximates each strip by a trapezium. For a smooth function the error of the composite rule is `−(T h² / 12) f''(ξ)` for some `ξ` in the interval: halving `h` divides the error by four. That "order 2" behaviour is what you measure in part 3; if you observe order 1, the samples are not evenly spaced or there is a discontinuity. Since `F` here is smooth and `F(0) = 0`, the rule behaves. For 81 samples (`h = 0.025 s`) the error is about 3 × 10⁻³ N·s, 0.06 % of the impulse: far below the noise. `quad`, SciPy's adaptive quadrature, evaluates the *model* wherever it likes and reaches 10⁻¹⁰ relative; it needs a function, not samples, which is the distinction between the two calls.

### 2.3 Fitting and the information in noise

`curve_fit` finds the `(a, τ)` that minimise the sum of squared residuals between the model and the samples. On noise-free samples it recovers the generating parameters to machine precision, which is the first test of any fitting code. On samples with noise of standard deviation 0.3 N (about 5 % of the peak) it recovers `a` to within ±1 N/s and `τ` to within ±0.005 s, and the fitted impulse to within 1–2 %: the fit *averages* the noise over 81 samples, so its impulse estimate is tighter than the single noisy integral. The initial guess matters: this model has a trivial local minimum at `a = 0`, so the starter derives `τ₀` from the time of the largest sample and `a₀` from the peak height, both of which the data provide.

### 2.4 Interpolating a table: linear in what?

The coarse ISA table gives pressure every 2 km. Between nodes, *linear* interpolation of `p` overestimates in a convex, decaying function (the chord lies above the curve), by up to 1.2 % here. In the isothermal layer `ln p` is exactly linear in `h` (P00.04 §2.2), so interpolating `ln p` linearly is *exact* there and nearly so in the troposphere: error 0.21 %. A cubic spline on `ln p` reaches 0.015 %. The lesson is general: interpolate the variable in which the physics is closest to linear. The scheme also decides what happens at the table edge; the module refuses to extrapolate because a straight line in `ln p` continued above 20 km misses the temperature rise in the stratosphere, and silent extrapolation is how tables lie.

### 2.5 Which submodule, and where you will build it yourself

| Question | Submodule and call | Where you build the method | What SciPy does *not* tell you |
|---|---|---|---|
| Impulse from samples | `integrate.trapezoid` | P03.04 Numerical Differentiation & Integration | The error, unless you measure convergence |
| Impulse of a model function | `integrate.quad` | P03.04 | Whether the function has a singularity it stepped over |
| Best-fit parameters | `optimize.curve_fit` | P03.14 Unconstrained Optimization | Whether the minimum is the global one |
| Peak of a curve | `optimize.minimize_scalar` | P03.14 | Anything outside the bounds you gave |
| Threshold crossing | `optimize.brentq` | P03.02 Root Finding | Which root, if there are several in the bracket |
| Value between table rows | `interpolate.interp1d`, `CubicSpline` | P03.03 Interpolation & Splines | How wrong it is between nodes |
| Forces from balance readings | `linalg.solve` | P03.08 Direct Linear Solvers | How much the noise was amplified (ask `cond`) |

The last column is the lesson. Every call returns a number; none returns a warranty. The derivation in each row of this section is how you buy one.

### 2.6 A three-cell balance as a linear system

A model on a sting is held by two vertical load cells, forward and aft, at `±ℓ` from a reference point, and one axial cell. Vertical force equilibrium: `r_fwd + r_aft = L`. Moment equilibrium about the reference point: `ℓ r_fwd − ℓ r_aft = M`. The axial cell reads drag plus a small interaction, 2 % of lift from a slight misalignment, which real balances have and calibrate out: `r_drag = D + 0.02 L`. Solving the first two for the cell readings,

```
r_fwd = L/2 + M/(2ℓ),    r_aft = L/2 − M/(2ℓ),    r_drag = 0.02 L + D,
```

so `r = B x` with `x = [L, D, M]ᵀ` and

```
B = [ 0.5    0    0.5/ℓ ]
    [ 0.5    0   −0.5/ℓ ]
    [ 0.02   1     0    ]
```

Recovering `x` from readings is `x = B⁻¹ r`, which `scipy.linalg.solve` does without forming the inverse. The *condition number* of `B` (ratio of its largest to smallest singular value) bounds how much relative error in `r` can be amplified in `x`: 6.7 for `ℓ = 0.15 m`, 67 for `ℓ = 0.015 m`. A short moment arm makes the moment hard to resolve, exactly as intuition says: two nearly equal readings must be subtracted.

Check by hand: readings `[50, 50, 2]` with `ℓ = 0.15`: `L = 100`, `M = 0.15 × (50 − 50) = 0`, `D = 2 − 0.02 × 100 = 0`.

## 3. Build the computational model

Open `starter/main.py`. The functions are grouped by submodule; `isa_pressure_exact` is provided (it is the P00.03 model).

### Step 1: the model and its analytic answers

```python
def thrust_model(t, a, tau):
    t = np.asarray(t, dtype=float)
    return np.where(t >= 0.0, a * t * np.exp(-t / tau), 0.0)

def analytic_impulse(a, tau, t_end):
    x = t_end / tau
    return a * tau**2 * (1.0 - np.exp(-x) * (1.0 + x))
```

`synthesize_thrust_data` samples the model on `np.linspace(0, t_end, n)` and adds `rng.normal(0, noise)` from a seeded generator; `noise = 0` returns the exact model so tests can compare with analytic values.

### Step 2: integrate

```python
from scipy import integrate

def total_impulse_samples(t, f):
    return float(integrate.trapezoid(f, t))

def total_impulse_model(a, tau, t_end):
    value, err = integrate.quad(lambda t: thrust_model(t, a, tau), 0.0, t_end)
    return float(value)
```

`trapezoid` takes `(y, x)` in that order, a classic source of silent error; the test that compares against the analytic impulse catches it. `trapezoid_convergence` loops over sample counts and returns `|I_n − I_exact|`.

### Step 3: optimize

```python
from scipy import optimize

def fit_thrust_model(t, f):
    i = int(np.argmax(f))
    tau0 = max(t[i], 1e-3)
    a0 = max(f[i], 1e-3) * np.e / tau0
    (a, tau), cov = optimize.curve_fit(thrust_model, t, f, p0=(a0, tau0))
    return float(a), float(tau)

def peak_thrust(a, tau, t_end):
    res = optimize.minimize_scalar(lambda t: -thrust_model(t, a, tau),
                                   bounds=(0.0, t_end), method="bounded")
    return float(res.x), float(-res.fun)
```

`minimize_scalar` minimises, so pass the negative thrust. `burnout_time` brackets the root between `τ` (the peak) and `t_end` and calls `brentq`; it raises `ValueError` if the threshold is above the peak or not yet reached, because a root finder given a bracket with no sign change will either fail or, worse, be asked to invent one.

### Step 4: interpolate

```python
from scipy import interpolate

def interpolate_pressure(h, kind="log-linear"):
    if np.any(h < ISA_TABLE_H[0]) or np.any(h > ISA_TABLE_H[-1]):
        raise ValueError("refusing to extrapolate")
    if kind == "linear":
        return interpolate.interp1d(ISA_TABLE_H, ISA_TABLE_P)(h)
    if kind == "log-linear":
        return np.exp(interpolate.interp1d(ISA_TABLE_H, np.log(ISA_TABLE_P))(h))
    if kind == "cubic":
        return np.exp(interpolate.CubicSpline(ISA_TABLE_H, np.log(ISA_TABLE_P))(h))
```

All three reproduce the table nodes exactly; they differ only between nodes, which is what the tests measure.

### Step 5: linalg

```python
from scipy import linalg

def balance_matrix(arm=0.15):
    return np.array([[0.5, 0.0, 0.5 / arm], [0.5, 0.0, -0.5 / arm], [0.02, 1.0, 0.0]])

def resolve_balance(readings, arm=0.15):
    return linalg.solve(balance_matrix(arm), np.asarray(readings, float))
```

`linalg.solve` factorises `B` and back-substitutes; it is both faster and more accurate than `inv(B) @ r`. `balance_condition_number` returns `np.linalg.cond(B)`.

### Step 6: the report

`main()` ties the four parts together on one synthetic dataset and prints, for each submodule, the SciPy result beside the analytic one:

```
integrate:
  trapezoid on noisy samples :   5.4166 N*s  -> class C
  quad on the model          :   5.3473 N*s
  analytic                   :   5.3473 N*s
  observed convergence rates : [1.98 1.99 2.   2.  ]  (trapezoid: 2)
optimize:
  fitted a = 58.015 N/s (true 60.0), tau = 0.3060 s (true 0.3)
  peak 6.532 N at t = 0.3060 s   (analytic a*tau/e = 6.622 N at 0.3 s)
```

Read it as a reviewer would: the noisy trapezoid is 1.3 % above the exact impulse, the fit is 3 % low on `a` and 2 % high on `τ` (the two parameters trade off along the ridge of the least-squares surface), and the peak inherits both. None of these differences is a bug; each is the noise doing what noise does, and the report puts the analytic value next to every number so that a real bug would stand out.

### Verify

- *Analytic:* `total_impulse_model` equals `analytic_impulse` to 10⁻⁸; `peak_thrust` returns `(τ, aτ/e)` to 10⁻⁶; `burnout_time` returns a point where the model equals the threshold to 10⁻⁸; `fit_thrust_model` on noise-free data recovers `(a, τ)` to 10⁻⁶ relative.
- *Limiting cases:* `analytic_impulse(a, τ, 0) = 0`; `analytic_impulse(a, τ, 50τ) → aτ²`; `thrust_model` is zero for `t < 0`; the balance with equal vertical readings gives zero moment.
- *Convergence:* trapezoid errors for `n = 11, 21, 41, 81, 161` fall by a factor 4 per doubling (observed rates 1.98, 1.99, 2.00, 2.00).
- *Conservation/linearity:* `analytic_impulse` is linear in `a`; `resolve_balance(r₁ + r₂) = resolve_balance(r₁) + resolve_balance(r₂)`; the residual `B x − r` is zero to 10⁻¹².
- *Monotonicity and ordering:* all three interpolants decrease with altitude; log-linear beats linear and cubic beats log-linear in maximum relative error; the motor-class letter never decreases with impulse.

### Compare

The trapezoid on 81 *noisy* samples gave 5.417 N·s against the exact 5.347: the 1.3 % difference is the noise, not the rule (the noise-free trapezoid error is 0.06 %). The fitted model's impulse is typically within 1–2 % because the fit averages the noise. Compare also `integrate.simpson` on the same noise-free samples: its error is 4 × 10⁻⁶ N·s at `n = 81` against the trapezoid's 3 × 10⁻³, and it converges at order 4. Simpson is "better" on smooth model data; on noisy measured data the noise dominates both, and the trapezoid's simplicity (and its honesty about not inventing curvature between samples) is why thrust-curve standards specify it. Phase 03 derives both rules and their error terms.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Two studies, both of numerical methods against a known truth.

*Study A: convergence.* Before running, predict the trapezoid error at `n = 321` from the error at `n = 161` (what factor, and why?). Run `trapezoid_convergence` for `n = 11, 21, 41, 81, 161, 321`, tabulate errors and the observed rate `log₂(e_n / e_{2n})`, and plot error against `h` on log–log axes with your P00.04 module. Then repeat with a thrust curve that *cuts off abruptly*: the same model for `t < 1.0 s` and exactly zero after (a motor whose nozzle is ejected), with the exact impulse from `analytic_impulse(a, τ, 1.0)`, and observe what the jump does to the rate.

*Study B: interpolation.* Predict which of linear, log-linear and cubic will have the smallest maximum error and where along the altitude range the maximum will sit. Evaluate all three on 2 001 points, record maximum relative error and its altitude, and plot the relative error curves. Then coarsen the table to 4 km spacing and repeat: by what factor does each scheme's error grow?

Raw-data format: `study, scheme_or_n, max_rel_error, location, observed_rate`.

## 5. Validate and improve

Study A: the predicted factor was 4 (second order); the observed is 3.99–4.00 for the smooth model. For the cut-off curve the observed rate falls to 1.4, 1.3, 1.2, 1.1, 1.0: first order, because the trapezium that straddles the jump has an error proportional to `h`, not `h²`, no matter how smooth the rest is. This is not a failure of the rule; it is the rule telling you the function is discontinuous, and the fix is to integrate the two pieces separately (or to place a sample exactly on the jump). Say which you did.

Study B: linear's maximum error sits near 19 km, where the chord over the steeply curving `p(h)` is worst relative to the small `p`; log-linear's sits near 9 km in the troposphere, where `ln p` is most curved; cubic's is smallest and sits near the tropopause, where the slope of `ln p` is discontinuous and no smooth spline can be exact. Doubling the table spacing to 4 km roughly quadruples the linear and log-linear errors (1.25 % → 5.1 %, 0.21 % → 0.77 %; both are second-order schemes between nodes) and multiplies the cubic error by about fifteen (0.015 % → 0.23 %), close to the sixteen a fourth-order scheme predicts; the cubic's maximum sits at the tropopause, where `ln p` has a slope discontinuity no smooth spline can follow.

Uncertainty: the noisy-impulse test accepts a 5 % error. Estimate your own: run `synthesize_thrust_data` with twenty seeds, compute the trapezoid impulse for each, and report mean and standard deviation. Expect a standard deviation near 0.07 N·s, about 1.3 %. Predict it first: the trapezoid is a weighted sum with weights `h` (and `h/2` at the ends), so independent noise of 0.3 N gives `σ_I = 0.3 × h × √(n − 1.5) ≈ 0.3 × 0.025 × 8.9 = 0.067 N·s`. If yours is much larger, check that your noise is not being added twice.

Improve: pass `sigma` to `curve_fit` so that it weights samples by their stated uncertainty, and use the returned covariance to put error bars on `(a, τ)`. Re-run and see whether the fitted impulse's spread across seeds shrinks.

## 6. Defend your work

Submit: `main.py` and test output; the two study tables and figures; the twenty-seed impulse statistics; and a half-page answer to the engineering question with a stated uncertainty on each of impulse, class, peak and burnout. Use `projects/_template/report-template.md`.

### Adversarial self-checks

- What if the sample rate doubles? The trapezoid error on the model falls by 4; the noise contribution to the impulse falls by √2 only. Which dominates at 40 Hz? At 400 Hz?
- What if `τ` doubles? The peak halves in height and moves to 0.6 s; the impulse to 2 s rises toward `a τ²` but the tail is cut off (`T/τ = 3.3`, bracket 0.86). Does the class change?
- Which assumption fails first on a real motor? The model shape: real curves have an ignition spike and a flat tail-off. Expect `curve_fit` to compromise and the fitted impulse to be biased. The trapezoid on the raw samples does not care about the shape, which is why it is the standard.
- What measurement could falsify the balance model? Apply a known pure moment (equal and opposite weights) and read a non-zero drag cell: the 2 % interaction coefficient is wrong or there is a second interaction.
- Does the result survive a change of interpolation scheme? For pressure at 5.5 km the three schemes differ by 1 %; a density used in a drag estimate inherits that. Would your drag conclusion change by 1 %?
- Which decision changes if the load-cell noise doubles? The class, if the impulse is within ~2 % of a boundary; otherwise none. Compute how far 5.35 N·s is from the C/D boundary at 10 N·s before answering.

## Check yourself

1. Which argument order does `scipy.integrate.trapezoid` take, and what happens if you swap them with evenly spaced `t`?
2. Why does `total_impulse_model` use `quad` while `total_impulse_samples` uses the trapezoid? Could you swap them?
3. State the exact peak time and peak value of `F = a t e^{−t/τ}` and the limit of the impulse as `t_end → ∞`.
4. The trapezoid errors at `n = 41` and `n = 81` are 1.26 × 10⁻² and 3.15 × 10⁻³ N·s. What is the observed order, and what would an order of 1 tell you?
5. Why is linear interpolation of `ln p` exact in the isothermal layer and not in the troposphere?
6. Give the condition number's meaning in one sentence, and say what happens to the recovered pitching moment when the balance arm shrinks by ten.
7. `brentq` is asked for a root on `[τ, t_end]` and the function has the same sign at both ends. What should the code do and why?
8. `curve_fit` returns `a = 0.0003, τ = 47`. What went wrong and what prevents it?
