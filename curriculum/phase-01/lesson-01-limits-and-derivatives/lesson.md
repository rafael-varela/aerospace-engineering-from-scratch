# P01.01  Limits, Continuity & Derivatives

> A derivative is a limit you can see on a log-log plot; trust the number only after you have watched it converge.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the derivative of a lift model by the chain rule and state the limit definition it rests on. 2. Implement and test forward, central, five-point and Richardson-extrapolated finite differences and measure their order of accuracy. 3. Locate a kink (continuous but not differentiable) in a piecewise model with one-sided limits. 4. Decide, from the error-versus-step curve, which step size to trust a numerical derivative at, and explain the two regimes of that curve. |
| **Prerequisites** | `P00.03` (NumPy arrays and broadcasting), `P00.07` (pytest with tolerances), `P00.09` (SymPy `diff`, `lambdify`, and the first look at a central difference). Skill check: you can evaluate tanh and cosh on a calculator, apply the chain rule to sin(3x), and read a log-log plot. |
| **Engineering question** | How much extra lift does a vertical gust produce on an aircraft wing, and how far past the linear range can that "lift per degree" number be trusted? |
| **Assumptions and limits** | The lift model C_L = C_Lmax tanh(a₀α/C_Lmax) is a smooth stand-in for a real section: it has the thin-airfoil slope 2π per radian at α = 0 and a plateau C_Lmax, but no sharp stall, no hysteresis and no Reynolds-number dependence. Finite differences assume the function is smooth at the point: at a kink they return a step-dependent average of the two one-sided slopes. In float64 arithmetic every finite difference is wrong below a step of about 10⁻⁵ because of rounding; the lesson quantifies this rather than hiding it. |
| **Analytical work** | The limit definition of the derivative; one-sided limits at the ISA tropopause; the chain rule applied to the lift model; the Taylor-based truncation error of each difference formula; the rounding error and the optimal step; hand values at α = 10°. |
| **Computational deliverable** | `starter/main.py`: `lift_coefficient`, `lift_curve_slope_exact`, `symbolic_lift_slope`, `lambdify_lift_slope`, `limit_sequence`, `one_sided_limits`, `is_continuous_at`, `difference_quotients`, four difference formulas plus `richardson_extrapolation`, `one_sided_derivatives`, `isa_temperature_two_layer`, `error_vs_step`, `observed_order`, `optimal_step_central`, `gust_lift_increment`. Tests prove: symbolic and numerical derivatives agree, each formula is exact on the polynomial degree it should be, observed orders are 1, 2 and 4, Richardson equals the five-point stencil, the ISA kink is detected, and rounding takes over below the optimal step. |
| **Practical activity** | **computational experiment**: map the error of four derivative formulas against step size over eleven decades, fit the slopes on each side of the minimum, and compare the minimum with the predicted optimal step. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `derivative_error_vs_h.png` with its table, the hand calculation at α = 10° against the code, and a short note on the two regimes of the error curve with fitted slopes. |
| **Assessment** | `assessments/phase-01/grading-rubric.md`; competencies 1 (reasoning), 2 (analytical), 3 (computational) and 4 (interpretation). Pass: the chain-rule derivation is complete; all public tests pass; the hand value at 10° agrees with the code to four figures; the error-curve note names both regimes with fitted slopes and the predicted optimum. |
| **References** | `@openstax-calculus-volume-1` §2.2 "The Limit of a Function", §2.3 "The Limit Laws", §2.4 "Continuity", §3.1 "Defining the Derivative", §3.2 "The Derivative as a Function", §3.3 "Differentiation Rules", §3.5 "Derivatives of Trigonometric Functions", §3.6 "The Chain Rule", §3.9 "Derivatives of Exponential and Logarithmic Functions" (all required-free). `@openstax-calculus-volume-1` §2.1 "A Preview of Calculus", §2.5 "The Precise Definition of a Limit", §3.7 "Derivatives of Inverse Functions", §3.8 "Implicit Differentiation" (supplementary-free). |
| **Extensions** | ★ Implement the complex-step derivative Im f(x + ih)/h, which has no subtraction and no rounding floor, and add it to the error plot. ★ Replace the tanh model with tabulated lift data from a NACA report (`P00.11`) and differentiate the table with `np.gradient`; explain why the result is noisier. |
| **Time** | unmeasured |

## 1. The engineering problem

An aircraft in cruise meets a vertical gust. For a moment the air arrives from slightly below, the angle of attack α (the angle between the wing chord and the oncoming air) increases, and the lift jumps. The structure, the passengers and the autopilot all care about the size of that jump. The quantity that converts "a bit more angle" into "this much more lift" is the **lift-curve slope**, dC_L/dα: the derivative of the lift coefficient with respect to angle of attack. Thin-airfoil theory (Phase 10) says the slope is 2π per radian, about 0.110 per degree, for small angles. Near stall it is much smaller; past stall it is negative. A gust analysis done with the wrong slope is wrong by the same factor.

The engineering question is therefore: given a lift curve C_L(α), what is its derivative at the flight condition, and how far can a linear estimate ΔL ≈ (dL/dα) Δα be pushed before it misleads? To answer it we need the derivative as a tool: as a definition (a limit), as a set of rules (chain rule and friends), as a symbolic object (SymPy) and as a number computed from samples (finite differences). The last one matters most in practice, because most lift curves in your career will arrive as tables from a wind tunnel or a CFD run, not as formulas, and a finite difference is the only derivative a table has.

Success criteria, written before we start: the symbolic derivative of the lift model and the hand derivative agree exactly; a central difference with step 10⁻⁴ rad reproduces it to 10⁻⁶ per radian; the observed orders of the forward, central and five-point formulas are 1.0, 2.0 and 4.0 to within 0.1 where truncation dominates; the step at which rounding takes over is predicted to within a factor of three; and the ISA temperature profile is found to be continuous but not differentiable at 11 km.

## 2. Understand and derive

### 2.1 Limits: what a function is heading for

Write lim_{x→a} f(x) = L for "f(x) can be made as close to L as you like by taking x close enough to a, without x being a". The value f(a) itself does not enter; the function need not even be defined there. The classic example is sin(x)/x at x = 0, undefined at 0 but heading for 1 from both sides:

| h | sin(h)/h |
|---|---|
| 10⁻¹ | 0.9983341665 |
| 10⁻² | 0.9999833334 |
| 10⁻³ | 0.9999998333 |
| 10⁻⁴ | 0.9999999983 |

The error 1 − sin(h)/h is h²/6 to leading order (you can read the 6 from the table: 0.00166 at h = 0.1), which is the first of many times in this phase that an error has a known power law. Limits can be taken from one side only: lim_{x→a⁺} approaches from above, lim_{x→a⁻} from below. The two-sided limit exists only when both one-sided limits exist and agree (`@openstax-calculus-volume-1` §2.2).

The limit laws (§2.3) say limits respect sums, products and quotients (quotients when the denominator's limit is not zero), so most limits are found by substitution. When substitution gives 0/0 you have an **indeterminate form** and must work harder: factor, rationalise, or use the squeeze theorem. The squeeze theorem proves sin(x)/x → 1: for 0 < x < π/2, cos x < sin(x)/x < 1, and both bounds go to 1. P01.02 adds l'Hôpital's rule, and P01.05 shows that a Taylor series answers all of these at once.

### 2.2 Continuity, and a kink at 11 km

f is **continuous at a** when lim_{x→a} f(x) exists and equals f(a): no jump, no hole, no blow-up. Most physical models are continuous by construction, but *piecewise* models are continuous only if someone made them so. The ISA from `P00.03` is the example. In the troposphere T(h) = 288.15 − 0.0065 h; from 11 km the temperature is held at 216.65 K. The two pieces meet: 288.15 − 0.0065 × 11000 = 216.65. T is continuous at 11 km. But its *slope* jumps from −0.0065 K/m to 0: the graph has a corner, a **kink**.

```
   T [K]
 288 ┤╲
     │ ╲   slope −0.0065 K/m
     │  ╲
     │   ╲
 217 ┤    ╲_____________  slope 0
     └────┴──────┴──────▶ h
          11 km   20 km
```

Continuity buys you the **intermediate value theorem** (§2.4): a continuous function that is 0.7780 at α = 8° and 1.0293 at α = 12° takes the value 1.0 somewhere in between. That is how you know the angle for a required C_L *exists* before you ask Newton's method (P01.02) to find it; for the lesson's model it is 11.44°.

### 2.3 The derivative as a limit

The **derivative** of f at x is the limit of the difference quotient (§3.1):

f′(x) = lim_{h→0} [f(x + h) − f(x)] / h.   (1)

It is the slope of the tangent line, the instantaneous rate of change, the sensitivity of output to input. The one-sided versions use h → 0⁺ and h → 0⁻; the derivative exists only if the two agree. At the ISA kink they do not (−0.0065 and 0), so T′(11 km) does not exist, even though T is continuous there. Continuity does not imply differentiability; differentiability does imply continuity (§3.2).

Here is (1) made numerical for the lift model at α = 8°, where the exact slope is 4.342675 per radian:

| h (rad) | [C_L(α+h) − C_L(α)]/h | error |
|---|---|---|
| 10⁻¹ | 3.3007 | 1.04 |
| 10⁻² | 4.2342 | 0.108 |
| 10⁻³ | 4.3318 | 0.0108 |
| 10⁻⁴ | 4.3416 | 0.00108 |
| 10⁻⁵ | 4.34257 | 0.000108 |

The error falls by exactly ten per decade: it is proportional to h. That is a property of this *formula*, not of the function, and §2.6 derives it.

### 2.4 The rules, with the aerospace quantities they differentiate

You will memorise these; the point here is to see each one attached to a quantity you will use (§3.3, §3.5, §3.6, §3.9).

- **Power rule.** Dynamic pressure q = ½ρV². dq/dV = ρV: at 50 m/s and sea level, 61.25 Pa per m/s. Relative to q itself, (dq/dV)/q = 2/V: one extra m/s adds 4 % to the loads at 50 m/s, 2 % at 100 m/s.
- **Product rule.** Lift L = q S C_L with both q(t) and C_L(t) varying in a manoeuvre: dL/dt = S (q̇ C_L + q Ċ_L).
- **Quotient rule.** Glide ratio C_L/C_D: d(C_L/C_D)/dα = (C_L′ C_D − C_L C_D′)/C_D². Its zero is the best-glide angle of attack (P01.02).
- **Chain rule.** Density along a climb: ρ(h(t)), so dρ/dt = (dρ/dh)(dh/dt). Nothing in physics is a function of one variable through only one path; the chain rule is the rule you use most.
- **Trigonometric.** The vertical force on a banked wing L cos φ: d/dφ = −L sin φ. At 60° bank, each extra degree removes sin 60° × (π/180) = 0.015 L from the vertical component, 3 % of the 0.5 L that is left.
- **Exponential and logarithm.** Isothermal pressure p = p₀ e^(−h/H): dp/dh = −p/H. With p = ρRT and H = RT/g₀ this is dp/dh = −ρ g₀, the hydrostatic equation from `P00.03`. The exponential is the only function proportional to its own derivative, which is why exponentials appear wherever a rate is proportional to an amount.

### 2.5 The lift model, differentiated by hand

The model is

C_L(α) = C_Lmax tanh(a₀ α / C_Lmax),   a₀ = 2π rad⁻¹,   C_Lmax = 1.4.   (2)

Why this shape: for small α, tanh u ≈ u, so C_L ≈ a₀α, the thin-airfoil result; for large α, tanh → 1 and C_L → C_Lmax. It is a cartoon of a real lift curve (no sharp stall, no post-stall drop), chosen because it is smooth everywhere and has a one-line derivative. Chain rule with u = a₀α/C_Lmax, du/dα = a₀/C_Lmax, and d(tanh u)/du = sech²u = 1/cosh²u:

dC_L/dα = C_Lmax · sech²(u) · a₀/C_Lmax = **a₀ / cosh²(a₀α/C_Lmax)**.   (3)

Check the limits: at α = 0, cosh 0 = 1 and the slope is a₀ = 6.2832 rad⁻¹ = 0.10966 deg⁻¹. As α → ∞, cosh → ∞ and the slope → 0: the plateau. Units: C_L is dimensionless, so dC_L/dα has units of 1/rad; a₀ has units 1/rad and cosh² is dimensionless. Good.

**Hand value at α = 10°.** α = 0.174533 rad. u = 6.28319 × 0.174533/1.4 = 0.78330. tanh u = 0.65460, so C_L = 1.4 × 0.65460 = 0.9164. sech²u = 1 − tanh²u = 1 − 0.42850 = 0.57150, so dC_L/dα = 6.28319 × 0.57150 = **3.591 rad⁻¹**, 57 % of the thin-airfoil value: at 10° this section has already lost 43 % of its lift-per-degree. The gust that would add 1.44 g from level flight at α = 0 adds 1.00 g from α = 8° (worked in §3).

### 2.6 Finite differences and their truncation error

Nobody evaluates (1) at h = 0. Instead we pick a small h and accept an error. How big? Expand f around x (P01.05 proves this Taylor expansion; for now take it as the polynomial that matches f and its derivatives at x):

f(x ± h) = f(x) ± h f′(x) + (h²/2) f″(x) ± (h³/6) f‴(x) + (h⁴/24) f⁽⁴⁾(x) ± (h⁵/120) f⁽⁵⁾(x) + …

**Forward difference.** [f(x+h) − f(x)]/h = f′ + (h/2) f″ + …: the error is (h/2) f″, **first order**. At α = 8°, C_L″ = −21.66 rad⁻², so at h = 10⁻² the predicted error is 0.108; the table in §2.3 shows 0.108.

**Central difference.** Subtract the two expansions: the f, f″ and f⁽⁴⁾ terms cancel, and

[f(x+h) − f(x−h)]/(2h) = f′ + (h²/6) f‴ + (h⁴/120) f⁽⁵⁾ + …   (4)

The error is (h²/6) f‴, **second order**: halving h divides the error by four. At α = 8°, C_L‴ = −12.85 rad⁻³, so at h = 10⁻² the predicted error is −2.14 × 10⁻⁴ (the estimate comes out low); the code measures −2.16 × 10⁻⁴. For a quadratic, f‴ = 0 and the central difference is *exact*.

**Five-point stencil.** Combine the central differences at h and 2h so the h² terms cancel: [−f(x+2h) + 8f(x+h) − 8f(x−h) + f(x−2h)]/(12h) = f′ − (h⁴/30) f⁽⁵⁾ + …, **fourth order**, exact for quartics.

**Richardson extrapolation.** The same idea as a recipe: if D(h) = f′ + c h² + …, then [4D(h/2) − D(h)]/3 = f′ + O(h⁴). Applied to the central difference it *is* the five-point stencil with step h/2 (expand both and compare); the code checks that identity to rounding. Richardson is worth knowing on its own because it works for any quantity with a known error power law: integrals in P01.03, ODE solvers in Phase 03.

### 2.7 Rounding error and the optimal step

Floating-point numbers carry about 16 significant figures: machine epsilon ε = 2⁻⁵² ≈ 2.2 × 10⁻¹⁶. Each f(x ± h) is rounded to about ε·f, and when the two nearly equal values are subtracted that rounding does not cancel, so the quotient carries an error of order ε·f/h. Truncation falls with h, rounding rises as 1/h, and the total

E(h) ≈ (h²/6) f‴ + ε f/h

is smallest where the two balance: h_opt = (3ε f/f‴)^(1/3). For the lift model at α = 8° (f = 0.778, f‴ = 12.85) that is 3.4 × 10⁻⁶ rad, and the smallest error you can reach is of order 10⁻¹⁰, not zero. "Make h as small as possible" is the single most common mistake with numerical derivatives; part 4 measures it.

### 2.8 Where this stops being valid

- At a kink, formula (4) returns the *average* of the two one-sided slopes, with no warning. At the ISA tropopause a central difference with any h gives −0.00325 K/m, a value the atmosphere never has. Use `one_sided_derivatives` when a model is piecewise.
- The linear estimate ΔL ≈ (dL/dα)Δα is a tangent line; it is only as good as the curve is straight. §3 shows it overshooting by 32 % for a 5.7° gust starting from 8°.
- Noisy data (a wind-tunnel table) has a different trade-off: differentiating amplifies noise by 1/h, so the best step is much larger than for a smooth function. Smoothing or fitting first is usual.

## 3. Build the computational model

Open `starter/main.py`. Each step is runnable on its own.

### Step 1: the lift model and its exact derivative

```python
import numpy as np

A0_THIN_AIRFOIL, CL_MAX = 2.0 * np.pi, 1.4

def lift_coefficient(alpha, a0=A0_THIN_AIRFOIL, cl_max=CL_MAX):
    return cl_max * np.tanh(a0 * np.asarray(alpha, dtype=float) / cl_max)      # (2)

def lift_curve_slope_exact(alpha, a0=A0_THIN_AIRFOIL, cl_max=CL_MAX):
    return a0 / np.cosh(a0 * np.asarray(alpha, dtype=float) / cl_max) ** 2     # (3)
```

`np.asarray(..., dtype=float)` lets a float, a list or an array through, as in `P00.03`.

### Step 2: the same derivative from SymPy

```python
import sympy as sp
alpha_s, a0_s, clmax_s = sp.symbols("alpha a0 C_Lmax", positive=True)

def symbolic_lift_slope():
    return sp.diff(clmax_s * sp.tanh(a0_s * alpha_s / clmax_s), alpha_s)

def lambdify_lift_slope():
    return sp.lambdify((alpha_s, a0_s, clmax_s), symbolic_lift_slope(), modules="numpy")
```

SymPy returns `a0*(1 - tanh(a0*alpha/C_Lmax)**2)`, which is (3) in the identity sech² = 1 − tanh². The test asks `sp.simplify(expr - hand) == 0`; remember from `P00.09` that `==` compares trees, not mathematics.

### Step 3: limits, one-sided

```python
def limit_sequence(f, a, hs, side="right"):
    hs = np.asarray(hs, dtype=float)
    if np.any(hs <= 0.0):
        raise ValueError("every step in hs must be positive")
    if side == "right":
        return f(a + hs)
    if side == "left":
        return f(a - hs)
    raise ValueError("side must be 'left' or 'right'")
```

`one_sided_limits` calls this with h₀·2⁻ᵏ for k = 0…11 and returns the last value on each side; `is_continuous_at` compares both with f(a) to a tolerance. The estimate is only as good as the last step (h₀·2⁻¹¹ ≈ 5 × 10⁻⁶), so the tolerance must be larger than that for a function like |x|. `difference_quotients(f, x, hs)` is definition (1) for a list of h, producing the table of §2.3.

### Step 4: the difference formulas

```python
def _check_h(h):
    if not np.isfinite(h) or h <= 0.0:
        raise ValueError("h must be a positive finite number")

def forward_difference(f, x, h):
    _check_h(h); x = np.asarray(x, dtype=float)
    return (f(x + h) - f(x)) / h

def central_difference(f, x, h):
    _check_h(h); x = np.asarray(x, dtype=float)
    return (f(x + h) - f(x - h)) / (2.0 * h)

def five_point_difference(f, x, h):
    _check_h(h); x = np.asarray(x, dtype=float)
    return (-f(x + 2*h) + 8*f(x + h) - 8*f(x - h) + f(x - 2*h)) / (12.0 * h)

def richardson_extrapolation(f, x, h):
    _check_h(h)
    return (4.0 * central_difference(f, x, 0.5 * h) - central_difference(f, x, h)) / 3.0
```

Every formula is vectorised in `x` for free because `f` is. The `ValueError` on h ≤ 0 is there because `h = 0` produces `nan` silently and `h < 0` produces a plausible but reflected answer, both of which would travel a long way before being noticed.

### Step 5: one-sided derivatives and the ISA kink

```python
def one_sided_derivatives(f, a, h):
    _check_h(h); a = float(a)
    left = 2.0 * backward_difference(f, a, 0.5 * h) - backward_difference(f, a, h)
    right = 2.0 * forward_difference(f, a, 0.5 * h) - forward_difference(f, a, h)
    return float(left), float(right)
```

Each side is a first-order formula improved by one Richardson step (2D(h/2) − D(h) cancels the h term). On a straight line both are exact, so at the kink the code returns (−0.0065, 0.0) to 10⁻⁹ with h = 10 m. `isa_temperature_two_layer` is the two-piece profile with `np.where`, exactly as in `P00.03`.

### Step 6: the error curve and its slope

```python
def error_vs_step(f, dfdx, x, hs):
    exact = float(dfdx(np.asarray(x, dtype=float)))
    return {name: np.array([abs(float(m(f, x, h)) - exact) for h in hs]) for name, m in METHODS.items()}

def observed_order(hs, errs):
    slope, _ = np.polyfit(np.log(hs), np.log(errs), 1)
    return float(slope)
```

`observed_order` is a least-squares straight line through (log h, log error). Fit it only in the truncation regime; include rounding-dominated points and the slope is meaningless. `optimal_step_central(f_scale, f3_scale)` returns (3ε f/f‴)^(1/3). `gust_lift_increment(slope, rho, v, s, u_gust)` is the sensitivity used as an engineer uses it: ΔL = ½ρV²S · slope · (u/V).

### Verify

Run `python -m pytest curriculum/phase-01/lesson-01-limits-and-derivatives`. The tests, in order of what they would catch:

1. **Analytic case.** The symbolic slope simplifies to a₀/cosh²; the lambdified version matches `lift_curve_slope_exact` to 10⁻¹³ and a central difference to 10⁻⁶.
2. **Exactness on polynomials.** Central is exact on 3x² − 2x + 1, five-point on x⁴ − x³ + x/2, and central is *not* exact on the quartic (a test that only checks successes proves little).
3. **Orders.** On eˣ with h = 0.02, 0.01, 0.005 the fitted slopes are 1.00, 2.00 and 4.00. Forward and backward errors have opposite signs and sum to O(h²).
4. **Identity.** Richardson with step h equals the five-point stencil with step h/2 to 10⁻¹².
5. **Limits and continuity.** sin(h)/h → 1 from both sides with error h²/6; a step function is caught, |x| passes, the ISA is continuous at 11 km but has one-sided derivatives (−0.0065, 0).
6. **Rounding.** The central-difference error of sin at x = 1 has its minimum between h = 10⁻⁷ and 10⁻⁴, and the errors at h = 10⁻¹ and 10⁻¹² are both more than 100× larger.
7. **Interface.** Every formula rejects h ≤ 0 and `nan`, and a (3, 4) array in gives a (3, 4) array out.

### Compare

Three independent routes to dC_L/dα at α = 8°: the hand formula (3) gives 4.342675; SymPy's `diff` gives the same expression in a different form; the central difference at h = 10⁻³ gives 4.342673, 2 × 10⁻⁶ low, which is the (h²/6)f‴ = 2.1 × 10⁻⁶ predicted. When the three agree to the predicted error, the derivative is right. When they do not, the disagreement tells you which one is wrong: if the hand and SymPy agree and the numeric does not, check h; if SymPy disagrees with the hand result, `sp.simplify` the difference before assuming a typo.

The gust answer, from `main()`: a 5 m/s vertical gust at 50 m/s is Δα = 0.1 rad = 5.7°. On a 10 kN light aircraft with S = 15 m², the linear estimate from α = 0 is ΔL = 14.43 kN (1.44 g), against the true change in the model of 13.53 kN: the tangent overshoots by 7 %. From α = 8° the linear estimate is 9.97 kN (1.00 g) against a true 7.58 kN: 32 % high, because 8° + 5.7° is close to the plateau. The derivative is a local tool; the second derivative (−21.66 rad⁻² at 8°) is what tells you how local.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Map the error of each formula against step size and find the two regimes.

1. Fix the lift model and α = 8°. For h = 10⁻¹, 10⁻², …, 10⁻¹², record `error_vs_step(lift_coefficient, lift_curve_slope_exact, a8, hs)`.
2. **Before running, write down**: the slope you expect on the left of the minimum for each formula (1, 2, 4), the slope on the right (−1 for all), and h_opt for the central difference from §2.7 (3.4 × 10⁻⁶ rad).
3. Plot the four error curves on log-log axes, add the rounding line ε·0.778/h, and save as `derivative_error_vs_h.png`.
4. Fit `observed_order` separately to h ∈ [10⁻⁴, 10⁻²] and to h ∈ [10⁻¹², 10⁻⁹].

Raw data format: CSV with columns `h, err_forward, err_central, err_five_point, err_richardson`.

What you should see (the author's run):

| h | forward | central | five-point | Richardson |
|---|---|---|---|---|
| 10⁻¹ | 1.04 | 3.4 × 10⁻² | 3.3 × 10⁻² | 3.1 × 10⁻³ |
| 10⁻² | 1.1 × 10⁻¹ | 2.2 × 10⁻⁴ | 5.6 × 10⁻⁶ | 3.5 × 10⁻⁷ |
| 10⁻³ | 1.1 × 10⁻² | 2.1 × 10⁻⁶ | 5.7 × 10⁻¹⁰ | 3.5 × 10⁻¹¹ |
| 10⁻⁴ | 1.1 × 10⁻³ | 2.1 × 10⁻⁸ | 6.8 × 10⁻¹³ | 1.4 × 10⁻¹² |
| 10⁻⁶ | 1.1 × 10⁻⁵ | 6.3 × 10⁻¹¹ | 6.3 × 10⁻¹¹ | 2.9 × 10⁻¹⁰ |
| 10⁻⁸ | 1.0 × 10⁻⁷ | 3.5 × 10⁻⁹ | 3.5 × 10⁻⁹ | 1.8 × 10⁻⁸ |
| 10⁻¹⁰ | 8.9 × 10⁻⁷ | 3.3 × 10⁻⁷ | 7.9 × 10⁻⁷ | 1.1 × 10⁻⁶ |
| 10⁻¹² | 7.3 × 10⁻⁵ | 7.3 × 10⁻⁵ | 8.2 × 10⁻⁵ | 7.3 × 10⁻⁵ |

Fitted orders on 10⁻² to 10⁻⁴: forward 1.00, central 2.00; five-point 4.00 on 10⁻² to 10⁻³ (it has already hit rounding by 10⁻⁴). On the right of the minimum every curve rises as roughly 1/h, and at h = 10⁻¹² all four formulas give the same wrong answer: there the subtraction has lost all 16 digits and the formula no longer matters.

## 5. Validate and improve

Put prediction and measurement side by side at α = 8°:

| quantity | predicted (§2.6–2.7) | measured | note |
|---|---|---|---|
| forward error at h = 10⁻² | (h/2)·21.66 = 0.108 | 0.108 | first order |
| central error at h = 10⁻² | (h²/6)·12.85 = 2.14 × 10⁻⁴ | 2.16 × 10⁻⁴ | the 0.7 % gap is the h⁴ term |
| five-point ratio 10⁻²→10⁻³ | 10⁴ | 9.9 × 10³ | fourth order |
| h_opt, central | 3.4 × 10⁻⁶ | minimum at 10⁻⁶ | within a factor of 3.4 |
| minimum central error | ~ε f/h_opt ≈ 5 × 10⁻¹¹ | 6.3 × 10⁻¹¹ | rounding floor |

Every difference has a cause you can name. The 0.7 % gap at h = 10⁻² in the central error is the next term in (4), (h⁴/120)f⁽⁵⁾ with f⁽⁵⁾ = −16 940 rad⁻⁵: it is −1.4 × 10⁻⁶, the same sign as the h² term, and adding it gives 2.156 × 10⁻⁴, the measured value to four figures. The measured minimum sits one decade to the *right* of h_opt in the table only because the table is sampled per decade; the formula for h_opt is an order-of-magnitude estimate, and that is all it claims. Rounding errors are not reproducible to better than a factor of two or three between machines and library versions, so the right-hand slope of −1 is a trend, not a measurement to three figures.

What the agreement does *not* show: that the model (2) describes a real wing. A measured lift curve from a NACA report (`P00.11`) has a slope near 0.1 per degree at small α for most sections, a C_Lmax between 1.2 and 1.8 depending on section and Reynolds number, and a stall far sharper than tanh. The numbers in this lesson validate the *derivative machinery*; validating the lift model is Phase 10's job.

Improvements worth making: for data on a grid, `np.gradient` implements the central difference with one-sided ends; for noisy data, fit a low-order polynomial over a window and differentiate the fit (Savitzky–Golay, in `scipy.signal`). Improvement worth *not* making: pushing h below 10⁻⁶ in the hope of more digits.

## 6. Defend your work

Submit: the completed `main.py` with all tests passing; the hand derivation of (3) with the value at 10°; `derivative_error_vs_h.png` with its CSV and the two fitted slopes per formula; the gust table from `main()`; and a one-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) saying what the lift model assumes, which derivative formula you would use for a wind-tunnel table and why, and where the linear gust estimate stops being safe.

### Adversarial self-checks

- **What if C_Lmax doubles?** The slope at α = 0 is unchanged (a₀), but the curve stays straight twice as far: at 10° the slope becomes 5.41 rad⁻¹ instead of 3.59. Check it with `lift_curve_slope_exact(np.radians(10), cl_max=2.8)`.
- **What if h goes to zero?** The formula divides by zero and returns `nan` or `inf`; your code raises `ValueError` instead. What if h goes to 1 rad? The forward difference of the lift model at 8° returns 0.62 rad⁻¹: not wrong in any floating-point sense, just the slope of a chord to a point 57° past the start of the plateau.
- **Which assumption fails first in the real world?** Smoothness of the lift curve. A real section stalls with a slope discontinuity, hysteresis between increasing and decreasing α, and Reynolds-number dependence; a finite difference across the stall returns an average of two regimes.
- **What measurement could falsify the model?** A measured C_L(α) whose slope at small α differs from 0.110 per degree by more than the thickness correction (thin-airfoil theory underestimates by a few per cent for a 12 % thick section) or whose plateau is not reached smoothly.
- **Does the result survive a change in step?** That is the entire content of part 4: it survives between 10⁻⁴ and 10⁻⁷ for the central difference and nowhere else to better than 10⁻⁸.
- **Which decision would change if the uncertainty doubled?** If the lift-curve slope were uncertain by 20 % instead of 10 %, the gust load factor would be uncertain by 0.3 g instead of 0.15 g at 1.44 g: the structural margin needed for a 5 m/s gust would double, and the certification gust speed would be met with a heavier wing.

## Check yourself

1. State the limit definition of the derivative and name the one thing about the point x = a that it does *not* use.
2. The ISA temperature is continuous at 11 km. Is it differentiable there? What are the two one-sided derivatives in K/m?
3. Differentiate C_L = C_Lmax tanh(a₀α/C_Lmax) and give the slope at α = 0 in per radian and per degree.
4. Why does the central difference have no h f″ term while the forward difference does?
5. A central-difference error falls from 2.2 × 10⁻⁴ at h = 10⁻² to 2.1 × 10⁻⁶ at h = 10⁻³. What is the observed order, and what is the error expected at h = 10⁻⁴?
6. Why does a derivative computed with h = 10⁻¹² come out worse than one with h = 10⁻⁶, and what constant sets the boundary?
7. Richardson extrapolation of the central difference with steps h and h/2 equals which other formula at which step?
8. A 5 m/s gust at 50 m/s changes α by about 0.1 rad. Why is the linear lift estimate 7 % high from α = 0 but 32 % high from α = 8°?
