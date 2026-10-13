# P01.08  Multivariable Functions & Partial Derivatives

> Hold everything else still, nudge one thing, and watch. Then add the nudges.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the partial derivatives of the lift of a finite wing with respect to angle of attack and airspeed, the total differential, the implicit derivative along a level curve and the chain-rule time derivative along a manoeuvre, and bound the error of the tangent-plane approximation. 2. Implement and test a lift-sensitivity module whose partials are computed three ways (analytic, central differences, SymPy) and agree. 3. Measure, by computational experiment, how the finite-difference error depends on the step and why one of the two partials shows no truncation error at all. 4. Decide, for a given gust, whether the linear sensitivity estimate is good enough or the full function must be evaluated. |
| **Prerequisites** | `P01.01`–`P01.02` (derivatives, linear approximation), `P01.05` (Taylor remainder), `P00.09` (SymPy `diff`, `lambdify`), `P00.03` (broadcasting). Skill check: you can differentiate $x^2 y^3$ with respect to $x$ treating $y$ as a constant, and you know that a central difference has error $\propto h^2$. |
| **Engineering question** | By how much does the lift of a light-aircraft wing change when a gust adds 5° of angle of attack or 10 m/s of airspeed, and can a first-order sensitivity estimate be trusted for that, or must the full lift function be evaluated? |
| **Assumptions and limits** | Incompressible flow (Mach below about 0.3); fixed wing geometry; lift coefficient depends on angle of attack only, through a lift-curve slope $a$ and a cubic roll-off toward stall; density fixed unless stated; quasi-steady (no unsteady aerodynamics). The cubic $C_L(\alpha)$ is a crude stand-in for stall: it has the right qualitative shape (a maximum, then a drop) but no claim to accuracy beyond about 12° for this wing. The model fails past stall, above Mach 0.3, and for a gust that arrives faster than the flow can adjust (a few chord-lengths of travel). |
| **Analytical work** | Partials of $L = \tfrac12\rho V^2 S\,C_L(\alpha)$; mixed partials and Clairaut's theorem; total differential and tangent plane with its second-order remainder; implicit differentiation of the level curve $L = W$; the chain rule for $L(\alpha(t), V(t))$; relative (logarithmic) sensitivities; a hand calculation at 4° and 55 m/s. |
| **Computational deliverable** | `starter/main.py`: `lift_coefficient`, `lift`, `lift_partials`, `partial_derivative`, `mixed_partial`, `lift_partials_sympy`, `lift_linearized`, `relative_sensitivities`, `level_curve_speed`, `level_curve_slope`, `lift_rate`, `stall_angle`, `tangent_plane_error_study`. Tests prove: $L \propto V^2$ and broadcasting; the linear model as the $k_3 = 0$ limit with sensitivities exactly 1 and 2; analytic partials equal central differences at second order and SymPy to $10^{-12}$; mixed partials agree (Clairaut); the tangent-plane error falls by 4 when the box halves and is exactly $\tfrac12\rho S C_L\,\Delta V^2$ along $V$; the level curve carries constant lift and its implicit slope matches a finite difference; the chain-rule rate matches a finite difference along a manoeuvre and vanishes along the level curve. |
| **Practical activity** | **computational experiment**: step-size study of the finite-difference partials ($h$ from $10^{-1}$ to $10^{-8}$) and a map of the tangent-plane error over gusts of $\pm 6°$ and $\pm 15$ m/s. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `lift_sensitivity.png`, the hand calculation at the reference point with its comparison to the code, the step-size table with both slopes, and the gust table with a decision on each row. |
| **Assessment** | `assessments/phase-01/grading-rubric.md`; competencies 1–4. Pass: all partials and the implicit derivative derived with the held-fixed variable named at each step; all public tests pass; the hand calculation agrees with the code to four figures; the gust decision is argued from the second-order remainder with a number. |
| **References** | `@openstax-calculus-volume-3` §4.1 "Functions of Several Variables" (required-free); §4.2 "Limits and Continuity" (supplementary-free, skim); §4.3 "Partial Derivatives" (required-free); §4.4 "Tangent Planes and Linear Approximations" (required-free); §4.5 "The Chain Rule" (required-free). |
| **Extensions** | ★ Add density as a third variable, $L(\alpha, V, \rho)$, with $\rho(h)$ from `P00.03`, and compute $\partial L/\partial h$ at 2 km by the chain rule. ★ Replace the cubic $C_L$ by a table interpolated from `@ntrs-naca-tr-824-summary-of-airfoil-data` and see what happens to the finite-difference partials at the table's knots. |
| **Time** | unmeasured |

## 1. The engineering problem

A four-seat light aircraft of 1540 kg cruises at 55 m/s with its wing at 4° angle of attack, producing 15.1 kN of lift to balance its weight. It flies into a gust. A vertical gust changes the angle of attack; a horizontal gust changes the airspeed; a real gust does both. The structures engineer needs the resulting lift increment to size the wing spar; the flight-dynamics engineer needs it to predict the bump the passengers feel. Both want a *sensitivity*: how many newtons per degree, how many newtons per metre per second. Both also need to know when that sensitivity is enough and when the gust is large enough that the lift must be recomputed from scratch.

Lift is a function of two variables, $L(\alpha, V)$, and a sensitivity to one of them while the other is held fixed is a **partial derivative**. The two partials together give the **total differential**, which is the linear estimate of the lift change for any small combined gust. The questions this lesson answers are: how to compute partials (three ways, so they can be checked against each other), how to combine them, how big the error of the linear estimate is, and how the same tools give the trim speed at every angle of attack (a **level curve** of $L$) and the rate of change of lift during a pull-up (the **chain rule**).

Success criteria, written before we start: the three ways of computing $\partial L/\partial\alpha$ and $\partial L/\partial V$ agree to $10^{-7}$ relative at the reference point; the error of the linear estimate over a box of gusts falls by a factor of 4 each time the box is halved; the lift along the computed trim curve is constant to $10^{-12}$; and the hand calculation agrees with the code to four figures.

## 2. Understand and derive

### 2.1 The lift function

Lift on a wing of area $S$ in air of density $\rho$ at speed $V$ is

$$L(\alpha, V) = \tfrac12\rho V^2 S\,C_L(\alpha) \quad [\mathrm{N}],$$   (1)

with $C_L$ the dimensionless **lift coefficient**, a property of the wing's shape and orientation. Units: (kg/m³)(m²/s²)(m²) = kg m/s² = N ✓. **Assumption 1 (incompressible, quasi-steady):** $C_L$ depends on the angle of attack $\alpha$ only. **Assumption 2 (model for $C_L$):** below stall a wing's lift coefficient is linear in $\alpha$ with **lift-curve slope** $a$ (about $2\pi$ per radian for a thin airfoil, reduced to $a \approx 5$ per radian for a finite wing of aspect ratio 7 by the downwash of its own tip vortices; Phase 10 derives both), measured from the **zero-lift angle** $\alpha_0$ (negative for a cambered wing). Near stall the lift rolls off. We model the whole curve as

$$C_L(\alpha) = a(\alpha - \alpha_0) - k_3(\alpha - \alpha_0)^3, \qquad \alpha \text{ in radians},$$   (2)

whose maximum (stall) is where $\mathrm{d}C_L/\mathrm{d}\alpha = a - 3k_3(\alpha - \alpha_0)^2 = 0$, at $\alpha_s = \alpha_0 + \sqrt{a/3k_3}$. The numbers for our wing: $S = 16$ m², $a = 5.0$ /rad, $\alpha_0 = -2°$, $k_3 = 12$ /rad³, $\rho = 1.225$ kg/m³. Then $\alpha_s = 19.4°$ and $C_{L,\max} = 1.24$, a little high in $\alpha$ for a real light-aircraft wing (which stalls nearer 16°) but with the right shape; §2.9 says what that costs.

**Hand calculation at the reference point** $\alpha = 4°$, $V = 55$ m/s: $\alpha - \alpha_0 = 6° = 0.10472$ rad; $C_L = 5.0 \times 0.10472 - 12 \times 0.10472^3 = 0.52360 - 0.01378 = 0.50982$; dynamic pressure $q = \tfrac12 \times 1.225 \times 55^2 = 1852.8$ Pa; $L = 1852.8 \times 16 \times 0.50982 = 15\,114$ N, which is 1541 kg of aircraft. Sanity: a Cessna-class four-seater at 55 m/s (107 kt) and 4° is right.

A function of two variables is a surface over the $(\alpha, V)$ plane (`@openstax-calculus-volume-3` §4.1). The curves of constant $L$ on that plane are its **level curves**; the one through $L = W$ (weight) is the **trim curve**, every combination of angle of attack and speed at which the wing carries exactly the aircraft. Solving (1) for $V$,

$$V(\alpha) = \sqrt{\frac{2L}{\rho S\,C_L(\alpha)}},$$   (3)

defined only where $C_L > 0$, and falling as $\alpha$ rises: more angle, less speed. At $C_{L,\max}$ it bottoms out at the **stall speed**, $35.2$ m/s for our wing at 15.1 kN.

### 2.2 Partial derivatives

The **partial derivative** of $L$ with respect to $\alpha$ is the ordinary derivative of $L$ along a line of constant $V$ (`@openstax-calculus-volume-3` §4.3):

$$\frac{\partial L}{\partial\alpha} = \lim_{h \to 0}\frac{L(\alpha + h, V) - L(\alpha - h, V)}{2h}.$$

Everything that does not contain $\alpha$ is a constant, so from (1) and (2),

$$L_\alpha \equiv \frac{\partial L}{\partial\alpha} = \tfrac12\rho V^2 S\,\big[a - 3k_3(\alpha - \alpha_0)^2\big] \quad [\mathrm{N/rad}],$$   (4)

$$L_V \equiv \frac{\partial L}{\partial V} = \rho V S\,C_L(\alpha) = \frac{2L}{V} \quad [\mathrm{N/(m/s)}].$$   (5)

(5) is worth a pause: because $L \propto V^2$, $L_V = 2L/V$ whatever $C_L$ is, so a 1 % increase in speed is always a 2 % increase in lift. The analogous **relative sensitivity** in $\alpha$ is $L_\alpha(\alpha - \alpha_0)/L = [a - 3k_3 d^2]/[a - k_3 d^2]$ with $d = \alpha - \alpha_0$: exactly 1 for the linear model, less as stall approaches, and zero at $C_{L,\max}$, where more angle of attack buys no lift at all.

At the reference point: $\mathrm{d}C_L/\mathrm{d}\alpha = 5.0 - 36 \times 0.010966 = 4.6052$ /rad, so $L_\alpha = 1852.8 \times 16 \times 4.6052 = 136\,522$ N/rad $= 2383$ N/deg, and $L_V = 2 \times 15\,114/55 = 549.6$ N/(m/s). Relative sensitivities: 0.946 in $\alpha$ (the cubic term is already costing 5 %) and 2.000 in $V$.

### 2.3 Second partials and Clairaut's theorem

Differentiate again, each way:

$$L_{\alpha\alpha} = \tfrac12\rho V^2 S\,(-6k_3 d), \qquad L_{VV} = \rho S\,C_L, \qquad L_{\alpha V} = \rho V S\,[a - 3k_3 d^2] = L_{V\alpha}.$$   (6)

The two **mixed partials** are equal: differentiating (4) in $V$ or (5) in $\alpha$ gives the same expression. That is Clairaut's theorem, true whenever the second partials are continuous (`@openstax-calculus-volume-3` §4.3), and it is a test: compute $L_{\alpha V}$ by finite differences in the two orders and they must agree. At the reference point $L_{\alpha\alpha} = -223\,518$ N/rad², $L_{VV} = 9.992$ N/(m/s)², $L_{\alpha V} = 4964$ N/(rad m/s). The sign of $L_{\alpha\alpha}$ is the lift curve bending over toward stall.

### 2.4 The total differential and the tangent plane

For a small combined change $(\Delta\alpha, \Delta V)$ the first-order estimate adds the two one-variable estimates (`@openstax-calculus-volume-3` §4.4):

$$\Delta L \approx \mathrm{d}L = L_\alpha\,\Delta\alpha + L_V\,\Delta V.$$   (7)

Geometrically, (7) is the **tangent plane** to the surface $L(\alpha, V)$ at the reference point, and the surface itself lies off the plane by the second-order remainder (two-variable Taylor's theorem, `P01.05` extended):

$$\Delta L - \mathrm{d}L = \tfrac12 L_{\alpha\alpha}\Delta\alpha^2 + L_{\alpha V}\Delta\alpha\,\Delta V + \tfrac12 L_{VV}\Delta V^2 + O(\Delta^3).$$   (8)

Halve both $\Delta$s and the remainder drops by 4: that is the convergence test. For a pure speed gust the remainder is *exactly* $\tfrac12 L_{VV}\Delta V^2 = \tfrac12\rho S C_L\Delta V^2$, with no $O(\Delta^3)$, because $L$ is a quadratic in $V$; for a pure $\alpha$ gust there is also a cubic term, $\tfrac16 L_{\alpha\alpha\alpha}\Delta\alpha^3$ with $L_{\alpha\alpha\alpha} = -6k_3 qS$.

**The gust, by hand.** $\Delta\alpha = +5° = 0.08727$ rad at fixed $V$: linear $\mathrm{d}L = 136\,522 \times 0.08727 = 11\,914$ N; second-order correction $\tfrac12(-223\,518)(0.08727)^2 = -851$ N; third-order $\tfrac16(-6 \times 12 \times 1852.8 \times 16)(0.08727)^3 = -236$ N; total $10\,827$ N against the exact $L(9°, 55) - L_0 = 10\,826$ N. The linear estimate is 10 % high. $\Delta V = +10$ m/s at fixed $\alpha$: linear $5496$ N, remainder $\tfrac12 \times 9.992 \times 100 = 500$ N, exact $5995$ N: the linear estimate is 8 % low. Both together: linear $17\,410$ N; the remainder (8) adds $-851 + 4964 \times 0.08727 \times 10 + 500 = 3981$ N plus the cubic $-236$ N, giving $21\,155$ N against the exact $21\,116$ N. The linear estimate is now 18 % low, and the cross term $L_{\alpha V}\Delta\alpha\Delta V$ is the biggest correction: lift from extra angle *at the higher speed* is worth more than either alone. That cross term is exactly what a one-variable analysis cannot see.

### 2.5 Implicit differentiation: the slope of the trim curve

Along the level curve $L(\alpha, V) = W$ the two variables are tied together. Differentiate the identity $L(\alpha, V(\alpha)) = W$ with the chain rule (`@openstax-calculus-volume-3` §4.5):

$$L_\alpha + L_V\frac{\mathrm{d}V}{\mathrm{d}\alpha} = 0 \quad\Longrightarrow\quad \frac{\mathrm{d}V}{\mathrm{d}\alpha} = -\frac{L_\alpha}{L_V}.$$   (9)

For our wing $-136\,522/549.6 = -248.4$ (m/s)/rad $= -4.34$ (m/s)/deg: each extra degree of angle of attack lets the aircraft fly 4.3 m/s slower at the same lift. The same number from (3) directly: $\mathrm{d}V/\mathrm{d}\alpha = -\tfrac12 V\,(\mathrm{d}C_L/\mathrm{d}\alpha)/C_L = -\tfrac12 \times 55 \times 4.6052/0.50982 = -248.4$ ✓. Two routes, one number, which is the test. Implicit differentiation is how you get a slope when you cannot, or would rather not, solve for $V(\alpha)$ explicitly.

### 2.6 The chain rule along a manoeuvre

In a pull-up both $\alpha$ and $V$ change with time. The rate of change of lift is

$$\frac{\mathrm{d}L}{\mathrm{d}t} = \frac{\partial L}{\partial\alpha}\frac{\mathrm{d}\alpha}{\mathrm{d}t} + \frac{\partial L}{\partial V}\frac{\mathrm{d}V}{\mathrm{d}t},$$   (10)

one term per path by which time reaches $L$. For $\dot\alpha = 2°$/s and $\dot V = -1$ m/s² at the reference point: $136\,522 \times 0.034907 - 549.6 \times 1 = 4766 - 550 = 4216$ N/s, so the load factor rises by $4216/15\,114 = 0.28\,g$ per second. Units: (N/rad)(rad/s) + (N/(m/s))(m/s²) = N/s ✓. Along the trim curve, $\dot V = (\mathrm{d}V/\mathrm{d}\alpha)\dot\alpha$ and (10) gives zero by (9): the chain rule and implicit differentiation are the same statement.

### 2.7 Finite differences in several variables

A partial derivative is a one-variable derivative, so the central difference of `P00.09` applies with only one coordinate perturbed:

$$\frac{\partial f}{\partial x_i} \approx \frac{f(\mathbf x + h\mathbf e_i) - f(\mathbf x - h\mathbf e_i)}{2h},$$   (11)

truncation error $\tfrac16 h^2\,\partial^3 f/\partial x_i^3$, rounding error about $\varepsilon|f|/h$. For $\partial L/\partial\alpha$, $L_{\alpha\alpha\alpha} = -6k_3 qS = -2.134 \times 10^6$ N/rad³, so at $h = 10^{-2}$ rad the truncation error is $\tfrac16 \times 10^{-4} \times 2.134 \times 10^6 = 35.6$ N/rad (relative $2.6 \times 10^{-4}$), and it falls 100-fold per decade of $h$ until rounding takes over near $h = 10^{-6}$. For $\partial L/\partial V$ the third derivative is zero ($L$ is quadratic in $V$) and the central difference is exact up to rounding at every $h$: the two columns of the step-size table look completely different, and part 4 asks you to say why before you run it.

### 2.8 How SymPy does it

`sp.diff(L, alpha)` treats every other symbol as a constant, which is precisely the definition of a partial derivative; the result is (4) in symbolic form, and substituting numbers gives the same $136\,521.6$ N/rad as the hand formula. The value of the symbolic route is not the number but the *structure*: `sp.simplify(sp.diff(L, V) - 2*L/V)` returns 0, proving (5) for every $C_L$, and `sp.simplify(sp.diff(L, alpha, V) - sp.diff(L, V, alpha))` returns 0, proving Clairaut for this $L$ without choosing a point.

### 2.9 Where the model stops being valid

- **Stall.** The cubic (2) is a caricature. A real $C_L(\alpha)$ has a sharper peak, hysteresis, and depends on Reynolds number. Partials computed past about 12° for this wing are properties of the caricature, not the aircraft.
- **Compressibility.** Above Mach 0.3 (about 100 m/s at sea level) $C_L$ grows with Mach by roughly $1/\sqrt{1 - M^2}$, adding a $V$-dependence to $C_L$ that (5) ignores.
- **Unsteady flow.** A sharp-edged gust takes several chord-lengths of travel (about 0.1 s here) before the lift reaches its new steady value; (7) is the *steady* increment.
- **Density.** $L_\rho = L/\rho$ is the third partial we have frozen; at 2000 m the ISA density is 1.006 kg/m³ and every lift number drops 18 %. The extension adds it.

## 3. Build the computational model

Open `starter/main.py`. The wing lives in a frozen `Wing` dataclass so every function takes the same parameter object; the tests build their own `Wing(k3=0.0)` to reach the linear model.

### Step 1: the function and its analytic partials

```python
def lift_coefficient(alpha, wing=DEFAULT_WING):
    d = np.asarray(alpha, dtype=float) - wing.alpha0
    return wing.a * d - wing.k3 * d**3                                          # (2)

def lift(alpha, V, wing=DEFAULT_WING):
    V = np.asarray(V, dtype=float)
    if np.any(V < 0.0):
        raise ValueError("airspeed must be non-negative")
    return 0.5 * wing.rho * V**2 * wing.S * lift_coefficient(alpha, wing)     # (1)

def lift_partials(alpha, V, wing=DEFAULT_WING):
    d = alpha - wing.alpha0
    dcl = wing.a - 3.0 * wing.k3 * d**2
    return 0.5 * wing.rho * V**2 * wing.S * dcl, wing.rho * V * wing.S * lift_coefficient(alpha, wing)   # (4), (5)
```

`alpha[:, None]` against `V[None, :]` broadcasts to a full grid, which is how `main()` draws the level curves.

### Step 2: finite-difference partials of any function

```python
def partial_derivative(f, x, i, h):
    e = np.zeros_like(x); e[i] = h                     # perturb one coordinate only
    return float((f(x + e) - f(x - e)) / (2.0 * h))   # (11)

def mixed_partial(f, x, i, j, h):
    return partial_derivative(lambda y: partial_derivative(f, y, i, h), x, j, h)
```

`f` takes a single array argument, so the same code differentiates $L(\alpha, V)$, $x^2 y$ or anything else; the tests use all three.

### Step 3: SymPy, linearization, sensitivities

`lift_partials_sympy` builds $L$ from symbols and returns `sp.diff(L, alpha)`, `sp.diff(L, V)` and the symbol table. `lift_linearized` is (7) added to $L$ at the reference point; `relative_sensitivities` returns $L_\alpha(\alpha - \alpha_0)/L$ and $L_V V/L$; `stall_angle` is $\alpha_0 + \sqrt{a/3k_3}$ and raises for $k_3 = 0$; `tangent_plane_error_study` evaluates $\max|L - L_{\text{lin}}|$ over an $11 \times 11$ grid in a box that shrinks with each scale.

### Step 4: level curve and chain rule

```python
def level_curve_speed(alpha, L_target, wing=DEFAULT_WING):
    cl = lift_coefficient(alpha, wing)
    if np.any(cl <= 0.0):
        raise ValueError("C_L must be positive on the level curve")
    return np.sqrt(2.0 * L_target / (wing.rho * wing.S * cl))                  # (3)

def level_curve_slope(alpha, V, wing=DEFAULT_WING):
    la, lv = lift_partials(alpha, V, wing)
    return -la / lv                                                             # (9)

def lift_rate(alpha, V, dalpha_dt, dV_dt, wing=DEFAULT_WING):
    la, lv = lift_partials(alpha, V, wing)
    return la * dalpha_dt + lv * dV_dt                                          # (10)
```

### Verify

Run `python -m pytest curriculum/phase-01/lesson-08-partial-derivatives` from the repo root. The tests, in order of how much they would hurt if wrong:

1. **The function.** $L(2V) = 4L(V)$; $L = 0$ at $\alpha_0$ and at $V = 0$; a $(5, 1) \times (1, 3)$ input gives a $(5, 3)$ output; negative speed raises.
2. **Three partials agree.** Analytic against central differences at second order in $\alpha$ (rates above 1.9 as $h$ halves) and to $10^{-10}$ in $V$ at any $h$; analytic against SymPy to $10^{-12}$; SymPy proves $L_V = 2L/V$ and Clairaut symbolically; the finite-difference mixed partials agree in both orders and with (6).
3. **Linearization.** The tangent-plane error over a shrinking box falls by $4.0 \pm 0.2$ per halving; along $V$ the error is exactly $\tfrac12\rho S C_L\Delta V^2$.
4. **Limits.** With $k_3 = 0$ the relative sensitivities are exactly 1 and 2 and `stall_angle` raises; at $\alpha_s$ the $\alpha$-sensitivity is zero and `argmax` of $C_L$ on a fine grid lands within 0.02° of it.
5. **Level curve and chain rule.** Lift along `level_curve_speed` is constant to $10^{-12}$ and the curve descends; the implicit slope matches a central difference of (3) to $10^{-7}$; `lift_rate` matches a finite difference of $L(\alpha(t), V(t))$ along a quadratic manoeuvre to $10^{-7}$ and vanishes along the trim curve.
6. **Interface.** `main()` writes the figure.

### Compare

`main()` prints the three routes side by side:

```
dL/dalpha = 136521.6 N/rad = 2382.8 N/deg;  dL/dV = 549.6 N/(m/s)
central differences h = 1e-02: dL/dalpha err 3.56e+01 N/rad, dL/dV err 2.17e-10 N/(m/s)
central differences h = 1e-03: dL/dalpha err 3.56e-01 N/rad, dL/dV err 1.22e-09 N/(m/s)
central differences h = 1e-04: dL/dalpha err 3.56e-03 N/rad, dL/dV err 1.70e-08 N/(m/s)
sympy: dL/dalpha = 136521.6, dL/dV = 549.6
mixed partials: L_alpha_V = 4964.410, L_V_alpha = 4964.410 N/(rad m/s)
```

The $\alpha$ column falls exactly 100-fold per decade, and $35.6$ N/rad at $h = 10^{-2}$ is the predicted $\tfrac16 h^2|L_{\alpha\alpha\alpha}|$ to three figures: the finite-difference error is not noise, it is a number you can compute in advance. The $V$ column is rounding only, growing as $h$ shrinks, as §2.7 said. SymPy's `diff` and the hand formula agree to all printed digits because they are the same algebra. The one established library in play is SymPy; there is no NumPy "partial derivative" to call, which is why `partial_derivative` is yours and why its error must be understood rather than trusted.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.**

**Experiment A: step-size study.** For $h = 10^{-1}, 10^{-2}, \dots, 10^{-8}$ (radians for $\alpha$, m/s for $V$), record `partial_derivative(f, x0, 0, h)` and `partial_derivative(f, x0, 1, h)` against the analytic values at the reference point. Prediction, written before running: the $\alpha$ error follows $\tfrac16 h^2|L_{\alpha\alpha\alpha}| = 3.56 \times 10^5 h^2$ N/rad down to a minimum near $h_{\text{opt}} = (3\varepsilon|L|/|L_{\alpha\alpha\alpha}|)^{1/3} \approx 2 \times 10^{-6}$ and then rises as $\varepsilon|L|/h$; the $V$ error has no $h^2$ branch and is at the rounding floor ($\sim 10^{-10}$) for every $h$ above $10^{-5}$. Plot both on log–log axes, fit the slopes on each side of the minimum, and save as `fd_partials_vs_h.png`.

**Experiment B: map of the tangent-plane error.** On a grid $\Delta\alpha \in [-6°, 6°]$, $\Delta V \in [-15, 15]$ m/s around the reference point, compute $L - L_{\text{lin}}$ and the second-order prediction (8). Prediction: the error is dominated by the cross term, so it is positive when $\Delta\alpha$ and $\Delta V$ have the same sign and negative otherwise, and (8) reproduces it to within the cubic term $\tfrac16 L_{\alpha\alpha\alpha}\Delta\alpha^3$ ($\le 410$ N at $6°$). Contour the error and mark the region where $|L - L_{\text{lin}}| < 2$ % of $L_0$: that is the region in which the sensitivity estimate is good enough.

Raw data format: two tables, $(h, \text{err}_\alpha, \text{err}_V)$ and $(\Delta\alpha, \Delta V, \Delta L_{\text{exact}}, \mathrm{d}L, \text{remainder})$, kept in your lesson notes.

## 5. Validate and improve

Hand (§2), code and the independent route side by side at the reference point:

| quantity | hand | code (analytic) | independent check | difference |
|---|---|---|---|---|
| $L$ (N) | 15 114 | 15 113.6 | | 0 |
| $L_\alpha$ (N/rad) | 136 522 | 136 521.6 | FD $h = 10^{-4}$: 136 521.62; SymPy 136 521.6 | $3 \times 10^{-8}$ |
| $L_V$ (N/(m/s)) | 549.6 | 549.58 | $2L/V = 549.58$; FD exact | $10^{-12}$ |
| $L_{\alpha V}$ (N/(rad m/s)) | 4964 | 4964.42 | FD both orders 4964.41 | $2.6 \times 10^{-6}$ |
| $\mathrm{d}V/\mathrm{d}\alpha$ ((m/s)/rad) | $-248.4$ | $-248.41$ | FD of (3): $-248.41$ | $10^{-8}$ |
| gust $+5°$: $\Delta L$ (N) | 10 827 (to 3rd order) | 10 826 (exact) | linear 11 914 | linear 10 % high |
| gust $+10$ m/s: $\Delta L$ (N) | 5995 (to 2nd order, exact) | 5995 | linear 5496 | linear 8 % low |
| both: $\Delta L$ (N) | 21 155 (to 3rd order) | 21 116 | linear 17 410 | linear 18 % low |

Every residual has a cause. The $10^{-8}$ on $L_\alpha$ is the $h^2$ truncation at $h = 10^{-4}$, predicted by §2.7. The $2.6 \times 10^{-6}$ on the mixed partial is the same truncation applied twice (the inner difference in $\alpha$ at $h = 10^{-3}$). The 39 N gap between the third-order hand estimate and the exact combined gust is the fourth-order cross terms, $\tfrac12 L_{\alpha\alpha V}\Delta\alpha^2\Delta V$ and friends, which (8) drops: $L$ is cubic in $\alpha$ and quadratic in $V$, so its Taylor series terminates at fifth order and you could, in principle, write the gust response exactly as a polynomial. The 10 %, 8 % and 18 % errors of the linear estimate are the engineering content: for a $1°$ or $1$ m/s gust the linear estimate is within 0.3 %, for a $5°$ and $10$ m/s gust it is wrong by a fifth, and the cross term is the main reason.

What is still unverified when everything agrees: the lift model. The three routes to $L_\alpha$ agree because they differentiate the same (2); a wind-tunnel $C_L(\alpha)$ for this wing at this Reynolds number would differ from the cubic by several percent above 10°, and *that* uncertainty, not the $10^{-8}$ of the arithmetic, is the precision of any gust load computed here. Phase 10 replaces (2) with a computed lift curve; `@ntrs-naca-tr-824-summary-of-airfoil-data` has measured ones.

Improvements worth making now: use the second-order estimate (8) rather than (7) for gusts above about 2°, since the Hessian costs three more numbers and removes most of the error; or, simplest of all, evaluate $L$ at the gust point directly, which is cheap for a formula but is not an option when "the function" is a two-hour CFD run, which is when sensitivities earn their keep. Improvements worth *not* making: pushing $h$ below $10^{-6}$ to "get more digits"; part 4 shows what happens.

## 6. Defend your work

Submit: the completed `starter/main.py` with all tests passing; `lift_sensitivity.png` and `fd_partials_vs_h.png`; the hand calculation with its comparison table; the two experiment tables with one paragraph each; and a one-page note (sections 1, 2, 5, 6 and 7 of `projects/_template/report-template.md`) stating for which gust sizes the linear sensitivity is adequate and why.

### Adversarial self-checks

- **What if the speed doubles?** $L$ quadruples, $L_\alpha$ quadruples, $L_V$ doubles, and the relative sensitivities do not change at all. Check with the code; if $L_V$ quadruples you have differentiated $V^2$ as $V^2$.
- **What if $V \to 0$? What if $\alpha \to \alpha_s$?** At $V = 0$ both partials vanish and the level curve (3) runs off to infinite speed as $C_L \to 0$; at $\alpha_s$, $L_\alpha = 0$, the implicit slope (9) is zero (the trim curve is horizontal: stall speed), and the relative $\alpha$-sensitivity is zero.
- **Which assumption fails first in the real world?** The cubic stall model, as soon as a gust takes $\alpha$ past about 12°; then quasi-steadiness for gusts shorter than a tenth of a second; compressibility only above 100 m/s.
- **What measurement could falsify the model?** A wind-tunnel lift curve: if the measured $\mathrm{d}C_L/\mathrm{d}\alpha$ at 4° is not within a few percent of 4.6 /rad, every sensitivity here is off by the same factor. A flight test with an accelerometer through a known gust gives $\Delta L$ directly.
- **Does the result survive a change in step size?** $L_\alpha$ by finite difference must agree to $10^{-6}$ for any $h$ between $10^{-5}$ and $10^{-3}$. If your answer depends on $h$ in that range, the function you are differentiating is not smooth (a table with kinks, a clipped $C_L$) and no finite-difference partial is meaningful there.
- **Which decision changes if the uncertainty doubled?** Suppose $a$ is known to ±10 % instead of ±5 %. The $+5°$ gust load is then $10\,826 \pm 1083$ N instead of ±541 N; a spar designed with a 1.5 safety factor on 15.1 kN has 7.5 kN of margin either way, so the structure does not change, but a stall-warning margin set at 2° above the predicted $\alpha$ for a given speed now has to be 3°. The decision that moves is the one with the smaller margin.

## Check yourself

1. Write $\partial L/\partial V$ for $L = \tfrac12\rho V^2 S C_L(\alpha)$ and say why it does not depend on the form of $C_L$.
2. What is held fixed when you compute $\partial L/\partial\alpha$, and what would be different about $\mathrm{d}L/\mathrm{d}\alpha$ along the trim curve?
3. State Clairaut's theorem for $L(\alpha, V)$ and the condition under which it holds.
4. The total differential $\mathrm{d}L = L_\alpha\Delta\alpha + L_V\Delta V$ misses one second-order term that the two single-variable estimates cannot see. Which one?
5. Along the level curve $L = W$, what is $\mathrm{d}V/\mathrm{d}\alpha$ in terms of the partials, and what is its sign below stall?
6. For a 1 % increase in airspeed, by what percentage does lift change? For a 1 % increase in $\alpha - \alpha_0$ in the linear model?
7. Why does the central difference of $\partial L/\partial V$ show no $h^2$ error while that of $\partial L/\partial\alpha$ does?
8. In a pull-up, $\alpha$ increases at 2°/s and $V$ decreases at 1 m/s². Write $\mathrm{d}L/\mathrm{d}t$ with its units and say which term dominates for the lesson's wing.
