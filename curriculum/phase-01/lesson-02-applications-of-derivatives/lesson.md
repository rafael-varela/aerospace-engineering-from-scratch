# P01.02  Applications of Derivatives

> The best you can do is where the derivative is zero; how much it matters is the second derivative; how fast it changes is the chain rule.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the best-glide and minimum-sink lift coefficients of a parabolic drag polar by setting a derivative to zero and confirm each with the second-derivative test. 2. Implement and test a glide-performance module: polar, L/D and its derivatives, Newton's method, linear sensitivity, related rates, and the best speed in a headwind. 3. Measure how the linearisation error of a sensitivity estimate grows with the size of the perturbation and show it is second order. 4. Decide, from the ground-range curve, how much faster to fly into a given headwind, and state when the parabolic-polar model stops deserving trust. |
| **Prerequisites** | `P01.01` (derivatives, chain rule, finite differences), `P00.03` (NumPy), `P00.06` (SciPy's `minimize_scalar`), `P00.07` (pytest). Skill check: you can apply the quotient rule to x/(a + bx²), solve a quadratic, and use √ and arctan on a calculator. |
| **Engineering question** | After an engine failure at 1000 m, how far can this aircraft glide, at what airspeed, how much faster should it fly into a 10 m/s headwind, and by how much does a 10 % increase in parasite drag shorten the answer? |
| **Assumptions and limits** | Steady, unaccelerated glide in still air with no thrust; the parabolic polar C_D = C_D0 + K C_L² with constant C_D0 and K (valid between about C_L = 0.2 and 90 % of C_Lmax; fails near stall, at high Mach, and for cambered sections whose minimum drag is not at C_L = 0); ISA troposphere density; angles small enough that cos γ ≈ 1 (here γ = 1.87°, cos γ = 0.99947). The model says nothing about the transition after the failure, the landing, or a wind that changes with height. |
| **Analytical work** | Force balance in a glide; the drag polar; maximisation of L/D by the quotient rule and the second-derivative test; l'Hôpital's rule at the ends of the polar; minimum sink; linear approximation of (L/D)max with respect to C_D0 and its second-order error; the related rate dV*/dt in the descent; the headwind optimum as the root of a transcendental equation; Newton's method and its quadratic convergence. |
| **Computational deliverable** | `starter/main.py`: `induced_drag_factor`, `drag_polar`, `lift_to_drag`, `ld_derivative`, `ld_second_derivative`, `best_glide_cl`, `max_lift_to_drag`, `min_sink_cl`, `glide_angle`, `airspeed`, `sink_rate`, `sink_rate_of_speed`, `glide_range`, `newton`, `best_glide_cl_newton`, `ld_max_sensitivity_cd0`, `linear_approximation`, `isa_density`, `isa_density_gradient`, `best_glide_speed`, `best_glide_speed_rate`, `ground_glide_ratio`, `headwind_best_speed`. Tests prove: the critical point is a maximum and beats every sampled C_L, induced drag equals parasite drag there, Newton reproduces it inside its basin and fails outside, derivatives match finite differences, the linearisation error is second order, the related rate matches a time difference, and the headwind optimum reduces to still-air best glide at zero wind. |
| **Practical activity** | **computational experiment**: a parameter sweep of (L/D)max against C_D0 and aspect ratio, of the best speed and ground range against headwind, and a measurement of the linearisation error against perturbation size. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `glide_polar.png`, the hand derivation of C_L* with the second-derivative test, the sweep tables, and a short note stating the speed-to-fly rule you extracted for headwinds and where the model fails. |
| **Assessment** | `assessments/phase-01/grading-rubric.md`; competencies 1–4. Pass: the derivation of (L/D)max is complete with the second-derivative test; all public tests pass; the hand values for the glider agree with the code to three figures; the headwind note gives a rule with a number and a stated range of validity. |
| **References** | `@openstax-calculus-volume-1` §3.4 "Derivatives as Rates of Change", §4.1 "Related Rates", §4.2 "Linear Approximations and Differentials", §4.3 "Maxima and Minima", §4.5 "Derivatives and the Shape of a Graph", §4.7 "Applied Optimization Problems", §4.8 "L'Hôpital's Rule", §4.9 "Newton's Method" (all required-free). `@openstax-calculus-volume-1` §4.4 "The Mean Value Theorem", §4.6 "Limits at Infinity and Asymptotes" (supplementary-free). |
| **Extensions** | ★ Add rising or sinking air (a vertical wind component) to the ground-range optimum and derive the "MacCready" speed-to-fly. ★ Replace the parabolic polar by a tabulated polar from a NACA report and redo the optimisation numerically; compare the optimum C_L. |
| **Time** | unmeasured |

## 1. The engineering problem

A 350 kg training glider (wing area 10.5 m², aspect ratio 20) is at 1000 m when it loses its lift source, or, the same problem in another aircraft, a light aeroplane's engine stops. From now on the only energy available is height. Every pilot is taught two numbers for their aircraft: the **best-glide speed**, which covers the most ground per metre of height lost, and the **minimum-sink speed**, which keeps the aircraft in the air for the longest time. Fly too fast and drag burns the height; fly too slowly and the wing works so hard that induced drag does the same. Somewhere in between is a best, and "best" is a word that means "a derivative is zero".

The engineering question is fourfold. How far can the glider go from 1000 m in still air? At what airspeed? How much faster should it fly into a 10 m/s headwind, where the still-air answer is wrong? And if the airframe is dirtier than the brochure says, with 10 % more parasite drag, how much range is lost? All four are derivative questions: an optimum (part 2.3), a sensitivity (2.6), a rate of change during the descent (2.7), and a root of an equation that has no closed form (2.8).

Success criteria, written down first: the best-glide lift coefficient found by Newton's method agrees with the closed form to 10⁻⁸; `scipy.optimize.minimize_scalar` finds the same optimum to 10⁻⁸; the linear sensitivity estimate is within 0.4 % of the exact value for a 10 % change in C_D0 and its error falls four-fold when the change halves; the headwind optimum at zero wind reproduces the still-air speed to 10⁻⁹; and the glide angle is small enough that the small-angle sink rate is within 0.2 % of the exact one.

## 2. Understand and derive

### 2.1 Forces in a steady glide

In a steady glide the aircraft moves at constant speed V along a straight path inclined at the **glide angle** γ below the horizontal. Three forces act: weight W down, lift L perpendicular to the path, drag D along the path, backwards.

```
                     L
                     ▲
                     │        path, angle γ below horizontal
            D ◄──────●──────────────────────────────▶ (direction of motion, down-slope)
                     │╲
                     │ ╲ γ
                     ▼  ╲
                     W   ╲
```

**Assumption 1 (steady):** no acceleration, so the forces balance along and across the path:

L = W cos γ,   D = W sin γ   ⟹   **tan γ = D/L = C_D/C_L**.   (1)

The glide angle depends only on the lift-to-drag ratio, not on weight: a heavier glider flies the same path faster. From L = ½ρV²S C_L = W cos γ,

**V = √(2 W cos γ / (ρ S C_L))**.   (2)

Units: √(N/((kg/m³) m²)) = √(kg m s⁻² m / kg) = m/s. Good. The **sink rate** is w = V sin γ (m/s, positive down), and the ground covered per metre of height lost is 1/tan γ = L/D.

**Assumption 2 (small angle):** when L/D is large, γ is small, sin γ ≈ tan γ ≈ γ and cos γ ≈ 1. For L/D = 30, γ = 1.87° and cos γ = 0.99947; the error is 0.05 %, and we use the small-angle forms freely for ranges while the code keeps the exact ones available.

### 2.2 The drag polar

**Assumption 3 (parabolic polar):** the drag coefficient is

**C_D = C_D0 + K C_L²**,   K = 1/(π e AR).   (3)

C_D0 is the **parasite drag** coefficient (skin friction and pressure drag at zero lift), K C_L² is the **induced drag** (the price of producing lift with a finite wing), AR = b²/S is the aspect ratio and e the Oswald efficiency factor (0.7–0.95 for real wings). Phase 10 derives (3); here it is the given model. For the glider: C_D0 = 0.015, e = 0.9, AR = 20, so K = 1/(π × 0.9 × 20) = 0.017684.

### 2.3 Maximum L/D: a critical point, confirmed

The glide ratio as a function of C_L is f(C_L) = C_L/(C_D0 + K C_L²). The **quotient rule** (P01.01) gives

f′(C_L) = [(C_D0 + K C_L²) − C_L · 2K C_L]/(C_D0 + K C_L²)² = **(C_D0 − K C_L²)/(C_D0 + K C_L²)²**.   (4)

A maximum of a smooth function on an open interval occurs only where the derivative is zero (`@openstax-calculus-volume-1` §4.3: a **critical point**). The denominator never vanishes, so f′ = 0 exactly when

**C_L* = √(C_D0/K)**,   and then   **(L/D)max = 1/(2√(C_D0 K))**.   (5)

In words: best glide is where induced drag equals parasite drag, K C_L*² = C_D0, and the total drag coefficient is 2C_D0. For the glider: C_L* = √(0.015/0.017684) = 0.9210, C_D = 0.030, (L/D)max = 1/(2 × 0.016287) = 30.70. The still-air glide from 1000 m covers 30.70 km.

Is it a maximum? Two ways to be sure. **Second-derivative test** (§4.5): differentiate (4) again,

f″(C_L) = −2K C_L (3C_D0 − K C_L²)/(C_D0 + K C_L²)³,   (6)

and at C_L*, where K C_L*² = C_D0, this is −K C_L*/(2C_D0²) = −36.19 for the glider: negative, so f is concave down there and the critical point is a maximum. **Behaviour at the ends** (§4.6, §4.8): f(0) = 0, and as C_L → ∞ the ratio C_L/(C_D0 + K C_L²) is of the form ∞/∞; **l'Hôpital's rule** replaces it by the ratio of derivatives, 1/(2K C_L) → 0. A positive function that is zero at both ends of (0, ∞) with a single critical point has its maximum there. Either argument suffices; engineers use both because each catches a different kind of algebra mistake.

Notice how flat the optimum is: L/D is within 1 % of its maximum for every C_L between 0.80 and 1.06. A pilot flying 5 % off the best-glide speed loses almost nothing in range; the first derivative is zero at the optimum, so only the second-order term is left. The flatness of optima is why sensitivities (2.6) matter more than exact optima in design.

### 2.4 Best-glide speed and sink

From (2) with cos γ ≈ 1 and sea-level density, V* = √(2 × 3432.3/(1.225 × 10.5 × 0.9210)) = 24.07 m/s (87 km/h). The exact form, with cos γ = 0.99947, gives 24.066 m/s. The sink rate is V/(L/D) = 0.784 m/s small-angle, 0.783 m/s exact. From 1000 m the glide lasts 21.3 minutes.

### 2.5 Minimum sink: a different derivative, a different optimum

Time in the air is set by the sink rate, w = V sin γ ≈ V C_D/C_L. Substituting V from (2):

w ≈ √(2W/(ρS)) · C_D / C_L^(3/2).

So the minimum-sink C_L minimises g(C_L) = (C_D0 + K C_L²)/C_L^(3/2). Differentiate (power and quotient rules):

g′ = [2K C_L · C_L^(3/2) − (C_D0 + K C_L²) · (3/2) C_L^(1/2)] / C_L³,

which vanishes when 2K C_L² = (3/2)(C_D0 + K C_L²), i.e. **C_L,ms = √(3C_D0/K) = √3 · C_L***.   (7)

For the glider: C_L,ms = 1.595, V = 18.28 m/s, sink 0.687 m/s, L/D = 26.59. Flying for endurance costs 13 % of the glide ratio and gains 14 % in time aloft (24.3 minutes from 1000 m instead of 21.3). The ratio of the two sink rates is 3^(3/4)/2 = 1.140, independent of the aircraft; the code reproduces 1.1401, the 0.03 % difference being the cos γ terms the small-angle derivation dropped.

### 2.6 Linear approximation: how much does dirt cost?

The **linear approximation** (§4.2) of a function near x₀ is the tangent line, f(x₀ + Δx) ≈ f(x₀) + f′(x₀)Δx. Apply it to (L/D)max = ½(C_D0 K)^(−1/2) as a function of C_D0:

d(L/D)max/dC_D0 = −¼ C_D0^(−3/2) K^(−1/2) = **−(L/D)max/(2C_D0)**,   (8)

a sensitivity of −1023 per unit C_D0, or −1.02 in L/D per 0.001 of C_D0. A 10 % increase in parasite drag (ΔC_D0 = 0.0015) predicts Δ(L/D) = −1.535, i.e. 29.16 against the exact 29.27: the linear estimate is 0.36 % low. The error is the next term, ½ f″ Δx² = (3/8)(L/D)max(ΔC_D0/C_D0)² = 0.115, which is what the code finds (0.106; the remainder is the third-order term). Halve the perturbation and the error falls four-fold: second order, and the test checks exactly that.

Engineering reading: each 0.001 of C_D0 (a mis-rigged fairing, a dirty wing, an open vent) costs 1 L/D point and 1 km of glide from 1000 m. A 30 % increase, which a badly iced wing can produce, costs 3.8 L/D points, and there the linear estimate is 3 % low: the tangent line is a local tool, as P01.01 warned.

### 2.7 Related rates: the speed to fly changes as you descend

The best-glide *lift coefficient* is fixed by (5), but the best-glide *speed* V* = √(2W/(ρ S C_L*)) depends on density, and density rises as the glider descends. Writing V* = C ρ^(−1/2) and applying the chain rule (§4.1) through ρ(h(t)):

**dV*/dt = (dV*/dρ)(dρ/dh)(dh/dt) = −(V*/(2ρ)) (dρ/dh) ḣ**.   (9)

In the ISA troposphere ρ = ρ₀(T/T₀)ⁿ with n = −g₀/(RL) − 1 = 4.2559, so dρ/dh = ρ n L/T. At 1000 m: ρ = 1.1116 kg/m³, dρ/dh = −1.092 × 10⁻⁴ kg/m⁴, V* = 25.27 m/s, and the glider descends at its own sink rate ḣ = −0.822 m/s. Then dV*/dt = −(25.27/2.223)(−1.092 × 10⁻⁴)(−0.822) = −1.02 × 10⁻³ m/s². Over the 20-minute glide that is −1.2 m/s in all: the pilot should bleed off about 0.12 m/s of airspeed per 100 m of descent. (Small, which is why the briefing gives one speed; but from 3000 m the best-glide speed is 27.9 m/s against 24.1 at sea level, a 16 % difference that does matter.)

### 2.8 Headwind: an optimum with no closed form

In a headwind V_w the ground speed is V − V_w, and the ground covered per metre of height is G(V) = (V − V_w)/w(V), where w(V) is the sink rate as a function of airspeed. With C_L = 2W/(ρSV²) in (3), the small-angle sink rate is

**w(V) = a V³ + b/V**,   a = ρ S C_D0/(2W),   b = 2 W K/(ρ S).   (10)

(Check the still-air optimum from (10): maximising V/w means minimising w/V = aV² + b/V², zero derivative at V⁴ = b/a, which is (2) with C_L* inserted: 24.07 m/s again.) Maximising G by the quotient rule, G′ = 0 gives

**F(V) ≡ w(V) − (V − V_w) w′(V) = 0**.   (11)

Geometrically: on the plot of sink against airspeed, the best speed is where the line from the point (V_w, 0) is tangent to the polar; `glide_polar.png` draws it. Equation (11) is a quartic in V with no useful closed form, so we solve it with **Newton's method** (§4.9):

V_{n+1} = V_n − F(V_n)/F′(V_n),   F′(V) = −(V − V_w) w″(V).

For a 10 m/s headwind: V = 27.81 m/s (C_L = 0.690), ground glide ratio 18.87 against 17.95 if the pilot stubbornly flies the still-air speed: 5 % more range for flying 3.7 m/s faster. In a 15 m/s wind, 31.0 m/s and 14.0; in a 20 m/s wind, 35.4 m/s and 10.2. The rule of thumb "add a third to a half of the headwind to the best-glide speed" comes straight from this table and is tested in part 4.

### 2.9 Newton's method: why it is fast and when it is not

Newton's method replaces f by its tangent line at xₙ and steps to where the line crosses zero. If f is smooth and f′(root) ≠ 0, the error eₙ₊₁ ≈ (f″/2f′) eₙ²: **quadratic convergence**, the number of correct digits doubles every step. For √2 from x₀ = 1 the errors are 0.41, 0.086, 0.0025, 2.1 × 10⁻⁶, 1.6 × 10⁻¹²: six steps to machine precision. The method has two failure modes, both of which the code must report rather than hide: f′(xₙ) = 0 (division by zero), and a start outside the **basin of attraction**. The basin for finding C_L* as the root of (4) ends at the inflection point of L/D, where f″ = 0; from (6) that is C_L = √(3C_D0/K), the minimum-sink C_L. Start Newton at C_L = 2.0 and it runs off to infinity; start at 1.4 and it converges to −C_L*, the L/D *minimum* in inverted flight, which is a true root of (4) and the wrong answer. Newton is fast, not safe: bracket first, or start close.

### 2.10 Where the model stops being valid

- **Steady glide.** Immediately after an engine failure the aircraft decelerates from cruise speed; the transition costs height that (1) does not include.
- **Parabolic polar.** Near stall, C_D rises faster than K C_L² and C_L stops increasing; the minimum-sink C_L of 1.6 may be at or beyond C_Lmax for a trainer. At low C_L, cambered sections have their minimum drag at C_L ≈ 0.2–0.3, not 0, which shifts both optima.
- **Constant density along the glide.** Handled by (9) for the speed; for the range, the glide ratio is independent of ρ and the error is nil.
- **Uniform wind.** Wind usually weakens near the ground; the headwind optimum then changes during the descent, and a single speed is a compromise.
- **No thermals or sink.** Rising air at 1 m/s doubles the endurance of this glider and changes the speed-to-fly entirely (extension).

## 3. Build the computational model

Open `starter/main.py`. The glider's data are in the `GLIDER` dictionary; every function takes its parameters explicitly so the tests can vary them.

### Step 1: polar, L/D and its two derivatives

```python
import numpy as np

def induced_drag_factor(aspect_ratio, oswald):
    if aspect_ratio <= 0.0 or oswald <= 0.0:
        raise ValueError("aspect ratio and Oswald factor must be positive")
    return 1.0 / (np.pi * oswald * aspect_ratio)                         # (3)

def drag_polar(cl, cd0, k):
    return cd0 + k * np.asarray(cl, dtype=float) ** 2

def lift_to_drag(cl, cd0, k):
    cl = np.asarray(cl, dtype=float)
    return cl / drag_polar(cl, cd0, k)

def ld_derivative(cl, cd0, k):
    cl = np.asarray(cl, dtype=float)
    return (cd0 - k * cl**2) / (cd0 + k * cl**2) ** 2                     # (4)

def ld_second_derivative(cl, cd0, k):
    cl = np.asarray(cl, dtype=float)
    return -2.0 * k * cl * (3.0 * cd0 - k * cl**2) / (cd0 + k * cl**2) ** 3   # (6)
```

Type the derivatives from your own derivation, then let the test compare them with a central difference of `lift_to_drag`: that is how you find a dropped factor of two.

### Step 2: the closed-form optima and the glide geometry

```python
def best_glide_cl(cd0, k):      return float(np.sqrt(cd0 / k))             # (5)
def max_lift_to_drag(cd0, k):   return float(1.0 / (2.0 * np.sqrt(cd0 * k)))
def min_sink_cl(cd0, k):        return float(np.sqrt(3.0 * cd0 / k))       # (7)

def glide_angle(cl, cd0, k):
    cl = np.asarray(cl, dtype=float)
    return np.arctan(drag_polar(cl, cd0, k) / cl)                          # (1), exact

def airspeed(cl, weight, rho, area, gamma=0.0):
    cl = np.asarray(cl, dtype=float)
    if np.any(cl <= 0.0) or weight <= 0.0 or rho <= 0.0 or area <= 0.0:
        raise ValueError("C_L, weight, density and area must be positive")
    return np.sqrt(2.0 * weight * np.cos(gamma) / (rho * area * cl))       # (2)

def sink_rate(cl, cd0, k, weight, rho, area):
    gamma = glide_angle(cl, cd0, k)
    return airspeed(cl, weight, rho, area, gamma) * np.sin(gamma)
```

`sink_rate_of_speed(v, ...)` is (10) and `glide_range(h, cd0, k)` is h(L/D)max. Keeping both the exact and the small-angle sink lets the test measure the small-angle error (0.08 % at best glide) instead of assuming it.

### Step 3: Newton's method, generic

```python
def newton(f, dfdx, x0, tol=1e-12, maxiter=50):
    x = float(x0)
    for n in range(1, maxiter + 1):
        d = dfdx(x)
        if d == 0.0 or not np.isfinite(d):
            raise ZeroDivisionError("derivative vanished in Newton's method")
        x_new = x - f(x) / d
        if not np.isfinite(x_new):
            raise RuntimeError("Newton's method diverged")
        if abs(x_new - x) <= tol * max(1.0, abs(x_new)):
            return x_new, n
        x = x_new
    raise RuntimeError(f"Newton's method did not converge in {maxiter} iterations")
```

The stopping rule is relative for large |x| and absolute near zero, a standard compromise. Returning the iteration count lets a test assert quadratic convergence (√2 in at most 7 steps). `best_glide_cl_newton` applies this to (4) with (6) as its derivative.

### Step 4: sensitivity and related rate

```python
def ld_max_sensitivity_cd0(cd0, k):
    return -max_lift_to_drag(cd0, k) / (2.0 * cd0)                          # (8)

def linear_approximation(f, dfdx, x0, dx):
    return float(f(x0) + dfdx(x0) * dx)

def best_glide_speed_rate(h, hdot, cd0, k, weight, area):
    v = float(best_glide_speed(h, cd0, k, weight, area))
    return -(v / (2.0 * isa_density(h))) * float(isa_density_gradient(h)) * hdot   # (9)
```

`isa_density` is the troposphere formula from `P00.03` and `isa_density_gradient` its derivative, dρ/dh = ρ n L/T, which the test checks against a 1 m central difference.

### Step 5: the headwind optimum

```python
def headwind_best_speed(v_wind, cd0, k, weight, rho, area):
    a = rho * area * cd0 / (2.0 * weight)
    b = 2.0 * weight * k / (rho * area)
    w  = lambda v: a * v**3 + b / v                                          # (10)
    w1 = lambda v: 3.0 * a * v**2 - b / v**2
    w2 = lambda v: 6.0 * a * v + 2.0 * b / v**3
    f  = lambda v: w(v) - (v - v_wind) * w1(v)                               # (11)
    df = lambda v: -(v - v_wind) * w2(v)
    v_still = float((b / a) ** 0.25)
    root, _ = newton(f, df, 1.2 * v_still + max(v_wind, 0.0))
    return root
```

The start is deliberately on the fast side of the optimum, where F is monotone, so Newton walks down to the root for any wind the glider can make headway against.

### Verify

Run `python -m pytest curriculum/phase-01/lesson-02-applications-of-derivatives`. The tests:

1. **Analytic optimum.** f′(C_L*) = 0 to 10⁻¹², f″(C_L*) < 0, K C_L*² = C_D0, and L/D at 500 sampled C_L values never exceeds (L/D)max.
2. **Newton.** From 0.3, 0.6 and 1.2 the root agrees with (5) to 10⁻¹⁰; (L/D)″ vanishes at the min-sink C_L; from 2.0 the method raises `RuntimeError`. On x² − 2 it converges in ≤ 7 iterations; on x² + 1 it reports failure; on a zero derivative it raises.
3. **Derivatives.** (4) and (6) match central differences to 10⁻⁷ and 10⁻⁶; dρ/dh matches a 1 m difference to 10⁻⁷; dV*/dt matches a 1 s difference to 10⁻⁶.
4. **Limiting cases and scaling.** C_D(0) = C_D0; L/D(0) = 0; V ∝ √W and ∝ 1/√ρ; glide range ∝ h and ∝ 1/√C_D0; zero wind gives the still-air speed to 10⁻⁹; ḣ = 0 gives dV*/dt = 0.
5. **Approximation quality.** Small-angle sink within 0.2 % of exact and always slightly above it; linearisation error falls with rate 2.0 ± 0.1 as the step halves.
6. **Optimality in wind.** The headwind speed beats every sampled speed in [0.6V, 1.6V] for V_w = −10, 0, +10 m/s, lies above the still-air speed in a headwind and below it, but above min-sink speed, in a tailwind.

### Compare

`main()` hands −L/D to `scipy.optimize.minimize_scalar(method="bounded")` on [0.05, 3]. It returns C_L = 0.92099403 after 11 function evaluations, 6 × 10⁻⁹ below √(C_D0/K); the L/D there differs from (L/D)max by 4 × 10⁻¹⁵. The library did not use a derivative at all (Brent's method brackets a minimum with parabolic interpolation), and it still got eight digits, because the optimum is flat: a 6 × 10⁻⁹ error in C_L costs 10⁻¹⁵ in L/D. That flatness is why a bracketing method is the safe default for a one-dimensional optimum and why Newton on the derivative is the fast one: Brent needs 11 evaluations and a bracket, Newton needs 5 and a good start. In Phase 09 you will meet problems with 50 design variables where neither luxury is affordable and gradients (P01.09) are the only way.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Three sweeps and one error measurement.

1. **Drag sweep.** (L/D)max for C_D0 = 0.010, 0.015, 0.020, 0.030 (K fixed): predict the scaling first (∝ C_D0^(−1/2)), then run: 37.6, 30.7, 26.6, 21.7.
2. **Span sweep.** (L/D)max for AR = 8, 12, 20, 30 with e = 0.9 and C_D0 fixed: 19.4, 23.8, 30.7, 37.6. Note the symmetry with the first sweep: (L/D)max depends only on the product C_D0 K.
3. **Wind sweep.** Best speed and ground glide ratio for V_w = 0, 5, 10, 15, 20 m/s: (24.1, 30.7), (25.6, 24.5), (27.8, 18.9), (31.0, 14.0), (35.4, 10.2). Compute the speed increase as a fraction of the wind: 0.30, 0.37, 0.46, 0.57. **Prediction to write down first:** is "add half the headwind" a good rule, and in which range?
4. **Linearisation error.** |linear − exact| for (L/D)max against ΔC_D0/C_D0 = 0.01, 0.02, 0.04, 0.1, 0.3 on log-log axes; fit the slope.

Raw data format: three CSV tables (`cd0, ld_max`; `ar, ld_max`; `v_wind, v_best, ground_ld`) and one (`frac, err`), with a one-line header naming the fixed parameters.

What you should find: the drag and span sweeps lie on the same power law; the wind sweep shows the speed increase growing from 30 % to 57 % of the wind, so "half the headwind" is right only near 15–20 m/s and overshoots at light winds (where the correct increase is a third); the linearisation error has slope 2.0 on the log-log plot up to 10 % and bends upward beyond, where the third-order term appears.

## 5. Validate and improve

| quantity (glider, sea level) | hand (§2) | code | difference | cause |
|---|---|---|---|---|
| C_L* | 0.9210 | 0.920994 | 0 | exact formula |
| (L/D)max | 30.70 | 30.6998 | 0 | exact formula |
| V* (m/s) | 24.07 (cos γ = 1) | 24.066 (exact) | 0.03 % | cos γ = 0.99947 |
| sink at V* (m/s) | 0.784 (V/(L/D)) | 0.783 | 0.08 % | cos^(3/2) γ in w = V sin γ |
| sink ratio w*/w_ms | 3^(3/4)/2 = 1.1398 | 1.1401 | 0.03 % | small-angle terms |
| Δ(L/D)max for +10 % C_D0 | −1.535 (linear) | −1.429 (exact) | 0.36 % of L/D | ½f″Δx² = 0.115 |
| dV*/dt at 1000 m (m/s²) | −1.02 × 10⁻³ | −1.0207 × 10⁻³ | 0 | chain rule is exact |
| best speed in 10 m/s headwind | 27.8 (graph tangent) | 27.808 (Newton) | reading accuracy | — |

Every difference has a named cause, and all are far smaller than the uncertainty in the inputs: nobody knows C_D0 of a real aircraft to better than 5 %, and the Oswald factor to better than 10 %, so (L/D)max, which goes as (C_D0 K)^(−1/2), is uncertain by about ½√(5² + 10²) ≈ 6 % and the still-air range from 1000 m by ±1.7 km. The sensitivity (8) is the right tool for that statement, and it is a derivative.

What the agreement does not show: that a real glider follows (3). Published polars for training gliders give (L/D)max of 28–35 at best-glide speeds near 25 m/s, so the numbers are realistic, but the parabolic fit is good only in the middle of the polar. Near stall the real sink rate rises much faster than (10), so the minimum-sink prediction of C_L = 1.6 is the least trustworthy number in this lesson. A measured polar (a flight test with a logged airspeed and vertical speed, or a wind-tunnel table) is what validates it.

Improvements worth making: replace `isa_density` with the full seven-layer `isa` from `P00.03`; let `headwind_best_speed` accept a vertical air velocity (rising air) and derive the MacCready speed-to-fly; bracket the Newton start with a coarse grid so the basin problem of §2.9 cannot arise. Improvement worth *not* making: trusting (7) beyond 90 % of C_Lmax.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the hand derivation of (5) and (7) with the second-derivative test; `glide_polar.png`; the four sweep tables; and a one-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) giving the glider's two speeds, the headwind rule you extracted with its range of validity, and the assumption you trust least.

### Adversarial self-checks

- **What if the weight doubles?** C_L*, (L/D)max and the range are unchanged; V* rises by √2 to 34.0 m/s and the sink rate with it, so the glide takes 15 minutes instead of 21. Check all four with the code. (This is why gliders carry water ballast: same range, faster.)
- **What if the speed goes to zero? To infinity?** w → ∞ at both ends of (10): b/V at low speed (induced drag) and aV³ at high speed (parasite). The optimum exists because both ends are bad.
- **Which assumption fails first in the real world?** The parabolic polar near stall, at the minimum-sink condition. The second is "steady": the first 10 s after an engine failure are not a glide.
- **What measurement could falsify the model?** A logged descent at several fixed airspeeds: the measured sink against V should be a V³ + b/V. If the minimum sink occurs at a lower C_L than √3 C_L*, the polar is not parabolic at high lift.
- **Does the result survive a change in solver or start?** Newton from 0.3, 0.6, 1.2 and Brent's method agree to 10⁻⁸; Newton from 2.0 does not converge; from 1.4 it converges to the wrong root. Say which starts you would allow in production code.
- **Which decision would change if the uncertainty doubled?** If C_D0 were known only to ±20 % instead of ±10 %, the range from 1000 m would be uncertain by ±3 km instead of ±1.5 km, and a landing field 29 km away would move from "reachable with margin" to "coin toss". The decision (turn back or go on) depends on the sensitivity, not on the nominal number.

## Check yourself

1. In a steady glide, which ratio sets the glide angle, and why does a heavier glider not descend more steeply?
2. At best glide, how do induced and parasite drag compare, and what is C_D in terms of C_D0?
3. State the second-derivative test and apply it in one line to confirm C_L* is a maximum of L/D.
4. What does l'Hôpital's rule say about the limit of C_L/(C_D0 + K C_L²) as C_L → ∞, and why does that help prove the maximum?
5. The minimum-sink C_L is √3 times the best-glide C_L. Is the minimum-sink *speed* higher or lower than the best-glide speed, and by what factor?
6. Write the linear approximation of (L/D)max for a change ΔC_D0 and give the size of the error term in terms of (L/D)max and ΔC_D0/C_D0.
7. In a related-rates problem, which rule links dV*/dt to dh/dt, and what sign does dV*/dt have in a descent?
8. Newton's method for √2 from x₀ = 1 has errors 0.41, 0.086, 0.0025, 2 × 10⁻⁶. What is this convergence called, and what would stop it working for the best-glide root started at C_L = 2?
