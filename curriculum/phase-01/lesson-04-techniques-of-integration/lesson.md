# P01.04  Techniques of Integration

> Every technique is a change of viewpoint that turns an integral you cannot do into one you can; the numerical rule is the viewpoint that always works, and its error bound is the price.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the roll moment of inertia of a tapered wing for three planforms using integration by parts, trigonometric substitution and a direct polynomial integral, and derive the take-off time under quadratic drag with partial fractions. 2. Implement and test the midpoint, trapezoid and Simpson rules with their error bounds, and show Simpson is exact for the linear-taper wing. 3. Measure how an end-point singularity reduces every rule to order 1.5 and how a substitution restores the full order (and, for a periodic integrand, far more). 4. Decide, from a cutoff sequence, whether an improper integral converges, and what the limiting value means for a wing of "infinite" span. |
| **Prerequisites** | `P01.03` (Riemann sums, the Fundamental Theorem, substitution), `P01.01` (orders of convergence), `P00.09` (SymPy `integrate`). Skill check: you can differentiate a product and a composite, you know sin² + cos² = 1 and the double-angle formulas, and you can write 1/((a − x)(a + x)) as a sum of two fractions. |
| **Engineering question** | What is the roll moment of inertia of a light aircraft's wing pair (which sets how quickly it responds to aileron), how does it depend on the planform, and how many integration panels does a numerical answer need to be trusted to four figures? |
| **Assumptions and limits** | The wing is a slender beam along the span with mass per unit length proportional to the local chord, λ = σ c(y), with constant areal density σ; the chordwise extent of the mass is ignored (it adds a term of order (c/s)² ≈ 3 % to I_xx); no fuel, engines or landing gear in the wing; the three planforms are idealisations (a real wing has a constant centre section and rounded tips). The numerical rules assume a smooth integrand; at the elliptic tip the chord has a square-root end and every rule slows down. The take-off model ignores rolling friction and the fall of propeller thrust with speed, both of which matter by 10–30 % in practice. |
| **Analytical work** | Three closed forms for I_xx (polynomial, trigonometric substitution, integration by parts twice); the improper integral for an exponentially tapered beam of unbounded span; the partial-fraction integral for the take-off time; the error laws of the trapezoid and Simpson rules and why Simpson is exact for cubics; hand values for the lesson's wing. |
| **Computational deliverable** | `starter/main.py`: `chord_linear`, `chord_elliptic`, `chord_exponential`, `wing_mass_linear`, `spanwise_centroid_linear`, `roll_inertia_linear`, `roll_inertia_elliptic`, `roll_inertia_exponential`, `roll_inertia_exponential_infinite`, `symbolic_roll_inertia_exponential`, `symbolic_roll_inertia_elliptic`, `midpoint`, `trapezoid`, `simpson`, `trapezoid_error_bound`, `simpson_error_bound`, `roll_inertia_numeric`, `roll_inertia_elliptic_substituted`, `improper_integral`, `takeoff_time`, `takeoff_time_numeric`. Tests prove: the closed forms agree with SymPy and with Simpson, Simpson is exact for cubics and reproduces the linear-taper inertia with two panels, trapezoid is second order and Simpson fourth on a smooth integrand, the error bounds hold and are tight to a factor 20, the tip singularity reduces Simpson to order 1.5 and the θ-substitution makes the trapezoid rule exact at n = 4, the improper integral converges monotonically to the closed-form limit while ∫dx/x does not, and the take-off time matches its numerical integral and both limits. |
| **Practical activity** | **computational experiment**: a convergence study of three rules on the elliptic planform with and without the trigonometric substitution; a comparison of measured errors with the textbook bounds; a cutoff study of the improper integral. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `beam_inertia_convergence.png` with its table, the three hand derivations, and a short note on which rule and which substitution you would use for each planform and why. |
| **Assessment** | `assessments/phase-01/grading-rubric.md`; competencies 1–4. Pass: all three closed forms derived with the technique named at each step; all public tests pass; the hand value for the lesson's wing agrees with the code to four figures; the convergence note explains the 1.5 and the exactness of the substituted trapezoid rule. |
| **References** | `@openstax-calculus-volume-2` §3.1 "Integration by Parts", §3.2 "Trigonometric Integrals", §3.3 "Trigonometric Substitution", §3.4 "Partial Fractions", §3.6 "Numerical Integration", §3.7 "Improper Integrals" (all required-free). `@openstax-calculus-volume-2` §3.5 "Other Strategies for Integration" (supplementary-free). |
| **Extensions** | ★ Add the chordwise contribution to I_xx (each strip is a thin plate of width c: add σ ∫ c³/12 dy) and show it is 3 % for this wing. ★ Add rolling friction μ(W − L) to the take-off model and show the integral is still a partial-fraction one, with a different v_t. |
| **Time** | unmeasured |

## 1. The engineering problem

Push the stick sideways and the ailerons produce a rolling moment; how fast the aircraft starts to roll is set by its **roll moment of inertia** I_xx = ∫ y² dm, the sum of every mass element times the square of its distance from the roll axis. For most aircraft the wings dominate I_xx, because their mass sits far from the centreline. A design study that compares planforms, a flight-dynamics model (Phase 15) and a simulator all need the number, and they need it before the wing is built, from the planform and a mass model.

The simplest honest mass model is a slender beam: at spanwise station y the wing has mass per unit length λ(y) = σ c(y), with σ an areal density (kg per m² of planform, about 10–15 kg/m² for a light aircraft wing) and c(y) the chord. Then I_xx = 2 ∫₀ˢ y² σ c(y) dy over the two half-wings of semispan s. For a **linearly tapered** chord this is a polynomial integral; for an **elliptic** planform it has a square root that needs a trigonometric substitution; for an **exponential** taper it needs integration by parts, twice; and for a wing of unbounded span it is an **improper integral** that converges. A fourth integral in the lesson, the time an aircraft needs to accelerate on the runway against drag that grows with the square of speed, needs **partial fractions**. These are the four techniques of the title, each attached to a question with a number as its answer.

The second half of the engineering question is numerical. A real planform is a table of chords, not a formula, so the integral will be done by a rule: midpoint, trapezoid or Simpson. Each has an error bound you can compute before you start, and each is defeated in the same way by a square-root tip. Knowing which rule, how many panels, and when to substitute first is the practical skill.

Success criteria, written first: the three closed forms agree with SymPy exactly and with Simpson's rule to four figures; Simpson reproduces the linear-taper inertia *exactly* with two panels (the integrand is a cubic); the measured errors lie below the textbook bounds and within a factor of three of them; the elliptic integral after substitution is exact with four trapezoid panels; the improper integral's cutoff sequence approaches the closed-form limit from below; and the take-off time's closed form matches its numerical integral to 10⁻⁸.

## 2. Understand and derive

### 2.1 The beam model and the lesson's wing

```
            y = 0 (roll axis)                   y = s (tip)
    ──────────┼───────────────────────────────────┤
              │  chord c(y)   ┌──────────────┐    │   mass per unit span  λ(y) = σ c(y)
              │   ▲           │   dm = σ c dy│    │   I_xx = 2 ∫₀ˢ y² σ c(y) dy  (both wings)
              │   │  c_r      │              │ c_t│
              │   ▼           └──────────────┘    │
    ──────────┼───────────────────────────────────┤
                              ◄── dy ──►  at distance y
```

**Assumption 1 (slender beam):** the chordwise spread of mass is ignored; the whole strip dm = σ c dy is at distance y. **Assumption 2 (constant σ):** the wing's structure weighs the same per square metre everywhere. The lesson's wing: span b = 11 m (s = 5.5 m), root chord c_r = 1.6 m, taper ratio τ = c_t/c_r = 0.6, σ = 12 kg/m². Planform area 14.08 m², wing mass 169 kg: realistic for a 1000 kg light aircraft. Units of I_xx: (kg/m²)(m)(m²)(m) = kg m². Good.

### 2.2 Linear taper: a polynomial, and the first check on Simpson

c(y) = c_r [1 − (1 − τ) y/s]. Then

I_xx = 2σ c_r ∫₀ˢ [y² − (1 − τ) y³/s] dy = 2σ c_r [s³/3 − (1 − τ) s³/4] = **2σ c_r s³ [1/3 − (1 − τ)/4]**.   (1)

For the lesson's wing: 2 × 12 × 1.6 × 166.375 × (0.33333 − 0.10000) = 38.4 × 166.375 × 0.23333 = **1490.7 kg m²**. Two sanity checks: a rectangular wing (τ = 1) gives 2σ c_r s³/3 = 2129.6, a triangular one (τ = 0) gives 2σ c_r s³/12 = 532.4, and ours sits between. The mass is 2σ ∫c dy = σ c_r s (1 + τ) = 169.0 kg, and the spanwise centroid of each half-wing is ȳ = ∫y c dy/∫c dy = s(1 + 2τ)/(3(1 + τ)) = 2.521 m, between a triangle's s/3 and a rectangle's s/2.

The integrand y² c(y) is a **cubic**. Simpson's rule (§2.6) is exact for every cubic, so `roll_inertia_numeric` with n = 2 returns 1490.72 to all sixteen digits. That is the first test, and it catches any slip in (1).

### 2.3 Elliptic planform: trigonometric substitution

c(y) = c_r √(1 − (y/s)²). The integrand y² √(1 − y²/s²) has no polynomial antiderivative. **Trigonometric substitution** (`@openstax-calculus-volume-2` §3.3): set y = s sin θ, so dy = s cos θ dθ, √(1 − y²/s²) = cos θ (for 0 ≤ θ ≤ π/2), and the limits become 0 and π/2:

∫₀ˢ y² √(1 − y²/s²) dy = s³ ∫₀^(π/2) sin²θ cos²θ dθ.

A **trigonometric integral** (§3.2): sin²θ cos²θ = ¼ sin²2θ = ⅛(1 − cos 4θ), whose integral over [0, π/2] is ⅛ × π/2 = π/16. Hence

**I_xx = 2σ c_r s³ π/16 = π σ c_r s³/8**.   (2)

For the lesson's root chord and span: π × 12 × 1.6 × 166.375/8 = **1254.4 kg m²**, 84 % of the tapered wing's, with 2 % less mass (165.9 kg): the ellipse carries its mass further inboard. SymPy, asked to integrate in y, runs for minutes; asked to integrate sin²cos² in θ it returns π/16 at once. The substitution is not only a hand technique.

### 2.4 Exponential taper: integration by parts, twice

c(y) = c_r e^(−y/ℓ), with ℓ = s/ln(1/τ) = 10.767 m so that the tip chord matches the tapered wing's. **Integration by parts** (§3.1), ∫u dv = uv − ∫v du, applied to ∫y² e^(−y/ℓ) dy with u = y², dv = e^(−y/ℓ) dy (so v = −ℓ e^(−y/ℓ)):

∫y² e^(−y/ℓ) dy = −ℓ y² e^(−y/ℓ) + 2ℓ ∫y e^(−y/ℓ) dy.

Parts again on the remaining integral (u = y, dv = e^(−y/ℓ) dy): ∫y e^(−y/ℓ) dy = −ℓ y e^(−y/ℓ) + ℓ ∫e^(−y/ℓ) dy = −ℓ y e^(−y/ℓ) − ℓ² e^(−y/ℓ). Collecting:

∫y² e^(−y/ℓ) dy = −e^(−y/ℓ) (ℓ y² + 2ℓ² y + 2ℓ³) + C.

Evaluate from 0 to s with the Fundamental Theorem:

**I_xx = 2σ c_r [2ℓ³ − e^(−s/ℓ) (s²ℓ + 2sℓ² + 2ℓ³)]**.   (3)

Numbers: e^(−s/ℓ) = τ = 0.6; s²ℓ = 325.70, 2sℓ² = 1275.18, 2ℓ³ = 2496.32, sum 4097.20; bracket 2496.32 − 0.6 × 4097.20 = 38.00; times 38.4: **1459.1 kg m²**. Between the ellipse and the trapezoid, as the planform is. SymPy's `integrate` reproduces (3) in a different arrangement; the test simplifies the difference to zero. The pattern of parts, "differentiate the polynomial, integrate the exponential, repeat until the polynomial is gone", is the one you will use for every Laplace transform in Phase 02.

### 2.5 Improper integral: a wing of unbounded span

Let the exponential wing extend to infinity. The mass 2σ c_r ∫₀^∞ e^(−y/ℓ) dy = 2σ c_r ℓ is finite, and so is the inertia. An **improper integral** (§3.7) is defined as a limit of cutoffs:

∫₀^∞ y² e^(−y/ℓ) dy = lim_{R→∞} ∫₀ᴿ y² e^(−y/ℓ) dy = lim_{R→∞} [2ℓ³ − e^(−R/ℓ)(R²ℓ + 2Rℓ² + 2ℓ³)] = 2ℓ³,

because e^(−R/ℓ) kills every power of R. So **I_xx(∞) = 4σ c_r ℓ³** = 95 859 kg m² for ℓ = 10.767 m: sixty-five times the real wing's, because the fictitious span is long and the inertia weights y². The fraction captured within R is 1 − e^(−R/ℓ)(R²/2ℓ² + R/ℓ + 1): 8 % at R = ℓ, 32 % at 2ℓ, 76 % at 4ℓ, 98.6 % at 8ℓ. The code builds the same sequence numerically and the test asks it to be monotone, to stay below the limit, and to reach it to 10⁻⁴ at 16ℓ. Contrast ∫₁^∞ dx/x: each decade of R adds ln 10 = 2.303 and the sequence never settles. The convergence test is not a formality; a model that produces a divergent integral is telling you it has no cutoff where physics has one.

### 2.6 The numerical rules and what they cost

For a smooth f on [a, b] with n panels of width h (§3.6):

- **Midpoint**: error ≈ (h²/24)[f′(b) − f′(a)], bound (b − a)³ max|f″|/(24 n²), from P01.03.
- **Trapezoid**: T_n = h[½f₀ + f₁ + … + f_{n−1} + ½fₙ]; error ≈ −(h²/12)[f′(b) − f′(a)], bound **(b − a)³ max|f″|/(12 n²)**. Opposite sign and twice the midpoint error, so the weighted combination (2M + T)/3 cancels the h² term; that combination *is*
- **Simpson**: S_n = (h/3)[f₀ + 4f₁ + 2f₂ + 4f₃ + … + 4f_{n−1} + fₙ], n even, bound **(b − a)⁵ max|f⁽⁴⁾|/(180 n⁴)**. Fourth order, and exact for every cubic because the error term carries f⁽⁴⁾.

For the exponential wing with the trapezoid rule, |f″| of y² e^(−y/ℓ) is largest at y = 0 where it equals 2, so the bound is (5.5)³ × 2 × 2σ c_r/(12 n²) = 16.6 kg m² at n = 8; the measured error is 7.4. At n = 32: bound 1.04, measured 0.47; at n = 128: 0.065 and 0.029. The bound is honoured and loose by a factor 2.2, because f″ is not 2 everywhere. Simpson at n = 8: bound 0.027, measured 0.018; at n = 32: 1.06 × 10⁻⁴ and 7.1 × 10⁻⁵. A bound you can compute *before* integrating tells you the n you need; a measured error afterwards tells you the bound was honest.

### 2.7 Partial fractions: the take-off roll

An aircraft of mass m accelerates under thrust T against aerodynamic drag k v² (k = ½ρ C_D S): m dv/dt = T − k v². The time to reach speed v is t = ∫₀ᵛ m dv′/(T − k v′²). Factor the denominator: T − k v² = (√T − √k v)(√T + √k v). **Partial fractions** (§3.4):

1/(T − k v²) = (1/(2√T)) [1/(√T − √k v) + 1/(√T + √k v)],

each piece a logarithm. Integrating from 0 to v,

**t(v) = (m/(2√(Tk))) ln[(√T + √k v)/(√T − √k v)] = (m/(2√(Tk))) ln[(v_t + v)/(v_t − v)]**,   (4)

with v_t = √(T/k) the speed at which drag equals thrust and acceleration stops. For a 1000 kg aircraft with T = 2500 N and k = 0.5 × 1.225 × 0.04 × 15 = 0.3675 kg/m: v_t = 82.5 m/s, m/(2√(Tk)) = 16.50 s, and the time to 30 m/s is 16.50 × ln(112.48/52.48) = 16.50 × 0.7624 = **12.58 s**, against 12.00 s without drag; to 40 m/s, 17.47 s against 16.00. The logarithm diverges as v → v_t: no finite time reaches the asymptotic speed. Simpson's rule on the integrand gives the same times to 10⁻⁸, which is the test, and it is also what you would do if thrust were a table T(v) with no closed form.

### 2.8 Where this stops being valid

- **Slender beam.** A strip of chord c has its own inertia about its centre, σ c · c²/12 per unit span; for this wing it adds about 3 % to I_xx (extension). Fuel tanks (dense, often outboard) can double the wing's contribution.
- **Constant σ.** Real wing structure is heavier at the root (bending loads); a λ ∝ c² law (Problem 7) is often closer and gives a smaller I_xx.
- **Idealised planforms.** A real wing has a constant-chord centre section, rounded tips and a fuselage that occupies the first metre. These are polynomial corrections to a polynomial integral; the technique is unchanged.
- **Rules on singular integrands.** The elliptic tip has c ∝ √(s − y): every rule is order 1.5 there (part 4). Substitute first.
- **Take-off.** Rolling friction μ(W − L) and a thrust that falls with speed change v_t and (4) by 10–30 %; the form of the integral survives, the numbers do not.

## 3. Build the computational model

Open `starter/main.py`. The wing's data are in `WING`.

### Step 1: planforms and closed forms

```python
import numpy as np

def chord_linear(y, c_r, taper, s):
    return c_r * (1.0 - (1.0 - taper) * np.asarray(y, dtype=float) / s)

def chord_elliptic(y, c_r, s):
    return c_r * np.sqrt(np.clip(1.0 - (np.asarray(y, dtype=float) / s) ** 2, 0.0, None))

def chord_exponential(y, c_r, ell):
    return c_r * np.exp(-np.asarray(y, dtype=float) / ell)

def roll_inertia_linear(sigma, c_r, taper, s):
    return 2.0 * sigma * c_r * s**3 * (1.0 / 3.0 - (1.0 - taper) / 4.0)             # (1)

def roll_inertia_elliptic(sigma, c_r, s):
    return np.pi * sigma * c_r * s**3 / 8.0                                            # (2)

def roll_inertia_exponential(sigma, c_r, ell, s):
    return 2.0 * sigma * c_r * (2.0 * ell**3 - np.exp(-s / ell) * (s**2 * ell + 2.0 * s * ell**2 + 2.0 * ell**3))   # (3)
```

The `np.clip` in the ellipse guards against a rounding-induced negative under the square root at y = s. `roll_inertia_exponential_infinite` is 4σ c_r ℓ³.

### Step 2: the three rules

```python
def trapezoid(f, a, b, n):
    if n < 1:
        raise ValueError("n must be at least 1")
    x = np.linspace(a, b, n + 1); fx = f(x); h = (b - a) / n
    return float(h * (0.5 * fx[0] + np.sum(fx[1:-1]) + 0.5 * fx[-1]))

def simpson(f, a, b, n):
    if n < 2 or n % 2:
        raise ValueError("Simpson's rule needs an even n >= 2")
    x = np.linspace(a, b, n + 1); fx = f(x); h = (b - a) / n
    return float(h / 3.0 * (fx[0] + fx[-1] + 4.0 * np.sum(fx[1:-1:2]) + 2.0 * np.sum(fx[2:-1:2])))
```

The slices `fx[1:-1:2]` (odd indices) and `fx[2:-1:2]` (even interior indices) are the 4s and 2s of Simpson's pattern. An odd n silently produces a wrong Simpson sum, so the code refuses it. `trapezoid_error_bound` and `simpson_error_bound` are the two formulas of §2.6; they take max|f″| or max|f⁽⁴⁾| as a number you supply, because bounding a derivative is your job, not the integrator's.

### Step 3: inertia by rule, with and without substitution

```python
def roll_inertia_numeric(chord, sigma, s, n, rule="simpson"):
    return 2.0 * sigma * RULES[rule](lambda y: y**2 * chord(y), 0.0, s, n)

def roll_inertia_elliptic_substituted(sigma, c_r, s, n):
    g = lambda th: np.sin(th) ** 2 * np.cos(th) ** 2
    return 2.0 * sigma * c_r * s**3 * trapezoid(g, 0.0, 0.5 * np.pi, n)
```

The substituted integrand is smooth and periodic with period π/2, and the trapezoid rule on a full period of a trigonometric polynomial is exact once n exceeds its degree: with sin²cos² = (1 − cos 4θ)/8 the rule is exact for n ≥ 2. The test asks n = 4 to agree with (2) to 10⁻¹⁴.

### Step 4: improper integral and take-off time

```python
def improper_integral(f, a, cutoffs, n=2000):
    return np.array([simpson(f, a, float(r), n) for r in cutoffs])

def takeoff_time(v, thrust, k, mass):
    vt = np.sqrt(thrust / k)
    if np.any(v >= vt) or np.any(v < 0.0):
        raise ValueError("speed must lie in [0, sqrt(T/k))")
    return mass / (2.0 * np.sqrt(thrust * k)) * np.log((vt + v) / (vt - v))          # (4)
```

`takeoff_time_numeric` is Simpson on m/(T − k v²). Refusing v ≥ v_t is physics, not pedantry: the integrand has a pole there and a Simpson sum across it returns a finite, wrong number.

### Verify

Run `python -m pytest curriculum/phase-01/lesson-04-techniques-of-integration`. The tests:

1. **Closed forms.** Simpson with n = 2 reproduces (1) to 10⁻¹³; (2) matches Simpson (n = 4000) to 10⁻⁵ and the substituted trapezoid (n = 4) to 10⁻¹⁴; (3) matches SymPy to 10⁻¹² and simplifies to the hand form; the s → ∞ and s → 0 limits of (3) are 4σ c_r ℓ³ and 0.
2. **Rules.** Simpson exact for a cubic with two panels; trapezoid and midpoint exact for a line; trapezoid not exact for a cubic; bad n and an unknown rule raise.
3. **Orders and bounds.** On eˣ, trapezoid rate 2.0 ± 0.05 and Simpson 4.0 ± 0.1; midpoint error is half the trapezoid's with opposite sign; on sin x the measured errors lie below the bounds and within a factor 20 of them; the Simpson bound falls 16-fold when n doubles.
4. **Singularity.** Simpson on the raw elliptic integrand has rates between 1.3 and 1.7 for n = 16 to 256; the substituted trapezoid rule is at 10⁻¹⁰ from n = 4.
5. **Improper.** The cutoff sequence is increasing, bounded by the limit, equal to (3) at each cutoff, and within 10⁻⁴ at R = 16ℓ; ∫dx/x grows by ln 10 per decade.
6. **Take-off.** (4) matches Simpson to 10⁻⁸ at three speeds; t(0) = 0; drag lengthens the roll; k → 0 gives m v/T; the time at 0.999 v_t exceeds 5× the time at 0.5 v_t; v ≥ v_t and bad parameters raise.

### Compare

`scipy.integrate.quad` on the elliptic integrand returns 1254.4379 with an estimated error of 10⁻¹¹, using adaptive subdivision near the tip; `scipy.integrate.simpson` on a uniform grid of 4001 points returns 1254.4350, the same 2 × 10⁻⁶ relative error as our own Simpson, because it is the same rule on the same grid. The library's adaptive routine is the right default for a one-off integral of a formula; our fixed-grid rules are what you use on tabulated data, where there is no formula to resample, and then the substitution has to be done on the *table* (interpolate c(y) at the θ-grid points). SymPy gives (2) and (3) exactly, in seconds for the exponential and for the θ-form of the ellipse, and not at all in reasonable time for the y-form of the ellipse: the computer algebra system needs the same hint a student does.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Three studies, each with a prediction written down first.

1. **Rules against the tip.** For n = 4, 8, …, 256 compute |I_xx − exact| for midpoint, trapezoid and Simpson on the raw elliptic integrand, and for the trapezoid rule on the θ-form. Predict the orders (2, 2, 4, and "exact") before running. Then run:

| n | midpoint | trapezoid | Simpson | trapezoid, θ |
|---|---|---|---|---|
| 4 | 57.6 | 218 | 103 | 2 × 10⁻¹³ |
| 16 | 8.26 | 28.8 | 11.8 | 2 × 10⁻¹³ |
| 64 | 1.06 | 3.65 | 1.44 | 2 × 10⁻¹³ |
| 256 | 0.134 | 0.458 | 0.179 | 2 × 10⁻¹³ |

All three raw rules fall by a factor of 2.8 = 2^(1.5) per doubling: the tip's √(s − y) defeats Simpson's fourth order as thoroughly as the trapezoid's second. The θ-form is exact at n = 4 (2 × 10⁻¹³ is rounding).

2. **Bounds against errors.** For the exponential wing, tabulate the trapezoid error and bound at n = 8, 32, 128 (§2.6) and the Simpson error and bound at n = 8, 32. Record the ratio bound/error: it should be constant (2.2 for the trapezoid, 1.5 for Simpson) because the bound uses max|f″| where the true error uses an average.

3. **The improper integral.** Compute I(R)/I(∞) for R/ℓ = 1, 2, 3, 4, 6, 8 numerically and compare with 1 − e^(−R/ℓ)(R²/2ℓ² + R/ℓ + 1): 0.0803, 0.3233, 0.5768, 0.7619, 0.9380, 0.9862. Then try f = 1/x from 1 with cutoffs 10, 100, 1000 and watch it not converge.

Raw data format: three CSV tables (`n, err_mid, err_trap, err_simp, err_theta`; `n, err, bound` per rule; `R_over_ell, I_frac`), each with a header naming the wing.

## 5. Validate and improve

| quantity (lesson's wing) | hand (§2) | code | difference | cause |
|---|---|---|---|---|
| I_xx linear | 1490.7 | 1490.72 (Simpson n = 2) | 0 | cubic integrand: exact |
| I_xx elliptic | 1254.4 | 1254.4379 (θ-trapezoid n = 4); 1254.41 (Simpson n = 1000) | 0; 0.03 | tip singularity, order 1.5 |
| I_xx exponential | 1459.1 | 1459.14 (Simpson n = 100) | 0 | smooth integrand: 10⁻⁸ |
| I_xx(∞), exponential | 95 859 | 95 857 (Simpson to 16ℓ) | 1.6 × 10⁻⁵ | the tail beyond 16ℓ, exactly e^(−16)(130) |
| trapezoid error, exponential, n = 32 | ≤ 1.04 (bound) | 0.465 | factor 2.2 | max|f″| at y = 0 against the average |
| take-off to 30 m/s | 12.58 s | 12.5756 s (Simpson n = 200) | 10⁻⁸ | — |

Every discrepancy has a named cause and the two that are not zero are predicted by the error laws. One more worth recording: Simpson with 2000 panels on ∫₁^1000 dx/x is 2 × 10⁻⁴ off, not 10⁻⁸, because h = 0.5 is coarse where 1/x changes fastest, near x = 1; a rule's accuracy depends on the integrand's scale of variation relative to h, and a wide interval with a sharp end is the classic trap. The fix is the same as for the tip: substitute x = eᵘ, which makes the integrand constant.

What the agreement does not show: that I_xx of a real aircraft is 1490 kg m². The beam model with σ = 12 kg/m² gives a wing contribution in the right range for a 1000 kg light aircraft (published I_xx values for such aircraft are 1000–1500 kg m² *in total*, fuselage and engine included, so our wing-only estimate is on the high side: real wings are lighter outboard than λ ∝ c says). A measurement, the bifilar pendulum or compound-pendulum swing done on prototypes, is what validates the number; Problem 7's λ ∝ c² law moves the estimate in the right direction.

Improvements worth making: implement Richardson extrapolation of the trapezoid rule (Romberg integration) and show it recovers fourth order on smooth integrands; add the chordwise strip inertia. Improvement worth *not* making: refining the raw elliptic integral to n = 10⁶ when four θ-panels are exact.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the hand derivations of (1)–(4) with the technique named at each step; `beam_inertia_convergence.png` with its three tables; and a one-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) stating which rule and which substitution you would use for a tabulated planform and what n the error bound demands for four figures.

### Adversarial self-checks

- **What if the span doubles?** I_xx ∝ s³ for every planform: eight times. The mass only doubles, so the radius of gyration √(I/m) doubles. What if σ doubles? I_xx doubles; the planform comparison is unchanged.
- **What if s → 0 or s → ∞?** (3) → 0 as s³ for small s (expand the exponential: the wing becomes a short rectangle) and → 4σ c_r ℓ³ for large s. Both limits are tested.
- **Which assumption fails first in the real world?** Constant σ: wing structure is heaviest at the root, and fuel is often outboard; the two errors have opposite signs and can be 20 % each.
- **What measurement could falsify the model?** A pendulum swing test of the aircraft about its roll axis: I_xx = (period/2π)² × (restoring moment per radian). If the measured total is below the wing-only beam estimate, the σ c(y) law is too heavy outboard.
- **Does the result survive a change in rule or in n?** Simpson, midpoint, trapezoid and the θ-form all converge to 1254.44 for the ellipse; only their rates differ. Simpson with n = 3 gives a `ValueError`, not a wrong number.
- **Which decision would change if the uncertainty doubled?** A simulator's roll time constant is I_xx/L_p; if σ were uncertain by 20 % instead of 10 %, the predicted roll response would be uncertain by 20 % and the aileron sizing would need a flight test to settle it. The integration technique contributes 10⁻⁶ to that uncertainty.

## Check yourself

1. Write the integration-by-parts formula and say which factor of y² e^(−y/ℓ) you differentiate and which you integrate, and why.
2. What substitution turns √(1 − y²/s²) into cos θ, and what do dy and the limits become?
3. Why is Simpson's rule exact for a cubic, and what does that imply for the linear-taper inertia?
4. The trapezoid error bound is (b − a)³ max|f″|/(12n²). By what factor does the bound fall when n doubles? When n quadruples?
5. Split 1/(T − k v²) into partial fractions and name the function that each piece integrates to.
6. Define ∫₀^∞ f dx as a limit, and say in one sentence why ∫₀^∞ y² e^(−y/ℓ) dy converges but ∫₁^∞ dx/x does not.
7. Three rules all show order 1.5 on the elliptic planform. What property of the integrand is responsible, and what fixes it?
8. The take-off time to v_t is infinite. Is that a flaw in the model or a correct statement about it? What does v_t mean physically?
