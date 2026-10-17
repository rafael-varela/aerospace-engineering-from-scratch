# P02.13  Calculus of Variations

> A function has a minimum where its derivative vanishes; a functional has one where its first variation vanishes, and the Euler–Lagrange equation is what "vanishes" means.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the Euler–Lagrange equation from the first variation of a functional, derive the Beltrami first integral for integrands with no explicit x dependence, and state Hamilton's principle. 2. Solve the brachistochrone: show the extremal is a cycloid, fit it through given end points, and give its travel time in closed form. 3. Implement and test, from scratch, a symbolic Euler–Lagrange operator, a path-time integrator for any speed field, the cycloid, Snell's invariant, a second-variation check, and a direct Ritz minimisation that approaches the exact minimum from above. 4. Decide, from the second-order behaviour of the functional and from the competing paths, that the cycloid is the minimum and not merely a stationary point, and say which assumption of the frictionless model fails first. |
| **Prerequisites** | `P02.06` (the pendulum equation, which Hamilton's principle reproduces), `P02.02` (second-order ODEs), `P01.09` (partial derivatives, the chain rule), `P01.04` (integration by parts: the whole method is one integration by parts), `P01.02` (optimisation of functions, the second-derivative test), `P00.09` (SymPy), `P00.06` (`scipy.optimize`), `P00.03`, `P00.07`. Skill check: you can integrate by parts, differentiate √(1 + y′²), and solve a transcendental equation by bracketing. |
| **Engineering question** | A frictionless parts chute must carry a component from one deck to a point 4.0 m along and 2.0 m below: what track profile delivers it fastest, how much faster is it than the straight ramp an engineer would draw first, and how do we know the answer is a minimum rather than just a stationary point? |
| **Assumptions and limits** | A point mass, frictionless, released from rest, in uniform gravity, constrained to a smooth planar track: speed depends only on depth, v = √(2gy). No rolling (a ball's rotational energy changes the effective g), no air drag, no track roughness. The Euler–Lagrange equation gives stationary paths; minimality is checked separately (second variation and competing paths), not assumed. The Ritz family (sines added to the ramp) cannot represent the cycloid's vertical start exactly, so it converges slowly. |
| **Analytical work** | First variation and the fundamental lemma; the Euler–Lagrange equation; Beltrami's identity; the shortest path as a check; the brachistochrone integrand, its first integral y(1 + y′²) = 2k, and the parametric cycloid; the closed-form time θ_f√(k/g); the ramp time; Snell's law as the first integral of a stratified Fermat problem; Hamilton's principle and the pendulum and spring–mass equations; the second variation as the reason the time increase is quadratic in the perturbation. |
| **Computational deliverable** | `starter/main.py`: `euler_lagrange`, `beltrami_first_integral`, `lagrangian_equation`, `path_time`, `gravity_speed`, `travel_time`, `straight_line_time`, `cycloid_parameter`, `brachistochrone`, `brachistochrone_time`, `snell_invariant`, `ritz_path`, `ritz_minimum_time`, `perturbed_time_increase`. Tests prove: the symbolic Euler–Lagrange operator gives y″ = 0 for the shortest path, y″ + ω²y = 0 for the oscillator and y″ = x for a forced case, and an exact extremal substituted back gives zero; the Beltrami integral of the brachistochrone integrand is y(1 + y′²) and the cycloid satisfies it; Hamilton's principle reproduces the pendulum and spring–mass equations; the path-time integrator is exact for the ramp and second-order accurate for the cycloid; the cycloid meets its end points to 10⁻¹⁰ and its numerical time matches θ_f√(k/g); the cycloid beats power-law paths and twenty random admissible paths; the time increase under a perturbation is second order on the cycloid and first order on the ramp; Snell's invariant is constant on the cycloid to 0.2 % and not on the ramp; the Ritz minimum lies above the cycloid time, improves with every added mode and comes within 2 %. |
| **Practical activity** | **computational experiment**: a Ritz convergence study, a second-variation measurement, and a sweep of end-point geometry showing when the fastest chute dips below its destination. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `brachistochrone.png`, the table of path times (cycloid, ramp, arc, Ritz), the perturbation table with its fitted exponents, and a short note stating why the cycloid is a minimum and which assumption you trust least. |
| **Assessment** | `assessments/phase-02/grading-rubric.md`; competencies 1–4. Pass: the derivation of the Euler–Lagrange equation and of the cycloid from Beltrami's identity is complete; all public tests pass; the hand values of k, θ_f and the two times agree with the code to four figures; the note gives the exponent of the perturbation law with a number. |
| **References** | `@mit-ocw-16-323-principles-of-optimal-control` Lecture 5: Calculus of variations (required-free; functionals, first variation, Euler–Lagrange); Lecture 2: Nonlinear optimization: constrained nonlinear optimization, Lagrange multipliers (required-free; the finite-dimensional warm-up); Lecture 6: Calculus of variations applied to optimal control (supplementary-free; the bridge to trajectory optimisation); `@mit-ocw-8-033-relativity` Session 8 "Variational Calculus" (supplementary-free; a short undergraduate introduction with the brachistochrone). Hamilton's principle for mechanical systems is original instruction in this lesson. |
| **Extensions** | ★ Add Coulomb friction μ to the chute (the speed at depth y becomes √(2g(y − μx)) for a path whose slope stays small) and re-run the Ritz minimisation: how does the optimum move? ★ Derive the catenary from the minimum of potential energy of a hanging cable with a length constraint (a Lagrange multiplier inside the integrand, Lecture 2). ★ Set up the minimum-time climb of an aircraft in the energy-height approximation as a variational problem and compare with the "climb at the speed of maximum specific excess power" rule. |
| **Time** | unmeasured |

## 1. The engineering problem

A workshop moves small parts between two deck levels by a gravity chute: the part starts from rest at the upper deck and slides on a smooth track to a point 4.0 m along and 2.0 m below. The obvious track is a straight ramp. The question is whether a curved track is faster, by how much, and what shape it should be. This is the **brachistochrone** ("shortest time") problem, posed by Johann Bernoulli in 1696, and it is the problem that created the calculus of variations: a question whose unknown is not a number but a *function*, the whole profile y(x), and whose answer is found by asking what happens to the travel time when the profile is nudged.

The same mathematics decides the shape of a hanging cable, the path of a light ray or a sound ray through a medium of varying speed, the equations of motion of every mechanical system (Hamilton's principle: the trajectory makes the action stationary, which is why the pendulum of `P02.06` can be derived without drawing a single force), and, in Phase 09, the minimum-time climb and minimum-fuel ascent of an aircraft or rocket. The 16.323 notes cited above go on to optimal control; this lesson stops at the point where the method is clear and the chute is solved.

Success criteria, written down first: the Euler–Lagrange operator, implemented symbolically, reproduces the three classical equations and accepts an exact extremal; the cycloid fitted through the end points reproduces them to 10⁻¹⁰ and its closed-form time matches the numerical path integral to 10⁻⁴; the straight ramp is slower by a predicted 25 %; a perturbation of amplitude ε changes the cycloid's time by an amount ∝ ε² (fitted exponent 2.0 ± 0.1) and the ramp's by ∝ ε (exponent 1.0 ± 0.15); Snell's invariant is constant along the cycloid to 0.2 %; and a direct minimisation over a six-parameter family lands within 1 % of the cycloid time from above.

## 2. Understand and derive

### 2.1 Functionals and the first variation

A **functional** assigns a number to a function: J[y] = ∫ₐᵇ F(x, y, y′) dx. The travel time of a chute, the length of a curve, the action of a mechanical system are all of this form. We want the y that makes J smallest among all functions with the same end values y(a), y(b). Take a candidate y and a **variation** η(x) with η(a) = η(b) = 0, and look at J[y + εη] as a function of the number ε (`@mit-ocw-16-323-principles-of-optimal-control` Lecture 5). If y is a minimum, this ordinary function of ε has a minimum at ε = 0, so its derivative there, the **first variation**, must vanish for every admissible η:

δJ = d/dε J[y + εη]|₀ = ∫ₐᵇ (F_y η + F_{y′} η′) dx = 0.   (1)

Integrate the second term by parts (`P01.04`): ∫ F_{y′} η′ dx = [F_{y′} η]ₐᵇ − ∫ (d/dx F_{y′}) η dx, and the boundary term dies because η vanishes at the ends. So

δJ = ∫ₐᵇ [F_y − d/dx F_{y′}] η dx = 0 for every η.

The **fundamental lemma**: if ∫ G η dx = 0 for every smooth η that vanishes at the ends, then G ≡ 0 (otherwise choose η to have the sign of G where G ≠ 0 and the integral is positive). Hence the **Euler–Lagrange equation**:

**F_y − d/dx (F_{y′}) = 0.**   (2)

It is a second-order ODE for y (the d/dx acts on F_{y′}(x, y, y′) through all three arguments). Its solutions are **extremals**: paths on which J is stationary. Whether a given extremal is a minimum, a maximum or neither is a separate question, exactly as a zero derivative is a separate question from a minimum in `P01.02`; part 2.7 addresses it.

**Check on the shortest path.** The length of a curve is ∫ √(1 + y′²) dx. F_y = 0 and F_{y′} = y′/√(1 + y′²), so (2) says d/dx[y′/√(1 + y′²)] = 0: y′ is constant, the path is a straight line. The SymPy operator returns −y″/(1 + y′²)^{3/2}, zero exactly when y″ = 0.

### 2.2 Beltrami's identity

When F has no explicit x dependence, F = F(y, y′), there is a first integral. Differentiate F − y′F_{y′} along a solution:

d/dx (F − y′F_{y′}) = F_y y′ + F_{y′} y″ − y″F_{y′} − y′ d/dx(F_{y′}) = y′[F_y − d/dx F_{y′}] = 0 by (2).

**F − y′ F_{y′} = C (constant) on every extremal.**   (3)

This is the conservation of energy of mechanics in disguise (with x as time, F as the Lagrangian, F − y′F_{y′} is minus the Hamiltonian), and it turns a second-order equation into a first-order one. `beltrami_first_integral` computes (3) and refuses an integrand that contains x, for which (3) is false.

### 2.3 The brachistochrone

Measure y *downward* from the start. A frictionless particle released from rest at y = 0 has, by energy conservation, speed v = √(2gy) at depth y, independent of the path. The time along an element ds = √(1 + y′²) dx is ds/v, so

**T[y] = ∫₀^{x_f} √[(1 + y′²)/(2gy)] dx.**   (4)

F has no x in it, so Beltrami applies. With F = √(1 + y′²)/√(2gy), F_{y′} = y′/[√(2gy)√(1 + y′²)], and

F − y′F_{y′} = [(1 + y′²) − y′²]/[√(2gy)√(1 + y′²)] = 1/√(2gy(1 + y′²)) = C,

**y (1 + y′²) = 2k** (writing 1/(2gC²) = 2k).   (5)

SymPy confirms it: `1/(2 g belt**2)` simplifies to y(1 + y′²). Solve (5) for y′ = √((2k − y)/y) and substitute y = k(1 − cos θ) = 2k sin²(θ/2): then 2k − y = 2k cos²(θ/2), y′ = cot(θ/2), and dx = dy/y′ = k sin θ dθ tan(θ/2) = 2k sin²(θ/2) dθ = k(1 − cos θ) dθ. Integrating with x(0) = 0:

**x = k(θ − sin θ), y = k(1 − cos θ): a cycloid,**   (6)

the curve traced by a point on the rim of a wheel of radius k rolling along the underside of the x axis. It starts vertically (y′ = cot(θ/2) → ∞ at θ = 0), which is the first surprise: the fastest chute begins with a free-fall drop, to buy speed early.

**Fitting the end point.** Dividing the two equations of (6), x_f/y_f = (θ_f − sin θ_f)/(1 − cos θ_f), a transcendental equation for θ_f on (0, 2π) that `cycloid_parameter` solves by bracketing (`P01.02`); then k = y_f/(1 − cos θ_f). For the chute, x_f/y_f = 2 gives **θ_f = 3.5084 rad (201.0°)** and **k = 1.0344 m**. Because θ_f > π, the track passes *below* its destination: the lowest point is y = 2k = 2.069 m at x = πk = 3.250 m, and the part climbs 7 cm in the last 0.75 m. (The end point sits exactly at the bottom when x_f/y_f = π/2; the test checks that case gives θ_f = π to 10⁻⁹.)

**The time.** Along the cycloid ds = √(dx² + dy²) = k√((1 − cos θ)² + sin²θ) dθ = 2k sin(θ/2) dθ and v = √(2gy) = 2√(gk) sin(θ/2), so ds/v = √(k/g) dθ, the same for every element:

**T = θ_f √(k/g) = 3.5084 × √(1.0344/9.80665) = 1.1394 s.**   (7)

The straight ramp of length L = √(x_f² + y_f²) = 4.472 m has constant acceleration g y_f/L along it, so T = √(2L/(g y_f/L)) = L√(2/(g y_f)) = **1.4281 s**: 25.3 % slower. A circular arc that starts vertically takes 1.1744 s, 3.1 % slower than the cycloid: close, because it has the right start, but wrong.

### 2.4 Snell's law is the same first integral

Fermat's principle says a ray takes the path of least time through a medium of speed c(x, y): T = ∫ ds/c. For a medium stratified in y, F = √(1 + y′²)/c(y) has no x, and Beltrami gives F − y′F_{y′} = 1/[c√(1 + y′²)] = constant. Since 1/√(1 + y′²) = dx/ds = sin θ, where θ is the angle of the path from the y axis (the direction of stratification),

**sin θ / c(y) = constant: Snell's law in differential form.**   (8)

The chute is a Fermat problem with c = √(2gy), so (8) holds on the cycloid: sin θ = dx/ds = sin(θ/2) in the parametrisation of (6), c = 2√(gk) sin(θ/2), and the ratio is 1/√(4gk) = 0.15699 s/m for the chute. `snell_invariant` evaluates sin θ/c at each segment midpoint of any polygonal path; on the cycloid it returns 0.15695 with a spread of 2 × 10⁻³ after the first segment (the chord direction is a poor estimate of the tangent where the curve is vertical), and on the ramp it varies by a factor of three. Any path that violates (8) is not an extremal; the test exploits that.

### 2.5 Hamilton's principle

For a mechanical system with coordinate q(t), kinetic energy T and potential energy V, define the **Lagrangian** L = T − V and the **action** S = ∫ L dt. **Hamilton's principle** states that the motion between two configurations at two times makes S stationary. The Euler–Lagrange equation (2) with x → t and y → q is then

**d/dt (∂L/∂q̇) − ∂L/∂q = 0,**   (9)

**Lagrange's equation** of motion. For the pendulum of `P02.06`, T = ½ml²θ̇², V = −mgl cos θ, and (9) gives ml²θ̈ + mgl sin θ = 0, i.e. θ̈ + (g/l) sin θ = 0, with no free-body diagram and no constraint force: the rod's tension never appears because the coordinate θ already respects the constraint. For a spring–mass, T = ½mq̇², V = ½kq², and (9) gives mq̈ + kq = 0. `lagrangian_equation` does (9) in SymPy and the tests check both. With several coordinates the same principle gives one equation per coordinate, which is how the attitude and flexible-body equations of Phase 05 are assembled; and with a constraint added through a multiplier (Lecture 2 of 16.323) the multiplier *is* the constraint force. The variational statement is not just elegant: it is the only practical way to derive equations for systems with dozens of coupled coordinates.

### 2.6 Direct minimisation: Ritz

Instead of solving (2), restrict y to a family that satisfies the end conditions, y(x) = y_f x/x_f + Σ_k a_k sin(kπx/x_f), and minimise T over the coefficients with an ordinary optimiser (`scipy.optimize.minimize`, Nelder–Mead, from a = 0: the ramp). Because the family is a subset of all admissible paths, the Ritz minimum is an **upper bound** on the true minimum and must decrease as modes are added. For the chute: 1.1836 s (1 mode, 3.9 % above the cycloid), 1.1684 (2, 2.5 %), 1.1545 (4, 1.3 %), 1.1496 (6, 0.9 %). The convergence is slow because the cycloid leaves the origin vertically (y ∝ x^{2/3}) and a finite sum of sines, whose slope at x = 0 is finite, cannot reproduce that; the Ritz path in the figure starts too flat and then overshoots. This is the general lesson of direct methods: they are only as good as the family, and the family should contain the known singular behaviour. Phase 09's trajectory optimisers are Ritz methods with better families.

### 2.7 Is it a minimum? The second variation

The Euler–Lagrange equation guarantees δJ = 0, not a minimum. Expand J[y + εη] = J[y] + ε δJ + ½ε² δ²J + …: on an extremal the linear term is gone and the change is ½ε²δ²J, *quadratic* in ε, with the sign of the **second variation** δ²J. On a non-extremal the change is linear in ε. `perturbed_time_increase` measures this with η = sin(πx/x_f): on the cycloid, ε = 0.1, 0.05, 0.025 m give ΔT = 9.46 × 10⁻⁴, 2.40 × 10⁻⁴, 6.06 × 10⁻⁵ s (ratios 3.94 and 3.97: exponent 2.0) and ΔT > 0 (δ²J > 0: a minimum along this direction); on the ramp the same perturbations give ΔT = −0.059, −0.032, −0.016 s (ratios 1.9 and 1.9: exponent 1.0, and negative, because sagging the ramp makes it faster). A positive second variation in one direction is necessary, not sufficient; for the brachistochrone the full argument (Jacobi's condition: no conjugate point before the end, which holds for θ_f < 2π) is classical and we take it as given. What we *can* check cheaply is that no competing path beats the cycloid: power-law profiles, the circular arc, twenty random admissible paths, and every Ritz family all come out slower, which is the test.

### 2.8 Where the model stops being valid

- **Frictionless.** Real chutes have μ ≈ 0.1–0.3; friction depends on the path (the normal force varies with curvature and slope), so the speed is no longer a function of depth alone and (5) fails. The vertical start is then wrong too: a part on a vertical wall does not slide at all.
- **Point mass.** A rolling ball or cylinder carries rotational kinetic energy: v² = 2gy/(1 + I/(mr²)), which rescales g by 5/7 for a solid sphere and leaves the *shape* of the optimum unchanged (the test on scaling g shows the time scales as g^{−1/2}).
- **Released from rest at the start.** With an initial speed v₀ the cycloid's cusp moves to a virtual start v₀²/(2g) above the real one.
- **Planar, smooth track.** No lateral forces, no bumps.
- **Stationary ≠ minimal.** Always check.

## 3. Build the computational model

Open `starter/main.py`. `chute` and `main` are given.

### Step 1: the symbolic operators

```python
def euler_lagrange(F, y, x):
    yp = sp.Derivative(y, x)
    return sp.simplify(sp.diff(F, y) - sp.diff(sp.diff(F, yp), x))       # (2)
```

`y` is an applied function `sp.Function("y")(x)`, and differentiating with respect to `sp.Derivative(y, x)` is legal in SymPy. `beltrami_first_integral` returns F − y′F_{y′} after substituting placeholder symbols to test whether F still contains x. `lagrangian_equation(T, V, q, t)` is (9) with L = T − V.

### Step 2: path times

```python
def path_time(x, y, speed):
    v = speed(x, y)
    ds = np.hypot(np.diff(x), np.diff(y))
    return float(np.sum(2.0 * ds / (v[1:] + v[:-1])))                       # exact when v^2 is linear in s
```

The weight 2/(v₁ + v₂) is the exact segment time when v² varies linearly along the segment, which is the case for v = √(2gy) with y linear on a chord; it also survives v₁ = 0 at the start. Guard against negative speeds and against a segment with zero speed at both ends. `gravity_speed(g)` returns the field √(2g·max(y, 0)); `travel_time` composes the two; `straight_line_time` is the closed form of 2.3.

### Step 3: the cycloid

`cycloid_parameter` brackets θ_f in (10⁻⁶, 2π − 10⁻⁶) with `brentq` on (θ − sin θ)/(1 − cos θ) − x_f/y_f and returns (k, θ_f). `brachistochrone` samples (6) uniformly in θ (uniform in *time*, by (7)); `brachistochrone_time` is (7); `snell_invariant` is (dx/ds)/v at segment midpoints.

### Step 4: Ritz and the second variation

`ritz_path` builds the family; `ritz_minimum_time` runs Nelder–Mead on `travel_time` with a penalty for any path that rises above the start (the particle could not get there); `perturbed_time_increase` returns T[y + ε sin(πx/x_f)] − T[y].

Running `main()` prints

```
shortest path: Euler-Lagrange gives -Derivative(y(x), (x, 2))/(Derivative(y(x), x)**2 + 1)**(3/2) = 0  (so y'' = 0: a straight line)
   so 1/(2 g C^2) = (Derivative(y(x), x)**2 + 1)*y(x) = const: y (1 + y'^2) = 2k, the cycloid's differential equation
pendulum from Hamilton's principle: g*sin(theta(t))/l + Derivative(theta(t), (t, 2)) = 0
spring-mass: k*q(t) + m*Derivative(q(t), (t, 2)) = 0
CHUTE from (0, 0) to (4.0, 2.0) m: cycloid k = 1.0344 m, theta_f = 3.5084 rad (201.0 deg), lowest point y = 2.069 m at x = 3.250 m
cycloid time: closed form 1.13943 s, path integral 1.13943 s
straight ramp: closed form 1.42809 s, path integral 1.42809 s (25.3 % slower)
circular arc (vertical start): 1.17444 s (3.1 % slower)
Snell invariant sin(theta)/v on the cycloid: mean 0.156951, spread 9.0e-03, prediction 1/sqrt(4 g k) = 0.156988
Ritz with 6 sine modes: T = 1.14960 s, excess over the cycloid 0.892 %, a1 = 0.9930 m
perturb the cycloid by eps = 0.1 m: T rises by 9.460e-04 s
```

and saves `brachistochrone.png`.

### Verify

- **Analytic cases:** y″ = 0 for the shortest path; y″ + ω²y = 0 for the oscillator; y″ = x for F = y′² + 2xy, with y = x³/6 substituted back giving exactly zero; the Beltrami integral of (4) equals y(1 + y′²) and the parametric cycloid satisfies it symbolically; the pendulum, spring–mass and free-particle equations from (9); the ramp time exact for any number of points; free fall (x_f = 0) gives √(2y_f/g); the end point at the cycloid's bottom gives θ_f = π.
- **Limiting cases and scaling:** doubling both end coordinates doubles k, keeps θ_f and multiplies the time by √2; quadrupling g halves the ramp time; a nearly vertical drop gives θ_f → 0 and the free-fall time; the speed field is zero (not NaN) above the start.
- **Convergence:** the cycloid's path-integral time exceeds the closed form by an amount that falls four-fold per doubling of points (second order); the perturbation law has exponent 2.0 ± 0.1 on the cycloid and 1.0 ± 0.15 on the ramp.
- **Optimality:** the cycloid beats five power-law paths and twenty random admissible paths; the Ritz minima for 1, 2, 4 modes lie strictly between the cycloid and the ramp and decrease monotonically, finishing within 2 %; Snell's invariant is constant on the cycloid to 0.2 % and varies by more than 50 % on the ramp.
- **Failure cases:** `beltrami_first_integral` raises for an integrand with explicit x; `path_time` raises for a negative speed or a stationary segment; `cycloid_parameter` raises for a non-positive end point; `straight_line_time` raises for y_f ≤ 0.

### Compare

There is no library brachistochrone; what SciPy offers is the general optimiser we used for Ritz (`scipy.optimize.minimize`) and, for problems with more structure, `scipy.integrate.solve_bvp` to solve the Euler–Lagrange equation as a boundary-value problem. The division of labour is the point: the calculus of variations turns "find the best function" into an ODE (2) or a first integral (3), which we solved exactly here; the optimiser attacks the functional directly and gives an upper bound that improves with the trial family. In Phase 09 both appear together: the necessary conditions of optimal control (`@mit-ocw-16-323-principles-of-optimal-control` Lecture 6 is the bridge) supply the structure, and a direct transcription supplies the numbers. Writing both for the chute, and watching the Ritz bound creep down toward 1.1394 s, is the calibration for believing either in a problem with no closed form.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Three studies; write each prediction first.

1. **Ritz convergence.** For the chute, record the Ritz time and the excess over (7) for 1, 2, 3, 4, 6, 8 modes. Predict the trend (monotone decrease) and whether the excess halves per added mode (it will not: the vertical start limits the rate). Then replace the sine family by y = y_f x/x_f + a₁ x^{2/3}(1 − x/x_f) + a₂ sin(πx/x_f) and see how much two parameters with the right singularity buy.
2. **Second variation.** Tabulate ΔT for ε = 0.2, 0.1, 0.05, 0.025, 0.0125 m on the cycloid, the ramp and the circular arc; fit the exponent for each. Predict 2, 1, 1.
3. **Geometry sweep.** For y_f = 2 m and x_f = 1, 2, π, 4, 6, 8 m: θ_f, k, whether the track dips below the destination (θ_f > π ⇔ x_f/y_f > π/2), and the ratio of the ramp time to the cycloid time. Predict that the ratio grows with x_f/y_f (the longer the flat run, the more the early drop pays).

Raw data format: three CSV tables with a one-line header naming the fixed parameters; `brachistochrone.png`.

What you should find: the sine-family excess falls from 3.9 % to 0.67 % at 8 modes while the x^{2/3} family reaches 0.04 % with two parameters; exponents 2.0, 1.0 and 1.0 (the arc gives 1.05 and 1.03 over successive halvings); dips for x_f = 4, 6, 8 m with x_f = π m exactly on the boundary; and a ramp-to-cycloid ratio of 1.03 at x_f = 1 m, 1.10 at 2 m, 1.25 at 4 m and 1.54 at 8 m.

## 5. Validate and improve

| quantity | hand (§2) | code | difference | cause |
|---|---|---|---|---|
| θ_f, k | 3.508 rad, 1.034 m | 3.5084 rad, 1.0344 m | 0 | same equation, bracketing to 10⁻¹⁴ |
| cycloid time | 1.1394 s (7) | 1.13943 s (path integral, 401 points) | 9 × 10⁻⁷ s | chords are slightly slower than arcs, O(h²) |
| ramp time | 1.4281 s | 1.42809 s | 0 | exact segment rule |
| ramp / cycloid | 1.253 | 1.2534 | 0 | — |
| Snell invariant | 0.15699 s/m | 0.15695 ± 0.0025 | 0.03 % mean | chord direction near the vertical start |
| ΔT exponent on cycloid | 2 | 1.98–2.0 | — | quadratic term of the expansion |
| ΔT exponent on ramp | 1 | 0.9–1.0 | — | linear term |
| Ritz, 6 modes | > 1.1394 s | 1.1496 s (+0.9 %) | upper bound, as required | the family cannot start vertically |

Every difference is named and every one is a discretisation or a family limitation, not a disagreement about the answer. What the agreement does not show is what a real chute does. A part on a real track has μ ≈ 0.2; the work against friction is roughly μmg x_f, which here is 40 % of the available mg y_f, so the final speed falls by about 23 %, and the vertical start becomes impossible (no normal force, no motion); the practical "fast chute" is a cycloid-like curve with the first 20 cm replaced by a steep straight section at the angle of friction, and the difference from the straight ramp shrinks from 25 % to perhaps 10 %. The model's value is that it tells us the *direction* of the improvement (drop early, run flat late) and the *bound* on it (25 % is the most a profile change can buy), which is what a designer needs before deciding whether a curved track is worth its cost.

Improvements worth making: Coulomb friction in `gravity_speed` (then the speed depends on the path and `path_time` must be marched along it); the x^{2/3} Ritz family; `solve_bvp` on the Euler–Lagrange equation of the friction problem, for which there is no first integral. Improvement worth *not* making: more sine modes. Eight modes buy 0.2 % over six, and the family is the wrong one.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivation of (2), (3), (5) and (6); `brachistochrone.png`; the three tables of part 4; and a one-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) giving the fastest-track time, the saving over the ramp, the evidence that the cycloid is a minimum, and the assumption you trust least.

### Adversarial self-checks

- **What if the drop doubles?** With y_f = 4 m and x_f = 4 m, x_f/y_f = 1 gives θ_f = 2.412 rad and no dip; the ramp is only 9.5 % slower than the cycloid. With y_f = 0.5 m the ratio is 8, θ_f = 4.97 rad, the track dips to 1.35 m, and the ramp is 98 % slower. The saving depends on the geometry, not on g. Check both with the code.
- **What if g → 0? x_f → 0?** Times scale as g^{−1/2}: with no gravity nothing moves. x_f → 0 gives θ_f → 0, k → ∞ and the free-fall time √(2y_f/g) = 0.639 s for 2 m.
- **Which assumption fails first in the real world?** Friction, at the vertical start, where the normal force is zero and a real part would not move at all.
- **What measurement could falsify the model?** Two tracks, a cycloid and a ramp, between the same points, a steel ball and a stopwatch (or two light gates): the predicted ratio of times is 1.253 for a sliding point and the same ratio for a rolling ball (g rescales by 5/7 for both). A measured ratio well below 1.25 means friction or rolling resistance differs between the tracks.
- **Does the result survive a change of solver or discretisation?** The closed form (7), the 401-point path integral, and the Beltrami check agree to 10⁻⁶; the Ritz bound is above all of them for every family; Nelder–Mead from the ramp and from a random start give the same coefficients to four figures (try it).
- **Which decision would change if the uncertainty doubled?** The decision is "curved track or straight?". If the friction coefficient were uncertain by ±0.1 instead of ±0.05, the predicted saving of a curved track could range from 5 % to 15 %, and the decision would hinge on a measured μ for the actual part and track material before any track is bent.

## Check yourself

1. Write the Euler–Lagrange equation for J = ∫ F(x, y, y′) dx and say what integration by parts contributed to its derivation.
2. What does the fundamental lemma of the calculus of variations assert, and why is it needed?
3. State Beltrami's identity and the condition on F under which it holds. What conserved quantity of mechanics is it?
4. Starting from y(1 + y′²) = 2k, show that y = k(1 − cos θ), x = k(θ − sin θ) satisfies it.
5. Why does the brachistochrone start vertically? Give the energy argument in one sentence.
6. A chute goes 6 m along and 2 m down. Does the fastest track dip below its destination? Give the criterion.
7. The travel time changes by 10⁻³ s when the cycloid is perturbed by 0.1 m. By about how much does it change for a 0.05 m perturbation, and what would the answer be for the straight ramp?
8. Write the Lagrangian of a simple pendulum and derive its equation of motion from Lagrange's equation in two lines.
