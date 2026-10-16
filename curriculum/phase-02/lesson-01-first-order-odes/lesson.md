# P02.01  First-Order ODEs

> A first-order differential equation tells you the slope at every point. Solving it is following the slope; the art is knowing when you can write the answer down and when you must step.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the vertical rocket equation from momentum conservation, solve it by separation of variables, and integrate it by parts to get the altitude. 2. Solve a first-order linear equation with an integrating factor and a separable nonlinear one (quadratic drag), and state where each solution stops being valid. 3. Implement and test Euler and RK4 steppers from scratch, measure their orders of convergence, and compare them with `scipy.integrate.solve_ivp`. 4. Read a direction field and a phase line: locate equilibria of an autonomous equation, classify their stability, and identify terminal velocity as one. 5. Decide how much altitude drag costs a sounding rocket and which modelling assumption limits the answer. |
| **Prerequisites** | `P01.03` (the Fundamental Theorem), `P01.04` (integration by parts, partial fractions, `∫ sec²`), `P01.01` (orders of convergence, Taylor remainder), `P00.03` (NumPy), `P00.06` (`solve_ivp`), `P00.07` (pytest). Skill check: you can integrate ln(a/u) by parts, write 1/(a² − v²) as partial fractions, and differentiate tan(α − βt). |
| **Engineering question** | A 250 kg sounding rocket burns 150 kg of propellant in 20 s at 2200 m/s exhaust velocity. How fast and how high is it at burnout, how high does it coast, and how much of the vacuum answer does atmospheric drag take away? |
| **Assumptions and limits** | Vertical flight only (no pitch, no wind); constant exhaust velocity and mass flow; constant g = 9.81 m/s² (0.7 % low at 45 km: an error we accept and quantify); exponential atmosphere ρ = 1.225 e^(−h/8500 m); constant drag coefficient C_D = 0.5, which is wrong at transonic speeds (C_D roughly doubles near Mach 1) and makes the drag loss the least trustworthy number here; rigid rocket, no thrust tail-off. The closed forms also assume constant mass (coast) or constant density (quadratic-drag coast) and are used as brackets, not answers. |
| **Analytical work** | Separation of variables for dm/dt = −ṁ and for the rocket equation; integration by parts for the burnout altitude; the integrating factor for dv/dt + (c/m)v = −g; separation and `∫ dv/(g + kv²)` for the quadratic-drag coast; the Euler step as a truncated Taylor series with its O(h²) local and O(h) global error; existence and uniqueness (the Lipschitz condition) and the singularity at m → 0; the phase line of dv/dt = −g − kv|v| and the stability test f′(v*) < 0. |
| **Computational deliverable** | `starter/main.py`: `burn_time`, `mass_at`, `tsiolkovsky_delta_v`, `vacuum_velocity`, `vacuum_altitude`, `apogee_from_burnout`, `air_density`, `rocket_rhs`, `euler_step`, `rk4_step`, `integrate_fixed_step`, `simulate_rocket`, `convergence_study`, `quadratic_drag_coast`, `quadratic_drag_velocity`, `linear_drag_velocity`, `terminal_velocity`, `equilibria`, `direction_field`. Tests prove: every closed form satisfies its differential equation by finite differences; the altitude is the integral of the velocity; Euler converges at order 1 and RK4 at order 4; the drag model agrees with `solve_ivp` to 10⁻⁴ m/s; the quadratic-drag closed form matches RK4 and reduces continuously to the vacuum case; the linear-drag solution approaches terminal velocity from any start; terminal velocity is the single, stable equilibrium of the autonomous equation; direction-field arrows are unit vectors along the slope. |
| **Practical activity** | **computational experiment**: a convergence study of Euler against RK4 on the vacuum burn, a drag-coefficient sweep, a burn-rate sweep that exposes the gravity-loss/drag-loss trade, and a scale-height sensitivity. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. This lesson models a rocket; it does not build or fire one (see `CONTRIBUTING.md`). |
| **Evidence** | Completed `main.py`, passing tests, `rocket_ascent.png`, the hand derivations of the burnout velocity and altitude, the four sweep tables with your written predictions, and a short note on the drag loss and its uncertainty. |
| **Assessment** | `assessments/phase-02/grading-rubric.md`; competencies 1–4. Pass: the derivation of the vacuum altitude is complete including the integration by parts; all public tests pass; your hand values for burnout agree with the code to three figures; the note states the apogee with drag, its sensitivity to the scale height, and the assumption you trust least. |
| **References** | `@openstax-calculus-volume-2` §4.1 "Basics of Differential Equations", §4.2 "Direction Fields and Numerical Methods", §4.3 "Separable Equations", §4.5 "First-Order Linear Equations" (required-free); §4.4 "The Logistic Equation" (supplementary-free, for the phase line). `@mit-ocw-18-03-differential-equations` Lecture 1 "The Geometrical View of y' = f(x,y)", Lecture 3 "Solving First-order Linear ODEs", Lecture 5 "First-order Autonomous ODEs" (required-free). |
| **Extensions** | ★ Replace the constant C_D by a Mach-dependent table (a transonic bump) and quantify how much the apogee moves. ★ Add a thrust tail-off over the last second of the burn and show the integrator's order is restored only if the tail-off is smooth. ★ Replace constant g by g₀(R/(R+h))² and compare the apogee. |
| **Time** | unmeasured |

## 1. The engineering problem

A sounding rocket lifts a 100 kg payload-and-structure (the **dry mass**) with 150 kg of propellant. The motor pushes 7.5 kg of gas out of the nozzle every second at 2200 m/s, so it burns for 20 s and produces 16.5 kN of **thrust**, 6.7 times the lift-off weight. Everything a rocket does is a rate of change: mass falls at a known rate, velocity changes at a rate set by thrust, drag and weight, and altitude changes at the rate of velocity. None of these quantities is known as a function of time in advance; what is known is the **differential equation**, the rule linking each rate to the current state. This lesson is about turning that rule into numbers.

The engineering question has three parts. How fast and how high is the rocket at **burnout**, the instant the propellant runs out? How much higher does it **coast** before gravity and drag stop it? And how much of the answer is lost to drag, which in a vacuum would be zero? The first part has a closed form, the famous rocket equation, if you ignore drag. The coast has a closed form if you freeze the air density. The full problem, with density falling by a factor of e every 8.5 km and drag growing with v², has none: you must step through it numerically, and then you must know how far to trust the steps.

Success criteria, written down first: the closed-form burnout velocity and altitude satisfy their differential equations to 10⁻⁷ by finite differences; the RK4 stepper reproduces them to 10⁻⁵ m/s with 400 steps and shows fourth-order convergence (error ÷ 16 per halving of the step) while Euler shows first order; the full drag model agrees with `solve_ivp` at rtol = 10⁻¹⁰ to within 10⁻⁴ m/s at burnout; and the frozen-density coast formula and the vacuum formula bracket the full-model apogee.

## 2. Understand and derive

### 2.1 What a first-order ODE is, and what solving it means

A **first-order ordinary differential equation** is a rule dy/dt = f(t, y): given the time and the current value, it tells you the slope. A **solution** is a function y(t) whose derivative is f(t, y(t)) at every t; it is a curve that is tangent to the slope field everywhere. There are infinitely many such curves, one through each starting point, so a problem is only well posed as an **initial value problem**: dy/dt = f(t, y), y(t₀) = y₀ (`@openstax-calculus-volume-2` §4.1).

The **existence and uniqueness theorem** says when exactly one curve passes through the start: if f and ∂f/∂y are continuous in a box around (t₀, y₀), the IVP has exactly one solution in some interval around t₀. The condition on ∂f/∂y (a **Lipschitz condition**) is what rules out splitting. Without it two solutions can share a point: dy/dt = 3y^(2/3), y(0) = 0 is solved both by y = 0 and by y = t³, because ∂f/∂y = 2y^(−1/3) blows up at y = 0. Every equation in this lesson is Lipschitz wherever we use it, and the theorem also tells us where to stop: the rocket equation has a 1/(m₀ − ṁt) term that blows up at t = m₀/ṁ = 33.3 s, so no solution exists past that time. Physically the rocket runs out of propellant at 20 s and switches to a different equation; mathematically the model must be stopped before the singularity, and the code refuses to evaluate the burn formulas past burnout.

### 2.2 Mass: the simplest separable equation

**Assumption 1 (constant mass flow):** the motor expels propellant at a constant rate ṁ. Then dm/dt = −ṁ. The equation is **separable**: the variables can be put on opposite sides, dm = −ṁ dt, and each side integrated (§4.3):

**m(t) = m₀ − ṁ t,   0 ≤ t ≤ t_b = (m₀ − m_dry)/ṁ**.   (1)

For our rocket t_b = 150/7.5 = 20.0 s. After t_b the mass is constant at m_dry: the solution is piecewise, and the pieces join continuously. Units: kg/s × s = kg. Trivial, but the pattern (separate, integrate, apply the initial condition, state the interval of validity) is the one every separable equation follows.

### 2.3 The rocket equation: momentum, then separation

Take the rocket plus the gas it will expel in the next dt as one system of mass m moving upward at v. **Assumption 2 (constant exhaust velocity):** the gas leaves at speed v_e relative to the rocket, so its velocity in the ground frame is v − v_e. **Assumption 3 (vertical flight, no drag for now):** the only external force is weight, −mg. Newton's second law for the system, F_ext = d(momentum)/dt, over the interval dt gives

−mg dt = [(m − ṁ dt)(v + dv) + ṁ dt (v − v_e)] − mv = m dv − ṁ v_e dt + O(dt²),

so

**dv/dt = ṁ v_e / m(t) − g = T/m(t) − g**,   T = ṁ v_e.   (2)

The thrust is the mass flow times the exhaust velocity: 7.5 × 2200 = 16 500 N. Units: (kg/s)(m/s)/kg = m/s². Good. At lift-off the acceleration is 16 500/250 − 9.81 = 56.2 m/s² (5.7 g); at burnout, with m = 100 kg, it is 155.2 m/s² (15.8 g). The rocket accelerates hardest at the end, when it is lightest, which is why payloads are qualified for the burnout load, not the lift-off one.

Equation (2) has no v on the right-hand side, so it is separable in the most trivial way: dv = (ṁv_e/(m₀ − ṁt) − g) dt, and ∫ ṁ dt/(m₀ − ṁt) = −ln(m₀ − ṁt) (substitute u = m₀ − ṁt). With v(0) = 0:

**v(t) = v_e ln(m₀/m(t)) − g t**.   (3)

At burnout the first term is **Tsiolkovsky's rocket equation**, Δv = v_e ln(m₀/m_f) = 2200 ln 2.5 = 2015.8 m/s, and the second is the **gravity loss** g t_b = 196.2 m/s. Burnout velocity in vacuum: **v_b = 1819.6 m/s**. Two things to notice. The velocity gain depends on the mass *ratio* through a logarithm: doubling the propellant does not double Δv. And the gravity loss depends only on how *long* the burn takes: a faster burn (larger ṁ, same propellant) loses less to gravity. Part 4 tests whether that is still true with drag.

### 2.4 Altitude: integrate again, by parts

Altitude is dh/dt = v, so h(t) = ∫₀ᵗ v dτ. The −gt term gives −½gt². For the logarithm, substitute u = m₀ − ṁτ, dτ = −du/ṁ, and integrate by parts (`P01.04`): ∫ ln(m₀/u) du = u ln(m₀/u) + u, because d/du[u ln(m₀/u)] = ln(m₀/u) − 1. Evaluating from u = m₀ to u = m(t):

**h(t) = v_e t − (v_e m(t)/ṁ) ln(m₀/m(t)) − ½ g t²**.   (4)

Check the limits: h(0) = 0 ✓. Differentiate back: with dm/dt = −ṁ, d/dt[m ln(m₀/m)] = −ṁ ln(m₀/m) + m(ṁ/m) = ṁ − ṁ ln(m₀/m), so dh/dt = v_e − (v_e/ṁ)(ṁ − ṁ ln(m₀/m)) − gt = v_e ln(m₀/m) − gt, which is (3) ✓ (the test does this check numerically). At burnout: h_b = 2200 × 20 − (2200 × 100/7.5) ln 2.5 − ½ × 9.81 × 400 = 44 000 − 26 878 − 1962 = **15 160 m**. Dimensions: m/s × s, (m/s × kg)/(kg/s), m/s² × s², all metres ✓.

### 2.5 Coast in vacuum: energy

After burnout, with no thrust and no drag, dv/dt = −g and the apogee follows from energy or from integrating twice: **h_apex = h_b + v_b²/(2g)** = 15 160 + 1819.6²/19.62 = 183.9 km. That is the vacuum answer and it is wildly optimistic, as part 3 will show, because this rocket leaves the dense atmosphere at Mach 3.5.

### 2.6 Drag makes the equation nonlinear

**Assumption 4 (quadratic drag):** D = ½ ρ(h) C_D A v|v|, opposing the motion. **Assumption 5 (exponential atmosphere):** ρ(h) = ρ₀ e^(−h/H) with ρ₀ = 1.225 kg/m³ and H = 8500 m, a one-parameter fit to the lower atmosphere that is within about 10 % up to 30 km and ignores the temperature structure. The velocity equation becomes

dv/dt = [T − ½ρ(h) C_D A v|v|]/m − g,   (5)

which is coupled to dh/dt = v through ρ(h) and is nonlinear in v. There is no closed form. This is the normal situation in engineering, and it is why the second half of this lesson is numerical. The reference area is A = π(0.15 m)² = 0.0707 m² and C_D = 0.5 is a rough supersonic value for a slender body; the real C_D of a sounding rocket rises sharply through Mach 1 and falls again above Mach 2, so the drag loss we compute is a model estimate, not a measurement.

### 2.7 Two coast models that do have closed forms

**Quadratic drag, constant density, constant mass** (the coast after burnout, with ρ frozen). Upward, dv/dt = −g − kv² with k = ρC_D A/(2m) (units 1/m). Separate: dv/(g + kv²) = −dt. With v = √(g/k) tan θ the integral is ∫ dv/(g + kv²) = (1/√(gk)) arctan(v√(k/g)), so

**v(t) = √(g/k) tan(arctan(v_b√(k/g)) − √(gk) t)**,   (6)

valid until v reaches zero at **t_apex = arctan(v_b√(k/g))/√(gk)**. The height gained follows from the chain rule, dh = v dt = v dv/(dv/dt) = −v dv/(g + kv²), integrated from v_b to 0: **Δh = ln(1 + kv_b²/g)/(2k)**. As k → 0, arctan x ≈ x and ln(1 + x) ≈ x recover t = v_b/g and Δh = v_b²/(2g): the vacuum case is the limit, and the test checks the limit is approached continuously. The tangent in (6) also shows a failure mode: the formula is for the *upward* leg only, because on the way down the drag changes sign. Never use a separable solution outside the interval where you derived it.

**Linear drag, constant mass** (a parachute-like regime, or a first model of any speed-proportional resistance): dv/dt + (c/m) v = −g. This is a **first-order linear equation** (§4.5), and the tool is the **integrating factor** μ(t) = e^(∫ c/m dt) = e^(ct/m). Multiply through: d/dt[v e^(ct/m)] = −g e^(ct/m). Integrate: v e^(ct/m) = −(mg/c) e^(ct/m) + C, and with v(0) = v₀,

**v(t) = (v₀ + mg/c) e^(−ct/m) − mg/c**.   (7)

Whatever v₀ is, v → −mg/c as t → ∞: the **terminal velocity**, approached with time constant m/c. For m = 100 kg and c = 2 kg/s the time constant is 50 s, the terminal velocity is −490.5 m/s, and from v₀ = 300 m/s the velocity after one time constant is −199.7 m/s. The integrating factor is the method for *every* first-order linear equation, including ones with time-varying coefficients, which is why 18.03 Lecture 3 spends a whole hour on it.

### 2.8 Autonomous equations and the phase line

For a falling body with quadratic drag the velocity equation is dv/dt = −g − kv|v| = f(v): no explicit t. Such an equation is **autonomous**, and its qualitative behaviour can be read without solving it (18.03 Lecture 5, `@openstax-calculus-volume-2` §4.4). Where f(v*) = 0 the solution is constant: an **equilibrium**. Here f = 0 at v* = −√(g/k), the terminal velocity, 212.9 m/s for the empty rocket at sea level. For v below v* (falling faster), f > 0 and v increases toward v*; for v above v*, f < 0 and v decreases toward v*. The equilibrium is **stable**, and the test for stability is the sign of the slope: f′(v*) = −2k|v*| < 0. The **phase line** is the v axis with arrows: all of them point at v*. Two warnings. The phase line says *where* solutions go, not *how fast*; and the rising leg of the same equation (v > 0, where f = −g − kv² < 0 everywhere) has no equilibrium at all: the rocket cannot hover.

### 2.9 Euler's method, and why we need better

When there is no closed form we follow the slope in small steps. The Taylor expansion y(t + h) = y(t) + h y′(t) + ½h² y″(t) + … with y′ = f gives **Euler's method** (§4.2): y_{n+1} = y_n + h f(t_n, y_n). The dropped term is ½h² y″, the **local error**, O(h²). Over a fixed interval there are (t_end − t₀)/h steps, so the **global error** is O(h): halve the step and the error halves. First order.

On the vacuum burn (3), Euler with 10 steps (h = 2 s) is 95.5 m/s off at burnout; 20 steps, 48.6; 40 steps, 24.5; 80, 12.3; 160, 6.17. The ratios are 1.96, 1.98, 1.99, 2.00: order one, and the error constant, error/h, settles to 49.3 m/s per second of step. That number is predictable: when the right-hand side depends on t alone, Euler's method is the left Riemann sum of ∫f dt, whose error is ½h[f(t_b) − f(0)] = ½h(155.2 − 56.2) = 49.5h. The Taylor remainder is a prediction, not a surprise.

The **classical Runge–Kutta method** (RK4) evaluates the slope four times per step, at the start, twice at the midpoint and once at the end, and combines them with weights 1, 2, 2, 1: y_{n+1} = y_n + (h/6)(k₁ + 2k₂ + 2k₃ + k₄). The weights are chosen so that the Taylor series of the step matches the exact one through h⁴; the local error is O(h⁵) and the global error O(h⁴): halve the step and the error falls sixteen-fold. On the same burn the RK4 errors are 2.23 × 10⁻³, 1.41 × 10⁻⁴, 8.82 × 10⁻⁶, 5.52 × 10⁻⁷, 3.45 × 10⁻⁸ m/s: ratios 15.8, 16.0, 16.0, 16.0. With 10 steps RK4 is already 40 000 times more accurate than Euler, at four times the cost per step. P03.05 derives the weights; here we use RK4 as the workhorse and verify its order rather than take it on trust.

One practical point the code must respect. At burnout the right-hand side jumps (thrust goes from 16.5 kN to zero). A stepper whose stages straddle the jump has an O(h) error in that one step, which drags the global order down to one however good the method. The cure is to integrate the burn and the coast as two separate problems with smooth right-hand sides, joined at t_b. `simulate_rocket` does exactly that, and the test checks that the switch-off lands on a grid point. `solve_ivp` solves the same problem with its `events` machinery; we do it by hand to see why it is needed.

### 2.10 Where the model stops being valid

- **Constant C_D.** The worst assumption. Through Mach 1 the real drag coefficient of a slender rocket roughly doubles, then falls; our rocket crosses 250–400 m/s between t = 4.2 s and 6.8 s, below 1.4 km altitude where the air is densest. Expect the drag loss to be wrong by tens of per cent.
- **Exponential atmosphere.** H = 8500 m is a compromise over the troposphere and stratosphere; part 4 shows ±1 km in H moves the apogee by +9/−6 km.
- **Constant g.** At 45 km, g is 0.7 % lower than at the surface, a 0.3 km effect on the apogee: negligible here, not negligible for an orbital launch.
- **Vertical flight.** A real sounding rocket is launched a few degrees off vertical and follows a curved path; its apogee is a little lower and its range is not zero.
- **Constant ṁ and v_e.** Solid motors have a thrust curve, not a step; part 4's extension shows what the tail-off does to the integrator.

## 3. Build the computational model

Open `starter/main.py`. `sounding_rocket()` is given; everything else is yours. The state vector is **y** = (h, v, m) and every right-hand side returns d**y**/dt in the same order.

### Step 1: the separable solutions

```python
def mass_at(t, m0, m_dry, mdot):
    t = np.asarray(t, dtype=float)
    if np.any(t < 0.0):
        raise ValueError("time must be non-negative")
    tb = burn_time(m0, m_dry, mdot)
    return np.where(t < tb, m0 - mdot * t, m_dry)                             # (1)

def vacuum_velocity(t, p):
    t = np.asarray(t, dtype=float)
    tb = burn_time(p["m0"], p["m_dry"], p["mdot"])
    if np.any(t < 0.0) or np.any(t > tb * (1.0 + 1e-12)):
        raise ValueError("valid only during the burn")
    m = p["m0"] - p["mdot"] * t
    return p["ve"] * np.log(p["m0"] / m) - p["g"] * t                         # (3)
```

`vacuum_altitude` is (4). Type it from your own derivation and let the test differentiate it numerically against `vacuum_velocity`: a dropped term in the integration by parts shows up as a residual of several m/s.

### Step 2: the right-hand side with drag

```python
def rocket_rhs(t, y, p, burning=True, drag=True):
    h, v, m = float(y[0]), float(y[1]), float(y[2])
    thrust = p["mdot"] * p["ve"] if burning else 0.0
    d = 0.5 * float(air_density(h)) * p["cd"] * p["area"] * v * abs(v) if drag else 0.0
    return np.array([v, (thrust - d) / m - p["g"], -p["mdot"] if burning else 0.0])   # (5)
```

`burning` is a flag, not a test on m, so that each phase is a smooth function (2.9). Write v|v| rather than v² so the drag opposes the motion on the way down too.

### Step 3: the steppers

```python
def euler_step(f, t, y, h):
    return y + h * f(t, y)

def rk4_step(f, t, y, h):
    k1 = f(t, y)
    k2 = f(t + 0.5 * h, y + 0.5 * h * k1)
    k3 = f(t + 0.5 * h, y + 0.5 * h * k2)
    k4 = f(t + h, y + h * k3)
    return y + h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
```

`integrate_fixed_step(f, t0, y0, t_end, n_steps, method)` allocates the grid with `np.linspace`, loops, and returns `(ts, ys)` with `ys[k]` the state at `ts[k]`. `simulate_rocket` calls it twice, once per phase, splitting `n_steps` in proportion to the phase durations, and concatenates. `convergence_study(method, ns, p)` integrates the drag-free burn with each step count and returns |v(t_b) − v_exact|.

### Step 4: the coast closed forms and the qualitative tools

`quadratic_drag_coast(v_b, k, g)` returns (t_apex, Δh) from 2.7 and must return the vacuum values for k = 0 rather than divide by zero. `linear_drag_velocity` is (7). `terminal_velocity` is √(2mg/(ρC_D A)). `equilibria(f, grid)` evaluates f on the grid, brackets each sign change, refines the root with `scipy.optimize.brentq`, and classifies it by the sign of a central-difference slope. `direction_field(f, t_grid, y_grid)` returns the unit vectors (1, f)/√(1 + f²) on the mesh for `plt.quiver`.

Running `main()` prints:

```
thrust 16.5 kN, T/W at lift-off 6.73, burn time 20.0 s
Tsiolkovsky delta-v 2015.8 m/s; gravity loss g t_b = 196.2 m/s
vacuum burnout: v_b = 1819.6 m/s, h_b = 15160 m, apogee 183.92 km
euler errors in v(t_b): [95.544 48.634 24.533 12.321  6.174]  rates [0.97 0.99 0.99 1.  ]
rk4   errors in v(t_b): [2.231e-03 1.408e-04 8.824e-06 5.519e-07 3.450e-08]  rates [3.99 4. 4. 4.]
with drag (RK4, 2400 steps): burnout v_b = 1180.9 m/s at h_b = 11570 m; drag costs 639 m/s
apogee with drag 45.56 km at t = 99.8 s
solve_ivp (RK45, rtol 1e-10): burnout v_b = 1180.9 m/s, difference from RK4 5.32e-09 m/s
coast with density frozen at 0.314 kg/m^3: k = 5.549e-05 1/m, apex after 52.7 s, height gained 19.68 km (vacuum: 120.4 s, 71.07 km)
full model gains 33.99 km in 79.8 s
terminal velocity of the empty rocket at sea level: 212.9 m/s
equilibria of dv/dt = -g - k v|v| on [-800, 800]: [(-420.447, True)]
```

### Verify

Run `python -m pytest curriculum/phase-02/lesson-01-first-order-odes`. The tests:

1. **Closed forms satisfy their equations.** Central differences of (3), (4), (6) and (7) match the right-hand sides to 10⁻⁶ relative; h(t_b) matches the trapezoid rule on v(t) to 10⁻⁷; burnout velocity equals Tsiolkovsky minus gt_b; two burns in sequence add their Δv.
2. **Limiting cases.** k → 0 in the quadratic-drag coast recovers v_b/g and v_b²/(2g) continuously; the linear-drag solution reaches −mg/c from v₀ = −500, 0 and 300 m/s after 20 time constants; c = 0 gives free fall; terminal velocity doubles when the mass quadruples.
3. **Steppers.** One Euler step of y′ = −2y has error ½(2h)² to within 30 %; one RK4 step is inside the first missing Taylor term; on the vacuum burn the global orders are 1.0 ± 0.1 and 4.0 ± 0.15; RK4 with 400 steps matches (3) and (4) to 10⁻⁵ m/s and 10⁻⁴ m.
4. **Oracle.** The full drag model agrees with `solve_ivp` (rtol 10⁻¹⁰) to 10⁻⁴ m/s at burnout and 10⁻² m along the coast; the switch-off sits on a grid point; drag lowers both burnout numbers; mass never drops below m_dry.
5. **Qualitative tools.** `equilibria` finds exactly one root of −g − kv|v|, at −v_terminal, stable; on the logistic equation it finds 0 (unstable) and 1 (stable); direction-field arrows are unit vectors with slope f and positive time component.

### Compare

`main()` hands `rocket_rhs` to `scipy.integrate.solve_ivp` with its default RK45 method and rtol = 10⁻¹⁰. It agrees with our fixed-step RK4 (h = 0.05 s) to 5 × 10⁻⁹ m/s at burnout and 10⁻⁵ m over the coast. RK45 is a Runge–Kutta pair (Dormand–Prince): it computes a fourth- and a fifth-order estimate from the same six slope evaluations, uses their difference as an error estimate, and *adapts* the step so that the estimate stays below the tolerance. On the coast, where the state changes slowly, it takes steps of several seconds; near lift-off, where the acceleration jumps, it takes small ones. The library is what you use in practice because of that adaptivity and because of `events` (to stop exactly at apogee or burnout) and dense output. We wrote our own to know what an "order of accuracy" means operationally, to be able to tell a stepping error from a modelling error, and to recognise the one failure mode the library shares with us: a discontinuous right-hand side costs any method its order unless the integration is split at the jump.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Four sweeps; write each prediction down before running.

1. **Convergence.** `convergence_study` for Euler and RK4 with n = 10, 20, 40, 80, 160 on the vacuum burn, log-log plot of error against n. Predict the slopes (−1 and −4) and the Euler error constant (≈ ½ t_b × typical d²v/dt², 2.9).
2. **Drag-coefficient sweep.** Apogee with drag for C_D = 0.3, 0.5, 0.8 (RK4, 3000 steps to 150 s): burnout velocities 1374, 1181, 979 m/s; apogees 75.1, 45.6, 27.0 km. Prediction first: is the apogee linear in C_D? (It is not; it falls faster than 1/C_D because higher drag also keeps the rocket lower, where the air is denser.)
3. **Burn-rate sweep.** Same propellant, ṁ = 5, 7.5, 15 kg/s (t_b = 30, 20, 10 s). In vacuum, burnout velocity 1722, 1820, 1918 m/s and apogee 172, 184, 196 km: the fast burn wins because gravity loss ∝ t_b. With drag: apogee 55.5, 45.6, 32.5 km, the *slow* burn wins because it reaches high speed only in thinner air. Predict which effect dominates *before* you run, then explain the reversal.
4. **Scale-height sensitivity.** Apogee for H = 7500, 8500, 9500 m: 55.0, 45.6, 39.1 km. A ±12 % uncertainty in the atmosphere model is a +21 %/−14 % uncertainty in the apogee.

Raw data format: four CSV tables (`n, err_euler, err_rk4`; `cd, v_b, apogee`; `mdot, t_b, v_b_vac, apogee_vac, v_b_drag, apogee_drag`; `H, apogee`), each with a one-line header naming the fixed parameters.

What you should find: the convergence lines are straight with the predicted slopes until the RK4 error reaches 10⁻⁸ m/s, where rounding (10⁻¹³ relative, times 1000 m/s, times a few hundred steps) begins to show; the drag sweep is strongly nonlinear; the burn-rate sweep reverses between vacuum and atmosphere; and the scale height, not the integrator, dominates the uncertainty of the apogee.

## 5. Validate and improve

| quantity | hand (§2) | code | difference | cause |
|---|---|---|---|---|
| t_b (s) | 20.0 | 20.0 | 0 | exact |
| Δv Tsiolkovsky (m/s) | 2015.8 | 2015.8 | 0 | exact |
| v_b vacuum (m/s) | 1819.6 | 1819.6 (RK4, 400 steps: within 10⁻⁵) | 0 | exact; RK4 error 3 × 10⁻⁸ at 160 steps |
| h_b vacuum (m) | 15 160 | 15 160 | 0 | exact |
| apogee vacuum (km) | 183.9 | 183.92 | 0 | exact |
| v_b with drag (m/s) | no closed form | 1180.9 (RK4) vs 1180.9 (`solve_ivp`) | 5 × 10⁻⁹ | both converged |
| apogee with drag (km) | 19.7 (ρ frozen at h_b) to 71.1 (vacuum coast from the drag burnout) | 34.0 gained, apogee 45.56 | — | brackets, as predicted |
| terminal velocity at ρ_b (m/s) | 420.4 | 420.447 (`equilibria`) | 0 | brentq to 10⁻¹² |

The numerical agreement is complete: every difference between closed form and stepper is explained by the stepper's order and rounding, and the two independent integrators agree to nine digits. What the agreement does *not* show is that the rocket flies this way. The drag loss (639 m/s at burnout, 138 km of apogee) rests on C_D = 0.5 at every Mach number and on an exponential atmosphere, and part 4 shows the apogee moving by ±20 % for a 12 % change in the scale height alone. A Mach-dependent C_D would change the burnout velocity by a comparable amount. The honest statement of the result is: apogee 46 km with the stated model, uncertain by at least ±10 km from the atmosphere and drag models, against a vacuum figure of 184 km that is a ceiling, not a prediction.

What validates such a model is a flight: a sounding rocket of this class carries a GPS or barometric altimeter and a telemetry link, and the measured apogee against the predicted one is the standard post-flight comparison. Without a flight, the next best check is the sensitivity table: a model whose answer swings by 20 % under a plausible parameter change should be reported with that swing.

Improvements worth making: read ρ(h) from the ISA layers of `P00.03` instead of one exponential; let C_D depend on Mach through a table; replace `simulate_rocket`'s two-phase split by a general event-detection step (bisect for the time at which m = m_dry inside a step, then restart). Improvement worth *not* making: a tighter tolerance on the integrator. It is already ten orders of magnitude below the modelling uncertainty.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the hand derivations of (3), (4), (6) and (7) with the interval of validity stated for each; `rocket_ascent.png`; the four sweep tables with your predictions written beside the results; and a one-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) giving the burnout state and apogee with drag, the bracket from the two closed forms, the scale-height sensitivity, and the assumption you trust least.

### Adversarial self-checks

- **What if the propellant mass doubles (m₀ = 400 kg, same ṁ)?** Burn time 40 s, Tsiolkovsky Δv = 2200 ln 4 = 3050 m/s, gravity loss 392 m/s. Check with the code that the vacuum burnout velocity is 2657 m/s and that the drag loss, 624 m/s, is *smaller* in proportion (23 % of the vacuum value against 35 % for the baseline), because the heavier rocket is slower in the dense air. (The singularity moves to m₀/ṁ = 53.3 s, still beyond t_b.)
- **What if the mass flow goes to infinity?** t_b → 0, gravity loss → 0, v_b → Tsiolkovsky's Δv, an impulsive burn. With drag, the same limit is the worst case: all the speed is reached at sea level. What does the burn-rate sweep in part 4 say about where the optimum lies?
- **Which assumption fails first in the real world?** Constant C_D through Mach 1. Second, the exponential atmosphere above 30 km.
- **What measurement could falsify the model?** A flight with an on-board accelerometer: the measured acceleration during the burn is (T − D)/m − g, and dividing out the known thrust and mass gives D(t), hence C_D(Mach). A transonic bump in that curve falsifies Assumption 4 directly.
- **Does the result survive a change of step or solver?** RK4 at 0.05 s and `solve_ivp` at rtol 10⁻¹⁰ agree to nine digits; Euler at 0.05 s does not (several hundred metres in apogee). Which decision would change if you used Euler with 120 steps and did not know its order?
- **Which decision changes if the uncertainty doubles?** If the apogee is uncertain by ±20 km instead of ±10 km, a payload that needs 40 km of altitude for its experiment moves from "probably" to "no margin", and the mission would carry more propellant or accept a 50 % chance of failing its objective.

## Check yourself

1. State the existence and uniqueness theorem for dy/dt = f(t, y) and give an equation that violates it at a point, with two solutions through that point.
2. Separate and solve dm/dt = −ṁ with m(0) = m₀, and state the interval on which the model applies.
3. Why does the burnout velocity depend on the mass *ratio* through a logarithm and not on the propellant mass linearly?
4. What is the gravity loss, what does it depend on, and why does a faster burn reduce it?
5. What is an integrating factor, and what is it for dv/dt + (c/m)v = −g?
6. The equation dv/dt = −g − kv|v| has one equilibrium. Where is it, is it stable, and what single calculation decides that?
7. Euler's method has local error O(h²) and global error O(h). Explain the difference in one sentence each, and give the corresponding orders for RK4.
8. The vacuum coast formula and the frozen-density coast formula gave 71.1 km and 19.7 km for a coast whose full-model answer is 34.0 km. Why must the truth lie between them, and which assumption makes the lower one too pessimistic?
