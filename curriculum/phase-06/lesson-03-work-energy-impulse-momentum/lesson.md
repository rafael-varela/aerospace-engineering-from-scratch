# P06.03  Work-Energy & Impulse-Momentum

> Work removes the energy over a distance; impulse removes the momentum over a time. The arresting gear has to do both, and the peak force is decided by how evenly it spreads the work.

| | |
|---|---|
| **Learning outcomes** | 1. Apply the work–energy theorem along a stopping stroke and the impulse–momentum theorem over its duration, and say what each one gives (a distance, a time, two different averages of the force). 2. Solve the work–energy equation in closed form for a force that depends on position *and* speed, F = F₀ + kx + cv², by writing it as a linear ODE for v² in x, and prove that the force is then monotonic along the stroke, so the peak is at engagement or at the stop. 3. Implement and test the stopping distance, peak force and stopping time from those two theorems, and check all three, plus ∫F dt = m v₀ and the energy balance, against `solve_ivp`. 4. Design the constant-force gear (efficiency 1), and show by a sweep how the mix of spring, damper and brake terms sets the peak deceleration and how much energy is stored rather than dissipated. 5. Decide how a gear set for one aircraft treats a heavier or faster one. |
| **Prerequisites** | `P04.02` (work, kinetic energy, the work–energy theorem), `P04.03` (momentum and impulse), `P06.02` (Newton's law along a path), `P02.01` (first-order linear ODEs and the integrating factor), `P03.02` (root finding), `P03.04` (numerical quadrature, integrable singularities), `P03.06` (`solve_ivp` events). Skill check: you can solve u′ + αu = βx + γ by an integrating factor. |
| **Engineering question** | A 20 t carrier jet (illustrative) catches the wire at 65 m/s. With a lumped cable-and-absorber force F = F₀ + kx + cv², how far does it run out, what is the peak deceleration, how long does it take, and how should the three terms be mixed to stop it within 100 m with the smallest peak? The answers: any of the designs below stops it in exactly 100 m, but the peak ranges from 2.154 g (constant force: c = 100 N·s²/m², k = 4225 N/m, or a pure brake) to 4.31 g (a pure spring), in times from 2.42 to 3.96 s; work–energy, impulse–momentum and `solve_ivp` agree to 10⁻⁹ m and 10⁻¹⁰ s. |
| **Assumptions and limits** | The aircraft is a particle moving in a straight line on the deck; its own thrust, wheel drag and the wind over the deck are ignored. The force law is an illustrative lump of the cable, its sheaves and its energy absorber, with constant coefficients: a real gear's force comes from cable geometry and a hydraulic engine whose valve is scheduled for each aircraft. The engagement is instantaneous (a real cable has wave dynamics and a rising force). The k x term is treated as storing energy, which a real gear must not hand back. No number here is data for any real arresting gear. |
| **Analytical work** | The work–energy theorem as d(½mv²)/dx = −F; for F = F₀ + kx + cv², a linear ODE in u = v² solved by an integrating factor; the stopping distance as its root; the force's own ODE F′ = k − (2c/m)F and the monotonic force; the constant-force design k = 2cF/m; the stopping time ∫dx/v with its integrable endpoint singularity removed by substitution; impulse = m v₀; the energy split between brake, damper and the stored k x term; stroke efficiency. |
| **Computational deliverable** | `starter/main.py`: `arrestor_force`, `v_squared`, `force_along_stroke`, `stopping_distance`, `peak_force`, `constant_force_design`, `stopping_time`, `momentum_average_force`, `energy_average_force`, `simulate_arrestment`, `energy_audit`, `efficiency` (`_phi1`, `_phi2` and `main` are given). Tests prove: the closed form solves its ODE for five force laws; the brake, spring and pure-damper limits (the last never stops); continuity as c → 0; more force stops sooner and more mass runs further; the force is monotonic and the peak at an end, never below either average; the constant-force design is constant with efficiency 1 and time 2x/v₀; a pure spring stops in a quarter period with efficiency ½; `solve_ivp` reproduces distance, time and peak to 10⁻⁷, its impulse equals m v₀ and its energy audit closes to 10⁻⁵. |
| **Practical activity** | **computational experiment**: five gear designs with the same 100 m runout compared on peak deceleration, time and energy split; the three routes (work–energy, impulse–momentum, time integration) compared; and an off-design sweep of aircraft mass and engagement speed. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `arrestment.png`, the derivations of the v²(x) closed form, the monotonic-force result and the constant-force design, the design and off-design tables, and a half-page note recommending a force mix. |
| **Assessment** | `assessments/phase-06/grading-rubric.md`; competencies 1–4. Pass: both theorems are stated and used for what each gives; the closed form is derived; all public tests pass; the three routes are compared with numbers; the note states the peak deceleration of the recommended design and how it behaves off design. |
| **References** | `@mit-ocw-16-07-dynamics` L9 "Linear impulse and momentum; collisions"; L12 "Work and energy"; L13 "Conservative internal forces and potential energy" (required-free). `@openstax-university-physics-vol1` Ch. 7: 7.1 Work; 7.3 Work-Energy Theorem; Ch. 8: 8.2 Conservative and Non-Conservative Forces; 8.3 Conservation of Energy; Ch. 9: 9.2 Impulse and Collisions (required-free). `@scipy-userguide` "Ordinary differential equations (solve_ivp)" (required-free: the terminal event that stops the run at v = 0). |
| **Extensions** | ★ Model the cable geometry: the hook pulls the middle of a cable spanning the deck, so the absorber's force T acts on the aircraft as 2T sin β with tan β = x/(span/2); find the absorber law T(payout, payout rate) that gives a constant force on the aircraft. ★ Add the aircraft's own thrust (pilots go to full power on touchdown in case the hook misses) and redo the energy audit. ★ Give the engagement a finite force rise time (a first-order lag on F) and find how much it raises the peak. |
| **Time** | unmeasured |

## 1. The engineering problem

A carrier deck gives a landing jet about a hundred metres to stop. The tailhook catches one of several steel cables stretched across the deck; the cable runs over sheaves to an energy absorber below deck, which pays out cable against a retarding force. In a couple of seconds the gear must take away the kinetic energy of an aircraft of twenty tonnes or more at well over 100 knots. Too little force and the aircraft runs off the end of the deck; too much and the hook, the airframe or the pilot is overloaded. The quantity that matters to all three is the peak deceleration, and for a fixed runout it depends entirely on *how* the force varies along the stroke.

The two theorems of this lesson answer complementary questions. The **work–energy** theorem relates the force to *distance*: the work ∫F dx equals the kinetic energy removed. The **impulse–momentum** theorem relates it to *time*: the impulse ∫F dt equals the momentum removed. Neither alone gives the peak force when the force depends on the speed, as a hydraulic damper's does. The trick of this lesson is that for a force of the form F₀ + kx + cv², the work–energy equation becomes a *linear* differential equation for v² as a function of x, which has a closed form. That gives the whole force profile along the stroke without ever integrating in time. Time integration with `solve_ivp` is then the independent check.

Success criteria, written first: the closed-form v²(x) solves its ODE; the stopping distance, the peak force and the stopping time from the two theorems agree with `solve_ivp` to its tolerance; the simulated impulse equals m v₀ and the simulated energy split adds up to the initial kinetic energy; the constant-force design has efficiency 1.

## 2. Understand and derive

### 2.1 The two theorems for a stopping stroke

The aircraft (mass m) engages at x = 0 with speed v₀ and is slowed by a force F ≥ 0 acting backward. Newton's law, m dv/dt = −F, can be read two ways.

**Work–energy** (OpenStax 7.3, Work-Energy Theorem; 16.07 L12). Multiply by v = dx/dt: m v dv/dt = d(½mv²)/dt = −F dx/dt, so along the path

d(½ m v²)/dx = −F,  and over the whole stroke  ∫₀^{x_s} F dx = ½ m v₀².   (1)

**Impulse–momentum** (OpenStax 9.2, Impulse and Collisions; 16.07 L9). Integrate in time:

∫₀^{t_s} F dt = m v₀.   (2)

From (1) the **distance-average force** is F̄_x = ½ m v₀²/x_s; from (2) the **time-average force** is F̄_t = m v₀/t_s. They are different averages of the same force, and the peak is at least as large as both. For the illustrative jet (m = 20 000 kg, v₀ = 65 m/s) the kinetic energy is 42.25 MJ and the momentum 1.300 MN·s. A **constant** force over x_s = 100 m is 422.5 kN, a deceleration of 21.125 m/s² = 2.154 g lasting 2x_s/v₀ = 3.077 s. For a fixed runout no force history can have a smaller peak than that, because the peak is never below the distance average.

### 2.2 A force that depends on position and speed

Lump the gear into

**F(x, v) = F₀ + k x + c v²,**   (3)

a constant drag F₀ (friction, a brake), a term that grows with payout (stretch and rising absorber pressure, written k x), and a hydraulic damper's c v² (fluid forced through an orifice gives a force close to proportional to the square of the flow rate). With u = v², (1) becomes

**u′ + (2c/m) u = −(2/m)(F₀ + k x),  u(0) = v₀²:**   (4)

a linear first-order ODE in x (`P02.01`), because v appears only as v². With α = 2c/m the integrating factor e^{αx} gives, for c > 0,

**u(x) = (v₀² − a) e^{−αx} + a + b x,  b = −k/c,  a = −F₀/c + k m/(2c²),**   (5)

and for c = 0 simply u = v₀² − (2/m)(F₀ x + ½ k x²). Check by substitution: u′ + αu = −α(v₀² − a)e^{−αx} + b + α(v₀² − a)e^{−αx} + αa + αbx = (b + αa) + αbx, and b + αa = −2F₀/m, αb = −2k/m ✓. The **stopping distance** x_s is the root of u(x_s) = 0, found with `brentq` (`P03.02`).

*A numerical warning.* For small c, a = km/(2c²) is huge and (5) subtracts two huge numbers. The code uses the equivalent form u = v₀² e^{−y} − (2F₀/m) x φ₁(y) − (2k/m) x² φ₂(y) with y = αx, φ₁ = (1 − e^{−y})/y and φ₂ = (e^{−y} − 1 + y)/y² (series below y = 0.01), which is accurate for every c including zero. The first version of the code used (5) as written; a design search that tried c = 10⁻⁶ N·s²/m² returned a stopping distance 26 % short, and the bracketing failed. The test `test_small_damper_is_continuous_with_no_damper` now guards it.

**Limits.** Pure brake (k = c = 0): x_s = m v₀²/(2F₀). Pure spring: x_s = v₀ √(m/k), a quarter of an oscillation. Pure damper (F₀ = k = 0): u = v₀² e^{−αx}, so v = v₀ e^{−cx/m}: the speed falls exponentially with distance and **never reaches zero**. A damper alone cannot stop an aircraft in a finite distance; it needs a brake or a spring term to finish the job.

### 2.3 Where is the peak? A one-line proof

Differentiate the force along the stroke, F(x) = F₀ + kx + c u(x), and use (4):

**F′(x) = k + c u′ = k − (2c/m) F.**   (6)

This is another linear ODE, with solution F(x) = k m/(2c) + (F(0) − k m/(2c)) e^{−αx}. So F relaxes monotonically from its engagement value F(0) = F₀ + c v₀² toward km/(2c). It never has an interior maximum: **the peak force is at engagement or at the stop**, max(F₀ + c v₀², F₀ + k x_s). Damper-heavy gears peak at engagement, spring-heavy gears at the end.

**The constant-force design.** F is constant along the whole stroke if F′ = 0, i.e. if **k = 2cF/m**, and then F = F₀ + c v₀² must equal the ideal ½ m v₀²/x_s. With F₀ = 0 and x_s = 100 m: **c = 422 500/65² = 100.0 N·s²/m², k = 2 × 100 × 422 500/20 000 = 4225 N/m.** The damper's falling share is exactly made up by the spring's rising share. (Substituting, k = m v₀²/(2x_s²): the constant-force spring term is always half the spring that would stop the aircraft alone in the same distance, m v₀²/x_s².)

### 2.4 The stopping time

From impulse–momentum the time is not given directly, but t = ∫dt = ∫₀^{x_s} dx/v(x) with v = √u from (5). Near the stop u ≈ |u′(x_s)|(x_s − x), so the integrand behaves like 1/√(x_s − x): infinite, but integrable (`P03.04`). The substitution x = x_s(1 − s²) turns it into a smooth integrand,

**t_s = ∫₀¹ 2 x_s s / √(u(x_s(1 − s²))) ds,**   (7)

which `scipy.integrate.quad` handles to 10⁻¹¹. With t_s, (2) gives the time-average force m v₀/t_s.

### 2.5 Where the energy goes

Along the stroke, (1) splits into three parts: the brake's work F₀ x_s (dissipated), the damper's ∫c v² dx = ∫c v³ dt (dissipated), and the k x term's ½ k x_s² (**stored**, like a spring). At the stop the stored part would push the aircraft back; a real gear must not return it. For a 100 m runout the stored share is k/(2 × 4225 N/m) of the energy: half for the constant-force design, all of it for a pure spring. A pure brake also gives a constant force and dissipates everything, but cannot adapt to speed (part 4).

### 2.6 Where the model stops

- **Engagement is not instantaneous.** The force on the hook rises over a short time as stress waves run along the cable; (3) jumps to F₀ + c v₀² at once. Damper-heavy designs, which peak at engagement, are the most affected.
- **Constant coefficients.** A real absorber's valve opening is scheduled with payout and set for each aircraft's weight; (3) is a fixed law.
- **Cable geometry** (Extension 1) turns the absorber force into a hook force that grows with the cable angle.
- **The aircraft is a particle** with no thrust, no wheel drag and no wind over the deck (which reduces the engagement speed relative to the deck in reality).

## 3. Build the computational model

Open `starter/main.py`. `_phi1`, `_phi2` and `main` are given.

### Step 1: the force and the work–energy closed form

`arrestor_force(x, v, F0, k, c)` is (3). `v_squared(x, v0, m, F0, k, c)` is (5) in the stable form of 2.2:

```python
y = 2.0 * c / m * x
return v0 ** 2 * np.exp(-y) - 2.0 * F0 / m * x * _phi1(y) - 2.0 * k / m * x ** 2 * _phi2(y)
```

`force_along_stroke` is F₀ + kx + c u(x).

### Step 2: distance, peak and the ideal design

`stopping_distance` brackets the root of u by doubling and calls `brentq`; it returns inf for a pure damper and uses the quadratic formula when c = 0. `peak_force` returns max(F(0), F(x_s)) and where it acts, from (6). `constant_force_design(v0, m, x_stop, F0)` returns (c, k) from 2.3.

### Step 3: time and the averages

`stopping_time` is (7) with `quad`. `momentum_average_force` and `energy_average_force` are m v₀/t_s and ½ m v₀²/x_s; `efficiency` is their distance average over the peak.

### Step 4: the independent check

`simulate_arrestment` integrates m ẍ = −F(x, ẋ) with `solve_ivp` (rtol = atol = 10⁻¹⁰, dense output) until a terminal event at v = 0. `energy_audit` samples the dense solution and integrates F dt (the impulse) and c v³ dt (the damper's dissipation) with the trapezoid rule, and adds F₀ x_s and ½ k x_s².

Running `main()` prints the ideal force, the constant-force design, the five designs, the three-way check with the energy split, the two averages, and the off-design table, and writes `arrestment.png` (deceleration along the stroke from work–energy; deceleration against time from `solve_ivp`; v(x) from the closed form with simulated points on top).

### Verify

Run `python -m pytest curriculum/phase-06/lesson-03-work-energy-impulse-momentum`. The tests:

1. **Closed form**: (4) is satisfied by central differences of `v_squared` for five force laws; u(0) = v₀², u(x_s) = 0.
2. **Limits**: brake m v₀²/(2F₀); spring v₀√(m/k); pure damper never stops (inf distance and time) with v = v₀ e^{−cx/m}.
3. **Small damper**: c = 10⁻⁹, 10⁻⁶, 10⁻³ give the spring-only distance to 10⁻⁵.
4. **Monotonicity in the parameters**: more F₀, k or c stops sooner; double mass runs further.
5. **Peak**: F(x) monotonic for five laws; the peak equals the sampled maximum, sits at an end, and is not below either average.
6. **Constant-force design**: runout as designed, F constant to 10⁻⁹, efficiency 1, t_s = 2x_s/v₀, the two averages equal; an impossible F₀ is refused.
7. **Pure spring**: t_s = (π/2)√(m/k), efficiency ½.
8. **Time integration** for five laws: stop distance and time to 10⁻⁷, impulse = m v₀ and energy balance to 10⁻⁵, peak to 10⁻⁶.
9. **Speed profile**: simulated v² against (5) along the stroke to 10⁻⁷.
10. `main` writes the figure.

### Compare

| design | x_s sim − WE (m) | t_s sim − IM (s) | F_peak sim − WE (N) | impulse/(m v₀) − 1 | energy balance − 1 |
|---|---|---|---|---|---|
| constant force | 9.2 × 10⁻¹⁰ | 3.4 × 10⁻¹¹ | 5.3 × 10⁻⁶ | 1.7 × 10⁻¹¹ | 3.1 × 10⁻⁸ |
| damper-heavy | 1.0 × 10⁻⁹ | 8.2 × 10⁻¹¹ | 0 | 5.6 × 10⁻⁸ | 1.9 × 10⁻⁷ |
| spring-heavy | −1.1 × 10⁻⁹ | 2.5 × 10⁻¹¹ | −7.1 × 10⁻⁶ | −1.0 × 10⁻⁸ | 4.2 × 10⁻⁹ |
| spring only | −1.9 × 10⁻⁹ | −2.2 × 10⁻¹² | −1.6 × 10⁻⁵ | −1.3 × 10⁻⁸ | −3.8 × 10⁻¹¹ |
| brake only | −1.4 × 10⁻¹⁴ | −1.2 × 10⁻¹³ | 0 | 4.4 × 10⁻¹⁶ | −2.2 × 10⁻¹⁶ |

Work–energy (WE) and impulse–momentum (IM) agree with time integration to the solver's tolerance. The impulse and energy columns are limited by the trapezoid rule on 4001 samples of the dense output (largest for the damper-heavy design, whose force falls steeply at engagement), not by the physics.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** All numbers are illustrative simulations of a lumped force law; nothing here is data for a real arresting gear. Write each prediction before running.

1. **Five ways to stop in 100 m.** Hold the runout at 100 m and change the mix: constant force (c = 100, k = 4225); damper-heavy (k = 2112.5, c chosen by `brentq` to give 100 m); spring-heavy (k = 6337.5, c chosen likewise); spring only (k = 8450, c = 0); brake only (F₀ = 422.5 kN). Predict the ranking of peak deceleration and where each peak acts. *Finding:*

| design | k (N/m) | c (N·s²/m²) | t_s (s) | F_peak (kN) | at | peak g | efficiency | energy stored in kx |
|---|---|---|---|---|---|---|---|---|
| constant force | 4225.0 | 100.000 | 3.0769 | 422.5 | everywhere | 2.154 | 1.000 | 50 % |
| damper-heavy | 2112.5 | 193.755 | 3.9566 | 818.6 | engagement | 4.174 | 0.516 | 25 % |
| spring-heavy | 6337.5 | 42.417 | 2.6681 | 633.7 | stop | 3.231 | 0.667 | 75 % |
| spring only | 8450.0 | 0 | 2.4166 | 845.0 | stop | 4.308 | 0.500 | 100 % |
| brake only (F₀ = 422.5 kN) | 0 | 0 | 3.0769 | 422.5 | everywhere | 2.154 | 1.000 | 0 % |

The damper-heavy design is the slowest (3.96 s) yet has nearly twice the ideal peak: its distance-average force is 422.5 kN like every design here, but its time-average is only 328.6 kN, because it spends a long time crawling at low speed and low force near the end.
2. **Three routes.** For each design compare the work–energy distance and peak, the impulse–momentum time, and `solve_ivp`. *Finding:* the Compare table; agreement to 10⁻⁹ m and 10⁻¹⁰ s.
3. **Off design.** Keep the constant-force gear set for 20 t at 65 m/s and arrest 16 t and 24 t aircraft and 60 and 70 m/s engagements. Compare with the pure brake set for the same case. Predict which gear adapts better to speed. *Finding:*

| mass (t) | v₀ (m/s) | c v² + kx gear: x_s (m) | peak g | peak at | brake gear: x_s (m) | g |
|---|---|---|---|---|---|---|
| 16 | 65 | 86.76 | 2.693 | engagement | 80.00 | 2.693 |
| 20 | 65 | 100.00 | 2.154 | — | 100.00 | 2.154 |
| 24 | 65 | 112.14 | 2.013 | stop | 120.00 | 1.795 |
| 20 | 60 | 94.24 | 2.030 | stop | 85.21 | 2.154 |
| 20 | 70 | 105.56 | 2.498 | engagement | 115.98 | 2.154 |
| 24 | 70 | 118.50 | 2.127 | stop | 139.17 | 1.795 |

The speed-sensitive gear absorbs speed errors: 70 m/s overruns by 5.6 m instead of 16 m, because the c v² term grows with v₀². Neither absorbs weight errors well (a 20 % heavier aircraft overruns by 12–20 m): that is why the real gear is set for each aircraft's weight before it lands.

Raw data format: CSV tables `design, F0_N, k_Npm, c_Ns2pm2, x_stop_m, t_stop_s, F_peak_N, x_peak_m, eff, stored_frac` (1); `design, dx_m, dt_s, dF_N, impulse_err, energy_err` (2); `mass_kg, v0_mps, x_stop_m, peak_g, peak_at, brake_x_m, brake_g` (3); plus `arrestment.png`.

## 5. Validate and improve

| quantity | analytic (§2) | `solve_ivp` | difference | cause |
|---|---|---|---|---|
| stopping distance, constant force (m) | 100.000 (root of (5)) | 100.000 | 9 × 10⁻¹⁰ | solver tolerance |
| stopping time, constant force (s) | 3.0769 (7); 2x_s/v₀ = 3.0769 | 3.0769 | 3 × 10⁻¹¹ | solver tolerance |
| peak force, damper-heavy (kN) | 818.6 at engagement (6) | 818.6 | 0 | — |
| impulse, all designs | m v₀ = 1.300 MN·s (2) | within 6 × 10⁻⁸ relative | 10⁻⁸ | trapezoid on samples |
| energy, all designs | ½ m v₀² = 42.25 MJ (1) | within 2 × 10⁻⁷ relative | 10⁻⁷ | trapezoid on samples |

The engineering answer: **for a 100 m runout the illustrative 20 t jet at 65 m/s cannot be stopped with less than 2.154 g, and a constant-force gear (or a brake) achieves it in 3.08 s; a damper-heavy mix doubles the peak to 4.17 g at engagement, a spring alone to 4.31 g at the end.** The work–energy closed form gives the whole force profile along the stroke, the impulse–momentum integral gives the time, and time integration confirms both to 10⁻⁹.

What to improve: the engagement transient (Extension 3), which most raises the peak of damper-heavy designs; the cable geometry (Extension 1); and a rule for not returning the stored ½kx² (a valve that locks at the stop). What not to improve: the solver tolerance.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (4)–(7) and of the constant-force design; `arrestment.png`; the three experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) recommending a force mix for the illustrative jet with its peak deceleration, time, energy split and off-design behaviour.

### Adversarial self-checks

- **What if the engagement speed doubles?** Kinetic energy quadruples. With the constant-force gear, 130 m/s is far off design: the c v² term starts at four times the force, so the peak is at engagement and close to four times the design value. A brake gear would need four times the runout.
- **What if c → 0 or k → 0?** c → 0 is the pure spring (the closed form must stay continuous: test 3). k → 0 with F₀ = 0 never stops: the gear needs a term that does not vanish as v → 0.
- **Which assumption fails first in the real world?** Instantaneous engagement: the real force rises over a short time, and the damper term's initial spike is smoothed and partly converted into a cable-wave load on the hook.
- **What measurement could falsify the model?** A recorded runout and deceleration history from a real arrestment: if the force were of form (3), F′ + (2c/m)F would be constant (= k) along the stroke; any other shape falsifies the lump.
- **Does the result survive a change of route?** Distance, time and peak agree by three routes to 10⁻⁹; the peak's location is proved by (6), not just observed.
- **Which decision changes if the uncertainty doubles?** If the aircraft's mass is known only to ±20 % instead of ±10 %, the constant-force gear's runout spans 86.8–112.1 m instead of 93.5–106.2 m: a 105 m deck limit, just exceeded at +10 % (106.2 m), is then exceeded by 7 m, and the gear must be set for the heaviest credible weight.

## Check yourself

1. What does the work–energy theorem give you about a stopping stroke that the impulse–momentum theorem does not, and vice versa?
2. Why does a force that depends on v² still give a closed-form v(x)?
3. Why can a damper with F = cv² alone never stop an aircraft in a finite distance?
4. Where along the stroke can the peak force of F₀ + kx + cv² occur, and why only there?
5. Which force history gives the smallest peak for a given runout, and what is its efficiency?
6. A gear stops a jet in 100 m and 4 s. Which average force follows from each theorem?
7. Why must the k x part of the force not return its energy at the end of the stroke?
8. Why is the real gear set for each aircraft's weight but not for its exact speed?
