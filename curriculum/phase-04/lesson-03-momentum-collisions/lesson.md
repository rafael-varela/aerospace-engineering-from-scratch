# P04.03  Momentum & Collisions

> A rocket pushes on nothing. It throws mass backward, and momentum conservation does the rest, one pellet at a time.

| | |
|---|---|
| **Learning outcomes** | 1. State the impulse–momentum theorem and the conservation of linear momentum for a system with no net external force, and apply them to one-dimensional collisions of any coefficient of restitution, deriving the kinetic energy lost as ½μ(1 − e²)v_rel². 2. Derive the velocity gain of a body that throws pellets backward from momentum conservation alone, show that in the limit of many small pellets it becomes the Tsiolkovsky equation Δv = u ln(m₀/m_f), and derive the 1/n error of the discrete version. 3. Explain why F = d(mv)/dt is wrong for a rocket, write the correct variable-mass equation m dv/dt = ṁu + F_ext, and derive the closed-form burnout velocity and altitude of a vertical drag-free ascent. 4. Integrate the rocket with `scipy.integrate.solve_ivp`, verify the closed forms, verify that the momentum of rocket plus exhaust is conserved in free space, and measure gravity and drag losses. 5. Compute the velocity change of a staged vehicle, find the best propellant split, and explain why a single stage cannot exceed u ln(1/ε). |
| **Prerequisites** | `P04.01` (Newton's second law, drag), `P04.02` (work, kinetic energy, conservative and non-conservative forces), `P01.07` (vectors), `P02.04` (first-order systems), `P03.05` (`solve_ivp`). Skill check: you can sum a geometric-looking series by comparing it with an integral, and you can write a three-variable system dz/dt, dv/dt, dm/dt. |
| **Engineering question** | A 10 t single-stage vehicle carries 7.5 t of propellant with an exhaust speed of 2942 m/s (specific impulse 300 s). What velocity can it reach in free space, how much of that is lost to gravity and drag in a vertical ascent, what thrust-to-weight ratio should it have, and how much more could the same propellant give if it were split into two stages? |
| **Assumptions and limits** | One-dimensional motion; constant exhaust speed relative to the rocket and constant mass flow; uniform gravity g₀ (no decrease with altitude: at 127 km g is 4 % smaller); vertical flight (no gravity turn); drag from a constant C_D with an exponential atmosphere of scale height 8.5 km (a crude model, adequate for the size of the drag loss, not for structural loads); no pressure thrust term (the nozzle is assumed ideally expanded); structure dropped instantly at staging. The rocket equation itself has no approximations beyond constant u. |
| **Analytical work** | Impulse–momentum; collisions with restitution e and their energy loss; the centre-of-mass frame; the pellet rocket and its sum; the Riemann-sum error u(m₀ − m_f)(1/m_f − 1/m₀)/(2n); the Tsiolkovsky equation; the variable-mass equation of motion and why d(mv)/dt fails; the burnout velocity u ln R − g t_b and altitude in closed form; the gravity loss g t_b; the staging equation and the single-stage limit u ln(1/ε). |
| **Computational deliverable** | `starter/main.py`: `impulse`, `centre_of_mass`, `collide_1d`, `collision_energy_loss`, `pellet_rocket_velocities`, `tsiolkovsky_dv`, `mass_ratio_for_dv`, `pellet_convergence`, `rocket_rhs`, `simulate_rocket`, `burnout_closed_form`, `velocity_losses`, `stack_dv`, `structure_mass`, `best_two_stage_split` (`RocketFlight` and `main` are given). Tests prove: momentum conserved in every collision to 10⁻¹² and kinetic energy lost exactly as derived, never gained; equal masses swap velocities and a ball bounces off a wall; the impulse of a half-sine kick; one pellet is a single momentum balance and each later throw gains more; the pellet error is positive, of slope −1.00 ± 0.02 in n and equal to the derived constant to 0.1 %; the mass ratio inverts Tsiolkovsky and doubling u square-roots it; in free space v(m) = u ln(m₀/m) throughout the burn and the momentum of rocket plus exhaust stays zero; the vertical ascent equals the closed-form burnout velocity and altitude to 10⁻⁹ for three mass flows; the gravity loss is g t_b and falls with thrust; drag in a thinning atmosphere costs under 10 % and a constant-density atmosphere five times more; a rocket that cannot lift off is refused; the right-hand side balances momentum with gravity as the only external force; the stack Δv identities; two stages beat one and beat u ln(1/ε), at an optimum with equal stage mass ratios. |
| **Practical activity** | **computational experiment**: convergence of the pellet rocket to the rocket equation, a free-space momentum audit of rocket and exhaust, a thrust-to-weight sweep separating gravity and drag losses, and a staging sweep over the propellant split. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. No rocket, model rocket motor, compressed-gas or water rocket is built or fired in this lesson. |
| **Evidence** | Completed `main.py`, passing tests, `rocket_momentum.png`, the derivations of the collision formulas, the pellet sum and its limit, the burnout state and the staging limit, the four experiment tables, and a half-page note answering the engineering question. |
| **Assessment** | `assessments/phase-04/grading-rubric.md`; competencies 1–4. Pass: the rocket equation is derived from momentum conservation (not quoted); all public tests pass; the momentum audit is reported; the note separates ideal Δv, gravity loss and drag loss with numbers. |
| **References** | `@openstax-university-physics-vol1` Ch. 9 "Linear Momentum and Collisions": 9.1 Linear Momentum, 9.2 Impulse and Collisions, 9.3 Conservation of Linear Momentum, 9.4 Types of Collisions, 9.6 Center of Mass, 9.7 Rocket Propulsion (required-free: the impulse–momentum theorem, elastic and inelastic collisions, and 9.7's derivation of the rocket equation, which this lesson reaches by a different route). `@leishman-intro-aerospace-flight-vehicles` chapter "Rocket Engines", https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/rocket-engines/ (supplementary-free: thrust, specific impulse and exhaust velocity of real engines). |
| **Extensions** | ★ Add the pressure thrust (p_e − p_a)A_e with an exponential atmosphere and show how the effective exhaust speed rises with altitude. ★ Fly a gravity turn in two dimensions (thrust along the velocity) and compare the gravity loss with the vertical ascent. ★ Find the optimal propellant split for three stages and compare the gain from the second and the third stage. |
| **Time** | unmeasured |

## 1. The engineering problem

Every launch vehicle, satellite thruster and attitude-control jet obeys one equation, Tsiolkovsky's: the velocity a rocket gains is its exhaust speed times the logarithm of the ratio of its full mass to its empty mass. It is the most consequential equation in astronautics because of the logarithm. To reach orbit, about 9.4 km/s including losses, with an exhaust speed of 2942 m/s (a kerosene engine at a specific impulse of 300 s) the vehicle must be 24.4 times heavier full than empty: 96 % propellant, with the tanks, engines, guidance and payload in the other 4 %. No structure is that light, which is why every orbital vehicle stages, and why the exhaust speed is worth any engineering effort (at 4400 m/s, hydrogen, the ratio drops to 8.5).

This lesson derives that equation the way a physicist should: not from a formula for thrust, but from conservation of momentum, applied to a cart that throws bricks off its back. Throw one 7.5 t brick at 2942 m/s from a 10 t cart and the cart moves off at 2207 m/s. Throw it as two 3.75 t bricks and it reaches 2868 m/s; as ten, 3768; as a thousand, 4075; as infinitely many, **4078 m/s = 2942 ln 4**. The rocket is the limit of the cart. Then the lesson adds what a real ascent costs: a 150 s vertical burn at a thrust-to-weight ratio of 1.5 loses 1471 m/s, 36 % of the ideal, to gravity alone, and the atmosphere takes another 59 m/s. And it shows that splitting the same propellant into two stages raises the ideal velocity of a 500 kg payload from 6243 to 8559 m/s.

On the way there are collisions, because docking is one: a 12 t spacecraft closing at 0.10 m/s on a 420 t station has 60 J of kinetic energy, and the docking mechanism must absorb 58.3 J of it, since after a perfectly inelastic collision the pair moves at only 2.78 mm/s.

Success criteria, written first: momentum conserved to 10⁻¹² in every collision; the pellet rocket converging to Tsiolkovsky at order 1/n with the derived error constant; the simulated burnout state equal to the closed form to 10⁻⁹; the momentum of rocket plus exhaust conserved in free space; and the staging gain and its limit computed and explained.

## 2. Understand and derive

### 2.1 Momentum and impulse

The **linear momentum** of a particle is p = mv (OpenStax 9.1). Newton's second law in its original form says the net external force equals the rate of change of momentum, F = dp/dt, and integrating over a time interval gives the **impulse–momentum theorem** (9.2):

**J ≡ ∫F dt = Δp.**   (1)

A short, hard kick and a long, soft push with the same impulse give the same change of momentum; what differs is the peak force. A half-sine force of peak 2000 N lasting 10 ms has impulse 2000 × 0.01 × 2/π = 12.7 N·s and sends a 0.45 kg ball from rest to 28.3 m/s. For a *system* of particles the internal forces come in third-law pairs and cancel in the sum, so the total momentum changes only through external forces, and with none **the total momentum is conserved** (9.3). The **centre of mass** x_cm = Σmᵢxᵢ/Σmᵢ (9.6) moves at the constant velocity P/M when no external force acts, whatever the parts do.

### 2.2 Collisions in one dimension

Two bodies m₁, m₂ with velocities v₁, v₂ collide head-on. Momentum is conserved: m₁v₁ + m₂v₂ = m₁v₁′ + m₂v₂′. One more equation is needed, and it is a property of the materials, not of mechanics: the **coefficient of restitution** e, the ratio of separation speed to approach speed, v₂′ − v₁′ = −e(v₂ − v₁), with e = 1 for an **elastic** collision (no kinetic energy lost), e = 0 for a **perfectly inelastic** one (the bodies stick), and anything between for a real one (9.4). Solving the pair,

**v₁′ = (P + m₂e(v₂ − v₁))/M,  v₂′ = (P + m₁e(v₁ − v₂))/M,**   (2)

with P = m₁v₁ + m₂v₂ and M = m₁ + m₂. The kinetic energy lost is easiest in the centre-of-mass frame, where the total momentum is zero and the two bodies approach with the relative speed v_rel = v₁ − v₂; the kinetic energy there is ½μv_rel² with the **reduced mass** μ = m₁m₂/M, and after the collision the relative speed is e v_rel, so

**ΔK = ½ μ (1 − e²) v_rel².**   (3)

The kinetic energy of the centre-of-mass motion, ½MV_cm², can never be lost (momentum fixes it); only the relative-motion part can. Checks: equal masses, elastic, the second at rest: the velocities swap. A ball on a wall (m₂ → ∞), elastic: it bounces back at the same speed.

**Docking.** A 12 t vehicle at 0.10 m/s into a 420 t station at rest, e = 0: v′ = 12 000 × 0.10/432 000 = 2.78 mm/s; μ = 11.67 t; ΔK = ½ × 11 667 × 0.01 = **58.3 J** of the 60 J. Nearly all the kinetic energy goes into the docking mechanism's dampers, because the station is so much heavier. If the mechanism springs back (e > 0) the vehicle rebounds, which is why docking systems are built to be as inelastic as possible.

### 2.3 The pellet rocket

A cart of mass m moving at v throws a pellet of mass dm backward at speed u *relative to the cart after the throw*. In the frame where the cart was moving at v, momentum before is mv; after, the cart (m − dm) moves at v + dv and the pellet at v + dv − u. Conservation:

m v = (m − dm)(v + dv) + dm (v + dv − u)  ⇒  **dv = u dm/m**,   (4)

exactly, with m the mass *before* the throw. Nothing about forces, thrust or pressure was used: only the momentum of the system before and after. Throwing n equal pellets of mass Δm = (m₀ − m_f)/n from rest:

**v_n = u Σᵢ Δm/mᵢ,  mᵢ = m₀ − iΔm, i = 0 … n − 1.**   (5)

With one pellet, v₁ = u(m₀ − m_f)/m₀; for the 10 t cart with 7.5 t of propellant at u = 2942 m/s, 2207 m/s. With more pellets each throw gains more than the last, because the cart is lighter when it throws (the tests check that the increments rise). The sum is a left Riemann sum of ∫dm/m, and as n → ∞,

**Δv = u ∫_{m_f}^{m₀} dm/m = u ln(m₀/m_f)**   (Tsiolkovsky).   (6)

Because 1/m decreases across each interval, the left sum always underestimates the integral, and the error of a Riemann sum of a smooth function is half the step times the difference of the end values (Problem 2):

**u ln(m₀/m_f) − v_n ≈ u (m₀ − m_f)(1/m_f − 1/m₀)/(2n),**   (7)

3310/n m/s here. The code confirms it: n × error = 3104 at n = 10, 3289 at 100, 3309.5 at 10 000. Many small throws beat a few large ones because each pellet is thrown from a lighter vehicle: that is the whole content of the logarithm.

### 2.4 The rocket equation of motion, and why F = d(mv)/dt is wrong

In the continuous limit with mass flow ṁ = −dm/dt > 0 and external force F_ext, applying (4) over a time dt gives the **variable-mass equation**

**m dv/dt = ṁ u + F_ext.**   (8)

The term ṁu is the **thrust**: the rate at which momentum is thrown backward relative to the rocket. For ṁ = 50 kg/s and u = 2942 m/s it is 147.1 kN. The exhaust speed is usually quoted as a **specific impulse** I_sp = u/g₀ in seconds (300 s here); it is the impulse per unit *weight* of propellant (Leishman, "Rocket Engines", gives real values).

A tempting shortcut is F_ext = d(mv)/dt = m dv/dt + v dm/dt. It is wrong: it says that a rocket coasting at constant speed (ṁ = 0) is fine, but that a rocket leaking propellant gently from a hole with zero relative speed (u = 0, which cannot change its velocity) would accelerate or decelerate depending on v, i.e. on the frame of the observer. Newton's second law for a *fixed* set of particles is F = dp/dt; the rocket is not a fixed set, because particles leave it. The correct bookkeeping includes the momentum carried away, and the code tracks it explicitly as a fourth state, the momentum of the exhaust, dP_ex/dt = ṁ(v − u). With F_ext = 0 the total m v + P_ex must stay constant; the tests check it.

### 2.5 Vertical ascent: the closed form and the gravity loss

With F_ext = −mg (vertical, no drag), constant ṁ and m(t) = m₀ − ṁt, (8) integrates directly: dv/dt = ṁu/(m₀ − ṁt) − g, so

**v(t) = u ln(m₀/m(t)) − g t,  v_b = u ln(m₀/m_f) − g t_b,  t_b = (m₀ − m_f)/ṁ**,   (9)

and integrating once more (with ∫ln(m₀/m)dt evaluated by the substitution m = m₀ − ṁt),

**z_b = u t_b − u (m_f/ṁ) ln(m₀/m_f) − ½ g t_b².**   (10)

For ṁ = 50 kg/s: t_b = 150 s, v_b = 4078.5 − 1471.0 = **2607.5 m/s**, z_b = 441 299 − 203 924 − 110 325 m = **127.05 km**. The term g t_b is the **gravity loss**: every second spent burning against gravity costs 9.81 m/s, whatever the engine. It falls as the thrust rises (shorter burn): 1839 m/s at a thrust-to-weight ratio of 1.2, 441 m/s at 5. It would vanish only with an infinite thrust, which drag forbids: the faster the rocket climbs through the dense lower atmosphere, the larger the drag loss and the **dynamic pressure** q = ½ρv² that the structure must carry. The thrust-to-weight ratio of a real first stage, about 1.2–1.5 at lift-off, is that compromise (Problem 4).

### 2.6 Staging

A single stage with structural coefficient ε = m_s/(m_s + m_p) and no payload can reach at most u ln((m_s + m_p)/m_s) = **u ln(1/ε)**, 6774 m/s for ε = 0.1 and u = 2942 m/s, below orbital velocity. The way out is to discard the empty tanks of a first stage before lighting a second: each stage then burns with only the stages above it as dead weight, and the velocities add,

**Δv = Σᵢ uᵢ ln(m₀ᵢ/m_fᵢ),**   (11)

with m₀ᵢ the mass of stage i and everything above it at ignition and m_fᵢ the same at burnout. With 500 kg of payload, 20 t of propellant and ε = 0.1: one stage gives 6243 m/s; two equal stages 7517 m/s; the best split (87 % of the propellant in the first stage) 8559 m/s, above the single-stage limit. For equal u and ε the optimum has equal mass ratios in both stages (a classical result; the test checks it numerically), which puts most of the propellant in the first stage because the second stage plus payload is light.

### 2.7 Where the model fails

- **Constant u**: a real engine's exhaust speed rises with altitude as the outside pressure falls (sea-level and vacuum I_sp differ by roughly 10–25 % depending on the engine; ★ extension).
- **Uniform g**: g falls 4 % by 127 km; the gravity loss is overestimated by a few per cent.
- **Vertical flight**: a real vehicle turns toward horizontal, so gravity acts less and less along the velocity and the gravity loss is far smaller than g t_b; orbit needs horizontal speed, which a vertical ascent never gives.
- **Drag model**: the exponential atmosphere with constant C_D gives the size of the drag loss and of max-q; C_D of a rocket changes strongly through Mach 1.
- **Instant staging**, no residual propellant, no engine start transients.

## 3. Build the computational model

Open `starter/main.py`. `RocketFlight` and `main` are given.

### Step 1: collisions and impulse

`collide_1d(m1, v1, m2, v2, e)` is (2), with a `ValueError` for e outside [0, 1]; `collision_energy_loss` is (3); `impulse` is `np.trapezoid(F, t)`; `centre_of_mass` uses `np.tensordot(m, x, axes=(0, 0)) / m.sum()` so that it works for 1-D and 2-D positions.

### Step 2: pellets and Tsiolkovsky

```python
def pellet_rocket_velocities(m0, mf, u, n):
    dm = (m0 - mf) / n
    m = m0 - dm * np.arange(n)                 # masses before each throw
    return np.concatenate([[0.0], np.cumsum(u * dm / m)])
```

`tsiolkovsky_dv` and `mass_ratio_for_dv` are (6) and its inverse; `pellet_convergence` returns u ln(m₀/m_f) − v_n for each n.

### Step 3: the continuous rocket

```python
def rocket_rhs(t, y, mdot, u, g, k_drag, scale_height=8500.0):
    z, v, m, _ = y
    drag = k_drag * math.exp(-max(z, 0.0) / scale_height) * v * abs(v)
    return np.array([v, (mdot * u - m * g - drag) / m, -mdot, mdot * (v - u)])
```

The fourth component is the exhaust's momentum: each kilogram leaves at v − u in the ground frame. `simulate_rocket` refuses a vehicle whose thrust does not exceed its weight on the pad, integrates to t_b with DOP853 at rtol 10⁻¹¹ and returns a `RocketFlight`; `burnout_closed_form` is (9)–(10); `velocity_losses` splits the ideal Δv into what was achieved, the gravity loss g t_b, and the remainder, which is the drag loss.

### Step 4: staging

`stack_dv(payload, stages)` applies (11) with each stage's dead weight being everything above it; `structure_mass(m_p, ε) = εm_p/(1 − ε)`; `best_two_stage_split` maximises the two-stage Δv over the first stage's propellant fraction with a bounded scalar search.

Running `main()` prints the docking and collision examples, the pellet table of 2.3, the vertical ascent against (9)–(10) (2607.474 m/s and 127.0509 km, both to all printed digits), the free-space momentum audit, the drag case, the thrust-to-weight sweep and the staging results, and writes `rocket_momentum.png` (the pellet error against n on log axes with the predicted 1/n line; the velocity history in free space, with gravity, and with gravity and drag; Δv against the propellant split).

### Verify

Run `python -m pytest curriculum/phase-04/lesson-03-momentum-collisions`. The tests:

1. **Collisions** (three mass pairs × four values of e): momentum conserved to 10⁻¹²; restitution law; kinetic energy lost equal to (3) and never negative; e = 0 bodies move together.
2. **Limits**: equal-mass elastic swap; a ball bouncing off a 10¹² kg wall; in the centre-of-mass frame e = 0 loses all the kinetic energy; e = 1.5 raises.
3. **Impulse and centre of mass**: the half-sine kick's 12.73 N·s and the 28.3 m/s it gives a 0.45 kg ball; centre of mass in 1-D and 2-D.
4. **Pellets**: one pellet is u(m₀ − m_f)/m₀; increments rise; the error is positive with slope −1.00 ± 0.02 and n × error equal to the constant of (7) to 0.1 %.
5. **Mass ratio**: inverts (6); Δv = 0 for no propellant; doubling u gives the square root of the mass ratio.
6. **Free space**: burnout equals (6) to 10⁻⁹; v(m) = u ln(m₀/m) throughout; momentum of rocket plus exhaust stays zero to 10⁻⁶ of m₀u; the exhaust's momentum is negative; m(t) linear.
7. **Vertical ascent** at T/W₀ = 1.2, 1.5, 3.6: burnout velocity and altitude equal (9)–(10) to 10⁻⁹; gravity loss g t_b; drag loss zero; the g → 0 limit.
8. **Losses**: burnout velocity rises with thrust and stays below the ideal; drag in an 8.5 km atmosphere costs between 0 and 10 % and exactly the difference from the drag-free burnout; a constant-density atmosphere costs more than five times as much.
9. **Lift-off**: T/W₀ = 0.9 raises under gravity but flies in free space.
10. **Right-hand side**: a = ṁu/m; d(mv)/dt + dP_ex/dt = 0 without gravity and −mg with it.
11. **Staging**: single-stage and two-stage Δv against hand values; ε = 0.1 structure; two stages beat one and beat u ln(1/ε); the optimum is a maximum with equal stage mass ratios to 0.1 %.
12. `main` writes the figure.

### Compare

The independent checks are the closed forms (6), (9), (10) and the conservation law. `solve_ivp` reproduces (9) and (10) to 10⁻⁹ relative and holds the total momentum of rocket and exhaust at 3 × 10⁻⁴ kg·m/s against a rocket momentum of 1.0 × 10⁷ kg·m/s: eleven orders of magnitude. The audit is worth its line of code: a wrong sign in the thrust, or the d(mv)/dt mistake of 2.4, conserves nothing, and the audit says so before any plot does. The pellet sum is a second, entirely algebraic route to (6), and its convergence at exactly the predicted rate is the evidence that the continuum limit has been taken correctly.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction down before running it.

1. **Pellets.** For m₀ = 10 t, m_f = 2.5 t, u = 2942 m/s, compute v_n for n = 1, 2, 4, 10, 100, 1000, 10 000 and plot the error on log axes. Predict v₁ by hand, the slope, and the constant of (7).
2. **Momentum audit.** In free space, track m v, P_ex and their sum through the burn. Predict the final rocket momentum and the exhaust momentum.
3. **Thrust-to-weight sweep.** For T/W₀ = 1.2, 1.5, 2, 3, 5: burn time, gravity loss and burnout speed without drag; then with drag (k = 0.144 kg/m, H = 8.5 km) and the maximum dynamic pressure. Predict the gravity loss from g t_b before running, and the direction in which drag and max-q move.
4. **Staging.** Two stages, 500 kg payload, 20 t of propellant, ε = 0.1: Δv against the first stage's share of the propellant from 5 % to 95 %. Predict that one stage cannot exceed 6774 m/s and that the optimum puts most of the propellant in the first stage.

Raw data format: CSV tables `n, v_n_mps, error_mps, n_times_error` (1); `t_s, m_kg, rocket_momentum, exhaust_momentum, total` (2); `TW0, burn_s, gravity_loss_mps, drag_loss_mps, burnout_mps, max_q_kPa` (3); `stage1_fraction, dv_mps` (4); plus `rocket_momentum.png`.

What you should find: (1) 2206.50, 2868.45, 3374.00, 3767.99, 4045.58, 4075.16, 4078.14 m/s against 4078.47 m/s; n × error approaching 3309.7; (2) the rocket ends with +1.0196 × 10⁷ kg·m/s and the exhaust with exactly −1.0196 × 10⁷ kg·m/s, total 3 × 10⁻⁴ kg·m/s; (3) gravity losses 1839, 1471, 1103, 736 and 441 m/s; with drag at T/W₀ = 1.5 the burnout speed is 2548 m/s (drag loss 59 m/s) and max-q 29.9 kPa at 58 s and 10.5 km; (4) 6243 m/s for one stage, 7517 m/s for equal stages, 8559 m/s at an 87.1 %/12.9 % split.

## 5. Validate and improve

| quantity | hand (§2) | code | difference | cause |
|---|---|---|---|---|
| docking velocity, energy absorbed | 2.78 mm/s, 58.3 J | 2.778 mm/s, 58.33 J | 0 | — |
| one-pellet velocity (m/s) | 2206.5 | 2206.50 | 0 | — |
| n × pellet error, n = 10 000 | 3309.7 | 3309.5 | 0.006 % | next term of the Riemann-sum error |
| Tsiolkovsky, R = 4 (m/s) | 4078.5 | 4078.471 (free-space ODE) | 10⁻⁹ | solver tolerance |
| burnout v at T/W₀ = 1.5 (m/s) | 2607.47 | 2607.474 | 10⁻⁹ | — |
| burnout altitude (km) | 127.051 | 127.0509 | 10⁻⁹ | — |
| total momentum, free space (kg·m/s) | 0 | 3.4 × 10⁻⁴ | 10⁻¹¹ of the rocket's | — |
| gravity loss, T/W₀ = 1.5 (m/s) | 1471.0 | 1471.0 | 0 | definition |
| drag loss (m/s) | small | 59.4 | — | no closed form |
| best two-stage Δv (m/s) | > 6774 | 8559 | — | — |

Everything with a closed form matches to the solver tolerance. The two numbers without one, the drag loss and the staging optimum, are checked by bounds (the drag loss is positive and below 10 %, the two-stage Δv exceeds the single-stage limit) and by an independent property (equal mass ratios at the optimum). The model's weakness is not numerical: a vertical ascent wastes most of its gravity loss, and the right improvement is a two-dimensional gravity-turn trajectory (★ extension), where the gravity loss of a real vehicle is typically 1–1.5 km/s on the way to orbit rather than g t_b. The second improvement is the altitude-dependent exhaust speed.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (2)–(4), (6), (7), (9), (10) and the single-stage limit; `rocket_momentum.png`; the four experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) answering the engineering question: the free-space Δv, the gravity and drag losses at your chosen thrust-to-weight ratio with the reason for the choice, and the staging gain.

### Adversarial self-checks

- **What if the propellant mass doubles?** The mass ratio does not double: with 15 t of propellant on the same 2.5 t dry mass, R = 7 and Δv = 2942 ln 7 = 5725 m/s, 40 % more for twice the propellant (and in practice the tanks grow too, so ε stays put and the gain is smaller still). The logarithm is the reason rockets stage instead of growing.
- **What if the exhaust speed goes to zero?** No thrust, no Δv, and (4) gives dv = 0 for each pellet: a leak changes nothing. This is the case that exposes the d(mv)/dt mistake, which would predict a velocity change from a leak.
- **Which assumption fails first in the real world?** Vertical flight: no launch vehicle flies straight up for 150 s, and the vertical gravity loss (1471 m/s) overstates a gravity turn's by a large factor. Next, constant u: sea-level and vacuum exhaust speeds differ by roughly 10–25 %.
- **What measurement could falsify the model?** A static-fire thrust measurement together with the measured mass-flow rate gives u = T/ṁ directly; a flight's measured burnout velocity below u ln R − g t_b by much more than the predicted drag loss would point at the thrust or the drag model.
- **Does the result survive a change of solver?** DOP853 matches the closed forms to 10⁻⁹; the pellet sum is a solver-free route to the same Δv; three routes, one number.
- **Which decision changes if the uncertainty doubles?** If the structural coefficient were uncertain between 0.08 and 0.12, the single-stage limit would span u ln(1/ε) = 6238 to 7431 m/s: a single stage still cannot reach orbit at any point of that range, so the decision to stage does not change; the best split moves by a few per cent.

## Check yourself

1. State the impulse–momentum theorem. A 2000 N peak half-sine force lasts 10 ms: what impulse does it deliver?
2. In a perfectly inelastic collision, which quantity is conserved and which is not? Where does the lost kinetic energy go?
3. A 12 t vehicle docks at 0.1 m/s with a 420 t station. Why is almost all of its kinetic energy absorbed by the docking mechanism?
4. Derive dv = u dm/m from momentum conservation for a single throw.
5. Why do many small pellets give more velocity than one large one with the same total mass?
6. What mass ratio does a rocket with u = 3000 m/s need for Δv = 3000 m/s? What Δv does a mass ratio of 10 give?
7. Why is F = d(mv)/dt wrong for a rocket?
8. What is the gravity loss of a 100 s vertical burn, and why can a higher thrust not reduce it to zero in the atmosphere?
