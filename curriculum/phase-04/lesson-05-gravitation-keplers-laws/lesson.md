# P04.05  Gravitation & Kepler's Laws

> One force law, F = GMm/r², and two conserved quantities, energy and angular momentum, give you every orbit: its size from the energy, its shape from the angular momentum, its period from the size alone. A simulation that does not conserve them is not simulating gravity.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the circular-orbit speed, the period (Kepler's third law), the escape speed and the specific orbital energy ε = −μ/2a from Newton's law of gravitation, and from energy and angular momentum conservation derive the vis-viva relation v² = μ(2/r − 1/a). 2. Recover the semi-major axis, eccentricity, periapsis and apoapsis of an orbit from one position and velocity, and decide from the sign of ε whether the object is bound. 3. Implement and test a leapfrog (kick-drift-kick) integrator for the two-body problem, show that it conserves angular momentum to rounding and keeps the energy error bounded, and measure its second-order convergence against a `scipy.integrate.solve_ivp` oracle. 4. Check Kepler's three laws on the simulated orbit (ellipse with the primary at a focus, constant areal velocity, T²/a³ = 4π²/μ) and reduce the full two-body problem to one body about the centre of mass. 5. Size a Hohmann transfer from low Earth orbit to geostationary orbit and explain its Δv in terms of orbital energy. |
| **Prerequisites** | `P03.07` (symplectic integrators: the leapfrog's bounded energy error), `P03.05` (RK4 and the two-body test problem, energy drift as a diagnostic), `P02.04` (second-order ODEs as first-order systems), `P01.07` (vectors, cross product), `P04.02` (work, potential energy), `P04.03` (momentum conservation, centre of mass), `P04.04` (angular momentum). Skill check: you can write ẍ = −μx/r³ as a four-component first-order system and compute r × v for a planar state. |
| **Engineering question** | For a satellite launched into a 400 km circular orbit, what are its speed and period, how much Δv does it take to reach geostationary orbit by a Hohmann transfer, and can a fixed-step simulation of its orbit be trusted for ten orbits? |
| **Assumptions and limits** | Two point masses (or spheres with spherically symmetric density); Newtonian gravity; no drag, no thrust, no third body. The Earth's oblateness (J₂) moves the orbit plane and the perigee by degrees per day in low orbit, drag decays orbits below about 500 km in months to years, and the Moon and Sun matter at geostationary altitude; general relativity is a 10⁻⁹ effect here (it matters for the clocks of `P04.13`, not for the orbit). Impulsive burns in the Hohmann transfer are an idealisation of finite burns. |
| **Analytical work** | Inverse-square law and g(h); circular speed and period; potential energy −GMm/r and escape speed; ε and h from a state; vis-viva from conservation at periapsis and apoapsis; Kepler's second law as conservation of angular momentum; the third law for ellipses; the centre-of-mass reduction with μ = G(M + m); the Hohmann Δv. |
| **Computational deliverable** | `starter/main.py`: the force law and its consequences, `orbital_elements`, a hand-written `leapfrog`, `propagate_leapfrog`, `propagate_solve_ivp` (oracle), `areal_velocity`, `period_from_trajectory`, `simulate_two_body`, `hohmann_transfer`. Tests (17) check the inverse-square law and Newton's third law, g(h) limits, v² / r = g, vis-viva limits (circular, parabolic), the element round trip, the leapfrog's order 2 on an oscillator, time reversibility, angular momentum to 10⁻¹², bounded energy error that scales as h², order 2 against DOP853, the closed orbit, the three laws, momentum and centre-of-mass conservation in the full problem, and the Hohmann numbers and limits. |
| **Practical activity** | **computational experiment**: a step-size and eccentricity study of the leapfrog against the `solve_ivp` oracle (position error, energy and angular-momentum conservation), a check of Kepler's third law from simulated periods over a decade of semi-major axes, and a Hohmann transfer sweep. |
| **Equipment** | None (a computer with the course's Python environment). |
| **Safety** | None beyond normal computer use. |
| **Evidence** | `two_body.png`; the step-size table (position error, energy error, angular-momentum error) with the observed order; the Kepler's-third-law table; the Hohmann table; a half-page note answering the engineering question. |
| **Assessment** | `assessments/phase-04/grading-rubric.md`; competencies 1–4. Pass: vis-viva and the Hohmann Δv derived by hand; all public tests pass; the 400 km orbit and LEO-to-GEO numbers reproduced to three figures; the note states the step that meets a stated energy tolerance and explains why the leapfrog's energy error does not grow. |
| **References** | `@openstax-university-physics-vol1` Ch. 13 "Gravitation": 13.1 "Newton's Law of Universal Gravitation", 13.2 "Gravitation Near Earth's Surface", 13.3 "Gravitational Potential Energy and Total Energy", 13.4 "Satellite Orbits and Energy", 13.5 "Kepler's Laws of Planetary Motion" (required-free; https://openstax.org/books/university-physics-volume-1/pages/13-4-satellite-orbits-and-energy). `@leishman-intro-aerospace-flight-vehicles` chapter "Astronautics, Space & Astrodynamics", from the book's contents at https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/ (supplementary-free; read for the aerospace framing of orbits and transfers). |
| **Extensions** | ★ Add the J₂ term of the Earth's oblateness to the acceleration and measure the regression of the node of a 400 km, 51.6° orbit. ★ Integrate a three-body problem (Earth, Moon, spacecraft) with the leapfrog and find a free-return trajectory by shooting. ★ Replace the fixed step by a time transformation (dt ∝ r^{3/2}) that keeps the leapfrog symplectic on very eccentric orbits. |
| **Time** | unmeasured |

## 1. The engineering problem

A launch vehicle drops a communications satellite into a 400 km circular parking orbit. The satellite must end up 35 786 km above the equator, where its orbital period equals one sidereal day and it hangs over one longitude. Three numbers decide whether the mission closes: the speed in the parking orbit (which the launcher had to provide), the velocity change Δv needed to get from there to geostationary orbit (which sets the propellant the satellite or an upper stage must carry), and the time the transfer takes (which sets how long the batteries must last before the solar arrays see the Sun). All three follow from one force law and two conservation laws, and this lesson derives them.

The second half of the problem is computational. Mission analysts do not stop at closed forms; they propagate the orbit numerically, with perturbations, for days or years. The integrator must not invent physics. Gravity conserves energy and angular momentum, so a numerical orbit whose energy drifts is spiralling in or out for no physical reason. `P03.05` showed RK4 losing energy slowly and Euler gaining it fast; `P03.07` built the symplectic integrators that do neither. Here you write the simplest of them, the leapfrog, by hand, and you hold it to the physics: angular momentum conserved to rounding, energy error bounded over ten orbits, second-order convergence of position against an independent solver.

Success criteria, written first: circular speed, period, escape speed and the LEO-to-GEO Hohmann Δv reproduced by hand and code to three figures; a leapfrog whose angular momentum error is below 10⁻¹² and whose energy error after ten orbits is no larger than after one; observed order 2.0 ± 0.1 in position against `solve_ivp`; T²/a³ from simulated periods equal to 4π²/μ within 10⁻⁴ over a factor of six in a.

## 2. Understand and derive

### 2.1 The law and the field

Newton's law of universal gravitation: two point masses M and m a distance r apart attract each other along the line joining them with force magnitude

**F = G M m / r²**,  G = 6.674 × 10⁻¹¹ N m² kg⁻².

The force on m is −GMm r̂/r², and the force on M is equal and opposite (Newton's third law, tested in the code). Dividing by m gives the acceleration of m, the **gravitational field** of M: **g(r) = −μ r̂ / r²** with **μ = GM**. For the Earth μ = 3.986004418 × 10¹⁴ m³/s², known to ten digits from satellite tracking, while G and M separately are known only to about five; that is why orbital mechanics is written in μ.

A uniform sphere attracts as if its mass were at its centre (OpenStax 13.2 explains why), so above the Earth's surface g(h) = μ/(R + h)². With R = 6371 km: g(0) = 9.820 m/s² (the measured 9.807 is lower because the Earth spins and bulges), and at 400 km g = 8.694 m/s², **88.5 % of the surface value**. Astronauts in orbit are not weightless because gravity is absent; they are falling freely along with their spacecraft.

### 2.2 Circular orbits and Kepler's third law

A satellite in a circular orbit of radius r has a centripetal acceleration v²/r supplied entirely by gravity: v²/r = μ/r², so

**v_c = √(μ/r)**,  **T = 2πr/v_c = 2π √(r³/μ)**.

At 400 km (r = 6771 km): **v_c = 7672.6 m/s, T = 92.41 min**. The period depends only on r, and T² ∝ r³: Kepler's third law for circles. Setting T to one sidereal day (86 164.1 s) and solving gives the geostationary radius r = (μT²/4π²)^{1/3} = **42 164.2 km** (altitude 35 793 km with the mean radius used here; 35 786 km above the equatorial radius), where v_c = 3074.7 m/s.

### 2.3 Energy: potential, total, escape

The work done by gravity in moving m from r₁ to r₂ is ∫ −GMm/r² dr = GMm(1/r₂ − 1/r₁), independent of the path, so gravity is conservative with potential energy **U(r) = −GMm/r**, zero at infinity. Per unit mass the **specific mechanical energy** of a satellite is

**ε = v²/2 − μ/r**,

and it is conserved along any unpowered trajectory. ε < 0 is bound; ε ≥ 0 escapes. Setting ε = 0 at radius r gives the **escape speed v_esc = √(2μ/r) = √2 v_c**: at 400 km, 10 850.7 m/s. For a circular orbit ε = μ/2r − μ/r = −μ/2r: the kinetic energy is half the magnitude of the potential energy, and the total is negative. At 400 km ε = −29.43 MJ/kg.

Notice the paradox that drives every orbit transfer: to go *up* to a larger orbit, ε must rise (toward zero), but the speed in the higher circular orbit is *lower*. Thrust adds energy; the satellite spends it climbing, and arrives slower.

### 2.4 Angular momentum and Kepler's second law

The gravitational force on m points along r, so the torque r × F about the centre is zero and the **specific angular momentum h = r × v** is conserved. Two consequences. The motion stays in the plane perpendicular to h. And in a short time dt the radius vector sweeps the thin triangle of area dA = ½|r × v dt|, so

**dA/dt = |h|/2 = constant**:

equal areas in equal times, **Kepler's second law**, which is conservation of angular momentum and nothing else. It holds for any central force, which is why a numerical method that conserves angular momentum exactly will satisfy the second law exactly, even with a wrong orbit.

### 2.5 Ellipses, the vis-viva relation and the orbital elements

Solving the equation of motion (OpenStax 13.5 states the result; `P03.02` used it through Kepler's equation) shows that bound orbits are ellipses with the centre of attraction at one focus: **Kepler's first law**. An ellipse is fixed by its semi-major axis a and eccentricity e; the closest point (periapsis; perigee for the Earth) is at r_p = a(1 − e) and the farthest (apoapsis) at r_a = a(1 + e).

Two conservation laws relate the speeds there. At both points the velocity is perpendicular to r, so h = r_p v_p = r_a v_a. Energy: ε = v_p²/2 − μ/r_p = v_a²/2 − μ/r_a. Eliminate v_a = v_p r_p/r_a from the energy equation:

v_p² (1 − r_p²/r_a²)/2 = μ(1/r_p − 1/r_a)  ⟹  v_p² = 2μ r_a / (r_p(r_p + r_a)).

Substitute back into ε: ε = μr_a/(r_p(r_p + r_a)) − μ/r_p = −μ/(r_p + r_a) = **−μ/(2a)**. The energy depends on the size of the orbit only, not on its shape. Then, since ε = v²/2 − μ/r holds everywhere on the orbit,

**v² = μ (2/r − 1/a)**,

the **vis-viva relation**: speed anywhere from the distance and the size alone. Its limits are the formulas above: r = a gives the circular speed, a → ∞ gives escape. And from h² = r_p² v_p² = μ a(1 − e²) we get the eccentricity back from a state:

**e = √(1 + 2εh²/μ²)**.

With ε and h computed from one position and velocity, a = −μ/2ε and e follow, and with them r_p and r_a. That is `orbital_elements`. The test orbit of this lesson, a = 7500 km and e = 0.1, has perigee 379 km and apogee 1879 km above the surface, period 6464.0 s, and speeds 8059.6 m/s at perigee and 6594.2 m/s at apogee; the code recovers a = 7500.0 km and e = 0.100000 from the perigee state.

Kepler's third law also holds for ellipses with a in place of r: **T = 2π√(a³/μ)** (Problem 2 sketches why). So T²/a³ = 4π²/μ = 9.9043 × 10⁻¹⁴ s²/m³ for every orbit about the Earth.

### 2.6 Two bodies, not one

The Earth also moves. With both masses free, the equations are m₁r̈₁ = Gm₁m₂(r₂ − r₁)/d³ and m₂r̈₂ = −Gm₁m₂(r₂ − r₁)/d³. Adding them: m₁r̈₁ + m₂r̈₂ = 0, so the total momentum is conserved and the **centre of mass** moves in a straight line at constant speed (`P04.03`). Subtracting the accelerations: the separation d = r₂ − r₁ obeys

**d̈ = −G(m₁ + m₂) d / |d|³**,

the one-body problem with μ = G(m₁ + m₂). Each body moves on an ellipse about the barycentre, scaled by m₂/(m₁ + m₂) and m₁/(m₁ + m₂). For the Earth and Moon (separation 384 400 km) the barycentre sits 4671 km from the Earth's centre, inside the Earth; the relative orbit has period 27.29 days (a circular orbit at the mean distance; the real sidereal month is slightly longer because the Moon's orbit is eccentric and perturbed by the Sun). For a satellite, m₂/m₁ ≈ 10⁻²¹, and the one-body problem with μ = GM_Earth is exact for practical purposes.

### 2.7 The Hohmann transfer

To move between two circular, coplanar orbits r₁ < r₂ with two impulsive burns, the cheapest way for ratios r₂/r₁ below 11.94 (where a bi-elliptic transfer through a very distant apoapsis starts to win) is the **Hohmann transfer**: an ellipse with periapsis r₁ and apoapsis r₂, so a_t = (r₁ + r₂)/2. Vis-viva gives the speeds on it, and the burns are the differences:

**Δv₁ = √(μ(2/r₁ − 1/a_t)) − √(μ/r₁)**,  **Δv₂ = √(μ/r₂) − √(μ(2/r₂ − 1/a_t))**,  **t = π√(a_t³/μ)**.

From 400 km to geostationary: a_t = 24 467.6 km, e_t = 0.7233, **Δv₁ = 2399.5 m/s, Δv₂ = 1457.2 m/s, total 3856.7 m/s, time of flight 5.29 h**. The energy bookkeeping: ε goes from −29.434 MJ/kg in the parking orbit to −8.145 MJ/kg on the transfer ellipse after the first burn; the second burn adds the rest to reach −μ/(2r_GEO). (A real geostationary transfer also removes the launch site's orbital inclination, which costs more; that is the reason launches from near the equator are prized.)

### 2.8 The leapfrog, and why it conserves what it does

The **leapfrog** (kick-drift-kick, also called velocity Verlet; `P03.07`) advances (r, v) by a step h as

v_{1/2} = v + (h/2) a(r);  r′ = r + h v_{1/2};  v′ = v_{1/2} + (h/2) a(r′).

One force evaluation per step (the last kick's a(r′) is the next step's first). It is second order, time-reversible (run it backward from (r′, −v′) and you land on (r, −v) to rounding: a test checks this), and symplectic, which `P03.07` showed means it exactly conserves a *modified* energy close to the true one, so the true energy error oscillates with an amplitude ∝ h² and does not grow. For a central force it does better still: each kick changes v along r, which leaves r × v unchanged, and each drift moves r along v, which also leaves r × v unchanged. So **the leapfrog conserves angular momentum exactly**, to rounding, at any step size, and therefore satisfies Kepler's second law exactly. Energy and the orbit's orientation (the perigee direction) are where its error shows.

### 2.9 Where the model fails

- **Oblateness.** The Earth's equatorial bulge (J₂ ≈ 1.08 × 10⁻³) makes the force non-central; h is no longer conserved, the orbit plane precesses (−5.0° per day at 400 km and 51.6° inclination, from dΩ/dt = −(3/2) n J₂ (R_eq/a)² cos i) and the perigee rotates. It is the first correction any real propagator adds.
- **Drag.** Below roughly 500 km the residual atmosphere (`P04.07`) removes energy; the orbit shrinks and, by the paradox of 2.3, the satellite *speeds up* as it decays.
- **Third bodies.** At geostationary altitude the Moon and Sun slowly tilt the orbit plane, by a fraction of a degree per year; station-keeping burns undo it.
- **The step.** A fixed step resolves perigee and apogee equally; on a very eccentric orbit the perigee passage needs far more steps than the rest (part 4). A variable step breaks the leapfrog's symplectic property unless it is done through a time transformation (★ extension).

## 3. Build the computational model

Open `starter/main.py`. The constants are given; every function body is yours.

### Step 1: the law and the circular-orbit formulas

```python
def gravity_acceleration(r, mu=MU_EARTH):
    r = np.asarray(r, dtype=float)
    rn = np.linalg.norm(r, axis=-1, keepdims=True)
    return -mu * r / rn ** 3
```

`axis=-1, keepdims=True` lets the same function take one vector or a stack of shape (n, 2). Then `gravitational_force` (the pair force, for the third-law test), `g_at_altitude`, `circular_speed`, `orbital_period`, `escape_speed`, `vis_viva` from 2.2–2.5.

### Step 2: energy, angular momentum, elements

`specific_energy(r, v)` and `specific_angular_momentum(r, v)` must accept stacks: use `np.sum(v * v, axis=-1)` and `r[..., 0] * v[..., 1] - r[..., 1] * v[..., 0]`. `orbital_elements` computes ε and h, then a = −μ/2ε, e = √(1 + 2εh²/μ²), p = h²/μ, r_p = p/(1 + e) and r_a = p/(1 − e) (infinite for e ≥ 1). `perigee_state(a, e)` puts the satellite at (a(1 − e), 0) moving along +y at the vis-viva perigee speed.

### Step 3: the leapfrog

```python
def leapfrog(accel, r0, v0, h, n_steps):
    r, v = np.array(r0, float), np.array(v0, float)
    R, V = np.empty((n_steps + 1, r.size)), np.empty((n_steps + 1, v.size))
    R[0], V[0] = r, v
    a = accel(r)
    for k in range(n_steps):
        v_half = v + 0.5 * h * a
        r = r + h * v_half
        a = accel(r)
        v = v_half + 0.5 * h * a
        R[k + 1], V[k + 1] = r, v
    return h * np.arange(n_steps + 1), R, V
```

`propagate_leapfrog(a, e, n_orbits, steps_per_orbit)` wraps it from perigee with h = T/steps_per_orbit. `propagate_solve_ivp` is the oracle: `two_body_rhs` and DOP853 at `rtol=1e-12`, sampled at fixed times. `energy_drift` returns (ε(t) − ε(0))/|ε(0)|.

### Step 4: Kepler's laws on the simulated orbit

`areal_velocity(R, dt)` is |r_k × r_{k+1}|/(2 dt), the area of each thin triangle per unit time. `period_from_trajectory(t, R)` finds the samples where |r| is a local minimum (perigee passages), refines each time by a parabola through the three samples around it, and averages the intervals; it raises `ValueError` if there are fewer than two passages. `keplers_third_law_ratio(a, T)` is T²/a³.

### Step 5: two free bodies and the transfer

`two_body_full_rhs` has eight state components; `simulate_two_body` starts the bodies in the centre-of-mass frame (r₁ = −m₂d/M, r₂ = m₁d/M and the same split for velocity) and returns positions, velocities, the centre of mass, the total momentum and the relative motion. `hohmann_transfer(r1, r2)` is 2.7, with absolute values so that it works downward too.

Running `main()` prints, in part,

```
g at the surface 9.820 m/s^2, at 400 km 8.694 m/s^2 (88.5 % of surface)
circular orbit at 400 km: v = 7672.6 m/s, T = 92.41 min, escape speed 10850.7 m/s = sqrt(2) x circular
geostationary radius 42164.2 km (altitude 35793.2 km), v = 3074.7 m/s
Hohmann 400 km -> GEO: dv1 = 2399.5 m/s, dv2 = 1457.2 m/s, total 3856.7 m/s, time of flight 5.29 h
leapfrog vs DOP853 (rtol 1e-12) after one orbit; energy and angular-momentum error over ten orbits
 steps/orbit  pos err 1 orbit (m)   max |dE/E| 10 orbits   max |dh/h| 10 orbits
         250            1.367e+04              6.760e-05              3.927e-15
         500            3.418e+03              1.688e-05              6.171e-15
        1000            8.546e+02              4.220e-06              1.122e-14
        2000            2.136e+02              1.055e-06              2.146e-14
observed order of the leapfrog position error: 2.00
```

and writes `two_body.png` (the orbit against the oracle, the energy error over ten orbits, the convergence plot).

### Verify

Run `python -m pytest curriculum/phase-04/lesson-05-gravitation-keplers-laws`. The 17 tests check, in order: the inverse-square law and its direction, on single vectors and stacks; Newton's third law for the pair force; g(0) = μ/R², g(R) = g(0)/4; v_c = √(μ/r), T = 2πr/v_c, v_esc = √2 v_c, v_c²/r = g(h); vis-viva at r = a, at perigee, at a → ∞; ε = −μ/2a and h = √(μa(1 − e²)) at perigee, elements recovered from any point of the orbit, and an escaping state flagged; the leapfrog's order 2 on an oscillator; time reversibility; angular momentum to 10⁻¹² and an energy error over ten orbits no larger than over one, falling four-fold when the step halves; order 2 against DOP853; the solve_ivp orbit closing after one period; the ellipse (r_p, r_a, sum of focal distances = 2a); constant areal velocity on an e = 0.5 orbit; T²/a³ from simulated periods; the full two-body problem (momentum, centre of mass, relative circle, barycentre radius); the Hohmann numbers, limits and reverse transfer; `main()` writing the figure.

### Compare

`scipy.integrate.solve_ivp` with DOP853 at rtol 10⁻¹² is the oracle: an eighth-order adaptive method whose error on one orbit is far below the leapfrog's. Two differences are worth noticing. First, the leapfrog's angular momentum error (10⁻¹⁴) is *smaller* than DOP853's at that tolerance; a cheap second-order method beats a superb eighth-order one on the quantity it is built to keep. Second, DOP853's energy error drifts slowly, like RK4's in `P03.05`, while the leapfrog's does not drift at all: after ten orbits at 1000 steps per orbit its maximum is 4.220 × 10⁻⁶, exactly the first orbit's maximum. For a one-orbit position, DOP853 wins easily; for a ten-thousand-orbit study of an orbit's long-term behaviour, the symplectic method's honesty about energy is worth more than its accuracy per step.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction down before running it.

1. **Step size.** Propagate the test orbit (a = 7500 km, e = 0.1) for ten orbits at 250, 500, 1000, 2000 steps per orbit. Record the one-orbit position error against the oracle, the maximum energy error over ten orbits and the maximum angular-momentum error. *Predict:* position error and energy error both fall by 4 per halving of the step; angular-momentum error stays at rounding level (it grows slightly with the number of steps, because rounding errors accumulate).
2. **Boundedness.** At 1000 steps per orbit, compare the maximum energy error in the first orbit with that over ten. *Predict:* equal. Repeat with RK4 from `P03.05` at the same cost (250 steps per orbit, four evaluations each) and predict a monotone negative drift.
3. **Eccentricity.** At 1000 steps per orbit, propagate a = 15 000 km with e = 0.1, 0.3 and 0.5 (the perigee stays above the surface for all three). *Predict* that the energy error grows steeply with e, because a fixed step spends most of its samples far from perigee where nothing happens.
4. **Kepler's third law.** For a = 7000, 10 000, 20 000 km and the geostationary radius, at e = 0.2 and 1000 steps per orbit, find the period from the simulated perigee passages and compute T²/a³. *Predict* 4π²/μ = 9.9043 × 10⁻¹⁴ s²/m³ to within the leapfrog's period error.
5. **Transfers.** Tabulate the Hohmann Δv₁, Δv₂, total and time of flight from 400 km to r₂ = 10 000, 20 000, 42 164, 100 000 and 384 400 km. *Predict* that the total Δv rises and then falls (as a function of r₂/r₁ it peaks at 15.58), and that the Moon's distance costs less total Δv from LEO than you might guess.

Raw data format: CSV files `steps_per_orbit, pos_err_1orbit_m, max_rel_energy_err_10, max_rel_h_err_10` (1–2), `e, steps_per_orbit, max_rel_energy_err` (3), `a_m, T_sim_s, T_kepler_s, T2_over_a3` (4), `r2_m, dv1_mps, dv2_mps, total_mps, tof_s` (5); plus `two_body.png`.

What you should find for items 1, 2 and 4 (from `main()`): the table in part 3, with observed order 2.00, energy errors in the ratio 4.00 per halving and angular-momentum errors of 4 × 10⁻¹⁵ to 2 × 10⁻¹⁴; at 1000 steps the energy error is 4.220 × 10⁻⁶ after one orbit and 4.220 × 10⁻⁶ after ten; the simulated periods 5828.60, 9952.16, 28 148.97 and 86 165.39 s give T²/a³ = 9.9046 × 10⁻¹⁴ in every case. Items 3 and 5 are yours.

## 5. Validate and improve

| quantity | theory (§2) | code | difference | cause |
|---|---|---|---|---|
| v_c, T at 400 km | 7672.6 m/s, 92.41 min | same | 0 | closed form both sides |
| LEO → GEO Hohmann | 3856.7 m/s, 5.29 h | same | 0 | closed form |
| elements from the perigee state | a = 7500 km, e = 0.1 | 7500.0 km, 0.100000 | < 10⁻¹⁰ | — |
| leapfrog position order | 2 | 2.00 | 0 | — |
| leapfrog energy error, 10 orbits vs 1 | equal (bounded) | 4.220 × 10⁻⁶ both | 0 | symplectic |
| leapfrog angular momentum | exact | ≤ 2 × 10⁻¹⁴ relative | rounding | — |
| areal velocity spread | 0 | 3.3 × 10⁻¹⁴ | rounding | exact h conservation |
| period of the test orbit | 6464.02 s | 6464.07 s | +0.048 s (7 × 10⁻⁶) | leapfrog phase error ∝ h² |
| T²/a³ | 9.9043 × 10⁻¹⁴ | 9.9046 × 10⁻¹⁴ | 3 × 10⁻⁵ | twice the period error |
| Earth–Moon barycentre | 4671 km from the centre | 4671 km | 0 | — |
| two-body momentum | conserved | 1.8 × 10⁻¹⁵ of the Moon's momentum | rounding | — |

Two rows depart from theory, and both are the leapfrog's one real error: it keeps the right energy (on average) and the right angular momentum, but its orbit runs slightly slow, a *phase* error. The simulated period is 0.048 s long at 1000 steps per orbit, which is 7.4 × 10⁻⁶, and T²/a³ inherits twice that relative error. Halve the step and the discrepancy falls four-fold: 0.012 s at 2000 steps per orbit (Problem 4 measures it on another orbit). So the right way to read the table is: conserved quantities exact or bounded, timing off by an amount you can predict and size. A position error of 855 m after one orbit at 1000 steps is almost entirely along-track: the satellite is on the right ellipse, a little late.

What to improve, and what not. For a mission timeline, 855 m after one orbit and growing linearly is not good enough; either refine the step (2000 steps: 214 m) or move to a higher-order symplectic method (Yoshida 4, `P03.07`). But before buying digits, compare with the model error: J₂ alone shifts the ground track of a 400 km, 51.6° satellite by about 30 km per orbit at the equator crossing (0.32° of nodal regression per orbit) relative to the point-mass ellipse. The first improvement to this model is physics (J₂, drag), not numerics.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of vis-viva, ε = −μ/2a, escape speed and the Hohmann Δv; `two_body.png`; the five experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) answering the engineering question: the parking-orbit speed and period, the LEO-to-GEO Δv split and time, and the leapfrog step that keeps the relative energy error below 10⁻⁶ over ten orbits of the test orbit, with the evidence (the measured error at one step and the h² law).

### Adversarial self-checks

- **What if the satellite's mass doubles?** Nothing in the orbit changes: the acceleration is μ/r², independent of m. The Δv is the same; the propellant mass doubles (the rocket equation of `P04.03` scales with the vehicle).
- **What if the step doubles?** Position and energy errors quadruple (13.7 km and 6.8 × 10⁻⁵ at 250 steps per orbit), angular momentum stays exact. If the step approaches the perigee-passage time scale the leapfrog loses the orbit altogether; the energy error then stops being small, the first warning.
- **Which assumption fails first in the real world?** For the 400 km parking orbit, J₂ (degrees per day of nodal regression), then drag (decay over months). For geostationary orbit, lunar and solar perturbations. A simulation that matches the two-body ellipse perfectly is still wrong about a real satellite by kilometres a day.
- **What measurement could falsify the model?** Radar or GNSS tracking of a satellite over a day: the point-mass prediction misses the observed node by about 5° per day at 51.6° inclination. Adding J₂ closes most of the gap; what remains is drag and higher harmonics.
- **Does the result survive a change of solver?** The leapfrog at 2000 steps and DOP853 agree to 214 m after one orbit and on ε and h to 10⁻⁶ and 10⁻¹⁴; the closed-form v_c, T and Δv need no solver at all and match both.
- **Which decision changes if the energy tolerance tightens tenfold?** The step: energy error ∝ h², so 10× tighter needs √10 = 3.2× more steps per orbit. The method choice changes only if the cost matters, in which case a fourth-order symplectic method is the answer, not RK4 (whose energy drifts).

## Check yourself

1. Write g at altitude h in terms of μ and R. What fraction of surface g acts on a satellite at 400 km?
2. Give the circular-orbit speed and the escape speed at radius r, and their ratio.
3. Write the specific orbital energy of an orbit with semi-major axis a. Is a satellite with ε = +1 MJ/kg bound?
4. State the vis-viva relation and its two limiting cases.
5. Kepler's second law is a consequence of which conservation law? Why does the leapfrog satisfy it exactly?
6. Two satellites orbit the Earth with the same semi-major axis, one circular and one with e = 0.5. Which has the longer period?
7. Why does a satellite that thrusts forward to raise its orbit end up moving more slowly?
8. The leapfrog's energy error after ten orbits equals its error after one. What would a monotone growth tell you about your code?
