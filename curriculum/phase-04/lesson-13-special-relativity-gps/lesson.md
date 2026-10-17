# P04.13  Special Relativity for GPS

> A GPS satellite's clock runs fast by 38 microseconds a day: slow by 7 because it moves, fast by 46 because it is high. Left uncorrected, that is 11 kilometres of position error by tomorrow. Relativity is not a curiosity in navigation; it is a line item in the error budget.

| | |
|---|---|
| **Learning outcomes** | 1. Derive time dilation, length contraction and the Lorentz transformation from the two postulates, show the invariance of the spacetime interval, and derive velocity addition, the relativistic Doppler shift and the energy–momentum relation from them. 2. State the weak-field gravitational time dilation ΔΦ/c² and say when it is valid. 3. Compute the clock-rate budget of a satellite on any circular orbit (velocity term, gravitational term, ground-clock terms), reproduce the GPS factory offset −4.465 × 10⁻¹⁰ and the net +38.6 µs/day, and check it against an independent form. 4. Implement and test the eccentricity correction, the Sagnac term and the accumulated pseudorange error, and decide for a given mission (LEO, GPS, GEO, Mars) which terms must be corrected to meet a stated position accuracy. |
| **Prerequisites** | `P04.05` (circular-orbit speed, vis-viva, the gravitational potential), `P04.02` (kinetic energy, work), `P03.02` (Kepler's equation, used for the eccentricity correction), `P02.16` (what a frequency offset does to a measured time), `P01.07` (cross products for the Sagnac geometry). Skill check: you can write v = √(GM/a) for a circular orbit and expand (1 − x)^−½ to first order in x. |
| **Engineering question** | A GPS satellite carries a 10.23 MHz atomic clock on a 26 562 km orbit. By how much, and in which direction, must its frequency be offset before launch so that a receiver on the ground sees 10.23 MHz, what periodic corrections remain, and how fast does the position error grow if any of them is left out? |
| **Assumptions and limits** | Special relativity: inertial frames, flat spacetime, speeds below c. Gravity enters only through the weak-field rate ΔΦ/c² (valid when GM/(rc²) ≪ 1: here 7 × 10⁻¹⁰, so the next-order terms are 10⁻¹⁹ and ignorable). Earth is a point mass plus J₂ for the geoid; the Sun's and Moon's potentials are the same for satellite and ground clock to the precision needed here. The ground clock sits on the geoid (sea level); a clock at altitude h gains gh/c² extra (about 0.94 ns/day per 100 m). The model fails for strong fields (neutron stars, black holes), for speeds where γ − 1 is not ≈ v²/2c² (above about 0.1c the series needs its next term), and for clocks whose own stability (10⁻¹³ per day for a GPS rubidium clock) is coarser than the effect being corrected. |
| **Analytical work** | The light clock and γ; the Lorentz transformation from linearity and the invariance of c; the interval; velocity addition; the Doppler factor and its transverse (second-order) form; E = γmc², p = γmv and E² − (pc)² = (mc²)²; the weak-field rate; the GPS budget by hand: −v²/2c², (GM/c²)(1/R − 1/a), the geoid terms, the net, the factory frequency; the crossover radius 3R/2; the eccentricity term from vis-viva and the orbit average of 1/r; the Sagnac term from the receiver's motion during signal flight. |
| **Computational deliverable** | `starter/main.py`: `gamma`, `gamma_minus_one`, `time_dilation_rate`, `length_contraction`, `lorentz_transform`, `lorentz_inverse`, `spacetime_interval`, `velocity_addition`, `relativistic_doppler`, `relativistic_energy`, `kinetic_energy`, `momentum`, `circular_orbit_speed`, `gravitational_rate`, `geoid_potential_ratio`, `gps_clock_budget` (a `ClockBudget` with the four terms), `gps_frequency_offset_oracle`, `factory_clock_frequency`, `crossover_radius`, `eccentricity_correction`, `sagnac_correction`, `max_sagnac`, `position_error_without_correction`, `pass_doppler`. Tests prove: γ(0.6c) = 1.25 and γ − 1 = v²/2c² + 3v⁴/8c⁴ to 10⁻¹² at 3.9 km/s; the interval is invariant under random boosts and the inverse boost round-trips; two boosts compose by velocity addition; c ⊕ v = c; the Doppler product f_app f_rec = f₀² and the transverse shift f₀/γ; E² − (pc)² = (mc²)²; the kinetic energy at 7.7 km/s is ½mv² to 10⁻⁹; the weak-field rate at R_E is −6.95 × 10⁻¹⁰ and follows 1/r; the velocity term lies in (−7.5, −6.9) µs/day, the gravitational term in (45, 46.3), the net in (38, 39) and equals 4.4647 × 10⁻¹⁰ to 0.02 % and the oracle to 10⁻¹²; the crossover is at 1.5 R_E; the eccentricity amplitude is 45.8 ns at e = 0.02 and averages to zero; the Sagnac term is zero overhead, antisymmetric east–west and 133 ns at the horizon; the error grows linearly to 11.6 km/day; the figure is written. |
| **Practical activity** | **computational experiment**: a sweep of the clock budget over altitude from 200 km to GEO (where does the net change sign; how do ISS, GPS, GEO compare), the accumulated error against time for each missing term, the eccentricity term over an orbit, and the Sagnac term over a sky pass. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `gps_relativity.png`, the hand-derived budget with each term to three figures, the four sweep tables, and a half-page note stating which corrections a 1 m, a 10 m and a 1 km navigation service must apply, with the numbers. |
| **Assessment** | `assessments/phase-04/grading-rubric.md`; competencies 1–4. Pass: the Lorentz transformation and the interval invariance are derived; the budget is computed by hand and agrees with the code to three figures; all public tests pass; the note names the terms a 1 m service must keep and gives the time after which each omitted term alone exceeds 1 m. |
| **References** | `@openstax-university-physics-vol3` Ch. 5 "Relativity": 5.1 Invariance of Physical Laws; 5.3 Time Dilation; 5.4 Length Contraction; 5.5 The Lorentz Transformation; 5.7 Doppler Effect for Light; 5.9 Relativistic Energy (required-free): the postulates, the derivations, the muon and Hafele–Keating examples (the chapter does not do GPS; this lesson does). `@openstax-university-physics-vol1` Ch. 13 "Gravitation": 13.7 Einstein's Theory of Gravity (required-free): gravitational time dilation, qualitatively. `@mit-ocw-8-033-relativity` Session 23 "General Relativity and Black Holes", topics "Gravitational Redshift" and "Application to the GPS System" (supplementary-free): the same budget as slides. `@ashby-relativity-in-gps` "Relativistic Effects on Satellite Clocks" ("Satellite orbits", "The eccentricity correction"), "Reference Frames and the Sagnac Effect", "GPS Coordinate Time and TAI" (supplementary-free): the authoritative treatment; link only, read for the derivation of the factory offset. |
| **Extensions** | ★ Add the Sun's and Moon's tidal potential differences between satellite and ground and show each is below 10⁻¹⁵. ★ Integrate the instantaneous rate −v²/2c² − GM/rc² around an e = 0.1 orbit with `solve_ivp` and recover −3GM/2ac² as the mean. ★ Compute the Sagnac term for a receiver at latitude 45° for a full sky pass from a real GPS almanac. |
| **Time** | unmeasured |

## 1. The engineering problem

A GNSS receiver works by timing. Each satellite broadcasts the time from its own clock; the receiver measures when the signal arrives, multiplies the delay by c, and gets a range, 1 m for every 3.3 ns. Four such ranges fix position and the receiver's own clock error. The whole system therefore rests on every satellite clock keeping the same time as the clocks on the ground, to a few nanoseconds. Two things stop it from doing so. The satellite moves at 3.87 km/s, and special relativity says a moving clock runs slow by v²/2c² = 8.3 × 10⁻¹¹, which is 7.2 µs a day. It also sits 20 200 km up, where the gravitational potential is less negative than at the ground, and general relativity says a clock higher in a potential runs fast by ΔΦ/c² = 5.3 × 10⁻¹⁰, which is 45.7 µs a day. The two have opposite signs and do not cancel: the net is +38.6 µs/day, and since 1 µs is 300 m, an uncorrected constellation would be 11.6 km wrong after one day and useless after an hour.

The fix is mostly done before launch: the 10.23 MHz clock is set to 10.229 999 995 43 MHz, so that seen from the ground it reads 10.23 MHz exactly. Two periodic terms remain, and the receiver must compute them. The orbit is slightly eccentric, so the speed and the height, and therefore both rate terms, vary around the orbit; the correction is −2√(GMa) e sin E / c², 46 ns peak-to-peak for e = 0.02, which is 14 m. And the Earth rotates under the signal while it is in flight, so a receiver that computes ranges in the rotating frame must add the Sagnac term, up to 133 ns (40 m) for a satellite on the horizon. Each of these is a number this lesson derives by hand and reproduces in code; together they are the simplest honest answer to "does relativity matter in engineering": yes, to the metre, every second, in every phone.

Success criteria, written first: γ and the Lorentz transformation derived from the postulates and verified by the interval invariance; the budget computed by hand to three figures and by code to the published value 4.4647 × 10⁻¹⁰ within 0.02 %; the sweep over altitude showing the sign change at 1.5 R_E; the periodic terms with their amplitudes; and a decision table saying which terms a 1 m, 10 m and 1 km service must keep.

## 2. Understand and derive

### 2.1 The two postulates and the light clock

Special relativity rests on two statements (OpenStax 5.1): the laws of physics are the same in every inertial frame, and the speed of light in vacuum is c in every inertial frame, whatever the motion of the source. The second is the strange one, and everything below follows from it.

A **light clock** is a pulse bouncing between two mirrors a distance L apart, perpendicular to the direction of motion: one tick is 2L/c in the clock's rest frame. Watch the same clock move past at speed v. In your frame the pulse travels a diagonal: in half a tick Δt/2 it goes a distance cΔt/2 while the mirror moves vΔt/2, so (cΔt/2)² = L² + (vΔt/2)², giving Δt = (2L/c)/√(1 − v²/c²). The moving clock's tick, measured by you, is longer by the **Lorentz factor**

**γ = 1/√(1 − v²/c²)**,   (1)

and the clock runs slow: **time dilation** (5.3). Its rate relative to yours is 1/γ − 1, which for v ≪ c is −v²/2c² (expand (1 − x)^½ ≈ 1 − x/2). Numbers: γ(0.6c) = 1.25; for a muon at 0.995c γ = 10.0, which is why a particle with a 2.2 µs lifetime reaches the ground from 15 km up (it covers 6.6 km per dilated lifetime, not 0.66). For a GPS satellite at 3873.8 m/s, γ − 1 = 8.349 × 10⁻¹¹, and the clock loses 8.349 × 10⁻¹¹ × 86 400 s = 7.21 µs per day. An airliner at 250 m/s loses 30 ns per day. The effect is tiny at any engineering speed, which is exactly why it only matters where the clocks are good enough to see it.

A numerical point the code enforces: at 4 km/s, `1/sqrt(1 - v**2/c**2) - 1` keeps six significant digits of a quantity whose true value is 8.349 × 10⁻¹¹, because the subtraction cancels ten of sixteen. The identity γ − 1 = β²/(√(1 − β²)(1 + √(1 − β²))), β = v/c, has no cancellation and is what `gamma_minus_one` uses; a test checks it against the series β²/2 + 3β⁴/8 to 10⁻¹².

### 2.2 Length contraction and the Lorentz transformation

A rod of rest length L₀ moving along its length is measured shorter: L = L₀/γ (5.4). The argument is the light clock turned sideways: a pulse sent along the rod and back takes 2γ²L/c in your frame (the forward leg takes L/(c − v), the return L/(c + v)), and must equal γ × 2L₀/c by time dilation, so L = L₀/γ. At GPS speed a 2 m satellite is contracted by 0.17 nm.

The general rule is the **Lorentz transformation** (5.5). Frame S′ moves at v along x. The transformation must be linear (a uniform motion must map to a uniform motion), must reduce to x′ = x − vt at low speed, and must keep x = ct ↔ x′ = ct′. Writing x′ = A(x − vt) and, by symmetry, x = A(x′ + vt′), and substituting a light ray x = ct, x′ = ct′ into both gives A² = 1/(1 − v²/c²): A = γ. Eliminating x′ then gives t′:

**t′ = γ(t − vx/c²),   x′ = γ(x − vt)**.   (2)

The inverse is the same with −v. Three consequences, each a test:

- The **interval** s² = c²Δt² − Δx² is the same in every frame (substitute (2); the cross terms cancel because γ²(1 − v²/c²) = 1). It is the invariant that replaces "time" and "distance" separately.
- Two boosts along x compose into one boost whose speed is the **velocity addition** w = (u + v)/(1 + uv/c²): never above c, and c ⊕ v = c. At u = v = 0.5c the result is 0.8c; at 100 and 200 m/s it is 300 m/s to one part in 10¹².
- Simultaneity is relative: events with Δt = 0 and Δx ≠ 0 in S have Δt′ = −γvΔx/c² ≠ 0 in S′. For GPS this is the Sagnac term in disguise (2.7).

### 2.3 The Doppler effect for light

A source of rest frequency f₀ moving at speed v is received at (5.7)

**f = f₀ / (γ(1 − (v/c) cos θ))**,   (3)

θ the angle between the source velocity and the line to the receiver, in the receiver's frame. Head-on (θ = 0): f = f₀√((1 + β)/(1 − β)); receding: f₀√((1 − β)/(1 + β)); their product is f₀². Transverse (θ = 90°): f = f₀/γ, the pure time-dilation shift, second order in v/c, with no classical counterpart. For the GPS L1 carrier at 1575.42 MHz a radial speed of 800 m/s gives a first-order shift of ±4204 Hz, which the receiver tracks; the second-order shift is −f₀(γ − 1)/γ = −0.132 Hz, constant over the pass, and it is the velocity term of the clock budget seen as a frequency instead of a time. The receiver never sees the two separately: the first-order term is antisymmetric over a symmetric pass and averages out, the second-order term does not, and that is what the factory offset removes.

### 2.4 Relativistic energy and momentum

With p = γmv and E = γmc² (5.9), direct algebra gives **E² − (pc)² = (mc²)²**, an invariant like the interval. The kinetic energy K = (γ − 1)mc² reduces to ½mv² for v ≪ c (one more term: ½mv²(1 + 3v²/4c²)); at 7.7 km/s the correction is 5 × 10⁻¹⁰, at 0.1c it is 0.76 %, at 0.9c K is 3.2 times ½mv². The code's `kinetic_energy` uses `gamma_minus_one` for the same cancellation reason as before. Where this matters in aerospace: the electrons and protons of the radiation belts (`P04.14`) are relativistic (a 1 MeV electron has γ = 2.96 and v = 0.941c), while the ions in an electric thruster's exhaust at 30 km/s have γ − 1 = 5 × 10⁻⁹ and are not.

### 2.5 Gravitational time dilation, weak field

General relativity adds (OpenStax Vol. 1, 13.7): a clock at rest where the gravitational potential is Φ runs, relative to a clock far away where Φ = 0, at rate 1 + Φ/c² (Φ negative, so slower). Between two clocks at Φ₁ and Φ₂ the rate difference is (Φ₂ − Φ₁)/c². This is the first term of an expansion in Φ/c², valid when |Φ|/c² ≪ 1; at the Earth's surface GM/(Rc²) = 6.95 × 10⁻¹⁰, so the next term is 10⁻¹⁹ and no clock can see it. A quick plausibility argument, not a derivation: a photon climbing a height h loses energy mgh with m = hf/c², so its frequency falls by gh/c²; a clock is a frequency, so a clock lower down is seen to tick slower by the same factor. That gives the useful local form gh/c² = 1.09 × 10⁻¹⁶ per metre, about 0.94 ns per day per 100 m: a clock on the second floor gains on the one in the basement, measurably with today's optical clocks.

### 2.6 The GPS budget, by hand

Take the satellite on a circular orbit of radius a = 26 561.75 km (period 43 082 s, half a sidereal day), the ground clock on the equator at R = 6378.137 km, GM = 3.986 004 418 × 10¹⁴ m³/s². Write the fractional rate of the satellite clock relative to the ground clock as the sum of four terms:

| term | expression | value | µs/day |
|---|---|---|---|
| velocity | −v²/2c², v = √(GM/a) = 3873.8 m/s | −8.3485 × 10⁻¹¹ | −7.213 |
| gravity | (GM/c²)(1/R − 1/a) | +5.2838 × 10⁻¹⁰ | +45.652 |
| ground rotation | +ω²R²/2c² (the ground clock moves at ωR = 465 m/s) | +1.2034 × 10⁻¹² | +0.104 |
| ground J₂ | +(GM/Rc²)(J₂/2) (the equatorial geoid is deeper than GM/R) | +3.764 × 10⁻¹³ | +0.033 |
| **net** | | **+4.4647 × 10⁻¹⁰** | **+38.575** |

The first two are the whole story to 0.4 %; the simple model (sphere, no rotation) gives 38.44 µs/day. The last two exist because the ground clock is not at rest in an inertial frame and not on a sphere: both make the *ground* clock slightly slower, so both add to the satellite's apparent gain. With them, the net is Ashby's published −4.4647 × 10⁻¹⁰ (sign: the satellite clock must be *slowed*). A second way to get the same number is the independent oracle the tests use: for a circular orbit v²/2 = GM/2a, so the satellite's rate relative to infinity is −v²/2c² − GM/ac² = −3GM/2ac², and the ground clock's is Φ₀/c² with Φ₀ = −(GM/R)(1 + J₂/2) − ω²R²/2 = −6.9693 × 10⁻¹⁰ c²; the difference is −2.5046 × 10⁻¹⁰ + 6.9693 × 10⁻¹⁰ = +4.4647 × 10⁻¹⁰. The two computations share no code and agree to round-off.

So the clock is set before launch to 10.23 MHz × (1 − 4.4647 × 10⁻¹⁰) = **10 229 999.995 43 Hz**, 4.57 mHz low, and the ground sees 10.23 MHz. Without it, the broadcast time drifts by 38.6 µs/day and the pseudorange by c × 4.4647 × 10⁻¹⁰ = 0.134 m every second: 8 m after a minute, 482 m after an hour, **11.56 km after a day**.

**Where the sign flips.** Keeping only the two main terms, −GM/2ac² + (GM/c²)(1/R − 1/a) = 0 when 1/R = 3/2a, i.e. **a = 3R/2**, altitude 3189 km: below it the velocity term wins and the satellite clock runs slow (the ISS at 400 km: −28.27 + 3.55 + 0.14 = −24.58 µs/day, the last term the ground clock's), above it gravity wins (GEO: −4.54 + 50.99 + 0.14 = +46.58 µs/day). The code finds the crossover by bisection as a check on the algebra.

### 2.7 The periodic terms

**Eccentricity.** On an elliptical orbit both v and r vary. From vis-viva, v²/2 = GM/r − GM/2a, so the instantaneous rate relative to infinity is −(v²/2 + GM/r)/c² = −(2GM/r − GM/2a)/c². Its orbit average is −3GM/2ac² because the time average of 1/r over a Kepler orbit is exactly 1/a (Problem 5 proves it); that average is what the factory offset removes. The part that varies, −2GM(1/r − 1/a)/c², integrates over time to a clock offset of

**Δt_ecc = −(2√(GMa)/c²) e sin E**,   (4)

E the eccentric anomaly (`P03.02`), which Ashby derives in "The eccentricity correction". For a = 26 562 km the amplitude is 2.29 × 10⁻⁶ e seconds: **45.8 ns at e = 0.02**, i.e. 13.7 m of range, 57 ns at e = 0.025. It is zero at perigee and apogee, negative just after perigee (the clock has just run slowest), and averages to zero over an orbit. The receiver computes it from the broadcast orbit; it is the one relativistic term that lives in receiver firmware.

**Sagnac.** The receiver computes the satellite position at the transmit time and its own position at the receive time, both in the Earth-fixed rotating frame. But the signal travels straight in the inertial frame, and during its flight of at least 67 ms the receiver moves ωR Δt ≈ 31 m or more eastward. The extra path is the component of that motion along the line of sight, (r_r − r_s)/ρ · (ω ẑ × r_r) × ρ/c, which simplifies to

**Δt_Sagnac = ω (x_s y_r − y_s x_r)/c² = ω ẑ·(r_s × r_r)/c²**,   (5)

equivalently 2ωA/c² with A the area of the triangle Earth-centre–satellite–receiver projected on the equatorial plane. It is zero for a satellite directly overhead, antisymmetric between a satellite to the east and one to the west, and largest for a satellite on the equatorial horizon of an equatorial receiver, where cos θ = R/a and Δt = ωaR sin θ/c² = **133.4 ns, 40 m**. Ignoring it is a 40 m error that depends on where in the sky each satellite is, so it does not even average out between satellites.

### 2.8 Where the model stops

- Strong fields: the weak-field rate is the first term in GM/rc²; the next is 10⁻¹⁹ here and irrelevant, but near a neutron star it is not.
- The Sun and Moon: their potentials differ between satellite and ground only through the tidal term, 6.8 × 10⁻¹⁶ (Moon) and 3.1 × 10⁻¹⁶ (Sun) at most, and periodic; a timing service at 10⁻¹⁶ would need them.
- Clock stability: a rubidium GPS clock drifts at 10⁻¹³ per day, so terms below 10⁻¹³ (the J₂ term of a satellite on an inclined orbit, the tidal terms) are below the noise and are absorbed in the daily clock-parameter upload instead of being modelled.
- Receiver motion: an aircraft at 250 m/s has its own velocity term (−30 ns/day) and its own Sagnac-like term; both are negligible for positioning, not for carrier-phase timing.
- The budget assumes a circular orbit; GPS orbits have e up to 0.02, and (4) is exactly the correction for that.

## 3. Build the computational model

Open `starter/main.py`. The constants are given; every function is yours.

### Step 1: kinematics

`gamma(v)` with a `ValueError` for |v| ≥ c; `gamma_minus_one(v)` in the cancellation-free form; `time_dilation_rate(v)` = 1/γ − 1 = −(γ − 1)/γ; `length_contraction`; `lorentz_transform(t, x, v)` returning (t′, x′) from (2), `lorentz_inverse` as the same with −v, `spacetime_interval`, `velocity_addition`, `relativistic_doppler(f0, v, theta)` from (3). All must accept NumPy arrays where the docstring says so.

### Step 2: dynamics

`relativistic_energy`, `kinetic_energy` (through `gamma_minus_one`), `momentum`.

### Step 3: gravity and the budget

`circular_orbit_speed(a, mu)`; `gravitational_rate(r, mu)` = −GM/rc²; `geoid_potential_ratio()` = Φ₀/c² with the J₂ and rotation terms; `gps_clock_budget(a, mu, R, omega, J2)` returning a `ClockBudget` whose `net` property sums the four terms and whose `us_per_day()` converts them; `gps_frequency_offset_oracle(a)` = −3GM/2ac² − Φ₀/c²; `factory_clock_frequency`; `crossover_radius` by `scipy.optimize.brentq` on the two-term balance.

### Step 4: the periodic terms and the error

`eccentricity_correction(e, E, a, mu)` from (4); `sagnac_correction(r_sat, r_rec, omega)` from (5); `max_sagnac`; `position_error_without_correction(t, offset)` = c × offset × t; `pass_doppler(f0, v_radial, v_total)` = first-order shift times the constant second-order factor.

Running `main()` prints

```
GPS orbit: a = 26561.75 km, v = 3873.8 m/s, gamma - 1 = 8.349e-11
  velocity         -8.3485e-11     -7.213 us/day
  gravity          +5.2838e-10    +45.652 us/day
  ground_rotation  +1.2034e-12     +0.104 us/day
  ground_j2        +3.7639e-13     +0.033 us/day
  net              +4.4647e-10    +38.575 us/day
oracle -3GM/(2ac^2) - Phi0/c^2 = +4.4647e-10
factory clock frequency: 10229999.99543 Hz  (nominal 10230000 Hz)
crossover radius (simple model): 9567.2 km = 1.5000 R_E, altitude 3189 km
  ISS (h = 400 km)   velocity  -28.27  gravity   +3.55  net  -24.58 us/day
  GPS                velocity   -7.21  gravity  +45.65  net  +38.58 us/day
  GEO                velocity   -4.54  gravity  +50.99  net  +46.58 us/day
eccentricity correction amplitude at e = 0.02: 45.8 ns
maximum Sagnac term: 133.4 ns
pseudorange error after one day without correction: 11.56 km
L1 Doppler at 800 m/s radial: +4204 Hz; second-order (time-dilation) shift: -0.132 Hz
```

and writes `gps_relativity.png`: the three rate terms against altitude with the GPS altitude marked; the accumulated error over 24 h; the eccentricity term over one orbit at e = 0.02; the Sagnac term for an equatorial receiver as the satellite moves from the western to the eastern horizon.

### Verify

Run `python -m pytest curriculum/phase-04/lesson-13-special-relativity-gps`. The tests:

1. **γ**: 1 at rest, 1.25 at 0.6c, 5/3 at 0.8c, monotone, raises at c; `gamma_minus_one` matches β²/2 + 3β⁴/8 to 10⁻¹² at 3873.8 m/s; the airliner loses 30 ns/day.
2. **Lorentz**: the interval is invariant under five random boosts of twenty random events; the inverse round-trips to 10⁻¹²; a light ray stays a light ray; v = 0 is the identity; two boosts compose by velocity addition; 0.5c ⊕ 0.5c = 0.8c, c ⊕ 0.9c = c, 100 ⊕ 200 = 300.
3. **Doppler**: approaching above f₀ above receding, product f₀², transverse f₀/γ, first-order limit to 10⁻¹¹, the pass is antisymmetric to first order with a −0.13 Hz constant.
4. **Dynamics**: E² − (pc)² = (mc²)² at four speeds; K = ½mv² at 7.7 km/s to 10⁻⁹ and larger than ½mv² at 0.9c.
5. **Gravity**: −6.9535 × 10⁻¹⁰ at R_E, zero at infinity, the 1/r law; the geoid ratio is GM/Rc² within 0.3 % and deeper.
6. **Budget**: the three ranges, the net within 0.02 % of 4.4647 × 10⁻¹⁰; equality with the oracle at three radii to 10⁻¹²; the terms scale with GM and the ground terms vanish with ω = J₂ = 0; the crossover at 1.5 R_E to 10⁻⁶, with ISS < 0 < GPS < GEO; the factory frequency 4.567 mHz low.
7. **Periodic terms**: zero at e = 0, 45.8 ns amplitude, zero at E = 0, negative at E = π/2, linear in e, zero mean; Sagnac zero overhead, positive for a western satellite, antisymmetric, linear in ω, 133.4 ns at the horizon.
8. **Error**: linear in t, 11.56 km at one day; the orbit speed is Keplerian with a 43 082 s period; `main` writes the figure.

### Compare

There is no library that computes a clock budget, but there are two independent checks built in. The oracle −3GM/2ac² − Φ₀/c² is a different formula (energy per unit mass rather than a sum of terms) and agrees to 10⁻¹². And the published number: Ashby's review gives the factory offset as −4.4647 × 10⁻¹⁰ and the eccentricity amplitude as 2.29 × 10⁻⁶ e; the MIT 8.033 slides present the same two-term budget. If your code gives 38.4 µs/day you have the two-term model (fine, say so); if it gives 38.6 you have the geoid terms; if it gives anything else, check the sign of the velocity term first, then whether you used the orbit radius or the altitude.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction first.

1. **Altitude sweep.** Tabulate the velocity, gravity and net terms (µs/day) for altitudes 200, 400, 1000, 3189, 10 000, 20 200, 35 786 km. Predict the sign change at 3189 km, −24.6 at the ISS and +46.6 at GEO, and that the gravity term saturates toward GM/Rc² × 86 400 = 60.1 µs/day as a → ∞ while the velocity term goes to zero.
2. **Error growth.** For each term omitted alone (velocity only, gravity only, ground terms only, eccentricity at e = 0.02, Sagnac at the horizon), tabulate the range error after 1 s, 1 min, 1 h, 1 day. Predict that the velocity term alone reaches 1 m after 40 s, the gravity term after 6.3 s, the ground terms after 35 min; that the eccentricity and Sagnac terms do not grow (they are bounded, 13.7 m and 40 m) but do not average out either.
3. **Eccentricity.** For e = 0.005, 0.01, 0.02, tabulate the amplitude in ns and metres, and the maximum *rate* of the correction, 2√(GMa) e n/c² with n the mean motion. Predict that the rate is 1.5 % of the factory offset at e = 0.02, so it cannot be absorbed in a constant.
4. **Sky pass.** For an equatorial receiver and a satellite moving from the western horizon through zenith to the eastern horizon, tabulate the Sagnac term every 15° of position angle. Predict an odd function with ±133 ns at the ends and zero at zenith, and that a receiver at latitude φ sees the equatorial-plane result scaled by roughly cos φ.
5. **Other planets.** Repeat 1 for Mars (GM = 4.2828 × 10¹³ m³/s², R = 3389.5 km): predict the crossover at 1.5 R regardless of the planet (the ratio is pure geometry), and a much smaller net at every altitude because GM/Rc² is nine times smaller.

Raw data format: CSV tables `altitude_km, velocity_us_day, gravity_us_day, net_us_day` (1, 5), `term, t_s, error_m` (2), `e, amplitude_ns, amplitude_m, max_rate` (3), `angle_deg, sagnac_ns` (4); plus `gps_relativity.png`.

What you should find: the sign change at 3189 km to the kilometre; the gravity term at GEO 50.99 µs/day against the 60.1 limit; the error from the velocity term alone 0.025 m/s so 1 m at 40 s; amplitudes 11.4, 22.9, 45.8 ns (3.4, 6.9, 13.7 m) and an eccentricity rate of 6.7 × 10⁻¹² s/s at e = 0.02, 1.5 % of the net; the Sagnac curve ω a R sin(angle)/c² shape; Mars crossover at 5084 km radius (1694 km altitude) and a net of −3.9 µs/day at 450 km.

## 5. Validate and improve

| quantity | hand (§2) | code | published | difference |
|---|---|---|---|---|
| γ − 1 at 3873.8 m/s | 8.349 × 10⁻¹¹ | 8.349 × 10⁻¹¹ | — | 0 |
| velocity term | −7.213 µs/day | −7.213 | ≈ −7 (commonly quoted) | 0 |
| gravity term | +45.652 µs/day | +45.652 | ≈ +45 (commonly quoted) | 0 |
| net, two-term model | +38.44 µs/day | +38.44 | "38 µs/day" | — |
| net, with geoid terms | +4.4647 × 10⁻¹⁰ | +4.4647 × 10⁻¹⁰ | −4.4647 × 10⁻¹⁰ (Ashby, as a frequency offset) | < 0.02 % |
| oracle −3GM/2ac² − Φ₀/c² | +4.4647 × 10⁻¹⁰ | same | — | 10⁻¹² |
| factory frequency | 10 229 999.995 43 Hz | same | 10.229 999 995 43 MHz | 0 |
| eccentricity amplitude, e = 0.02 | 45.8 ns | 45.8 ns | 2.29 × 10⁻⁶ e s (Ashby) | 0 |
| Sagnac maximum | 133.4 ns | 133.4 ns | ω a R sin θ/c² (§2.7) | 0 |
| error after one day | 11.56 km | 11.56 km | "about 10 km" in popular accounts | rounding |

Every number is reproduced. The one honest disagreement is with the popular "10 km per day", which is the 38 µs/day figure rounded down and multiplied by c: 38.6 × 10⁻⁶ × 299 792 458 = 11.57 km. Quote the number you computed and say how.

What to improve: the two-term budget is the physics; the geoid terms are the first thing to add (they are 0.4 % and a 1 m service needs them within 35 minutes); the eccentricity term is next and must be in every receiver; the Sagnac term is a frame-choice artefact that disappears if you compute in the inertial frame, and is 40 m if you do not. What not to improve here: nothing below 10⁻¹³ (tidal terms, inclination-dependent J₂ terms), which the satellite clock's own drift hides and the daily upload absorbs.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of γ, of (2) and the interval invariance, of velocity addition, of (3) and its transverse limit, of the budget table with each term to three figures, of the crossover radius and of (4) from vis-viva; `gps_relativity.png`; the five sweep tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6, 7) stating for a 1 m, a 10 m and a 1 km navigation service which of the five terms must be applied and the time after which each omitted term alone exceeds the budget.

### Adversarial self-checks

- **What if the orbit radius doubles?** The velocity term halves (v² ∝ 1/a), the gravity term grows toward its 60.1 µs/day ceiling, the net grows; at a = 53 124 km the net is +49.4 µs/day. The crossover does not move: it is 1.5 R.
- **What if v → 0?** γ → 1, the velocity term vanishes, and the budget is the pure gravitational redshift between two static clocks; the Doppler shift reduces to f₀; the Sagnac term is unchanged because it comes from the *receiver's* motion.
- **Which assumption fails first in the real world?** The circular orbit: e = 0.02 produces a 46 ns periodic term that no constant offset removes, which is why the eccentricity correction is in the receiver. Next, the clock's own drift (10⁻¹³/day), which hides every term below it.
- **What measurement could falsify the model?** A satellite launched with an uncorrected clock: its broadcast time would gain 38.6 µs/day against ground time, and the pseudoranges would grow at 0.134 m/s. That measurement was made on the first GPS satellite with the offset switchable, and the drift matched.
- **Does the result survive a change in input?** Using the altitude (20 200 km) instead of the radius in the velocity term gives −9.5 µs/day instead of −7.2; using the polar radius instead of the equatorial changes the gravity term by 0.2 µs/day. The oracle catches the first; nothing catches the second, which is why the geoid, not a radius, defines ground time.
- **Which decision changes if the receiver's required accuracy loosens from 1 m to 1 km?** The factory offset still matters (1 km is reached in 2 h); the eccentricity and Sagnac terms (14 m, 40 m) no longer do; the geoid terms (0.14 µs/day, 1 km after 24 days) matter only if the clock parameters are not refreshed daily.

## Check yourself

1. State the two postulates and write the Lorentz factor. What is γ at 0.6c, and what does it do to a clock?
2. A moving clock's rate is 1/γ − 1. To second order in v/c, what is that, and what is its value for a GPS satellite in µs/day?
3. Write the weak-field gravitational rate between two static clocks at potentials Φ₁ and Φ₂. Which runs faster?
4. Give the two main GPS terms with signs and magnitudes, and the net. Why do they not cancel at GPS altitude, and at what altitude would they?
5. Why is the eccentricity correction not absorbed in the factory offset, and what is its amplitude for e = 0.02?
6. What is the Sagnac term physically, when is it zero, and what is its maximum for GPS?
7. Write the relativistic Doppler formula and say which part of it survives when the pass is symmetric.
8. Write E² − (pc)² = (mc²)². What is the kinetic energy of a 1 kg mass at 7.7 km/s, relative to ½mv²?
