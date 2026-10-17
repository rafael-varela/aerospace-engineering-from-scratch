# P04.09  Magnetism & Induction

> A current makes a field, a field pushes on a current, and a changing field makes a current. Three sentences, one constant μ₀, and a small satellite can turn itself around by pushing on the Earth.

| | |
|---|---|
| **Learning outcomes** | 1. Compute the Lorentz force on a moving charge and on a current element, and derive the gyroradius and cyclotron frequency of a charged particle in a uniform field. 2. Derive the field on the axis of a circular loop and of a finite solenoid from the Biot–Savart law, implement the law as a sum over wire segments, and show that the sum converges at second order to the closed forms and agrees with a `scipy.integrate.quad` oracle off the axis. 3. Model the Earth's field as a centred dipole, evaluate it along a circular orbit, and compute the torque m × B a magnetorquer coil can produce, including the direction about which it can produce none. 4. Apply Faraday's and Lenz's laws to a rotating coil and a motor, derive the RL time constant, and simulate a one-axis B-dot detumble whose unsaturated phase decays with the closed-form time constant I/(kB²). |
| **Prerequisites** | `P04.08` (current, resistance, Kirchhoff, the RC time constant), `P01.07` (vectors, the cross product and its right-hand rule), `P03.04` (quadrature: the Biot–Savart integral is one), `P03.05` and `P02.04` (`solve_ivp` on a first-order system), `P00.14` (reporting a number with its uncertainty). Skill check: you can compute a × b by hand and say which way it points. |
| **Engineering question** | A 3U CubeSat (I = 0.04 kg·m² about the spin axis) separates from its deployer tumbling at 10°/s and must be detumbled by an air-core coil of 200 turns on an 8 cm face, drawing 0.12 A, in a 500 km orbit where the Earth's field is 25–49 μT: what torque can the coil produce, about which axis can it produce none, how long does the detumble take, and how fast can the coil current switch? The answers: m = 0.154 A·m², up to 5.7 μN·m, never about the local field direction, about 71 minutes (4250 s) to reach 0.1°/s with a 30 μT field, and 0.18 ms for the current to settle. |
| **Assumptions and limits** | Magnetostatics (steady or slowly varying currents, no radiation); thin wires; vacuum permeability everywhere (no iron: a ferromagnetic torquer rod multiplies the moment by hundreds and is not modelled); the Earth's field as an untilted, centred dipole, correct to 10–20 % in magnitude at low altitude and wrong in detail (the real dipole is tilted 11°, offset from the centre, and the field has non-dipole parts); the detumble about one axis perpendicular to a fixed field, which is the planar case of a three-axis problem whose field rotates once per orbit. Eddy currents, hysteresis and the inductance of a flat coil are estimated, not computed. |
| **Analytical work** | The Lorentz force; the gyroradius mv/(qB) and cyclotron frequency qB/m; the Biot–Savart integral on the axis of a loop and its dipole limit; the finite solenoid by integrating loop fields; μ₀nI for the infinite solenoid; the dipole field 3(m·r̂)r̂ − m and its 2:1 pole-to-equator ratio; the torque on a rectangular loop from the force on its sides; Faraday's law for a rotating coil; the RL equation and its time constant; the B-dot law and its time constant I/(kB²). |
| **Computational deliverable** | `starter/main.py`: `lorentz_force`, `gyroradius`, `cyclotron_frequency`, `wire_field`, `loop_field`, `loop_field_axis`, `loop_field_quad`, `straight_wire_field`, `solenoid_field_axis`, `solenoid_inductance`, `dipole_field`, `earth_dipole_field`, `field_along_orbit`, `magnetic_moment`, `torque_on_dipole`, `coil_resistance`, `faraday_emf`, `generator_emf`, `motor_back_emf`, `rl_current`, `rl_simulate`, `detumble_time_constant`, `saturation_rate`, `bdot_detumble`, `momentum_dump_time` (`loop_points`, `orbit_positions`, the `Torquer` dataclass and `main` are given). Tests prove: F ⊥ v and F ⊥ B with magnitude qvB sin θ; r·ω_c = v; the segment sum on the loop axis converges at order 2.0 ± 0.2 to μ₀IR²/(2(R² + z²)^{3/2}) and matches a `quad` oracle off the axis to 10⁻⁴; B_z is even and B_ρ odd in z; a long polyline gives μ₀I/(2πd) to 10⁻³; the solenoid gives μ₀nI at the centre and half at the end and equals the `quad` integral of loop fields to 10⁻⁸; the dipole has the 2:1 ratio, 1/r³ and zero divergence and curl; the orbit field is constant on an equatorial orbit and has a 2:1 range on a polar one; m × B has no component along B or m; the emf equals −dΦ/dt and peaks at NABω; the RL closed form matches `solve_ivp` to 10⁻⁹ A with the 63.2 % point at τ; the detumble decays as exp(−t/T_d) when unsaturated and linearly when saturated, monotonically throughout. |
| **Practical activity** | **computational experiment**: a convergence study of the Biot–Savart sum, a field map of a loop, the dipole field along orbits of three inclinations, the torque available around an orbit, the RL rise of the coil current against `solve_ivp`, and a detumble parameter sweep over gain, coil limit and field strength. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `magnetorquer.png`, the derivations of the loop axis field, the solenoid, the loop torque and the RL and B-dot time constants, the four experiment tables, and a half-page note sizing the torquer for a stated detumble requirement with the evidence. |
| **Assessment** | `assessments/phase-04/grading-rubric.md`; competencies 1–4. Pass: the on-axis loop field and the torque on a loop are derived from the Lorentz force; all public tests pass; the convergence order is reported and the off-axis comparison with `quad` explained; the note states the torque, the uncontrollable axis and the detumble time with their basis. |
| **References** | `@openstax-university-physics-vol2` Ch. 11 "Magnetic Forces and Fields": 11.2 Magnetic Fields and Lines, 11.3 Motion of a Charged Particle in a Magnetic Field, 11.4 Magnetic Force on a Current-Carrying Conductor, 11.5 Force and Torque on a Current Loop (required-free): the Lorentz force, the gyroradius, and τ = m × B; 11.6 The Hall Effect (supplementary-free): how a magnetometer chip measures B. Ch. 12 "Sources of Magnetic Fields": 12.2 Magnetic Field Due to a Thin Straight Wire, 12.5 Ampère's Law, 12.6 Solenoids and Toroids (required-free): the Biot–Savart law and its three classic results; 12.7 Magnetism in Matter (supplementary-free): why an iron core multiplies the moment. Ch. 13 "Electromagnetic Induction": 13.1 Faraday's Law, 13.2 Lenz's Law, 13.3 Motional Emf, 13.6 Electric Generators and Back Emf (required-free); 13.5 Eddy Currents (supplementary-free). Ch. 14 "Inductance": 14.2 Self-Inductance and Inductors, 14.4 RL Circuits (required-free). |
| **Extensions** | ★ Replace the fixed field in `bdot_detumble` by the orbit field of `field_along_orbit` rotating in the body frame and integrate the three-axis Euler equations with `solve_ivp`; show the detumble still converges and find the residual rate about the field direction. ★ Add a tilted dipole (11°) and compare the field along a sun-synchronous orbit over one day with the untilted model. ★ Compute the mutual inductance of two coaxial loops by integrating `loop_field_quad` over one loop's area and check Neumann's reciprocity M₁₂ = M₂₁. |
| **Time** | unmeasured |

## 1. The engineering problem

Every CubeSat leaves its deployer tumbling. Before the solar panels can point at the Sun, the antenna at the ground and the camera at anything, the tumble must be removed, and most small satellites remove it with the cheapest actuator in spaceflight: a coil of wire. A current I in N turns of area A has a magnetic moment m = NIA, and in the Earth's field B it feels a torque τ = m × B. The field in low orbit is a few tens of microtesla, so a 0.15 A·m² coil gives a few micronewton-metres, a torque that would be laughable on the ground and is the standard tool in orbit, where there is nothing else to push against and nothing to stop a small torque from acting for an hour. The same physics, in reverse, is the magnetometer that measures the field (the Hall effect), the brushless motor in a reaction wheel (force on a current, back-emf from the spinning magnet) and the eddy-current damper on a deployable boom.

This lesson builds the torquer from the Lorentz force up. Part 2 derives the field of a current from the Biot–Savart law, first for the cases with closed forms (the loop axis, the straight wire, the solenoid), then the Earth's field as a dipole, then the torque on the coil, then Faraday's law and the RL time constant that decides how fast the coil can be switched. Part 3 implements the law as a sum over wire segments and checks that sum against the closed forms and against an adaptive quadrature. Part 4 runs the design: the field along three orbits, the torque available around the orbit, the one-axis detumble with a saturating coil.

Success criteria, written first: the segment sum converges to the on-axis closed form at second order and agrees with `quad` off the axis to 10⁻⁴; the solenoid formula reproduces μ₀nI at the centre of a long coil; the dipole model gives 31 μT at the surface equator and 62 μT at the pole and 25–49 μT around a 500 km orbit; the torque never has a component along B; the RL rise reaches 63.2 % at t = L/R; the detumble follows exp(−t/T_d) with T_d = I/(kB²) once the coil is unsaturated.

## 2. Understand and derive

### 2.1 The Lorentz force

A charge q moving at velocity v in electric field E and magnetic field B feels

**F = q(E + v × B).**   (1)

The magnetic term is perpendicular to v, so it does no work: a magnetic field changes the direction of a charge's motion, never its speed. In a uniform B with v ⊥ B the motion is a circle: qvB = mv²/r gives the **gyroradius** r = mv/(qB) and the **cyclotron frequency** ω_c = qB/m, independent of speed (OpenStax 11.3). A 1 MeV-class proton at 10⁶ m/s in the 30 μT of low orbit feels 4.8 × 10⁻¹⁸ N, circles with r = 348 m at 457 Hz; this is why the radiation belts are belts, charged particles spiralling along field lines, and why a magnetic field can trap them.

Summed over the charges in a wire, (1) becomes the force on a current element, dF = I dl × B, and on a straight length L in a uniform field, F = IL × B (11.4): 1 A in 1 m across 30 μT is 30 μN. Dimensionally, [B] = N/(A·m) = T; the Earth's surface field is 3–6 × 10⁻⁵ T, a fridge magnet 10⁻² T, an MRI machine 1–3 T.

### 2.2 The Biot–Savart law and its three closed forms

A current element I dl at position r′ produces at r the field

**dB = (μ₀/4π) I dl × (r − r′)/|r − r′|³,  μ₀ = 4π × 10⁻⁷ T·m/A,**   (2)

and the field of a circuit is the integral of (2) around it (12.2–12.6 build the same results). Three integrals are worth doing by hand because every numerical Biot–Savart code must reproduce them.

**The axis of a circular loop.** Radius R in the plane z = 0, current I anticlockwise seen from +z, field point (0, 0, z). Every element has |r − r′| = √(R² + z²) and dl ⊥ (r − r′); the transverse components of dl × (r − r′) cancel around the loop and the axial component is dl √(R² + z²) × R/√(R² + z²) = R dl. Integrating dl over 2πR:

**B_z(z) = μ₀IR² / (2(R² + z²)^{3/2}).**   (3)

At the centre B = μ₀I/(2R) (a 4 cm loop with 0.12 A: 1.885 μT); far away, B ≈ μ₀IR²/(2|z|³) = μ₀m/(2π|z|³) with **m = IπR²**: the field on the axis of a dipole. At z = 10R the exact (3) is 1.857 nT and the dipole formula 1.885 nT, 1.5 % high, the leading correction being −(3/2)(R/z)².

**The straight wire.** Integrating (2) along an infinite wire at distance d gives B = μ₀I/(2πd), azimuthal by the right-hand rule, the same result Ampère's law gives in one line (12.5): 1 A at 5 cm is 4.0 μT.

**The solenoid.** A coil of n turns per metre, radius R, length L, is a stack of loops; integrating (3) over the windings from −L/2 to L/2 with the field point at z,

**B_z(z) = (μ₀nI/2)[ (z + L/2)/√((z + L/2)² + R²) − (z − L/2)/√((z − L/2)² + R²) ].**   (4)

For L ≫ R the bracket is 2 at the centre, B = μ₀nI (Ampère's law again, 12.6), and 1 at either end: the end field is half the centre field, because half the windings are missing. A 10 cm coil of 2000 turns/m and 5 mm radius at 1 A: 2.5008 mT at the centre against μ₀nI = 2.5133 mT, 1.2551 mT at the end. The flux linkage NΦ = N·μ₀nI·πR² per ampere is the **self-inductance** L_ind = μ₀N²A/L (14.2), 39 μH for 200 turns of that coil; it is what resists a change of current.

### 2.3 The Earth as a dipole

Far from any current the field of a loop is that of a point dipole,

**B(r) = (μ₀/4π)[3(m·r̂)r̂ − m]/r³,**   (5)

and the Earth's main field is, to 10–20 %, a dipole of moment 8.0 × 10²² A·m² at the centre, pointing roughly toward the geographic south pole (so that field lines leave near the south pole and enter near the north; the "north" end of a compass points along B). The lesson's model takes it untilted and centred, with the surface equatorial field B₀ = 31 μT as the input (5 then gives m = B₀R_E³/(μ₀/4π)). At the equator r̂ ⊥ m and |B| = (μ₀/4π)m/r³ = B₀; at the pole r̂ ∥ m and |B| = 2B₀ = 62 μT: the **2:1 pole-to-equator ratio** is the first thing a magnetometer on a polar orbit sees. The field falls as 1/r³: at 500 km the equatorial value is 31 × (6371/6871)³ = 24.7 μT. In general |B| = B₀(R_E/r)³√(1 + 3 sin²λ) with λ the magnetic latitude.

Along a circular orbit of inclination i the magnetic latitude reaches ±i, so an equatorial orbit sees a constant 24.7 μT, the ISS inclination (51.6°) sees 24.7–41.7 μT (mean 33.7), and a sun-synchronous orbit (98°) sees 24.7–49.1 μT (mean 37.9). (5) is divergence-free and curl-free away from the origin, as a vacuum field must be; the tests check both by finite differences.

Where the model fails: the real dipole is tilted 11° and offset by about 500 km from the centre, the South Atlantic Anomaly weakens the field by a third over Brazil, and the non-dipole part is 10 % at the surface. For a torque *budget* the model is adequate; for attitude *determination* from a magnetometer it is not, and spacecraft carry a full spherical-harmonic model.

### 2.4 Torque on a coil

Take a rectangular loop of sides a and b, current I, in a uniform B. The forces on opposite sides are equal and opposite (IL × B with L reversed), so the net force is zero; but the two forces on the sides of length a act along lines separated by b sin θ, θ the angle between the loop normal and B, and form a couple of torque IaB × b sin θ = IAB sin θ. In vector form, with the **magnetic moment** m = NIA n̂,

**τ = m × B.**   (6)

(11.5; Problem 2 does the derivation for any orientation.) Two consequences run the whole lesson. The torque is at most mB, reached when m ⊥ B; and **τ·B = 0 always**: a coil, or any set of coils, cannot produce torque about the local field direction. A magnetorquer system is therefore under-actuated at every instant and relies on B rotating in the body frame (once per orbit, as the satellite moves through the dipole field) to reach the missing axis over time.

The lesson's torquer is 200 turns on an 8 cm square (A = 6.4 × 10⁻³ m²) of 0.25 mm copper wire: m = 200 × 0.12 × 0.0064 = 0.154 A·m², resistance ρ(4N·side)/(πd²/4) = 21.9 Ω, power I²R = 0.32 W, and torque up to 0.154 × 41.7 μT = 6.4 μN·m at the high-latitude end of the ISS orbit. With the coil axis along the dipole axis the torque around the orbit ranges from 0 (over the equator, where B ∥ z) to 5.69 μN·m.

### 2.5 Faraday, Lenz and the RL time constant

A changing flux Φ = ∫B·dA through a circuit drives an emf

**ℰ = −dΦ/dt,**   (7)

with the sign (Lenz's law, 13.2) saying that the induced current's own field opposes the change. A coil of N turns and area A rotating at ω in a uniform B has Φ = NAB cos ωt and ℰ = NABω sin ωt (13.6): the generator. The torquer coil spun at 50 Hz in 0.5 T would give 201 V peak. The same equation, in a motor, is the **back-emf** e = k_e ω that rises with speed until it nearly cancels the supply voltage: a brushless reaction-wheel motor at no load spins at V/k_e, and its current, hence torque, falls as it speeds up. A conductor moving at v through B sees the motional emf vBL (13.3); a 650 m tether at orbital speed across 30 μT sees 148 V (Problem 3).

The coil's own flux resists a change of its own current: the self-inductance gives e = −L di/dt, and a step V applied to a coil of resistance R obeys L di/dt = V − Ri,

**i(t) = (V/R)(1 − e^{−t/τ}),  τ = L/R,**   (8)

(14.4): the current reaches 63.2 % of V/R at t = τ and its initial slope is V/L. The torquer's inductance is not a solenoid's (it is a flat coil), and the solenoid formula with the side as the effective length gives an order of magnitude, 4 mH; with 21.9 Ω that is τ ≈ 0.18 ms, so the coil current can be switched thousands of times per second and the controller, not the coil, sets the bandwidth. When the current is switched *off*, (8) runs in reverse and the energy ½Li² must go somewhere: without a flyback diode it appears as a voltage spike L di/dt across the switch.

### 2.6 One-axis detumble: the B-dot law

The simplest detumble controller commands a moment proportional to the rate of change of the measured field in the body frame, m = −k dB_body/dt. For a body spinning at ω in a fixed field, dB_body/dt = −ω × B, so m = k ω × B and the torque is

**τ = m × B = k(ω × B) × B = k[(ω·B)B − B²ω].**   (9)

With ω ⊥ B this is −kB²ω: a torque against the spin, proportional to it. Then I dω/dt = −kB²ω and

**ω(t) = ω₀ e^{−t/T_d},  T_d = I/(kB²).**   (10)

The component of ω along B is untouched, as 2.4 said it must be. A real coil has a maximum moment m_max, so while kωB > m_max the moment sits at m_max, the torque is m_max B and the rate falls linearly at m_max B/I until ω_s = m_max/(kB), after which (10) takes over. For the CubeSat (I = 0.04 kg·m², ω₀ = 10°/s, B = 30 μT, k = 5 × 10⁴ A·m²·s/T, m_max = 0.154 A·m²): the saturated phase lasts 626 s and ends at 5.87°/s, then T_d = 889 s; the rate is 3.86°/s at 1000 s, 1.25 at 2000, 0.13 at 4000 and 0.014 at 6000. The absolute lower bound for any controller limited to τ_max = m_max B is Iω₀/τ_max = 1515 s; B-dot takes longer because it throttles back once the moment falls below the limit. In three dimensions the field rotates once per orbit, the B-dot torque reaches each axis in turn, and the detumble takes a few orbits.

### 2.7 Where the models stop

- Biot–Savart as used here assumes thin wires and steady current; at the coil's switching frequency the quasi-static approximation holds easily (the wavelength at 5 kHz is 60 km), and it fails only for antennas, the subject of `P04.11`.
- No magnetic materials: a permalloy rod inside a torquer coil raises m by a factor of hundreds (12.7), with hysteresis and saturation that this lesson does not model.
- The dipole Earth: 10–20 % on magnitude, worse on direction; the detumble time scales as 1/B², so a 20 % field error is a 40 % time error.
- The one-axis detumble is the planar case; the three-axis one with a rotating field is an extension.
- The flat coil's inductance is an order of magnitude estimate; measure it.

## 3. Build the computational model

Open `starter/main.py`. `loop_points`, `orbit_positions`, the `Torquer` dataclass and `main` are given.

### Step 1: the Lorentz force and its circle

`lorentz_force(q, v, E, B)` is (1) with `np.cross`. `gyroradius` and `cyclotron_frequency` follow from qvB = mv²/r.

### Step 2: Biot–Savart as a sum over segments

`wire_field(points, I, r, closed)` discretises (2) with the midpoint rule: for each segment, dl is the vector between consecutive points and r′ its midpoint; sum (μ₀I/4π) dl × (r − r′)/|r − r′|³. Vectorise over segments and over field points (shape (m, n, 3)):

```python
dl = P[1:] - P[:-1]; mid = 0.5 * (P[1:] + P[:-1])
R = r[:, None, :] - mid[None, :, :]
B = MU0 * I / (4 * np.pi) * np.sum(np.cross(dl, R) / np.linalg.norm(R, axis=-1)[..., None] ** 3, axis=1)
```

`loop_field` applies it to the n-gon from `loop_points`; `loop_field_axis` is (3); `loop_field_quad(rho, z, R, I)` integrates the exact integrand over the loop angle with `scipy.integrate.quad` (write the integrand yourself: with the field point at (ρ, 0, z) and the source at (R cos φ′, R sin φ′, 0), B_z ∝ ∫R(R − ρ cos φ′)/d³ dφ′ and B_ρ ∝ ∫Rz cos φ′/d³ dφ′). `straight_wire_field` is μ₀I/(2πd); `solenoid_field_axis` is (4); `solenoid_inductance` is μ₀N²A/L.

### Step 3: the Earth's field

`dipole_field(r, m)` is (5), accepting one position or a stack. `earth_dipole_field` calls it with m = −M ẑ. `field_along_orbit(altitude, inclination)` evaluates it at `orbit_positions`.

### Step 4: the coil

`magnetic_moment`, `torque_on_dipole` (6), `coil_resistance`. The `Torquer` dataclass collects them.

### Step 5: induction

`faraday_emf(t, flux)` is −`np.gradient`; `generator_emf` is NABω sin ωt; `motor_back_emf` is k_e ω; `rl_current` is (8) with an initial current; `rl_simulate` integrates L di/dt = V(t) − Ri with `solve_ivp` for any source.

### Step 6: the detumble

`detumble_time_constant` and `saturation_rate` are the two constants of 2.6; `bdot_detumble` integrates I dω/dt = −min(k|ω|B, m_max) B sign(ω) with `solve_ivp` (pass `max_step` so the kink at saturation is resolved); `momentum_dump_time` is Iω₀/τ_max.

Running `main()` prints

```
proton at 1e+06 m/s across 30 uT: force 4.81e-18 N, gyroradius 348 m, cyclotron frequency 457 Hz

loop R = 4 cm, I = 0.12 A: centre field 1.885 uT, far field at z = 10 R: 1.857 nT (dipole: 1.885 nT)
segments   |B_z(0,0,R/2) - closed form| / closed form
       8   8.54e-02
      16   2.08e-02
      32   5.15e-03
      64   1.29e-03
     128   3.21e-04
observed orders: [2.04 2.01 2.   2.  ]
off axis (rho, z) = (R/2, 0.3 R): quad (B_rho, B_z) = (0.4916, 1.8108) uT, 360-gon (0.4917, 1.8109) uT, B_phi = 1.2e-23 T
wire 40 m long in 5 mm segments, 1 A, at 5 cm: 4.0000 uT against mu0 I / 2 pi d = 4.0000 uT
solenoid 2000 turns/m, L = 10 cm, R = 5 mm, 1 A: centre 2.5008 mT, end 1.2551 mT, mu0 n I = 2.5133 mT

Earth dipole: m = 8.016e+22 A m^2; surface equator 31.0 uT, pole 62.0 uT
  500 km, i =   0.0 deg: |B| min 24.71, max 24.71, mean 24.71 uT
  500 km, i =  51.6 deg: |B| min 24.71, max 41.67, mean 33.71 uT
  500 km, i =  98.0 deg: |B| min 24.71, max 49.07, mean 37.86 uT

torquer: N = 200, side 8 cm, wire 0.25 mm, I = 0.12 A: m = 0.154 A m^2, R = 21.9 ohm, P = 0.32 W, L ~ 4.0 mH, tau = L/R ~ 0.18 ms
torque with m along z, 500 km, i = 51.6: max 5.69 uN m, min 0.00 uN m; component along B: 4.0e-22 N m
```

followed by the induction lines (peak emf 201.062 V against NABω = 201.062 V; RL at t = τ: 63.2 %, `solve_ivp` against the closed form 2.3 × 10⁻¹² A) and the detumble table, and writes `magnetorquer.png` (the loop's field lines; |B| around three orbits; the RL rise; the detumble with its exponential tail).

### Verify

Run `python -m pytest curriculum/phase-04/lesson-09-magnetism-induction`. The tests:

1. **Lorentz**: F ⊥ v, F ⊥ B, |F| = |q|vB sin θ on random vectors; zero for v ∥ B; the E term adds; r ω_c = v and the scalings of r with v and B.
2. **Loop axis**: μ₀I/(2R) at the centre; the dipole limit at z = 40R to 0.2 %; evenness in z.
3. **Segment sum**: purely axial on the axis; relative error against (3) falling with observed order 1.8–2.2 from 16 to 128 segments; below 10⁻³ at 128.
4. **Off axis**: B_ρ and B_z of a 720-gon against an independent `quad` oracle to 10⁻⁴, B_φ = 0; `loop_field_quad` against the same oracle to 10⁻⁹.
5. **Symmetry**: B_z even and B_ρ odd in z; invariance under rotation about the axis; linearity in I and sign reversal.
6. **Wire**: a 40 m polyline in 5 mm segments at 5 cm against μ₀I/(2πd) to 10⁻³, azimuthal direction by the right-hand rule.
7. **Solenoid**: μ₀nI at the centre and half at the end of a 1 m coil to 10⁻⁴; a 2 cm coil against the `quad` integral of (3) to 10⁻⁸; inductance ∝ N² and A/L.
8. **Dipole**: B₀ at the equator, 2B₀ at the pole, directions, 1/8 at 2R_E, 20–50 μT around a 500 km orbit; divergence and curl zero by central differences at 1 m; batch shape.
9. **Orbit field**: constant on the equatorial orbit with the (R_E/r)³ value; max/min = 2 on the polar orbit; the orbit radius.
10. **Torque**: τ ⊥ B and τ ⊥ m on random vectors, |τ| ≤ mB, zero for m ∥ B; NIA; coil resistance ∝ N and 1/d²; the `Torquer` power and moment in the CubeSat range.
11. **Induction**: −dΦ/dt against NABω sin ωt; a linearly rising flux gives a constant negative emf; k_e ω.
12. **RL**: i(τ) = 0.632 V/R, initial slope V/L, closed form against `solve_ivp` to 10⁻⁹ A, `rl_simulate` against the closed form, decay from an initial current.
13. **Detumble**: with an unsaturated coil ω = ω₀ exp(−t/T_d) to 10⁻⁶; with a saturated coil a linear fall at m_max B/I, monotone throughout, ending below 0.1 ω₀.
14. **main** writes the figure.

### Compare

`scipy.integrate.quad` is the oracle for the off-axis field: it integrates the exact Biot–Savart integrand adaptively to 10⁻¹⁰ and the 360-gon sum agrees with it to 10⁻⁴ (0.4917 against 0.4916 μT in B_ρ, 1.8109 against 1.8108 μT in B_z). The sum is what you use on a real coil, which is not a circle but a polyline of whatever shape the winding took; `quad` is what you use when a one-dimensional integrand exists and you want all the digits. `solve_ivp` (RK45 at rtol 10⁻¹⁰) agrees with the RL closed form to 2 × 10⁻¹² A, and is the tool for the detumble, whose saturating right-hand side has no closed form across the kink; the lesson passes `max_step` so the step never jumps over it, and checks the two phases separately against (10) and the linear slope.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction first.

1. **Convergence of the segment sum.** Relative error of `loop_field` against (3) at z = R/2 and at z = 0.05R for 8–512 segments. Predict order 2 at both, and predict whether moving the point toward the loop's plane matters (the midpoint rule's error scales with (segment length / distance to the wire)², and on the axis that distance is √(R² + z²), never less than R).
2. **Field map.** `loop_field` on a 121 × 121 grid of the x–z plane, |B| contours and field lines; compare |B| on the axis at z = R, 2R, 3R with (3) and with the dipole formula. Predict how far out the dipole formula is within 20 %.
3. **Field along orbits.** |B| and its components around orbits at 400, 500 and 800 km for i = 0°, 51.6°, 98°; mean and range. Predict the mean to scale as (R_E/r)³ and the range to grow with inclination up to 2:1 at 90°.
4. **Torque around the orbit.** For the 51.6° orbit, the torque from `torque_on_dipole` with m along each body axis in turn (the body held fixed in the dipole frame), and the component of each along B. Predict zero along B always, and that at least one of the three axes has a torque minimum of exactly zero somewhere on the orbit.
5. **Detumble sweep.** `bdot_detumble` at k = 2 × 10⁴, 5 × 10⁴, 10⁵; m_max = 0.05, 0.15, 0.5 A·m²; B = 25, 35, 45 μT; record the time to 0.1°/s. Predict the time to scale as 1/B² when unsaturated and 1/B when saturated, and to be bounded below by Iω₀/(m_max B).

Raw data format: CSV tables `n_seg, z_over_R, rel_error` (1); `z_over_R, B_exact, B_sum, B_dipole` (2); `altitude_km, inclination_deg, B_min_uT, B_max_uT, B_mean_uT` (3); `axis, u_deg, tau_uNm, tau_dot_B` (4); `k, m_max, B_uT, t_to_0p1_degps` (5); plus `magnetorquer.png`.

What you should find: orders 2.00 at both heights, with the error at z = 0.05R only 1.6 times the z = R/2 value at equal n (5.0 × 10⁻⁴ against 3.2 × 10⁻⁴ at 128 segments; 3.1 × 10⁻⁵ against 2.0 × 10⁻⁵ at 512), because the distance to the wire barely changes; the dipole formula 183 % high at z = R, 40 % at 2R, 17 % at 3R and 1.5 % at 10R; means 24.7, 33.7, 37.9 μT at 500 km with the ranges in the table above; the torque with m along z vanishing over the equator; detumble times from 934 s (k = 10⁵, m_max = 0.5, B = 45 μT) to 15 340 s, 4.3 hours (k = 2 × 10⁴, m_max = 0.05, B = 25 μT), the lower bound Iω₀/(m_max B) never violated; at k = 2 × 10⁴ raising m_max from 0.15 to 0.5 changes nothing (14 737 s at 25 μT both ways), because the coil is no longer saturated and only the gain matters.

## 5. Validate and improve

| quantity | theory (§2) | code | difference | cause |
|---|---|---|---|---|
| loop centre field, R = 4 cm, 0.12 A | μ₀I/2R = 1.885 μT | 1.885 μT | 0 | — |
| far field at z = 10R | dipole 1.885 nT | 1.857 nT | 1.5 % | −(3/2)(R/z)² correction |
| segment-sum order | 2 | 2.04, 2.01, 2.00, 2.00 | 0.04 | pre-asymptotic at 8 segments |
| off-axis B_ρ, B_z | `quad` 0.4916, 1.8108 μT | 0.4917, 1.8109 μT | 2 × 10⁻⁴ | 360 segments |
| straight wire, 1 A, 5 cm | 4.0000 μT | 4.0000 μT | < 10⁻⁴ | 5 mm segments |
| solenoid centre | μ₀nI = 2.5133 mT | 2.5008 mT | 0.5 % | finite length: L/R = 20 |
| solenoid end | half the centre | 1.2551 mT | 0.4 % of half | same |
| dipole pole/equator | 2 | 2.000 | 0 | — |
| 500 km equatorial field | 31 × (6371/6871)³ = 24.71 μT | 24.71 μT | 0 | — |
| polar-orbit range | 2:1 | 49.07/24.71 = 1.99 at 98° | — | 98° is not 90° |
| generator peak emf | NABω = 201.06 V | 201.06 V | 0 | — |
| RL at t = τ | 63.2 % | 63.2 % | 0 | — |
| detumble, unsaturated | exp(−t/T_d), T_d = 889 s | matches to 10⁻⁶ | — | — |
| detumble, saturated | slope m_max B/I = 0.396°/s per min | matches to 10⁻⁵ | — | — |

Every row with a closed form is met. The three that depart (far field, solenoid, polar range) depart by exactly the amount the next term of the theory predicts, which is the check that matters: a discrepancy you can compute is a validation, one you cannot is a bug or a wrong model. The engineering answer is that the torquer gives **0.154 A·m², up to 5.7 μN·m with the coil axis along the dipole axis on the ISS orbit, nothing about B, and a detumble from 10°/s to 0.1°/s in about 71 minutes (4250 s)** at 30 μT, for 0.32 W.

What to improve: the detumble with the real rotating field (Extension 1), which makes all three axes reachable and lengthens the time to a few orbits; the flat coil's inductance by numerical integration of the flux of `loop_field` through its own area, or better, a measurement; and the dipole by a tilted, offset one. What not to improve: the Biot–Savart sum beyond 128 segments, since the coil's real winding geometry is uncertain at the per-cent level, and the detumble to better than 10 %, since the field is uncertain by that much.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (3), (4), (6), (8) and (10); `magnetorquer.png`; the five experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) sizing a torquer coil for a stated satellite, orbit and detumble requirement, with the torque, the uncontrollable axis, the power and the time, and the evidence for each.

### Adversarial self-checks

- **What if the current doubles?** The moment, the torque and the saturated detumble slope double; the power quadruples (0.32 → 1.26 W, the real constraint on a CubeSat bus); the unsaturated time constant I/(kB²) does not change at all, because it depends on the gain, not the coil, until the coil saturates.
- **What if the orbit altitude goes to 2000 km?** B falls as 1/r³, 24.7 → 13.7 μT at the equator, and the detumble time grows as 1/B² (unsaturated) to 1/B (saturated): roughly three times longer. At geostationary altitude (r = 6.6 R_E) the field is 0.1 μT and magnetorquers are useless; that is why GEO satellites carry thrusters and wheels.
- **Which assumption fails first in the real world?** The fixed field of the one-axis model: in orbit B rotates in the body frame once per orbit (and faster, if the body tumbles), so the B-dot controller sees a changing field even with no tumble and must distinguish the two; next the dipole magnitude (10–20 %); then the air core, if a rod is used.
- **What measurement could falsify the model?** A three-axis magnetometer log over one orbit: the dipole predicts a 2:1 magnitude range on a polar orbit and 24.7 μT at the equator at 500 km; a 30 % departure is the tilt and the anomaly, a 300 % departure is a wrong model or a magnetised satellite. On the bench, the coil's field on its axis at a known distance against (3), with the ambient field subtracted.
- **Does the result survive a change of discretisation or solver?** The field at 128 segments is within 3.2 × 10⁻⁴ of the closed form and at 512 within 2 × 10⁻⁵; the agreement between the segment sum and `quad` is the discretisation check and the agreement between `solve_ivp` and (10) the solver check.
- **Which decision changes if the field uncertainty doubles?** From ±20 % to ±40 % on B the unsaturated detumble time (∝ 1/B²) spans a factor of about 5 instead of 2; the gain and the coil would be sized for the weak-field case, or the schedule would allow a day. The uncontrollable axis is not uncertain at all: τ·B = 0 is exact.

## Check yourself

1. Write the Lorentz force and explain in one sentence why a magnetic field never changes a particle's speed.
2. State the field on the axis of a circular loop at its centre and far away; what quantity replaces IR² in the far-field form?
3. Why is the field at the end of a long solenoid half the field at its centre?
4. The Earth's field at the magnetic equator is 31 μT at the surface. What is it at the pole, and at 500 km over the equator?
5. A magnetorquer coil has m = 0.15 A·m² and the local field is 40 μT. What is the largest torque it can produce, and about which axis can it produce none?
6. A coil's flux is increasing. Which way does the induced current flow, and what law says so?
7. What is the time constant of a coil of 4 mH and 22 Ω, and what fraction of the final current flows after one time constant?
8. In the B-dot law, how does the detumble time scale with the field strength when the coil is unsaturated, and when it is saturated?
