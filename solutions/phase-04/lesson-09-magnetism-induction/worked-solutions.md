# P04.09  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` (and the probe scripts that call it) and given to the figures shown.

## 1. A charged particle in the Earth's field (s = 15: v = 2.5 × 10⁶ m/s)

**Assumptions.** Uniform field, v ⊥ B, non-relativistic (v/c < 0.01).

(a) F = evB = 1.602 × 10⁻¹⁹ × 2.5 × 10⁶ × 3 × 10⁻⁵ = **1.20 × 10⁻¹⁷ N**. r = m_p v/(eB) = **870 m**. T = 2πm_p/(eB) = **2.19 ms** (457 Hz, independent of v).

(b) Electron: r = m_e v/(eB) = **0.474 m**, f = eB/(2πm_e) = **840 kHz**. With B out of the page, F = qv × B points to the right of v for a positive charge: the proton circles **clockwise** seen with B toward you, the electron **anticlockwise**.

(c) F = ILB = 1 × 1 × 3 × 10⁻⁵ = **30 μN**, perpendicular to both the wire and the field, direction given by L × B with L along the current.

**Check.** r·ω_c = v: 870 m × 2π/2.19 ms = 2.5 × 10⁶ m/s. The mass ratio 1836 equals the radius ratio 870/0.474.

**Common mistakes.** Using the kinetic energy instead of the momentum in r; forgetting that the period does not depend on speed; getting the sense of rotation from the electron's charge wrong.

## 2. The loop axis, the solenoid and the torque on a loop

(a) Source at r′ = R(cos φ′, sin φ′, 0), dl = R(−sin φ′, cos φ′, 0)dφ′, field point (0, 0, z): r − r′ = (−R cos φ′, −R sin φ′, z), |r − r′| = √(R² + z²) for every φ′. The cross product dl × (r − r′) = R dφ′ (z cos φ′, z sin φ′, R sin²φ′ + R cos²φ′) = R dφ′ (z cos φ′, z sin φ′, R). The x and y components integrate cos φ′ and sin φ′ over a period: zero. The z component integrates to 2πR². So **B_z = (μ₀I/4π) 2πR²/(R² + z²)^{3/2} = μ₀IR²/(2(R² + z²)^{3/2})**. For |z| ≫ R, (R² + z²)^{−3/2} = |z|⁻³(1 + R²/z²)^{−3/2} = |z|⁻³(1 − (3/2)R²/z² + …), so B = (μ₀/2π)(IπR²)/|z|³ × **(1 − (3/2)(R/z)² + …)**. At z = 10R the correction is −1.5 %, which is exactly the 1.857 against 1.885 nT of the lesson.

(b) A slice dz′ at z′ carries current nI dz′: dB = μ₀nI R² dz′/(2(R² + (z − z′)²)^{3/2}). Substitute z − z′ = R tan u, dz′ = −R sec²u du: the integrand becomes (μ₀nI/2) cos u du, which integrates to (μ₀nI/2)[sin u] between the limits, sin u = (z − z′)/√((z − z′)² + R²). With z′ = ∓L/2 this is (4). Centre, z = 0: (μ₀nI/2) × 2 (L/2)/√(L²/4 + R²) → **μ₀nI** as L/R → ∞. End, z = L/2: (μ₀nI/2)[L/√(L² + R²) − 0] → **μ₀nI/2**.

(c) Put B along x, the rotation axis along y, the sides of length a parallel to y. The sides of length b lie in the x–z plane at y = ±a/2 with opposite currents; I L × B for them is along ±y, equal and opposite and on the same line (the loop's centre line): no force, no torque. The sides of length a carry currents ±I ŷ: F = ±IaB(ŷ × x̂) = ∓IaB ẑ, equal and opposite, net force zero. Their lines of action are separated, measured perpendicular to ẑ, by b sin θ. The couple is **τ = IaB · b sin θ = IAB sin θ**, and it turns n̂ toward B, i.e. it points along n̂ × B. With m = IA n̂ this is **τ = m × B**. A coil of N turns multiplies by N.

**Check.** Units: A·m² × T = A·m² × N/(A·m) = N·m. θ = 0: no torque (stable equilibrium); θ = 90°: maximum IAB.

**Common mistakes.** Integrating the magnitude |dB| instead of its components; dropping the minus sign of dz′ in the substitution and getting the limits reversed; computing the lever arm as b instead of b sin θ.

## 3. An electrodynamic tether (s = 15: L = 650 m)

**Assumptions.** Straight tether perpendicular to v and B; constant field; the ionospheric circuit is a fixed 100 Ω.

(a) **ℰ = vBL = 7600 × 3 × 10⁻⁵ × 650 = 148.2 V**; **I = 1.482 A**; **P = ℰI = 219.6 W**.

(b) **F = ILB = 1.482 × 650 × 3 × 10⁻⁵ = 0.0289 N**. The emf drives the current along v × B; the force on that current, IL × B with L ∥ v × B, is along (v × B) × B = −B²v (for v ⊥ B): **opposite to v**, as Lenz's law requires (the induced current opposes the motion that causes it).

(c) a = 0.0289/100 = **2.89 × 10⁻⁴ m/s²**, Δv per 5600 s orbit = **1.62 m/s**. A 100 kg satellite at 500 km needs roughly 100–200 m/s to deorbit (an estimate, not computed here), so about 100 orbits, a week: yes, a useful deorbit device that needs no propellant.

(d) Fv = (vBL/R)(LB)v = (vBL)²/R = ℰ²/R = P: the drag power becomes heat in the circuit (or useful power if a load replaces part of R). A battery driving current against the emf reverses I and therefore F: the tether becomes a thruster, paid for in electrical power.

**Check.** 0.0289 N × 7600 m/s = 219.6 W = ℰI.

**Common mistakes.** Using the orbital speed relative to the field as 7600 m/s without noting that the ionosphere co-rotates (a few per cent correction); forgetting that the circuit must close, without which no current flows and there is no force.

## 4. Size a torquer for a 1U CubeSat (s = 15: m_req = 0.25 A·m²)

(a) 400 km polar orbit: **|B| min 25.8 μT, max 51.6 μT, mean 39.8 μT**. Torque with m ⊥ B: **12.9 μN·m max, 6.46 μN·m min**.

(b) Side 9 cm, each design set to give 0.25 A·m²:

| N | d (mm) | I (A) | R (Ω) | P (W) | L (mH) | L/R (ms) |
|---|---|---|---|---|---|---|
| 100 | 0.20 | 0.3086 | 19.3 | 1.834 | 1.13 | 0.059 |
| 100 | 0.25 | 0.3086 | 12.3 | 1.174 | 1.13 | 0.092 |
| 100 | 0.30 | 0.3086 | 8.6 | 0.815 | 1.13 | 0.132 |
| 150 | 0.20 | 0.2058 | 28.9 | 1.223 | 2.54 | 0.088 |
| 150 | 0.25 | 0.2058 | 18.5 | 0.782 | 2.54 | 0.138 |
| 150 | 0.30 | 0.2058 | 12.8 | 0.543 | 2.54 | 0.198 |
| 200 | 0.20 | 0.1543 | 38.5 | 0.917 | 4.52 | 0.117 |
| 200 | 0.25 | 0.1543 | 24.6 | 0.587 | 4.52 | 0.184 |
| 200 | 0.30 | 0.1543 | 17.1 | 0.408 | 4.52 | 0.264 |
| 300 | 0.20 | 0.1029 | 57.8 | 0.611 | 10.18 | 0.176 |
| 300 | 0.25 | 0.1029 | 37.0 | 0.391 | 10.18 | 0.275 |
| 300 | 0.30 | 0.1029 | 25.7 | 0.272 | 10.18 | 0.397 |

Three designs meet 0.5 W: 200 turns of 0.30 mm (0.408 W, 72 m of wire, 45.6 g of copper, 2.64 V), 300 of 0.25 mm (0.391 W, 47.5 g, 3.81 V) and 300 of 0.30 mm (0.272 W, 68.4 g, 2.64 V). A reasonable choice is **200 turns of 0.30 mm**: lightest, and its 2.64 V drive fits a 3.3 V bus; the 300-turn 0.30 mm coil saves 0.14 W for 23 g more copper. Power at fixed moment scales as 1/(N d²) through R ∝ N/d² and I ∝ 1/N, so wire mass (∝ N d²) buys power directly. All time constants are below 0.4 ms.

(c) At the mean field 39.8 μT, the saturation rate m_max/(kB) = **18.0°/s > 5°/s**: the coil is **not** saturated at the start. **T_d = I/(kB²) = 63.2 s**; the rate is 4.23°/s at 10 s, 1.02°/s at 100 s, and falls below 0.1°/s at **247.5 s**, consistent with T_d ln 50 = 247 s. The lower bound Iω₀/(m_max B) = **17.5 s**: the unsaturated B-dot law uses only a fraction of the coil, and a higher gain would be faster.

**Check.** Time to 0.1°/s from (10): 63.2 × ln(5/0.1) = 247 s, as simulated.

**Common mistakes.** Computing power from the current needed by one design and the resistance of another; using the 2:1 polar-orbit maximum field for the detumble estimate rather than the mean; forgetting that the planar model ignores the field's rotation.

## 5. A Helmholtz pair (R = 0.1 m, 1 A)

(a) From `wire_field` with two 720-gons (deviation B(z)/B(0) − 1):

| d | B(0) (μT) | at z = 0.1R | at z = 0.2R |
|---|---|---|---|
| 0.8R | 10.058 | −4.09 × 10⁻³ | −1.71 × 10⁻² |
| R | 8.992 | −1.14 × 10⁻⁴ | −1.76 × 10⁻³ |
| 1.2R | 7.923 | +3.44 × 10⁻³ | +1.23 × 10⁻² |

(b) Let f(z) = (R² + z²)^{−3/2}; B(z) ∝ f(z − d/2) + f(z + d/2), even in z, so all odd derivatives vanish at 0. f″(z) = 3(4z² − R²)(R² + z²)^{−7/2} vanishes at z = R/2: with d = R the two loops sit at the inflection points of each other's field, B″(0) = 0, and the leading deviation is fourth order. Expanding to fourth order gives **B(z)/B(0) − 1 = −(144/125)(z/R)⁴**: −1.15 × 10⁻⁴ at 0.1R, matching −1.14 × 10⁻⁴; at 0.2R the formula gives −1.84 × 10⁻³ against −1.76 × 10⁻³, the sixth-order term showing. For d < R the deviation is negative and second order (the field peaks at the centre); for d > R it is positive and second order (the field has a shallow dip at the centre).

(c) **B(0) = (4/5)^{3/2} μ₀I/R = 8.9918 μT**; the 720-gon sum gives 8.9919 μT. To test a magnetometer, place it at the centre, null the Earth's field with three orthogonal pairs, then step the current in each pair to known values and fit the magnetometer's scale factors, offsets and axis misalignments; the uniformity of 10⁻⁴ over ±1 cm means the sensor's size does not matter.

**Common mistakes.** Driving the two loops in opposite senses (that is the anti-Helmholtz gradient coil); expecting second-order uniformity to vanish for any d.

## 6. Where a magnetorquer is the wrong tool

(a) At GEO, r = 42 164 km, the dipole gives **0.107 μT**, 230 times weaker than at 500 km. The lesson's 0.154 A·m² coil gives at most **1.6 × 10⁻⁸ N·m**, while solar radiation pressure on a 20 m² array with a 0.5 m centre-of-pressure offset gives 4.56 × 10⁻⁶ N/m² × 20 m² × 0.5 m = **4.6 × 10⁻⁵ N·m**, three thousand times more. Verdict: magnetorquers cannot even hold a GEO satellite against the disturbances; use wheels with thrusters for momentum dumping.

(b) On the coil axis at 10 cm the dipole estimate is μ₀m/(2πz³) = **30.8 μT** (15.4 μT in the coil plane), as large as the Earth's 25–50 μT. Flight software must either sample the magnetometer only while the torquer is off (duty cycling, with a settling time of several L/R plus eddy currents in the structure) or subtract a calibrated model of the coil field. The dipole estimate is rough at 10 cm (z ≈ 2.5 times the coil half-size), but the conclusion does not depend on a factor of two.

(c) Every coil's torque is mᵢ × B, perpendicular to B, so the sum is perpendicular to B: the component of the desired torque along B is unreachable at that instant. Over an orbit the direction of B in the body frame rotates (twice per orbit in latitude on a polar orbit, as the dipole field turns from horizontal to vertical), so a direction unreachable now is reachable later, and the time-average controllability is complete.

(d) The orbital motion through the dipole field (B rotates in inertial space at about twice the orbital rate on a polar orbit, so a non-tumbling satellite still sees dB/dt) and the noise of the magnetometer (plus the torquer's own field from (b)). At 400 km the period is 5544 s, the orbital rate 0.065°/s and the rotation of B about twice that, 0.13°/s; B-dot cannot drive the body rate much below the rate at which B itself rotates, and the controller is usually switched to a different law there.

## Check-yourself answers

1. F = q(E + v × B); the magnetic part is perpendicular to v, so it does no work and cannot change |v|.
2. μ₀I/(2R) at the centre; μ₀IR²/(2|z|³) = μ₀m/(2π|z|³) far away: the magnetic moment m = IπR² is the only property of the loop the far field remembers.
3. The end of a long solenoid sees windings on one side only; by superposition two half-solenoids joined end to end give the full μ₀nI, so each gives half.
4. 62 μT at the pole (twice the equatorial value); 31 × (6371/6871)³ = 24.7 μT at 500 km over the equator.
5. mB = 0.15 × 40 × 10⁻⁶ = 6 μN·m, with m ⊥ B; none about the local field direction.
6. In the sense whose field opposes the increase (it reduces the flux): Lenz's law, the minus sign in Faraday's law.
7. τ = 4 × 10⁻³/22 = 0.18 ms; 1 − e⁻¹ = 63.2 % of the final current.
8. Unsaturated: T_d = I/(kB²), so the time scales as 1/B²; saturated: the torque is m_max B, so the time scales as 1/B.
