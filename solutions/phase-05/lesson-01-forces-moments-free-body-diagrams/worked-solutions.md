# P05.01  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` (and the probe script that calls it) and given to the figures shown.

## 1. Moments in three dimensions (s = 15: F = (−1000, 3000, 1500) N)

**Assumptions.** Rigid body; the force acts at r = (0, −1.2, 0.3) m; moments about the origin.

(a) M = r × F = (r_y F_z − r_z F_y, r_z F_x − r_x F_z, r_x F_y − r_y F_x) = ((−1.2)(1500) − (0.3)(3000), (0.3)(−1000) − 0, 0 − (−1.2)(−1000)) = **(−2700, −300, −1200) N·m**, |M| = **2969.8 N·m**. `moment_about_point` returns the same vector.

(b) About the x axis: M_x = x̂ · M = **−2700 N·m**. About the axis along (1, 0, 1)/√2: M_u = (−2700 − 1200)/√2 = **−2757.7 N·m**.

(c) M · F = 0 and M · r = 1.3 × 10⁻¹⁴ (round-off): a cross product is perpendicular to both factors, so the moment lies in the plane normal to the plane containing r and F.

**Check.** Dimensional: m × N = N·m. The x-axis moment comes only from F_y and F_z (an x force is parallel to the axis and contributes nothing): −1.2 × 1500 − 0.3 × 3000 = −2700.

**Common mistakes.** Swapping the order (F × r reverses every sign); forgetting to normalise the axis vector (gives √2 times the answer).

## 2. Equivalent force systems

(a) Let P′ = P + λû. Then (r − P′) × F = (r − P) × F − λ û × F, and û · (û × F) = 0 because a triple product with a repeated vector vanishes. So **û · ((r − P′) × F) = û · ((r − P) × F)** for every λ.

(b) M_B = Σ (r_i − B) × F_i = Σ ((r_i − A) + (A − B)) × F_i = M_A + (A − B) × Σ F_i = **M_A + (A − B) × R**. For a couple, R = F − F = 0, so M_B = M_A about every B.

(c) In the plane z = 0, R = (R_x, R_y, 0) and M_O = (0, 0, M). A single force R placed at a point P has moment P × R about O; we need P × R = M_O. Take p = R × M_O/|R|². Expanding, p × R = (R × M_O) × R/|R|² = (R · R) M_O/|R|² − (M_O · R) R/|R|² = M_O, since M_O · R = 0 for a planar system. Every P = p + λR also works (sliding along the line), so the system is **one force R on the line through p along R**, with zero pitch. If R = 0 the system is a **pure couple** M_O, the same about every point, and has no line of action.

**Check.** The test `test_central_axis_carries_a_moment_parallel_to_the_resultant` places R at p for random planar systems and recovers M_O to 10⁻¹⁰.

**Common mistakes.** Treating the moment about an axis as the magnitude of the moment vector; forgetting that (5) requires the same R on both sides.

## 3. A rocket on its launch rail (s = 15: m = 13.0 kg)

**Assumptions.** Rigid rocket at rest; buttons frictionless (normal force only); stop takes only an axial force; g = 9.80665 m/s².

FBD: weight W = 127.49 N downward at the CG; button forces N₁, N₂ across the rail at 0.3 m and 1.1 m; stop force S along the rail at the aft end. Three unknowns, three equations.

Axes: a = (cos 80°, sin 80°) up the rail, n = (−sin 80°, cos 80°) across it. W has components W_a = −W sin 80° = −125.55 N and W_n = −W cos 80° = −22.14 N.

- Along the rail: **S = 125.55 N**.
- Moments about the stop (it drops out): the weight at (0.9 a + 0.05 n) has moment 0.9 W_n − 0.05 W_a = −19.92 + 6.28 = −13.65 N·m. The buttons' moments are 0.3 N₁ and 1.1 N₂, so 0.3 N₁ + 1.1 N₂ = 13.65.
- Across the rail: N₁ + N₂ = 22.14.

Solving: **N₂ = 8.76 N, N₁ = 13.38 N**. `solve_fbd_2d` returns (13.381, 8.756, 125.550) N with a residual of 3 × 10⁻¹⁴ N.

(d) The forward button carries less than the aft one here because the CG's offset from the rail line makes the axial weight component W_a act with a lever of 0.05 m, a moment of 6.28 N·m that opposes the moment of W_n. Without the offset the buttons would carry 5.53 N (aft) and 16.60 N (forward): the 5 cm offset reverses which button is more loaded. Small offsets matter whenever a large force (here the along-rail weight, six times the across-rail one) acts with a small lever.

**Check.** N₁ + N₂ = 22.14 N = W cos 80°; S = W sin 80°.

**Common mistakes.** Putting the weight on the rail line (dropping the 0.05 m offset); resolving the weight with sin and cos swapped for an angle measured from the horizontal.

## 4. A raked main gear (s = 15: m = 900 kg)

`landing_load_case(900)` gives V = 11 915.1 N and, at D/V = 0.25, D = 2 978.8 N. With C = (0.10, 0):

| D/V | A_x (N) | A_y (N) | T (N) |
|---|---|---|---|
| 0 | 2 383.0 | −10 326.4 | −2 864.0 |
| 0.25 | −1 787.3 | −15 092.4 | 5 728.1 |
| 0.5 | −5 957.5 | −19 858.5 | 14 320.1 |

`brace_force_by_moments` gives the same three brace forces.

(b) With the contact 0.10 m ahead of the trunnion, V has a moment about A of 0.10 V = 1 191.5 N·m (anticlockwise seen from starboard: it tends to swing the foot of the leg forward about A). The brace must supply −1 191.5 N·m, and with its lever of 0.4160 m that is T = −2 864.0 N: **compression**. A brace on a raked leg is a strut on every landing, even with no drag.

(c) The brace is idle when the ground force passes through A: the line through C with direction (−D, V) reaches x = 0 at height 1.2 m when D/V = 0.10/1.2 = **0.0833**. The code gives T = 9 × 10⁻¹⁴ N there.

(d) For the vertical leg of the lesson at the same mass and D/V = 0.25, T = 8 592.1 N (tension); the raked leg's 5 728.1 N is a third less, because the rake makes V oppose the drag's moment. But the raked brace reverses sign over the load range (compression at low drag, tension at high), so it must be sized as a strut (buckling) and checked for fatigue under reversal. The choice needs the full set of load cases, the buckling and fatigue allowables and the space available.

**Check.** At D/V = 0.25: T × d = 5 728.1 × 0.4160 = 2 383 N·m = D ℓ − 0.10 V = 3 574.5 − 1 191.5.

**Common mistakes.** Keeping C at the origin of the old geometry; measuring the rake with the wrong sign (aft instead of ahead).

## 5. Wing lift distributions

(a) With y = s sin u, dy = s cos u du and √(1 − y²/s²) = cos u:
W = w₀ s ∫₀^{π/2} cos²u du = **π w₀ s/4**, and ∫ y w dy = w₀ s² ∫₀^{π/2} sin u cos²u du = w₀ s²/3. So **ȳ = (w₀s²/3)/(πw₀s/4) = 4s/(3π) = 0.4244 s**.

(b) Split the chord into a rectangle of height λc_r (centroid s/2) and a triangle of height (1 − λ)c_r falling to zero at the tip (centroid s/3). Areas λc_r s and (1 − λ)c_r s/2:
ȳ = [λ s²/2 + (1 − λ)s²/6]/[λ s + (1 − λ)s/2] = **s(1 + 2λ)/(3(1 + λ))**. For λ = 0.5 and s = 5 m, ȳ = 2.2222 m.

(c) Schrenk distribution on s = 5 m, 5000 N:

| points | W (N) | ȳ (m) |
|---|---|---|
| 41 | 4 996.30 | 2.16972 |
| 161 | 4 999.54 | 2.17186 |
| 641 | 4 999.94 | 2.17211 |
| 2561 | 4 999.99 | 2.17214 |

The average of the closed-form centroids is (2.12207 + 2.22222)/2 = **2.17214 m**, matched at 2561 points. (The total converges from below because the trapezoidal rule cuts off the elliptic part's steep tip.)

(d) Root bending moment = L ȳ: elliptic 10 610.3 N·m, chord-proportional 11 111.1 N·m, Schrenk **10 860.7 N·m**. The spar is sized by this moment (`P05.10`); a distribution that puts the lift further outboard costs spar weight, and choosing the elliptic one when the real wing behaves like Schrenk would underestimate the root moment by 2.3 %.

**Common mistakes.** Averaging w₀ values instead of distributions scaled to the same total; forgetting that the chord-proportional distribution's total must also be set to 5000 N.

## 6. What the free-body diagram left out

(a) **Spin-up and spring-back**: the drag rises to a peak while the tyre skids up to speed, then the stored strain energy in the leg and brace throws the wheel forward: D reverses and the brace goes into compression. **Dynamic overshoot of V**: the oleo strut's response can make the peak vertical load exceed a quasi-static estimate; V does not load this brace (it passes through A) but it loads the trunnion and, on a raked leg, the brace too.

(b) The brace is no longer a two-force member: it is a beam loaded at three points (two pins and the actuator attachment). Its end forces are not along the pin line, so at B the leg sees a force with two unknown components, and the leg alone has four unknowns. You must take the brace's own free-body diagram too (a frame, `P05.04`) to close the problem.

(c) A trunnion that resists a moment M_A takes part of the drag moment: ΣM_A gives T d = D ℓ − M_A, so **the brace force goes down**, by an amount that depends on the relative stiffness of the joint and the brace. The problem becomes statically indeterminate (`P05.02`).

(d) Taking moments about C is legitimate, but then the trunnion reactions have levers (A_x has a lever of 1.2 m about C) and must be kept, and the brace force's lever must be measured from C, not from A. The two likely mistakes are dropping A_x from the moment equation, and reusing the lever d measured from A.

## Check-yourself answers

1. M = F d = 2000 × 0.3 = 600 N·m; sliding the force along its line changes nothing, because the perpendicular distance is unchanged.
2. Moving the point along the axis changes the lever by a multiple of û, whose cross product with F is perpendicular to û and drops out of the dot product.
3. M_B = M_A + (A − B) × R; once R = ΣF = 0, the moment sum is the same about every point, so you may take moments wherever most unknowns vanish.
4. Two equal, opposite, parallel forces; their resultant is zero, so by the transport rule their moment does not depend on the point.
5. Uniform: at L/2 with W = w₀L. Triangular: at 2L/3 from the zero end, W = w₀L/2. Elliptic: at 4s/(3π) = 0.4244 s, W = πw₀s/4.
6. With forces only at its two pins, equilibrium of the brace needs those two forces equal, opposite and collinear: along the pin line.
7. The trunnion reactions act at A and have no lever about A. What remains are the drag (lever ℓ), the vertical load (lever zero when C is below A) and the brace (lever d).
8. When the matrix is singular: the three unknown forces are concurrent or all parallel. The supports then cannot resist a moment about the meeting point (or a force across the parallel lines): the body is improperly constrained.
