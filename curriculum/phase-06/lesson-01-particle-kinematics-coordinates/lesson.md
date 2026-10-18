# P06.01  Particle Kinematics in Many Coordinates

> One motion, three descriptions. Pick the coordinates the measurement comes in, and remember that their unit vectors turn.

| | |
|---|---|
| **Learning outcomes** | 1. Derive velocity and acceleration in polar and cylindrical coordinates from the time derivatives of the turning unit vectors, and name each term (radial, centripetal, transverse, Coriolis-like 2ṙθ̇). 2. Describe the same motion in path (normal–tangential) coordinates: speed, tangential acceleration, normal acceleration, curvature and the tangent–normal–binormal triad, and convert between Cartesian, cylindrical and path descriptions. 3. Implement and test the conversions, and estimate velocity and acceleration from sampled tracking data by central differences, by a cylindrical and a Cartesian route. 4. Show that the estimates converge to the closed-form trajectory at second order, and quantify how measurement noise sets an optimal difference step. 5. Decide what a tracking record can and cannot tell you about an aircraft's acceleration. |
| **Prerequisites** | `P04.01` (kinematics and Newton's laws), `P01.06` (parametric and polar curves), `P01.07` (vectors, dot and cross products), `P01.15` (rotation matrices: a change of basis), `P03.04` (numerical differentiation, truncation error), `P02.16` (random variables: the standard deviation of a sum), `P00.14` (measurement uncertainty). Skill check: you can differentiate r(t) cos θ(t) twice by the product and chain rules. |
| **Engineering question** | A ground radar reports an aircraft's horizontal range, azimuth and height ten times a second while the aircraft flies a climbing turn. From those samples alone, what are its speed, its normal acceleration and the radius of its turn, and how good are those numbers when the radar's azimuth is noisy? The answers: 200.25 m/s, 20.0 m/s² (2.04 g) and 2005 m from noise-free samples to 0.1 %; with 1 mrad of synthetic azimuth noise, no better than about 1.1 m/s² RMS (5 %), at a difference step near 6 s. |
| **Assumptions and limits** | A particle: the aircraft's attitude and size are ignored. A flat, non-rotating Earth over the 12 km of the example (lesson `P06.04` adds rotation). Uniform sampling. The noise is synthetic, Gaussian and independent from sample to sample; real radar errors are correlated and include biases. Cylindrical coordinates about the radar fail on the axis (r → 0: the azimuth rate blows up as an aircraft passes overhead). |
| **Analytical work** | ė_r = θ̇ e_θ and ė_θ = −θ̇ e_r; the polar and cylindrical velocity and acceleration; the inverse conversion from Cartesian rates; a_t = v·a/\|v\|, a_n = \|v × a\|/\|v\|, κ = \|v × a\|/\|v\|³; the helix curvature R/(R² + b²); why straight-line flight has zero polar acceleration; the truncation and noise errors of a central second difference and the optimal step. |
| **Computational deliverable** | `starter/main.py`: `polar_kinematics`, `cylindrical_kinematics`, `cylindrical_to_cartesian`, `cartesian_to_cylindrical`, `path_kinematics`, `helix_turn`, `helix_curvature`, `central_differences`, `track_from_radar`, `track_cartesian`, `optimal_step` (`radar_samples`, `observed_order`, `accel_error_vs_stride` and `main` are given). Tests prove: the polar components against direct differentiation; uniform circular motion is purely centripetal; straight-line flight gives a_r = a_θ = 0 to 10⁻¹⁰; Cartesian ↔ cylindrical round trips; the conversion matrix is a rotation; a_t² + a_n² = \|a\|² with an orthonormal triad and no binormal acceleration; the helix closed form is self-consistent; finite differences of radar samples converge at order 2.0 ± 0.1 by both routes; the azimuth is unwrapped through ±π; synthetic noise makes small steps worse and the optimum falls within a factor 2 of the model. |
| **Practical activity** | **computational experiment**: a convergence study of the two differentiation routes against the closed-form helix, then a step-size sweep on a synthetic noisy radar record, compared with an error model, and a sample-rate comparison. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `radar_track.png`, the derivation of the polar acceleration and of the optimal step, the convergence and noise tables, and a half-page note answering the engineering question with its uncertainty. |
| **Assessment** | `assessments/phase-06/grading-rubric.md`; competencies 1–4. Pass: the polar acceleration is derived from the unit-vector derivatives (not quoted); all public tests pass; the observed orders are reported for both routes; the note gives the normal acceleration with an uncertainty from the noise study and says which step was used and why. |
| **References** | `@mit-ocw-16-07-dynamics` L5 "Other coordinate systems": Polar Coordinates (r − θ), Velocity vector, Acceleration vector, Change of basis, Cylindrical Coordinates (r − θ − z); L6 "Intrinsic coordinates": Intrinsic coordinates: Tangential, Normal and Binormal components (required-free: the unit-vector derivatives and the path components, with an aircraft helix example). `@openstax-university-physics-vol1` 4.4 Uniform and Nonuniform Circular Motion: Centripetal Acceleration; Total Acceleration during Circular Motion (required-free). `@mit-ocw-2-003sc-engineering-dynamics` Lecture 4: Movement of a Particle in Circular Motion w/ Polar Coordinates (supplementary-free, video). `@openstax-calculus-volume-3` 3.4 Motion in Space and `@openstax-calculus-volume-2` 7.3 Polar Coordinates (supplementary-free: the calculus of space curves and polar geometry). |
| **Extensions** | ★ Replace the three-point stencil by a least-squares fit of a quadratic over a sliding window (a Savitzky–Golay filter) and show that, unlike the stencil, it benefits from a higher sample rate. ★ Add spherical coordinates (slant range, azimuth, elevation), which is what a real radar measures, and repeat the convergence study. ★ Estimate the torsion of the helix from the third derivative and compare with b/(R² + b²). |
| **Time** | unmeasured |

## 1. The engineering problem

A ground radar does not measure x, y and z. It measures where the target is *from the radar*: a range and an angle (and, for a three-dimensional radar, a height or an elevation). An air-traffic or range-safety system then has to answer questions that are posed in other terms. How fast is the aircraft going? Is it turning, and how hard? Will it pass inside a keep-out circle? A pilot and a structures engineer ask the same questions in yet another frame: they care about the acceleration along the flight path, which changes the speed, and the acceleration across it, which bends the path and loads the wing. That cross-path acceleration divided by g is most of the load factor of lesson `P06.02`.

So one motion has at least three useful descriptions: **Cartesian** (fixed axes, simple calculus), **cylindrical** (the coordinates the radar reports), and **path** or **normal–tangential** coordinates (the ones the aircraft feels). The difficulty is that the unit vectors of the second and third descriptions turn as the aircraft moves. Differentiate the position without accounting for that and you get the wrong acceleration, sometimes by more than the acceleration itself.

The test case is a fighter-sized aircraft (illustrative numbers) flying a climbing turn: radius 2000 m about a vertical axis 8 km east and 6 km north of the radar, turn rate 0.1 rad/s (200 m/s horizontally), climbing at 10 m/s from 3000 m. Success criteria, written first: the cylindrical and path descriptions reproduce the closed-form helix exactly; differences of noise-free radar samples converge to it at second order by both routes; the normal acceleration and turn radius are recovered to 0.1 %; and the noise study finds the best difference step and the error it leaves, in agreement with a stated error model.

## 2. Understand and derive

### 2.1 Turning unit vectors

In the horizontal plane put the radar at the origin, x east, y north, and describe the aircraft by its horizontal range r and azimuth θ (measured from east toward north). Define the unit vectors

e_r = (cos θ, sin θ),  e_θ = (−sin θ, cos θ).

They depend on θ, and θ depends on time. Differentiating with the chain rule:

**ė_r = θ̇ e_θ,  ė_θ = −θ̇ e_r.**   (1)

Each turns at the rate θ̇, and its derivative is perpendicular to it (a unit vector can only rotate). This is the whole content of curvilinear kinematics; the rest is the product rule (16.07 L5, "Polar Coordinates").

### 2.2 Polar and cylindrical velocity and acceleration

The position is **r** = r e_r. Differentiate once, using (1):

**v = ṙ e_r + r θ̇ e_θ.**   (2)

Differentiate again; each term of (2) gives two terms:

d(ṙ e_r)/dt = r̈ e_r + ṙ θ̇ e_θ,  d(r θ̇ e_θ)/dt = ṙ θ̇ e_θ + r θ̈ e_θ − r θ̇² e_r,

so

**a = (r̈ − r θ̇²) e_r + (r θ̈ + 2 ṙ θ̇) e_θ.**   (3)

Name the four terms. r̈ is the acceleration of the range itself. −rθ̇² is the **centripetal** term: a particle at constant r going round at θ̇ accelerates inward at rθ̇². rθ̈ is the angular acceleration term. 2ṙθ̇ is the term students forget: one ṙθ̇ comes from the radial velocity being rotated, the other from the transverse velocity rθ̇ growing because r grows. (In lesson `P06.04` the same structure reappears as the Coriolis acceleration.) Units check: every term is m/s².

**Cylindrical coordinates** add a height z along a fixed axis e_z, which does not turn, so the extra components are simply ż and z̈:

v = (ṙ, rθ̇, ż),  a = (r̈ − rθ̇², rθ̈ + 2ṙθ̇, z̈)   in (e_r, e_θ, e_z).   (4)

**Converting** components to Cartesian is a rotation by θ about z (`P01.15`): the rows of the matrix are e_r, e_θ, e_z, it is orthonormal with determinant +1, and its transpose converts back (16.07 L5, "Change of basis"). The **inverse problem**, getting the cylindrical rates from Cartesian v and a, follows from projecting on e_r and e_θ and inverting (2) and (3):

r = √(x² + y²), θ = atan2(y, x), ṙ = v·e_r, θ̇ = (v·e_θ)/r, r̈ = a·e_r + rθ̇², θ̈ = (a·e_θ − 2ṙθ̇)/r.   (5)

Every one of these divides by r: on the axis (an aircraft directly over the radar) the azimuth and its rates are undefined.

**A worked example: straight and level flight.** An aircraft passes the radar on a straight line at constant velocity. Its range first falls, then rises; its azimuth swings fastest at closest approach. Neither r nor θ is linear in time. Yet a = 0, so by (3) r̈ = rθ̇² and rθ̈ = −2ṙθ̇ at every instant. The second identity says d(r²θ̇)/dt = r(rθ̈ + 2ṙθ̇) = 0: **r²θ̇ is constant on a straight line**, the geometric fact behind Kepler's second law. For an aircraft at 200 m/s passing at a miss distance d, θ̇ peaks at V/d at closest approach and |θ̈| peaks at (3√3/8)(V/d)²: 0.040 rad/s and 0.0010 rad/s² at 5 km, but 4.0 rad/s and 10.4 rad/s² at 50 m. That is why tracking antennas cannot follow an overhead pass, and why cylindrical coordinates about the radar are a poor choice for a target near the axis.

### 2.3 Path coordinates

The path description uses the velocity itself as the first axis. With speed V = |v| and unit tangent e_t = v/V, differentiate v = V e_t: a = V̇ e_t + V ė_t. The unit tangent turns as the path curves; its rate is V κ e_n, where κ is the **curvature** (1/ρ, ρ the radius of curvature) and e_n the **principal normal**, pointing to the centre of curvature (16.07 L6). So

**a = V̇ e_t + (V²/ρ) e_n,  a_t = V̇ = v·a/V,  a_n = V²/ρ = |v × a|/V,  κ = |v × a|/V³.**   (6)

The third vector, the **binormal** e_b = e_t × e_n, is perpendicular to the plane in which the path is momentarily turning (the osculating plane); the acceleration has no component along it. Because (6) is built from v and a alone, it does not care which basis they were written in: compute a_t and a_n from Cartesian or from cylindrical components and you must get the same numbers, and a_t² + a_n² = |a|² always (the tests check both).

OpenStax 4.4 (Centripetal Acceleration; Total Acceleration during Circular Motion) is the special case of a circle: a_n = V²/r toward the centre and a_t = dV/dt along the path.

### 2.4 The climbing turn in closed form

The aircraft's path is a helix about the vertical line through the turn centre (c_x, c_y):

x = c_x + R cos ωt,  y = c_y + R sin ωt,  z = z₀ + w t,

with R = 2000 m, ω = 0.1 rad/s and climb rate w = 10 m/s. Differentiating by hand: v = (−Rω sin ωt, Rω cos ωt, w), a = (−Rω² cos ωt, −Rω² sin ωt, 0). The speed is constant, V = √((Rω)² + w²) = **200.250 m/s**, so a_t = 0. The acceleration is horizontal, Rω² = **20.000 m/s² (2.039 g)**, toward the turn axis. From (6), with the pitch parameter b = w/ω = 100 m, the curvature is

**κ = R/(R² + b²) = 4.9875 × 10⁻⁴ m⁻¹,  ρ = 2005.00 m:**   (7)

slightly larger than R because a climbing path is less curved than its circular shadow. The normal acceleration V²κ = 20.000 m/s² equals |a|, as it must when a_t = 0.

**The same instant in three descriptions.** At t = 12 s the aircraft is at (8724.7, 7864.1, 3120.0) m with v = (−186.408, 72.472, 10.000) m/s and a = (−7.2472, −18.6408, 0) m/s². The radar sees r = 11 745.8 m, θ = 42.030°, ṙ = −89.941 m/s, θ̇ = 0.015208 rad/s, r̈ = −15.1468 m/s², θ̈ = −5.328 × 10⁻⁴ rad/s². The cylindrical components are v = (−89.941, 178.635, 10.000) and a = (−17.864, −8.994, 0). Look at the terms: rθ̇² = 2.717 m/s² and 2ṙθ̇ = −2.736 m/s², each about 14 % of the whole acceleration. Drop either and the radar would report a wrong turn. The path description gives V = 200.250 m/s, a_t = 0 (8.9 × 10⁻¹⁶ in floating point), a_n = 20.000 m/s² and ρ = 2005.00 m, identical to (7).

### 2.5 Derivatives from samples, and what noise does

The radar gives samples, not functions. A central difference over a step h (a stride of k samples at spacing Δt, h = kΔt) estimates

f′ ≈ (f₊ − f₋)/(2h),  f″ ≈ (f₊ − 2f₀ + f₋)/h²,

both with truncation errors of order h² (`P03.04`): (h²/6) f‴ and (h²/12) f⁗. Now add independent noise of standard deviation σ to every sample. The second difference combines three samples with weights 1, −2, 1, so its noise has standard deviation √(1 + 4 + 1) σ/h² = √6 σ/h² (`P02.16`). The total error is roughly

E(h) ≈ (h²/12)|f⁗| + √6 σ/h²,   (8)

which falls and then rises. Setting dE/dh = 0 puts the two terms equal at

**h* = (12√6 σ/|f⁗|)^{1/4},  E(h*) = 2√(√6 σ |f⁗|/12).**   (9)

For the turn, a Cartesian coordinate is R cos ωt plus a straight line, so |f⁗| = Rω⁴ = 0.2 m/s⁴. A 1 mrad azimuth error at the mean range of 10.05 km is 10.0 m across the line of sight; with 5 m in range, the error per horizontal axis is about √((5² + 10²)/2) = 7.93 m. Then (9) gives h* = 5.84 s and an error bound of 1.14 m/s² (the two terms added). For the whole vector, add the parts in quadrature at h = 6 s: truncation 0.60 m/s², noise 0.54 m/s² on each horizontal axis and 0.34 m/s² vertically (σ_z = 5 m), together 1.03 m/s². Either way: about 5 % of the turn acceleration, however clever the arithmetic, unless more samples are used per estimate.

### 2.6 Where the model stops

- **Overhead passes**: the cylindrical description is singular on the axis (2.2). Track in Cartesian or spherical coordinates near it.
- **Correlated noise and bias**: (8) assumes independent errors. A slowly wandering bias does not average out and is not reduced by any step.
- **Manoeuvre onset**: a central difference with h = 6 s smears a manoeuvre that starts within ±6 s of the estimate. A wide step trades noise for lag.
- **Earth curvature and rotation**: over 12 km the flat-Earth error is small; over hundreds of kilometres it is not (`P06.04`).

## 3. Build the computational model

Open `starter/main.py`. `radar_samples`, `observed_order`, `accel_error_vs_stride` and `main` are given.

### Step 1: polar and cylindrical components

`polar_kinematics(r, rd, rdd, th, thd, thdd)` returns (2) and (3) as (v_r, v_θ, a_r, a_θ). `cylindrical_kinematics` adds (ż, z̈) and stacks the results with shape (…, 3). Everything must accept arrays.

```python
return rd, r * thd, rdd - r * thd ** 2, r * thdd + 2.0 * rd * thd
```

### Step 2: conversions

`cylindrical_to_cartesian(comp, th)` applies the rotation: x = c_r cos θ − c_θ sin θ, y = c_r sin θ + c_θ cos θ, z = c_z. `cartesian_to_cylindrical(pos, vel, acc)` is (5).

### Step 3: path coordinates

`path_kinematics(vel, acc)` returns a dict with speed, a_t, a_n, κ, ρ (inf on a straight line) and e_t, e_n, e_b from (6). Get e_n from the part of a perpendicular to e_t, normalised, and e_b = e_t × e_n.

### Step 4: the closed-form trajectory

`helix_turn(t, ...)` returns position, velocity and acceleration of 2.4, differentiated by hand; `helix_curvature(R, omega, climb)` is (7).

### Step 5: differences of samples, two routes

`central_differences(t, f, k)` returns (t[k:−k], f′, f″) with stride k. `track_from_radar(t, r, th, z, k)` differences r, θ and z, applies `cylindrical_kinematics` and converts to Cartesian: the **cylindrical route**. `track_cartesian` converts the samples to x, y, z first and differences those: the **Cartesian route**. The two are not identical (differencing does not commute with the nonlinear map from (r, θ) to (x, y)), but both are O(h²) from the truth. `radar_samples` unwraps the azimuth with `np.unwrap`: a jump from +π to −π would otherwise look like a 2π/h azimuth rate.

### Step 6: the error model

`optimal_step(sigma, f4)` is (9).

Running `main()` prints the closed form, the instant at t = 12 s in three descriptions, the convergence table and the noise sweep, and writes `radar_track.png` (plan view of the radar and the track; r(t) and θ(t) as the radar sees them; the error against step for the noise-free and noisy cases on log axes).

### Verify

Run `python -m pytest curriculum/phase-06/lesson-01-particle-kinematics-coordinates`. The tests:

1. **Polar components** against central differences of x = r cos θ, y = r sin θ for an arbitrary r(t), θ(t).
2. **Uniform circular motion**: v_θ = Rω, a_r = −Rω², a_θ = 0.
3. **Straight-line flight**: r and θ vary strongly, yet a_r and a_θ vanish to 10⁻¹⁰ and the speed is preserved.
4. **Round trips** Cartesian → cylindrical → Cartesian for 20 random states, speed the same in both bases; the conversion matrix is a proper rotation.
5. **Path coordinates**: a_t² + a_n² = |a|², orthonormal triad, no binormal component, a = a_t e_t + a_n e_n, κV² = a_n; a straight line has ρ = ∞; a counterclockwise circle has e_b up.
6. **Helix**: differences of the closed form agree with its derivatives to their truncation error; κ = R/(R² + b²); zero climb gives a circle.
7. **Sampled radar**: both routes converge at order 2.0 ± 0.1 in velocity and acceleration; they agree with each other; the turn radius and a_n are recovered to 10⁻³; a track crossing azimuth ±π is handled.
8. **Noise**: a seeded synthetic record gives an error at h = 0.5 s more than 20 times the best, rising again at large h, with the best step within a factor 2 of (9); the noise generator is reproducible.
9. `main` writes the figure.

### Compare

The independent check is the closed-form helix: the radar samples are generated from it, differenced, and compared with its exact derivatives. The noise-free errors, RMS over 120 s, are

| h (s) | cylindrical route (m/s²) | Cartesian route (m/s²) |
|---|---|---|
| 4.0 | 3.111 × 10⁻¹ | 2.653 × 10⁻¹ |
| 2.0 | 7.814 × 10⁻² | 6.658 × 10⁻² |
| 1.0 | 1.953 × 10⁻² | 1.666 × 10⁻² |
| 0.5 | 4.878 × 10⁻³ | 4.166 × 10⁻³ |

Observed orders 1.993, 2.000, 2.001 (cylindrical) and 1.994, 1.999, 2.000 (Cartesian): second order, as the truncation term predicts. At h = 1.0 s the cylindrical route's average normal acceleration is 19.984 m/s² and its radius of curvature 1999.8 m, both within 0.3 % of the closed form, with |a_t| below 0.012 m/s². The cylindrical route is 17 % less accurate here because r(t) and θ(t) have larger fourth derivatives than x(t) and y(t); the radar's own coordinates are not always the best ones to difference in.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** All data here are synthetic: generated from the closed-form helix, with synthetic Gaussian noise where stated. Nothing is a radar measurement. Write each prediction down before running.

1. **Convergence.** Difference the noise-free samples at h = 4, 2, 1 and 0.5 s by both routes. Predict the observed order and which route is more accurate. *Finding:* the table above; order 2.00 for both, the Cartesian route about 15 % more accurate.
2. **Step against noise.** Add synthetic noise σ_r = 5 m, σ_θ = 1 mrad, σ_z = 5 m (seed 1) and sweep the step from 0.5 s to 14 s. Predict the best step from (9) before running. *Finding:* RMS errors (m/s²) of 125.1, 30.2, 7.86, 3.49, 1.90, 1.32, 1.09, 1.30, 1.86 and 3.38 at h = 0.5, 1, 2, 3, 4, 5, 6, 8, 10 and 14 s: the best is 1.086 m/s² at h = 6 s (5.4 % of the turn acceleration); the model predicts h* = 5.84 s and 1.03 m/s² (quadrature sum, 2.5).
3. **A different noise level.** Double σ_θ to 2 mrad (Problem 4). Predict how far h* moves: by (9), only by the fourth root of the noise ratio.
4. **Sample faster.** Repeat experiment 2 at 20 Hz instead of 10 Hz with the same noise. Predict the change. *Finding (seed 15):* best error 1.084 m/s² at h = 6 s at 20 Hz against 1.087 m/s² at 10 Hz: no improvement, because the three-point stencil uses three samples whatever the rate. More samples help only if the estimator uses them (Extension 1).

Raw data format: CSV tables `h_s, err_cyl_mps2, err_cart_mps2` (1) and `h_s, rms_err_mps2` (2–4), each with a header line stating that the data are SYNTHETIC and giving the seed; plus `radar_track.png`.

## 5. Validate and improve

| quantity | closed form (§2) | code | difference | cause |
|---|---|---|---|---|
| speed (m/s) | 200.250 | 200.250 (`path_kinematics`) | < 10⁻¹² | — |
| radius of curvature (m) | 2005.00 (7) | 2005.00 from exact v, a; 1999.8 from samples, h = 1 s | 0; 0.26 % | truncation, O(h²) |
| normal acceleration (m/s²) | 20.000 | 20.000; 19.984 from samples, h = 1 s | 0; 0.08 % | truncation |
| acceleration error, h = 0.5 s (m/s²) | (h²/12) f⁗ ≈ 4 × 10⁻³ | 4.9 × 10⁻³ (cyl.), 4.2 × 10⁻³ (Cart.) | factor 1.2 | f⁗ of r, θ larger than of x, y |
| best step with noise (s) | 5.84 (9) | 6.0 (sweep) | 3 % | one seeded record, flat minimum |
| error at the best step (m/s²) | 1.03 (quadrature, h = 6 s) | 1.09 RMS | 6 % | noise added in (r, θ), not x, y; one record |

The engineering answer: **from noise-free samples the aircraft's speed is 200.25 m/s, its normal acceleration 20.0 m/s² (2.04 g) and its turn radius 2005 m, recovered to better than 0.3 % with a 1 s step. With 1 mrad of azimuth noise at 10 km, a three-point estimate cannot do better than about 1.1 m/s² RMS (5 %), at a step near 6 s.** The noise, not the coordinates, limits the answer.

What to improve: use all the samples in a window (a least-squares fit) instead of three; model the aircraft's motion and filter (a Kalman filter, `P18.03`, in Phase 18: Estimation, Navigation, Guidance & Autonomy); and track in Cartesian coordinates near the radar's axis. What not to improve: the truncation error, which is already a thousand times smaller than the noise error at the useful steps.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (1)–(3), (6), (7) and (9); `radar_track.png`; the experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) giving the speed, normal acceleration and radius of the turn, with the uncertainty set by the noise study and the step you chose.

### Adversarial self-checks

- **What if the noise doubles?** By (9) the best step grows by 2^{1/4} = 1.19 and the error by √2 = 1.41. Problem 4 checks this with a seeded record.
- **What if the turn tightens?** f⁗ = Rω⁴ grows fast with ω, so the best step shrinks and the error grows: twice the normal acceleration (about 4 g) at the same speed halves R and doubles ω, multiplying f⁗ by 8, the minimum error by √8 = 2.8 and shrinking h* by 8^{1/4} = 1.7.
- **Which assumption fails first in the real world?** Independent noise. Real radar angle errors include slow bias and multipath that a difference does not see as noise; a bias in θ̇ appears directly as a false transverse acceleration 2ṙθ̇.
- **What measurement could falsify the model?** An aircraft that reports its own acceleration (an onboard inertial unit, recorded) flown through the radar's coverage: the radar's normal acceleration should match within the band predicted by (9).
- **Does the result survive a change of route or step?** The two routes agree to O(h²); halving the step multiplies the noise-free error by 1/4, as the table shows.
- **Which decision changes if the uncertainty doubles?** A range-safety rule that flags turns above 2 g would, with a 2.2 m/s² (0.2 g) error, have to flag at 1.8 g to be sure of catching every real 2 g turn, and would raise more false alarms.

## Check yourself

1. Write ė_r and ė_θ in terms of θ̇. Why must each be perpendicular to its unit vector?
2. Name the four terms of the polar acceleration and give the units of each.
3. An aircraft flies a straight line at constant speed past a radar. Which of r̈, θ̈, a_r, a_θ are zero?
4. What are a_t and a_n for an aircraft at constant speed in a level circle of radius 1500 m at 150 m/s?
5. Why is the radius of curvature of a climbing helix larger than its radius R?
6. Why does a central second difference of noisy samples get worse as the step shrinks?
7. By what factor does the best difference step change if the sample noise grows by 16?
8. Why are cylindrical coordinates about the radar a bad choice for an aircraft passing directly overhead?
