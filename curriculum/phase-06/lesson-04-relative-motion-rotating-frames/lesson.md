# P06.04  Relative Motion & Rotating Frames

> Nothing pushes the shell sideways. The ground turns under it, and in the ground's frame that shows up as two extra accelerations you must add by hand.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the transport theorem (d/dt)_inertial = (d/dt)_rotating + ω × and from it the five-term acceleration a_P = a_O + α × r + ω × (ω × r) + 2ω × v_rel + a_rel, and name each term. 2. Write Newton's law in the Earth-fixed rotating frame with Coriolis and centrifugal terms, and convert states between the rotating and inertial frames. 3. Implement and test the integration of a long-range shot in both frames and show that they land in the same place, with the rotating frame's Jacobi integral conserved. 4. Derive the first-order Coriolis deflection on a flat Earth, −ω × (v₀t² + gt³/3), and the Foucault pendulum's exact precession −Ω sin λ, and check both by simulation. 5. Decide when Earth's rotation matters for an aerospace calculation and which of its two effects is which. |
| **Prerequisites** | `P06.01` (turning unit vectors, path coordinates), `P01.15` (rotation matrices, change of basis), `P01.07` (cross products), `P04.05` (gravitation, μ/r²), `P04.06` (the pendulum), `P02.02` (second-order linear ODEs), `P03.06` (adaptive integrators, events). Skill check: you can compute ω × (ω × r) for ω along z and say which way it points. |
| **Engineering question** | A shell leaves the muzzle at 1000 m/s and 45° elevation at 49° N on a spherical Earth without air. How far, and which way, does the Earth's rotation move its impact point when it is fired north, east, south and west, how much of that is the Coriolis term, and does a rotating-frame integration agree with an inertial one? And why does a 67 m pendulum in Paris turn its swing plane once every 31.8 hours? The answers: rotation shifts the impact by 0.6–1.2 km, mostly to the right (the Coriolis cross-range alone is 590, 839, 1064 and 820 m for N, E, S, W, within 4 % of the flat-Earth first-order formula), while the centrifugal term moves it by up to 340 m; the two integrations agree to 10⁻⁴ m; the swing plane turns at −Ω sin λ = −11.33 °/h, to 2.5 × 10⁻⁷ in simulation. |
| **Assumptions and limits** | A spherical Earth of mean radius 6371 km with point-mass gravity, turning at constant Ω about a fixed axis (no oblateness, no orbital motion, no precession). No atmosphere: real shells lose most of their speed to drag, which shortens the flight and the deflection. Launch elevation is measured from the radial (geocentric) vertical, not the plumb line. The Foucault pendulum is linearised (small amplitude, the Ω² terms absorbed into g) and undamped. |
| **Analytical work** | The transport theorem; the five-term acceleration; Newton's law in the rotating frame; frame conversions; the Jacobi integral; the local ENU components of Ω; the first-order deflection on a flat Earth and its impact-time correction (Eötvös); the right-hand rule for deflection in each hemisphere; the Foucault equations, their exact solution from the non-rotating view and the precession rate Ω sin λ. |
| **Computational deliverable** | `starter/main.py`: `five_term_acceleration`, `rotating_frame_acceleration`, `inertial_to_rotating`, `rotating_to_inertial`, `enu_basis`, `launch_state`, `simulate_rotating`, `simulate_inertial`, `first_order_deflection`, `first_order_impact_shift`, `flat_earth_rotating`, `foucault_rhs`, `foucault_closed_form` (`rot_z`, `gravity`, `_ground`, `range_and_cross`, `omega_enu`, `simulate_foucault`, `swing_plane_rate` and `main` are given). Tests prove: the five-term sum equals a brute-force second derivative of an inertial position for an accelerating, spinning frame; its limits; Coriolis does no work; frame conversions round-trip; ENU is orthonormal and right-handed; rotating and inertial integrations land within 10⁻³ m with the Jacobi integral conserved to 10⁻¹⁰; a non-rotating Earth gives no cross-range; the flat-Earth integration matches the first-order formula to 2 %; deflection is right in the north and left in the south and linear in Ω; a northward shot at the equator lands west by ΩT²v_up/3; the spherical Coriolis cross-range is within 6 % above the flat formula; the Foucault closed form solves its equations, matches the integration to 10⁻⁶ m and precesses at −Ω sin λ to 10⁻⁴, with no precession at the equator. |
| **Practical activity** | **computational experiment**: the same shot in two frames; the impact shifts split into their Coriolis and centrifugal parts and compared with the flat-Earth formulas for four azimuths; a latitude sweep from 60° S to 80° N; and a one-hour Foucault run against its closed form. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `rotating_earth.png`, the derivations of the transport theorem, the five-term acceleration, the first-order deflection and the Foucault solution, the experiment tables, and a half-page note answering the engineering question. |
| **Assessment** | `assessments/phase-06/grading-rubric.md`; competencies 1–4. Pass: the five-term acceleration is derived from the transport theorem; all public tests pass; the rotating-frame and inertial-frame impacts are compared with numbers; the note separates the Coriolis and centrifugal effects and states the sign rule. |
| **References** | `@mit-ocw-16-07-dynamics` L7 "Relative motion using translating axes"; L8 "Relative motion using rotating axes": Time derivative of a vector in a rotating frame: Coriolis' theorem; Velocity vector; Acceleration vector; Acceleration in Rotating Cartesian Coordinates (required-free: the transport theorem and the five terms). `@libretexts-cline-variational-principles` Ch. 12 "Non-inertial Reference Frames": 12.3 Rotating Reference Frame; 12.5 Newton's Law of Motion in a Non-Inertial Frame; 12.7 Centrifugal Force; 12.8 Coriolis Force; 12.10 Effective gravitational force near the surface of the Earth; 12.11 Free Motion on the Earth; 12.13 Foucault pendulum (required-free: the Earth applications at a more advanced level). `@mit-ocw-2-003sc-engineering-dynamics` Lecture 3: Motion of Center of Mass; Acceleration in Rotating Ref. Frames (supplementary-free, video). `@mit-ocw-16-61-aerospace-dynamics` lecture 2 "Coriolis "Demystified"" (supplementary-free, handwritten). `@openstax-university-physics-vol1` 6.3 Centripetal Force, subsection Inertial Forces and Noninertial (Accelerated) Frames: The Coriolis Force (supplementary-free, qualitative). |
| **Extensions** | ★ Add quadratic air drag with an exponential atmosphere and redo the four shots: how much of the deflection survives? ★ Measure the elevation from the plumb line (gravity plus centrifugal) instead of the radius, and show that most of the "centrifugal" shift disappears. ★ Integrate an aircraft flying a great circle at 250 m/s for an hour in the rotating frame and find the bank angle the Coriolis term requires to hold the track. |
| **Time** | unmeasured |

## 1. The engineering problem

Every aircraft and rocket is guided in a frame attached to the Earth: runways, radar sites, targets and maps all turn with it. The Earth turns at Ω = 7.292 × 10⁻⁵ rad/s, once per sidereal day. That is slow, but a long flight is long: in the 146 s that a shell fired at 1000 m/s spends in the air, the Earth turns through 0.6°, and the ground at 49° N moves 44 km eastward. Newton's laws hold in an inertial frame, not in the ground's. Write them in the ground's frame and two extra accelerations appear, the **Coriolis** and the **centrifugal** terms, which are not forces at all but the price of describing motion from a turning platform.

Getting them right matters wherever the time of flight or the required accuracy is large: artillery and ballistic-missile aiming, an inertial navigation system (which must subtract the Coriolis term from its accelerometer readings, `P18.07` in Phase 18: Estimation, Navigation, Guidance & Autonomy), the large-scale winds (the Coriolis term turns air moving toward a low-pressure area into a rotating storm), and launch-vehicle trajectories (the Earth's eastward surface speed is free velocity for an eastward launch).

This lesson derives the five-term acceleration, then answers its question two independent ways: by integrating in the rotating Earth-fixed frame with the fictitious terms, and by integrating in an inertial frame where they do not exist and rotating the answer back. Success criteria, written first: the two integrations land in the same place to the solver tolerance; the rotating run conserves its Jacobi integral; the flat-Earth Coriolis integration matches the first-order formula; and the Foucault simulation matches its closed form and precesses at −Ω sin λ.

## 2. Understand and derive

### 2.1 The transport theorem

Let a frame B rotate at angular velocity ω relative to an inertial frame I. A vector A has components in B that change at the rate (dA/dt)_B. But B's own axes turn, and a unit vector fixed in B changes, seen from I, at ω × (unit vector) (the same fact as ė_r = θ̇ e_θ in `P06.01`, now in three dimensions). So

**(dA/dt)_I = (dA/dt)_B + ω × A**   (1)

for any vector A (16.07 L8, "Time derivative of a vector in a rotating frame: Coriolis' theorem").

### 2.2 The five-term acceleration

Let the moving frame's origin O be at R_O, and the point P at r relative to O. The inertial position is R_O + r. Apply (1) once to r: v_P = v_O + v_rel + ω × r, where v_rel = (dr/dt)_B. Apply (1) again to v_rel + ω × r, using (dω/dt)_I = (dω/dt)_B = α:

**a_P = a_O + α × r + ω × (ω × r) + 2 ω × v_rel + a_rel.**   (2)

Five terms: **origin** a_O; **Euler** α × r (the frame's spin changing); **centripetal** ω × (ω × r), pointing toward the axis with size ω²r_⊥; **Coriolis** 2ω × v_rel, from two sources that each contribute ω × v_rel (the relative velocity being rotated, and ω × r changing because r moves); and the **relative** acceleration a_rel. Units: each is m/s².

**On the Earth.** Put O at the Earth's centre (a_O = 0 for our purposes) and ω = Ω e_z constant (α = 0). Newton's law in the inertial frame is a_P = f, the gravitational force per unit mass. Solve (2) for what an observer on the ground sees:

**a_rel = f − 2 Ω × v_rel − Ω × (Ω × r).**   (3)

The ground observer sees two accelerations that no force causes: the **Coriolis acceleration** −2Ω × v_rel, perpendicular to the velocity and to the axis, so it never does work (a test checks v · a_Cor = 0); and the **centrifugal acceleration** −Ω × (Ω × r), away from the axis, a static field that depends only on position.

**Sizes.** An airliner flying north at 250 m/s at 49° N feels a Coriolis acceleration 2ΩV sin λ = **0.0275 m/s² (2.8 mg)**, pointing east (to its right); the centrifugal term at the surface is Ω²R cos λ = **0.0222 m/s²** away from the axis. Both are a few thousandths of g. The Coriolis term matters because it acts in one direction for the whole flight.

### 2.3 Converting between frames

The frames coincide at t = 0 and the rotating one has turned by Ωt since. With R_z the rotation about the spin axis (`P01.15`):

**r_R = R_z(−Ωt) r_I,  v_R = R_z(−Ωt) v_I − Ω × r_R,**   (4)

the second from (1). At launch the gun is carried by the ground, so a shell fired at v₀ relative to the ground has inertial velocity v₀ + Ω × r₀: at 49° N the ground moves east at ΩR cos λ = 305 m/s.

**The Jacobi integral.** Dot (3) with v_rel: the Coriolis term vanishes, f = −∇(−μ/r) and the centrifugal term is −∇(−½|Ω × r|²). So

**J = ½|v_rel|² − μ/|r| − ½|Ω × r|²**   (5)

is conserved in the rotating frame: an energy-like check of the rotating-frame integration.

### 2.4 The local frame and the first-order deflection

At latitude λ use east, north, up (ENU) axes. The Earth's rotation has components **Ω = Ω (0, cos λ, sin λ)**: a horizontal part pointing north and a vertical part. On a flat Earth with uniform g = (0, 0, −g) (the centrifugal term absorbed into g), (3) is r̈ = g − 2Ω × ṙ. Write r = r₀(t) + r₁(t) with r₀ the ordinary parabola, v₀ + g t its velocity, and r₁ of first order in Ω. Then r̈₁ = −2Ω × (v₀ + g t), and integrating twice from rest,

**r₁(t) = −Ω × (v₀ t² + g t³/3).**   (6)

**Reading (6).** For horizontal motion only the vertical part of Ω acts sideways: −2Ω sin λ e_z × v points to the **right of the velocity in the northern hemisphere and to the left in the southern**, with size 2ΩV sin λ. The horizontal part Ω cos λ acts on vertical motion: a rising shell is pushed west, a falling one east, and for a full trajectory the net is a westward shift −ΩT² cos λ v_up/3. At the equator a northward shot (no sideways sin λ term) therefore lands **west** by ΩT²v_up/3.

**The impact point.** At the vacuum flight time T = 2v_up/g, (6) also has a vertical part r₁,z. A shell that is still above (or below) the ground at T lands later (earlier) by δT = r₁,z/v_up, moving along the horizontal velocity meanwhile. To first order the impact shift is r₁,h(T) + v_h δT. The vertical part is the **Eötvös effect**: firing east, the Coriolis acceleration has an upward component, the shell stays up longer and flies further; firing west, the opposite.

For the lesson's shot (1000 m/s at 45°, 49° N, T = 144.2 s on the flat Earth), (6) and its correction give a cross-range of **574.8 m** fired north, **1043.8 m** fired south and **809.3 m** east or west, all to the right, and a down-range change of **+469 m** east and **−469 m** west. North and south differ because the westward drift from the vertical motion is to the left going north and to the right going south.

### 2.5 The Foucault pendulum

A pendulum of length L has small-swing frequency ω₀ = √(g/L). In the horizontal plane (x east, y north), only the vertical component Ω_z = Ω sin λ of the rotation acts sideways on the horizontal swing (the horizontal component changes the string tension, to this order), so (3) linearised is

**ẍ = −ω₀² x + 2Ω_z ẏ,  ÿ = −ω₀² y − 2Ω_z ẋ.**   (7)

Now look from axes that do not turn about the vertical: r = R_z(−Ω_z t) s. In those axes the Coriolis term disappears and a centrifugal one appears, so **s̈ = −(ω₀² + Ω_z²) s**: an isotropic oscillator at W = √(ω₀² + Ω_z²). Released from rest on the floor at (a, 0), it starts with velocity (0, Ω_z a) in the non-turning axes, so s = (a cos Wt, (Ω_z a/W) sin Wt): a very thin fixed ellipse. Seen from the floor,

**r(t) = R_z(−Ω_z t) s(t):  the swing plane turns clockwise (in the north) at exactly Ω sin λ.**   (8)

For the 67 m pendulum at 48.85° N (the latitude of Paris; the length of the famous Panthéon pendulum): period 16.42 s, precession **−5.4909 × 10⁻⁵ rad/s = −11.326 °/h, one turn in 31.79 h**; the ellipse's minor semi-axis is only Ω_z a/W = 0.43 mm for a 3 m swing. At the equator Ω_z = 0 and there is no precession; at the poles it turns once a sidereal day.

### 2.6 Where the model stops

- **No air.** A real shell at 1000 m/s loses most of its speed to drag; the deflection then follows the real (shorter, steeper) flight. The physics of (3) is unchanged; only the trajectory is.
- **Gravity and the vertical.** Point-mass gravity on a sphere is a model. The real "vertical" is the plumb line, set by gravity plus the centrifugal term (and the Earth's oblateness); a gun is laid relative to it. Here the elevation is from the radius, which makes the centrifugal term look like a separate effect (Extension 2).
- **Small-angle Foucault.** Real pendulums have elliptical swings started by imperfect release, and an elliptical swing precesses on its own (the Airy precession) unless the amplitude is small or the bob is driven carefully.

## 3. Build the computational model

Open `starter/main.py`. `rot_z`, `gravity`, `_ground`, `range_and_cross`, `omega_enu`, `simulate_foucault`, `swing_plane_rate` and `main` are given.

### Step 1: the five terms

`five_term_acceleration(a_O, alpha, omega, r, v_rel, a_rel)` returns the sum (2) and a dict of the five terms. `rotating_frame_acceleration(r, v, f, omega)` is (3).

### Step 2: frames and the launch

`inertial_to_rotating` and `rotating_to_inertial` are (4) and its inverse. `enu_basis(lat, lon)` returns the rows E, N, U in Earth-fixed components; `launch_state(lat, lon, speed, elevation, azimuth)` returns r₀ = R U and v₀ along (cos el (sin az E + cos az N) + sin el U).

### Step 3: two integrations

`simulate_rotating` integrates (3) with point-mass gravity (DOP853, rtol 10⁻¹¹) to the ground event |r| = R; its `coriolis` and `centrifugal` switches exist only for the experiment that separates the two terms. `simulate_inertial` converts the launch state with (4), integrates r̈ = f with nothing else, and rotates the impact point back by the Earth's turn during the flight:

```python
rI0, vI0 = rotating_to_inertial(r0, v0, 0.0, Omega)       # the gun moves with the ground
...                                                       # integrate r'' = -mu r / |r|^3 to |r| = R
rR, vR = inertial_to_rotating(r_hit, v_hit, t_hit, Omega) # where is that on the turned Earth?
```

### Step 4: flat Earth and the closed form

`flat_earth_rotating(v0_enu, lat)` integrates r̈ = g − 2Ω × ṙ in ENU to z = 0. `first_order_deflection` is (6); `first_order_impact_shift` adds the impact-time correction of 2.4.

### Step 5: Foucault

`foucault_rhs` is (7); `foucault_closed_form` is (8).

Running `main()` prints the airliner's terms, the two-frame comparison, the split of the impact shift for four azimuths, the latitude sweep and the Foucault run, and writes `rotating_earth.png` (Coriolis cross-range of a northward shot against latitude with the flat formula; the Foucault bob's first 150 s; the swing plane turning over an hour).

### Verify

Run `python -m pytest curriculum/phase-06/lesson-04-relative-motion-rotating-frames`. The tests:

1. **Five terms** against a brute-force second difference of R_O(t) + R_z(θ(t)) r(t) for a frame that accelerates and spins up: agreement to 10⁻⁵.
2. **Limits**: no rotation gives a_O + a_rel; a point fixed in a turning frame has only ω²r toward the axis.
3. **Coriolis** is perpendicular to v and to the spin axis.
4. **Conversions** round-trip; at t = 0 the velocities differ by Ω × r.
5. **ENU** is orthonormal with E × N = U; `launch_state` gives the speed, elevation and azimuth asked for.
6. **Two frames**: the impact points agree within 10⁻³ m and the flight times to 10⁻⁹; the Jacobi integral (5) holds to 10⁻¹⁰.
7. **No rotation**: no cross-range.
8. **Flat Earth**: integration against (6) with its correction to 2 % for N, E and S shots.
9. **Hemispheres**: right in the north, left in the south, size 2ΩV sin|λ| t²/2; linear in Ω.
10. **Equator**: a northward shot lands west by ΩT²v_up/3 exactly.
11. **Sphere against flat**: the spherical Coriolis cross-range is above the flat value and within 6 % of it.
12. **Foucault**: (8) satisfies (7) and the release condition; integration matches it to 10⁻⁶ m; the precession is −Ω sin λ to 10⁻⁴; none at the equator; +Ω at the south pole; no rotation, no precession.
13. `main` writes the figure.

### Compare

The independent check is the inertial integration. For the four shots:

| azimuth | flight (s) | range (km) | \|rotating − inertial\| (m) | flight-time difference (s) |
|---|---|---|---|---|
| north | 146.170 | 102.621 | 1.6 × 10⁻⁵ | 2.0 × 10⁻⁸ |
| east | 147.196 | 103.273 | 5.3 × 10⁻⁵ | 7.1 × 10⁻⁸ |
| south | 146.170 | 102.975 | 7.4 × 10⁻⁵ | 9.9 × 10⁻⁸ |
| west | 145.157 | 102.328 | 7.2 × 10⁻⁵ | 9.9 × 10⁻⁸ |

Two completely different calculations, one with the fictitious terms and one without, land within a tenth of a millimetre after 100 km. A sign error in the Coriolis term, a factor of 2 missing, or a forgotten Ω × r₀ in the launch velocity would each show up here as hundreds of metres. On a non-rotating Earth the same shot ranges 102.649 km in 145.950 s; on a flat Earth with uniform g, 101.972 km in 144.210 s.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** All results are simulations of a vacuum shot on a model Earth, not firing data. Write each prediction down before running it.

1. **Two frames.** Run the four shots both ways. Predict the size of the disagreement. *Finding:* the table above, 10⁻⁵–10⁻⁴ m.
2. **Splitting the rotation's effect.** For each azimuth, compare the full rotating run with runs that switch off the Coriolis term or the centrifugal term, with a run on a non-rotating Earth, and with the flat-Earth formulas. Predict the sign of each cross-range from the right-hand rule. *Finding* (m; cross-range positive to the right of the line of fire):

| azimuth | Coriolis: down | Coriolis: cross | centrifugal: down | centrifugal: cross | flat 1st order: down | cross | flat sim: down | cross | no rotation at all: down | cross |
|---|---|---|---|---|---|---|---|---|---|---|
| north | −0.4 | 590.3 | −27.6 | 0.8 | 0.0 | 574.8 | −0.3 | 574.8 | −28.0 | 590.3 |
| east | 468.0 | 839.1 | 156.7 | 184.3 | 469.0 | 809.3 | 464.6 | 818.6 | 624.3 | 1018.3 |
| south | −11.8 | 1063.5 | 338.4 | 4.2 | 0.0 | 1043.8 | −11.5 | 1043.7 | 326.6 | 1063.5 |
| west | −477.1 | 819.9 | 155.7 | −174.3 | −469.0 | 809.3 | −473.2 | 800.0 | −320.8 | 640.7 |

Every Coriolis cross-range is to the right; it is largest going south and smallest going north, as 2.4 predicts. The spherical Coriolis cross-ranges are 1.3–3.7 % above the flat-Earth first-order formula; most of that difference is the 1.4 % longer flight on the sphere (gravity weakens with height and the ground curves away), which enters as T². The centrifugal term is a static pull of 0.0168 m/s² toward the equator and 0.0146 m/s² outward at 49° N; over 146 s the horizontal part alone moves the shell about 179 m south, which is why it shifts the east and west shots sideways and the north and south shots along range. A gun laid by a plumb line would include this field in its "gravity" (Extension 2).
3. **Latitude.** The Coriolis cross-range of a northward shot from 60° S to 80° N. Predict where it changes sign. *Finding:* −1128.0, −861.3, −360.6, 236.8, 590.3, 767.4 and 1013.4 m at −60°, −30°, 0°, 30°, 49°, 60° and 80°, against −1107.4, −845.7, −357.4, 226.6, 574.8, 749.9 and 994.0 m from the flat formula. It changes sign not at the equator but between 0° and 30° N, because the westward drift of the vertical motion (Ω cos λ) is still there when the sideways term (Ω sin λ) vanishes.
4. **Foucault.** A 67 m pendulum at 48.85° N, 3 m amplitude, released from rest, integrated for one hour. Predict the swing plane's turn. *Finding:* −11.326° in one hour; the fitted rate −5.490885 × 10⁻⁵ rad/s against −5.490886 × 10⁻⁵ predicted (2.5 × 10⁻⁷ relative); the integration matches (8) to 7.2 × 10⁻⁸ m.

Raw data format: CSV tables `azimuth, t_flight_s, range_km, diff_rot_inertial_m` (1); `azimuth, cor_down_m, cor_cross_m, cen_down_m, cen_cross_m, flat1_down_m, flat1_cross_m, flatsim_down_m, flatsim_cross_m, norot_down_m, norot_cross_m` (2); `lat_deg, cross_sphere_m, cross_flat_m` (3); `t_s, x_m, y_m` (4, every 50th sample); plus `rotating_earth.png`.

## 5. Validate and improve

| quantity | analytic (§2) | simulation | difference | cause |
|---|---|---|---|---|
| impact point, rotating vs inertial (m) | identical | 1.6–7.4 × 10⁻⁵ apart | 10⁻⁴ | solver tolerance |
| Jacobi integral (5) | constant | constant to 10⁻¹⁰ (test) | 10⁻¹⁰ | solver tolerance |
| Coriolis cross-range, north, 49° N (m) | 574.8 (6) flat | 574.8 flat sim; 590.3 sphere | 0; 2.7 % | first order; longer spherical flight |
| equator, northward: west shift | −ΩT²v_up/3 = −357.4 m | sphere −360.6 m | 0.9 % | sphere |
| Foucault precession (rad/s) | −5.490886 × 10⁻⁵ (8) | −5.490885 × 10⁻⁵ | 2.5 × 10⁻⁷ | turning-point fit |
| Foucault trajectory (m) | (8) | max error 7.2 × 10⁻⁸ | 10⁻⁷ | solver tolerance |

The engineering answer: **Earth's rotation moves the impact of a 100 km vacuum shot from 49° N by 0.6–1.2 km. The Coriolis term moves it to the right of the line of fire by 590 m (north) to 1064 m (south), and along range by ±470 m east and west, within 4 % of the flat-Earth first-order formula; the centrifugal term, a static field, moves it by up to 340 m and belongs with gravity once the gun is laid by a plumb line. Rotating-frame and inertial-frame integrations agree to 10⁻⁴ m. The Foucault pendulum's plane turns at Ω sin λ: 11.33 °/h, one turn in 31.8 h in Paris.**

What to improve: drag (Extension 1), which changes the trajectory and so the size of the deflection; the plumb-line vertical (Extension 2); an oblate Earth for longer ranges. What not to improve: the solver, which already agrees with the independent frame to 10⁻⁴ m.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (1), (2), (3), (5), (6) and (8); `rotating_earth.png`; the four experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) answering the engineering question and stating the sign rule.

### Adversarial self-checks

- **What if the muzzle speed doubles?** The flight time doubles, the range quadruples, and (6) scales as v₀T² ~ v₀³: the deflection grows about eightfold (first order). At 1600 m/s the flat formula gives 2.35 km of cross-range on a 261 km flat-Earth range, where Earth curvature is no longer small.
- **What if Ω → 0?** Both fictitious terms vanish and the cross-range goes to zero (test 7); the first-order deflection is linear in Ω (test 9).
- **Which assumption fails first in the real world?** No air: a real shell's flight is shorter and lower, and its deflection is set by its real flight time.
- **What measurement could falsify the model?** A Foucault pendulum with a long, carefully started swing: if the swing plane did not turn at Ω sin λ (once friction and Airy precession are accounted for), equation (3) would be wrong. Historically it did, which is why the pendulum was the first laboratory proof of the Earth's rotation.
- **Does the result survive a change of frame?** Yes: that is the main check, to 10⁻⁴ m.
- **Which decision changes if the uncertainty doubles?** A gunner who ignores the 590 m cross-range at 100 km misses by more than a target's size; if the muzzle speed is uncertain by 2 % the range error (about 4 km) dwarfs the Coriolis correction, and correcting for Coriolis first is pointless.

## Check yourself

1. State the transport theorem and say what each side's derivative means.
2. Name the five terms of the acceleration of a point seen from a moving frame.
3. Why does the Coriolis acceleration do no work?
4. Which way is an aircraft flying south in the southern hemisphere deflected, relative to its track?
5. Why is a northward shot at the equator deflected at all, and which way?
6. Why must the inertial integration start with v₀ + Ω × r₀ rather than v₀?
7. How long does a Foucault pendulum at 30° N take to turn its swing plane once?
8. What is the Eötvös effect on a shell fired east?
