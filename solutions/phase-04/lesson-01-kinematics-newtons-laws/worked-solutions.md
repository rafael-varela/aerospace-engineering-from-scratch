# P04.01  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for **s = 15**; numerical results are from `solution.py` (ρ = 1.20 kg/m³, g = 9.80665 m/s², ping-pong k = 0.1396 m⁻¹, v_t = 8.38 m/s) and given to the figures shown. SYNTHETIC data are generated from the model and are never measurements.

## 1. A toss by hand and by code (s = 15: v₀ = 6.5 m/s, θ = 40°)

**Assumptions.** Point mass, launch and landing at the same height, still air, no spin; (a) no drag.

(a) T = 2v₀ sin θ/g = 2 × 6.5 × 0.6428/9.807 = **0.8521 s**; R = v₀² sin 80°/g = 42.25 × 0.9848/9.807 = **4.243 m**; H = (6.5 × 0.6428)²/(2 × 9.807) = **0.890 m**. `vacuum_state(0.8521, 6.5, 40)` = (4.243, −4 × 10⁻¹⁶, 4.979, −4.178): z = 0 and x = R, exact.

(b) `simulate_projectile(6.5, 40.0, 0.1396)`: **range 3.005 m (70.8 % of vacuum)**, time of flight **0.771 s** (90.5 %), apex **0.732 m** (82.3 %), impact speed **4.67 m/s** (71.8 % of v₀), impact angle **48.6°** below the horizontal (launch 40°).

(c) The range changed most (−29 %), then the impact speed (−28 %, which is 48 % of the kinetic energy lost to drag), the apex (−18 %), the time (−9 %). The time is shorter because the ball never gets as high: the vertical speed is eroded on the way up, so the apex (0.73 m) is reached sooner than the vacuum apex (0.89 m) would be, and the descent from a lower height, even slowed by drag, takes less time than the vacuum descent from the higher one; the lower apex wins over the slower fall.

**Check.** Dimensional: v₀²/g = 4.31 m sets the scale of R and 2v₀/g = 1.33 s of T; both results are a fraction below 1 of these. Limiting: `simulate_projectile(6.5, 40, 0.0)` returns (a) to 10⁻⁸.

**Common mistakes.** Using sin θ instead of sin 2θ in the range; reporting the impact angle above the horizontal; expecting the time of flight to grow with drag.

## 2. The closed forms and their expansions

(a) Newton's second law along the downward axis with drag opposing the motion: m dv/dt = mg − ½ρC_D A v² = mg(1 − v²/v_t²) with v_t² = 2mg/(ρC_D A). Separate: dv/(1 − v²/v_t²) = g dt. The left side integrates to v_t artanh(v/v_t) (differentiate to check: v_t · (1/v_t)/(1 − v²/v_t²)). So v_t artanh(v/v_t) = gt + C, and **"from rest" sets C = 0** (v = 0 at t = 0). Inverting, **v = v_t tanh(gt/v_t)**. Then s = ∫₀ᵗ v dt′ = v_t ∫₀ᵗ tanh(gt′/v_t) dt′ = (v_t²/g) ln cosh(gt/v_t), with s(0) = 0 (release point as origin): **s = (v_t²/g) ln cosh(gt/v_t)**.

(b) Let u = gt/v_t and a = gh/v_t². Then h = (v_t²/g) ln cosh u gives a = u²/2 − u⁴/12 + O(u⁶). In vacuum u₀² = 2a (from a = u₀²/2). Write u = u₀(1 + ε): u²/2 = a(1 + 2ε + ε²), u⁴/12 = (4a²/12)(1 + 4ε + …) = a²/3 + O(a²ε). To first order in ε and a: a = a + 2aε − a²/3, so **ε = a/6** and t = t₀(1 + gh/(6v_t²)). At h = 2 m: a = 0.2793, ε = 0.04654, t ≈ 0.63866 × 1.04654 = **0.6684 s**; `vertical_drop_time(2.0, 8.38)` = **0.6688 s**, 0.06 % higher, from the u⁶ term that a = 0.28 is large enough to feel.

(c) dv/dt = −kv² ⇒ −(1/v²) dv = k dt ⇒ d(1/v)/dt = k ⇒ 1/v = 1/v₀ + kt ⇒ **v = v₀/(1 + kv₀t)**. x = ∫₀ᵗ v dt′ = ln(1 + kv₀t)/k. Setting v = v₀/2: 1 + kv₀t = 2, so x = **ln 2/k**, independent of v₀ (4.96 m for the ping-pong ball).

(d) From ẍ = −k|v|ẋ: ẋ(t) = ẋ₀ − k∫₀ᵗ|v|ẋ dt′. While ẋ > 0 the integrand is positive, so **ẋ is strictly decreasing**; it cannot reach zero in finite time (near ẋ = 0 the equation is ẍ ≈ −k|ż|ẋ, an exponential decay), so ẋ stays positive throughout.

For the impact angle, use the angle ψ of the velocity above the horizontal as the independent variable. With u = ẋ = v cos ψ, the component of (1) normal to the velocity gives v ψ̇ = −g cos ψ (drag has no normal component), so ψ falls monotonically from θ at launch to −α at impact. Then du/dψ = u̇/ψ̇ = (−kvu)/(−g cos ψ/v) = kv²u/(g cos ψ) = ku³/(g cos³ψ) > 0: **u is an increasing function of ψ**, i.e. the horizontal speed is larger at every point of the ascent than at the point of the descent with the mirror angle, u(ψ) > u(−ψ) for ψ > 0. The height change is dz = ż dt = v sin ψ · (−v dψ/(g cos ψ)) = −(u²/g) tan ψ sec²ψ dψ, so equal launch and landing heights require

  ∫ from −α to θ of u(ψ)² tan ψ sec²ψ dψ = 0.

Suppose α ≤ θ. Split the interval into [−α, α] and [α, θ]. On [α, θ] the integrand is non-negative. On [−α, α] pair ψ with −ψ: tan ψ sec²ψ is odd, so the pair contributes [u(ψ)² − u(−ψ)²] tan ψ sec²ψ > 0 for 0 < ψ ≤ α. The integral is then strictly positive, a contradiction. Hence **α > θ**: the descent is steeper than the ascent. (With k = 0, u is constant, the pairs cancel and α = θ, the vacuum symmetry.) The simulation agrees: 48.6° against 40° in Problem 1, 52.8° against 45° in the lesson.

**Check.** (a) and (c) differentiate back to their ODEs; (b) at v_t = 40 m/s, h = 1 m the expansion agrees with the code to 10⁻⁵ (a test).

**Common mistakes.** Forgetting the "from rest" constant; expanding ln cosh only to u² (which gives no correction at all); treating (6) as valid when gravity acts (it is not: |v| includes ż).

## 3. A supply drop under a parachute (s = 15: m = 21.5 kg)

**Assumptions.** Quadratic drag with a constant C_D on the canopy's projected area; the parachute is open from release and the system falls as one point mass; still air.

(a) A = π(1.2)²/4 = 1.131 m²; k = 1.20 × 1.5 × 1.131/(2 × 21.5) = **0.04734 m⁻¹**; v_t = √(9.807/0.04734) = **14.39 m/s**; t = (v_t/g) arccosh(exp(gh/v_t²)) with gh/v_t² = 1.4206: **3.080 s**; impact speed v_t tanh(gt/v_t) = **13.97 m/s**, **97 % of v_t**. Vacuum: t₀ = 2.474 s, impact 24.26 m/s. The parachute adds 0.61 s (25 %) and removes 42 % of the impact speed (67 % of the impact energy).

(b) A = π(0.4)²/4 = 0.1257 m², C_D = 1.0: k = 0.003507 m⁻¹, v_t = **52.9 m/s**, t = **2.517 s** (+1.8 % on vacuum), impact **23.0 m/s** (44 % of v_t). The time is the integral of 1/v over the whole fall and most of the fall happens before the speed is large (the first half of the height takes 71 % of the time in vacuum), where drag is weak; the impact speed is the end state, where drag has had the whole fall to act. The parachute's time effect is modest (25 %) while its speed effect (−42 %) is what matters for the cargo.

(c) `simulate_projectile(10, 0, 0.04734, z0=30, max_time=100)`: lands **16.15 m** downrange after **3.179 s**, with impact velocity (**1.76, −13.81**) m/s, speed 13.92 m/s. Vacuum: 10 × 2.474 = **24.74 m**. Pure horizontal decay over 3.179 s: ln(1 + 0.04734 × 10 × 3.179)/0.04734 = **19.40 m**. Ordering: vacuum (no drag on x) > pure decay (drag on x computed with |v| = ẋ only) > coupled (drag on x computed with the full |v|, which the 14 m/s vertical speed dominates, so the horizontal speed is killed faster: from 10 to 1.76 m/s). The flight takes 0.1 s longer than the pure drop because the drag on the vertical motion is larger when |v| includes ẋ.

(d) In the air frame nothing changes: the package lands 16.15 m east of the release point *in the air*. The air moves 4 m/s south for 3.179 s, 12.7 m, so in the ground frame the package lands 16.2 m east and **12.7 m south** of the release point (which is also where the drone's ground track would have put the still-air landing point moved by the wind). Relative motion: v_ground = v_air + v_wind, integrated over the same flight time.

**Check.** Energy at impact in (a): ½ × 21.5 × 13.97² = 2.1 kJ against 6.3 kJ in vacuum. Dimensional: 1/k = 21 m in (a), the length scale of the drag, comparable to the 30 m fall: strong drag, as the 97 % of v_t confirms.

**Common mistakes.** Using the package's area instead of the canopy's; taking the impact speed as v_t (it is 97 % of it: close, but the question asks for the fraction); adding the wind's displacement in the wrong direction (from the north means toward the south).

## 4. Designing the drop experiment (s = 15)

(a) With u_meas = 0.41/f only (Type A ignored), E_n = (t − t₀)/u:

| h (m) | excess (ms) | E_n at 30 fps (u = 13.6 ms) | 60 fps (6.8 ms) | 120 fps (3.4 ms) |
|---|---|---|---|---|
| 0.6 | 4.9 | 0.36 | 0.72 | 1.44 |
| 0.9 | 9.0 | 0.66 | 1.33 | 2.65 |
| 1.2 | 13.9 | 1.02 | 2.05 | 4.09 |
| 1.5 | 19.5 | 1.43 | 2.86 | 5.73 |
| 1.8 | 25.7 | 1.89 | 3.77 | 7.55 |

At 30 fps the vacuum model is never rejected even in the most optimistic case; at 60 fps from 1.2 m up; at 120 fps from 0.9 m up.

(b) SYNTHETIC drops, seeds 0–19 (`make_synthetic_drop_data`, labelled synthetic; design only):

| fps | n | mean fitted v_t (m/s) | sd over seeds | mean reported u_v | u_CD/C_D = 2u_v/v_t |
|---|---|---|---|---|---|
| 30 | 5 | 8.52 | 0.65 | 1.76 | 42 % |
| 30 | 10 | 8.36 | 0.28 | 1.59 | 38 % |
| 60 | 5 | 8.36 | 0.32 | 0.84 | 20 % |
| 60 | 10 | 8.36 | 0.21 | 0.81 | 19 % |
| 120 | 5 | 8.33 | 0.15 | 0.41 | 10 % |
| 120 | 10 | 8.37 | 0.10 | 0.40 | 10 % |

The reported u_v is *larger* than the scatter over seeds by a factor 2–6, because the frame term is entered as a fully systematic Type B uncertainty at every height while the generator's quantisation phase is random per drop and partly averages out: the budget is conservative, which is the right side to err on. The means are unbiased (8.38 true).

(c) **60 fps with ten drops per height** gives ±19 % on C_D, the cheapest combination meeting ±20 %. Doubling the drops from 5 to 10 at 60 fps moves u_v from 0.84 to 0.81 m/s (the Type A part shrinks, the frame term does not), while doubling the frame rate halves u_v outright; frame rate first, drops second, and 120 fps if the phone has it.

(d) With a +15 ms release bias at 60 fps, n = 10: the fitted v_t drops to **6.25 ± 0.38 m/s** (true 8.38): the bias is absorbed as extra drag, because a constant late start looks like a larger 1/v_t² term, especially at the low heights where the genuine excess is small; C_D would come out as 0.90, 80 % high, with a reported uncertainty that does not cover it. The bean bag reveals it: its fit should return a very large v_t (unbiased SYNTHETIC runs give a median of 17.8 m/s with u ≈ 9 m/s, i.e. "undetermined, large"), but with the bias it returns **8.2 ± 0.25 m/s**, a terminal speed for a 100 g bag that no drag model allows. That is what a control object is for.

**Check.** The 1/v_t² scaling of (5): from 30 to 60 fps the frame term halves and u_v halves (1.76 → 0.84), as the linearised fit predicts.

**Common mistakes.** Treating the frame term as Type A and dividing it by √n; reading the SYNTHETIC scatter as the real experiment's uncertainty (the real one adds release scatter and the frame-rule ambiguity); reporting C_D without an uncertainty.

## 5. A UAV survey square (s = 15: W = 8 m/s from 040°)

(a) `heading_for_track(25, track, 8, 40)`:

| track | heading | ground speed (m/s) | leg time (s) |
|---|---|---|---|
| 000° | 011.9° | 18.34 | 109.1 |
| 090° | 075.8° | 19.10 | 104.7 |
| 180° | 168.1° | 30.59 | 65.4 |
| 270° | 284.2° | 29.38 | 68.1 |

Total **347.3 s** against 320 s in still air (+8.5 %).

(b) Along the wind: t = d/(V − W) + d/(V + W) = d(V + W + V − W)/(V² − W²) = 2dV/(V² − W²) = (2d/V) · **1/(1 − W²/V²)**. Across the wind: the heading must cancel the full crosswind, leaving ground speed √(V² − W²) both ways, t = 2d/√(V² − W²) = (2d/V) · **1/√(1 − W²/V²)**. Both factors exceed 1 for 0 < W < V, and 1/(1 − x) > 1/√(1 − x) for 0 < x < 1, so **the along-wind round trip is worse**: here 1.114 against 1.056.

(c) The 270° leg sees a crosswind component W sin(40° − 270°) = W sin 130° = 0.766W. It cannot be held when 0.766W > 25 m/s, i.e. **W > 32.6 m/s**; `heading_for_track` raises `ValueError` ("crosswind component exceeds the true airspeed"). Well before that, at W = 25 m/s the 000° leg (ground speed V cos WCA − W cos 40°) has a ground speed near zero and the leg time diverges.

(d) 1200/347.3 = 3.46 → **3 squares** in the wind; 1200/320 = 3.75 → **3** in still air: the wind costs the margin, not a lap, with this battery.

**Check.** Each heading substituted into `wind_triangle` returns the wanted track to 10⁻⁸° (a test). The slow legs (000°, 090°) are the ones with a headwind component from 040°, as they should be.

**Common mistakes.** Adding the wind to the track direction rather than the heading; using the wind's "to" direction; computing the round-trip time as 2d/V_mean (the harmonic mean, not the arithmetic mean, of the ground speeds governs).

## 6. Where the model stops

(a) **No spin** fails: backspin gives a Magnus force upward, of the order of the drag at a few hundred rpm; the ball flies *longer* than predicted (the model under-predicts the range). Measurement: mark the ball with a pen and read the spin rate from the video; or compare tosses with the ball resting on the fingertips (little spin) against flicked tosses.
(b) **Constant C_D** fails: at Re = 70 × 0.0427/1.5 × 10⁻⁵ = 2 × 10⁵ a smooth sphere is at the drag crisis, and a dimpled golf ball has C_D ≈ 0.25 at that speed, half the sphere value; with C_D = 0.5 the model under-predicts the carry badly, and the ball's backspin lift (a) makes it worse. Measurement: a launch monitor's measured carry against the simulation at the measured launch speed.
(c) **Sphere with C_D = 0.5, Re > 10³** fails: a feather is not a sphere, its Reynolds number is of order 10²–10³ where the drag is partly viscous (linear in v) and its C_D on any sensible area is far above 0.5; the model with sphere parameters predicts a fall time far too short. Measurement: the fall time (a second or more from 2 m) against the 0.64 s prediction; a terminal speed well below 2 m/s.
(d) **Still air** fails: a draught of 0.5 m/s is 6 % of v_t and gives the ball a horizontal drift of several centimetres in a 0.67 s fall, and an up- or down-draught of the same size changes the fall time by about 6 % of the drag excess (2 ms), which is below the frame term; the visible symptom is a landing point that wanders. Measurement: the horizontal landing scatter against the still-air repeatability; close the door and repeat.
(e) **ρ** fails: at 2500 m the ISA density is 0.957 kg/m³, 20 % below the 1.20 in the code, so k is 20 % too large in the code and the predicted drag excess 20 % too large; the fitted C_D (which absorbs ρ) would come out 20 % low. Measurement: barometric pressure and temperature on the day; the fit should use the local ρ.
(f) A steel bearing: m = 7800 × π(0.008)³/6 = 2.09 g, k = 1.20 × 0.5 × π(0.004)²/(2 × 0.00209) = 0.0072 m⁻¹, v_t = 36.9 m/s, drag excess at 2 m = **1.5 ms**, against a frame term of 3.4 ms even at 120 fps: the camera cannot see it, and modelling the drag changes the prediction by less than the uncertainty. **Do not bother**: state the 1.5 ms as the model's uncertainty and use free fall. (The ball bearing is, however, hard and small; it is outside the project's permitted objects.)

## Check-yourself answers

1. ΣF = ma with F = W + F_D: the weight mg downward, exerted by the Earth; the drag ½ρC_D A v² opposite to v, exerted by the air. Nothing else acts after release.
2. k = ρC_D A/(2m), units 1/m. 1/k is the drag length: over a distance ln 2/k = 0.69/k the speed of a body moving through air with no other force halves, whatever its speed.
3. v_t = √(g/k) = √(9.807/0.14) = 8.37 m/s. A vacuum drop from 2 m reaches √(2gh) = 6.26 m/s, 75 % of v_t: the ball is deep in the drag regime.
4. v = v_t tanh(gt/v_t), s = (v_t²/g) ln cosh(gt/v_t). For gt ≪ v_t, tanh u ≈ u and ln cosh u ≈ u²/2, so v ≈ gt and s ≈ ½gt².
5. t ≈ t₀(1 + gh/(6v_t²)): the correction grows with h and shrinks as 1/v_t², so a heavier or smaller body (larger v_t) shows less of it.
6. ẋ decreases monotonically (drag always has a component against it), so at impact the horizontal speed is smaller than at launch while the ball must still come down; the velocity vector is therefore tilted more steeply. A lower launch angle keeps the ball in the fast, low part of the flight for less time at high |v| and trades apex height for horizontal speed; the optimum moves from 45° to 42.6° at 6 m/s and 37° at 15 m/s.
7. v_ground = (0, 50) + (−20, 0) = (−20, 50): ground speed 53.9 m/s, track 338.2°. To make good 000°: sin WCA = (20/50) sin 90° = 0.4, WCA = 23.6°, heading 023.6°, ground speed 50 cos 23.6° − 20 cos 90° = 45.8 m/s.
8. No: the drag excess at 1.2 m is 13.9 ms and the frame term at 30 fps is 13.6 ms, so E_n ≈ 1.0 before any Type A scatter is added; the two models cannot be told apart at that height and frame rate.
