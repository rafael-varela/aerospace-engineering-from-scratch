# P06.01  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe script that calls it, given to the figures shown. All noisy data are synthetic.

## 1. One radar plot, three descriptions (s = 15: r = 13.0 km)

**Assumptions.** Flat Earth, level flight, the radar's numbers exact.

(a) v_r = ṙ = **−150.0 m/s**; v_θ = rθ̇ = 13 000 × 0.010 = **130.0 m/s**; speed **198.49 m/s**.

(b) The four terms: r̈ = 2.0, rθ̇² = 1.3, rθ̈ = −1.3, 2ṙθ̇ = −3.0 m/s². So a_r = 2.0 − 1.3 = **0.70 m/s²** and a_θ = −1.3 − 3.0 = **−4.30 m/s²**. The 2ṙθ̇ term dominates: forget it and a_θ would be −1.3 instead of −4.3 m/s².

(c) At θ = 30°: v = **(−194.90, 37.58) m/s**, a = **(2.756, −3.374) m/s²** (east, north).

(d) a_t = **−3.345 m/s²**, a_n = **2.791 m/s²** (a_n/g = **0.285**), ρ = **14 117 m**; check √(3.345² + 2.791²) = 4.357 m/s² = |a| ✓. The aircraft is **slowing down** (a_t < 0) and **turning left** (counterclockwise): (v × a)_z = (−194.90)(−3.374) − (37.58)(2.756) = +554 m²/s³ > 0.

**Independent check.** Polar and path descriptions give the same |a| (4.357 m/s²), as they must: |a|² = a_r² + a_θ² = 0.49 + 18.49 = 18.98.

**Common mistakes.** Using θ in degrees in rθ̇ (θ̇ is already in rad/s); dropping 2ṙθ̇; reading the sign of a_θ as the turn direction (a_θ is about the radar, not about the path).

## 2. Derivations

(a) e_r = (cos θ, sin θ). By the chain rule, ė_r = (−sin θ, cos θ) θ̇ = θ̇ e_θ, and ė_θ = (−cos θ, −sin θ) θ̇ = −θ̇ e_r. Product rule on **r** = r e_r: v = ṙ e_r + r ė_r = ṙ e_r + rθ̇ e_θ. Product rule again: a = r̈ e_r + ṙ θ̇ e_θ + (ṙθ̇ + rθ̈) e_θ + rθ̇ ė_θ = (r̈ − rθ̇²) e_r + (rθ̈ + 2ṙθ̇) e_θ.

(b) a_θ = rθ̈ + 2ṙθ̇ = (1/r) d(r²θ̇)/dt. On a straight line at constant velocity a = 0, so a_θ = 0 and **r²θ̇ is constant** (it equals V d, twice the rate at which the radius vector sweeps area). Put the track on y = d, x = Vt: θ = atan(Vt/d) (measured from the perpendicular), θ̇ = Vd/(d² + V²t²), largest at t = 0: **θ̇_max = V/d**, at closest approach. θ̈ = −2V³dt/(d² + V²t²)²; dθ̈/dt = 0 gives d² + V²t² = 4V²t², i.e. **Vt = ±d/√3**, where |θ̈| = 2V³d(d/(√3V))/((4/3)d²)² = **(3√3/8)(V/d)²**. It peaks when the aircraft is d/√3 before and after closest approach, at 30° from the perpendicular. For V = 200 m/s: 0.0010 rad/s² at d = 5 km and 10.39 rad/s² at d = 50 m (probe, to four figures).

(c) v = (−Rω sin ωt, Rω cos ωt, w), a = (−Rω² cos ωt, −Rω² sin ωt, 0). v × a = (wRω² sin ωt, −wRω² cos ωt, R²ω³), so |v × a| = Rω²√(w² + R²ω²) and |v|³ = (R²ω² + w²)^{3/2}. Then κ = Rω²/(R²ω² + w²) = **R/(R² + b²)** with b = w/ω. As w → 0, κ → 1/R (a circle). As R → 0, κ → 0 (a vertical straight line). For the lesson's helix: 4.9875 × 10⁻⁴ m⁻¹, ρ = 2005.00 m.

**Common mistakes.** Writing ė_r = θ̇ e_r (the derivative of a unit vector is perpendicular to it); forgetting that |θ̈| is zero at closest approach, where θ̇ is largest.

## 3. A tracking antenna for a sounding rocket

**Assumptions.** Flat ground; the rocket a particle on a vertical line; constant acceleration from rest.

(a) z = at²/2, R = √(D² + z²), e = atan(z/D). Differentiating, ė = (ż/D)/(1 + z²/D²) = **D ż/(D² + z²)**.

(b) ė(t) = D a t/(D² + a²t⁴/4). Setting dė/dt = 0: (D² + a²t⁴/4) − a²t⁴ = 0, so D² = 3a²t⁴/4 and **z = D/√3 = 1154.70 m**, at **t = (4D²/(3a²))^{1/4} = 8.774 s**, elevation **e = 30.000°**. There **ė_max = (3/4) a t/D = 0.09871 rad/s = 5.655 °/s**. (Numerical maximum on a 10⁻⁴ s grid: 0.098706 rad/s at 8.7738 s.)

(c) At t = 8.774 s: R = 2309.40 m, Ṙ = zż/R = 131.61 m/s, R̈ = 37.500 m/s², and ë = 0 (ė is at its maximum), so that a_R = R̈ − Rė² = **15.000 m/s²** and a_e = Rë + 2Ṙė = **25.981 m/s²**. These are a sin 30° = 15.000 and a cos 30° = 25.981 m/s². The rocket's acceleration is one vertical vector of 30 m/s²; its polar components are just its projections on e_R (30° above the horizon) and e_e (perpendicular to it). Because ë = 0 there, a_e = 2Ṙė = 2 × 131.61 × 0.09871 exactly: the whole transverse acceleration is the 2Ṙė term.

(d) ë = d/dt [D a t/(D² + a²t⁴/4)] is largest at lift-off: **ë(0) = a/D = 0.0150 rad/s²**, falling to zero at t = 8.774 s and then negative. (The probe's largest |ë| on the grid was 0.0150 rad/s² at the first time step.)

**Independent check.** a_R² + a_e² = 15.0² + 25.98² = 900 = a² ✓.

**Common mistakes.** Maximising e instead of ė; treating ë = 0 as meaning no transverse acceleration.

## 4. Your own noisy radar (s = 15)

(a) RMS errors (m/s²) at h = 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14 s: 30.47, 7.58, 3.35, 1.90, 1.30, **1.087**, 1.12, 1.30, 1.86, 2.57, 3.38. Best step **6 s**, error **1.087 m/s²**. `optimal_step` with the per-axis noise 7.93 m and f⁗ = Rω⁴ = 0.2 m/s⁴ gives h* = **5.84 s** and the bound 1.14 m/s²; the minimum is flat, so 6 s and 5.84 s are the same answer.

(b) σ_θ = 2 mrad: per-axis noise √((5² + 20.09²)/2) = 14.64 m, a factor 1.846. Predicted: h* × 1.846^{1/4} = × 1.166, so 6.81 s; error × 1.846^{1/2} = × 1.36. Measured: best step **7 s**, error **1.427 m/s²**, a factor **1.31** on the error. The prediction is a little high because the vertical noise and the truncation error did not double.

(c) At 20 Hz the best error is **1.084 m/s²** at 6 s: no better. The stencil uses only three samples, at t − h, t and t + h; doubling the rate at a fixed step changes nothing it uses. To gain from more samples, fit a curve through all the samples in the window.

**Common mistakes.** Comparing a single seeded record's optimum with the model to more figures than the flat minimum allows; differencing at stride 1 because the data are there.

## 5. Is that 6 g real?

(a) At h = 0.5 s the lesson's noise study gives an RMS error of **125 m/s², about 12.8 g**. A 6 g reading is well inside the noise: it is **not evidence** of anything.

(b) Use h ≈ 6 s (error about 1.1 m/s², 0.11 g). The price is time resolution: a central difference over ±6 s blends data from 12 s; a manoeuvre that starts suddenly shows up spread over that window, and a real-time display (which cannot use future samples) would lag by about 6 s.

(c) (1) A slowly drifting azimuth bias adds δθ̇ to θ̇ and so a false centripetal term ≈ 2rθ̇ δθ̇ and a false 2ṙ δθ̇ term; at the lesson's t = 12 s state a drift of 0.1 mrad/s gives 2rθ̇δθ̇ = 0.036 m/s² and 2ṙδθ̇ = −0.018 m/s²: small, but it does not average out with a longer step. (2) An azimuth jump of 2π that is not unwrapped, or a pass close to the radar's axis where θ̇ ≈ V/d and θ̈ ≈ 0.65(V/d)² become enormous, puts huge spurious values into (3).

(d) An onboard inertial record of the aircraft's own accelerations, time-aligned with the radar (or a second radar at another site, which sees different geometry and independent noise).

## Check-yourself answers

1. ė_r = θ̇ e_θ, ė_θ = −θ̇ e_r. A unit vector has constant length, so d(e·e)/dt = 2e·ė = 0: its derivative is perpendicular to it.
2. r̈ (rate of change of range rate), −rθ̇² (centripetal), rθ̈ (angular acceleration), 2ṙθ̇ (the coupling of radial and transverse motion); all in m/s².
3. a_r and a_θ are zero; r̈ and θ̈ are not (r̈ = rθ̇² and θ̈ = −2ṙθ̇/r).
4. a_t = 0, a_n = 150²/1500 = 15.0 m/s² toward the centre.
5. Its tangent tilts upward, so the path turns less sharply per metre travelled: κ = R/(R² + b²) < 1/R.
6. The noise in the second difference is √6 σ/h²: it grows as h shrinks, while the truncation error falls only as h².
7. By 16^{1/4} = 2: the best step doubles (and the best error grows by 16^{1/2} = 4).
8. The axis r = 0 is a singularity: the azimuth and its rates are undefined there, and close to it θ̇ ≈ V/d and θ̈ ≈ 0.65 (V/d)² become huge, so tiny position errors produce enormous apparent accelerations.
