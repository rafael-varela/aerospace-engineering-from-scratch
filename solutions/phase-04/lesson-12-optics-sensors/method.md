# P04.12  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
n₂ = 1.44 + s/1000; θ_c = arcsin(n₂/n₁). At the end face sin θ_a = n₁ sin θ_core with θ_core = 90° − θ_c, so sin θ_a = n₁ cos θ_c = √(n₁² − n₂²). Check: `refraction_angle(1.0, n1, theta_a)` and 90° minus it equals θ_c. R = ((1.76 − 1)/(1.76 + 1))²; two surfaces pass (1 − R)².

## 2.
(a) Edge of a detector of width np at x = np/2 = f tan(FOV/2); one pixel near the axis: dα = dx/f. Since arctan u < u for u > 0, 2 arctan(np/2f) < np/f; Taylor arctan u = u − u³/3. (b) x_i = x₀ + e_i + r_i with Var e = σ², r uniform on [−½, ½] (variance 1/12) when x₀ + e_i is spread over many pixels relative to the grid; centroid = mean of N independent x_i; variances add. (c) V_out = V_ex[R(1 + x)/(R(1 + x) + R) − ½] = V_ex x/(2(2 + x)) → V_ex x/4. Full bridge: one side (1 + x)/2, the other (1 − x)/2, difference x. (d) k = γ/(γ − 1): q_c = p[kε + k(k − 1)ε²/2 + …] = ½γpM²[1 + M²/4 + …] = ½ρV²(1 + M²/4 + …).

## 3.
R(20 °C) from `thermistor_resistance(293.15)`. dV/dT = V_s R_top/(R_top + R_T)² × dR_T/dT; for fixed R_T, R_top/(R_top + R_T)² is maximised at R_top = R_T. Numerically: central difference of `divider_voltage(3.3, R_top, thermistor_resistance(T ± 0.01))`; K per count = (3.3/4096)/|dV/dT|. Power V²/R_T at 20 °C; dissipation constant = P/(K per count).

## 4.
`StarTracker` as in the hint; `.scale/ARCSEC`, `np.degrees(.fov)`, `fov_solid_angle(fov, fov)`, `.n_stars`. `airy_radius(600e-9, 1.4)`, divided by 3.45 μm. `star_electrons(5.5, D, 0.05)`; rms of `simulate_centroid_errors` for each window; per-star = rms × scale; cross-boresight = `attitude_error_cross_boresight(rms, scale, n_stars)`. Compare each factor with the baseline: scale ratio, electron ratio, star-count ratio.

## 5.
`quarter_bridge`, linear V_ex·GF·ε/4, `full_bridge`; linearity error (linear − exact)/exact. Gains 1.6/V_bridge. Full-scale strain range of the ADC = 3.3 V divided by (gain × V_bridge/ε); divide by 4096. `bits_needed(3000e-6, 1e-6)`. For (d) think about everything between the gauge and the counts that is not the ADC.

## 6.
For each decision find the lesson result: focus → pixel locking (2.7) and the read-noise table (2.6); 9 × 9 window → the same table at 10 e⁻; one star → a single direction fixes two axes, not roll, and gives no averaging over n. Compute with `StarTracker()`, `star_electrons` and `simulate_centroid_errors` at 10 e⁻ for the alternatives. The falsification test: the lesson's adversarial self-checks.
