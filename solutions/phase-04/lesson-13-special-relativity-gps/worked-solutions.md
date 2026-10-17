# P04.13  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for $s = 15$; numerical results are from `solution.py` (and the probe script used to write this file) and given to the figures shown.

## 1. An airliner's clock ($s = 15$: v = 230 m/s)

**Assumptions.** Non-rotating spherical Earth, ground clock at rest on it, constant speed and altitude for 10 h = 36 000 s; weak field.

(a) `time_dilation_rate(230)` = −v²/2c² = −2.943 × 10⁻¹³; over 36 000 s: **−10.59 ns**.

(b) (GM/c²)(1/R − 1/(R + 11 km)) = 1.197 × 10⁻¹²; over 36 000 s: **+43.10 ns**. The local form gh/c² gives 43.21 ns; the 0.26 % difference is the h/R ≈ 0.17 % fall of g over 11 km and the use of the standard g instead of GM/R² (9.798 m/s² at the equator).

(c) Gravity wins: net **+32.5 ns**. A clock stable to 10⁻¹³ accumulates 3.6 ns of its own error over 10 h, so the net is about nine times the clock noise and is measurable; the velocity term alone is three times the noise. This is the scale of the Hafele–Keating experiment discussed in OpenStax 5.3.

**Check.** Both terms are second order in small quantities (v/c and gh/c²); a rough order: (230/3 × 10⁸)²/2 ≈ 3 × 10⁻¹³, 10⁵/10¹⁷ ≈ 10⁻¹². Consistent.

**Common mistakes.** Using R + h in the velocity term (it has no R); using 1/R − 1/h; forgetting that eastward flight on a rotating Earth adds to the speed in the inertial frame (outside this problem's model).

## 2. The Lorentz transformation and what follows

(a) Put x = ct, x′ = ct′ into x′ = A(x − vt): ct′ = A(c − v)t. Into x = A(x′ + vt′): ct = A(c + v)t′. Multiply: c²tt′ = A²(c² − v²)tt′, so **A² = 1/(1 − v²/c²), A = γ** (positive so that the boost is continuous with the identity at v = 0). Now x = γ(x′ + vt′) with x′ = γ(x − vt): x = γ²(x − vt) + γvt′, so γvt′ = x(1 − γ²) + γ²vt. Since 1 − γ² = −γ²v²/c², t′ = γ(t − vx/c²). ∎

(b) c²Δt′² − Δx′² = γ²[(cΔt − βΔx)² − (Δx − βcΔt)²] = γ²[(c²Δt² − Δx²)(1 − β²)] + γ²[−2βcΔtΔx + 2βcΔtΔx] = c²Δt² − Δx², using γ²(1 − β²) = 1. ∎

(c) In (ct, x) the boost is L(β) = γ[[1, −β], [−β, 1]]. L(β₂)L(β₁) = γ₁γ₂[[1 + β₁β₂, −(β₁ + β₂)], [−(β₁ + β₂), 1 + β₁β₂]] = γ₁γ₂(1 + β₁β₂)[[1, −β_w], [−β_w, 1]] with **β_w = (β₁ + β₂)/(1 + β₁β₂)**. The prefactor must be γ_w: 1 − β_w² = (1 − β₁²)(1 − β₂²)/(1 + β₁β₂)², so γ_w = γ₁γ₂(1 + β₁β₂). ✓. For |β₁|, |β₂| < 1: 1 − β_w = (1 − β₁)(1 − β₂)/(1 + β₁β₂) > 0 and 1 + β_w = (1 + β₁)(1 + β₂)/(1 + β₁β₂) > 0 (the denominator is positive), so **|w| < c**. The test `test_two_collinear_boosts_compose_by_velocity_addition` checks this numerically to 10⁻¹².

(d) θ = 0: f = f₀/(γ(1 − β)) = f₀√((1 + β)/(1 − β)) (higher). θ = π: f₀√((1 − β)/(1 + β)) (lower). θ = π/2: **f = f₀/γ**, a pure time-dilation shift, second order in β. Product: f₀²·√((1 + β)/(1 − β))·√((1 − β)/(1 + β)) = **f₀²**. For small β the head-on case is f₀(1 + β + β²/2 + …), whose first-order part is the classical Doppler shift.

**Common mistakes.** Taking A = −γ; mixing the sign of v between the forward and inverse transformations; in (d), using the angle in the source frame (the formula changes).

## 3. A navigation constellation at Mars ($s = 15$: h = 450 km)

**Assumptions.** Spherical, non-rotating Mars; weak field; circular orbits; ground clock on the surface.

(a) a = 3839.5 km, v = √(GM/a) = 3339.8 m/s. **Velocity −5.36 µs/day, gravity +1.42 µs/day, net −3.94 µs/day.** At this altitude the velocity term dominates: the satellite clock runs slow.

(b) `crossover_radius(4.2828e13, 3389.5e3)` = **5084.25 km** = 1.5 R, altitude **1694.75 km**. The balance GM/2a = GM(1/R − 1/a) divides through by GM, leaving 1/R = 3/(2a): pure geometry, the same for every spherical body.

(c) Kepler: a = (GM T²/4π²)^{1/3} = **20 448 km** for T = 88 775 s. There v = 1447 m/s; terms −1.01 and +10.13 µs/day, **net +9.13 µs/day**. The 450 km satellite (net negative) must be *sped up* before launch, by 4.56 × 10⁻¹¹; the areostationary one (net positive) *slowed*, by 1.056 × 10⁻¹⁰.

**Check.** Every Mars term is smaller than its Earth counterpart at the same R-scaled radius by GM_Mars R_E/(GM_E R_Mars) = 0.20: GEO-like Earth net at 6.6 R is +46.6, Mars areostationary at 6.03 R is +9.1 ≈ 0.2 × 46. ✓

**Common mistakes.** Leaving ω and J₂ at their Earth values in `gps_clock_budget`; using the altitude instead of the radius.

## 4. The periodic terms ($s = 15$: e = 0.025)

(a) Amplitude 2√(GMa) e/c² = **57.2 ns**, **17.2 m** of range. n = √(GM/a³) = 1.4584 × 10⁻⁴ rad/s; maximum rate **8.35 × 10⁻¹² s/s**, which is **1.87 %** of the factory offset 4.4647 × 10⁻¹⁰ (1.5 % at e = 0.02).

(b) Sagnac term (ns) for an equatorial receiver; position angle p from the western (−90°) to the eastern (+90°) horizon:

| p (°) | −90 | −75 | −60 | −45 | −30 | −15 | 0 | 15 | 30 | 45 | 60 | 75 | 90 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Δt (ns) | 133.43 | 128.89 | 115.56 | 94.35 | 66.72 | 34.54 | 0 | −34.54 | −66.72 | −94.35 | −115.56 | −128.89 | −133.43 |

The curve is 133.43 sin(−p) ns: odd, zero overhead, ±40 m at the horizons.

(c) Neither grows: the eccentricity term is periodic with the orbit (worst 17.2 m here), the Sagnac term depends only on geometry (worst 40 m). Both are bounded, so they cannot be "averaged out" by waiting, but neither runs away; that is why they are computed in the receiver rather than removed by a frequency offset.

**Check.** The amplitude is linear in e: 45.8 ns × 0.025/0.02 = 57.2 ns ✓.

**Common mistakes.** Using the mean anomaly instead of the eccentric anomaly (the difference is O(e), so the amplitude is right but the phase is not); reporting the Sagnac term with the wrong sign convention (positive when the receiver moves away from the satellite).

## 5. Why the factory offset uses −3GM/2a

(a) T = 2π/n. ⟨1/r⟩ = (1/T)∫₀ᵀ dt/r = (n/2π)∫₀^{2π} [(1 − e cos E)/n] dE / [a(1 − e cos E)] = (1/2π)∫₀^{2π} dE/a = **1/a**, exactly, for every e. ∎

(b) Vis-viva: v²/2 = GM/r − GM/2a, so −(v²/2 + GM/r)/c² = −(2GM/r − GM/2a)/c². Average with (a): −(2GM/a − GM/2a)/c² = **−3GM/(2ac²)**. The remainder is −2GM(1/r − 1/a)/c². Its time integral: ∫(1/r − 1/a)dt = ∫[1/(a(1 − e cos E)) − 1/a](1 − e cos E)/n dE = (1/(an))∫ e cos E dE = e sin E/(an). With an = √(GM/a), −2GM e sin E/(c² √(GM/a)) = **−2√(GMa) e sin E/c²**, equation (4). ∎

(c) e = 0.1, `quad` over E of rate(E)(1 − e cos E)/2π: mean **−2.5046 × 10⁻¹⁰**, equal to −3GM/(2ac²) to the printed digits (relative difference 0 in double precision). Instantaneous rates relative to a clock far away: perigee **−24.85 µs/day**, apogee **−19.02 µs/day**, mean **−21.64 µs/day**. The ground clock's −60.2 µs/day (Φ₀/c²) is the same at every point, so relative to the ground these are +35.4, +41.2 and +38.6 µs/day.

**Common mistakes.** Averaging over E or over true anomaly instead of over time (the weights differ); forgetting that the averaged factor −3/2 includes both terms.

## 6. Which corrections does a service need?

Range-error rates and bounds (from the lesson's budget): velocity 0.0250 m/s, gravity 0.158 m/s, geoid terms 4.7 × 10⁻⁴ m/s (secular); eccentricity 13.7 m at e = 0.02, Sagnac 40 m (bounded).

| term | 1 m service | 10 m service | 1 km service |
|---|---|---|---|
| velocity | yes (1 m in 40 s) | yes (400 s) | yes (11.1 h < 1 day) |
| gravity | yes (6.3 s) | yes (63 s) | yes (1.75 h) |
| geoid (rotation + J₂) | yes (35 min) | yes (5.9 h) | no (24 days; the daily upload absorbs it) |
| eccentricity, e = 0.02 | yes (13.7 m) | yes (13.7 m > 10 m) | no |
| Sagnac | yes (40 m) | yes (40 m) | no |

The term never to omit is the gravitational one: it is the largest (0.16 m of range every second), it reaches even the 1 km requirement in under two hours, and it is the one that sets the sign of the net offset.

**Common mistakes.** Treating the periodic terms as secular (they do not grow); forgetting that a once-daily clock upload hides any constant rate below 1 km/day of range, which is exactly what the J₂ and rotation terms are for the coarse service.

## Check-yourself answers

1. The laws of physics are the same in every inertial frame; c is the same in every inertial frame. γ = 1/√(1 − v²/c²); γ(0.6c) = 1.25, so a clock moving at 0.6c ticks once for every 1.25 ticks of yours.
2. −v²/2c²; for GPS (v = 3873.8 m/s) −8.35 × 10⁻¹¹, i.e. **−7.21 µs/day**.
3. Rate difference (Φ₂ − Φ₁)/c²; the clock at the higher (less negative) potential runs faster.
4. Velocity −7.21 µs/day, gravity +45.65 µs/day, net +38.4 (two terms) or +38.6 µs/day with the geoid terms. They do not cancel because GPS is far above the crossover a = 3R/2 (altitude 3189 km), where they would.
5. It varies around the orbit (zero mean), so no constant frequency offset removes it; the amplitude at e = 0.02 is 45.8 ns (13.7 m).
6. The receiver rotates with the Earth while the signal is in flight, adding a path ω ẑ·(r_s × r_r)/c; zero when the satellite is overhead (or in the receiver's meridian plane); at most 133.4 ns (40 m) for GPS.
7. f = f₀/(γ(1 − β cos θ)); on a symmetric pass the first-order part is antisymmetric and averages out, the second-order part f₀(1/γ − 1) (−0.132 Hz at L1) survives.
8. K = (γ − 1)mc² = 2.9645 × 10⁷ J, which is ½mv² × (1 + 4.9 × 10⁻¹⁰): relativity changes it in the tenth digit.
