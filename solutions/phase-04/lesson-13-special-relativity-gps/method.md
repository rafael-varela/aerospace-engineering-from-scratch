# P04.13  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. v = 200 + 2s; ΔT_v = `time_dilation_rate(v)` × 36 000 s. 2. ΔT_g = (GM/c²)(1/R − 1/(R + h)) × 36 000 s; local form gh/c² × 36 000 s; compare. 3. Net = sum; compare 10⁻¹³ × 36 000 s with the net.

## 2.
(a) x = ct, x′ = ct′ in both: ct′ = A(c − v)t and ct = A(c + v)t′; multiply: c² = A²(c² − v²). Then substitute x′ into x = A(x′ + vt′) and solve for t′. (b) Substitute and use γ²(1 − β²) = 1. (c) Matrix L(β) = γ[[1, −β], [−β, 1]]; L(β₂)L(β₁) = γ₁γ₂(1 + β₁β₂)[[1, −β_w], [−β_w, 1]] with β_w = (β₁ + β₂)/(1 + β₁β₂); check the prefactor is γ_w. c − w ∝ (1 − β₁)(1 − β₂) > 0. (d) Substitute θ and simplify with γ = 1/√((1 − β)(1 + β)).

## 3.
(a) `gps_clock_budget(R + h, mu, R, omega=0, J2=0).us_per_day()`. (b) `crossover_radius(mu, R)`; algebra gives a = 3R/2 independent of GM. (c) a = (GM T²/4π²)^{1/3}; budget at that a; sign of net decides slow (positive net) or fast (negative net) adjustment.

## 4.
(a) A = |`eccentricity_correction(e, π/2)`|; A·c in metres; n = √(GM/a³), rate = A·n; divide by `gps_clock_budget().net`. (b) Loop p from −90° to 90° in 15° steps with `sagnac_correction`. (c) Both are bounded and periodic; the worst case is the amplitude times c.

## 5.
(a) ⟨1/r⟩ = (1/T)∫(1/r)dt = (n/2π)∫(1/(a(1 − e cos E)))(1 − e cos E)/n dE. (b) Rate = −(2GM/r − GM/2a)/c²; average with (a); write the remainder as −2GM(1/r − 1/a)/c², integrate in time using dt and r in E, and use √(GMa)/... to reach −2√(GMa)e sin E/c². (c) `quad` over E of rate(E)(1 − e cos E) divided by 2π; evaluate at E = 0 and π; convert to µs/day.

## 6.
For each term: secular rate in m/s = |fractional rate| × c; time to requirement = requirement / rate; compare with the daily upload interval. Periodic terms: compare the amplitude in metres with the requirement.
