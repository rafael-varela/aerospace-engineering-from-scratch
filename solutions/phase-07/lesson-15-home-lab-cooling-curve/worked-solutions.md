# P07.15  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe scripts that call it, given to the figures shown. Results from the generator are **SYNTHETIC**, not measurements.

## 1. Predict your own bar (s = 15: L = 55 mm)

(a) V = 27.869 cm³, A = 54.02 cm², L_c = V/A = **5.159 mm**, ρcV = **67.42 J/K**.

(b) At θ = 14 K in 21.0 °C air: h_conv = **5.639**, h_rad = **0.558 W/(m² K)**, **τ = 2014 s (33.6 min)**, **Bi = 1.91 × 10⁻⁴**.

(c) Longer, slightly. τ = (ρc/h)(V/A) and V/A = DL/(4L + 2D) grows with L toward D/4 = 6.35 mm, because the ends are a smaller share of the area: 5.064 mm at 50 mm, 5.159 mm at 55 mm. (h_conv also falls a little, since the ends, which have the larger h, are a smaller share.)

**Check.** τ ∝ L_c/h: 1968 × (5.159/5.064) × (6.22/6.197) ≈ 2012 s, consistent with 2014 s.

## 2. Derivations

(a) ρcV dθ/dt = −hAθ ⇒ dθ/θ = −dt/τ ⇒ **θ = θ₀e^(−t/τ), τ = ρcV/(hA)**. The heat hAθ leaving the surface must first cross the metal; its conduction resistance over a length of order L_c is L_c/(kA), in series with the surface resistance 1/(hA). The temperature drop splits in proportion to the resistances: internal/external = (L_c/kA)/(1/hA) = hL_c/k = **Bi**. For Bi ≪ 1 the internal drop is negligible.

(b) Try T_m − T∞ = a e^(−t/τ) + b e^(−t/τ_p). The b term solves the homogeneous equation; substituting the a term: −(τ_p/τ)a + a = θ₀ ⇒ **a = θ₀τ/(τ − τ_p)**; T_m(0) = T∞ ⇒ b = −a. So T_m = T∞ + θ₀τ/(τ − τ_p)[e^(−t/τ) − e^(−t/τ_p)]. For t ≫ τ_p the second exponential is gone and the reading is an exponential with **the same τ** and amplitude θ₀τ/(τ − τ_p) (1.01 θ₀ for τ_p = 20 s, τ = 1800 s).

(c) With T = T_m + d, T∞ = T_m − d: (T² + T∞²)(T + T∞) = (2T_m² + 2d²)(2T_m) = 4T_m³(1 + d²/T_m²), so 4εσT_m³ is low by the factor 1/(1 + (d/T_m)²), with d/T_m = (T − T∞)/(T + T∞).

| θ (K) | exact h_rad (W/(m² K)) | 4εσT_m³ |
|---|---|---|
| 5 | 0.53294 | 0.53290 |
| 14 | 0.55782 | 0.55752 |
| 28 | 0.59854 | 0.59731 |

The approximation is good to 0.2 % over the lab's range.

(d) For small λ, J₀ ≈ 1 − λ²/4 and J₁ ≈ λ/2, so λ(λ/2) ≈ Bi(1 − λ²/4) ⇒ **λ₁² ≈ 2Bi** (to first order in Bi). Then C₁ = 2J₁/(λ(J₀² + J₁²)) → 1, J₀(λ₁r/R) → 1 everywhere, and the higher eigenvalues are at least 3.83 (the zeros of J₁), so their terms die as e^(−14.7Fo). What remains is **exp(−2Bi Fo) = exp(−2ht/(ρcR))**, which is exp(−t/τ) with V/A = R/2 for an infinite cylinder. The tests confirm λ₁² = 2 × 10⁻⁴ at Bi = 10⁻⁴.

## 3. The bar in a thermal-vacuum chamber

(a) Bare (ε = 0.09): h_rad = **0.558 W/(m² K)**, τ = **21 961 s = 6.10 h**, 11 times the 32.8 min in room air. Without air, a bare aluminium part barely cools.

(b) Exact 50 → 30 °C: **25 309 s**; linearised τ ln(29/9) = **25 696 s** (1.5 % long, because h_rad at 14 K under-represents the hotter early part).

(c) Black tape (ε = 0.9): h_rad = 5.578 W/(m² K), τ = **2196 s (0.61 h)**; exact 50 → 30 °C **2531 s**, linearised **2570 s**.

(d) With T_w = 0: ρcV dT/dt = −εσAT⁴ ⇒ ∫dT/T⁴ = −(εσA/ρcV)t ⇒ **t = ρcV/(3εσA)(1/T₂³ − 1/T₁³)**. From 323.15 to 303.15 K: **5009 s** bare, **501 s** black. A real chamber's shroud at liquid-nitrogen temperature is close to this limit.

**Check.** All times scale as 1/ε, as they must for pure radiation.

## 4. Your own SYNTHETIC run (s = 15: seed 730)

SYNTHETIC; room readings 21.3, 21.3, 21.2 °C (mean 21.27 °C); first readings 21.4, 36.1, 43.0, 46.1 °C (the probe lag).

(a)

| fit | τ (s) | u(τ) (s, statistical) | χ²/dof | other |
|---|---|---|---|---|
| log-linear, T∞ = room | 1774.7 | 1.0 | 36.4 | — |
| nonlinear, T∞ fixed | 1770.9 | 1.1 | 35.3 | — |
| nonlinear, T∞ free | 1484.0 | 2.9 | 2.24 | T∞ = 23.11 ± 0.02 °C |
| lagged, whole record | 1763.3 | 1.0 | 37.6 | τ_p = 18.2 ± 0.06 s |

The statistical uncertainties are far too small, because χ²/dof ≫ 1: the model, not the noise, limits them.

(b) Prediction at 21.27 °C: 1968 ± 197 s; **E_n = 0.49**: consistent.

(c) Windows: **1615 ± 3, 1862 ± 6, 2121 ± 13 s** (predicted 1851, 2121, 2399 s at the windows' mean excesses 18.7, 9.8, 5.4 K). h is not constant; it falls as the bar cools, and the rise of τ (31 %) matches the correlation's (30 %).

(d) **f_conv = 1.149 ± 0.001**, χ² = 199 for 231 dof, with ε fixed at 0.09. The truth is f_conv = 1.12 with ε = 0.12: both raise h, and over this small temperature range radiation is nearly linear in θ while convection goes as θ^(1/4), so the data distinguish them only weakly; a smaller ε is compensated by a larger f. To separate them, change ε deliberately (the black-taped bar) while keeping the room the same.

**Common mistakes.** Reporting 1774.7 ± 1.0 s as if the fit were good; trusting the free-T∞ fit because its χ² is lowest.

## 5. When does the lumped model fail?

| case | Bi | τ (s) | centre at t = τ | corner at t = τ | lumped |
|---|---|---|---|---|---|
| aluminium, air (h 6.9) | 2.09 × 10⁻⁴ | 1775.4 | 0.3681 | 0.3678 | 0.3679 |
| aluminium, water bath (h 500) | 1.52 × 10⁻² | 24.5 | 0.3806 | 0.3599 | 0.3679 |
| 304, air | 2.53 × 10⁻³ | 2348.4 | 0.3700 | 0.3665 | 0.3679 |
| 304, water bath | 0.184 | 32.4 | 0.5110 | 0.2853 | 0.3679 |

(b) Aluminium in air: uniform to 3 × 10⁻⁴. Aluminium in water: about 2 % spread, still lumped to a few per cent. Steel in air: lumped. Steel in water: Bi = 0.18, the centre is 40 % above the lumped value and the corner 22 % below: not lumped.

(c) The warm-up is the water-bath row. An aluminium bar reaches the water temperature within a few τ ≈ 25 s; a steel one needs longer and has a cooler core for the first minutes. Leave either in the water for 5 minutes (more than ten water-bath time constants) so that it starts the cooling run uniform.

## 6. Why is the measured τ 20 % short?

(a) 20 % short means about 0.8 × 1968 = 1574 s. Convection 30 % stronger gives **1546 s**; ε = 0.3 gives **1628 s**. Both are plausible: a 30 % draught effect is common in a room with any air movement; ε = 0.3 needs a heavily oxidised or anodised bar (Table 10.1's upper range). Tell them apart with windows (draughts make the windowed τ scatter and break the smooth θ^(1/4) trend) and with the black-tape comparison or a cardboard draught shield.

(b) T∞ 1 K too high: **τ = 1615 s** instead of 1775 s for seed 730, **shorter by 9 %**, and the end of the semi-log plot bends down.

(c) A 2 g tip adds **1.6 %** and a 5 g tip **4.1 %** to ρcV = 61.29 J/K: τ *longer*. The stem acting as a fin draws extra heat from the bar: τ *shorter*. The two partly cancel; neither explains 20 %.

(d) Ranking: (1) draughts (largest and commonest; test: repeat inside a cardboard enclosure open at the bottom); (2) a wrong T∞ (test: three room readings, fit with T∞ fixed); (3) a higher emittance (test: black-tape bar alongside, or a bar of known finish); (4) probe effects (test: repeat with the thermistor, which is much lighter, under foil tape). A χ²/dof of 60 with a smooth, curved residual points first to the variable h and a T∞ error; a jagged residual points to draughts.

## Check-yourself answers

1. ρcV dT/dt = −hA(T − T∞); θ = θ₀e^(−t/τ) with τ = ρcV/(hA).
2. Bi = hL_c/k, the ratio of the internal conduction resistance to the surface resistance. For the lesson's bar it is about 2 × 10⁻⁴: the internal temperature differences are two ten-thousandths of the surface-to-air difference.
3. τ depends on the dimensions to the first power and they are measured to 0.2 %; the correlation is uncertain by about 10 % and h enters τ directly.
4. To avoid scalds: the hard limit is **60 °C**, and tap-hot water can be hotter than that, so it is measured, and cold water added if needed, before the bar goes in.
5. The first 5τ_p (about 2 min). After that the lag only scales the amplitude by τ/(τ − τ_p); the decay rate, and so τ, is unchanged.
6. One exponential does not describe the data to the reading noise: h is changing (natural convection weakens as the bar cools), or T∞ is wrong. Fit windows, check T∞, and fit the variable-h model.
7. T∞ and τ are strongly correlated; when the model is too simple, T∞ absorbs the curvature, moving away from the true room temperature (1.8 K in the lesson's SYNTHETIC case) and biasing τ by 16 %.
8. So that the comparison is a test: a prediction adjusted after seeing the data can always be made to agree, and then it tests nothing.
