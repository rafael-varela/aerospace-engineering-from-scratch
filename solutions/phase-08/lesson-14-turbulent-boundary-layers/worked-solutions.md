# P08.14  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and probe scripts that call it, given to the figures shown. Skin-friction formulas are **correlations**; every velocity profile used here is **synthetic**.

## 1. The fuselage (s = 15: L_f = 6.5 m)

**Assumptions.** Smooth flat plate, zero pressure gradient, turbulent from the nose, incompressible (M = 0.18).

(a) Re_L = 60 × 6.5/1.8630 × 10⁻⁵ = **2.09 × 10⁷**. C_F: power law **0.002541**, Prandtl–Schlichting **0.002676**, Kármán–Schoenherr **0.002609**, ITTC **0.002649**: a **5.3 %** spread. The power law is the lowest because Re_L is above the range (Re < 10⁷) its fit was made for.

(b) White's local c_f = **0.002307**; u_τ = 60√(0.002307/2) = **2.038 m/s**; τ_w = ρu_τ² = **3.78 Pa**; ν/u_τ = **9.14 μm**; δ = **82.8 mm** (1/7th law); δ⁺ = **9050**.

(c) First cell at y⁺ = 1: **9.1 μm**; the log layer starts at y⁺ = 30, **0.27 mm** above the skin. The whole inner layer is thinner than a coat of paint's roughness tolerance would suggest, which is why surface finish matters for turbulent drag.

(d) Laminar C_F = 1.328/√Re_L = 2.90 × 10⁻⁴: the turbulent value is **9.2 times** larger (Prandtl–Schlichting).

**Check.** u_τ/U = 3.4 %, typical of a turbulent layer; δ/L = 1.3 %.

**Common mistakes.** Using the trailing-edge c_f as the average; using ν/U instead of ν/u_τ for y⁺.

## 2. Derivations

(a) Inner law: U = u_τ f(y⁺), so y dU/dy = u_τ y⁺f′(y⁺). Defect law: U_e − U = u_τ g(y/δ), so y dU/dy = −u_τ (y/δ)g′(y/δ). In the overlap both hold for any y⁺ and y/δ (their ratio δ⁺ is large and arbitrary), so each side must equal a constant, written u_τ/κ. Then dU/dy = u_τ/(κy) and **u⁺ = ln(y⁺)/κ + B**. Spalding: with k = κu⁺, e^k − 1 − k − k²/2 − k³/6 = k⁴/24 + …, so for small u⁺, y⁺ = u⁺ + O(u⁺⁴): **u⁺ = y⁺**. For large u⁺ the exponential dominates: y⁺ ≈ e^(−κB)e^(κu⁺), i.e. **u⁺ ≈ ln(y⁺)/κ + B** (the code's u⁺ differs from the log law by less than 0.05 for y⁺ ≥ 10³).

(b) With η = y/δ: δ*/δ = ∫₀¹(1 − η^(1/7))dη = 1 − 7/8 = **1/8**; θ/δ = ∫₀¹(η^(1/7) − η^(2/7))dη = 7/8 − 7/9 = **7/72**; **H = (1/8)/(7/72) = 9/7**. Momentum integral: (7/72)dδ/dx = 0.0225(ν/(Uδ))^(1/4) ⇒ δ^(1/4)dδ = (72/7)(0.0225)(ν/U)^(1/4)dx ⇒ (4/5)δ^(5/4) = 0.2314(ν/U)^(1/4)x ⇒ **δ/x = (0.2893)^(0.8)Re_x^(−1/5) = 0.3707Re_x^(−1/5)**. Then c_f = 2τ_w/(ρU²) = 0.045(ν/(Uδ))^(1/4) = 0.045 × 0.3707^(−1/4)Re_x^(−1/5) = **0.0577Re_x^(−1/5)**. Averaging x^(−1/5) over 0..L gives (1/L)∫x^(−1/5)dx = (5/4)L^(−1/5): **C_F = (5/4)c_f(L)** (0.074/0.0592 = 1.25 for the fitted constants).

(c) Per unit span, the drag of the plate is ½ρU²[turbulent drag of 0..L, minus turbulent drag of 0..x_tr, plus laminar drag of 0..x_tr]. Dividing by ½ρU²L: **C_F = C_F,t(Re_L) − (x_tr/L)[C_F,t(Re_tr) − C_F,l(Re_tr)]**, with x_tr/L = Re_tr/Re_L. Assumption: the turbulent layer behind x_tr is the same as if it had been turbulent from the leading edge, i.e. it has the momentum thickness of a fully turbulent layer at x_tr, not the much smaller laminar θ.

**Common mistakes.** Integrating the 1/7th law to the wall for c_f (its slope is infinite there; that is why Blasius's pipe law is needed).

## 3. A rake on the wing (synthetic profile, s = 15)

At x = 1.2 m, Re_x = 3.86 × 10⁶, White's c_f = 0.002981, so the profile was built with **u_τ = 2.3165 m/s** and δ = 21.42 mm (Π = 0.55, 2 % noise, seed 15). The lowest point is at y⁺ = 6.2.

(a) Default window: **11 points**, u_τ = **2.3309 m/s (+0.62 %)**; τ_w = **4.94 Pa** against the 4.88 Pa that built the profile (+1.2 %, twice the u_τ error).

(b)

| window / constants | points | u_τ (m/s) | error |
|---|---|---|---|
| y⁺ ≥ 30, y/δ ≤ 0.1 | 8 | 2.3068 | −0.42 % |
| y⁺ ≥ 30, y/δ ≤ 0.2 (default) | 11 | 2.3309 | +0.62 % |
| y⁺ ≥ 30, y/δ ≤ 0.3 | 13 | 2.3380 | +0.93 % |
| y⁺ ≥ 50, y/δ ≤ 0.2 | 9 | 2.3418 | +1.09 % |
| κ = 0.384, B = 4.17 | 11 | 2.3326 | +0.69 % |

(c) The window moves the answer most (−0.4 to +1.1 %): a wider window takes in the wake (above the log law) or drops points, and the noise in a handful of points matters. The constants matter little because the two log laws cross at **y⁺ = 152**, inside the window (y⁺ 30–620 here), so they disagree by little over the points fitted. **No** for y⁺ > 500: at y⁺ = 1000 the two laws give u⁺ = 21.85 and 22.16 (1.4 % apart), so the fitted u_τ would move by about that much, and the wake would contaminate the window as well.

**Common mistakes.** Fixing the window in millimetres before knowing u_τ; reading the fit's agreement as a measurement of κ (the profile was built with κ = 0.41).

## 4. Correlations and the mixed plate (s = 15: Re_L = 4.5 × 10⁶)

(a) Integrated: power **0.003456**, Prandtl–Schlichting **0.003424**, Kármán–Schoenherr **0.003354**, ITTC **0.003464**. Integrals of the local laws: power **0.003452**, Schultz-Grunow **0.003387**, White **0.003543**. The full spread is 5.6 % (0.003354 to 0.003543). Within one family the agreement is closer: the power law's local and integrated forms agree to 0.1 %, Schultz-Grunow's local law against Prandtl–Schlichting to 1.1 %. Choosing a correlation is a ±3 % decision.

(b)

| Re_tr | Prandtl (P–S) | Prandtl (power) | momentum match | x₀/L |
|---|---|---|---|---|
| 3 × 10⁵ | 0.00321 | 0.00322 | 0.00325 | 0.044 |
| 5 × 10⁵ | 0.00307 | 0.00307 | 0.00315 | 0.080 |
| 10⁶ | 0.00273 | 0.00271 | 0.00289 | 0.175 |
| 2 × 10⁶ | 0.00209 | 0.00207 | 0.00232 | 0.371 |

Fully turbulent 0.00342, fully laminar 0.000626.

(c) Prandtl's subtraction gives the turbulent part behind x_tr the thick layer (large θ, low c_f) it would have had if turbulent from the leading edge; the momentum match starts it with the thin laminar θ at x_tr, so its c_f is higher just behind transition. The later the transition, the larger the difference between the two θ's and the larger the gap: 1 % at Re_tr = 3 × 10⁵, 11 % at 2 × 10⁶. For a laminar-flow wing, where transition is far back, trust the **momentum match**: it conserves the momentum deficit across transition, which is the physics the momentum integral says matters; Prandtl's formula is a convenient fit for early transition.

**Common mistakes.** Mixing a correlation family between the turbulent terms (use the same one in both).

## 5. A friction law with no correlation in it

(a) c_f from the log-wake law (×10⁻³):

| Re_x | Π = 0.4 | Π = 0.55 | Π = 0.7 | κ = 0.384, B = 4.17 | White |
|---|---|---|---|---|---|
| 10⁶ | 3.624 | 3.486 | 3.350 | 3.394 | 3.759 |
| 4.83 × 10⁶ | 2.832 | 2.739 | 2.647 | 2.645 | 2.876 |
| 10⁷ | 2.546 | 2.467 | 2.388 | 2.374 | 2.570 |

Π = 0.4 comes within 1–4 % of White; Π = 0.55 is 4–7 % below; the alternative κ, B lower c_f by 2.6–3.8 %.

(b)

| δ⁺ | Ue⁺ | c_f | Re_θ | H |
|---|---|---|---|---|
| 10³ | 24.53 | 3.32 × 10⁻³ | 2.76 × 10³ | 1.408 |
| 10⁴ | 30.15 | 2.20 × 10⁻³ | 2.95 × 10⁴ | 1.285 |
| 10⁵ | 35.76 | 1.56 × 10⁻³ | 3.08 × 10⁵ | 1.227 |

H falls towards 1 as Re grows: the turbulent profile gets fuller (the viscous sublayer is a shrinking fraction of the layer), and every decade of δ⁺ adds ln(10)/κ = 5.6 to Ue⁺.

(c) κ, B and Π (and how Π varies with Re_θ at low Re). Compare with published DNS profiles: Schlatter & Örlü's zero-pressure-gradient boundary layers (Re_θ 670–4060, `@kth-schlatter-orlu-tbl-dns-data`) for Π and c_f, and Lee & Moser's channel profiles (`@lee-moser-2015-channel-dns-data`) for κ and B in the log layer. Both are free to download under their citation terms, with no license to redistribute, so the learner downloads them and nothing is committed.

## 6. What do you tell the designer?

(a) About **90 N, range 83–153 N**. The central value uses transition at the e⁹ point of `P08.12` (83 N by Prandtl's subtraction, 95 N by the momentum match, which is more defensible). The upper end, 150–153 N, is transition at Re_x = 5 × 10⁵ (0.155 m), what insects, rain, propeller wash or a rough leading edge can produce; the lower end is the e⁹ estimate itself, because a production wing is unlikely to stay laminar longer than a smooth flat plate in quiet air. The all-laminar 29.7 N is not a credible bound.

(b) (1) The pressure gradient of the real aerofoil: a favourable gradient ahead of the peak delays transition (lower drag), the adverse recovery thickens the turbulent layer and adds pressure (form) drag (higher). (2) Surface roughness, steps and rivets: raise c_f of the turbulent layer and trip transition early (higher). Also: wing–fuselage interference and the propeller slipstream (higher).

(c) The transition line in flight (infrared thermography or hot films): it fixes the term that moves the answer by a factor of two; the choice of correlation moves it by 3 %.

## Check-yourself answers

1. u_τ = √(τ_w/ρ); y⁺ = yu_τ/ν; u⁺ = U/u_τ. They are the velocity and length scales of the inner layer, set by the wall shear and viscosity alone.
2. Viscous sublayer u⁺ = y⁺ (y⁺ ≲ 5), buffer layer (5–30), log layer u⁺ = ln y⁺/κ + B (y⁺ ≳ 30 up to about 0.2δ), then the wake region.
3. Choose the u_τ that makes the measured points in the log region lie on the log law; it gives the wall shear without measuring it, but it assumes the log law and its constants.
4. Turbulent mixing carries momentum towards the wall, so the profile is fuller and the wall gradient steeper: C_F is 5.6 times the laminar value at the wing's Re_L and 9.2 times at the fuselage's, and it falls only as about Re^(−1/5).
5. 0.0592Re_x^(−1/5) local and 0.074Re_L^(−1/5) average; the average is 5/4 of the end value.
6. On the wing: 166 N turbulent from the nose (Prandtl–Schlichting; 163–168 N over four correlations) against 29.7 N laminar, 5.6 times.
7. It treats the turbulent layer behind transition as if it had grown from the leading edge (too thick), so it under-estimates the drag when transition is late (11 % at Re_tr = 2 × 10⁶ here).
8. Mass and x-momentum, in integral form across the layer (the momentum integral) and in the averaged differential form with a turbulent (Reynolds) stress that the log law and the correlations model empirically.
