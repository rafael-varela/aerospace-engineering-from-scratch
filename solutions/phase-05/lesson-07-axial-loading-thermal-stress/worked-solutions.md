# P05.07  Worked solutions

Each solution states assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe scripts that call it, given to the figures shown.

## 1. A drag-brace tube (s = 15: P = 35 kN)

**Assumptions.** Uniform compression, no buckling, typical 4130 (E = 205 GPa, α = 12.2 × 10⁻⁶ K⁻¹, Fty = 435 MPa).

(a) A = π(25.4² − 22.2²)/4 = **119.63 mm²**. σ = −35 000/119.63 = **−292.6 MPa**; δ = −0.8563 mm (**0.856 mm shorter**); |σ|/Fty = **0.673**.

(b) σ = −EαΔT = −205 000 × 12.2 × 10⁻⁶ × 30 = **−75.0 MPa**, whatever the tube's size.

(c) ΔT = Fty/(Eα) = 435/(205 000 × 12.2 × 10⁻⁶) = **173.9 K**.

**Check.** δ/L = −1.427 × 10⁻³ = σ/E.

**Common mistakes.** Using the outside diameter alone for the area; adding (b) to (a) as if the loaded tube were also restrained (the end load cannot be applied to a tube whose ends are fixed to rigid walls).

## 2. Four closed forms

(a) With A(x) = A₀ + (A₁ − A₀)x/L, δ = (P/E)∫₀ᴸ dx/A(x) = (P/E)·(L/(A₁ − A₀))·ln(A₁/A₀) = **PL ln(A₀/A₁)/(E(A₀ − A₁))**. As A₁ → A₀, ln(A₀/A₁) = −ln(1 + (A₁ − A₀)/A₀) ≈ (A₀ − A₁)/A₀, so δ → PL/(EA₀).

(b) Junction displacement u (positive toward segment 2). Segment 1 (fixed at its far end) elongates by u, segment 2 by −u. Hooke's law with thermal strain: elongation = NL/(EA) + αLΔT, so N₁ = k₁(u − α₁L₁ΔT), N₂ = k₂(−u − α₂L₂ΔT). The junction is pulled back by N₁ and forward by N₂ and P: N₂ − N₁ + P = 0. Substituting: −k₂u − k₂α₂L₂ΔT − k₁u + k₁α₁L₁ΔT + P = 0, so **u = (P + k₁α₁L₁ΔT − k₂α₂L₂ΔT)/(k₁ + k₂)**. Every quantity is linear in the pair (P, ΔT), so the response to both is the sum of the responses to each: superposition.

(c) Compliance c = ∫₀ᵗ ds/(Eπ(u² − a²)) with u = D_w/2 + s tan α, a = d/2, du = tan α ds. With 1/(u² − a²) = (1/2a)(1/(u − a) − 1/(u + a)): c = (1/(Eπ tan α · 2a)) [ln((u − a)/(u + a))] from u₀ = D_w/2 to u₁ = D_w/2 + t tan α, i.e. c = (1/(πEd tan α)) ln[(D₁ − d)(D_w + d)/((D₁ + d)(D_w − d))] with D₁ = D_w + 2t tan α. **k = 1/c** is (7).

(d) From (8), F_m = F_b − P; into (9): (F_b − F_i)/k_b = (F_i − F_b + P)/k_m ⇒ F_b(1/k_b + 1/k_m) = F_i(1/k_b + 1/k_m) + P/k_m ⇒ **F_b = F_i + P k_b/(k_b + k_m) = F_i + CP**, and **F_m = F_i − (1 − C)P**. F_m = 0 at **P_sep = F_i/(1 − C)**. Thermal: the bolt and members are connected end to end through the head and nut, so any change of preload ΔF is common to both (same force: series); the free misfit δ_th = (Σαᵢtᵢ − α_b l)ΔT is absorbed as ΔF/k_b + ΔF/k_m, so **ΔF = δ_th k_b k_m/(k_b + k_m)**. Under an external load the head moves the same distance whether you follow the bolt or the members (same displacement: parallel), and the load divides between them.

**Check.** The tests compare (2), (5), (7), (10) and (12) with FE or quadrature to 10⁻⁸ or better.

**Common mistakes.** Sign errors in the junction equilibrium (draw the free-body diagram of the junction); treating bolt and members as parallel for the thermal case, which gives the sum instead of the series stiffness.

## 3. A liquid-oxygen feed line (s = 15: L = 2.0 m)

**Assumptions.** Uniform temperature; rigid mounts; no buckling of the line; constant typical 2024-T3 properties (E = 73.1 GPa, α = 23.2 × 10⁻⁶ K⁻¹, Fty = 345 MPa).

(a) A = π(50.8² − 48.3²)/4 = 194.58 mm². ΔT = −203 K. σ = −EαΔT = **+344.3 MPa (tension)**, a force of **66.99 kN**, practically equal to the typical yield strength (345 MPa). The free shrinkage would be **9.419 mm**.

(b) Line stiffness k_l = EA/L = **7112.0 N/mm**. Target stress 115.0 MPa ⇒ F = 22 377 N. Series springs: F = δ/(1/k_l + 1/k_bel) ⇒ 1/k_bel = δ/F − 1/k_l ⇒ **k_bel ≤ 3567.3 N/mm**. The line stretches F/k_l = 3.146 mm of the 9.419 mm, so the bellows must absorb **6.273 mm** of stroke.

(c) `bar_fe_solve([0, 2000, 2001], None, [7112.0 * 2000, 3567.3 * 1], [0, 0, 0], [0, 2], eps0=[23.2e-6 * -203, 0])` gives N = 22 377.0 N and a line stress of **115.00 MPa**.

(d) The expansion coefficient of aluminium falls substantially at cryogenic temperatures, so the real thermal strain over the 20 °C to −183 °C drop is smaller than α(20 °C) × 203 K: the stress in (a) is overestimated. The modulus rises somewhat at low temperature, which raises it. Aluminium alloys also get stronger in the cold, so the margin is larger than the room-temperature Fty suggests. A real design integrates α(T) over the temperature range.

**Check.** Rigid-bellows limit: k_bel → ∞ returns (a); soft limit k_bel → 0 gives zero stress and the bellows takes the full 9.419 mm.

**Common mistakes.** Treating the stress as compressive (cooling a restrained bar puts it in tension); forgetting that the bellows stroke is the free shrinkage minus the line's own stretch.

## 4. Can torque control meet the requirements? (seed 15)

(a) Thermal changes: −4003 N cold, +2668 N hot (from `preload_at`).

| torque (N·m) | separate below 18 kN when cold | bolt above Fty at limit when hot |
|---|---|---|
| 28.8 | 22.79 % | 0.73 % |
| 32.0 | 4.33 % | 3.46 % |
| 35.2 | 0.32 % | 11.48 % |
| 38.4 | 0.01 % | 28.30 % |

Raising the torque trades one failure for the other; no torque gets both below 1 %.

(b) Cold, clamped at ultimate: F_i − 4.00 kN − (1 − C) × 18 kN ≥ 0 ⇒ **F_i ≥ 16.45 kN**. Hot, below yield at limit: F_i + 2.67 kN + C × 12 kN ≤ 880 × 36.61 N ⇒ **F_i ≤ 25.85 kN**. The window is ±22 % about its centre (21.15 kN), while the 1st–99th percentile preload range under torque control is ±37 % or so: the distribution cannot fit, whichever way it is shifted. Shifting the torque slides the distribution along the window, which is why one fraction falls as the other rises.

(c) With a standard deviation of 0.01: **0.00 % and 0.00 %**; the preload range is **17.9–22.6 kN**, inside the window. A preload method with ±5 % scatter (one standard deviation) is needed: tightening to a measured bolt stretch (mechanical or ultrasonic), a load-indicating washer, or a torque-angle method that reads the stiffness of the joint as it tightens.

**Check.** The lesson's seed 507 gives 4.19 % and 3.55 % at 32 N·m; the differences between seeds (0.14 and 0.09 percentage points) are about twice the sampling standard error of 0.06 percentage points, plausible for two draws.

**Common mistakes.** Drawing a new random K for the hot and cold cases of the same joint (a joint has one preload; the temperatures change it deterministically).

## 5. A rotating compressor blade (s = 15: 15 000 rpm, ω = 1570.8 rad/s)

**Assumptions.** Uniform bar, radial direction only; no bending, no aerodynamic load; rigid disk.

(a) At radius r the outer part has centrifugal force N(r) = ∫_r^{r₀+L} ρAω²r′ dr′ = ρAω²((r₀ + L)² − r²)/2. At the root, r = r₀: N = ρAω²(r₀L + L²/2), and **σ_root = ρω²(r₀L + L²/2) = 4430 × 1570.8² × (0.25 × 0.1 + 0.005) = 327.9 MPa**, independent of A (doubling the area doubles the load too). Integrating N/(EA) from r₀ to r₀ + L gives **δ = ρω²(r₀L²/2 + L³/3)/E = 152.08 μm**.

(b) With consistent nodal loads:

| elements | tip stretch (μm) | root-element stress (MPa) |
|---|---|---|
| 2 | 152.0805 | 255.05 |
| 4 | 152.0805 | 292.62 |
| 8 | 152.0805 | 310.55 |
| 16 | 152.0805 | 319.31 |
| 32 | 152.0805 | 323.63 |

The tip stretch is exact for every mesh, to round-off: for 1-D bar elements with consistent loads the nodal displacements are exact (the element shape functions contain the exact solution's behaviour at the nodes, a property called nodal superconvergence). The root-element stress is the average stress over the first element, not the root value, and its error halves with each halving of h (72.9, 35.3, 17.4, 8.6, 4.3 MPa): **first order**.

(c) Free thermal growth αΔTL = 8.6 × 10⁻⁶ × 150 × 100 mm = **129.0 μm**; total tip-radius increase ≈ **281 μm** (152.1 + 129.0), before the disk's own growth, which is usually larger.

**Check.** Total load f.sum() = 16 395.9 N (for A = 50 mm²) equals ρAω²(r₀L + L²/2), and the FE root reaction balances it.

**Common mistakes.** Using r₀ alone instead of integrating over the blade; reading the first element's stress as the root stress and calling the 4 MPa difference an error in the FE solver.

## 6. Joint design decisions

(a) The separation load and the clamping force at ultimate both scale with F_i, so a high preload buys margin against separation, and keeps the bolt's alternating load at C times the external alternating load; preload is also lost to embedding and (for dissimilar metals) to cold, so the installed value must start well above the minimum needed. The ceiling is the bolt's yield, with the torque scatter and the hot case taken into account, which is why preloads well above half the yield strength are common.

(b) All-shank grip: k_b = 286.0 kN/mm, **C = 0.331** (against 0.308). A stiffer bolt takes a larger share of the alternating load, which is worse for bolt fatigue (the stiffness ratio is the subject of the Fastener Design Manual's "Fatigue-Resistant Bolts" section); and it increases the series stiffness in (12), so the cold preload loss grows (4.30 kN against 4.00 kN).

(c) At a flange near a corner or a free edge, on a thin lug, or where bolts are close together, the cone is cut off or overlaps its neighbour's: the effective member area is smaller, k_m falls and C rises. The model then overestimates the members' share and underestimates the bolt's alternating load.

(d) Measure the bolt's stretch (with a micrometer across the bolt ends, or ultrasonically from the time of flight along it), or use a strain-gauged bolt or a load-indicating washer. The preload is F_i = k_b × stretch; with the known torque the joint's own nut factor follows, K = T/(F_i d), and a set of such measurements gives its scatter, which is the input the Monte Carlo study needs.

## Check-yourself answers

1. δ = 10 000 × 1000/(205 000 × 50) = 0.976 mm.
2. Restrained means zero total strain, so σ/E = −αΔT: the stress depends only on E, α and ΔT; length and area cancel.
3. Compatibility: the deformations of the connected parts must fit together (here, the total length cannot change).
4. C = k_b/(k_b + k_m), the bolt's share of an external load. Stiffer plates raise k_m and lower C.
5. F_b = 20 + 0.3 × 10 = 23 kN; F_m = 20 − 0.7 × 10 = 13 kN; separation at 20/0.7 = 28.6 kN.
6. It falls: aluminium contracts more than titanium, so the plates shrink away from the bolt head and nut and the bolt relaxes.
7. While the joint is clamped, the bolt sees only C times the alternating load (about a third for the lesson joint); the members' unloading absorbs the rest.
8. A minimum clamp sets a lower bound on preload and a maximum stress an upper bound; the window between them for the lesson joint is ±22 %, while ±12.5 % per standard deviation spreads the preloads over about ±37 % (three standard deviations), so a few per cent of joints fall outside on each side.
