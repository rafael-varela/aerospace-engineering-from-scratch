# P05.18  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe scripts that call it, given to the figures shown. Material values are typical handbook values. Results from `make_synthetic_cantilever_data` are SYNTHETIC, not measurements.

## 1. Predict before you load (s = 15: m = 125 g)

**Assumptions.** Rigid clamp, small deflection, E = 69 ± 2 GPa (typical), g = 9.81 m/s².

(a) I = 0.020 × 0.002³/12 = **1.333 × 10⁻¹¹ m⁴ = 13.33 mm⁴**; **EI = 0.920 N·m²**.

(b) P = 0.125 × 9.81 = 1.226 N. At the load: Pa³/(3EI) = 1.226 × 0.015625/2.760 = **6.94 mm**. At x = 150 mm: P x²(3a − x)/(6EI) = **3.00 mm**. At x = 300 mm (on the straight tail): P a²(3x − a)/(6EI) = **9.02 mm**.

(c) σ = 6Pa/(bh²) = 6 × 1.226 × 0.25/(0.020 × 4 × 10⁻⁶) = **23.0 MPa**. Safe load: P_max = 72.5 × 10⁶ × 0.020 × 4 × 10⁻⁶/(6 × 0.25) = 3.867 N = **394 g**. 125 g is safe.

(d) Terms: E 2.90 %, b 0.25 %, h 3 × 0.02/2.00 = 3.00 %, a 3 × 1/250 = 1.20 %. Total √(8.40 + 0.06 + 9.00 + 1.44) = **4.35 %, δ = 6.94 ± 0.30 mm**. `predict_deflection` gives the same, 6.942 ± 0.302 mm. Here the thickness term is the largest: a 2 mm strip is more sensitive to a 0.02 mm error than a 3 mm one.

(e) q = 2700 × 9.81 × 0.020 × 0.002 = 1.059 N/m; w_q(a) with L = 0.3 m is **0.91 mm**. It is present in both the unloaded and the loaded reading, so their difference does not contain it.

**Check.** Units: N × m³/(N·m²) = m. The tail is straight: (9.02 − 6.94)/50 mm = 0.0416, equal to the slope at the load Pa²/(2EI) = 0.0417.

**Common mistakes.** Using h = 2 mm as 0.02 m; forgetting the tail beyond the load and applying the cubic at x = 300 mm (which gives a wrong, smaller value); quoting the safe load in newtons as grams.

## 2. Derivations

(a) E I w′ = P(ax − x²/2) + C₁, E I w = P(ax²/2 − x³/6) + C₁x + C₂; w′(0) = 0 gives C₁ = 0, w(0) = 0 gives C₂ = 0. **w = P x²(3a − x)/(6EI)** for x ≤ a. At a: w(a) = Pa³/(3EI), w′(a) = Pa²/(2EI). For x > a, w″ = 0, so w = w(a) + w′(a)(x − a) = **P a²(3x − a)/(6EI)**, a line of slope Pa²/(2EI).

(b) ln δ = ln 4 + ln P + 3 ln a − ln E − ln b − 3 ln h, so dδ/δ = dP/P + 3 da/a − dE/E − db/b − 3 dh/h, and for independent inputs **(u_δ/δ)² = (u_E/E)² + (u_b/b)² + (3u_h/h)² + (3u_a/a)²** (+ (u_P/P)²). At a fixed mark x ≤ a: ln g = 2 ln x + ln(3a − x) − ln 6, so d ln g/d ln a = a × 3/(3a − x) = **3a/(3a − x)**, which is **1.5 at x = a**. It is smaller than 3 because two of the three powers of a in a³/3 come from the reading point moving outward with the load; at a fixed mark only the moment arm changes. (The mark's own position then enters with d ln g/d ln x = 2 − x/(3a − x) = 1.5 at x = a.)

(c) w(a) = Pa³/(3EI) + θ₀a = Pa³/(3EI) + Pa²/k_θ. Dividing by P: **s = a³/(3EI) (1 + 3EI/(k_θ a))**. The rigid-clamp formula gives E_apparent = a³/(3sI) = **E/(1 + 3EI/(k_θ a))**. For EI = 3.881 N·m², a = 0.29 m, k_θ = 2000 N·m/rad: 3EI/(k_θ a) = **0.0201**, so the clamp makes E look **2.0 % low**. This matches the SYNTHETIC generator: 0.042 of 2.118 mm/N.

(d) Readings: unloaded r₀ + w_q(x), loaded r₀ + w_q(x) + w_P(x), where r₀ is the ruler's offset. The difference is **w_P(x) exactly**, whatever w_q looks like. This relies on **linearity** (small deflections, linear elastic material): the response to self-weight plus P equals the sum of the responses. At large deflection it holds only approximately.

**Common mistakes.** In (b), treating the reading at a fixed mark as if it moved with a; in (c), adding θ₀ instead of θ₀a to the deflection.

## 3. A thin steel rule (s = 15: m = 50 g)

**Assumptions.** E = 193 GPa and F_ty = 215 MPa (typical, annealed 304); vertical load; rigid clamp.

(a) I = 0.025 × 0.0008³/12 = 1.067 mm⁴, EI = 0.2059 N·m², P = 0.4905 N. **Linear 10.98 mm, elastica 10.96 mm, linear high by 0.21 %, δ/a = 0.046, σ = 44.1 MPa.**

(b) P_max = 107.5 × 10⁶ × 0.025 × 0.0008²/(6 × 0.24) = 1.194 N = **122 g**. There: **linear 26.74 mm, elastica 26.40 mm, linear high by 1.26 %, δ/a = 0.110**. Linear theory is 1 % high at **108 g**, where the linear δ/a is **0.099**, the same ≈ 0.10 as the aluminium strip of the lesson (whose 1 % point is at 1398 g, computed only). The agreement is exact, not a coincidence: the ratio depends only on P a²/(EI).

(c) √(5.0² + 0.2² + 3.75² + 1.25²) = **6.4 %**. Against the prediction, a 1.3 % effect is invisible. But E, h and b are scale factors: they multiply every deflection by the same number and cannot bend the line. The nonlinearity bends it: at the safe load the elastica is **0.33 mm** below the straight line through the low-load points. With a ruler (0.42 mm on each difference) that is below one standard uncertainty for a single point, so with a ruler it is not detectable; with a dial gauge (0.014 mm) it is twenty-four standard uncertainties and clearly detectable as curvature.

(d) **Safety glasses** (a thin hardened rule can snap or slip from the clamp and spring), the cut or worn end taped, and the clamp re-checked before each load because a thin rule slips more easily.

**Check.** At 50 g the elastica result is within 0.3 % of linear, as it must be at δ/a = 0.05.

**Common mistakes.** Taking E for spring steel as 200 GPa and the yield as that of hardened steel (unsafe if the rule is annealed); expecting a 1 % effect to show against a 6 % prediction band.

## 4. Reduce a SYNTHETIC dataset (s = 15)

All results SYNTHETIC.

(a) **Slope 2.185 ± 0.056 mm/N, intercept +0.10 ± 0.26 mm, χ² = 7.2 for 7 dof. E = 66.1 ± 2.3 GPa; E_n = 0.47** against 69 ± 2 GPa.

(b) **EI = 3.934 ± 0.408 N·m², θ₀ = 0.71 ± 4.70 mrad, χ² = 1.6 for 5 dof.**

(c) EI_true = 67.5 × 10⁹ × 0.02502 × 0.00301³/12 = 3.838 N·m². True compliance 2.118 + 0.042 = **2.160 mm/N**; apparent E through the nominal dimensions **66.9 GPa**. Your 66.1 ± 2.3 GPa differs from it by 0.33 standard uncertainties: consistent. (The seed-518 fit in the lesson, 69.2 ± 2.4 GPa, is also consistent, 1 u away on the other side.)

(d) **u(θ₀) = 0.16 mrad.** The generator's θ₀ at 700 g is 0.7 × 9.81 × 0.29/2000 = **1.00 mrad**, about six standard uncertainties: detected. With the ruler's 4.7 mrad it was not.

**Check.** χ²/dof near 1 in (a) and (b) shows the assumed 0.42 mm per difference matches the generator's noise, as it should.

**Common mistakes.** Using 0.30 mm instead of √2 × 0.30 mm for a difference of two readings; calling the SYNTHETIC result a measurement.

## 5. A 1 % measurement of E

(a)

| term | value |
|---|---|
| slope (ruler) | 2.67 % |
| thickness 3u_h/h, u_h = 0.02 mm | 2.00 % |
| load position 3u_a/a, u_a = 1 mm | 1.03 % |
| width u_b/b, u_b = 0.05 mm | 0.20 % |
| **total** | **3.50 %** |

(b) With the dial gauge the slope term is 0.09 % and the total 2.26 %, dominated by thickness. **u_h = 0.005 mm** (a micrometer, or more caliper readings along the strip) gives 0.50 %; with u_a still 1 mm the total is 1.17 %, so the load position must improve too: **u_a = 0.5 mm** (a filed notch, a ruler read under load) gives 0.52 %. Total **0.75 %**.

(c) The 2 % clamp bias is **outside** the 1 % target and twice the whole budget: it must be removed, not ignored. Fit the profile at one load with the dial gauge: u(θ₀) = 0.16 mrad, so the correction θ₀a carries 0.16 × 10⁻³ × 0.29 = 0.046 mm on δ = 14.4 mm at 700 g, **0.32 %**. Subtract θ₀a from each load-point deflection (scaled by P), or use k_θ from the fit; the total becomes √(0.75² + 0.32²) = **0.81 %**. Better still, clamp onto a steel block or between two steel bars to make k_θ much larger.

(d) At a fixed mark at x = a, the load position enters with sensitivity 1.5 (0.52 % for u_a = 1 mm), but the mark's own position adds another 1.5 u_x/x: with both at 1 mm, √2 × 0.52 = 0.73 % against 1.03 %. A modest gain, and only if the mark is located well; a dial gauge on a fixed stand is in effect a fixed mark.

**Common mistakes.** Improving the width measurement (0.2 %) first; treating the clamp bias as a random term and adding it in quadrature.

## 6. Judgment: a surprising modulus

A model review, in outline. **E_apparent/E = 56/69 = 0.812.** **Candidate 1, load position measured from the wrong reference.** E scales as (a_true/a_used)³; the result needs a_true = **270.5 mm**, an error of 19.5 mm. They measured "from the end of the strip under the clamp" and about 20 mm is under the clamp, so a = 290 mm from the end is 270 mm from the clamp edge: E_apparent = 69 × (270/290)³ = **55.7 GPa**. This explains the whole result with an error the description makes almost certain. **Candidate 2, thickness.** E scales as (h_true/h_used)³; it needs h = **2.80 mm** against 3.00 used, a 7 % error. Possible only if the nominal "3 mm" from the label was used without a caliper, and then only if the bar is at the bottom of its tolerance. **Candidate 3, a compliant clamp.** It needs 3EI/(k_θ a) = 0.232, k_θ = **173 N·m/rad**, a clamp that lets the strip rotate 11.5 mrad at 700 g (3.3 mm at the load point). Possible with a loose clamp or a soft table; it would be visible as the strip lifting at the clamp's edge. A wrong alloy cannot do it: aluminium alloys all have E within a few per cent of 69 GPa. **Cheapest tests.** Re-measure a from the clamp edge (no cost, one minute): if it reads about 270 mm, cause 1 is confirmed, and recomputing gives E near 69 GPa. Measure h with a caliper at five points: confirms or excludes cause 2. Then a profile fit for θ₀: 11.5 mrad against a ruler's 4.7 mrad resolution is 2.4 u, marginal; a dial gauge at the clamp edge settles it. **Verdict.** The conclusion is not acceptable. E_n = 2.3 is a disagreement that demands a cause, and the measurement description already contains one that fits; "the handbook is wrong" is the last explanation to reach for, after the measurement's own references have been checked. The report should state the corrected E, or withdraw the result until the test is done.

## Check-yourself answers

1. **δ = Pa³/(3EI)**. Doubling h multiplies I by 8, so δ falls to **1/8**.
2. It enters cubed (I = bh³/12), so its relative uncertainty counts three times: **3 × 0.02/3.00 = 2.0 %**.
3. σ = 6Pa/(bh²) = 6 × 3.924 × 0.29/(0.025 × 9 × 10⁻⁶) = **30.3 MPa**.
4. The difference cancels the self-weight sag, the ruler's zero offset and any slow drift of the clamp, table or camera between loads (superposition, and the two readings taken seconds apart).
5. It adds **θ₀a = Pa²/k_θ** at the load point. That is proportional to P, like the beam term, so the line's slope rises and reads as a lower E; only the shape along the strip (x against x²(3a − x)) separates them.
6. At **δ/a ≈ 0.10**; linear theory over-predicts (the true deflection is smaller, because the load's lever arm shortens as the strip bends).
7. E_n = 0.2/(2√(2.4² + 2²)) = **0.03**: consistent. You may conclude only that the measurement does not contradict the typical value at this uncertainty, not that E is 69.2 GPa or that the clamp is rigid.
8. So that the comparison is a test: a prediction written after seeing the data can be adjusted, consciously or not, to agree with it. The date shows the order.
