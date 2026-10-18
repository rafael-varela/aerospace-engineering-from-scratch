# P10.23  Worked solutions

Each solution states assumptions, units, the derivation, the result, an independent check, and the common mistakes. Seeded problems are worked for s = 15. Numbers come from `solution.py` and a probe script that calls it. Every number derived from the generated datasets is **SYNTHETIC**: it comes from `make_synthetic_bench_data`, whose "truth" is a set of stated assumptions, not a measurement.

## 1. How small are the numbers?

(a), (b), (c) with ρ = 1.204 kg/m³, ρ_w = 998 kg/m³, S = 0.02 m²:

| V [m/s] | q [Pa] | water column [mm] | reading at 10° [mm] | Re | Re/2 × 10⁶ | lift at CL 0.5 [g] | drag at CD 0.05 [g] |
|---|---|---|---|---|---|---|---|
| 3 | 5.42 | 0.554 | 3.19 | 1.996 × 10⁴ | 1/100 | 5.52 | 0.55 |
| 4 | 9.63 | 0.985 | 5.67 | 2.661 × 10⁴ | 1/75 | 9.82 | 0.98 |
| 5 | 15.05 | 1.539 | 8.86 | 3.326 × 10⁴ | 1/60 | 15.34 | 1.53 |

(d) At 4 m/s: lift 9.8 steps of a 1 g scale and 98 of a 0.1 g scale; drag, through a ratio-2 bell-crank, 2.0 steps of a 1 g scale and 196 of a 0.01 g scale. **Drag cannot be measured usefully at the minimum-cost tier** (two resolution steps, and the sting's own drag is about the same size); lift can, to about 10 % at the larger angles.

**Common mistakes.** Forgetting that the inclined reading is longer than the column (that is the point of inclining it); using g = 9.81 to convert grams-force without the factor 1000.

## 2. The instruments

(a) The two legs are at pressures p₀ + q (pitot) and p₀ (static). The liquid surfaces differ in vertical height h, so q = (ρ_w − ρ)g h (the air column's weight is subtracted). Along a leg inclined at θ the meniscus moves Δl = h/sin θ, so **q = (ρ_w − ρ)g Δl sin θ**. Then ∂q/∂θ = q cot θ, and the incline's relative contribution is **(cot θ) u_θ = 5.67 × 0.00349 = 1.98 %** for θ = 10° ± 0.2°. A shallow incline stretches the reading (a 1 mm column becomes 5.7 mm), cutting the reading-resolution term; but cot θ grows as θ falls, so the incline must be measured more carefully (rise over run with a ruler gives ±0.1–0.2°).

(b) Moments about the pivot: D l_D = ΔF l_s, so **ΔF = D l_D/l_s = R D**. Calibrating with known horizontal loads measures R directly, including the real lever arms (where the drag acts, where the pad touches), pivot friction (shown by up/down hysteresis) and any flexibility of the arm, which a ruler measurement of l_D and l_s does not see.

(c) C = F/(qS): u_C² = (u_F/(qS))² + (C u_q/q)² + (C u_S/S)². The relative form divides by C, which is zero at zero lift; the absolute form stays finite and says correctly that at α = 0 the uncertainty is the force resolution divided by qS.

## 3. A small drone wing

(a) Re = 1.225 × 12 × 0.20/1.789 × 10⁻⁵ = **1.64 × 10⁵**: six times the bench, twelve times below Ladson's lowest.

(b)

| α | x_sep/c (Thwaites) | x_tr/c (Michel) | first |
|---|---|---|---|
| 0° | 0.626 | 0.817 | laminar separation |
| 2° | 0.451 | 0.664 | laminar separation |
| 4° | 0.282 | 0.490 | laminar separation |

No: even at the UAV's Re the laminar layer reaches separation before Michel's criterion is met (at 2 × 10⁶ transition comes first up to 4°). In reality the separated shear layer usually transitions and reattaches as a laminar separation bubble at this Re.

(c) McGhee et al., NASA TM-4062, the Eppler 387 at Reynolds numbers of 60 000 to 460 000 in the Langley Low-Turbulence Pressure Tunnel (`@ntrs-nasa-tm-4062-eppler-387-low-re`).

(d) Between the two: cd several times Ladson's but well below the bench's, a lift slope closer to the high-Re value, and drag that depends strongly on bubble behaviour, which varies from section to section. That is why small-UAV designers choose sections tested at their Re (the UIUC data, `P10.14`) rather than scaling NACA data down.

## 4. A synthetic dataset of your own (s = 15: seed 215; SYNTHETIC)

(a) SYNTHETIC results: q = **9.06 ± 0.29 Pa**, V = **3.91 ± 0.06 m/s**, Re = **2.55 × 10⁴**; R = **2.044 ± 0.003**; sting tares lift −0.03 g, drag **1.07 ± 0.05 g** on the drag scale; lift slope **0.0620 ± 0.0005 per degree (fit) ± 0.0021 (with q)**; zero-lift at α_set = **−0.600 ± 0.025°**; CL_max = **0.559 ± 0.028 at α_set = 10°**.

(b) Against the generator's stated truth: R = 2.04 (within 2u ✓); flow angle +0.6° in the jet, i.e. zero lift at α_set = −0.6° (✓); slope 0.0596 per degree at small angles (✓ with the q uncertainty; ✗ with the fit's u alone, 4u away). The fit's uncertainty sees only the scatter about the line; the manometer's q error scales every CL by the same factor and is invisible in the scatter. Here it is partly a planted bias: the incline was recorded as 10.0° while the generator used 10.15°, so q is underestimated by 1.5 % and CL overestimated. CL_max: truth 0.539 at 11.9° (true angle), measured 0.559 ± 0.028 at 10.6° true (✓).

(c) Against the inviscid finite-span band (0.0731–0.0822 per degree, midpoint 0.0777, u = 0.0026 as a uniform band): **E_n = 2.33**, a disagreement. Against Ladson's 0.1035 ± 0.0005 per degree at 2.0 × 10⁶: **E_n = 9.8**, a clear disagreement. Both say the same thing: the bench section does not behave as inviscid theory or high-Re data predict, which is the Reynolds-number effect the lesson is about (and, for this synthetic set, what the generator assumed).

## 5. The uncertainty budget (seed 230; SYNTHETIC)

(a) Relative standard uncertainties:

| α_set | CL | force term u_L/L | q term u_q/q | of which manometer | blockage band 2ε | unsteadiness (1.5 %/√2) |
|---|---|---|---|---|---|---|
| 0° | 0.037 | 12.7 % (absolute 0.005 in CL) | 4.2 % | 3.7 % | 1.7 % | 1.1 % |
| 4° | 0.267 | 2.4 % | 4.6 % | 3.7 % | 2.5 % | 1.1 % |
| 10° | 0.529 | 2.6 % | 5.3 % | 3.7 % | 3.7 % | 1.1 % |

CD at 0°: u_D/D = 4.0 %; at 4°: 2.6 %.

(b) For CL at 4° the q term dominates, and within it the manometer (3.7 %): a better airspeed measurement (a digital differential-pressure sensor, enhanced tier, or a longer averaging of more manometer readings with a better-measured incline) helps most. For CD at 0° the force term (4.0 %) and the q term (4.2 %) are comparable; a finer drag measurement and a better q help equally; the sting-drag tare (1.05 g of a 2.4 g reading) must be measured well, because it is subtracted.

(c) 1 g scales, two runs: lift resolution per reading pair √2 × 1/√12 g = 0.41 g → 0.0040 N, averaged over two runs 0.0028 N, so u(CL) from resolution = 0.0028/(qS = 0.190 N) = **0.015**, plus the q term (4.6 % of 0.27 = 0.012): **u(CL) ≈ 0.019** at 4°. Drag: the same 0.41 g on the scale divided by R = 2: **u(CD) ≈ 0.0075** from resolution alone, about 18 % of CD = 0.042. CD is barely measurable; the sting tare would be uncertain by the same amount.

## 6. Reynolds number or rig?

(a) Blockage: in a closed section it raises q at the model, so an uncorrected CL is too high; in an open jet the sign can reverse (lowering CL). Jet non-uniformity and swirl from the fan: the model sees a lower and varying local q, and a twist-like spanwise variation of α (lowers the slope). End-plate leakage: a finite effective aspect ratio (lowers the slope). Balance friction or a compliant sting that lets α change under load (either sign; lowers the slope if the model weathercocks nose-down). Flow angularity shifts the curve but not the slope.

(b) Run the same sweep at two or three fan speeds (for example 3, 4 and 5 m/s: Re = 2.0, 2.7 and 3.3 × 10⁴) without touching the rig. Rig effects that scale with q (blockage, end plates, non-uniformity) give the same CL(α) at every speed; a Reynolds-number effect makes the slope and CL_max rise and CD fall with speed. Turning the model upside down and repeating isolates the flow angle (the curves shift in opposite directions).

(c) In a closed section the walls constrain the flow around the model, which speeds it up (positive blockage). In an open jet the free boundary can expand, the jet deflects, and the correction is smaller and can have the opposite sign; the 1/4 rule is a closed-section estimate. Without a measurement of the jet's response, applying it could make the result worse; carrying ±ε on velocity (±2ε on q) is honest.

(d) "At Re ≈ 2.6 × 10⁴ the foam NACA 0012 section between end plates gave a lift slope of 0.057 ± 0.001 per degree (scatter only; ± 0.0024 including the airspeed uncertainty), about 55 % of the published value at 2 × 10⁶ (Ladson, NASA TM-4074) and 74 % of inviscid finite-span theory. A repeat at three fan speeds is needed to confirm that the difference is a Reynolds-number effect and not a rig effect; the blockage and jet non-uniformity are carried as uncertainties of up to 4 % on q." (If the data were the synthetic set, the first sentence begins "In the SYNTHETIC dataset…".)

## Check-yourself answers

1. About 9.6 Pa at 4 m/s: a water column of under 1 mm; inclining the tube at 10° stretches it to about 5.7 mm along the tube.
2. Re ≈ 2.7 × 10⁴ on the bench against 2.0–18.9 × 10⁶ in Ladson's tunnel: about 75 times lower than the lowest.
3. The lift reading falls by L/g (the model pushes up on the air and the air pushes up on the model); the drag scale rises by R D/g; tares are the same readings with the sting alone (and with the model at zero lift) and are subtracted.
4. At low Re the laminar boundary layer separates before it can become turbulent (laminar separation, sometimes a long bubble), so the effective shape loses camber and suction: lower lift slope, earlier and softer stall, and drag several times the high-Re value.
5. A written, dated prediction cannot be adjusted after seeing the data; it makes the comparison a test.
6. Fan guarded on both sides with intact finger-proof grilles, at most 50 W, on a stable surface; no propellers, drones or unguarded blades; no cutting tools near the running fan; hot-wire cutting only with an adult and ventilation, otherwise a craft knife with the fan unplugged; unplug the fan before adjusting the model.
7. The manometer's uncertainty of q multiplies every CL by the same factor; it does not add scatter about the fitted line, so the fit's standard error misses it.
