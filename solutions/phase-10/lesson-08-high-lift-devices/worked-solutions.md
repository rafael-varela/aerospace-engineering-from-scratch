# P10.08  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and a probe script that calls it, given to the figures shown.

## 1. An aileron (s = 15: E = 0.35, δ = 10°)

**Assumptions.** Thin-airfoil flap theory, sharp hinge, no gap.

(a) cos θ_h = 2(0.35) − 1 = −0.30, **θ_h = 107.46°** (1.8755 rad), sin θ_h = 0.9539. τ = 1 − (1.8755 − 0.9539)/π = **0.7067**. Δcl = 2π × 0.7067 × 0.17453 = **0.7749**; Δα_L0 = −0.7067 × 10° = **−7.067°**; Δcm_c/4 = −½ × 0.9539 × 1.30 × 0.17453 = **−0.1082**.

(b) (4/π)√0.35 = 0.7533, 6.6 % above τ. The approximation is within 2 % only for small flaps: τ(0.10) = 0.3958 against 0.4026 (+1.7 %), τ(0.05) = 0.2823 against 0.2847 (+0.8 %), so for E up to about 0.1.

(c) `flap_increments(0.35, 0.17453)` and the quadrature of `flap_slope` agree with (a) to 10⁻⁴.

**Common mistakes.** Using cos θ_h = 1 − 2E (that puts the hinge at x = E, from the leading edge).

## 2. Derivations

(a) The slope is 0 for θ < θ_h and −δ for θ > θ_h. A0 = α − (1/π)∫_{θ_h}^{π}(−δ) dθ = **α + δ(π − θ_h)/π**. Aₙ = (2/π)∫_{θ_h}^{π}(−δ) cos nθ dθ = (2δ/π) sin(nθ_h)/n. So cl = 2π(A0 + A1/2) = 2πα + 2δ(π − θ_h) + 2δ sin θ_h = 2πα + 2πδ[1 − (θ_h − sin θ_h)/π]: **Δcl = 2πτδ**, and zero lift at α = −τδ: **Δα_L0 = −τδ**. cm_c/4 = (π/4)(A2 − A1) = (δ/4)(sin 2θ_h − 2 sin θ_h) = (δ/2) sin θ_h (cos θ_h − 1): **Δcm_c/4 = −½ sin θ_h (1 − cos θ_h) δ**.

(b) g(θ) = sin θ (1 − cos θ), g′ = cos θ − cos²θ + sin²θ = cos θ − cos 2θ = 0 ⇒ cos θ = cos 2θ ⇒ θ = 2π/3 (the root in (0, π)). Then cos θ_h = −½ = 2E − 1 ⇒ **E = 0.25**, and g = (√3/2)(3/2) = 3√3/4, so |Δcm|/δ = **3√3/8 = 0.6495**. The code gives −0.6495 at 0.25 and smaller magnitudes either side (−0.6400 at 0.20, −0.6416 at 0.30).

(c) Let θ_h = π − φ. cos θ_h = −cos φ = 2E − 1 ⇒ cos φ = 1 − 2E ⇒ φ² ≈ 4E, φ ≈ 2√E. θ_h − sin θ_h = π − φ − sin φ ≈ π − 2φ + φ³/6, so τ = (2φ − φ³/6)/π ≈ 2φ/π = **(4/π)√E**. The test at E = 10⁻⁶ confirms it to 10⁻³ relative.

(d) As E → 1 the hinge moves to the leading edge (x_h → 0), the "flap" is the whole plate at angle α + δ, and (3) becomes −∫₀¹ Δcp x dx, the moment about the leading edge per q c². For a flat plate the load acts at the quarter chord with cl = 2π(α + δ), so the moment is −cl/4 and **c_h → −(π/2)(α + δ)**. The code at E = 0.999 agrees to 2 %.

**Common mistakes.** Dropping the n in the denominator of Aₙ; differentiating with respect to E instead of θ_h without the chain rule (same answer, more algebra).

## 3. The elevator's hinge moment

**Assumptions.** 2-D flat-plate tail, no gap, no balance, no downwash correction beyond the stated α_t.

(a) c_hα = **−0.687** /rad, c_hδ = **−0.988** /rad (E = 0.35). At α_t = −2°, δ = −10°: c_h = **+0.196** (≈ −0.687 × (−0.0349) − 0.988 × (−0.1745)).

(b) q = ½ × 1.225 × 30² = 551.3 Pa; c_f = 0.35 × 0.571 = 0.200 m; H = 0.196 × 551.3 × 0.200² × 3.5 = **15.2 N·m, trailing edge down**: the airflow tries to push the up-elevator back towards neutral, and the pilot (through the linkage) must hold it.

(c) Both derivatives negative means the surface floats into the wind and resists deflection: stable but heavy as q grows. A horn or a set-back hinge reduces |c_hδ| so stick forces stay manageable. If c_hδ became positive (over-balance), a small deflection would grow on its own: the control would snatch.

**Common mistakes.** Using c (tail chord) instead of c_f in H; forgetting the span.

## 4. Slat deflection (s = 15: δ_s = 25°)

(a) α = 10°:

| δ_s | main peak cp | slat peak cp | total cl |
|---|---|---|---|
| 20° | −1.486 | −7.013 | 1.549 |
| 25° | −1.586 | −4.626 | 1.535 |
| 30° | −1.711 | −3.022 | 1.521 |

(b) cl_max by the criterion: δ_s = 25°: **cl_max = 1.728 at α = 11.53°, limited by the slat** (20°: 1.564 at 10.12°; 30°: 1.916 at 13.15°; slat-limited in all three).

(c) Turning the slat nose-down lowers its own incidence, so its peak falls (−7.0 → −3.0 at 10°) and it reaches the critical value at a higher α; the main element's peak rises a little (the slat carries less load, so it shields less) but stays far below critical. The best deflection balances the two peaks; here the slat limits throughout, so more deflection would help until the main element took over or the geometry overlapped.

**Common mistakes.** Rotating the slat about its leading edge (it changes the gap); using the main element's peak alone.

## 5. Your own trainer (s = 15: m = 775 kg, flapped fraction 0.65)

(a) W = 7602.8 N. Clean: CL_max = 1.440, **V_s = 26.80 m/s, margin 1.119**. Slotted flap: Δcl_max = 1.272, CL_max = 0.9(1.6 + 0.65 × 1.272) = **2.184, V_s = 21.76 m/s, margin 1.379**.

(b) Needed CL_max = 2 × 7602.8/(1.225 × 12 × (30/1.3)²) = **1.942**. The flapped wing meets it (2.184); the clean wing does not.

(c) With flaps, 1.3 × 21.8 ≈ **28.3 m/s** is the minimum; 30 m/s keeps 1.38 V_s and is a sensible choice. Test first the clean cl_max (the calibration of the whole criterion) and the 3-D factor 0.9, best by measuring the stall speed in flight.

**Check.** V_s ∝ √m: √(775/750) × 21.69 m/s would be 22.05 m/s with the lesson's flapped fraction; the larger flapped fraction brings it to 21.76 ✓.

**Common mistakes.** Using the mass in kg in W; applying the 0.9 factor twice.

## 6. Trusting the cl_max criterion

(a) The criterion only watches the *leading-edge* peak. Behind a sharp hinge at 30–40° the inviscid flow turns a corner with a singular suction peak followed by a steep adverse gradient; a real boundary layer separates there, so the flap loses effectiveness and the section stalls from the hinge, not from the nose. The model does not see the separation, so its plain-flap Δcl_max keeps rising with δ.

(b) Defence: the flap's boundary layer starts at the flap's own stagnation point in the high-energy slot flow (Smith's fresh-boundary-layer effect), and the critical peak calibrated on the main element's nose at the main chord's Reynolds number does not apply to it. Attack: the lesson has no alternative criterion for the flap, so its exclusion is unverified; a 30° flap can and does stall. The honest position: exclude it, report its peak, and flag the result as uncalibrated for flap stall.

(c) Use TR-938's single-slotted-flap section data: for a configuration close to the model's, compare the measured cl_max increment with the criterion's, and fit a separate critical peak for the flap. Better still, a pressure-tapped two-element test, comparing the measured main-element and flap peaks at the measured stall angle with the panel peaks.

## Check-yourself answers

1. cos θ_h = 2(0.3) − 1 = −0.4, θ_h = 113.6°.
2. τ = −dα_L0/dδ = 1 − (θ_h − sin θ_h)/π: the angle-of-attack change one radian of flap is worth. τ → (4/π)√E for a small flap and τ → 1 for a full-chord flap.
3. Its extra load sits aft of the quarter chord. The effect per degree is largest for a 25 % flap.
4. Both negative: lift behind the hinge pushes the trailing edge up, so a free flap floats trailing edge up.
5. Slat effect: the slat's circulation lowers the velocity at the main element's nose and so its suction peak. Circulation effect: a downstream element raises the velocity at the upstream element's trailing edge and so its circulation.
6. It unloads the main element as much as it lifts itself, but by lowering the main element's suction peak it lets the main element reach a higher angle before separating.
7. V_s ∝ 1/√CL_max: it falls to 1/√1.44 = 1/1.2 of its value, about 17 % lower.
8. Any two of: separation on the flap or behind a hinge; the slat's lower Reynolds number; confluent boundary layers and off-surface recovery; the assumed clean cl_max.
