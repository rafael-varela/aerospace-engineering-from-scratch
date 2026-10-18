# P10.08  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. θ_h = arccos(2E − 1). 2. τ, Δcl = 2πτδ, Δα_L0 = −τδ, Δcm = −½ sin θ_h (1 − cos θ_h) δ. 3. (4/π)√E and the ratio; scan E for the 2 % range with `flap_effectiveness`. 4. `flap_increments(E, δ)`; `thin_airfoil(flap_slope(E, δ), 0.0, m=40000)`.

## 2.
(a) A0 = α − (1/π)∫ slope dθ; Aₙ = (2/π)∫ slope cos nθ dθ over [θ_h, π]; cl and cm from A0, A1, A2. (b) d/dθ[sin θ (1 − cos θ)] = cos θ − cos 2θ = 0 → θ_h = 2π/3 → E from cos θ_h = 2E − 1. (c) θ_h = π − φ with cos φ = 1 − 2E, so φ ≈ 2√E; θ_h − sin θ_h ≈ π − 2φ + φ³/6. (d) E → 1: x_h → 0, the integral (3) becomes the leading-edge moment of the flat-plate load at α + δ.

## 3.
1. `hinge_moment(0.35, ±0.01, 0)` and `(0, ±0.01)` for the derivatives. 2. `hinge_moment(0.35, −2°, −10°)`. 3. q = ½ρV², c_f = 0.35 × 0.571, H = c_h q c_f² × 3.5. 4. Sign: positive is TE-down.

## 4.
1. `slat_configuration(slat_deflection=δ_s)`; `element_peaks(panel_basis(config), 10.0)`. 2. `critical_peak(naca4("2412", 100))`, `clmax_by_peak(config, cp_crit)`. 3. Repeat for 20° and 30°; tabulate.

## 5.
1. W = (700 + 5s) g. 2. Clean CL_max = 0.9 × 1.6; flapped with `wing_clmax(1.6, Δcl_max, 0.9, fraction)` where Δcl_max comes from `clmax_by_peak(slotted_flap_configuration(30.0), cp_crit, [0])`. 3. `stall_speed`; margin; needed CL_max = 2W/(ρS(30/1.3)²).

## 6.
(a) Separation behind the hinge; the inviscid peak at the hinge corner. (b) Fresh boundary layer, slot jet; but no calibration exists in the lesson. (c) TR-938 slotted-flap data (cl_max with and without flap) or a pressure-tapped test; compare the measured main-element peak at stall with cp_crit.
