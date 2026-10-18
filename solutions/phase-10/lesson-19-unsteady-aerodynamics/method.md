# P10.19  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. k = 2πf (c/2)/V for (a) and (b). 2. For (c): Ω = V_tip/R, V = 0.75 V_tip, b = c/2, k = Ωb/V. 3. `theodorsen(k)`: |C| and `np.angle`. 4. For (c), `theodorsen_lift(k, 1, 0, -0.5)`, add the two parts, divide the magnitude by 2π.

## 2.
(a) Write cl_c(s) = 2π[α₃/₄(0)φ(s) + ∫₀ˢ φ(s − σ) α₃/₄′(σ) dσ]; start the input at −∞ so the first term disappears; substitute τ = s − σ. (b) Transform each term; multiply by p = ik; take the limits. (c) Integrate ik∫φe^{−ikτ}dτ by parts with u = φ, dv = ik e^{−ikτ}dτ. (d) Set 1 − ½e^{−βs} = 0.9; solve for s; then for β given s = 12.7.

## 3.
1. k from Problem 1(c). 2. Quasi-steady: `quasi_steady_lift(k, α̂, 0, -0.5)`; Theodorsen: sum of `theodorsen_lift`. 3. Phase of each relative to α̂, in degrees (= degrees of azimuth). 4. Argue (b) and (c) from what a phase error does in a feedback loop and from the assumptions list.

## 4.
1. k and a from your seed. 2. `theodorsen_lift` and `quasi_steady_lift` at α̂ = 1. 3. `simulate_harmonic_pitch(k, radians(2), a)["total"]`; compare with the exact and the Jones-based C. 4. `lift_loop` and `loop_area`; `aero_work_per_cycle(-1j α̂, ĉl)`. 5. `aero_work_per_cycle(1, theodorsen_moment(k, 1, 0, a))`.

## 5.
1. Gust shape w(s) = (w₀/2)(1 − cos(πs/H)) for 0 ≤ s ≤ 2H (s in semi-chords: H in chords times 2). 2. `gust_response`, take the maximum, divide by 2π w₀/V. 3. Δn = ρ V CL_α w₀/(2 W/S) times the ratio.

## 6.
(a) Loop `pitch_flutter_boundary(a)`; x_pivot/c = (1 + a)/2. (b) Write cm about the pivot as apparent-mass terms plus (a + ½)/2 × circulatory cl; for a ≥ −½ show the circulatory work is negative or zero. (c) Combine the two single-degree-of-freedom results with the role of coupling. (d) Name dynamic stall.
