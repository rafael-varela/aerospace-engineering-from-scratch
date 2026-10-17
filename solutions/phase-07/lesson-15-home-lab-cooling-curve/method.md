# P07.15  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. `bar_geometry(0.0254, L)`. 2. ρcV with Table A.1. 3. `predict_tau(14.0, 294.15, L=L)`. 4. Compare V/A = DL/(4L + 2D) for the two lengths.

## 2.
(a) First law: ρcV dT/dt = −hA(T − T∞); separate and integrate. Conduction resistance L_c/(kA) in series with 1/(hA): the internal drop is the fraction (L_c/k)/(1/h) = Bi of the external one. (b) Substitute the trial solution; a = θ₀τ/(τ − τ_p); b = −a; for t ≫ τ_p the b term is gone. (c) Write T = T_m + d, T∞ = T_m − d; expand (T² + T∞²)(T + T∞) = 4T_m³ + 4T_m d². Evaluate both with `h_rad_linear` and the 4εσT_m³ formula. (d) J₀ ≈ 1 − λ²/4, J₁ ≈ λ/2 ⇒ λ²/2 ≈ Bi; C₁ → 1, other terms decay fast.

## 3.
(a) h_rad = `h_rad_linear(0.09, T∞ + 14, T∞)`; τ = ρcV/(h_rad A). (b) `quad` of ρcV/(εσA(T⁴ − T_w⁴)) from 303.15 to 323.15 K; linearised τ ln((50 − 21)/(30 − 21)). (c) Repeat with ε = 0.9. (d) With T_w = 0: dt = −ρcV dT/(εσAT⁴), integrate.

## 4.
1. Generate, convert to kelvin, room mean from the ambient rows. 2. The three fits as in `main()`. 3. `predict_tau_uncertainty` at the measured T∞; `en_number`. 4. Windows by masks on t; compare the trend with `predict_tau` at each window's mid excess. 5. `fit_variable_h`; compare with `truth`; explain the f–ε trade-off by what each does to h.

## 5.
1. For each material and h: `biot_number(h, k, Lc)` and `time_constant`. 2. α = k/(ρc); `short_cylinder(0, 0, [τ], ...)` and `short_cylinder(1, 1, [τ], ...)`. 3. The warm-up is the water-bath case; a few τ of the water-bath case brings the bar to the water temperature.

## 6.
(a) `predict_tau(14, 294.15, f_conv=1.3)` and `(eps=0.3)`; compare with 0.8 × 1968 s; look at the windowed τ and at a black-tape comparison. (b) `fit_loglinear` with T∞ + 1 K. (c) m_p c_p/(ρcV); a fin draws heat out of the bar. (d) Rank by size and plausibility; one test each.
