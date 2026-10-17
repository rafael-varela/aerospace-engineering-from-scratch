# P10.16  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. D = 0.40 + s/100; cosine grid on [0, 3]; r = `ellipsoid_radius`. 2. `body_volume`, (π/6)LD², `wetted_area`. 3. q = ½ρ55²; `normal_force_and_moment(x, r, α, q, 1.0)`; repeat at another x_ref. 4. `lamb_coefficients(L/D)`, `spheroid_munk_moment`; over-prediction = 1/(k₂ − k₁) − 1.

## 2.
(a) Momentum per unit length ρSw, w = Vα; force per length V d(ρSw)/dx. (b) ∫S′(x_ref − x)dx = [S(x_ref − x)] + ∫S dx; S(0) = 0; closed body S(L) = 0. (c) S = πR²x/L ⇒ S′ constant; N = 2qαπR²; centre of pressure ∫x dN/N. Cone: S′ ∝ x ⇒ centre at 2L/3. (d) Expand artanh e = e + e³/3 + e⁵/5 for small e; a₀, b₀ → 2/3, k → ½. For e → 1, artanh e → ∞ but (1 − e²) artanh e → 0: a₀ → 0, b₀ → 1.

## 3.
1. D = 15 m; Vol = (π/6)LD²; `lamb_coefficients(4)`. 2. q = ½ρ20²; M_sb = 2qα·Vol, M_ex = `spheroid_munk_moment`. 3. S_f = 2q·Vol(k₂ − k₁)/(q·3.0·25). 4. m_ex = k₂ρVol, m_sb = ρVol.

## 4.
1. L = 7 + s/20; grid; `body_volume`, `wetted_area`, f = L/1.1. 2. `fuselage_pitch_derivatives(x, r, 2.3)` and with `viscous_cut=True`; multiply the first by k₂ − k₁ at the new fineness. 3. `trainer_neutral_point(Cm, CN)` for both. 4. ΔC_D = C_f(55L/ν)·FF(f)·S_wet/12.

## 5.
1. `best_fineness(4.131, V)` for three V. 2. In a probe: `module.body_form_factor = lambda f: 1 + 1.5/f**1.5 + 7/f**3`; `best_fineness`, `spheroid_drag_area` at three f; restore the original. 3. Compare the optimum locations and the drag at 6.36.

## 6.
(a) Separation/thickening on the tail cone; wing upwash ahead of the wing. (b) `allen_perkins_cn` with and without η at 10° and 20° (S_base = π·0.55², A_plan from the profile); root of (total − 2 × potential). (c) Circular sections (the cross-flow of a flat-sided section has a larger apparent mass). (d) A fuselage-alone and wing–body pitching-moment test.
