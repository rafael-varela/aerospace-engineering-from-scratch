# P05.08  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. T = `shaft_torque(P, 6000)`. 2. J = `polar_moment(50, 47)`. 3. τ_o = T × 25/J, τ_i = T × 23.5/J. 4. φ = `twist_angle(T, 2000, G, J)` in degrees. 5. Ratio Fsy/τ_o. 6. r_m = 24.25 mm; τ = T/(2πr_m²t); compare.

## 2.
(a) γ = r dφ/dx (compatibility); τ = Gγ (Hooke); T = ∫τr dA (equilibrium) ⇒ dφ/dx = T/GJ. (b) Axial force balance on a wall element ⇒ q constant; ∮q p ds = 2A_m q = T. Work ½T dθ/dx per unit length = ∮(τ²/2G)t ds. (c) r_o⁴ − r_i⁴ = 4r_m³t + r_m t³ ⇒ J_exact = 2πr_m³t(1 + t²/(4r_m²)). (d) J_open = 2πr_m t³/3; ratio J_closed/J_open.

## 3.
(a) Vertices R(cos u, sin u) for u from π/2 to 3π/2 (n + 1 points); thicknesses: n × 1.2 for the nose, 2.0 for the closing web; `bredt_batho(T, verts, t, G)`. (b) A_m = πR²/2; J = 4A_m²/(πR/1.2 + 2R/2.0). (c) `open_section_J([πR, 2R], [1.2, 2.0])`; ratio.

## 4.
(a) T_rated = `shaft_torque(P, 20900)`; limit × 1.25; ultimate × 1.5. For k = 0 and 0.8: `size_hollow_shaft` with Fsy at limit and Fsu at ultimate; `size_for_twist` with 5°; take the maximum. (b) φ_x = 32 T_lim L/(πG d_u⁴(1 − k⁴)). (c) Compare φ_x with the requirement.

## 5.
(a) `prandtl_fd(10, 1, h)` for three h; error against `rectangle_torsion_constant(10, 1)`; order log₂ of error ratios. (b) `prandtl_fd(2, 1, t/k, hole=(1.8, 0.8))` for k = 2, 4, 8, 16; order from successive differences (no exact value). (c) `bredt_batho(1, [[0,0],[1.9,0],[1.9,0.9],[0,0.9]], 0.1, 1)`; ratio.

## 6.
(a) Transient over-torque in acceleration, governor overshoot, a stalled rotor. (b) Torsional buckling of a thin tube; whirl (bending critical speed) between bearings. (c) The bay becomes an open section locally; ring doubler or frame closes it. (d) Stress concentration × alternating stress at the spline root.
