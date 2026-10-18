# P06.01  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. v_r = ṙ, v_θ = rθ̇, V = √(v_r² + v_θ²). 2. a_r = r̈ − rθ̇², a_θ = rθ̈ + 2ṙθ̇; list r̈, rθ̇², rθ̈, 2ṙθ̇. 3. `cylindrical_to_cartesian(np.array([[v_r, v_θ, 0]]), θ)` and the same for a. 4. `path_kinematics(v, a)`: a_t, a_n, ρ; a_n/g; the sign of a_t; the sign of (v × a)_z.

## 2.
(a) ė_r = (−sin θ, cos θ)θ̇ = θ̇ e_θ (chain rule), similarly ė_θ; then differentiate r e_r and (2) by the product rule, substituting (1). (b) a_θ = 0 ⇒ d(r²θ̇)/dt = 0. With x = Vt, y = d: θ = atan(Vt/d), θ̇ = Vd/(d² + V²t²), maximum at t = 0; θ̈ = −2V³dt/(d² + V²t²)², set its derivative to zero ⇒ V²t² = d²/3; substitute. (c) v = (−Rω sin, Rω cos, w), a = (−Rω² cos, −Rω² sin, 0); |v × a| = Rω²√(R²ω² + w²); |v|³ = (R²ω² + w²)^{3/2}; simplify with b = w/ω; take the limits.

## 3.
(a) R = √(D² + z²), e = atan(z/D), differentiate e. (b) ė(t) = D a t/(D² + a²t⁴/4); dė/dt = 0 ⇒ D² = 3a²t⁴/4 ⇒ z = D/√3; then t, e and ė. (c) Ṙ = zż/R, R̈ = (ż² + z z̈)/R − (zż)²/R³, ë from the derivative of ė (or numerically); compare with a projected on e_R = (cos e, sin e) and e_e = (−sin e, cos e). (d) ë(t) from (b); look for its largest magnitude for t ≥ 0.

## 4.
1. `t = np.arange(0, 120.05, 0.1)`, `helix_turn(t)`, `accel_error_vs_stride(t, pos, acc, strides, 5, 1e-3, 5, seed=s)`. 2. Per-axis σ from the mean of `np.hypot(pos[:,0], pos[:,1])`; `optimal_step(σ, R ω⁴)`. 3. Repeat at 2 mrad; predicted factors from (9). 4. Repeat with `dt = 0.05` and strides doubled.

## 5.
(a) Lesson experiment 2 at h = 0.5 s, divided by g. (b) The best step from experiment 2; a central difference over h uses data from t − h to t + h. (c) A drift δθ̇ enters the centripetal term (2rθ̇ δθ̇ to first order) and 2ṙ δθ̇; an un-unwrapped azimuth jump or an overhead pass (r → 0) gives huge spurious rates. (d) An onboard inertial record of the same aircraft, time-aligned.
