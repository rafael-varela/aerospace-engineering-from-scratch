# P06.12  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. ω₀ = √(g/l); T₀ = 2π/ω₀; E = −ω₀² cos a; k² from its definition. 2. `pendulum_period(E, w0)`; the series; central difference of `phase_area` with dE = 10⁻⁶ω₀². 3. `phase_area(E, w0)`; p_max = √(2(E + ω₀²)); separatrix 2ω₀; speed 2ω₀l.

## 2.
(a) dH = Σ(q̇ dp + p dq̇) − Σ(∂L/∂q dq + ∂L/∂q̇ dq̇); cancel the dq̇ terms with (1); replace ∂L/∂q by ṗ (Lagrange); compare with dH = Σ(∂H/∂q dq + ∂H/∂p dp). (b) A(E) = 2∫ p(q, E) dq between turning points; differentiate (the boundary terms vanish because p = 0 there); ∂p/∂E = 1/p; dq/p = dt. (c) div v = ∂q̇/∂q + ∂ṗ/∂p = H_pq − H_qp = 0; with damping div v = −c; Reynolds transport: dA/dt = ∫ div v dA. (d) Write each Jacobian; multiply; compare determinants.

## 3.
(a) `legendre_transform(L, [th], [thd], [p])`; compare with T + V; argue from ∂L/∂t and from the work done by the hoop's drive. (b) ∂H/∂θ = 0 with p = 0: sin θ(g/R − Ω² cos θ) = 0; stability from the sign of ∂²H/∂θ² at p = 0. (c) Equilibria θ = 0, π and cos θ₀ = g/(RΩ²) if Ω > Ω_c; frequencies √(H_θθ/(mR²)). (d) `solve_ivp` of θ̈ = Ω² sin θ cos θ − (g/R) sin θ; evaluate both quantities per unit mass along the run.

## 4.
(a) `pendulum_hamiltonian` at the centre; `polygon_area`; compare the largest energy on the boundary with ω₀². (b) `flow_blob` for 5 s; ratio. (c) Loops over steppers and h with `iterate`. (d) Ratios of successive RK4 changes; log₂.

## 5.
(a) δa = 2a²v_pδv/μ; exact: energy of the perturbed state, a₁ = −μ/(2E₁). (b) `along_track_drift(a, δa, n, e)`. (c) Two runs; d · v/|v|. (d) `radial_action(E, h)`; central difference in E.

## 6.
(a) Separate volume from the spread along one direction. (b) Identify the non-symplectic step; propose a symplectic or much tighter integration and a test. (c) Compare the size of a structure-preserving benefit over one orbit with a plain accuracy requirement. (d) Name the force, say what it does to a and therefore to the phase drift.
