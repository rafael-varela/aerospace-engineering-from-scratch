# P10.05  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15. Numbers come from `solution.py` and the probe scripts that call it, given to the figures shown.

## 1. One panel by hand

**Assumptions.** Panel along the x-axis, so panel and global axes coincide; clockwise sheet γ = 1.

(a) For a uniform sheet, u = β/(2π) and v = −ln(r₁/r₂)/(2π).
- (0.5, 0.5): the ends are seen at 45° and 135°, β = π/2, r₁ = r₂: **(u, v) = (0.25, 0)**.
- (0.5, −0.5): β = −π/2: **(−0.25, 0)**.
- (2, 0): β = 0, r₁ = 2, r₂ = 1: **(0, −ln 2/(2π)) = (0, −0.1103)**.

(b) A point vortex Γ = 1 at (0.5, 0) gives (0.3183, 0), (−0.3183, 0) and (0, −0.1061). The difference is largest close to the panel (0.25 against 0.318 at a distance of half the panel length): near a panel its length matters; far away it looks like a point vortex (a test checks this at 360 panel lengths to 10⁻⁵).

(c) `linear_vortex_influence` returns the same values. At (0.5, 0.5) the a part (strength falling from 1 to 0) gives (0.125, −0.0342) and the b part (rising) gives (0.125, +0.0342). The point is above the midpoint, so the two halves contribute equally to u, the velocity along the panel; the a part has more vorticity on the left, whose clockwise swirl pushes the point down, and the b part the mirror image.

**Common mistakes.** Reversing the sign of v for a clockwise vortex; measuring β as the angle of the point instead of the angle subtended.

## 2. The influence integrals and the system

(a) With s = (s − ξ) + ξ: ∫ s η/D ds = η∫(s − ξ)/D ds + ξ∫η/D ds = η · ½ln D|₀^L + ξβ = **η ln(r₂/r₁) + ξβ**. Similarly ∫ s(ξ − s)/D ds = −∫(s − ξ)²/D ds + ξ∫(ξ − s)/D ds = −∫(1 − η²/D) ds + ξ ln(r₁/r₂) = **−L + ηβ + ξ ln(r₁/r₂)**.

(b) A clockwise element γ ds at (s, 0) induces du = γη/(2πD), dv = −γ(ξ − s)/(2πD). Substituting γ = g_a(1 − s/L) + g_b s/L: u = [g_a(I₁ − I₃/L) + g_b I₃/L]/2π and v = −[g_a(I₂ − I₄/L) + g_b I₄/L]/2π. On the panel (η = 0, 0 < ξ < L) the subtended angle jumps from +π (above) to −π (below): the tangential velocity jumps by the sheet strength, so a collocation point on the panel must say which side it is on. The code takes the outside (η → 0⁻ for a counter-clockwise surface).

(c) N panels have N + 1 node strengths but only N midpoints for flow tangency: one equation short. Any multiple of the "circulation-only" solution could be added; the Kutta condition picks one. γ₀ and γ_N are minus the surface velocities along the node direction at the upper and lower trailing-edge nodes. Those directions point forward on the upper surface and aft on the lower, so equal speeds leaving the edge mean equal and opposite γ: **γ₀ + γ_N = 0**.

(d) Across a vortex sheet the tangential velocity jumps by γ. If the fluid inside is at rest, the outside velocity is the whole jump: |u_t| = |γ| (u_t = −γ along the node direction for clockwise γ). The circulation is ∮γ ds; with γ linear on each panel the trapezoid rule is exact: Γ = Σ L_j(γ_j + γ_{j+1})/2, and c_l = ρVΓ/(½ρV²c) = **2Γ/(Vc)**.

**Common mistakes.** Using the inside limit for the self-influence; counting the closed trailing edge as one node and losing the Kutta row.

## 3. A "Kutta condition" on a circle

(a) Unit circle, polygon from θ = 0 round to 2π:

| N | Γ at α = 0 | Γ at α = 5° | 4πVa sin α |
|---|---|---|---|
| 40 | −7 × 10⁻¹⁶ | 1.09415 | 1.09523 |
| 80 | −4 × 10⁻¹⁵ | 1.09496 | 1.09523 |
| 160 | 5 × 10⁻¹⁵ | 1.09516 | 1.09523 |

The circulation converges to **4πVa sin α**, with pressure drag of 10⁻¹⁶.

(b) γ₀ + γ_N = 0 at the node θ = 0 forces equal and opposite surface velocities there, which on a smooth circle means zero: a stagnation point at the rear. By `P10.03`, u_θ(0) = −2V sin(−α) − Γ/(2πa) = 0 gives Γ = 4πVa sin α. The method did exactly what it was told.

(c) The Kutta condition is an input, not a result: on a sharp-edged airfoil it is the physically right input, but on a body without a sharp edge the panel method will produce whatever circulation the "trailing-edge" node implies. For a section with a rounded or blunt trailing edge, the inviscid lift depends on where you put that node; the real circulation is set by where the boundary layer separates.

## 4. Your own exact benchmark (s = 15: μ = −0.065 + 0.045i, τ = 12.5°)

Exact at α = 3°: c_l = **0.652351**, c_m,c/4 = **−0.074441**.

(a)

| N | c_l (Γ) | c_l (Cp) | c_m,c/4 | max ΔCp 5–95 % |
|---|---|---|---|---|
| 20 | 7.98 × 10⁻³ | 1.27 × 10⁻² | 6.32 × 10⁻⁴ | 1.11 × 10⁻² |
| 40 | 1.51 × 10⁻³ | 2.30 × 10⁻³ | 6.37 × 10⁻⁴ | 9.58 × 10⁻³ |
| 80 | 6.61 × 10⁻⁴ | 8.96 × 10⁻⁴ | 1.01 × 10⁻⁵ | 3.74 × 10⁻³ |
| 160 | 1.71 × 10⁻⁴ | 2.28 × 10⁻⁴ | 4.08 × 10⁻⁶ | 1.15 × 10⁻³ |
| 320 | 4.34 × 10⁻⁵ | 5.73 × 10⁻⁵ | 1.16 × 10⁻⁶ | 2.91 × 10⁻⁴ |

Orders: c_l(Γ) 2.40, 1.19, 1.95, 1.98; c_l(Cp) 2.46, 1.36, 1.98, 1.99; c_m −0.01, 5.98, 1.30, 1.82; Cp 0.22, 1.36, 1.70, 1.98.

(b) Relative c_l errors: 1.22 %, 0.23 %, 0.10 %, 0.026 %, 0.0067 %. **0.1 % needs N ≈ 80; 0.01 % needs N ≈ 260** (from 0.026 % at 160 with order 2: 160 × √2.6).

(c) At small N the panels do not resolve the leading-edge suction peak and the error is not yet in its asymptotic form (the next term in the error expansion is as large as the first), so successive ratios wander. Once the nose is resolved the second-order term dominates and the order settles near 2. The c_m error, a small difference of large moments, is the last to settle.

**Common mistakes.** Quoting an order from the coarsest pair; comparing Cp at points that are not the same surface point (panel midpoints are off the surface by O(h²), fine, but must be matched by circle angle, not by x).

## 5. Three NACA sections (N = 160)

(a)

| section | slope / 2π | α_L0 | c_m,c/4 at 4° | min Cp at 4° (x/c) |
|---|---|---|---|---|
| 0012 | 1.1001 | 0.000° | −0.0055 | −1.532 (0.012) |
| 2412 | 1.0988 | −2.154° | −0.0612 | −1.439 (0.014) |
| 4412 | 1.0972 | −4.305° | −0.1170 | −1.372 (0.023) |

(b) Thin-airfoil theory: α_L0 = −2.077° (2412), −4.154° (4412); c_m,c/4 = −0.0531, −0.1062. α_L0 differs by 4 %, c_m by 10–15 %: thickness changes the moment more than the zero-lift angle, because it moves the aerodynamic centre while the zero-lift angle stays mostly a camber property.

(c) For the symmetric 0012 at c_l = 0.483, c_m,c/4 = −0.0055 means dc_m/dc_l ≈ −0.0114, so **x_ac ≈ 0.261c** in inviscid flow: thickness moves the aerodynamic centre aft of the quarter chord. (Viscous effects move it back towards 0.25c.)

## 6. Where not to trust the panel Cp

(a) At a trailing edge of angle τ the exact velocity behaves like (distance)^{π/(2π − τ) − 1 + …}; with the Kutta circulation the leading term cancels and the velocity falls to zero as distance to a small power (0.03 for τ = 10°). Over a polygon panel next to the edge the exact Cp changes from nearly 1 to its neighbours' values; no finite panel resolves that, so the largest whole-surface error stays O(1) (0.47 at N = 640). The error is confined to the last panel or two, so integrated quantities still converge.

(b) Stop the boundary layer a little ahead of the edge and evaluate the trailing-edge quantities there; `P10.06` uses x/c = 0.99. The choice changes the predicted drag by a fraction of a per cent (its panel-count study: 0.006146, 0.006171, 0.006182 for 100, 200, 400 panels), far below the model's other errors.

(c) From largest to smallest for c_l of the trainer at cruise: (1) **the finite wing**: an elliptic wing of AR 8.33 has a lift slope 19 % below its section's (2π AR/(AR + 2), `P10.09`, `P10.10`); (2) **the boundary layer**: the measured NACA 0012 slope is 12–17 % below the exact inviscid one (`P10.03`), and one displacement-thickness iteration takes 2.1 % off the 2412's cruise c_l (`P10.06`); (3) **panel count**: 0.014 % at N = 160 for the 2412 (§4 of the lesson); (4) **the trailing-edge corner**: a local Cp error that hardly touches c_l.

## Check-yourself answers

1. N + 1 node values of γ; N tangency equations at the panel midpoints plus the Kutta condition γ₀ + γ_N = 0.
2. That the speeds leaving the trailing edge on the two surfaces are equal (the flow leaves the edge smoothly); without it the system is singular.
3. From the circulation (c_l = 2Γ/(Vc)) and from integrating Cp over the panels; they agree as N grows, and the pressure drag tends to zero.
4. About 2 for c_l, and for Cp away from the trailing edge once the nose is resolved.
5. Because the exact surface velocity near a finite-angle trailing edge varies over a region smaller than any panel.
6. It tends to thin-airfoil theory: for 2 % camber at α = 4° the excess lift is 11.3 % at 12 % thickness and 1.2 % at 1.5 %, falling in proportion to the thickness.
7. Karman–Trefftz: the panel method is run on its exact surface and compared with the exact c_l, c_m and Cp from `P10.03`.
