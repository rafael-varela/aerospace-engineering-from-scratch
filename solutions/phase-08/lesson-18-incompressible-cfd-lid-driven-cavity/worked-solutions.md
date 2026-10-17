# P08.18  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Everything is dimensionless (lengths in L, velocities in U, ψ in UL, ω in U/L). Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe scripts that call it, given to the figures shown. Timings are from one laptop core and will differ on yours.

## 1. The wall and the cell

(a) N = 33: h = 1/32 = **0.03125**, U h Re = **3.125**, lid term −2U/h = **−64**. Cell Reynolds numbers: **6.25** (17), 3.125 (33), **1.5625** (65), **0.78125** (129). The 17 and 33 grids are outside U h Re ≲ 2 at Re = 100; they still converge (the steady problem is mildly convective), but the coarse-grid profiles are visibly off (`cavity.png`).

(b) On y = 0, ψ = 0 for all x so ψ_xx = 0, and ω = −ψ_yy = −2x²(1 − x)²; at x = 0.5, **ω = −0.125 exactly**. Thom's estimate: **−0.10986** (17), **−0.11731** (33), **−0.12112** (65); errors 0.01514, 0.00769, 0.00388; observed orders **0.98, 0.99**: first order, as the Taylor series predicts.

(c) The O(h) error is confined to the wall nodes, one cell thick; it enters the interior equations only through the stencils of the first interior row, so its effect on the interior solution is a perturbation of size O(h) acting over a layer of width O(h), which is O(h²). The tests and the convergence study confirm order 2 in the interior.

**Check.** The lid term in (a) is exactly what `wall_vorticity(zeros, U=1)` returns on the lid (a test checks it to 10⁻¹⁴).

**Common mistakes.** Using N instead of N − 1 for 1/h; evaluating Thom at a corner node.

## 2. Derivations

(a) Write the momentum equations as u_t-free forms: u u_x + v u_y = −p_x + ν∇²u and u v_x + v v_y = −p_y + ν∇²v (ν = 1/Re). Apply ∂/∂x to the second and ∂/∂y to the first and subtract. Pressure: −p_yx + p_xy = 0, because mixed partial derivatives commute; that is why it disappears. Viscous terms: ν∇²(v_x − u_y) = ν∇²ω. Convective terms: ∂_x(u v_x + v v_y) − ∂_y(u u_x + v u_y) = u ω_x + v ω_y + (u_x v_x + v_x v_y − u_y u_x − v_y u_y). The bracket is v_x(u_x + v_y) − u_y(u_x + v_y) = ω(u_x + v_y) = 0 by continuity. Result: **u ω_x + v ω_y = ν∇²ω**.

(b) u_x + v_y = ψ_yx − ψ_xy = 0. ω = v_x − u_y = −ψ_xx − ψ_yy, so **∇²ψ = −ω**.

(c) On the lid y = 1, expand one cell into the fluid: ψ(x, 1 − h) = ψ − hψ_y + (h²/2)ψ_yy − (h³/6)ψ_yyy + …, all evaluated on the wall. There ψ = 0, ψ_y = u = U, and since ψ is constant along the wall ψ_xx = 0, so ψ_yy = −ω_w. Hence ψ₁ = −hU − (h²/2)ω_w + O(h³), and **ω_w = −2ψ₁/h² − 2U/h + O(h)**: dividing the O(h³) remainder by h² leaves O(h). On a fixed wall the same steps with zero wall speed give ω_w = −2ψ₁/h² (on the bottom wall the expansion is in +y and ψ_y = 0, with the same result).

(d) On the full grid D_x = D ⊗ I and D_y = I ⊗ D, with D the 1-D centred difference. Then D_x D_y = D ⊗ D = D_y D_x, so D_x(D_yψ) + D_y(−D_xψ) = 0 at every interior node. The wall velocities set by `velocities` (u = U on the lid, zero elsewhere, v = 0) enter the divergence only through u at the side walls and v at the bottom and lid, where they equal the differenced values because ψ = 0 along those walls. Flux: with u_0 = 0, u_{N−1} = U and u_j = (ψ_{j+1} − ψ_{j−1})/(2h), the trapezoid sum h[U/2 + Σ_{j=1}^{N−2}(ψ_{j+1} − ψ_{j−1})/(2h)] telescopes to **hU/2 + (ψ_{N−1} + ψ_{N−2} − ψ_1 − ψ_0)/2 = hU/2 + (ψ_{N−2} − ψ_1)/2**. Near the lid ψ_{N−2} = −hU − (h²/2)ω_lid + …, near the bottom ψ_1 = −(h²/2)ω_b + …, so the sum is (h²/4)(ω_b − ω_lid) + O(h³) = **O(h²)**. The code gives 0.0088, 0.0019, 0.00045 on 17, 33, 65 (order 2.2, 2.1).

**Common mistakes.** Forgetting that continuity is what kills the ω(u_x + v_y) term (in compressible flow it stays); a sign error in the lid direction in (c).

## 3. Grid convergence at Re = 65 (s = 15)

(a)

| N | ψ_min | x_c | y_c | u(0.5, 0.5) |
|---|---|---|---|---|
| 17 | −0.095400 | 0.6057 | 0.7609 | −0.18604 |
| 33 | −0.100129 | 0.5952 | 0.7563 | −0.20096 |
| 65 | −0.101268 | 0.5933 | 0.7550 | −0.20463 |
| 129 | −0.101540 | 0.5929 | 0.7546 | −0.20552 |

(b) On 33/65/129: u(0.5, 0.5) order **2.05**, extrapolate **−0.20580**, GCI **0.17 %**; ψ_min order **2.07**, extrapolate **−0.101625**, GCI **0.11 %**. On 17/33/65: orders **2.02** (u) and **2.05** (ψ_min).

(c) The centre moves **5.7 × 10⁻⁴** (less than a tenth of the 129 grid's spacing) between 65 and 129. The 65 grid's ψ_min is 0.35 % from the extrapolate: **yes**, adequate for 1 %.

**Check.** The Re = 100 study of the lesson gives the same orders (2.04–2.07) and a slightly larger GCI; a lower Re is a smoother flow.

**Common mistakes.** Taking the order from ψ_min located at the minimum node (not refined): the node position jumps between grids and the "order" is noise.

## 4. A faster stream over the cavity

(a) From the Stokes start at Re = 1000 on 65 × 65 the residuals go 1049, 360, 623, 6040, … and the call raises "Newton did not converge". Newton's method converges quadratically only inside a basin round the solution; the Stokes field (a symmetric vortex) is far from the Re = 1000 flow, whose vortex has moved and whose wall layers are thin, so the linearisation is poor and the steps overshoot.

(b) 65 × 65 by continuation:

| Re | iterations | ψ_min | centre |
|---|---|---|---|
| 100 (Stokes start) | 4 | −0.102768 | (0.617, 0.739) |
| 400 | 5 | −0.109052 | (0.558, 0.610) |
| 1000 | 6 | −0.105252 | (0.535, 0.571) |

(c) 129 × 129: Re = 400 in 5 iterations, ψ_min −0.112776 at (0.555, 0.607); Re = 1000 in 7, **−0.115473** at (0.532, 0.567). 33 × 33: Re = 400 converges (−0.094132 at (0.574, 0.631)); Re = 1000 **does not converge**: the cell Reynolds number is U h Re = 1000/32 = 31, far outside U h Re ≲ 2, and the central-difference steady problem has no solution near the start (or no physical one).

(d) −0.10525 (65) against −0.11547 (129): a **9 %** change in one refinement, so neither is converged; a 257 grid (or a better convection scheme) is needed to say anything to a few per cent. The Re = 100 study (GCI 0.2–0.4 % on 129) says nothing about Re = 1000: the grid needed grows with Re because the wall layers thin like Re^(−1/2) and the cell Reynolds number grows like Re. A grid study belongs to one flow, not to a code.

**Common mistakes.** Taking a smaller continuation step as a cure for a grid that is too coarse.

## 5. What each route costs (Re = 100)

(a)

| N | Newton its | Newton time | march steps | dt | march time |
|---|---|---|---|---|---|
| 17 | 4 | 0.023 s | 2463 | 0.0160 | 0.22 s |
| 33 | 4 | 0.05 s | 2461 | 0.0160 | 0.38 s |
| 65 | 4 | 0.25 s | 8041 | 0.00488 | 5.3 s |

(b) The limits are h²Re/4 (diffusion) and 2/(Re U²) = 0.02 (convection–diffusion, central differences). h²Re/4 = 0.098 (17), 0.024 (33), 0.0061 (65). On 17 and 33 the smaller is 0.02, so dt = 0.8 × 0.02 = 0.016 on both and the step count is the same; on 65 diffusion takes over, dt = 0.8 × 0.0061 = 0.0049, and the steps grow by 0.016/0.0049 = 3.3. From there on dt ∝ h².

(c) On 129, dt = 0.8 × (1/128)² × 100/4 = 0.0012, a quarter of the 65 value, so about **32 000 steps**, assuming the flow needs the same physical time (about 39 time units) to settle. Each Poisson solve has four times the unknowns and costs about four to five times more with a sparse LU, so about **100 s**, fifty times Newton's 2 s.

(d) On 33 × 33, max |ψ_march − ψ_Newton| = **5.0 × 10⁻¹¹**. The march stops when max|∂ω/∂t| < 10⁻⁹, which bounds its distance from the steady state by roughly that rate times the slowest decay time (a few time units), about 10⁻⁹–10⁻⁸ in ω and less in ψ. The agreement is 10⁶ times below the discretisation error (GCI ~10⁻³ in ψ), so the two routes solve the same discrete equations.

**Common mistakes.** Comparing the two routes' answers with each other and calling that a validation: it verifies the algebra (the Jacobian), not the physics or the grid.

## 6. Corners, published data and what to claim

(a) Next to the top-right corner, ω[N − 2, N − 1] = **−17.4, −39.8, −85.9, −178.6** on 17, 33, 65, 129: it doubles with each refinement, ω ∝ 1/h. Mid-lid ω = **−8.46, −6.99, −6.66, −6.59**: it converges. The lid velocity jumps from U to 0 at the corner, so over one cell the shear is about U/h; the exact solution has infinite vorticity at the corner, and the grid just samples it closer each time. Quantities away from the corner converge at order 2 regardless; the singularity is local.

(b) Against. The code's own error on the 129 grid is about 0.3 % of the lid speed (GCI 0.28 % for the centre u), it converges at order 2.07, and its Richardson-extrapolated profiles differ from the table by up to 0.005 (u) and 0.009 (v): the same size as the 129 − Ghia difference. A code whose error falls as h² cannot be getting worse; what is not falling is the difference from a table that has its own discretisation error, of order half a per cent. The 65 grid's smaller difference is a coincidence of two errors of opposite sign.

(c) Can: show that the code solves the 2-D steady laminar equations correctly and how fine a grid is needed at low Re; give the structure of a lid-driven recirculation. Cannot: a wheel well at Re = 10⁶ is turbulent, three-dimensional and unsteady; the "lid" is a shear layer from an incoming boundary layer, not a wall; the opening and the gear inside change the geometry; and central differences on a uniform grid would need U h Re ≲ 2, about 10⁶ cells in each direction. Recommendation: use the code as a verified building block; for the well, a RANS or scale-resolving solver with a turbulence model (`P08.15`), validated against a cavity-flow experiment at matching Re and Mach number.

## Check-yourself answers

1. Because u_x + v_y = ψ_yx − ψ_xy = 0 for any smooth ψ. No fluid crosses a wall, so the wall is a streamline, and ψ is constant along a streamline (taken as zero).
2. The pressure gradient: the curl of a gradient is zero. It removes the pressure, the unknown with no evolution equation of its own, and with it the need for a pressure–velocity coupling.
3. ω_lid = −2ψ₁/h² − 2U/h, from a Taylor series one cell into the fluid; first-order accurate at the wall.
4. The two centred differences are Kronecker products acting on different indices, so they commute: D_x D_yψ − D_y D_xψ = 0 exactly, up to round-off.
5. p = ln[(0.20719 − 0.20099)/(0.20867 − 0.20719)]/ln 2 = ln 4.19/ln 2 = 2.07. The grids are in the asymptotic range and the scheme achieves its design order, so Richardson extrapolation is justified.
6. U h/ν (U h Re in dimensionless form): the ratio of convection to diffusion across one cell. When it exceeds about 2, central differences of the convection term produce oscillations and the steady discrete problem can lose its solution.
7. Because the benchmark has its own error (here a 129 grid solution), and once your error falls below it the difference levels off. Tell them apart by your own convergence study: if your solution converges at the design order to an extrapolate that differs from the benchmark by the plateau, the plateau is the benchmark's.
8. Good: a simple geometry with all the hard features (recirculation, wall layers, corner singularities); one parameter (Re); published high-quality solutions. Poor: real wells are 3-D, turbulent, unsteady, and driven by a free shear layer rather than a sliding wall.
