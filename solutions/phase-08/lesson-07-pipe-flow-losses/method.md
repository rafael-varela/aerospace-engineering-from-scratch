# P08.07  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. `reynolds(1.204, V, 0.10, 1.81e-5)`; compare with 2300 and 4000. 2. `friction_colebrook(Re, 0.15e-3/0.10)` (returns the iterates: count them) and `friction_haaland`. 3. `darcy_weisbach_dp(f, 8, 0.1, rho, V)`; Q = VπD²/4; power = Q dp; `minor_loss_dp(K_total, rho, V)` with K_total = 4 × 0.3 + 0.5 + 1.0. 4. `friction_colebrook(Re, 0)` and (2 log₁₀(3.7/(ε/D)))⁻².

## 2.
(a) Replace f, V and Re; simplify to 128μLQ/(πD⁴). (b) g′ = 1 + 2b/((a + bx) ln 10) > 0, g″ = −2b²/((a + bx)² ln 10) < 0; Newton x_{k+1} = x_k − g/g′; argue with the tangent of a concave function lying above it. Fixed point: φ(x) = −2 log₁₀(a + bx), rate |φ′(x*)| = 2b/((a + bx*) ln 10). (c) dp = f(L/D)ρV²/2 with V = 4ṁ/(ρπD²); f constant → D⁻⁵; f ∝ Re^{−1/4} ∝ D^{1/4} → D^{−4.75}. Check: `segment_loss` at D = 10, 20, 40 mm with ε = 0 and with ε/D fixed, slope of log dp against log D.

## 3.
(a) `segment_loss(3.68, RHO_LOX, MU_LOX, 0.025, 1.5, 1.5e-6, K)`; compare Re with the fuel line's and f with the smooth value. (b) `injector_area_for_dp(3.68, RHO_LOX, 0.7, 0.4e6)`, then `feed_line_budget` with head 0.6 m, 3 g. (c) Fix each tank pressure; for p_c = 2.1 MPa solve `feed_line_budget(m, ...)["p_tank"] = p_tank` for m with `brentq`, once for each propellant; form O/F.

## 4.
1. `line_diameter_for_loss(mdot, RHO_RP1, MU_RP1, 2.5, ROUGHNESS_DRAWN, K, 0.1 * 2e6)` and with 0.05. 2. ID = OD − 2 × wall; `segment_loss` for each ID; loss = friction + minor; express as a percentage of 2 MPa. 3. Think about mass, the propellant trapped in the line, and the line's response time.

## 5.
(a) Loop x ← −2 log₁₀(a + bx) until |x − x*| < 10⁻¹² x*, with x* from Newton; `friction_colebrook(Re, e, x0=7.0)` and the default start; errors |x_k − x*|; ratios e_{k+1}/e_k². (b) `segment_loss(mdot, 1004, 0.9e-3, 0.006, 1.5, 1.5e-6, 0)`; 128μLQ/(πD⁴) with Q = ṁ/ρ; repeat for the larger flows and look at Re.

## 6.
Four pairs of `feed_line_budget` runs: μ × 0.7 and × 1.3; roughness × 10 and × 30; every K × 0.5 and × 1.5; Cd = 0.665 and 0.735 with the injector area kept at its design value. Subtract the baseline tank pressure; rank by the largest change.
