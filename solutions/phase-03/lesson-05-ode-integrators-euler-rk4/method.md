# P03.05  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. f(t, y) = −y + t. Euler: y₁ = 1 + h(−1). Heun: k₁ = −1, k₂ = f(h, 1 + hk₁); y₁ = 1 + h(k₁ + k₂)/2. RK4: k₁ = −1, k₂ = f(h/2, 1 + hk₁/2), k₃ = f(h/2, 1 + hk₂/2), k₄ = f(h, 1 + hk₃); y₁ = 1 + h(k₁ + 2k₂ + 2k₃ + k₄)/6. Compare with y(h) = h − 1 + 2e^{−h}. 2. `global_error` at h, h/2, h/4; ratios e(h)/e(h/2) and e(h/2)/e(h/4); `observed_order`. 3. The error is Chᵖ + Dh^{p+1} + …; the ratio is 2ᵖ(1 + (D/C)h)/(1 + (D/C)h/2) and tends to 2ᵖ only as h → 0.

## 2.
(a) Local: y(t + h) − [y + hf(t, y)] = h²y″(ξ)/2 by Taylor with remainder. Global: let eₙ = |y(tₙ) − yₙ|; yₙ₊₁ − y(tₙ₊₁) = [yₙ − y(tₙ)] + h[f(tₙ, yₙ) − f(tₙ, y(tₙ))] − h²y″/2, so eₙ₊₁ ≤ (1 + hL)eₙ + h²M/2; iterate from e₀ = 0: eₙ ≤ (h²M/2)Σ(1 + hL)ᵏ = (hM/2L)[(1 + hL)ⁿ − 1] ≤ (hM/2L)(e^{LT} − 1). (b) Heun: k₁ = λy, k₂ = λ(y + hλy) = λy(1 + z); y₁ = y + (h/2)(k₁ + k₂) = y(1 + z + z²/2). RK4: k₁ = λy, k₂ = λy(1 + z/2), k₃ = λy(1 + z/2 + z²/4), k₄ = λy(1 + z + z²/2 + z³/4); combine with weights 1, 2, 2, 1 over 6. Heun's limit: 1 + z + z²/2 = 1 at z = −2 and is > −1 everywhere on the negative axis. RK4: solve |R(−x)| = 1 numerically (`stability_limit_real(4)`). Euler on iω: |1 + ihω| = √(1 + h²ω²) > 1. (c) f = g(t): k₁ = g(t), k₂ = k₃ = g(t + h/2), k₄ = g(t + h); the step adds (h/6)[g(t) + 4g(t + h/2) + g(t + h)]. (d) n = T/h steps each contributing ≤ Ch^{p+1}, propagated with a factor ≤ e^{L(T − tₙ)}: total ≤ e^{LT}CTh^p. For y′ = −ky (k > 0) the propagation factor is e^{−k(T−tₙ)} < 1 and early errors are damped.

## 3.
(a) Formulas as given; the no-drag apogee is v₀²/2g. (b) `integrate(rhs, (0, 14), [0, v₀], h, method)`; z_exact(t) = (v_t²/g) ln[cos(gt/v_t − arctan(v₀/v_t))/cos(arctan(v₀/v_t))]; `observed_order`. (c) Integrate to 16 s; first index i with v < 0; t_a ≈ tᵢ₋₁ − vᵢ₋₁Δt/(vᵢ − vᵢ₋₁); z_a ≈ zᵢ₋₁ + vᵢ₋₁δ + ½ā δ² with δ = t_a − tᵢ₋₁ and ā = (vᵢ − vᵢ₋₁)/Δt. Compare ratios of successive errors with 2, 4, 16. (d) The required h for 0.1 m from the RK4 errors and (2); evaluations = 4 × 15/h.

## 4.
(a) `orbit_period(a)`; perigee altitude a(1 − e) − R_E; speeds |v| of `kepler_state` at t = 0 and t = T/2 (check against √(μ(1 ± e)/(a(1 ∓ e)))). (b) `orbit_position_error` for the three methods and three step counts; `observed_order(T/n, errs)`. (c) `step_for_tolerance(T/800, e_800, 4, 1.0)`, n = ⌈T/h⌉, then a run at n. (d) Errors at 1, 2, 4 orbits; the growth law is close to linear; the one-day error at 800 steps is e₁ × (86 400/T) × (growth factor); then (2) with order 4; check with `orbit_position_error(a, e, 86400/T, n, "rk4")`; evaluations 4 × n × 86 400/T. (e) `energy_drift(propagate_orbit(a, e, 5, 400, m)[1])[-1]`.

## 5.
(a) `max_stable_step(p, -lam)` for p = 1, 2, 4; runs with 0.95 and 1.05 of it to t = 10; |y(10)| against cos 10. (b) steps = ⌈10/(0.95 h_max)⌉, times stages; RK4 error at that step against the exact formula; two runs at 0.01 and 0.005 give the order-4 constant, then (2) for 10⁻⁶. (c) `solve_ivp(f, (0, 10), [0.0], method=..., rtol=1e-6, atol=1e-9).nfev`. (d) Euler on y′ = iy: |R|² = 1 + h², n = 2π/h steps per period; RK4: `abs(stability_function(1j*h, 4))` raised to the n.

## 6.
(a) Error after 13.4 orbits at 800 steps ≈ 13.4 × 4.7 mm × (settling factor measured in Problem 4 or the lesson's part 4); (2) with p = 4 for 10 m; add margin because the one-day order is slightly above 4 over this range; evaluations = 4 × steps × 13.4 × 2000; time = evaluations / 2 × 10⁶. (b) Same with p = 2 (Heun) and p = 1 (Euler) from 4.8 km and 3.8 × 10⁶ m. (c) Fixed step: the perigee speed ratio (1 + e)/(1 − e) and the fourth-power error law. (d) rtol 10⁻³ on a 7.5 × 10⁶ m position is kilometres; compare with the Kepler closed form or with DOP853 at rtol 10⁻¹³.

## 7.
(a) y″ = d(f(y))/dt = f′y′ = f′f; y‴ = (f′f)′ = f″f·f + f′·f′f. (b) Taylor of f about y with increment c₂hf. (c) h¹: b₁ + b₂ = 1; h²: b₂c₂ = ½; h³: the method gives b₂c₂²f″f²/2 and no f′²f term, the exact solution needs f″f²/6 + f′²f/6; the second cannot be matched. (d) Induction over stages: kᵢ = λy Pᵢ(z) with Pᵢ a polynomial of degree ≤ i − 1; the step y Σbᵢz Pᵢ(z) has degree ≤ s. Order p means the step agrees with the exact e^z y through zᵖ; a polynomial of degree ≤ s = p that agrees with e^z through zᵖ is the truncated series.
