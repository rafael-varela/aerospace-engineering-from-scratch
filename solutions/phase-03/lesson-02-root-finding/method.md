# P03.02  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. M in radians; E₀ = M + e sin M. 2. For k = 0, 1: f = E_k − e sin E_k − M, f′ = 1 − e cos E_k, E_{k+1} = E_k − f/f′. 3. Residual f(E₂). 4. Error estimate: e₂ ≈ (f″/2f′) e₁² with e₁ ≈ |E₂ − E₁| and f″ = e sin E. 5. `solve_kepler(M, e)`: steps = `iterates.size − 1`; ν from `true_anomaly`. 6. `bisection_steps_needed(M − e, M + e, 1e-14)`.

## 2.
(a) 0 = f(x_k) + (x* − x_k) f′(x_k) + ½(x* − x_k)² f″(ξ); divide by f′(x_k), rearrange to x* − x_{k+1}; assume f′(x*) ≠ 0 and f′ continuous so f′(x_k) stays away from zero. (b) Substitute f = (x − x*)²h into f/f′ and expand in (x − x*); for h ≡ 1 the step is exactly (x − x*)/2. (c) x_{k+1} − x* = g(x_k) − g(x*) = g′(ξ)(x_k − x*); converges if |g′| ≤ ρ < 1 near x*. (d) The linear interpolant p through the two points satisfies f(x) − p(x) = ½f″(ξ)(x − x_{k−1})(x − x_k); evaluate at x* where f = 0, and expand p(x*) about the secant root x_{k+1} as p′·(x* − x_{k+1}) with p′ ≈ f′; then substitute the power-law ansatz and match exponents.

## 3.
(a) Differentiate; use sinh H ≥ H to get a point where f > 0; f(0) < 0 gives the bracket. (b) `newton(f, df, H0)` with f = e sinh H − H − M_h, df = e cosh H − 1; compare `iterates.size` for the two starts; `brentq` on the bracket. (c) ν = 2 atan(√((e + 1)/(e − 1)) tanh(H/2)); ν_∞ = arccos(−1/e).

## 4.
(a) `solve_kepler(M, e, method)` for the four methods; steps and `n_evals`. (b) Reference `brentq(..., xtol=1e-15)`; `errors_from_iterates`; `fitted_convergence_order` (catch the "too few errors" case for Newton); `observed_linear_rate` for fixed point against e|cos E*|. (c) Digits gained per evaluation = log₁₀(e₀/e_final)/evaluations; re-weight Newton's evaluations by (1 + cost ratio)/2.

## 5.
(a) g(t) = ‖`orbit_state(a, e, t)[0]`‖ − R; `bisection(g, 0, T/2, tol=1e-6)` and `secant(g, 1000, 2000, tol=1e-6)`; solve counts are `n_evals`. (b) E = arccos((1 − R/a)/e); M = E − e sin E; t = M/n. (c) By symmetry the second crossing is at T − t; on [0, T], g(0) and g(T) have the same sign (both at perigee), so no sign change.

## 6.
(a) Evaluations per step for each method against the budget; Newton's start E₀ = M + e sin M or π; reasons to reject bisection (count), fixed point (rate near 1 at high e), secant (needs two starts and can diverge). (b) Step-size test at a few ulp of E; cap near 10–12; on cap, bracket [M − e, M + e] and bisect (or Brent). (c) f′ → 1 − e near perigee, overshoot from the standard start; switch the start to π or use a cubic initial guess, count the extra evaluations. (d) |f(E)| = |f′(ξ)| × error with f′ ∈ [1 − e, 1 + e]; convert a residual tolerance to an error tolerance at e = 0.95.

## 7.
(a) g(E₁) − g(E₂) = e(sin E₁ − sin E₂) = e cos ξ (E₁ − E₂). (b) Uniqueness by contradiction with the contraction inequality; existence: |E_{k+1} − E_k| ≤ eᵏ|E₁ − E₀| so the sequence is Cauchy; the limit is a fixed point by continuity; the error bound by induction. (c) k ≥ log(tol/|E₀ − E*|)/log e; compare with `solve_kepler(M, e, "fixed_point", tol=1e-12)`; the true rate is e|cos E*| ≤ e.
