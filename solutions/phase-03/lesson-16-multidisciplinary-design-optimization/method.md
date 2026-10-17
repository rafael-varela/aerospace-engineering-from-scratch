# P03.16  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. L = W₀; D from (1) with q = 2909 Pa; φ = D/L; W_fuel = L(e^{κφ} − 1); W_TO = L + W_fuel; M_root = 3.75 W_TO b/(3π). 2. h = 0.12 S/b; A = M_root/(hσ_a); W_wing = 2700 g (2/3) A b + 60 S. 3. k = (∂M_root/∂W_wing)(∂W_wing/∂M_root) from `aero_mission_partials` and `structures_partials` at the converged point (or your hand formulas of problem 2); sweeps n with kⁿ e₀ ≤ 10⁻⁶ W_wing.

## 2.
(a) M_root = n b W_TO/(3π), W_TO = L e^{κφ}, L = W₀ + W_wing: ∂W_TO/∂L = e^{κφ}(1 + Lκ ∂φ/∂L) with ∂φ/∂L = −q(f_rest + c_dw S)/L² + 1/(qπeb²). ∂W_wing/∂M_root = (2/3)ρ_m g b²/(τSσ_a). Multiply; note the only W₀-dependence is inside L, φ and e^{κφ}. (b) y₂⁺ = G₂(G₁(y₂* + e)) = y₂* + G₂'G₁' e + O(e²). (c) dF/dx = ∂F/∂x − (∂F/∂y)(∂R/∂y)⁻¹(∂R/∂x); define ψᵀ = (∂F/∂y)(∂R/∂y)⁻¹, i.e. (∂R/∂y)ᵀψ = (∂F/∂y)ᵀ. Direct: n_x solves; adjoint: n_F solves; break-even at n_x = n_F. (d) w'' = 2σ_a/(Eh); w(b/2) = σ_a (b/2)²/(Eh) at ultimate; divide by 1.5; divide by b/2; substitute h = τS/b; set equal to δ_max and solve for b²/S.

## 3.
(a) W_b = gRW/(η(L/D)e) = kW; W = W₀ + kW ⇒ W = W₀/(1 − k). k = 1: the battery needed to carry the aircraft weighs as much as the aircraft. (b) Error after n sweeps from W₀: e_n = kⁿ e₀ with e₀ = W − W₀ = kW; solve kⁿ · k/(1 − k) ≤ 10⁻⁶. Newton: one step (linear residual). (c) Total: d/dR [W₀/(1 − k)] = W₀ (dk/dR)/(1 − k)² with dk/dR = k/R. Partial: d/dR [kW] at fixed W = kW/R. Ratio 1/(1 − k). (d) k = 1 ⇒ R = η(L/D)e/g.

## 4.
(a) `mda_gauss_seidel(x0, p, tol=1e-11)`, `mda_newton(x0, p)`, `coupling_strength`; ratios of successive residuals; Newton order from three residuals. (b) `total_derivatives_adjoint`, `total_derivatives_fd(x0, p, 1e-6)`, `functions_partials(x0, y, p)[0][0, 0]`. (c) `mdf_optimize(p, x0)`, `idf_optimize(p, x0)`; constraints from `functions_of_interest(..)[1]` and `design_constraints_x_only`. (d) `dF` from the adjoint at x*, `J = design_constraints_x_only_jac(x*, p)[0]`; `lam = lstsq(J[:, None], dF[0])`; saving = 0.1 λ/g.

## 5.
(a) `fd_convergence_study(x0, p, hs, mda_tol=..., mda_solver=...)`; `np.polyfit(log h, log err, 1)` over the falling part. (b) Write y_GS(x) = y(x) + ε(x) with ε smooth and O(tol); the central difference of F(x, y_GS) is dF/dx + (∂F/∂y) dε/dx + O(h²); dε/dx is O(tol/scale), independent of h. (c) Read the error at h ≈ 10⁻⁷ (the nearest tabulated step) in each column; compare with the gradient accuracy SLSQP needs to see descent (a few per cent of the gradient).

## 6.
(a) Measure k by two or three Gauss–Seidel sweeps; Newton needs coupling derivatives the codes do not give, so use Gauss–Seidel with Aitken or relaxation if k < 1; at k = 0.9 expect ~130 sweeps for 10⁻⁶ and consider a quasi-Newton (Broyden) update of ∂R/∂y from the sweep history. (b) FD: 2 × 150 MDAs per gradient; direct: 150 coupled linear solves; adjoint: 13 (12 constraints + objective). Without discipline partials only FD is possible; the adjoint needs each code's partials (or a differentiable replacement). (c) MDF: every iterate consistent, needs a fast MDA; IDF: disciplines in parallel, no inner loop, but 4 + coupling-variable sets of equality constraints and larger optimiser problem; with expensive codes and no derivatives the MDA count decides. (d) Order by sensitivity of the objective to the assumption: wing weight law near the gate limit, drag polar at the new C_L, mission/fuel-fraction model; each needs an independent estimate (FE run, CFD or wind tunnel, flight data).
