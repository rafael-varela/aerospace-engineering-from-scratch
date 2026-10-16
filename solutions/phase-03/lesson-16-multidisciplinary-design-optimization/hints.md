# P03.16  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. With W_wing = 0 the lift is W₀ exactly; q = ½ρV² = 2909 Pa and κ = 2.574 are fixed numbers. For (b), after the first sweep the error in W_wing shrinks by the factor k every sweep, so the count is ln(target/initial error)/ln k.
2. (a) Write W_TO = L e^{κφ} with φ = D/L and differentiate L e^{κφ} with respect to L; (4) is linear in M_root. (b) Substitute y₂ = y₂* + e into G₂(G₁(y₂)) and expand to first order. (c) Insert dy/dx = −(∂R/∂y)⁻¹ ∂R/∂x into (7) and move the inverse onto ∂F/∂y by transposition. (d) Curvature κ_b = M/(EI) = σ_a/(E h/2) when the outer-fibre stress is σ_a everywhere; constant curvature means a parabola.
3. The battery weight is proportional to the total weight: W_b = k W with k = gR/(η (L/D) e), so the fixed point is linear and solvable in one line. Newton on a linear residual is exact in one step. For (c), differentiate W = W₀/(1 − k) with respect to R, remembering that k ∝ R.
4. Everything is a call: `mda_gauss_seidel`, `mda_newton`, `coupling_strength`, `total_derivatives_adjoint`, `total_derivatives_fd`, `functions_partials`, `mdf_optimize`, `idf_optimize`. For (d), at the optimum ∇F = λ ∇g_defl with the deflection constraint the only active one; solve for λ by least squares with `np.linalg.lstsq` on the one active gradient row.
5. Central differences cancel errors that are *odd* in h and errors that are the *same* on both sides. The Gauss–Seidel stopping error is nearly the same at x + h and x − h, so what survives is its derivative, not its value over h.
6. Think in counts: an FD gradient costs 2 n_x MDAs, the direct method n_x linear solves, the adjoint n_F; and a commercial code with no derivatives can only be differentiated by finite differences of its outputs, which is a very different accuracy from a coupled adjoint.
