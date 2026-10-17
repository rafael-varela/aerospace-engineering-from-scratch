# P07.11  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. T from s; ΔG° = g_H2 + ½ g_O2 − g_H2O; K = exp(−ΔG°/RT). 2. Solve α √(α/2)/(1 − α) · (1/(1 + α/2))^½ = K with `brentq` on (0, 1). 3. Full solve; compare x_H2O, x_H2, x_O2; list OH, H, O.

## 2.
(a) L = G/RT − Σ π_i(Σ a_ij n_j − b_i); ∂/∂n_j; the derivative identity. (b) Σ ν_j × (a), atom balance kills the π terms; substitute μ_j. (c) Substitute Δln n_j into the element and total equations; collect π and Δln n. (d) Hessian of Σ n_k ln(n_k/n) is positive semidefinite (it is the ideal mixing term); linear constraints keep the feasible set convex.

## 3.
1. `equilibrium_flame_temperature("Jet-A(g)", phi, "air", 30e5, 298.15, 800)` for both φ; x_NO. 2. `equilibrium_tp` at fixed T for φ = 0.396. 3. K_p from `equilibrium_constant`; x_N₂, x_O₂ from the solve. 4. Kinetics (Zeldovich) and residence time.

## 4.
1. O/F, p_c from s; φ; three solves; `complete_flame_temperature`. 2. Loop over O/F; √(T_c/M̄); explain with light H₂ in the products.

## 5.
1. Two solvers; table; relative differences. 2. `iterations` at 3000 and 1000 K; the 2 × 2 element block when one species holds both elements. 3. `rng.uniform` starting log-amounts; `np.allclose`.

## 6.
1. Write the energy balance with `equilibrium_tp(b0, T, 1e5, species)` and `brentq`. 2. Fewer species = fewer ways to absorb energy = hotter; condensed carbon for rich hydrocarbons. 3. Frozen versus shifting equilibrium; residence time.
