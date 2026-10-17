# P05.09  Method

How to attack each problem, without the final numbers. Read `hints.md` first; read `worked-solutions.md` only after you have an answer of your own.

## 1. Simply supported equipment beam

1. Free-body diagram of the whole beam: the two reactions up, the uniform load as its resultant 1.5L kN at L/2, the point load at L/3.
2. ΣM about x = 0 gives the roller reaction; ΣF gives the pin reaction. Check the sum.
3. On 0 ≤ x < L/3, V = R_A − 1500x; on the right segment subtract a further 4000 N. Evaluate V just left and just right of the point load.
4. If V changes sign by jumping across the point load, M is largest there; if V is still positive just to the right, solve V = 0 in the right segment instead. Which case you are in depends on L (work out the value of L at which V just to the right of the load is zero).
5. M at the maximum by summing moments of the left free body about the cut. Confirm with the code; `extreme_moments` returns ((x, M_max), (x, M_min)).

## 2. Derivations

1. (a) Slice: vertical forces V + w dx − (V + dV) = 0; moments about the right face, dropping dx². For the jumps, integrate across an interval of width ε that contains only the point load and let ε → 0; for the couple, the moment equation picks up −M₀.
2. (b) Write the subtracted load as w_b + k⟨x − b⟩ and show that it equals the original at b in value and slope. Integrate each bracket term with ∫⟨x − a⟩ⁿ dx = ⟨x − a⟩ⁿ⁺¹/(n + 1), with V = M = 0 at x = 0⁻.
3. (c) Reactions by resultant and centroid (the resultant acts at 2L/3). V(x) = R_A − w₀x²/(2L); set it to zero. M(x) = R_A x − w₀x³/(6L); evaluate at the root.
4. (d) Force: ∫₀ˢ √(1 − (x/s)²) dx with x = s sin θ. Moment: x dx = −(s²/2) du with u = 1 − (x/s)². Centre: moment over force.

## 3. Wing on trestles

1. Loads: uniform structure weight 80g over 6.7 m, engine 150g at 1.9 m, both down. Reactions from ΣM about the pin and ΣF.
2. Use `shear_moment(..., left=True)` and the default right-continuous values at each of x₁, 1.9 m and x₂.
3. The overhangs (0 to x₁ and x₂ to the tip) act as cantilevers: hogging moments at the trestles. Between the trestles the engine produces sagging. `extreme_moments` gives both extremes.
4. Sweep x₂ on a fine grid, recompute reactions and extremes for each, and keep the position that minimises max(|M_max|, |M_min|). The optimum balances the sagging and hogging peaks.
5. Compare with the in-flight root moment from the lesson. Think about what else a wing on two points can suffer: local bearing loads at the trestle pads, stability on the trestles, and the position of the centre of gravity relative to them.

## 4. Your own pull-up

1. w₀ from the lesson's formula with your n. Closed-form root shear and moment from `exact_lift_root`. The centre of lift is the ratio; inspect whether n cancels.
2. For each panel count, the lift alone on a fixed root: solve reactions, evaluate V and M at x = 0. Relative error against the closed form; observed order from successive pairs.
3. With relief: `wc.loads(400)` plus reactions. For the independent check, replace the panels by `DistributedFunction(0, s, wc.lift)` plus `wc.inertia_loads()`, solve reactions and call `sections_check(0.0, ...)`.
4. Evaluate V with `left=True` and without at `wc.engine_x`; the difference against n·m_engine·g.

## 5. Burn order

1. Build a function root_M(mass, inboard fuel, outboard fuel): `WingCase(mass=mass, fuel_mass=0.0).loads(400)` gives lift, structure and engine; add one `DistributedLoad` per non-empty tank with intensity −m_fuel·n·g/length; fixed root reactions; M at x = 0.
2. Evaluate the four cases. Explain the zero-fuel case: lift falls by the weight of the burned fuel times n, at the centre of lift, while the relief lost is the same weight at the tanks' centroids. Compare the two lever arms.
3. Evaluate root_M at two masses with empty tanks, fit the straight line and solve for the mass that matches the full-fuel case. Check by evaluating it.
4. The recommendation weighs the root moment against what else outboard fuel does: roll inertia, fuel system complexity and the moment at the outboard end of the inboard tank.

## 6. Judgment

1. (a) The difference between the two-fitting load and the point load is a self-equilibrated set: zero shear and moment outside the fittings; inside, a triangular moment correction whose peak follows from statics of two half loads. Compare with M at the engine from the lesson's table.
2. (b) Compare the 10-panel root moment with w₀s²/3 and check the sign of the error. A coarse model is conservative only if its error has the safe sign.
3. (c) No lift: the root reaction must carry all the inertia loads. Moment about the root of each: structure at its midpoint, fuel at the tank's centroid, engine at its station. Sagging or hogging from the sign; which cap is compressed follows.
4. (d) The lift-alone root moment is lift times centre-of-lift distance; a 5 % change in the distance is a 5 % change in the moment. Compare with the 0.7 m engine move in lesson part 4.
