# P08.15  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Re_b uses the full gap as the length; `duct_flow` returns everything else. Dean's τ_w is c_f × ½ρV².
2. (a) Integrate the averaged equation from the wall and use symmetry on the centreline. (b) Multiply the root by (1 + √…)/(1 + √…). (c) Constant stress plus P = ε gives two equations for k and ε once ν_t = C_μk²/ε.
3. ∫ ln y dy = y ln y − y. For (b), the relation is implicit in Re_τ: solve Re_b(Re_τ) = target with a root finder.
4. Change one constant at a time, and fit the log law only where it exists (30 ≤ y⁺ ≤ 0.2Re_τ).
5. In a pipe the pressure force on the end faces of a cylinder of radius r is πr²Δp, and the shear acts on 2πrL. The bulk velocity is the flow rate over πR².
6. Sum a geometric series for the cell count. For (b), compare 30 with 0.1–0.2 Re_τ.
