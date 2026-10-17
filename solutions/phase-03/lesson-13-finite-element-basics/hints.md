# P03.13  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Each element contributes (EA/h)[[1, −1], [−1, 1]] to the 2 × 2 block of its two nodes; the middle node receives both. Strike the first row and column for the constraint. The consistent load of a uniform q over an element is qh/2 at each node, and a uniform load is a case where the uniform bar's nodal exactness applies.
2. Multiply by v, integrate over the bar, integrate the derivative term by parts; the boundary term at L is the tip force, the one at 0 vanishes because v(0) = 0. In the parent coordinate dx = (h/2)dξ and N' = ∓1/h. For (c), write vᵀKv as an integral of a square. For (d), differentiate Π with respect to u.
3. σ = E(u' − αΔT) inside ∫σv' dx gives ∫EAu'v' dx − ∫EAαΔT v' dx: the second integral is a load. For a hat function v' = ∓1/h on the element, so the element load is EAαΔT times [−1, +1]. With both ends fixed every node stays put.
4. P = m_tipΩ²L; use `solve_rotor_blade` and `exact_rotating_bar`; the element stress is the exact element average. Linear extrapolation from element centres is exact for a stress that varies linearly and over-shoots for the quadratic one here.
5. Lumped: `np.diag(M.sum(axis=1))`. Fit p from two consecutive errors with `observed_order`. Mode k has k − ½ half-wavelengths on the bar; count elements per half-wave at the n where the ratio is within 1 %.
6. EA(x_mid)u' = P in the last element gives σ_e = P/A(x_mid) exactly; the tip is h/2 further out. The natural boundary condition says what the tip stress is without any solve.
