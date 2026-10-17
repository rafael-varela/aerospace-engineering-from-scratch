# P02.14  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. The seam's normal is at a + 90° from the axis; the normal stress across the seam is σ_y′ of (9) with θ = a, not σ_x′. For the invariants, I₂ in plane stress collapses to σ_xσ_y − τ_xy² and I₃ to zero. The true maximum shear uses σ₃ = 0.
2. For (2), fix i: the only non-zero terms have {j, k} = {m, n} = the other two indices. For Rodrigues, write v = (v · n)n + v_⊥ and note that e′_i · e_j is the j-component of the rotated basis vector. For (5), the step n_q = Q_jq n′_j is the inverse of (4), which is Qᵀ because Q is orthogonal.
3. The 45° gauge reads ε_x′ at θ = 45°: cos 2θ = 0, sin 2θ = 1, so ε_b = (ε_x + ε_y)/2 + ε_xy. Plane stress means σ₃₃ = 0, and the inverse law ε₃₃ = [σ₃₃ − ν(σ₁₁ + σ₂₂)]/E combined with the plane-stress formulas gives ε₃₃ in terms of ε_x + ε_y alone.
4. The seam line that carries the largest normal stress is perpendicular to σ₁, so its angle is θ_p + 90°. Reversing the torque reflects the stress state about the axis: the seam at a sees what the seam at −a saw. σ_vm² for plane stress is a quadratic in τ_xy; solve it for τ_xy and convert to T.
5. The tensor has 16 non-zero entries, not 4: fill every index order of each symmetric pair. The sign of C̄₁₁₁₂ depends on whether you rotate the axes by +θ or −θ; the problem fixes −θ so that the fibres sit at +θ. Engineering shear strain γ_xy = 2ε_xy is why the third column of Q̄ is doubled.
6. Tresca is never less conservative than von Mises; their ratio is between √3/2 and 1. Yield is not what fails a riveted seam: fatigue at the holes does, and that depends on the local stress range, not the average.
