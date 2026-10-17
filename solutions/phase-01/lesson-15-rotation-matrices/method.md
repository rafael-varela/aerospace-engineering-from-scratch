# P01.15  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. Convert the three angles to radians. 2. Write R₃(ψ), R₂(θ), R₁(φ) numerically. 3. Form R₂R₃, then R₁(R₂R₃). 4. For (b), sum the squares of row 1; dot row 1 with row 3; expand the determinant along row 1 with the three 2×2 minors. 5. For (c), θ = −arcsin C₁₃, ψ = atan2(C₁₂, C₁₁), φ = atan2(C₂₃, C₃₃); compare with `euler_321_from_dcm(C)`. 6. For (d), g times the third column; `axis_angle_from_dcm(C)`.

## 2.
(a) Transpose of a product reverses order; CᵀC = I for each factor; determinant of a product is the product of determinants. For the transpose, (Cᵀ)ᵀCᵀ = CCᵀ = I and det Cᵀ = det C. (b) Write the dot product as uᵀv and insert CᵀC = I. Take u = v for the length. (c) Substitute sin θ = ±1, cos θ = 0 into rows 2 and 3 of the full matrix; group the terms as sin φ cos ψ ∓ cos φ sin ψ and sin φ sin ψ ± cos φ cos ψ; recognise sin(φ ∓ ψ) and cos(φ ∓ ψ). State that only the combination φ ∓ ψ is determined. (d) C₁₁ and C₁₂ share the positive factor cos θ, so atan2 sees (sin ψ, cos ψ) up to a positive scale and returns ψ in the right quadrant; same for the pair (C₂₃, C₃₃) with cos θ > 0.

## 3.
(a) `dcm_from_axis_angle([1, 1, 1], 120 deg)`; the matrix is a permutation; read each row as "tracker axis *i* in body coordinates". (b) C_b/i = C_b/t C_t/i = C_t/bᵀ C_t/i. The other order would apply the mounting rotation in the inertial frame, where it does not act. (c) Thrust in inertial = C_b/iᵀ (1, 0, 0) = first column of C_b/i read as a column, i.e. first row of C_b/i. (d) `euler_321_from_dcm(C_b/i)`; compare with (40°, 20°, 10°) and the mounting angle; conclude.

## 4.
(a) Loop over N, call `integration_errors` twice (with and without `renormalize`); fit `np.polyfit(np.log(dt), np.log(err), 1)`. (b) Compute the prediction for each N and the ratio. (c) Build `step = dcm_from_axis_angle(axis, rate * dt)`, multiply it N times, measure `orthogonality_error` and `misalignment_angle` against `dcm_from_axis_angle(axis, rate * T)`. Explain that for a fixed axis the exact steps commute and compose to the exact total. (d) Solve ωT(ωΔt)²/3 = 0.1° in radians for Δt; f = 1/Δt.

## 5.
(a) Write **v** = **n**(**n**·**v**) + **v**⊥ with **v**⊥ = **v** − **n**(**n**·**v**); in the plane perpendicular to **n** use the orthogonal pair (**v**⊥, **n** × **v**); the components of the fixed vector in a frame turned by *a* are cos *a* **v**⊥ − sin *a* (**n** × **v**) for the passive sense (plus for active); assemble, replace **n**(**n**·**v**) by **nn**ᵀ**v** and **n** × **v** by [**n**×]**v**; test with **n** = e₃. (b) tr(**nn**ᵀ) = 1, tr[**n**×] = 0; [**n**×] is antisymmetric, **nn**ᵀ and I are symmetric. (c) d(arccos x)/dx = −1/√(1 − x²), infinite at x = 1; set x = 1 − ε and use arccos(1 − ε) ≈ √(2ε). (d) Substitute cos *a* = −1, sin *a* = 0; the result is even in **n**.

## 6.
(a) δ = 2π × 10 / 100; orthogonality error δ²; actual angle atan δ in degrees. Decide from whether δ is small. (b) δ_max = √(3 × 0.1° in rad / (ωT)); f = ω/δ_max; comment; propose exact Rodrigues steps or quaternion integration. (c) List the functions that take or return angles versus those that work on the matrix. (d) Earth rate × 2 h, or wing flex, with a number.
