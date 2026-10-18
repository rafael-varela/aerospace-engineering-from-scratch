# P06.15  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. A, I from the tube formulas; EI = E I; ρA; mass ρAL; μ. 2. `natural_frequencies(EI, rhoA, L, "clamped-free", 3)`. 3. Same with `m_tip=0.2`. 4. `rayleigh_tip_mass_estimate`; percentage against (c).

## 2.
(a) Slice of length dx: V′ = q − ρAẅ, M′ = V, M = EIw″. (b) φ = C₁ cosh + C₂ sinh + C₃ cos + C₄ sin; φ(0) = φ′(0) = 0; φ″(L) = φ‴(L) = 0; determinant; simplify. (c) Newton for the mass: m_t ẅ(L) = shear force from the beam; with w = φ cos ωt; redo (b) with the modified fourth row. (d) Integrate ∫φⱼφᵢ⁗ by parts twice; boundary terms vanish by φ(0) = φ′(0) = 0 and φ″(L) = φ‴(L) = 0; subtract the swapped identity. (e) Apply the stencil to sin(jπxᵢ/L); collect cos terms; check the ghost-node rows.

## 3.
(a) φ″ = φ‴ = 0 at both ends; determinant; rigid translation and rotation have β = 0. (b) `beta_L("clamped-clamped", 3)` and (3). (c) `fem_matrices(L, EI, rhoA, 60)`; `eigh(K, M)`; the first two eigenvalues are round-off around zero. (d) Zeros of the mode-1 shape φ = cosh + cos − σ(sinh + sin), σ = (cosh βL − cos βL)/(sinh βL − sin βL); for rotation, look where φ′ = 0 (the antinodes of deflection).

## 4.
(a) Loop over element counts; `fem_frequencies(..., 5)` against `natural_frequencies(..., 5)`; note the first count at which each mode is below 10⁻⁴. (b) Loop the tank position with `point_masses=[(x, m)]`; arg-min and arg-max of f₂ over the outboard half; zero of `mode_shape` for mode 2. (c) ω² ≈ (φᵀKφ)/(φᵀMφ + m φ(x)²).

## 5.
(a) Half-wavelength of mode n ≈ L/(n − ½); set it equal to ten depths. (b) Stiffness (skins), mass distribution, root fixity; directions. (c) Read experiment 1's lumped column at 8 elements. (d) Separate scaling laws and modelling choices from absolute numbers for the real wing.
