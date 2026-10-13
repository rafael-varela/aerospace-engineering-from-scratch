# P01.12  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15. Every number below is reproduced by `solution.py` at the stated resolution.

## 1. Divergence, curl and flux by hand

**Derivation.** (a) ∇ · **v** = s + (−1) + (1 − s) = 0 for every s. ∇ × **v**: each component of **v** depends only on its own coordinate, so every mixed partial is zero: (0, 0, 0). ∇ · **F** = ∂(y²)/∂x + ∂(z²)/∂y + ∂(x²)/∂z = 0. ∇ × **F** = (∂(x²)/∂y − ∂(z²)/∂z, ∂(y²)/∂z − ∂(x²)/∂x, ∂(z²)/∂x − ∂(y²)/∂y) = (−2z, −2x, −2y). SymPy agrees.

(b) Both fields are divergence-free, so by the divergence theorem their flux through *any* closed surface is zero: 0 and 0 m³/s. `box_flux` at n = 3 returns 0.0 for both (for **v**, which is linear, the midpoint rule is exact; for **F**, each component is constant along its own coordinate, so the two opposite faces cancel exactly).

(c) **G** = (15x², 0, 0): flux through x = 1 is 15 × 1 × (area 1) = 15; through x = 0 it is 0; the y and z faces carry nothing. Volume integral: ∫₀¹ 30x dx = 15. `divergence_theorem_check` at n = 2 gives (15.0, 15.0, 0.0): the integrand 30x is linear, the midpoint rule is exact for it, and the residual is zero, not merely small.

**Result.** 0; (0, 0, 0); 0; (−2z, −2x, −2y); fluxes 0, 0; 15 and 15, residual 0.

**Check.** Dimensions: with coordinates in m and the coefficients in 1/s, **v** is in m/s and ∇ · **v** in 1/s. ✓ ∇ · (∇ × **F**) = ∂(−2z)/∂x + ∂(−2x)/∂y + ∂(−2y)/∂z = 0, as it must. ✓

**Common mistakes.** Computing ∇ × **F** with the wrong cyclic order (sign flips); believing a non-zero field must have non-zero flux through a cube.

## 2. Green for a rectangle, and the continuity equation

**Derivation.** (a) Q term: for each fixed y, the fundamental theorem gives ∫_a^b ∂Q/∂x dx = Q(b, y) − Q(a, y). Integrate over y from c to d: ∬ ∂Q/∂x dA = ∫_c^d Q(b, y) dy − ∫_c^d Q(a, y) dy. Counter-clockwise, the right side x = b is traversed with y increasing (dy > 0), giving the first term as ∫ Q dy on that side; the left side x = a is traversed with y decreasing, giving the second term with its minus sign; on the top and bottom dy = 0 and Q dy contributes nothing. Hence ∬ ∂Q/∂x dA = ∮ Q dy. P term: ∬ ∂P/∂y dA = ∫_a^b [P(x, d) − P(x, c)] dx. The bottom y = c is traversed with x increasing, contributing ∫_a^b P(x, c) dx to ∮ P dx; the top with x decreasing, contributing −∫_a^b P(x, d) dx. So ∮ P dx = −∬ ∂P/∂y dA. Subtracting the two results: ∮ (P dx + Q dy) = ∬ (∂Q/∂x − ∂P/∂y) dA. Orientation enters exactly twice: in the sign of dy on the vertical sides and of dx on the horizontal ones.

(b) Mass in a fixed V: M(t) = ∭_V ρ dV. Mass crosses S at the rate ∯_S ρ**v** · **n** dA (outward positive), so dM/dt = −∯ ρ**v** · **n** dA. *Divergence theorem step:* the right side equals −∭_V ∇ · (ρ**v**) dV. Because V is fixed, dM/dt = ∭ ∂ρ/∂t dV. Hence ∭_V [∂ρ/∂t + ∇ · (ρ**v**)] dV = 0. *Arbitrary-volume step:* this holds for every V; if the integrand were positive at some point it would be positive in a small ball around it (continuity), and the integral over that ball would be positive, a contradiction. So ∂ρ/∂t + ∇ · (ρ**v**) = 0 pointwise.

(c) Take the control volume bounded by two cross-sections A₁, A₂ and the tube wall between them. On the wall **v** is tangent to the surface, so **v** · **n** = 0 and the wall carries no flux. Steady flow: ∂ρ/∂t = 0, so ∯ ρ**v** · **n** dA = 0, i.e. −ρ₁v₁A₁ + ρ₂v₂A₂ = 0 (inflow counted negative at A₁ where **v** opposes **n**). Thus ρvA is the same at every section: the mass flow rate.

**Check.** Apply (a) to P = −y, Q = x on [0, 1]²: left side = circulation of a rigid rotation = 2 × area = 2; right side = ∬ (1 + 1) dA = 2. ✓

**Common mistakes.** Losing the sign of the P term; taking d/dt inside a *moving* volume (that needs the Reynolds transport theorem, Phase 05).

## 3. Weighing the Earth with a flux

**Assumptions.** Uniform density; spherical; Newtonian gravity.

**Derivation.** (a) **g** = c**r** with c = −(4/3)πGρ; ∇ · **g** = c(1 + 1 + 1) = −4πGρ. Divergence theorem: ∯ **g** · **n** dA = ∭ (−4πGρ) dV = −4πGρV = −4πG M_enclosed.

(b) ∇ · **g** = −4π × 6.674 × 10⁻¹¹ × 5515 = −4.625 × 10⁻⁶ s⁻². Through a 1 km cube (10⁹ m³): flux = −4625.3 m³/s². `box_flux` with the cube at x from 3000 to 3001 km, n = 4: −4625.318 m³/s² (exact: the field is linear). Divide by −4πG: 5.515 × 10¹² kg, which is ρ × 10⁹ m³ = 5.515 × 10¹² kg. ✓

(c) At the surface g = (4/3)πGρR, so ρ = 3g/(4πGR) = 3 × 9.82/(4π × 6.674 × 10⁻¹¹ × 6.371 × 10⁶) = 5514 kg/m³, 2.04 times the density of crustal rock. The interior must be much denser than the crust: iron-nickel core at 10 000–13 000 kg/m³.

(d) Uniform model: g(r) = (4/3)πGρr = 9.82 × 3371/6371 = 5.20 m/s² at r = 3371 km. The real value, about 10 m/s², is higher because the mass is concentrated toward the centre: at r = 3371 km the real Earth still encloses most of its mass (the core alone is a third of it) within a sphere whose radius is 53 % of R, so GM(r)/r² is larger than the uniform model's. The failed assumption is uniform density.

**Result.** ∇ · **g** = −4.63 × 10⁻⁶ s⁻²; flux −4625 m³/s², mass 5.515 × 10¹² kg; ρ̄ = 5514 kg/m³; g(3371 km) = 5.20 m/s² in the uniform model.

**Check.** Dimensions of ∇ · **g**: (m/s²)/m = s⁻²; of G ρ: (m³ kg⁻¹ s⁻²)(kg m⁻³) = s⁻². ✓ The mean density 5514 kg/m³ against the accepted 5515 kg/m³ is the classic result of Cavendish's measurement of G.

**Common mistakes.** Forgetting the factor 3 in the divergence of **r**; using R in km.

## 4. Residuals and where they come from

**Result.** (a) Smooth field, residuals for n = 2, 4, 8, 16, 32, 64: 3.978 × 10⁻⁴, 9.936 × 10⁻⁵, 2.484 × 10⁻⁵, 6.209 × 10⁻⁶, 1.552 × 10⁻⁶, 3.880 × 10⁻⁷. Fitted slope −2.00. `tplquad` gives 0.2384032 for the volume integral; `box_volume_integral` at n = 64 gives 0.2384036; the difference 4 × 10⁻⁷ is the midpoint rule's error at n = 64, consistent with the residual column.

(b) (x²y, yz², z sin x): residual 0.0 at n = 2, 4, 8 (rounding). (x³, 0, 0): 1.60 × 10⁻³, 4.00 × 10⁻⁴, 1.00 × 10⁻⁴: second order again. Proof for the x component: fix a cell column at (y_j, z_k) with cell centres x_i = x₀ + (i + ½)Δx. The volume sum contributes Σ_i ∂_x F_x(x_i, y_j, z_k) Δx · Δy Δz; the two x faces contribute [F_x(x₁, y_j, z_k) − F_x(x₀, y_j, z_k)] Δy Δz. The two agree exactly iff the midpoint rule integrates ∂_x F_x exactly over [x₀, x₁], which it does for any function at most linear in x (the midpoint rule is exact for polynomials of degree ≤ 1). If F_x is at most quadratic in x, ∂_x F_x is at most linear: exact. For x³, ∂_x F_x = 3x² is quadratic and the midpoint rule's error per cell is (Δx³/24) f″ = (Δx³/24) × 6 = Δx³/4; summed over the n cells of the column it is Δx²(x₁ − x₀)/4 = (0.4/n)² × 0.1 = 0.016/n² per unit of yz area. Multiplying by the x-face area 0.4 × 1 = 0.4 m² gives 0.0064/n²: 1.6 × 10⁻³ at n = 2, 4.0 × 10⁻⁴ at n = 4, 1.0 × 10⁻⁴ at n = 8, exactly the measured residuals.

(c) Rankine vortex on the 6 m square, errors against Γ = 300 for n = 25, 50, 100, 200, 400, 800. Loop: 5.09 × 10⁻², 1.27 × 10⁻², 3.18 × 10⁻³, 7.96 × 10⁻⁴, 1.99 × 10⁻⁴, 4.97 × 10⁻⁵: slope −2.00, clean. Area integral: 4.20, 4.20, 1.91, 0.230, 0.115, 0.029 m²/s: slope about −1.3 by least squares but the sequence is irregular (no change from n = 25 to 50, a jump of 8× from 100 to 200). The vorticity is a step at r = a; a cell is counted fully in or fully out, so the error is ω_z times the mis-assigned area along the circle, which depends on how the grid lines happen to fall relative to r = 1.5 m. First order on average, noisy in detail. The loop lies entirely outside the core, where the field is smooth, so it is untouched.

**Check.** The expected per-cell error formula in (b) reproduces the measured cubic residual to three figures, which is as strong a check of the midpoint rule's error constant as one can make.

**Common mistakes.** Reporting the quadratic field's residual as "below tolerance" rather than recognising it as exactly zero; fitting one slope through the noisy area column and quoting it to three figures.

## 5. Stokes on two surfaces and one closed surface

**Result.** (a) R = 1.3 m, **F** = (y eˣ, z sin y, x cos z):

| n | rim | cap | disc | cap − disc |
|---|---|---|---|---|
| 10 | −6.4781 | −6.5397 | −6.5060 | −3.37 × 10⁻² |
| 20 | −6.5040 | −6.5194 | −6.5110 | −8.41 × 10⁻³ |
| 40 | −6.5105 | −6.5144 | −6.5123 | −2.10 × 10⁻³ |
| 80 | −6.5122 | −6.5131 | −6.5126 | −5.25 × 10⁻⁴ |

(b) cap − disc falls by a factor 4.01, 4.00, 4.00: second order. In the continuum it is exactly zero because cap ∪ (disc with reversed normal) is a closed surface and ∯ (∇ × **F**) · **n** dA = ∭ ∇ · (∇ × **F**) dV = 0.

(c) On the hemisphere **n** = **r**/R, so (0, 0, 2ω) · **n** = 2ωz/R. ∬ z dA = ∫₀^{2π}∫₀^{π/2} R cos s · R² sin s ds dt = 2πR³ [sin² s/2]₀^{π/2} = πR³. Flux = (2ω/R)(πR³) = 2πωR², equal to the rim circulation ωR · 2πR. ✓

(d) With s and t swapped, **r**_s × **r**_t reverses and the cap flux would come out as +6.513 instead of −6.513. You detect it without knowing the answer because the disc and rim still give −6.51: three independent computations of the same quantity must agree, and the odd one out is the wrong one. A cheaper check: the flux of the position field **r** through the cap must be positive (outward).

**Check.** All three columns converge to −6.5127 (extrapolating the disc with Richardson: −6.5126 + (−6.5126 + 6.5123)/3 ≈ −6.5127).

**Common mistakes.** Using the lower hemisphere with an outward normal (it points −z at the pole and gives the negative of the rim circulation); comparing cap and disc at *different* n and calling the difference a bug.

## 6. When "flux in equals flux out" fails

(a) Steady, 40 m/s (Mach 0.12): ρ uniform to 0.7 % (½M² compressibility), ∂ρ/∂t = 0. Both volume and mass flux balance to within 1 %. Sign of trouble: a face-by-face imbalance larger than the rake's 1–2 % uncertainty.
(b) Start-up: ∂ρ/∂t ≠ 0 inside the volume (pressure waves settle in a few tunnel lengths over the speed of sound, i.e. tens of milliseconds, but the velocity field takes seconds to establish), and the unsteady term is responsible. Neither flux balances instant by instant; a measurable sign is the test-section static pressure still drifting.
(c) Mach 0.8: steady, so mass flux balances, but ρ varies by about 25 % between faces, so *volume* flux does not; ∇ · **v** ≠ 0 while ∇ · (ρ**v**) = 0. Sign: a pitot rake calibrated for incompressible flow shows a volume imbalance of order 20 % that disappears when density from the static-pressure and temperature probes is included.
(d) Leak: the control surface is not closed and impermeable; 2 % of the mass flow leaves through a face the rake does not cover. Mass flux through the instrumented faces does not balance, by 2 %, steadily. Sign: the imbalance is constant in time and on one specific face; sealing the hatch removes it.

**Common mistakes.** Blaming compressibility in (b) or unsteadiness in (c).

## 7. Volume and centroid from a surface alone

**Derivation.** (a) ∇ · (**r**/3) = 1, so ∯ (**r**/3) · **n** dA = ∭ 1 dV = V. ∇ · (x²/2, 0, 0) = x, so ∯ (x²/2) n_x dA = ∭ x dV = V x_c.

(b) Ellipsoid a = 2, b = 1, c = 0.5 m, exact V = (4/3)π abc = 4.18879 m³. `surface_flux` of **r**/3: 4.2061 (n = 10), 4.1931 (20), 4.1899 (40), 4.1891 (80); errors 1.73 × 10⁻², 4.31 × 10⁻³, 1.08 × 10⁻³, 2.69 × 10⁻⁴: order 2.00.

(c) Shifted to (1, 0, 0): V = 4.1891 (unchanged, as it must be) and ∯ (x²/2) n_x dA / V = 0.99994 m at n = 80: the centroid is at x = 1 to 6 × 10⁻⁵.

(d) For a flat triangle the normal area vector is **n** dA → **A** = ½(**p**₂ − **p**₁) × (**p**₃ − **p**₁), and **r** · **A** is the same for every **r** in the triangle's plane (the plane's equation), equal to **p**₁ · **A**. So the triangle contributes (1/3) **p**₁ · **A** = (1/6) **p**₁ · [(**p**₂ − **p**₁) × (**p**₃ − **p**₁)] = (1/6) **p**₁ · (**p**₂ × **p**₃), the signed volume of the tetrahedron with the origin. Reversing one triangle's vertex order flips the sign of its term: the computed volume is off by twice that tetrahedron's volume, which is why mesh checkers test for consistently oriented facets.

**Result.** V = 4.189 m³ (exact 4.18879), order 2; x_c = 1.0000 m; triangle term (1/6) **p**₁ · (**p**₂ × **p**₃).

**Check.** For a sphere a = b = c = R the parametrisation is `sphere_surface` of `P01.11` and the flux of **r**/3 is (1/3) · R · 4πR² = (4/3)πR³. ✓ The centroid flux is independent of the origin only when the surface is closed; shifting the ellipsoid changed ∯ (x²/2) n_x dA from 0 to V, exactly as V x_c requires.

**Common mistakes.** Forgetting the 1/3; using an outward normal on some facets and inward on others (the signed volumes then partly cancel).

## Check-yourself answers

1. Divergence: 1/s, the fractional rate of volume change of a fluid parcel (outward flux per unit volume). Curl: 1/s, twice the local angular velocity (circulation per unit area).
2. 6 s⁻¹, along the rotation axis.
3. The outward flux of a field through a closed surface (surface integral) equals the integral of its divergence over the enclosed volume (volume integral); the normal points out of the volume.
4. Symbolic: ∇ · (∇ × **F**) = 0, so the divergence theorem gives zero. Stokes: a closed surface has no boundary curve, so the circulation around its (empty) rim is zero.
5. Circulation measures the vorticity enclosed by the loop, not the speed along it, and outside the core the vorticity is zero.
6. (x, y, z) and (x², y², z²) make it exact (each component at most quadratic in its own coordinate); (x³, 0, 0) and (y eˣ, 0, 0) do not.
7. ρ**v** is divergence-free in steady flow; the mass flow rate ρvA is constant along the nozzle, not the volume flow rate.
8. ∇ · (∇ × **F**) = 0, combined with the divergence theorem applied to the closed surface formed by the two caps.
