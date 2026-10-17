# P02.12  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for $s = 15$; numerical results are from `solution.py` and given to the figures shown.

## 1. From a harmonic function to an analytic one ($s = 15$)

**Derivation.** (a) u = x³ − 3xy² + 15x: u_xx = 6x, u_yy = −6x, sum 0 ✓.

(b) v_y = u_x = 3x² − 3y² + 15 ⇒ v = 3x²y − y³ + 15y + g(x). Then v_x = 6xy + g′(x) must equal −u_y = 6xy, so g′ = 0 and g is a constant (take 0). Since (x + iy)³ = x³ − 3xy² + i(3x²y − y³), **f = z³ + 15z**.

(c) f′(z) = 3z² + 15; at z = 1 + i, (1 + i)² = 2i, so f′ = **15 + 6i**. `complex_derivative(lambda z: z**3 + 15*z, np.array([1+1j]))` returns 15 + 6i to 10⁻⁹; `cauchy_riemann_residual` returns 1.3 × 10⁻⁹.

(d) u − iv = conj(f) = z̄³ + 15z̄: a function of z̄, whose Cauchy–Riemann residual is not zero (for z̄ alone it is 2). Conjugating v flips the sign of v_y and v_x, so u_x = v_y becomes u_x = −v_y, which fails except where u_x = 0.

**Check.** v is also harmonic: v_xx = 6y, v_yy = −6y ✓. Level curves: ∇u·∇v = u_x v_x + u_y v_y = (3x² − 3y² + 15)(6xy) + (−6xy)(3x² − 3y² + 15) = 0 ✓.

**Common mistakes.** Dropping g(x); forgetting the +15 in f′.

## 2. The Kutta condition and the Kutta–Joukowski theorem

(a) With s = Re^{−iβ} in w = U[e^{−iα} − R²e^{iα}/s²] + iΓ/(2πs): the doublet term is −Ue^{iα}e^{2iβ}, the vortex term iΓe^{iβ}/(2πR). Multiply by e^{−iβ}: U[e^{−i(α+β)} − e^{i(α+β)}] + iΓ/(2πR) = −2iU sin(α + β) + iΓ/(2πR). Setting this to zero: **Γ = 4πUR sin(α + β)**. The sign is positive for α + β > 0 with the potential's vortex term written +iΓ/(2π) log s, which is a *clockwise* circulation (tangential speed Γ/(2πr) in the −θ direction). With the other convention (−iΓ/2π log s) the formula reads Γ = −4πUR sin(α + β) and positive lift corresponds to negative Γ; either is consistent, but mixing them gives a stagnation point on the wrong side.

(b) On the circle s = Re^{iθ}, the same algebra with e^{−iθ} in place of e^{iβ} gives 2iU sin(θ − α) + iΓ/(2πR) = 0, so sin(θ − α) = −Γ/(4πUR) = −sin(α + β). Solutions: θ − α = −(α + β), i.e. **θ = −β** (the trailing edge), and θ − α = π + (α + β), i.e. **θ = π + 2α + β** (the front stagnation point). At α = −β: θ = −β and θ = π − β, mirror images in the line through the centre at angle −β: the flow is symmetric about the chord line of the circle and there is no circulation.

(c) F_x + iF_y = −∮ p n ds. On a counter-clockwise contour the outward normal times the arc length is n ds = −i dz (rotate the tangent dz clockwise by 90°). So F_x + iF_y = i∮ p dz, and conjugating, F_x − iF_y = −i∮ p dz̄. Bernoulli: p = p_∞ + ½ρU² − ½ρ|w|²; the constant integrates to zero round a closed curve; so F_x − iF_y = (iρ/2)∮ |w|² dz̄ = (iρ/2)∮ w w̄ dz̄. On a streamline the flow is tangent to the contour, so w dz = (u − iv)(dx + i dy) = u dx + v dy + i(u dy − v dx) is real (u dy − v dx = 0), and therefore w̄ dz̄ = conj(w dz) = w dz. Hence w w̄ dz̄ = w² dz and **F_x − iF_y = (iρ/2)∮ w² dz**.

(d) w = Ue^{−iα} + iΓ/(2πz) + a₂/z² + …, so w² = U²e^{−2iα} + iΓUe^{−iα}/(πz) + O(z⁻²). The contour can be deformed to a large circle (w² is analytic outside the body); the only term with a non-zero integral is the 1/z term, with residue iΓUe^{−iα}/π. F_x − iF_y = (iρ/2)(2πi)(iΓUe^{−iα}/π) = −iρUΓe^{−iα} = −iρUΓ(cos α − i sin α) = −ρUΓ sin α − iρUΓ cos α. So F_x = −ρUΓ sin α, F_y = ρUΓ cos α: the force is perpendicular to the free stream with magnitude **ρUΓ** and no component along it (**D′ = 0**). The shape of the body never entered: only the far-field circulation, which is a contour integral and therefore the same for every contour enclosing the body.

**Common mistakes.** Using n ds = +i dz (inward normal); forgetting that |w|² is not analytic and must be converted to w² via the streamline property; treating the O(z⁻²) terms as contributing.

## 3. A wake vortex near the ground ($s = 15$: Γ = 350 m²/s, h = 25 m)

**Derivation.** (a) Clockwise vortex Γ at ih, image counter-clockwise at −ih: F = (iΓ/2π)[log(z − ih) − log(z + ih)], w = (iΓ/2π)[1/(z − ih) − 1/(z + ih)] = (iΓ/2π)(2ih)/(z² + h²) = **−Γh/(π(z² + h²))**.

(b) On the ground z = x is real, so w = −Γh/(π(x² + h²)) is real: v = 0 everywhere on the runway (the ground is a streamline) and **u = −Γh/(π(x² + h²))**, a flow toward −x for a clockwise vortex (the vortex's lower side moves to the left, and the image doubles it). Maximum |u| = **Γ/(πh) = 350/(π × 25) = 4.46 m/s** at x = 0, directly below the vortex; the code's scan of x ∈ [−20, 20] m finds 4.456 m/s at x = 0. (The velocity directly under a single vortex at height h would be Γ/(2πh); the image doubles it.)

(c) The real vortex is carried by the velocity of its image evaluated at ih: w = −(iΓ/2π)/(ih + ih) = −Γ/(4πh), real, so the drift is horizontal at **Γ/(4πh) = 1.11 m/s**, toward −x for clockwise Γ (and toward +x for the counter-clockwise partner of a real vortex pair: the two vortices of a wake separate as they approach the ground, which is why a crosswind of about 1 m/s can hold one of them over the runway).

(d) `circulation_and_flux(w, circle_path(25j, 5.0))` returns (−350.0, 0): the counter-clockwise contour integral of a clockwise vortex is −Γ. `circle_path(0, 60.0)` enclosing both returns (0, 0): the circulations cancel and neither vortex is a source.

**Check.** Dimensions: Γ/h = m²/s ÷ m = m/s ✓. Limit h → ∞: ground velocity → 0 and the drift → 0 ✓.

**Common mistakes.** Taking the image with the same sign (then the ground is not a streamline); including a vortex's own (infinite) self-induced velocity in its motion.

## 4. Your own section ($s = 15$: x₀ = −0.021875 m, y₀ = 0.0125 m)

| quantity | value |
|---|---|
| R | 0.27216 m |
| β | 2.632° |
| chord | 1.0065 m |
| thickness / chord | 10.5 % |
| camber / chord | 2.3 % |
| Γ at 4° | 19.751 m²/s |
| L′ | 1209.7 N/m |
| C_l at 4° | 0.7849 |

(c) Lift slope from 0° and 8°: **6.745 per radian**; 2π(1 + 0.77 × 0.1045) = 6.789 (0.6 % above the code). Zero-lift angle **−2.63°**.

(d) Pressure-integrated lift 1209.7 N/m against ρUΓ = 1209.7 N/m (agreement to 10⁻⁵); drag 10⁻⁴ N/m, below 10⁻⁷ of the lift.

**Check.** C_l = 2π sin(α + β)(4R/chord) = 2π × sin 6.632° × 1.0816 = 0.785 ✓. Thinner and less cambered than the lesson section, so both the slope and C_l(4°) are lower ✓.

**Common mistakes.** Using chord = 4c instead of the measured chord (0.7 % error here); mixing degrees and radians in the slope.

## 5. How fast is too fast? ($s = 15$: α = 3.75°)

(a) Γ = **26.975 m²/s**, L′ = **1652.2 N/m**, C_l = **1.0700**; pressure-integrated lift 1652.2 N/m.

(b) Maximum surface speed **78.05 m/s** (1.561U) at z = (−0.472 + 0.041i) m, on the upper surface 3.6 % of chord behind the leading edge; minimum C_p = 1 − 1.561² = **−1.437**.

(c) M_local = 78.05/340 = **0.23**. Since |w| ∝ U, the local Mach number reaches 0.3 when 1.561U = 102 m/s, i.e. at **U = 65.3 m/s**. Above that, the incompressible C_p is in error by more than about 5 % at the suction peak (the Prandtl–Glauert factor 1/√(1 − 0.3²) = 1.048).

(d) At ζ = c(1 + 10⁻⁶) the code returns 44.717 m/s, the same as the l'Hôpital limit at c itself to five figures: both w_ζ (8.9 × 10⁻⁵ m/s) and dz/dζ (2.0 × 10⁻⁶) vanish *linearly* with the distance from c, so their ratio is well-conditioned, provided Γ is exactly the Kutta value. Perturb Γ by 1 % and w_ζ(c) becomes about 0.01 × 2U sin(α + β) = 0.15 m/s, which divided by 2 × 10⁻⁶ gives 8 × 10⁴ m/s: the "finite trailing-edge velocity" exists only for one value of Γ, and a numerical Kutta condition must be imposed exactly (as a constraint), not approximately. Panel methods do precisely that, adding one equation to the linear system.

**Common mistakes.** Reporting the maximum speed at a non-surface point; forgetting that C_p min depends on U² but M_local on U.

## 6. What the model may and may not say

(a) **Viscosity**: boundary-layer displacement de-cambers the section and moves the effective Kutta point; lift too high by 5–15 % at 4°, and drag (zero in the model) is C_d ≈ 0.007 in reality. **Compressibility**: at U = 50 m/s, M_local ≤ 0.23, lift low by about 2 % (Prandtl–Glauert 1/√(1 − M²) with M = 0.15 free stream); unaffected for practical purposes. **Finite span**: the model is per unit span of an infinite wing; a wing of aspect ratio 8 has a lift slope about 20 % lower (Phase 10), so for a *wing* the number is too high by that much. **Finite trailing-edge angle**: a real 10–15° edge angle lowers the lift slope by about 1–2 % relative to the cusp; small. Stall is a viscous effect that is "unaffected" at 4° and decisive above 12°.

(b) The pressure recovery on the rear 30 % of the upper surface: measured C_p stays more negative than the potential curve (the boundary layer cannot follow the steep adverse gradient and thickens), and C_p at the trailing edge does not reach the model's value. At 4° this is a few per cent of the suction peak; it is the first visible sign of viscosity and the quantity a boundary-layer code must predict.

(c) No. C_l = 2.0 at 12° is above the maximum lift of any plain 12 %-thick section (measured C_l,max ≈ 1.4–1.6 at Reynolds numbers of a few million); the section has either stalled or is about to. Use measured section data (a wind-tunnel polar for a similar thickness and camber) or a viscous airfoil code, and treat 12° as near C_l,max, not on the linear slope.

(d) "The five-figure agreement between the pressure integral, ρUΓ and the contour integral shows that we have solved the potential-flow equations correctly: that is *verification*. Whether potential flow is the right model for this wing at this Reynolds number is a question about the equations, not about our arithmetic, and only a measured pressure distribution or a measured lift curve can answer it: that is *validation*. The model is verified; it is validated only over the front half of the section, at low incidence, to a few per cent."

## Check-yourself answers

1. u_x = v_y, u_y = −v_x. Then u_xx = v_yx = v_xy = −u_yy, so u_xx + u_yy = 0.
2. φ is the velocity potential ((u, v) = ∇φ), ψ the stream function (u = ψ_y, v = −ψ_x). Incompressibility makes φ harmonic and ψ exist; irrotationality makes ψ harmonic and φ exist; together they are the Cauchy–Riemann equations for F = φ + iψ. dF/dz = u − iv.
3. F = U[z e^{−iα} + R²e^{iα}/z] + (iΓ/2π) log z: the first term is the free stream, the second the doublet that makes the circle a streamline, the third the vortex that gives circulation (and lift).
4. It maps |ζ| = c onto the segment [−2c, 2c] traversed twice (a flat plate of chord 4c). At ζ = ±c the derivative vanishes and angles are doubled: a smooth curve through c acquires a cusp.
5. The rear stagnation point must lie at the trailing-edge point of the circle, so that the mapped flow leaves the sharp edge smoothly; Γ = 4πUR sin(α + β).
6. w_z = w_ζ/(dz/dζ) is 0/0 there: the Kutta circulation makes w_ζ vanish linearly at c, as dz/dζ does, and l'Hôpital gives the finite ratio (dw_ζ/dζ)/(d²z/dζ²).
7. F_x − iF_y = (iρ/2)∮ w² dz; on a large circle w = Ue^{−iα} + iΓ/(2πz) + …, so w² has residue iΓUe^{−iα}/π at the 1/z term; the residue theorem gives F_x − iF_y = −iρUΓe^{−iα}, a force of magnitude ρUΓ perpendicular to the stream: L′ = ρUΓ, D′ = 0.
8. C_l = 2π sin(8°) × 1.10 = 0.962; zero-lift angle −3°.
