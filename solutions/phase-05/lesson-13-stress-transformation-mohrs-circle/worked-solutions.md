# P05.13  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe script that calls it, given to the figures shown. Material values are typical, not design allowables. The stress states in Problems 3, 4 and 6 are illustrative and the rosette readings in Problem 5 are synthetic: none is a measurement.

## 1. Mohr's circle by hand (s = 15: σ_x = 70, σ_y = −20, τ_xy = −30 MPa)

**Assumptions.** Plane stress; tension positive; θ anticlockwise from x; τ of the x face plotted upward.

(a) X = (70, −30), Y = (−20, +30). C = **25.00 MPa**, R = √(45² + 30²) = **54.08 MPa**. σ₁ = **79.08 MPa**, σ₂ = **−29.08 MPa**. θ_p = ½ atan2(−60, 90) = **−16.85°**: the σ₁ direction is 16.85° *clockwise* from x. On the circle X lies below the σ axis at 33.7° (= 2 × 16.85°) from the σ₁ point.

(b) Maximum in-plane shear **54.08 MPa** on planes at θ_p − 45° = **−61.85°** (and at +28.15°), carrying the normal stress C = **25.00 MPa**. Because σ₁ > 0 > σ₂, the third principal stress (0) lies between them, so the largest circle is the in-plane one: the absolute maximum shear is also **54.08 MPa**.

(c) At θ = 30°: σ_x′ = 25 + 45 cos 60° − 30 sin 60° = **21.52 MPa**, σ_y′ = **28.48 MPa**, τ_x′y′ = −45 sin 60° − 30 cos 60° = **−53.97 MPa**. The point (21.52, −53.97) is 60° clockwise from X, close to the bottom of the circle.

(d) `principal_stresses_2d(70, -20, -30)` gives (79.08, −29.08, −16.85°); `eigvalsh` of [[70, −30], [−30, −20]] gives 79.083 and −29.083 MPa.

**Independent check.** σ_x′ + σ_y′ = 50.00 = σ_x + σ_y, and (σ_x′ − C)² + τ_x′y′² = 3.48² + 53.97² = 2925 = R².

**Common mistakes.** Plotting Y at (σ_y, +τ_xy) so that X and Y are not opposite; rotating the wrong way on the circle; forgetting the factor 2 between the circle angle and the body angle.

## 2. Derivations

(a) The wedge (lesson, 2.1): faces OB (the −x face, length dA cos θ) and OA (the −y face, length dA sin θ), cut face AB of length dA with outward normal n = (cos θ, sin θ). Forces on OB: (−σ_x, −τ_xy) dA cos θ; on OA: (−τ_xy, −σ_y) dA sin θ. Along n:
σ_x′ dA = σ_x cos²θ dA + τ_xy sin θ cos θ dA + τ_xy sin θ cos θ dA + σ_y sin²θ dA,
so σ_x′ = σ_x cos²θ + σ_y sin²θ + 2τ_xy sin θ cos θ, which with the double-angle identities is (1). Along the face direction (−sin θ, cos θ):
τ_x′y′ = −σ_x sin θ cos θ + σ_y sin θ cos θ + τ_xy(cos²θ − sin²θ) = −(σ_x − σ_y)/2 · sin 2θ + τ_xy cos 2θ.
Only **equilibrium** is used; no compatibility and no constitutive law, so (1) holds for any material.

(b) With (σ_x − σ_y)/2 = R cos φ₀, τ_xy = R sin φ₀: σ_x′ − C = R cos(φ₀ − 2θ), τ_x′y′ = R sin(φ₀ − 2θ); squaring and adding gives (3). dσ_x′/dθ = −(σ_x − σ_y) sin 2θ + 2τ_xy cos 2θ = 2τ_x′y′, which is zero (extreme normal stress) exactly where the shear is zero: **tan 2θ_p = 2τ_xy/(σ_x − σ_y)**. tan 2θ repeats every 90° in θ, so the two principal directions are perpendicular. dτ_x′y′/dθ = −(σ_x − σ_y) cos 2θ − 2τ_xy sin 2θ = 0 gives tan 2θ_s = −(σ_x − σ_y)/(2τ_xy) = −1/tan 2θ_p, so 2θ_s = 2θ_p ± 90°: **θ_s = θ_p ± 45°**.

(c) Expanding det(σ − sI) along the first row:
(σ_x − s)[(σ_y − s)(σ_z − s) − τ_yz²] − τ_xy[τ_xy(σ_z − s) − τ_yzτ_xz] + τ_xz[τ_xyτ_yz − τ_xz(σ_y − s)]
= −s³ + I₁s² − I₂s + I₃, with I₁, I₂, I₃ as in the lesson. Put s = p + x, p = I₁/3: the deviator has zero trace, so x³ − J₂x − J₃ = 0 with J₂ = ½ s:s and J₃ = det s. With x = r cos α: r³cos³α − J₂ r cos α − J₃ = 0. Choosing r² = 4J₂/3 makes the first two terms (r³/4)(4cos³α − 3cos α) = (r³/4) cos 3α, so **cos 3α = 4J₃/r³ = (J₃/2)(3/J₂)^{3/2}**, and the three roots are α, α − 2π/3, α − 4π/3. A real symmetric matrix has three real eigenvalues; a depressed cubic has three real roots exactly when 4J₂³ ≥ 27J₃², which is |cos 3α| ≤ 1.

(d) The principal stresses are σ₁, σ₂ and 0; with σ₁ ≥ σ₂ ≥ 0 the largest circle runs from 0 to σ₁, radius **σ₁/2**. It acts on the two planes whose normals bisect the σ₁ direction and the skin normal z (at ±45° to the skin surface, containing the σ₂ direction): planes that slant through the thickness.

**Common mistakes.** Using areas dA instead of dA cos θ and dA sin θ; writing tan θ_p instead of tan 2θ_p; dropping the τ² terms from I₂.

## 3. A welded titanium pressure tank (s = 15: p = 3.5 MPa)

**Assumptions.** Thin-walled membrane; uniform shear from the twist; x axial, y hoop.

(a) σ_axial = pr/(2t) = **145.83 MPa**, σ_hoop = pr/t = **291.67 MPa**. r/t = 83.3, and p is 1.2 % of the hoop stress: the thin-wall model is good to about that order.

(b) Weld at 35°, normal at 125°: **σ_n = 224.90 MPa, |τ| = 75.36 MPa**. Without the twist: σ_n = 243.69 MPa, |τ| = 68.52 MPa. (The twist's sign relative to the helix decides whether it adds to or subtracts from the weld shear; with this sense it adds 7 MPa of shear and removes 19 MPa of normal stress.)

(c) σ₁ = **294.36 MPa** at **82.33°** from the axis, σ₂ = **143.14 MPa** (at 172.33°). A weld laid along either principal direction, at **82.33°** or **172.33°** (−7.67°) from the axis, would carry no shear.

(d) Outer surface: principal stresses 294.36, 143.14, 0: τ_abs = **147.18 MPa**. Inner surface: σ_r = −3.5 MPa: τ_abs = (294.36 + 3.5)/2 = **148.93 MPa**, 1.2 % more. The radial stress matters at the level of the thin-wall approximation itself, not more.

(e) Compare the **normal stress across the weld** with the weld's tensile allowable (with its weld efficiency), and check the weld shear separately. A longitudinal seam (normal along the hoop direction) carries the full hoop stress, 291.67 MPa here; the 35° helix carries 224.90 MPa across it, 23 % less, at the cost of 75 MPa of shear along it. Whether the trade pays depends on the weld's tensile and shear allowables, which must both be checked; spiral-welded tubes rely on the reduced normal stress. σ₁ is 0.334 of the typical F_ty of Ti-6Al-4V, so the parent plate is not critical; the weld and its heat-affected zone are.

**Independent check.** σ_n(125°) + σ_n(35°) = σ_x + σ_y = 437.50 MPa, and the principal stresses agree with `eigvalsh`.

**Common mistakes.** Putting the weld angle itself into (1) instead of the angle of its normal; using the in-plane maximum shear for the inner surface.

## 4. A 3-D state at a lug (s = 15: σ_x = 150, τ_yz = 30 MPa; synthetic)

(a) I₁ = **120.00 MPa**, I₂ = **−9625.00 MPa²**, I₃ = **−349 500 MPa³**. (7) gives σ₁ = **165.414 MPa**, σ₂ = **28.562 MPa**, σ₃ = **−73.976 MPa**; `eigh` agrees to 1.5 × 10⁻⁸ Pa.

(b) n₁ = (0.9509, 0.3063, −0.0436), n₂ = (−0.2706, 0.8916, 0.3631), n₃ = (0.1501, −0.3335, 0.9307), each parallel to `eigh`'s eigenvector (|dot| = 1 to 12 digits), and N Nᵀ = I to round-off.

(c) τ_abs = (165.414 + 73.976)/2 = **119.69 MPa**, on planes bisecting n₁ and n₃ with normal stress (σ₁ + σ₃)/2 = **45.72 MPa**.

(d) On n = (1, 1, 1)/√3: σ_n = **76.67 MPa**, τ_n = **94.02 MPa**, inside the Mohr region. It is **not** the octahedral plane: that plane is equally inclined to the *principal* directions, and carries σ_n = I₁/3 = 40.00 MPa. In the model's axes (1, 1, 1)/√3 picks up the off-diagonal terms: σ_n = (σ_x + σ_y + σ_z + 2τ_xy + 2τ_yz + 2τ_xz)/3.

(e) Rotated by 25° about z: σ_x′ = 163.04, σ_y′ = 16.96, σ_z′ = −60.00, τ_x′y′ = −17.04, τ_y′z′ = 35.64, τ_x′z′ = −5.45 MPa. Invariants: 120.000 MPa, −9625.000 MPa², −349 500.0 MPa³, unchanged. σ_z is unchanged because the rotation is about z.

**Common mistakes.** Taking (1, 1, 1)/√3 as octahedral in any axes; sorting the eigenvalues from `eigh` the wrong way (it returns them ascending).

## 5. A delta rosette (s = 15: 950, 200, −75 µε; synthetic readings)

(a) With (8) at 0°, 60° and 120°: ε₀ = ε_x; ε₆₀ = ε_x/4 + 3ε_y/4 + (√3/4)γ; ε₁₂₀ = ε_x/4 + 3ε_y/4 − (√3/4)γ. Adding the last two: ε_y = (2(ε₆₀ + ε₁₂₀) − ε₀)/3; subtracting: γ = 2(ε₆₀ − ε₁₂₀)/√3. Numbers: **ε_x = 950.0, ε_y = −233.3, γ_xy = 317.5 µε**, by `rosette_strains` and by the closed forms.

(b) σ_x = **71.62 MPa**, σ_y = **6.58 MPa**, τ_xy = **8.73 MPa**; σ₁ = **72.77 MPa**, σ₂ = **5.43 MPa** at **7.51°** to the 0° gauge.

(c) Principal strains **970.9 µε** and **−254.3 µε** at 7.51° (the same angle: isotropy). E/(1 − ν²)(ε₁ + νε₂) = **72.77 MPa** = σ₁. ε_z = −ν(ε_x + ε_y)/(1 − ν) = **−353.0 µε**: the skin thins, although no stress acts through its thickness.

(d) 10 µε noise, seed 15, 20 000 trials: delta rosette std(σ₁) = **0.78 MPa**, std(θ_p) = **0.38°**; rectangular 0/45/90 at the same strain state: **0.79 MPa**, **0.56°**. The magnitudes are equally well determined; the delta rosette, with its gauges spread evenly over 180°, determines the angle better here. (With 20 000 trials the standard deviations are known to about 0.5 %, so the difference in angle is real.)

**Common mistakes.** Using γ in place of γ/2 in the principal-strain radius; applying uniaxial σ = Eε to each gauge.

## 6. What drove this crack?

(a) Pressure alone: σ_x = 37.5 MPa (axial), σ_y = 75 MPa (hoop), no shear; σ₁ is **hoop** (θ_p = 90°), so a pressurisation fatigue crack runs **along the fuselage axis** (longitudinally).

(b) A crack at 25° to the axis needs σ₁ at 115° (equivalently −65°). τ = (σ_x − σ_y)/2 · tan 2θ_p = (−18.75) tan(−130°) = **−22.3 MPa** (magnitude 22.3 MPa); then σ₁ = 85.42 MPa, σ₂ = 27.08 MPa. Fuselage torsion of 100–200 kN·m gives τ = T/(2πr²t) = **2.5–5.0 MPa** in this skin, a factor of 4–9 too small. (Transverse fuselage shear is largest at the sides and near zero at the crown, `P05.11`.)

(c) With +40 MPa of bending, σ_x = 77.5 > σ_y = 75 MPa. On Mohr's circle X is then to the right of the centre, so the σ₁ direction is within ±45° of x for *any* τ (equation (4) with a positive denominator). The crack, perpendicular to σ₁, would run within 45° of the hoop direction, at least 45° from the axis. **No shear can produce a crack at 25° to the axis** in this state.

(d) The global membrane stresses at the crown do not explain a crack at 25°, so the colleague's explanation is doubtful as stated. Possibilities: (1) a **local stress field** at the rivet row (a stress concentration at the hole, a lap joint's secondary bending, load transfer through the fastener), which can rotate σ₁ locally; (2) **residual stresses** from riveting or cold working; (3) the edge of a **repair doubler**, antenna or cut-out nearby; (4) **multiple-site damage**, small cracks at adjacent holes linking along the rivet row, so the path follows the row rather than σ₁; (5) a load case other than pressure plus bending (a ground or door-surround load). Ask for: fractography (beach marks show the growth direction and origin), the structural drawing and repair history of that bay, an inspection of neighbouring holes, and a strain-gauge survey under a pressure test. Until then, treat the crack as unexplained and do not attribute it to routine pressurisation fatigue.

**Common mistakes.** Taking the crack to run *parallel* to σ₁; assuming any angle can be produced by some shear.

## Check-yourself answers

1. Moment equilibrium of the element about its centre: the shear forces on the x faces form a couple that only the shear on the y faces can balance, so τ_yx = τ_xy (with no body couples).
2. 50°, twice the body angle (in the opposite sense with the lesson's plotting convention).
3. C = 40 MPa, R = √(20² + 15²) = 25 MPa, σ₁ = 65 MPa, σ₂ = 15 MPa.
4. Both in-plane principal stresses are tensile, and the third (out of plane) is zero; the largest Mohr circle runs from 0 to σ₁, larger than the in-plane circle from σ₂ to σ₁. It acts on planes at 45° between the σ₁ direction and the skin normal, slanting through the thickness.
5. I₁ = trace, I₂ = sum of the principal 2 × 2 minors, I₃ = determinant. They are the coefficients of det(σ − sI), whose roots are the principal stresses, which are physical and so cannot depend on the axes.
6. σ₁ = +τ, σ₂ = −τ (and 0 out of plane), at ±45° to the shear axes.
7. At a double root cos 3α = ±1, where arccos has infinite slope: a round-off error ε in its argument becomes an error of order √ε in α. `eigh` works on the matrix with orthogonal transformations, which never amplify round-off that way; symmetric eigenvalues are perfectly conditioned.
8. γ_xy = 2ε₄₅ − ε₀ − ε₉₀ = 800 − 500 − 100 = 200 µε.
