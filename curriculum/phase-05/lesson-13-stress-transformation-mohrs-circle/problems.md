# P05.13  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical. Draw the stress element (and, where it helps, Mohr's circle) before writing any equation. Tension is positive; θ is anticlockwise from x.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

Material data (typical, not design allowables; publicly viewable ASM Aerospace Specification Metals data sheets): 2024-T3 E = 73.1 GPa, ν = 0.33; Ti-6Al-4V (annealed) F_ty = 880 MPa.

## 1. (Warm-up: Mohr's circle by hand)

A point in a wing skin carries σ_x = 40 + 2s MPa, σ_y = −20 MPa, τ_xy = −30 MPa.

(a) Draw Mohr's circle to scale by hand. Read off the centre, the radius, σ₁, σ₂ and the angle from x to the σ₁ direction (say in which sense you rotate).
(b) Find the maximum in-plane shear, the angle of its plane, and the normal stress on that plane. Is the absolute maximum shear the same, and why?
(c) Compute σ_x′, σ_y′ and τ_x′y′ on axes rotated by 30° with (1), and mark the point on your circle.
(d) Confirm (a)–(c) with your code and with `numpy.linalg.eigvalsh`.

Answer form: a sketch, five numbers; three numbers and a sentence; three numbers; a check.

## 2. (Derivation: the wedge, the circle and the cubic)

(a) Derive both equations (1) of the lesson from the equilibrium of a wedge, with a labelled free-body diagram. Which of equilibrium, compatibility and the constitutive law did you use?
(b) Eliminate θ from (1) to obtain (3), Mohr's circle. Then show that the angle 2θ_p of the principal direction satisfies tan 2θ_p = 2τ_xy/(σ_x − σ_y), that the two principal directions are perpendicular, and that the maximum in-plane shear planes lie at 45° to them.
(c) Show that the characteristic equation det(σ − sI) = 0 of a 3 × 3 symmetric tensor is s³ − I₁s² + I₂s − I₃ = 0 with the invariants of the lesson. Then remove the mean stress p and substitute x = r cos α into x³ − J₂x − J₃ = 0 to derive (7). Why is the argument of the arccos never outside [−1, 1] in exact arithmetic?
(d) For plane stress with σ₁ ≥ σ₂ ≥ 0, show that the absolute maximum shear is σ₁/2 and find the planes on which it acts.

Answer form: four derivations, a named ingredient, two sentences.

## 3. (Unfamiliar scenario: a welded titanium pressure tank)

A helium pressurant tank for a launcher upper stage is a Ti-6Al-4V cylinder of radius r = 0.25 m and wall t = 3.0 mm, made by rolling a plate and welding a seam along a helix. Its operating pressure is p = 2 + s/10 MPa. A misaligned mounting strap also twists it, giving a uniform wall shear of 20 MPa (illustrative).

(a) Compute the axial and hoop stresses. Is the thin-wall model reasonable (compare r/t, and p with the hoop stress)?
(b) The weld runs at 35° to the tank axis. Compute the normal stress across the weld and the shear along it, with and without the twist.
(c) Find the principal stresses and their directions. At which angles could the weld have been placed to carry no shear?
(d) Find the absolute maximum shear at the outer surface (radial stress 0) and at the inner surface (radial stress −p). Does the radial stress matter?
(e) The weld's strength is often lower than the parent plate's. Which stress on the weld would you compare with a weld allowable, and why might a helical seam be preferred to a longitudinal one?

Answer form: two stresses and a sentence; four stresses; two stresses, an angle and two angles; two stresses and a sentence; a paragraph.

## 4. (Numerical, uses your code: a 3-D state at a lug)

A finite-element model of a wing–fuselage lug reports, at one element, the (synthetic, illustrative) stress tensor in MPa:

σ_x = 120 + 2s, σ_y = 30, σ_z = −60, τ_xy = 45, τ_yz = 15 + s, τ_xz = −20.

(a) Compute the three invariants and the principal stresses by `principal_stresses_3d`. Compare with `numpy.linalg.eigh` and report the largest difference.
(b) Compute the three principal directions by `principal_direction`, check that they are mutually perpendicular, and compare them with `eigh`'s eigenvectors.
(c) Compute the absolute maximum shear and the normal stress on its plane.
(d) Compute the traction (σ_n, τ_n) on the plane with normal (1, 1, 1)/√3 in the model's axes, and check that it lies in the Mohr region. Is this the octahedral plane? Why or why not?
(e) Rotate the tensor by 25° about z with `rotation_z` and `transform_tensor`. Show that the invariants are unchanged.

Answer form: three invariants, three stresses and a difference; three vectors and two checks; two stresses; two stresses, a check and a sentence; a check.

## 5. (Numerical, uses your code: a delta rosette)

A delta (0°/60°/120°) rosette on a 2024-T3 skin panel gives the synthetic readings ε₀ = 800 + 10s, ε₆₀ = 200, ε₁₂₀ = −150 + 5s µε.

(a) Find ε_x, ε_y and γ_xy with `rosette_strains`. Derive the closed forms ε_x = ε₀, ε_y = (2(ε₆₀ + ε₁₂₀) − ε₀)/3, γ_xy = 2(ε₆₀ − ε₁₂₀)/√3 and check them.
(b) Compute σ_x, σ_y, τ_xy and the principal stresses and angle.
(c) Find the principal strains and show that applying (9) to them gives the same σ₁. What is the out-of-plane strain ε_z?
(d) Add independent Gaussian noise of 10 µε (1σ) to each reading (seed s, 20 000 trials) and report the standard deviations of σ₁ and θ_p. Repeat for a 0/45/90 rosette at the same strain state with the same noise. Which rosette pins down the angle better here?

Answer form: three strains and a derivation; five stresses and an angle; two strains, a check and a strain; four numbers and a sentence.

## 6. (Judgment: what drove this crack?)

A routine inspection finds a 12 mm fatigue crack growing from a rivet hole in the crown skin of a pressurised fuselage (r = 2.0 m, t = 1.6 mm, Δp = 60 kPa). The crack runs at **25° to the fuselage axis**. A colleague says it is a normal pressurisation fatigue crack.

(a) With pressure alone, which way does σ₁ act and which way should a fatigue crack (perpendicular to σ₁) run?
(b) Fatigue cracks grow perpendicular to the largest tensile stress. What shear stress τ_xy would turn σ₁ so that the crack would run at 25° to the axis, with pressure alone? Fuselage torsion of 100–200 kN·m (illustrative) gives what shear in this skin?
(c) With the crown also carrying a fuselage-bending stress of +40 MPa, can any shear produce a crack at 25° to the axis? Explain using Mohr's circle.
(d) What would you conclude, and what would you ask for before agreeing or disagreeing with your colleague? List at least three possibilities and the evidence that would separate them.

Answer form: a direction and a sentence; two stresses; a yes or no with an argument; a judgment with a list.
