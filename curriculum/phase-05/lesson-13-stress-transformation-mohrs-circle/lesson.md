# P05.13  Stress Transformation & Mohr's Circle

> The stress at a point does not change when you turn your head; only its components do. Turning the axes until the shear disappears tells you how the material is really loaded.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the plane-stress transformation equations from the equilibrium of a wedge, and show that they are the 2 × 2 block of the tensor rotation σ′ = RσRᵀ. 2. Construct Mohr's circle by hand and read from it the principal stresses, their directions, the maximum in-plane shear and the stresses on any plane. 3. Explain why the absolute maximum shear of a plane-stress state uses the third principal stress, zero, and compute it. 4. Find the principal stresses of a 3-D state from the invariants by the closed-form (trigonometric) solution of the characteristic cubic, check them against `numpy.linalg.eigh`, find the principal directions, and show that every plane's traction lies in the region bounded by the three Mohr circles. 5. Reduce a strain-gauge rosette's three readings to principal stresses, and estimate how gauge noise and misalignment propagate into them. |
| **Prerequisites** | `P02.14` (the stress tensor, its transformation law and the eigenproblem, in index notation), `P01.16` (eigenvalues of a symmetric matrix), `P01.15` (direction-cosine matrices), `P05.06` (Hooke's law with Poisson's ratio, the pressurised cylinder), `P05.08` and `P05.11` (shear stresses from torsion and transverse shear), `P05.10` (bending stress), `P00.13`. Skill check: given a 2 × 2 symmetric matrix you can find its eigenvalues by hand, and you can write the hoop and axial stresses of a thin pressurised cylinder. |
| **Engineering question** | A fuselage side panel 2 m from the axis, with 1.6 mm skin at a 60 kPa cabin pressure differential, carries an axial stress of 37.5 MPa, a hoop stress of 75 MPa and an (illustrative) shear of 30 MPa from fuselage shear and torsion. What is the largest normal stress, which way does it act, and what is the largest shear stress at the point? The answers: 91.6 MPa at 61.0° from the fuselage axis; a maximum in-plane shear of 35.4 MPa, but an absolute maximum of 45.8 MPa on planes inclined through the skin's thickness. |
| **Assumptions and limits** | Stress is a continuous field and is uniform over the small element considered; the stress tensor is symmetric (moment equilibrium of the element, no body couples); plane stress in thin skins (σ_z = τ_xz = τ_yz = 0); linear isotropic elasticity for the strain-gauge reduction only (the transformation itself holds for any material, elastic or not); the gauge reads the average strain under its grid. |
| **Analytical work** | Wedge equilibrium and the transformation equations; Mohr's circle and its angle-doubling; principal stresses and directions, maximum in-plane shear at 45°; the third principal stress of plane stress and the absolute maximum shear; σ′ = RσRᵀ; invariants and the characteristic cubic; the trigonometric solution of the cubic; principal directions by cross products; the traction on a plane and the three-circle Mohr region; strain transformation and the rosette; plane-stress Hooke's law. |
| **Computational deliverable** | `starter/main.py`: `transform_stress_2d`, `mohr_circle_2d`, `principal_stresses_2d`, `max_inplane_shear`, `plane_stress_abs_max_shear`, `transform_tensor`, `invariants`, `principal_stresses_3d`, `principal_direction`, `absolute_max_shear`, `traction_on_plane`, `mohr_circles_3d`, `in_mohr_region`, `rosette_strains`, `plane_stress_from_strain`, `fuselage_skin_stress` (`stress_tensor`, `rotation_z` and `main` are given). Tests prove: the 2-D formulas equal RσRᵀ and preserve the invariants; rotated stress points lie on Mohr's circle and a 180° rotation is the identity; principal planes carry no shear and hold the extreme normal stresses; maximum in-plane shear is at 45° to them; pure shear, uniaxial and hydrostatic limits; the absolute maximum exceeds the in-plane one when the principal stresses share a sign; the closed-form cubic matches `eigh` on 200 random tensors and on repeated roots; invariants survive rotation and give the characteristic polynomial; principal directions match `eigh`'s eigenvectors; every plane's traction lies in the Mohr region and the absolute maximum shear is reached at 45° between σ₁ and σ₃; the rosette inverts the strain transformation; Hooke's law round-trips. |
| **Practical activity** | **computational experiment**: a shear sweep on the fuselage side panel showing when the absolute maximum shear leaves the plane of the skin; a Monte Carlo comparison of the closed-form cubic with `eigh` and its loss of accuracy near repeated roots; 20 000 random planes plotted in the 3-D Mohr region; the propagation of gauge noise and misalignment through a rosette reduction. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `mohr_circle.png`, the wedge derivation and a hand-drawn Mohr's circle for the side panel, the four experiment tables, and a half-page note on the side panel's principal stresses and maximum shear. |
| **Assessment** | `assessments/phase-05/grading-rubric.md`; competencies 1–4. Pass: the transformation equations are derived from wedge equilibrium; Mohr's circle is drawn by hand and agrees with the code; all public tests pass; the absolute maximum shear is reported with the third principal stress named; the cubic is checked against `eigh` and its accuracy near repeated roots is explained. |
| **References** | `@mit-ocw-3-11-mechanics-of-materials` Roylance module "Tensor Transformations" (PDF "Transformation of Stresses and Strains", https://ocw.mit.edu/courses/3-11-mechanics-of-materials-fall-1999/resources/mit3_11f99_trans/): Introduction; Direct approach; Mohr's circle; General approach; Principal stresses and planes in three dimensions (required-free). `@mit-ocw-16-001-unified-materials-structures` lecture notes "Stress Equilibrium" and "Transformation of Stress Components" (required-free). `@tudelft-alderliesten-aerospace-structures-materials` Ch. 1 "Material physics & properties": 1.7 Bi-axial loading (required-free, read online, link only). `@nasa-tm-x-73305-astronautic-structures-manual` Vol. I, Section A3 "COMBINED STRESSES": 3.1.0 Combined Stresses (https://ntrs.nasa.gov/citations/19760071125) (supplementary-free): the design-manual form of the same formulas. The interactive Mohr's-circle tools that the 3.11 module page once linked no longer run in current browsers; this lesson's code replaces them. |
| **Extensions** | ★ Add the Newton refinement of each root of the cubic and show that one step restores full precision near a double root. ★ Find the principal directions of a state with a double root as an eigenspace (any vector in a plane) and handle it in `principal_direction`. ★ Write the strain-gauge reduction for an orthotropic composite skin, where principal strain and principal stress directions differ, and find the angle between them for a [0/90] laminate (preview of Phase 06). |
| **Time** | unmeasured |

## 1. The engineering problem

A pressurised fuselage is a thin cylinder. Cutting it along its length and across it (`P05.06`) gives a hoop stress Δp r/t and an axial stress Δp r/(2t); fuselage bending adds axial stress (`P05.10`), and fuselage shear and torsion add shear (`P05.08`, `P05.11`). At a point in the skin there are therefore three stress components in the axes we happened to choose, the fuselage's axial and circumferential directions. None of them is necessarily the largest normal stress at the point, and none is the largest shear.

That matters because cracks and yield do not care about our axes. A fatigue crack grows perpendicular to the largest tensile stress; a ductile metal yields when the largest shear reaches a limit (`P05.14`). A spirally wound tube, a skin splice at an angle, a fastener row cut diagonally, or a strain gauge glued at 45° all need the stress on some other plane. And an analyst handed a finite-element stress tensor with six components needs the three numbers that describe it independently of the axes: the principal stresses.

This lesson builds that tool. Its mathematics, the rotation of a symmetric tensor, was done in `P02.14`; here it becomes mechanics: where the formulas come from (equilibrium of a small wedge), how to see them (Mohr's circle), how to get the 3-D answer without a library (the closed-form cubic) and how to check it with one (`numpy.linalg.eigh`), and how to get the stress at a point on a real structure from three strain-gauge readings.

Success criteria, written first: the 2-D formulas agree with RσRᵀ to round-off for any angle; the principal stresses from Mohr's circle, from the closed-form cubic and from `eigh` agree to 10⁻⁹ of the stress scale (10⁻⁷ near repeated roots, for a reason section 5 explains); every random plane's traction falls inside the three-circle Mohr region; the rosette reduction reproduces the strain state it was generated from to 10⁻¹⁰.

## 2. Understand and derive

### 2.1 Stress on an inclined plane: the wedge

Take x and y in the plane of the skin, with x along the fuselage axis and y around it. The stress components on the element are σ_x, σ_y and τ_xy, tension positive; τ_xy acts on the face whose outward normal is +x, in the +y direction. Moment equilibrium of the element about its centre requires τ_yx = τ_xy: the tensor is symmetric (Roylance, "Tensor Transformations", Introduction; 16.001 "Stress Equilibrium").

Cut the element by a plane whose outward normal n makes an angle θ with x, anticlockwise. The wedge left behind has three faces: the cut face of length dA (per unit thickness), the −x face of length dA cos θ and the −y face of length dA sin θ.

```
        y (hoop)
        ^               n  (outward normal of the cut face AB, at θ from x)
        |              /
        B             /
        |\           /
 σ_x <--| \      σ_x′ (along n) and τ_x′y′ (along AB) act on AB, length dA
        |  \
 τ_xy ↓ |   \
        O----A------> x (axial)

   OB: the −x face, length dA cos θ, carrying σ_x (toward −x) and τ_xy (toward −y)
   OA: the −y face, length dA sin θ, carrying σ_y (toward −y) and τ_xy (toward −x)
```

The forces on the −x face are −σ_x dA cos θ along x and −τ_xy dA cos θ along y; on the −y face, −τ_xy dA sin θ along x and −σ_y dA sin θ along y. The cut face carries σ_x′ dA along n = (cos θ, sin θ) and τ_x′y′ dA along the face direction (−sin θ, cos θ). Summing forces along n and along the face, and using cos²θ = (1 + cos 2θ)/2, sin²θ = (1 − cos 2θ)/2, 2 sin θ cos θ = sin 2θ:

**σ_x′ = (σ_x + σ_y)/2 + (σ_x − σ_y)/2 · cos 2θ + τ_xy sin 2θ,**
**τ_x′y′ = −(σ_x − σ_y)/2 · sin 2θ + τ_xy cos 2θ,**   (1)

and σ_y′ is σ_x′ at θ + 90°. These are the transformation equations (Roylance, "Direct approach"). They came from equilibrium alone, so they hold for any material, elastic or plastic, and for any stress, however it arose. In matrix form, with R the 2 × 2 rotation whose rows are the new axes,

**σ′ = R σ Rᵀ,**   (2)

which is the transformation law of `P02.14` (Roylance, "General approach"); the first test checks that (1) and (2) agree. Adding σ_x′ and σ_y′ gives σ_x + σ_y at every angle: the trace is invariant.

### 2.2 Mohr's circle

Write C = (σ_x + σ_y)/2 and R = √(((σ_x − σ_y)/2)² + τ_xy²). Squaring and adding the two equations (1) eliminates θ:

**(σ_x′ − C)² + τ_x′y′² = R².**   (3)

Every pair (σ, τ) on any plane through the point lies on one circle of centre (C, 0) and radius R: Mohr's circle (Roylance, "Mohr's circle"). Writing (σ_x − σ_y)/2 = R cos φ₀ and τ_xy = R sin φ₀, (1) becomes σ_x′ − C = R cos(φ₀ − 2θ), τ_x′y′ = R sin(φ₀ − 2θ). So:

- The x face is the point X = (σ_x, τ_xy); the y face is Y = (σ_y, −τ_xy), diametrically opposite.
- Rotating the plane by θ in the body moves the point by **2θ** around the circle. With τ plotted upward as in this lesson's figure, the motion is in the opposite sense (clockwise for an anticlockwise θ). Textbooks differ on this plotting choice; whichever you use, state it, and check one angle with (1).
- The circle crosses the σ axis at the **principal stresses** σ₁,₂ = C ± R, on planes with **no shear**. The direction of σ₁ is at

**θ_p = ½ atan2(2τ_xy, σ_x − σ_y)**   (4)

from x (the two-argument arctangent picks the σ₁ branch; the one-argument tan 2θ_p = 2τ_xy/(σ_x − σ_y) gives either principal direction).
- The top and bottom of the circle are the **maximum in-plane shear** R, at 45° in the body from the principal directions, on planes that also carry the normal stress C.

```
         τ (x face plotted upward)
         ^
         |                 X (σx, τxy)
         |            .  ´ |  ` .
         |         .       |  2θp  .
   ------+-------σ2--------C--------σ1-----> σ
         |         .       |       .
         |            .  _ | _  .
         |   Y (σy, −τxy)
```

**Worked example: the side panel.** σ_x = 37.5 MPa (axial), σ_y = 75 MPa (hoop), τ_xy = 30 MPa (illustrative, from fuselage shear and torsion). C = 56.25 MPa, R = √(18.75² + 30²) = 35.38 MPa. σ₁ = **91.63 MPa**, σ₂ = **20.87 MPa**, θ_p = ½ atan2(60, −37.5) = **61.00°** from the fuselage axis: the largest tension is closer to hoop than to axial, tipped by the shear. A fatigue crack from a rivet hole here would run at 61° − 90° = −29° to the axis. The maximum in-plane shear is **35.38 MPa** on planes at θ_s = θ_p − 45° = 16.00°, with 56.25 MPa normal stress on them. On a splice line whose normal is at 30°, (1) gives σ = 72.86 MPa and τ = 31.24 MPa; at 45°, 86.25 MPa and 18.75 MPa.

### 2.3 Plane stress has three principal stresses

A thin skin has σ_z ≈ 0 (the cabin air pushes on the inner surface with Δp = 0.06 MPa, 0.16 % of the smaller membrane stress). That zero is a principal stress too: the z face carries no shear. The three principal stresses of the side panel are therefore 91.63, 20.87 and 0 MPa, and the largest shear at the point is

**τ_abs = (σ_max − σ_min)/2,**   (5)

the radius of the largest of three Mohr circles (one through σ₁ and σ₂, one through σ₂ and 0, one through σ₁ and 0). For the side panel, τ_abs = 91.63/2 = **45.81 MPa**, on planes at 45° between the in-plane σ₁ direction and the skin normal: planes that cut *through the thickness* at a slant. The in-plane R, 35.38 MPa, is smaller. Whenever σ₁ and σ₂ have the same sign, as in nearly every pressurised skin, the in-plane circle is not the largest and (5) must use the zero.

The fuselage crown shows how large the difference can be. With a fuselage-bending stress of +40 MPa (illustrative) and 8 MPa of torsional shear, σ_x = 77.5, σ_y = 75, τ_xy = 8 MPa: σ₁ = 84.35, σ₂ = 68.15 MPa at 40.56°, an in-plane maximum shear of only **8.10 MPa**, but an absolute maximum of **42.17 MPa**, 5.2 times larger. A yield check with the in-plane value would be wrong by that factor. Taking the inner-surface radial stress −Δp instead of zero changes the absolute maximum to 42.20 MPa.

### 2.4 Three dimensions: the eigenproblem

In 3-D the stress tensor σ is a symmetric 3 × 3 matrix, and on new axes whose unit vectors are the rows of R it becomes RσRᵀ. The traction on a plane of unit normal n is **t = σn** (Cauchy's formula), with normal component σ_n = n·t and shear magnitude τ_n = √(|t|² − σ_n²). A **principal plane** carries no shear, so t is parallel to n:

**σ n = s n,  det(σ − sI) = 0,  s³ − I₁s² + I₂s − I₃ = 0,**   (6)

with I₁ = σ_x + σ_y + σ_z, I₂ = σ_xσ_y + σ_yσ_z + σ_zσ_x − τ_xy² − τ_yz² − τ_xz², I₃ = det σ (Roylance, "Principal stresses and planes in three dimensions"). Because the principal stresses do not depend on the axes, neither do I₁, I₂ and I₃: they are the **invariants**. The principal stresses are the eigenvalues of σ and the principal directions its eigenvectors, which for a symmetric matrix are real and mutually perpendicular (`P01.16`).

**The cubic in closed form.** Remove the mean stress p = I₁/3: the deviator s = σ − pI has invariants J₂ = ½ s:s and J₃ = det s, and its eigenvalues satisfy x³ − J₂x − J₃ = 0. Substituting x = r cos α with r = 2√(J₂/3) and using cos 3α = 4cos³α − 3cos α gives

**cos 3α = (J₃/2)(3/J₂)^{3/2},  σ_k = p + r cos(α − 2πk/3), k = 0, 1, 2,**   (7)

with α = ⅓ arccos(…). All three roots are real because a symmetric matrix has real eigenvalues; that is why the argument of the arccos is always in [−1, 1] (round-off can push it a hair outside, so the code clips it). This needs no iteration and no library. The independent check is `numpy.linalg.eigh`, a numerical symmetric eigensolver built on entirely different arithmetic.

**Principal directions.** Each row of σ − σ_k I is perpendicular to n_k (because (σ − σ_k I)n_k = 0), so the cross product of two rows is parallel to n_k. Taking the largest of the three cross products avoids a pair of nearly parallel rows. For a repeated root every vector in a plane is principal, the rows are all parallel, and the cross products vanish: there is then no unique direction to find.

**Three Mohr circles.** Draw circles on [σ₂, σ₁], [σ₃, σ₂] and [σ₃, σ₁]. Every plane through the point gives a (σ_n, τ_n) inside the big circle and outside the two small ones. The absolute maximum shear (σ₁ − σ₃)/2 sits on planes whose normal bisects n₁ and n₃, with normal stress (σ₁ + σ₃)/2.

**Worked example: a fitting (synthetic stress state).** σ_x = 180, σ_y = −40, σ_z = 90, τ_xy = 60, τ_yz = 30, τ_xz = −25 MPa, an illustrative state of the size a finite-element model reports near a lug, not a measured one. I₁ = 230 MPa, I₂ = 275 MPa², I₃ = −1.199 × 10⁶ MPa³. (7) gives σ₁ = **198.040191**, σ₂ = **95.413423**, σ₃ = **−63.453615 MPa**; `eigh` agrees to 7.5 × 10⁻⁸ Pa. The σ₁ direction is (0.9617, 0.2221, −0.1608), within 16° of x. τ_abs = **130.75 MPa**. Rotating the axes by 30° about z changes every component but leaves I₁ = 230 MPa and I₂ = 275 MPa².

### 2.5 From strain gauges to stress

Strain is also a symmetric tensor, with the tensor shear ε_xy = γ_xy/2, so it obeys (1) with γ/2 in place of τ. A gauge at angle θ reads

**ε_θ = ε_x cos²θ + ε_y sin²θ + γ_xy sin θ cos θ.**   (8)

Three gauges at known angles (a **rosette**) give three equations for ε_x, ε_y, γ_xy. For a rectangular rosette at 0°, 45°, 90°: ε_x = ε₀, ε_y = ε₉₀, γ_xy = 2ε₄₅ − ε₀ − ε₉₀. Plane-stress Hooke's law (`P05.06`, with σ_z = 0) then gives the stresses:

**σ_x = E/(1 − ν²)(ε_x + νε_y),  σ_y = E/(1 − ν²)(ε_y + νε_x),  τ_xy = Gγ_xy,  G = E/(2(1 + ν)).**   (9)

For an isotropic material the principal strain and principal stress directions coincide (a test checks it), so the angle can be read from either. **Worked example (synthetic readings).** A 0/45/90 rosette on a 2024-T3 panel (typical E = 73.1 GPa, ν = 0.33, from the publicly viewable ASM Aerospace Specification Metals data sheet; not design allowables) reads 1200, 900 and −300 µε. Then ε_x = 1200, ε_y = −300, γ_xy = 900 µε; σ_x = 90.3, σ_y = 7.9, τ_xy = 24.7 MPa; σ₁ = **97.2 MPa**, σ₂ = **1.0 MPa** at **15.5°** to the 0° gauge. The readings are invented for the example; they are not a test record.

### 2.6 Where the model stops

- The transformation is exact at a point. A real gauge averages strain over its grid (often 3–6 mm), so near a hole or a fillet it under-reads the peak.
- Plane stress fails near free edges, at joints and in thick fittings, where σ_z and the transverse shears are not small; then the full 3-D tensor of 2.4 is needed.
- The rosette reduction (9) assumes linear isotropic elasticity. Past yield, or in a composite, the strains still transform by (8) but the stresses no longer follow (9).
- The principal stresses describe the load; whether the material survives it is the failure theory's question (`P05.14`).

## 3. Build the computational model

Open `starter/main.py`. `stress_tensor`, `rotation_z` and `main` are given.

### Step 1: the 2-D tools

`transform_stress_2d(sx, sy, txy, theta)` is (1) and must accept an array of angles. `mohr_circle_2d` returns (C, R). `principal_stresses_2d` returns σ₁ ≥ σ₂ and θ_p by (4). `max_inplane_shear` returns R and θ_p − 45°. `plane_stress_abs_max_shear` is (5) with the zero included.

### Step 2: the tensor tools

`transform_tensor(sigma, R)` is RσRᵀ. `invariants` returns (I₁, I₂, I₃). `fuselage_skin_stress(dp, r, t, sigma_bending, tau)` returns (σ_x, σ_y, τ_xy) with x axial and y hoop.

### Step 3: the closed-form cubic

`principal_stresses_3d(sigma)` is (7): compute p, the deviator, J₂ and J₃; if J₂ is zero to round-off the state is hydrostatic and all three roots are p; otherwise clip the arccos argument to [−1, 1] and return the roots sorted from largest to smallest. Do not call `eigh` or `eigvalsh` here: they are the check.

### Step 4: directions, tractions and the Mohr region

`principal_direction(sigma, s)` takes the largest cross product of two rows of σ − sI and normalises it. `absolute_max_shear` is (σ_max − σ_min)/2. `traction_on_plane(sigma, n)` normalises n first and returns (σ_n, τ_n). `mohr_circles_3d` returns the three (centre, radius) pairs; `in_mohr_region(sn, tn, principal)` tests that the point is inside the big circle and outside the two small ones, with a small tolerance scaled by the big radius.

### Step 5: the rosette

`rosette_strains(eps, angles_deg)` builds the 3 × 3 system (8) and solves it. `plane_stress_from_strain` is (9).

### Verify

Run `python3 -m pytest` in the lesson folder. The 15 tests check, in order:

1. (1) equals RσRᵀ for random states and angles, and both 2-D invariants (σ_x′ + σ_y′ and σ_x′σ_y′ − τ_x′y′²) are unchanged;
2. rotated x and y faces both satisfy (3); θ = 180° returns the same stresses and θ = 90° swaps σ_x and σ_y and reverses τ;
3. on the principal planes τ = 0, θ_p points at σ₁, σ₁ and σ₂ are the maximum and minimum of σ_x′ over a fine sweep of angles, and they equal the eigenvalues of the 2 × 2 matrix;
4. the maximum in-plane shear is at θ_p − 45°, with normal stress C;
5. pure shear gives σ₁,₂ = ±τ at 45°; uniaxial tension gives θ_p = 0 and R = σ/2; equal biaxial tension gives a circle of zero radius but an absolute maximum shear of σ/2;
6. the absolute maximum shear exceeds R when σ₁ and σ₂ share a sign and equals R when they do not;
7. the fuselage membrane stresses Δp r/t and Δp r/(2t);
8. the cubic matches `eigh` on 200 random tensors to 10⁻⁹ of the scale;
9. hydrostatic, uniaxial and double-root states, including a randomly rotated double root, to 10⁻⁷ of the scale;
10. invariants survive a random rotation, each eigenvalue satisfies (6), and `rotation_z` reproduces (1) inside the 3-D tensor;
11. principal directions are parallel to `eigh`'s eigenvectors, satisfy σn = sn, and are mutually perpendicular;
12. the σ₁ plane carries σ₁ and no shear; 500 random planes lie in the Mohr region and none exceeds τ_abs; the plane bisecting n₁ and n₃ carries τ_abs with normal stress (σ₁ + σ₃)/2; an impossible traction lies outside;
13. the rosette inverts (8) for rectangular and delta rosettes;
14. Hooke's law (9) inverts the compliance, uniaxial stress gives σ_y = 0, and principal strain and stress directions coincide;
15. `main` writes the figure.

### Compare

Run `python3 main.py`. It prints the side-panel and crown results of 2.2 and 2.3 with an `eigh` check, the fitting of 2.4 with the cubic and `eigh` side by side, and the rosette of 2.5, and writes `mohr_circle.png`: the side panel's circle with the rotated x face marked every 15° (the dots step around the circle by 30°, the angle doubling made visible), and 1500 random planes of the fitting inside its three circles.

## 4. Experiment at home, or investigate real data

This is a computational experiment. All stress states are illustrative, and the rosette readings are synthetic: none of the numbers below is a measurement.

**Experiment 1: when does the largest shear leave the skin?** Sweep the side panel's shear from 0 to 60 MPa (σ_x = 37.5, σ_y = 75 MPa). Record θ_p, σ₁, σ₂, R and τ_abs.

| τ_xy (MPa) | θ_p (°) | σ₁ (MPa) | σ₂ (MPa) | R in-plane (MPa) | τ_abs (MPa) |
|---|---|---|---|---|---|
| 0 | 90.00 | 75.00 | 37.50 | 18.75 | 37.50 |
| 10 | 75.96 | 77.50 | 35.00 | 21.25 | 38.75 |
| 20 | 66.58 | 83.66 | 28.84 | 27.41 | 41.83 |
| 30 | 61.00 | 91.63 | 20.87 | 35.38 | 45.81 |
| 40 | 57.56 | 100.43 | 12.07 | 44.18 | 50.21 |
| 50 | 55.28 | 109.65 | 2.85 | 53.40 | 54.83 |
| 60 | 53.68 | 119.11 | −6.61 | 62.86 | 62.86 |

Until σ₂ changes sign, at τ_xy = √(σ_xσ_y) = 53.03 MPa (where σ₁σ₂ = σ_xσ_y − τ_xy² = 0), the largest shear is out of plane. With pressure alone (τ = 0) it is twice the in-plane value. As the shear grows, θ_p swings from hoop (90°) toward 45°, the pure-shear direction.

**Experiment 2: the closed-form cubic against `eigh`.** On 10 000 random symmetric tensors (seed 2024, components of order 100 MPa), the largest difference between (7) and `eigvalsh`, relative to the largest principal stress, is 6.1 × 10⁻¹⁴ and the median 5.5 × 10⁻¹⁶. Now approach a double root: diag(100, 100(1 + δ), −50) MPa rotated by a fixed random rotation.

| δ | cubic error / 100 MPa | `eigh` error / 100 MPa |
|---|---|---|
| 10⁻¹ | 1.5 × 10⁻¹⁶ | 1.5 × 10⁻¹⁶ |
| 10⁻³ | 8.2 × 10⁻¹³ | 3.0 × 10⁻¹⁶ |
| 10⁻⁵ | 1.1 × 10⁻¹⁰ | 7.5 × 10⁻¹⁶ |
| 10⁻⁷ | 5.3 × 10⁻¹⁰ | 1.5 × 10⁻¹⁶ |
| 10⁻⁹ | 1.7 × 10⁻⁸ | 1.5 × 10⁻¹⁶ |
| 10⁻¹¹ | 1.4 × 10⁻⁸ | 4.5 × 10⁻¹⁶ |
| 0 | 3.0 × 10⁻⁸ | 4.5 × 10⁻¹⁶ |

The cubic's error grows to about √ε_machine = 1.5 × 10⁻⁸ and stops there; `eigh` stays at round-off. Section 5 explains why.

**Experiment 3: random planes in the Mohr region.** For the fitting of 2.4, draw 20 000 random unit normals (seed 13) and compute (σ_n, τ_n) on each. All 100.00 % fall in the three-circle region. The largest τ_n found is 130.74 MPa against τ_abs = 130.75 MPa, and only 1.01 % of planes come within 1 % of it: the absolute maximum lives on a narrow set of planes, so a search over a few hand-picked planes can miss it. The plane bisecting n₁ and n₃ carries τ_n = 130.75 MPa with σ_n = 67.29 MPa = (σ₁ + σ₃)/2. On the octahedral plane (equally inclined to the three principal directions) σ_n = I₁/3 = 76.67 MPa and τ_n = 107.57 MPa, which is √(2J₂)/√3 (the quantity behind von Mises' criterion in `P05.14`).

**Experiment 4: rosette error propagation.** Take the synthetic readings of 2.5 (1200, 900, −300 µε; nominal σ₁ = 97.17 MPa at 15.48°). Add independent Gaussian noise to each reading (seed 42, 20 000 trials; the sampling error of each standard deviation is under 0.5 %):

| gauge noise (µε, 1σ) | std of σ₁ (MPa) | std of θ_p (°) |
|---|---|---|
| 5 | 0.36 | 0.18 |
| 20 | 1.45 | 0.73 |

Then misalign the rosette. If the whole rosette is rotated by δ, the readings are those of a correct rosette on rotated axes: σ₁ is unchanged (97.17 MPa) and θ_p is wrong by exactly δ. If only the 45° grid is off by δ, the reduction mixes the wrong strains:

| 45° gauge misaligned by | σ₁ (MPa) | error | θ_p (°) |
|---|---|---|---|
| 1° | 96.44 | −0.75 % | 14.73 |
| 2° | 95.73 | −1.48 % | 13.93 |
| 5° | 93.76 | −3.51 % | 11.32 |

A rigid misalignment of the whole rosette only costs the angle; a defect inside the rosette costs magnitude. Stacked foil rosettes are made so that their grids are positioned relative to each other to within a fraction of a degree, which is why they are preferred to three single gauges glued separately.

## 5. Validate and improve

**The independent check.** `numpy.linalg.eigh` uses Householder reduction to a tridiagonal matrix and an iterative eigenvalue algorithm, nothing in common with (7). Agreement to 10⁻¹³ on 10 000 random tensors, together with the invariant and eigenvector tests, means the cubic, the invariants and the directions are right. The 2-D formulas are checked against RσRᵀ built independently in the test, and the rosette against the forward strain transformation.

**Why the cubic loses half its digits near a double root.** At a double root cos 3α = ±1, where the slope of arccos is infinite: arccos(1 − η) ≈ √(2η). The argument (J₃/2)(3/J₂)^{3/2} is computed with a round-off error of order ε_machine ≈ 2.2 × 10⁻¹⁶, so α is wrong by about √ε ≈ 1.5 × 10⁻⁸, and the roots by r times that. Experiment 2 shows exactly this. The physics is not ill-conditioned (the eigenvalues of a symmetric matrix are perfectly conditioned, and `eigh` keeps 10⁻¹⁶); the formula is. Two remedies: refine each root with one Newton step on the cubic (extension 1), or use `eigh` in production code and keep (7) as the transparent, library-free check. In engineering terms 1.5 × 10⁻⁸ of 100 MPa is 1.5 Pa, harmless for a stress; the lesson's tolerance of 10⁻⁷ for repeated roots is the honest statement of what the formula delivers.

**Principal directions at repeated roots.** When two principal stresses are equal, any direction in their plane is principal. `principal_direction` then returns a vector built from round-off, and `eigh` returns *some* orthonormal pair. Both are correct; neither is unique. Code that tracks "the σ₁ direction" from step to step of an analysis must expect it to jump when two principal stresses cross.

**Sign and plotting conventions.** The commonest error with Mohr's circle is the sense of rotation. Check one angle with (1) every time: for the side panel, θ = 30° must give σ_x′ = 72.86 MPa, τ_x′y′ = 31.24 MPa.

**Improve.** Add the Newton refinement; handle repeated roots in `principal_direction` by returning a basis of the eigenspace; extend the rosette reduction to any number of gauges by least squares, which gives a residual to detect a faulty gauge.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the wedge derivation of (1); a hand-drawn Mohr's circle of the side panel with σ₁, σ₂, θ_p and the maximum in-plane shear marked and checked against the code; `mohr_circle.png`; the four experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) giving the side panel's principal stresses, their directions, and the absolute maximum shear, and saying which assumption matters most.

### Adversarial self-checks

- **What if the cabin pressure doubles?** Every stress doubles if the shear doubles with it (linearity); the principal directions do not change, and τ_abs becomes 91.6 MPa. If only the pressure doubles and the shear stays at 30 MPa, the state is that of the τ = 15 MPa row scaled by two, and θ_p moves toward hoop.
- **What if τ_xy goes to zero or to infinity?** At zero, the axes are already principal (θ_p = 90° for hoop > axial) and τ_abs = σ_hoop/2. As the shear dominates, σ₁ → τ at 45° and σ₂ → −τ, and the in-plane and absolute maxima coincide.
- **Which assumption fails first in the real world?** Uniform stress over the point: at a rivet hole the stress concentrates by a factor of about 3, and the gauge reads an average. Then plane stress, near a frame attachment.
- **What measurement could falsify the model?** A four-gauge rosette (0/45/90 plus a 135° check gauge): (8) predicts the fourth reading from the other three. A disagreement larger than the gauge noise means a faulty gauge, a strain gradient under the rosette, or a misaligned grid.
- **Does the result survive a change of method?** Mohr's circle, (7) and `eigh` agree to 10⁻⁹ of the scale away from repeated roots and to 10⁻⁷ near them; the rosette inverts its own forward model to 10⁻¹⁰.
- **Which decision changes if the uncertainty doubles?** If the side-panel shear is uncertain by ±10 MPa rather than ±5 MPa, σ₁ spans 83.7–100.4 MPa instead of 87.5–96.0 MPa, and θ_p 57.6°–66.6° instead of 59.1°–63.4°; a fatigue analysis sensitive to the crack direction at a rivet row would need a measurement, which is what the rosette is for.

## Check yourself

1. Why is τ_yx = τ_xy, and which equilibrium condition gives it?
2. On Mohr's circle, how far around the circle does the point move when the plane rotates by 25° in the body?
3. A state has σ_x = 60 MPa, σ_y = 20 MPa, τ_xy = 15 MPa. What are C, R, σ₁ and σ₂?
4. Why is the absolute maximum shear in a pressurised skin larger than the maximum in-plane shear, and on which planes does it act?
5. Name the three invariants of the stress tensor and say why they cannot depend on the axes.
6. In pure shear τ, what are the principal stresses and their directions?
7. Why does the closed-form cubic lose accuracy near a double root while `eigh` does not?
8. A rectangular rosette reads ε₀ = 500, ε₄₅ = 400, ε₉₀ = 100 µε. What is γ_xy?
