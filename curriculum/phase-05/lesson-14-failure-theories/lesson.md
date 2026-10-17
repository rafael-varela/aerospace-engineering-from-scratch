# P05.14  Failure Theories

> A tensile test gives you one number, and the part carries six stress components. A failure theory is the rule that turns six into one, so the number from the test can decide.

| | |
|---|---|
| **Learning outcomes** | 1. Compute the principal stresses of a 3-D stress state as eigenvalues and of a plane state from Mohr's circle, and turn them into the von Mises and Tresca equivalent stresses. 2. Derive the von Mises criterion from the distortion energy, show that hydrostatic stress does not enter it, and prove that Tresca/von Mises lies between 1 and 2/√3 = 1.155. 3. Report the yield and ultimate margins of safety of an aerospace part at limit load with the factors 1.0 on yield and 1.5 on ultimate, and decide which one governs from the ratio F_tu/F_ty. 4. Apply the Coulomb–Mohr criterion to a brittle material whose compressive strength is several times its tensile strength, and say why a ductile criterion would be wrong there. 5. Implement and test all of this, with an equilibrium check that the section stresses integrate back to the applied loads and a root-finding check that each criterion is homogeneous in the load. |
| **Prerequisites** | `P05.13` (stress transformation, principal stresses, Mohr's circle), `P05.10` (flexure formula), `P05.08` (torsion of a circular tube), `P05.07` (axial stress), `P05.06` (yield strength, ultimate strength, the tensile test), `P01.16` (eigenvalues of a symmetric matrix), `P02.14` (the stress tensor and its rotation), `P03.10` (root finding with `brentq`). Skill check: given σ_x, σ_y and τ_xy you can find the in-plane principal stresses by hand. |
| **Engineering question** | The 4130 steel tube (50 mm outside diameter, 3 mm wall) that carries a light aircraft's main landing-gear leg into the fuselage sees, at the critical limit-load case, a bending moment of 1.2 kN·m, a torque of 0.8 kN·m and an axial compression of 20 kN. Where on the tube is the stress state worst, what are the yield and ultimate margins of safety by von Mises and by Tresca, and how much do the two criteria disagree? The answers: the compression side; von Mises 322 MPa and Tresca 332 MPa (3 % apart); MS_yield = +0.35 and +0.31, MS_ultimate = +0.39 and +0.35; yield governs, because 4130's F_tu/F_ty = 1.54 is above 1.5. |
| **Assumptions and limits** | Linear elastic stress analysis up to the first yield; isotropic, homogeneous material with equal yield strength in tension and compression (true for these metals to a few per cent; not for cast iron, ceramics or composites); stresses proportional to load, so that the ultimate-load stress is 1.5 times the limit-load stress (beyond yield this is only a conservative approximation, since plasticity redistributes stress); static loading (no fatigue, which governs most real landing-gear parts, and no stress concentration at the weld or the bolt holes, which a real check must add); material strengths are typical values, not the A- or B-basis allowables a certified design must use. |
| **Analytical work** | Principal stresses of a plane state from Mohr's circle; von Mises from the deviatoric stress, √(3J₂), and from the distortion energy; Tresca as σ₁ − σ₃ with the out-of-plane zero; the ratio bounds 1 ≤ σ_T/σ_vM ≤ 2/√3; the stresses in a tube under bending, torsion and axial load and the location of the critical point; the margin of safety with aerospace factors; the Coulomb–Mohr envelope and its three closed forms. |
| **Computational deliverable** | `starter/main.py`: `principal_stresses`, `principal_stresses_plane`, `von_mises`, `von_mises_principal`, `tresca`, `max_shear_stress`, `strain_energy_density`, `volumetric_energy_density`, `distortion_energy_density`, `margin_of_safety`, `margins`, `load_factor_to_yield`, `coulomb_mohr_factor`, `rankine_factor`, `tube_section`, `tube_stress_at`, `critical_point`, `tube_resultants`, `sigma_tau_envelope` (`plane_stress_tensor`, the `Material` table and `main` are given). Tests prove: principal stresses equal `numpy.linalg.eigh` and Mohr's circle; three routes to von Mises agree to 10⁻¹⁰; hydrostatic stress gives zero for both criteria; the uniaxial and pure-shear limits; Tresca/von Mises in [1, 2/√3] on 200 random states; invariance under rotation; distortion energy = (1 + ν)σ_vM²/(3E); `brentq` load factor = the homogeneous closed form; the tube stresses integrate back to N, M and T to 10⁻⁴; Coulomb–Mohr reproduces S_t, S_c and S_tS_c/(S_t + S_c) and collapses to Tresca when S_t = S_c. |
| **Practical activity** | **computational experiment**: a sweep of the stress ratio τ/σ showing where Tresca and von Mises part company; the same tube in four aerospace alloys, finding which margin governs; a brittle cast part under five stress states by Coulomb–Mohr and by Rankine; the failure envelopes plotted with the tube's limit and ultimate states on them. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `failure_theories.png`, the derivation of von Mises from the distortion energy and of the ratio bounds, the three experiment tables, and a half-page margin-of-safety note for the gear tube in the format of part 6. |
| **Assessment** | `assessments/phase-05/grading-rubric.md`; competencies 1–4. Pass: von Mises is derived from the distortion energy and Tresca from the largest Mohr circle; all public tests pass; the margins are reported with the factor, the allowable and the criterion named; the note says which margin governs and why. |
| **References** | `@mit-ocw-3-11-mechanics-of-materials` Roylance, module "Yield and Plastic Flow": Introduction; Multiaxial stress states; Effect of hydrostatic pressure (required-free): the two criteria and why pressure does not cause yield. `@nasa-tm-x-73305-astronautic-structures-manual` Vol. I, Section A1 "STRESS AND STRAIN": 1.4.0 Theories of Failure (1.4.1 Elastic Failure; 1.4.2 Interaction Curves); Section A3 "COMBINED STRESSES": 3.2.0 Stress Ratios, Interaction Curves, and Factor of Safety (supplementary-free): how aerospace practice writes margins of safety. `@tudelft-alderliesten-aerospace-structures-materials` Ch. 8 "Considering strength & stiffness": 8.2 Structural performance; 8.3 Selecting the appropriate criterion; Ch. 9 "Design & certification": 9.2 Safety, regulations and specifications (required-free; read online, link only): where the factor 1.5 comes from. |
| **Extensions** | ★ Add a stress-concentration factor K_t at a weld toe and redo the margins: which criterion now governs, and what does that say about where the tube would actually crack? ★ Implement the interaction-curve form R_b² + R_s² = 1 of the Astronautic Structures Manual for bending plus torsion of the tube and compare it with your von Mises envelope. ★ Plot the von Mises cylinder and the Tresca hexagonal prism in principal-stress space (σ₁, σ₂, σ₃) and show that the hydrostatic line is their common axis. |
| **Time** | unmeasured |

## 1. The engineering problem

A landing-gear leg on a light aircraft is a cantilever. At touchdown the wheel is pushed up and back, and the leg carries those loads to the fuselage as a bending moment, a torque (because the wheel axle is offset from the leg's axis) and an axial force. Where the leg enters the fuselage, a 4130 steel tube takes all three at once. Each one alone is a lesson you have done: the flexure formula of `P05.10` gives the bending stress, the torsion formula of `P05.08` the shear stress, `P05.07` the axial stress. What you have not done is decide whether the *combination* is safe. The tensile test of `P05.06` says 4130 yields at about 435 MPa in uniaxial tension. The tube's most loaded point has a normal stress of 289 MPa and a shear stress of 81 MPa acting together. Is that more or less dangerous than 435 MPa of pure tension? Neither number alone tells you.

A **failure theory** answers by mapping the multi-axial stress state to a single **equivalent stress** that can be compared with the uniaxial strength. For ductile metals, which fail by yielding (atomic planes slipping over one another), the two standard theories are **Tresca** (yield when the largest shear stress reaches its value in the tensile test) and **von Mises** (yield when the energy stored in changing the material's *shape* reaches its value in the tensile test). For brittle materials, which crack before they yield and are much stronger in compression than in tension, **Coulomb–Mohr** is the usual choice. Aerospace practice then wraps the equivalent stress in a **margin of safety**: the factor by which the load could grow before the criterion is met, minus one, after applying a factor of safety of 1.0 on yield at the **limit load** (the largest load expected in service) and 1.5 on rupture at the **ultimate load** (1.5 × limit).

Success criteria, written first: the principal stresses agree with `eigh`; three routes to von Mises agree to 10⁻¹⁰; the section stresses integrate back to the applied loads; the tube's margins come out positive with the governing one identified and explained; the brittle check reproduces its closed forms. If the answer were wrong in the unconservative direction, the gear would take a permanent bend on a firm landing; wrong in the conservative direction, the tube would be heavier than it needs to be on an aircraft where every kilogram is fuel or payload.

## 2. Understand and derive

### 2.1 Principal stresses, again

At any point the stress state is a symmetric 3 × 3 tensor σ (`P02.14`, `P05.13`). Its eigenvalues are the **principal stresses**, ordered σ₁ ≥ σ₂ ≥ σ₃, acting on three mutually perpendicular planes that carry no shear. For a **plane stress** state (σ_x, σ_y, τ_xy; everything with a z index zero), one principal stress is the zero out of plane and the other two come from Mohr's circle with centre c = (σ_x + σ_y)/2 and radius R = √(((σ_x − σ_y)/2)² + τ_xy²):

**σ_a, σ_b = c ± R.**   (1)

The critical point on the gear tube has σ_x = −289.4 MPa (bending plus axial), σ_y = 0 and τ_xy = 81.4 MPa (torsion): c = −144.7, R = 166.0, so σ_a = 21.3 MPa, σ_b = −310.8 MPa. Sorting with the out-of-plane zero gives **σ₁ = 21.3, σ₂ = 0, σ₃ = −310.8 MPa**. A compressive normal stress plus a shear produces a small *tensile* principal stress, a fact that matters for brittle materials (2.5).

### 2.2 Tresca: the largest shear stress

Tension test a ductile bar to its yield stress F_ty. Mohr's circle for uniaxial tension runs from 0 to F_ty, so the largest shear stress at yield is F_ty/2, on planes at 45° to the axis, which is where the slip lines (Lüders bands) of a yielding steel bar appear. Tresca's hypothesis is that yield starts in *any* stress state when the largest shear stress reaches that same value. The largest of the three Mohr circles has diameter σ₁ − σ₃, so the criterion is

**σ_T = σ₁ − σ₃ = F_ty.**   (2)

Two warnings. First, use all three principal stresses: equal biaxial tension σ_x = σ_y = 200 MPa has equal in-plane principal stresses (in-plane Mohr circle of zero radius) but σ₃ = 0 out of plane, so σ_T = 200 MPa, not zero. Second, Tresca's equivalent stress is a function of the *differences* of principal stresses, so adding the same pressure to all three changes nothing: a hydrostatic state, however large, never yields. Bridgman's experiments on metals under pressures of thousands of atmospheres confirmed this (Roylance, "Effect of hydrostatic pressure").

### 2.3 von Mises: the distortion energy

The elastic strain energy density of a stress state is u = ½ σ : ε, with Hooke's law ε = ((1 + ν)σ − ν tr(σ) I)/E. Split the stress into a mean (hydrostatic) part p = tr(σ)/3 times I and a **deviatoric** remainder s = σ − pI. The mean part changes the volume and the deviatoric part changes the shape, and because tr(s) = 0 the cross terms in the energy vanish, so u splits cleanly into a volumetric part

**u_v = (1 − 2ν)(tr σ)² / (6E)**   (3)

and a distortion part u_d = u − u_v. Writing u_d in principal stresses and collecting terms,

**u_d = ((1 + ν)/(6E)) [(σ₁ − σ₂)² + (σ₂ − σ₃)² + (σ₃ − σ₁)²] = ((1 + ν)/(3E)) σ_vM²,**   (4)

with the **von Mises stress**

**σ_vM = √(½[(σ₁ − σ₂)² + (σ₂ − σ₃)² + (σ₃ − σ₁)²]) = √(3J₂),  J₂ = ½ s : s.**   (5)

In uniaxial tension σ_vM = σ, so yield in a tensile test happens at u_d = (1 + ν)F_ty²/(3E), and the von Mises criterion says yield happens in any state when the distortion energy reaches that value: **σ_vM = F_ty**. In components, without finding principal stresses,

σ_vM = √(½[(σ_x − σ_y)² + (σ_y − σ_z)² + (σ_z − σ_x)²] + 3(τ_xy² + τ_yz² + τ_zx²)),

and for the shaft-like state σ_x = σ, τ_xy = τ it collapses to **σ_vM = √(σ² + 3τ²)**, against Tresca's √(σ² + 4τ²). At the gear-tube point: σ_vM = √(289.4² + 3 × 81.4²) = **322.0 MPa**; σ_T = √(289.4² + 4 × 81.4²) = **332.1 MPa**. The distortion energy there is 217.4 kJ/m³, and the test checks that (4) holds on random tensors. For 4130 the yield value of u_d is (1.29)(435 MPa)²/(3 × 205 GPa) = 396.9 kJ/m³.

Like Tresca, von Mises depends only on the deviatoric stress, so hydrostatic pressure does not enter. Unlike Tresca, it is smooth: a circle-like ellipse in the σ₁–σ₂ plane rather than a hexagon, which is why plasticity theory and finite-element codes prefer it.

### 2.4 How far apart can they be?

Both criteria agree in uniaxial tension (σ_T = σ_vM = σ) and in equal biaxial tension. They differ most in **pure shear**, where σ₁ = τ, σ₃ = −τ: σ_T = 2τ, σ_vM = √3 τ. A shaft in pure torsion therefore yields at τ = F_ty/2 by Tresca and at τ = F_ty/√3 = 0.577F_ty by von Mises (217.5 against 251.1 MPa for 4130). In general, with σ₂ between σ₁ and σ₃, write σ₂ = σ₃ + f(σ₁ − σ₃) with 0 ≤ f ≤ 1; then σ_vM = (σ₁ − σ₃)√(1 − f + f²), whose minimum over f is √3/2 at f = ½, so

**1 ≤ σ_T / σ_vM ≤ 2/√3 = 1.1547.**   (6)

Tresca is never less conservative than von Mises, and never more than 15.5 % more so. Torsion tests on thin-walled tubes put real ductile metals closer to von Mises, which is why it is the default in aerospace stress analysis; Tresca survives where a quick, conservative, hand-checkable bound is worth more than 15 %.

### 2.5 Brittle materials: Coulomb–Mohr

Cast iron, ceramics and many castings fail by cracking, in tension at an ultimate strength S_t and in compression at a much larger S_c (gray cast iron, typically S_t ≈ 214 MPa and S_c ≈ 752 MPa). Mohr's idea: draw the tensile-test circle (from 0 to S_t) and the compression-test circle (from −S_c to 0); a state fails when its largest circle touches the straight line tangent to both. Working out where the tangent line sits gives

**σ₁/S_t − σ₃/S_c = 1/n,**   (7)

with n the factor of safety and σ₁ ≥ σ₂ ≥ σ₃ including the out-of-plane zero. Three closed forms follow: uniaxial tension fails at S_t; uniaxial compression at S_c; pure shear at **τ = S_tS_c/(S_t + S_c)**, 166.6 MPa for the cast iron, well below S_t. When both in-plane principal stresses are tensile the line would allow σ₁ > S_t, which no brittle material does, so a **tension cut-off** σ₁ ≤ S_t is added: 1/n = max(σ₁/S_t − σ₃/S_c, σ₁/S_t). With S_t = S_c = S, (7) is σ₁ − σ₃ = S: Tresca. The cruder **Rankine** (maximum normal stress) criterion ignores the interaction between σ₁ and σ₃; in pure shear it gives n = S_t/τ, 28 % more optimistic than Coulomb–Mohr for the cast iron.

### 2.6 Stresses in the gear tube and where they are worst

For the tube (outside diameter d_o = 50 mm, wall t = 3 mm): A = π(d_o² − d_i²)/4 = 443.0 mm², I = π(d_o⁴ − d_i⁴)/64 = 122 812 mm⁴, J = 2I = 245 624 mm⁴, outer radius c = 25 mm. At a point at angle θ around the circumference (θ = 0 on the tension side of the bending) and radius r,

**σ_x = P/A + M r cos θ / I,  τ = T r / J.**   (8)

With M = 1.2 kN·m, T = 0.8 kN·m and P = −20 kN: bending ±244.3 MPa at the outer fibre, axial −45.2 MPa, torsional shear 81.4 MPa. The shear is the same all round the outer surface; the normal stress is 199.1 MPa at θ = 0 and −289.4 MPa at θ = 180°. Both criteria grow with |σ_x|, so the **critical point is the compression side** (σ_vM = 322.0 MPa there against 244.0 MPa on the tension side, and 288.2 MPa on the inner surface at θ = 180°). If the axial load reversed to tension, the critical point would move to θ = 0. An equilibrium check that every analysis of this kind should carry: integrate σ_x and τ over the annulus and you must get back N = −20 kN, M = 1.2 kN·m, T = 0.8 kN·m. `tube_resultants` does it with a midpoint rule and recovers all three to four decimal places.

### 2.7 Margins of safety, aerospace style

The airworthiness rules for transport and light aircraft (summarised in Alderliesten 9.2) ask two things of a structure: at **limit load** no detrimental permanent deformation, and at **ultimate load = 1.5 × limit** no rupture for at least three seconds. Written as margins with the equivalent stress at limit load σ_eq,

**MS_yield = F_ty / (1.0 × σ_eq) − 1,  MS_ult = F_tu / (1.5 × σ_eq) − 1,**   (9)

and the part passes when both are ≥ 0 (Astronautic Structures Manual A3.2.0 uses the same form). A margin of +0.35 means the load could rise 35 % before the criterion is met. Treating the ultimate-load stress as 1.5σ_eq assumes the stress is still proportional to load, which is conservative for a ductile part that has yielded locally. Which margin governs depends only on the material: MS_ult < MS_yield exactly when **F_tu/F_ty < 1.5**. Typical room-temperature values used in this lesson (not design allowables; certified designs use the lower A- or B-basis values of MMPDS, and Roylance's "Material Properties" appendix tabulates handbook values of the same kind):

| material | E [GPa] | ν | F_ty [MPa] | F_tu [MPa] | F_tu/F_ty |
|---|---|---|---|---|---|
| 2024-T3 | 73.1 | 0.33 | 345 | 483 | 1.40 |
| 7075-T6 | 71.7 | 0.33 | 503 | 572 | 1.14 |
| Ti-6Al-4V (annealed) | 113.8 | 0.342 | 880 | 950 | 1.08 |
| 4130 (normalised) | 205 | 0.29 | 435 | 670 | 1.54 |

For the 4130 tube, **MS_yield = 435/322.0 − 1 = +0.351** (von Mises) and **+0.310** (Tresca); **MS_ult = 670/(1.5 × 322.0) − 1 = +0.387** and **+0.345**. Yield governs, by a little, because 1.54 > 1.5. The load factor to first yield, found by `brentq` on the criterion, is 1.351 (von Mises) and 1.310 (Tresca): exactly 1 + MS_yield, because both criteria are homogeneous of degree one in the stress.

### 2.8 Where the models stop

- Yield criteria describe the *onset* of yielding of an isotropic metal. They do not predict fatigue cracking under repeated loads, which is what usually limits landing gear, nor creep, nor buckling (`P05.15`), which can occur at stresses far below yield.
- Stress concentrations at welds, holes and fillets multiply the local stress by 2–3; the nominal stresses used here are the input to that check, not a replacement for it.
- Rolled plate and extrusions are mildly anisotropic, and some alloys yield at different stresses in tension and compression (F_cy ≠ F_ty); then the equivalent stress needs the matching allowable.
- Composites need their own criteria (maximum stress, Tsai–Wu), because their strength depends on direction.
- Coulomb–Mohr is fitted to two tests; between them it is an interpolation, and under triaxial compression it predicts strengthening that only some brittle materials show.

## 3. Build the computational model

Open `starter/main.py`. `plane_stress_tensor`, the `Material` table and `main` are given.

### Step 1: principal stresses two ways

`principal_stresses(sigma)` is `np.linalg.eigvalsh` of the symmetrised tensor, sorted descending. `principal_stresses_plane(sx, sy, txy)` is (1). They are independent routes to the same numbers, and the first test compares them.

### Step 2: the criteria

`von_mises(sigma)` from components (no eigenvalues); `von_mises_principal(p)` from (5); `tresca(sigma)` as σ₁ − σ₃ from `principal_stresses`, so the out-of-plane zero of a plane state is included automatically; `max_shear_stress` is half of it.

### Step 3: the energy route

`strain_energy_density` is ½ σ : ε with isotropic Hooke's law; `volumetric_energy_density` is (3); `distortion_energy_density` is their difference. The test checks (4), which is the von Mises criterion arrived at from energy, not from a formula.

### Step 4: margins and load factors

`margin_of_safety(allowable, applied, factor)` is (9); `margins(sigma, mat)` returns both criteria and both margins as a dictionary. `load_factor_to_yield` brackets and solves λ with σ_eq(λσ) = F_y by `brentq`:

```python
g = lambda lam: f(lam * sigma_unit) - F_y
hi = 1.0
while g(hi) < 0.0:
    hi *= 2.0
return brentq(g, 0.0, hi, xtol=1e-14, rtol=1e-14)
```

It should equal F_y/σ_eq(σ); if it does not, your criterion is not homogeneous, which means it is wrong.

### Step 5: brittle criteria

`coulomb_mohr_factor(sigma, S_t, S_c)` is (7) with the tension cut-off; `rankine_factor` the maximum normal stress.

### Step 6: the tube

`tube_section` gives A, I, J, c; `tube_stress_at(theta, r, M, T, P, d_o, t)` is (8) as a tensor in local (axial, hoop, radial) coordinates; `critical_point` scans the outer surface; `tube_resultants` integrates σ_x, σ_x r cos θ and τ r over the annulus on a midpoint grid. `sigma_tau_envelope` traces each criterion's boundary in the (σ, τ) plane by scaling unit states on rays.

Running `main()` prints

```
4130 tube 50 x 3 mm: A = 443.0 mm^2, I = 122812 mm^4, J = 245624 mm^4
limit load M = 1.2 kN m, T = 0.8 kN m, P = -20 kN: bending 244.3 MPa, axial -45.2 MPa, torsional shear 81.4 MPa
critical point at theta = 180 deg: sigma = -289.4 MPa, tau = 81.4 MPa
tension side (theta = 0): sigma = 199.1 MPa, von Mises 244.0 MPa
principal stresses at the critical point: 21.3, 0.0, -310.8 MPa
von Mises 322.0 MPa, Tresca 332.1 MPa (ratio 1.031)
  MS yield (FS 1.0, F_ty 435):   von Mises +0.351   Tresca +0.310
  MS ultimate (FS 1.5, F_tu 670): von Mises +0.387   Tresca +0.345
resultants from the integrated stresses: N = -20.0000 kN, M = 1.2000 kN m, T = 0.8000 kN m
distortion energy 217.43 kJ/m^3 = (1+nu) sigma_vm^2 / 3E = 217.43 kJ/m^3
```

followed by the three experiment tables of part 4, and writes `failure_theories.png` (the von Mises, Tresca and Coulomb–Mohr envelopes in the σ–τ plane with the tube's limit and ultimate states marked; the equivalent stress around the tube's circumference).

### Verify

Run `python -m pytest curriculum/phase-05/lesson-14-failure-theories`. The tests:

1. **Principal stresses**: equal to `eigh`, sorted, trace preserved; Mohr's circle (1) against `eigvalsh` with both invariants.
2. **Three routes to von Mises**: components, principal stresses and √(3J₂) agree to 10⁻¹⁰.
3. **Hydrostatic**: both criteria zero for pI and unchanged when pI is added to any state.
4. **Limits**: uniaxial σ gives σ for both; pure shear τ gives √3τ and 2τ; equal biaxial tension gives σ for Tresca because of the out-of-plane zero.
5. **Ratio bounds** (6) on 200 random tensors; invariance under random rotations QσQᵀ.
6. **Distortion energy** (4) on random states; zero for hydrostatic stress.
7. **Homogeneity**: `brentq` load factor equals F_y/σ_eq; σ_eq(3σ) = 3σ_eq(σ), σ_eq(−σ) = σ_eq(σ).
8. **Margins**: the definition (9), the sign of the stress irrelevant, Tresca the conservative one in shear; the material table is sane.
9. **Tube**: A, J = 2I, the thin-wall limit I → πr³t; the integrated stresses return N, M, T to 10⁻⁴ (equilibrium); the critical point moves with the sign of the axial load.
10. **Brittle**: Coulomb–Mohr's three closed forms, the tension cut-off, never more optimistic than Rankine in plane stress, Tresca when S_t = S_c; the envelopes lie on their criteria.
11. **main** writes the figure.

### Compare

`numpy.linalg.eigh` is the established tool for the principal stresses and the oracle for Mohr's circle; agreement is to round-off. For the margins themselves there is no library to compare with, so the independent checks are internal: the energy route to von Mises, the root-finding route to the load factor, and the equilibrium of the section stresses. A finite-element code would report the same von Mises value at the tube surface (it is the default contour plot of every structural solver) and would differ only where the beam formulas (8) stop being exact, near the fuselage fitting.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction first.

1. **Where the criteria disagree.** For σ_x = 1, τ_xy = r with r from 0 to ∞, tabulate σ_vM, σ_T and their ratio. Predict the ratio at r = 0 and as r → ∞ from part 2.4, and guess where it is halfway.
2. **Same tube, four alloys.** Keep the loads and the tube, change the material through the table of 2.7. Predict, from F_tu/F_ty alone, which alloys will be ultimate-critical. Predict whether 2024-T3 passes.
3. **A brittle part.** For a cast part with S_t = 214 MPa and S_c = 752 MPa (typical gray cast iron), compute the factor of safety by Coulomb–Mohr and by Rankine for uniaxial tension, uniaxial compression, pure shear, equal biaxial tension and equal tension–compression, each at a 100 MPa scale. Predict which states the two criteria agree on.
4. **Envelopes.** Plot the three envelopes in the σ–τ plane and mark the tube at limit and at ultimate load. Predict whether the ultimate-load point lies inside the 4130 *yield* envelope (it need not: the ultimate check uses F_tu).

Raw data format: CSV tables `tau_over_sigma, vm_over_sigma, tresca_over_sigma, ratio` (1); `material, F_ty, F_tu, MS_y_vm, MS_u_vm, MS_y_tr, MS_u_tr, governs` (2); `state, n_coulomb_mohr, n_rankine` (3); plus `failure_theories.png`.

What you should find:

```
tau/sigma   von Mises/sigma   Tresca/sigma   ratio
     0.00            1.0000         1.0000   1.0000
     0.25            1.0897         1.1180   1.0260
     0.50            1.3229         1.4142   1.0690
     1.00            2.0000         2.2361   1.1180
     2.00            3.6056         4.1231   1.1435
      inf            1.7321         2.0000   1.1547
```

The tube's state (τ/σ = −0.28) sits near the 0.25 row: the 3 % disagreement is small because bending dominates. A shaft in nearly pure torsion would see the full 15.5 %.

```
material     F_ty  F_tu   MS_y(vM)  MS_u(vM)  MS_y(Tr)  MS_u(Tr)  governs
2024-T3       345   483    +0.072    +0.000    +0.039    -0.030  ultimate
7075-T6       503   572    +0.562    +0.184    +0.515    +0.148  ultimate
Ti-6Al-4V     880   950    +1.733    +0.967    +1.650    +0.907  ultimate
4130          435   670    +0.351    +0.387    +0.310    +0.345  yield
```

Only 4130, with F_tu/F_ty = 1.54, is yield-critical. 2024-T3 passes by von Mises with an ultimate margin of +0.0001 (zero to three places) and *fails* by Tresca (−0.030): the choice of criterion decides the answer, which is exactly when you must say which one you used and why. The high-strength alloys gain margin from their strength but lose relative margin to the 1.5 factor, because their ultimate strength is barely above their yield.

```
  uniaxial tension       Coulomb-Mohr n = 2.140   Rankine n = 2.140
  uniaxial compression   Coulomb-Mohr n = 7.520   Rankine n = 7.520
  pure shear             Coulomb-Mohr n = 1.666   Rankine n = 2.140
  biaxial tension        Coulomb-Mohr n = 2.140   Rankine n = 2.140
  tension + compression  Coulomb-Mohr n = 1.666   Rankine n = 2.140
```

They agree whenever only one principal stress is non-zero or both are tensile, and part company in the second and fourth quadrants, where a compressive principal stress helps open a crack in the perpendicular direction.

## 5. Validate and improve

| quantity | theory (§2) | code | difference | cause |
|---|---|---|---|---|
| principal stresses, plane | Mohr's circle (1) | `eigvalsh` | < 10⁻¹² relative | — |
| von Mises, three routes | (5), components, √(3J₂) | same | < 10⁻¹⁰ | round-off |
| distortion energy | (1 + ν)σ_vM²/3E = 217.43 kJ/m³ | 217.43 kJ/m³ | 0 | — |
| Tresca/von Mises | 1 to 1.1547 | 1.031 at the tube; bounds held on 200 random states | — | — |
| pure shear yield, 4130 | F_ty/√3 = 251.1, F_ty/2 = 217.5 MPa | same | 0 | — |
| tube resultants | N = −20, M = 1.2, T = 0.8 | −20.0000, 1.2000, 0.8000 | < 10⁻⁴ | midpoint rule, 40 × 360 cells |
| load factor to yield | 1 + MS_y = 1.351 (vM), 1.310 (Tr) | `brentq` 1.351, 1.310 | < 10⁻¹² | — |
| Coulomb–Mohr, pure shear | S_tS_c/(S_t + S_c) = 166.6 MPa | n = 1.666 at 100 MPa | 0 | — |

Every row with a closed form is met to round-off or to the integration error, which is the check that the code computes what the theory says. What the code cannot check is the theory: whether the tube really yields at σ_vM = 435 MPa depends on the material, its heat treatment, the weld next to the critical point, and whether the stresses near the fuselage fitting follow the beam formulas. The engineering answer is **MS_yield = +0.35, MS_ultimate = +0.39 by von Mises (+0.31 and +0.35 by Tresca), yield governing, on typical properties**. With a B-basis allowable perhaps 10–15 % below the typical F_ty, the yield margin would shrink to about +0.15 to +0.22: still positive, and a reason to state which allowable a margin uses every time you write one.

What to improve: a stress-concentration factor at the weld (Extension 1), the real allowables, and a fatigue check under the landing spectrum. What not to improve: the integration grid of `tube_resultants` or the scan resolution of `critical_point`, which are already far finer than the load uncertainty.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivation of (4)–(5) from the energy split and of the bounds (6); `failure_theories.png`; the three tables of part 4; and a half-page margin note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) for the gear tube: load case, critical point, stress components, equivalent stress by both criteria, allowables with their source, factors of safety, margins, the governing case and the main limitation.

### Adversarial self-checks

- **What if the torque doubles?** τ goes to 162.8 MPa and σ_vM at the critical point to 404.1 MPa; MS_yield falls from +0.35 to +0.08. The ratio Tresca/von Mises grows toward 1.155 because the state moves toward pure shear: the more torsion, the more the choice of criterion matters.
- **What if the axial load goes to zero?** The tension and compression sides become equally loaded at σ = ±244.3 MPa and σ_vM = 282.1 MPa. The compressive axial load *cost* 40 MPa of equivalent stress at the critical point; a tensile one would have moved the critical point to the other side.
- **Which assumption fails first in the real world?** The nominal stress: the tube is welded to a fitting, and the weld toe has a stress concentration and a heat-affected zone with lower strength. Next, the static treatment of a landing load, which is dynamic and repeated thousands of times.
- **What measurement could falsify the model?** Strain-gauge rosettes on the tube's compression side under a known static load: they give ε in three directions, hence σ_x and τ, which can be compared with (8) directly; and a proof-load test to limit load followed by a measurement of permanent set.
- **Does the result survive a change of method?** The margin by Tresca and by von Mises differs by 0.04 here; principal stresses by Mohr and by `eigh` agree; the load factor by root finding agrees with the closed form. The result is robust to the method and sensitive only to the inputs: the loads, the wall thickness (a 2.5 mm wall would bring MS_yield down to +0.16) and the allowable.
- **Which decision changes if the allowable were 15 % lower?** For 4130 nothing changes sign (MS_yield = +0.15). For 2024-T3 the already-zero ultimate margin becomes clearly negative, and the alloy is ruled out for this tube.

## Check yourself

1. Why does a hydrostatic stress state not cause yield by either Tresca or von Mises?
2. A thin-walled shaft is in pure torsion at τ = 200 MPa. What are its Tresca and von Mises equivalent stresses?
3. A plane stress state has σ_x = σ_y = 150 MPa and no shear. What is the Tresca equivalent stress, and why is it not zero?
4. Write the margin of safety at ultimate load for an equivalent limit stress of 300 MPa and F_tu = 480 MPa.
5. Which margin governs, yield or ultimate, for a material with F_ty = 400 MPa and F_tu = 560 MPa, and how can you tell without computing either?
6. What is the largest possible ratio of the Tresca to the von Mises stress, and in what stress state is it reached?
7. For a brittle material with S_t = 200 MPa and S_c = 600 MPa, what shear stress causes failure by Coulomb–Mohr?
8. Name one failure mode of a landing-gear tube that no yield criterion predicts.
