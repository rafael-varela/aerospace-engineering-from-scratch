# P05.14  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31). Material properties are the typical values of the lesson's table (2.7), not design allowables.

## 1. (Warm-up: a wing-skin panel)

A 7075-T6 lower wing skin panel at limit load carries σ_x = (120 + 5s) MPa spanwise, σ_y = −40 MPa chordwise and τ_xy = 60 MPa (plane stress).

(a) Find the in-plane principal stresses by Mohr's circle and the three principal stresses including the out-of-plane one. Check with `principal_stresses`.
(b) Find the absolute maximum shear stress, the Tresca stress and the von Mises stress.
(c) Report MS_yield and MS_ultimate by both criteria with the aerospace factors.

Answer form: three principal stresses, three stresses, four margins.

## 2. (Derivation: from energy to von Mises, and the two bounds)

(a) Starting from u = ½ σ : ε and isotropic Hooke's law, show that the volumetric part of the energy is (1 − 2ν)(tr σ)²/(6E) and that the remainder is u_d = ((1 + ν)/(6E))[(σ₁ − σ₂)² + (σ₂ − σ₃)² + (σ₃ − σ₁)²]. Conclude that equating u_d to its value in a tensile test at yield gives σ_vM = F_ty.
(b) Prove 1 ≤ σ_T/σ_vM ≤ 2/√3 for every stress state, and name the states that reach each bound.
(c) Derive the Coulomb–Mohr failure shear stress τ_f = S_tS_c/(S_t + S_c) from (7), and show that the criterion becomes Tresca when S_t = S_c.

Answer form: three derivations.

## 3. (Unfamiliar scenario: a titanium engine-mount lug)

A finite-element model of a Ti-6Al-4V engine-mount lug reports, at its worst element at limit load, the full stress tensor (MPa, in the element's x, y, z axes)

```
σ = [[180,  60,   0],
     [ 60, -40,  30],
     [  0,  30,  90]]
```

(a) Find the three principal stresses, the von Mises and Tresca stresses, and both margins at yield and ultimate.
(b) The lug sits in a pressurised bay; adding a cabin-pressure-like hydrostatic compression of 100 MPa to every normal stress (a deliberately large number) changes which numbers above, and which not? Show it with the code.
(c) Compute the total, volumetric and distortion strain energy densities. What fraction of the stored energy drives yield?

Answer form: numbers with units; one sentence for (b) backed by numbers.

## 4. (Code: sizing the gear tube's wall)

The landing-gear loads of the lesson are revised upward by a factor (1 + s/100): M = 1.2(1 + s/100) kN·m, T = 0.8(1 + s/100) kN·m, P = −20(1 + s/100) kN. The outside diameter stays 50 mm.

(a) With `critical_point` and `scipy.optimize.brentq`, find the smallest wall thickness t that gives MS_yield = +0.15 by von Mises, and the same by Tresca.
(b) Check MS_ultimate at that thickness. Which governs?
(c) Report the mass per metre of each design (density 7850 kg/m³) and the mass penalty of designing to Tresca.

Answer form: two thicknesses to 0.01 mm, one margin, two masses per metre.

## 5. (A brittle casting)

A cast-iron gearbox mounting foot (typical S_t = 214 MPa, S_c = 752 MPa) carries σ_x = −(100 + 2s) MPa and τ_xy = 80 MPa.

(a) Compute the principal stresses and the factor of safety by Coulomb–Mohr and by Rankine.
(b) The load reverses so that σ_x = +(100 + 2s) MPa with the same shear. Repeat.
(c) Explain why the two load directions give such different factors, and why using Tresca with F = S_t would be wrong for this part in both cases.

Answer form: four factors of safety and a short argument.

## 6. (Judgment: a zero margin)

The lesson's experiment 2 shows the gear tube made of 2024-T3 with MS_ultimate = +0.000 by von Mises and −0.030 by Tresca. A colleague writes in a design review: "von Mises is the more accurate criterion for aluminium, so the 2024-T3 tube passes; we can save the steel's weight."

Write a half-page response. Address: what the margin is computed against (typical versus allowable properties); whether the 1.5σ assumption at ultimate is conservative or not for a ductile tube; which features the nominal stress analysis left out (weld or fitting, stress concentration, fatigue); and whether a zero margin can ever be accepted. End with a recommendation and the one piece of evidence you would ask for first.

Answer form: a structured argument of about 250 words.
