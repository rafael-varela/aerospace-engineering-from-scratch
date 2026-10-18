# P10.16  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact (closed form) and which are numerical. Draw the body, the cross-flow and the sign of dS/dx before writing any equation.

Several problems use your personal seed **s**: the day of the month you were born (1–31). ISA sea level: ρ = 1.225 kg/m³, ν = 1.461 × 10⁻⁵ m²/s.

## 1. (Warm-up: an external fuel tank)

A drop tank is a prolate spheroid of length 3.0 m and diameter D = 0.40 + s/100 m, carried at 55 m/s and α = 4°.

(a) Compute its volume (closed form and `body_volume`) and its wetted area.
(b) With slender-body theory, find the normal force and the moment about a point 1.0 m behind its nose. Does the moment depend on that point?
(c) Correct the moment with Lamb's coefficients. By what percentage does slender-body theory over-predict it?

Answer form: two volumes and an area; a force, a moment and a sentence; a moment and a percentage.

## 2. (Derivation: slender-body theory)

(a) Derive dN/dx = 2qα dS/dx from the momentum of the cross-flow carried by each section's apparent mass.
(b) Integrate by parts to obtain N = 2qαS_base and M = 2qα(Vol − S_base(L − x_ref)) about x_ref. Show that for a closed body M is independent of x_ref.
(c) Show that the paraboloid r = R√(x/L) has uniform loading, C_Nα = 2 on the base and its centre of pressure at L/2. Where is the centre of pressure of a cone (r ∝ x)?
(d) Show from (5) of the lesson that k₁ = k₂ = ½ for a sphere (take the limit e → 0) and that (6) reduces to (4) as L/D → ∞.

Answer form: four derivations.

## 3. (Unfamiliar scenario: an airship hull)

An airship envelope is a prolate spheroid 60 m long with fineness ratio 4, flying at 20 m/s at α = 5°.

(a) Compute its volume, Lamb's k₁ and k₂, and the pitching moment by slender-body theory and by the exact solution.
(b) The fins must cancel the moment's derivative dM/dα. With a fin lift slope of 3.0 /rad (on fin area, including interference) at 25 m behind the centre of volume, what fin area is needed for neutral stability? Comment on the size.
(c) Compute the transverse apparent mass exactly and by slender-body theory, and compare with the mass of air the hull displaces. Why does apparent mass matter far more for an airship than for an aeroplane?

Answer form: a volume, two coefficients and two moments; an area and a sentence; three masses and a sentence.

## 4. (Numerical, uses your code: a longer fuselage)

Stretch the trainer's fuselage to L = 7 + s/20 m, keeping the same diameter and the same shape scaled along its length (`trainer_fuselage_radius(x * 7.0 / L)`).

(a) Compute its volume, wetted area and fineness.
(b) Compute Cm_α about the wing quarter chord (x = 2.3 m) with the spheroid-corrected slender body and with the cut loading, and the neutral point for each.
(c) Compute the fuselage's ΔC_D at 55 m/s.
(d) In one paragraph: was the stretch good or bad for stability and drag? What would you do to the tail?

Answer form: three numbers; two Cm_α and two neutral points; a ΔC_D; a paragraph.

## 5. (Numerical, uses your code: how reliable is the optimum fineness?)

(a) Find the fineness of least drag for a spheroid of volume 4.131 m³ at 30, 55 and 90 m/s. Why does speed hardly matter?
(b) Another empirical body form factor in use is FF = 1 + 1.5/f^{1.5} + 7/f³. Replace `body_form_factor` by it (in a probe, not in your module) and recompute the optimum at 55 m/s and D/q at L/D = 3, 6.36 and 10. Compare with the lesson's correlation.
(c) What can a designer conclude from the two correlations, and what can she not?

Answer form: three numbers and a sentence; a number and two rows of three; two sentences.

## 6. (Judgment: which neutral point would you sign off?)

(a) The two loading models of §4 bracket the neutral-point shift between 0.048 c and 0.093 c. Name the physical effect that makes the truth lie inside the bracket, and the one effect left out of both that would move it further forward.
(b) Using Allen–Perkins, by how much does C_N on the maximum section exceed the potential value at 10° and at 20°? At which angle does the viscous term exceed the potential term? Does that matter for the neutral point at cruise?
(c) The real fuselage has a flat-sided cabin and a canopy. Which assumption of slender-body theory do they violate, and in which direction would you expect them to change Cm_α?
(d) Propose one wind-tunnel test that would settle the bracket, and say what you would measure.

Answer form: two effects; two percentages, an angle and a sentence; an assumption and a direction; a test.
