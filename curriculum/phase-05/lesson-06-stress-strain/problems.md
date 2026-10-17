# P05.06  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical. Material values are the typical values of `MATERIALS` in `starter/main.py`; every tensile record in this set is **synthetic** and must be labelled so in anything you submit.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: a control-rod end in tension)

A 7075-T6 rod of diameter d = 6 + s/10 mm and length 300 mm carries an axial tension of 10 kN.

(a) Find the engineering stress, the strain, the elongation and the change of diameter.
(b) Find the true stress and true strain, and the ratio of true to engineering stress. Is the difference worth keeping?
(c) What tension would bring the rod to its typical yield strength?

Answer form: four numbers with units, two numbers and a ratio, one force.

## 2. (Derivation: true measures, Considère and the elastic constants)

(a) Starting from constant volume, derive σ_t = σ(1 + ε) and ε_t = ln(1 + ε). State the assumption that fails after the maximum load.
(b) Derive Considère's criterion dσ_t/dε_t = σ_t from dF = 0. For a power-law material σ_t = K ε_tᵐ, show that the uniform true strain is ε_t = m. The lesson's synthetic 2024-T3 material is Ramberg–Osgood with n = 9; compare its uniform true strain from `considere_point` with 1/n and explain the difference.
(c) Show that the 0.2 % offset yield of the Ramberg–Osgood curve (6) is σ_y, and find the 0.05 % offset (proof) stress in terms of σ_y and n.
(d) A square element in pure shear τ has principal stresses +τ and −τ on planes at 45°. Using only the uniaxial Hooke's law and Poisson's ratio for those principal stresses, compute the change of length of the element's diagonals and show that τ = Gγ with G = E/2(1 + ν).

Answer form: four derivations, one comparison of two numbers.

## 3. (Unfamiliar scenario: a titanium pressurant sphere)

A spherical Ti-6Al-4V helium tank has mid-surface radius R = 300 mm and wall t = 1.2 mm, and is pressurised to p = 1.5 + s/30 MPa.

(a) By a free-body diagram of a hemisphere, show that the wall stress is σ = pR/2t in every tangential direction, and evaluate it.
(b) With `hooke_strain`, find the tangential strain, the thickness strain, the growth of the radius, the change of wall thickness and the increase of internal volume (ΔV/V ≈ 3ε for small strains).
(c) What would the radius growth be if you ignored the Poisson effect of the second tangential stress? By what percentage would that overestimate it?
(d) What fraction of the typical yield strength is the wall stress? Name one reason a real tank is not designed to this ratio.

Answer form: a derivation and a stress, five numbers with units, a number and a percentage, a ratio and a sentence.

## 4. (Numerical, uses your code: reduce a titanium coupon)

Generate a **synthetic** Ti-6Al-4V record with `synthetic_tensile_test(E=113.8e3, sigma_y=880.0, n=20.0, eps_fracture=0.14, drop=0.10, noise_stress=2.0, noise_strain=2e-5, seed=s)`.

(a) Reduce it with `analyze_tensile(..., window=21)` and tabulate E, Fty, Ftu, ε_u, elongation, toughness and resilience.
(b) Compare each property with the generator's inputs and with `considere_point(113.8e3, 880.0, 20.0)`. Explain every difference larger than the noise.
(c) Repeat (a) with the modulus windows 5–20 %, 10–40 % and 20–60 %. Which property is most sensitive, and why is this titanium coupon less sensitive to the window than the 2024-T3 one of the lesson?
(d) Reduce 100 seeds and report the mean and standard deviation of E and Fty. Is your own seed's result within 2σ of the mean?

Answer form: a 7-row table, a comparison table, a 3-row table and a sentence, four numbers and a yes/no.

## 5. (Numerical and derivation: diagnose a crosshead-strain record)

A laboratory sends the **synthetic** 2024-T3 record `synthetic_tensile_test(**COUPON, compliance=c)` with c = (3 + s/5) × 10⁻⁶ per MPa, but does not tell you c.

(a) Reduce it with `analyze_tensile`. What modulus do you get, and why should it alarm you?
(b) Assuming the true modulus is the clean-record value 72.49 GPa, estimate c from (7) and compare with the value used.
(c) Show algebraically that the offset yield with the fitted modulus does not depend on c, and confirm it numerically.
(d) Compute the offset yield the laboratory would report if it used the handbook E = 73.1 GPa with this record. By how much is it wrong, and in which direction?

Answer form: one modulus and a sentence, two numbers, a derivation and a number, a number and a direction.

## 6. (Judgment: what the numbers can and cannot tell you)

Answer each in two or three sentences, with a number where one is asked for.

(a) A certificate for a lot of 2024-T3 sheet reports E = 52 GPa and Fty = 300 MPa. List the two questions you ask the laboratory before rejecting the lot, and say which of the two numbers you would distrust first.
(b) A crashworthy seat must absorb the kinetic energy of 90 kg moving at 4 m/s by plastic deformation of a metal strap. Estimate the toughness of each alloy in `MATERIALS` as (Fty + Ftu)/2 × elongation and find the minimum strap mass for each. Which alloy wins per kilogram, and why is the modulus of resilience the wrong property for this job?
(c) Why does a stress engineer use an allowable rather than the typical value, even when the typical value comes from thousands of tests?
(d) Name one manufacturing process for which the true stress–strain curve beyond the uniform elongation matters, and explain why the engineering curve cannot provide it.

Answer form: four short paragraphs, four masses.
