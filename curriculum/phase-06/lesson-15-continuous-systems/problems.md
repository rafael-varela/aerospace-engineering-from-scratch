# P06.15  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31). Material data are typical, not design allowables; structural numbers are illustrative.

## 1. (Warm-up: an instrument mast)

A spacecraft's magnetometer mast is a carbon-fibre tube (E = 100 GPa, ρ = 1600 kg/m³), 60 mm outside and 56 mm inside diameter, L = (1 + s/10) m long, clamped at the spacecraft and carrying a 0.2 kg sensor at its tip.

(a) Its A, I, EI, ρA, own mass and the mass ratio μ = m_tip/(ρAL).
(b) The first three bending frequencies without the sensor, from (3) and (4).
(c) The first three frequencies with the sensor, from the tip-mass equation (5) with `natural_frequencies(..., m_tip=0.2)`.
(d) The SDOF estimate of `P06.13` for f₁ with the sensor, and its error in per cent.

Answer form: six numbers; three frequencies; three frequencies; a frequency and a percentage.

## 2. (Derivation: the beam toolkit)

(a) Derive EI w⁗ + ρA ẅ = 0 from the equilibrium of a slice with d'Alembert's inertia force, starting from the results of `P05.12`.
(b) Derive the cantilever frequency equation (4) and its mode shape from the general solution of φ⁗ = β⁴φ.
(c) Derive the tip-mass shear condition EI φ‴(L) = −m_t ω² φ(L) from a free-body diagram of the mass, and show that it leads to (5).
(d) Prove the orthogonality relation (6) for a cantilever by integrating by parts twice. Where is each boundary condition used?
(e) Show that sin(jπxᵢ/L) is an exact eigenvector of the pinned finite-difference matrix of 2.4 and derive its eigenvalue (7).

Answer form: five derivations.

## 3. (Unfamiliar scenario: a rocket in flight)

In flight a slender launch vehicle is a free–free beam. Model one as uniform, L = (30 + s/2) m, EI = 5 × 10⁸ N·m², ρA = 800 kg/m (illustrative).

(a) Write the free–free boundary conditions and show that the frequency equation is cos βL cosh βL = 1, the same as for the clamped–clamped beam. What extra solutions does the free–free beam have?
(b) The first three flexible frequencies from `beta_L("clamped-clamped", 3)`.
(c) Check them with the finite-element model with no supports at all (all degrees of freedom free; `scipy.linalg.eigh(K, M)` on the full matrices). What are the first two eigenvalues, and why?
(d) The vehicle's guidance uses a rate gyro. Where along the vehicle (as x/L) should it be mounted so that it senses as little as possible of the first bending mode's *deflection*? Is that the right criterion for a rate gyro, which senses rotation? Where would you put it then?

Answer form: conditions and a sentence; three frequencies; two eigenvalues, three frequencies and a reason; two positions and a sentence.

## 4. (Numerical, uses your code: a model fit for purpose)

Use the lesson's wing spar.

(a) How many cubic elements does `fem_frequencies` need to give each of the first five frequencies to 10⁻⁴ relative? Tabulate the error of each mode for 3, 4, 6, 8, 10, 12, 16 and 20 elements, and state a rule in elements per mode.
(b) A tank of (20 + s) kg can go anywhere between the root and the tip. Sweep its position on a 0.05 m grid (100 elements). Where must it go to lower f₂ most, and where to leave f₂ unchanged? Compare the second position with the node of the second mode from `mode_shape`.
(c) Explain both results with the Rayleigh quotient of `P06.14`.

Answer form: a table and a rule; two positions with frequencies and a comparison; a paragraph.

## 5. (Judgment: is Euler–Bernoulli enough?)

(a) For the lesson's spar (depth about 0.18 m on a 5 m span), estimate the mode number at which the half-wavelength of a bending mode falls to about ten depths. Why is that a rough limit of Euler–Bernoulli theory?
(b) A ground vibration test measures the first bending frequency of the real wing at 8.5 Hz, against the model's 6.22 Hz. Name three model assumptions that could explain the difference and the direction of each.
(c) A colleague proposes a lumped-mass FE model with 8 elements "because it is simpler". What error should they expect in f₁ and f₃ (experiment 1), and in which direction?
(d) Which result of this lesson would you trust enough to use in a design decision, and which would you not?

Answer form: four labelled items.
