# P03.13  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact, which are numerical, and what the numerical error is.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: a two-element bar by hand)

A steel tie rod (E = 200 GPa, A = 200 mm², L = 1 m) is fixed at one end and loaded at the other by P = 10 + s kN. Model it with two equal elements.

(a) Write the two element stiffness matrices, assemble the 3 × 3 global matrix, apply the constraint and solve for the nodal displacements by hand.
(b) Give the stress in each element and the reaction. Check against σ = P/A and u(L) = PL/(EA).
(c) Now add a uniform body force q = 5 kN/m (no tip load). Write the consistent load vector and solve again. Compare the tip displacement with the exact qL²/(2EA). Why does it agree exactly?

Answer form: the matrices, three displacements, two stresses, one reaction, one sentence.

## 2. (Derivation: the weak form and what follows from it)

(a) Starting from −(EAu')' = b with u(0) = 0 and EAu'(L) = P, derive the weak form (2) and state which boundary condition is built into the space and which appears in the equation.
(b) With hat functions, derive the element stiffness matrix (4) and the consistent mass matrix (ρAh/6)[[2, 1], [1, 2]] by direct integration in the parent coordinate.
(c) Prove that the assembled K is symmetric and positive semi-definite, that the vector of ones is in its null space, and that fixing any one node makes it positive definite.
(d) Show that Ku = f is the stationarity condition of Π(u) = ½uᵀKu − fᵀu, that the stationary point is a minimum, and that at the minimum ½uᵀKu = ½fᵀu.

Answer form: four derivations.

## 3. (Unfamiliar scenario: thermal strain)

A 2 m aluminium bar (E = 70 GPa, A = 600 mm², α = 23 × 10⁻⁶ K⁻¹) is heated uniformly by ΔT = 50 K. The stress–strain law is σ = E(u' − αΔT).

(a) Insert this law into the weak form and show that the thermal strain appears as an equivalent element load f_e = EAαΔT [−1, 1]ᵀ, with no change to K.
(b) Using the lesson's code with four elements, solve the bar with both ends fixed. Give the nodal displacements, the mechanical strain u', the total stress σ and the reactions. Check against σ = −EαΔT.
(c) Solve again with the far end free. Give the tip displacement and the stress, and check against u(L) = αΔTL.
(d) A bar with both ends fixed and a *linearly varying* ΔT(x) from 0 to 50 K: without running anything, say what the stress distribution is and why the four-element model gets it exactly.

Answer form: a derivation, two short tables, one paragraph.

## 4. (Numerical, uses your code: your own rotor)

Repeat the spar case with Ω = 380 + s rpm and m_tip = 1.5 + s/10 kg, all else as in `rotor_blade`.

(a) Give the tip load, the spar's own centrifugal load and the exact root stress, tip extension and hub reaction from (5).
(b) Solve on n = 8 elements. Report the tip displacement, the root element stress, the exact stress at the first element's centre and at the root, and the reaction. Which numbers are exact and which carry a discretisation error?
(c) Extrapolate the first two element stresses linearly to x = 0 and compare with the exact root stress. Is the extrapolation better or worse than the element value, and why?
(d) Run `convergence_study` on n = 4, 8, 16, 32 and report the observed orders.

Answer form: five numbers for (a), a short table for (b), two numbers for (c), two orders for (d).

## 5. (Numerical, uses your code: consistent against lumped mass)

For the spar (not spinning), compute the first three axial frequencies on n = 4, 8, 16, 32 with the consistent mass matrix and with a lumped one (the diagonal of row sums of M).

(a) Tabulate the ratios ω_FE/ω_exact for both.
(b) Show that the consistent ratios are above 1 and decrease with n while the lumped ones are below 1 and increase. Fit the first-mode error of each to n^{−p} and report p.
(c) How many elements per half-wavelength does each model need for 1 % accuracy on a given mode? Use the third mode to answer.
(d) The blade spins at 400 rpm = 6.67 Hz. Is any axial frequency near a rotor harmonic kΩ for k ≤ 10? Say why the question, for a real blade, must be asked of the bending modes instead.

Answer form: two 4 × 3 tables, two exponents, one number, one paragraph.

## 6. (Judgment: the tapered strut and the limits of the bar)

A landing-gear drag strut tapers linearly with c = 0.3 + s/100 from A₀ = 800 mm² over L = 1.5 m and carries P = 50 kN.

(a) How many uniform elements put the last element's stress within 2 % of the exact tip stress? Within 0.5 %? Tabulate n = 4 … 128.
(b) Show (from EAu' = P and the mid-element area) that the last element's stress equals the exact stress at the element's centre to rounding, and explain why the tip error therefore falls only as 1/n.
(c) Propose a recovery of the tip stress that is exact for this problem without any refinement, and say what it relies on.
(d) The strut carries bending as well as axial load in service, and its ends are pinned in clevises. List, in order of importance, the three assumptions of the bar model that fail for the real part, with one measurable consequence each.

Answer form: a table, a derivation, one proposal, one ranked list.
