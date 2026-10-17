# P06.14  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31). All structural numbers are illustrative.

## 1. (Warm-up: a two-mass instrument stack)

An instrument of m₂ = 300 kg sits on a k₂ = 2 × 10⁸ N/m bracket on a platform of m₁ = (200 + 5s) kg, which is joined to the spacecraft by k₁ = 5 × 10⁸ N/m (fixed base).

(a) Write M and K. Solve det(K − ω²M) = 0 by hand for both frequencies (Hz).
(b) Find both mode shapes (x₂/x₁), mass-normalise them, and verify φ₁ᵀMφ₂ = 0 numerically.
(c) Compute the participation factors and effective masses for base excitation and check that they sum to the total mass.
(d) Check (a)–(c) with `modes` and `participation`.

Answer form: two matrices and two frequencies; two ratios, two vectors and a number; two factors, two masses and a sum; a check.

## 2. (Derivation: the modal toolkit)

(a) Prove the orthogonality relations (4) from (2) and the symmetry of M and K. What happens when two frequencies are equal?
(b) Show that for mass-normalised modes the effective masses Γᵢ² sum to rᵀMr. (Hint: expand r in the modes.)
(c) Derive the apparent mass (12) for proportional damping: start from (10), write the absolute acceleration of the masses, and sum the inertia forces.
(d) Show that with C = aM + bK the modal damping ratios are (8). Sketch ζ(ω) and say which modes are over- and under-damped relative to the two fitted ones.

Answer form: four derivations and a sketch.

## 3. (Unfamiliar scenario: the stack after separation)

After separation from the launch vehicle the lesson's stack floats free in orbit: no spring k₁ to the interface.

(a) Write the free–free stiffness matrix (only k₂ and k₃). Show that it is singular and find its null vector. What motion is it?
(b) Solve K φ = ω² M φ (`scipy.linalg.eigh` works with a singular K). Report the three frequencies and shapes. Compare the flexible frequencies with the fixed-base ones of the lesson, and explain why they are not simply the higher two of those.
(c) Show that the flexible modes carry zero linear momentum (Σmᵢφᵢ = 0) and explain why from orthogonality with the rigid mode.
(d) A thruster pulse excites the stack. Which modes can it excite, and what happens to the rigid-body mode?

Answer form: a matrix, a vector and a sentence; three frequencies, three shapes and an explanation; a check and a reason; a sentence.

## 4. (Numerical, uses your code: thrust onset)

At engine ignition the interface acceleration jumps from 0 to a₀ = (1 + s/10) g and stays there (a step). Use the lesson's stack with Rayleigh damping of 2 % in modes 1 and 2.

(a) In relative coordinates, M z̈ + C ż + K z = −M r a₀. What are the static spring forces (the forces each spring carries when the stack has settled)?
(b) Compute the response for 0.1 s with `direct_integration` and with `modal_integration` using 1, 2 and 3 modes. Report the peak force in each spring and the dynamic amplification factor (peak over static) for each.
(c) An SDOF system under a step overshoots by a factor of 2 (undamped). Why can a spring in the stack exceed 2?
(d) What is the payload's peak absolute acceleration in g?

Answer form: three forces; a table of peaks and factors for four models; a sentence; a number.

## 5. (Numerical, uses your code: a damper that is not proportional)

A viscous damper of c_d = 10⁴(1 + s/5) N·s/m is fitted between the service and payload modules of the lesson's stack, with no other damping.

(a) Compute the diagonal modal damping ratios from ΦᵀCΦ and the exact damping ratios and damped frequencies from the eigenvalues of the state matrix [[0, I], [−M⁻¹K, −M⁻¹C]].
(b) Compare the payload driving-point receptance peak near f₁ from `frf_direct` with that from `frf_modal` using the diagonal ratios.
(c) With this damper, what payload acceleration does a 1 g interface sine at f₁ produce, and what interface force?
(d) The damper barely damps mode 1. Use the mode shapes to explain why, and suggest where a damper would help mode 1 more.

Answer form: two sets of three ratios and three frequencies; two peaks and a ratio; an acceleration and a force; an explanation.

## 6. (Judgment: modal models in a test campaign)

Give a short argued answer to each.

(a) A colleague's 50-mode finite-element model and the lesson's 3-mass model both predict f₁ within 2 % of the test. Which would you use to set the notch, and why?
(b) The measured interface force at f₁ is 40 % below the prediction. Name two explanations that the model can test and one that it cannot.
(c) Why does a launch-vehicle user's guide ask for the spacecraft's fundamental frequency and its effective masses rather than its full mass and stiffness matrices?
(d) A test engineer proposes to skip the low-level sweep because "we already have the modal model". What do you answer?

Answer form: four labelled items.
