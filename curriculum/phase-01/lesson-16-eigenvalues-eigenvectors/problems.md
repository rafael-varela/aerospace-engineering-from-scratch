# P01.16  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: a 2×2 tensor by hand)

A flat plate has the in-plane inertia tensor (kg m²)

$$
A = \begin{pmatrix} 20 + s & -4 \\ -4 & 10 + s/2 \end{pmatrix}.
$$

(a) Write the characteristic polynomial from the trace and determinant and solve it. (b) Find a unit eigenvector for each eigenvalue and the angle each makes with the x axis. (c) Verify the two eigenvectors are perpendicular and that their sum and product match the trace and determinant. (d) Confirm with `jacobi_eigen` and with `np.linalg.eigh`.

Answer form: two eigenvalues, two angles, two checks.

## 2. (Derivation: the spectral theorem and the triangle inequality)

(a) Prove that the eigenvalues of a real symmetric matrix are real.
(b) Prove that eigenvectors belonging to distinct eigenvalues of a symmetric matrix are orthogonal.
(c) Prove that for any unit vector **n**, **n**ᵀ I **n** = Σ m_k |**r**_⊥,k|², the moment of inertia about the axis **n** through the reference point.
(d) Prove I_xx + I_yy ≥ I_zz from the definition, state when equality holds, and explain why the inequality then holds for the principal moments in any frame.

Answer form: four short proofs.

## 3. (Unfamiliar scenario: natural frequencies of a two-mass chain)

Two masses m₁ = 2 kg and m₂ = 1 kg are joined in a line by three identical springs of stiffness k = 1000 + 100s N/m, the outer springs anchored to walls. The equations of motion are M **ẍ** + K **x** = 0 with M = diag(m₁, m₂) and K = k [[2, −1], [−1, 2]].

(a) Show that the substitution **x** = M^{−1/2} **y** turns the problem into the *symmetric* eigenvalue problem S **y** = ω² **y** with S = M^{−1/2} K M^{−1/2}. Why is the symmetric form preferable to M⁻¹K, which has the same eigenvalues?
(b) Compute S, its eigenvalues with `jacobi_eigen`, and the two natural frequencies in Hz.
(c) Recover the mode shapes **x** = M^{−1/2} **y** and normalise each so that x₁ = 1. In which mode do the masses move in opposite directions?
(d) Verify K **x** = ω² M **x** for each mode, and verify that the two mode shapes are M-orthogonal: **x**₁ᵀ M **x**₂ = 0.

Answer form: a 2×2 matrix, two frequencies, two mode shapes, two checks.

## 4. (Numerical, uses your code: a payload is added)

A 6 kg payload is bolted to the satellite at (0.40, −0.30, s/50) m in body axes.

(a) Recompute the centre of mass, the inertia tensor and the principal moments with your code.
(b) By how many degrees has the major axis moved relative to its direction without the payload? How far is it now from body z?
(c) Run `power_iteration` on the new tensor and compare the iteration count with log(10⁻¹²)/log(λ₂/λ₁).
(d) Is a spin about body z closer to or further from the major principal axis than before? Give the moment of inertia about body z (Rayleigh quotient) and the major moment.

Answer form: a tensor, three moments, two angles, two iteration counts, one judgment with numbers.

## 5. (Judgment: spin stabilisation)

(a) The satellite of the lesson was drawn with body z as its intended spin axis, but the major principal axis is 19° away. Describe, using **H** = I **ω**, what happens if it is spun about body z anyway, and give the angle between **H** and **ω** for that spin.
(b) The minor axis (moment 16.26 kg m²) has the highest spin rate for a given angular momentum. Why would a designer still refuse to spin about it, and which real event demonstrates the point?
(c) Name a condition under which `spin_axis_stability` refuses to answer, and say what physics replaces the rule there.
(d) If the propellant tank drains from 16 kg to 4 kg during the mission, does the major axis stay the major axis? Compute both cases with your code and comment on the margin.

Answer form: four short paragraphs, each with a number.

## 6. (Hard, uses your code: eigenvalues and stability of a damped system)

A one-degree-of-freedom system has ẍ + 2ζωẋ + ω²x = 0 with ω = 2 rad/s and ζ = s/100. Written as d**u**/dt = A**u** with **u** = (x, ẋ), A = [[0, 1], [−ω², −2ζω]].

(a) Find the eigenvalues of A by hand and show that their real part is −ζω and their imaginary part ±ω√(1 − ζ²).
(b) Compute `linear_system_decay_rates(A)` and confirm. Then try `power_iteration(A)` and explain its behaviour from the eigenvalue magnitudes.
(c) Compute **u**(t) = exp(At) **u**(0) at t = 5 s for **u**(0) = (1, 0) with `scipy.linalg.expm` and show that |**u**(5)| is below e^{−5ζω} times a constant you can justify.
(d) Explain why `matrix_function_symmetric` cannot be used for this A, and what property of A the Jacobi method needs that A lacks.

Answer form: a derivation, two numerical checks, two short explanations.
