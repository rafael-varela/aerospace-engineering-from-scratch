# P09.13  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact, which come from your code, and which are simulations of a representative model.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: a second-order sensor in two realisations)

An accelerometer's mechanical element behaves as ÿ + 2ζω_n ẏ + ω_n² y = ω_n² u, with ω_n = 2 + s/10 rad/s (a deliberately slow, scaled example) and ζ = 0.3.

(a) Write the controllable and the observable canonical forms by hand. (b) Evaluate G(jω_n) by hand from the TF and with `ss_eval` from both realisations. Explain why |G(jω_n)| = 1/(2ζ). (c) Write a third realisation with the physical states (y, ẏ) and show by (5) of the lesson which T relates it to the controllable form.

Answer form: two sets of (A, B, C, D); one complex number with a one-line reason; one matrix T.

## 2. (Derivation: invariance, duality, Faddeev–LeVerrier, the aircraft TF)

(a) Prove that the transfer function and the eigenvalues are unchanged by x = Tz. (b) Derive the controllable canonical form (3) for a biproper G (deg num = deg den) using the internal variable w, and show where D = b₀ comes from. (c) Show that the observable form (Aᵀ, Cᵀ, Bᵀ, D) has the same transfer function. Explain why this argument needs G to be scalar. (d) Derive the recursion N_k = AN_{k−1} + c_{k−1}I by expanding (sI − A) adj(sI − A) = det(sI − A)I, and show that the s⁰ coefficient is the Cayley–Hamilton theorem. (e) Derive θ/δ and h/δ of the lesson's model by hand and give their DC behaviour (type, high-frequency roll-off).

Answer form: five derivations.

## 3. (Numerical, uses your code: a hidden mode)

Take G(s) = (s + 1)/((s + 1)(s + 2)), and do not cancel.

(a) Build its controllable canonical form and give A, C and the eigenvalues. (b) Simulate the unit step for 5 s with dt = 10 ms and compare with the one-state model ẋ = −2x + u, y = x. Give the largest difference and y(5 s). (c) Now start the 2-state model from x₀ = (1, 1) with zero input and give y at 0, 1 and 5 s. Which mode does the output show? (d) Compute the ranks of [C; CA] and [B, AB]. Say in one sentence what P09.14 will call this realisation.

Answer form: two matrices and two eigenvalues; two numbers; three numbers and one sentence; two ranks and one sentence.

## 4. (Unfamiliar scenario: a satellite with a flexible solar array)

A satellite hub (inertia J₁ = 3 kg m²) carries a solar array (J₂ = 1 + s/10 kg m²) through a torsional spring k = 2 N m/rad with damping c = 0.02 N m s/rad. The reaction wheel torque u acts on the hub, and the star tracker measures the hub angle θ₁. These are representative numbers.

(a) Write the equations of motion and the state-space model with x = (θ₁, θ₂, ω₁, ω₂). (b) Get the TF θ₁/u with `ss_to_tf`, and give its poles and zeros. Identify the rigid-body poles, the flexible-mode resonance and the antiresonance (the zero), and check them against √(k(J₁ + J₂)/(J₁J₂)) and √(k/J₂). (c) Simulate a 0.01 N m torque step for 30 s and compare θ₁(30 s) with SciPy and with the rigid-body result u t²/(2(J₁ + J₂)). (d) Can `modal_form` be applied to this A? Report cond(V) and explain.

Answer form: two equations and four matrices; a TF with four poles and two zeros; three numbers; one number and one sentence.

## 5. (Judgment: which coordinates, which model)

(a) A colleague stores a 10th-order flexible-aircraft model in companion form with coefficients printed to five significant figures. Use the lesson's sensitivity table to say what can go wrong, and recommend a better representation. (b) Give two situations in which the transfer function of the lesson's aircraft is enough, and two in which you need the state-space model. (c) The pitch-altitude model predicts a steady pitch rate under a held elevator. What prediction would a 30-second simulation make for altitude, and why should nobody believe it? Name the missing state and the derivative that would fix it. (d) A vertical gust adds an angle w_g (rad) to the aerodynamic angle of attack, so the force and moment terms see α + w_g while the kinematics (rows 3 and 4) see the inertial α. Write the new column of B for w_g, derive the TF from w_g to h, and say which part differs from h/δ.

Answer form: four short paragraphs with specific numbers or equations.

## 6. (Hard: how accurate is Faddeev–LeVerrier?)

(a) For A = Q diag(1, 2, …, n) Qᵀ with Q a random orthogonal matrix (seed 1315, from the QR factorisation of a standard normal matrix), compute the characteristic polynomial by `faddeev_leverrier` and by `np.poly`. Give the largest relative coefficient error against the exact integer coefficients for n = 5, 10, 15, 20. (b) Repeat with A = diag(1, …, n) (Q = I). Explain why the error is zero there up to n = 15. (c) For diag(1, …, n), compute the roots of the Faddeev–LeVerrier coefficients with `np.roots` and give the largest root error for n = 5, 10, 15, 20. Explain why a polynomial with *exact* coefficients can still have inaccurate computed roots (Wilkinson's polynomial). (d) What do (a)–(c) say about using `ss_to_tf` as a step in a design pipeline for high-order models?

Answer form: two four-entry tables, one four-entry table, two explanations, one recommendation.
