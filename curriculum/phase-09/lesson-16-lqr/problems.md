# P09.16  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact (closed form) and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: the scalar regulator by hand)

A first-order unstable mode ẋ = a x + b u with a = s/10 s⁻¹, b = 1, cost J = ∫(q x² + r u²) dt with q = 1, r = 0.1.

(a) Write the scalar Riccati equation 2aP − b²P²/r + q = 0 and solve it for the stabilising root P > 0. Give P, the gain K = bP/r and the closed-loop pole a − bK, and show the pole equals −√(a² + b²q/r). (b) Take the limit r → ∞ (expensive control) and r → 0 (cheap control) of the pole and of K. Explain in one sentence each why the expensive-control pole is −|a| and not 0. (c) Check (a) with `care_hamiltonian` and `care_newton_kleinman` on the 1 × 1 matrices and with `scipy.linalg.solve_continuous_are`.

Answer form: three numbers with a one-line derivation, two limits with one sentence each, three agreeing numbers.

## 2. (Derivation: why the Hamiltonian works, and why LQR has margins)

(a) For J = ∫(xᵀQx + uᵀRu) dt subject to ẋ = Ax + Bu, introduce a costate λ and show that the optimality conditions give u = −½R⁻¹Bᵀλ and the linear system d/dt(x, λ/2) = H (x, λ/2) with H = [[A, −BR⁻¹Bᵀ], [−Q, −Aᵀ]]. (b) Show that if (x, λ/2) evolves in an invariant subspace spanned by the columns of [X₁; X₂], then λ/2 = Px with P = X₂X₁⁻¹, and that substituting into the costate equation gives the CARE. (c) Show that the eigenvalues of H come in pairs λ, −λ (use J H J⁻¹ = −Hᵀ with J = [[0, I], [−I, 0]]), so the stable half exists whenever no eigenvalue lies on the imaginary axis. (d) For a single input, derive the return-difference identity (1 + L(−jω))R(1 + L(jω)) = R + G(−jω)ᵀQG(jω) with L = K(jωI − A)⁻¹B and G = (jωI − A)⁻¹B, conclude |1 + L(jω)| ≥ 1, and deduce from the Nyquist plot that the gain may be multiplied by any k ∈ (½, ∞) and the phase margin is at least 60°.

Answer form: four derivations, each ending in a boxed equation or inequality.

## 3. (Unfamiliar scenario: roll attitude of a small aircraft)

Roll dynamics of a light UAV, about one axis and decoupled: φ̇ = p, ṗ = L_p p + L_δ δ_a, with roll damping L_p = −2 s⁻¹ and aileron effectiveness L_δ = (15 + s) s⁻² (representative numbers).

(a) Use Bryson's rule with φ_max = 10°, p_max = 30°/s, δ_max = 15° to write Q and R. (b) Compute K, the closed-loop poles, their damping ratio and natural frequency with `lqr`. (c) For a 10° initial bank error, give the peak aileron deflection from `simulate` and compare it with K₁ × 10°. (d) Compute the phase margin and crossover with `loop_transfer` and `phase_margin`, and confirm the Kalman inequality numerically. (e) Repeat (b)–(d) with R divided by 10. What did the cheaper control buy and what did it cost?

Answer form: two matrices, a gain with poles, a peak deflection, a margin and crossover, and a two-row comparison table with one sentence.

## 4. (Numerical, uses your code: your satellite)

Your satellite has J = (5 + s) kg·m² about the axis of interest and the lesson's wheel (0.02 N·m, 0.2 N·m·s).

(a) With the lesson's Bryson weights (5°, 0.5°/s, 0.02 N·m), give K, the closed-loop poles and the hand formulas k₁ = √(q₁/r), k₂ = √(q₂/r + 2Jk₁). (b) For a 5° error, give the peak torque, the 2 % settling time of θ and the peak change in wheel momentum. (c) For a 15° error with the 0.02 N·m limit, give the time spent at the limit and the 2 % settling time; compare with the unlimited loop. (d) How many Newton–Kleinman iterations does `care_newton_kleinman` take from the Bass gain, and how far is its P from the Hamiltonian P?

Answer form: a gain, two poles, three numbers, a four-entry comparison, an iteration count and a relative difference.

## 5. (Judgment: what the guarantee does not cover)

The lesson's satellite design has a phase margin of 77° and an infinite upper gain margin. (a) The guaranteed margins assume the full state is measured without delay and the loop is broken at the plant input. List three things on a real spacecraft that violate those assumptions. (b) The rate is in practice estimated by a star tracker and gyro through an observer (`P09.15`). What happens to the guarantee (name the known result), and what would you check instead? (c) The delay margin of the design is PM/ω_c. Compute it and say whether a 0.25 s attitude-determination latency is a concern. (d) A flexible solar array has a lightly damped mode at 0.8 rad/s. Is the rigid-body LQR design still trustworthy? Give one quantitative reason.

Answer form: four short paragraphs with specific numbers or mechanisms.

## 6. (Hard: the ½ gain margin is tight, and the cost of a wrong gain)

(a) For the scalar plant ẋ = a x + b u (a > 0) with q fixed, show that as r → ∞ the LQR gain tends to K = 2a/b and that the loop with gain kK is stable if and only if k > ½. So the lower gain margin ½ cannot be improved in general. (b) For the lesson's satellite, compute J(kK) = x₀ᵀX(k)x₀ with `cost_of_gain` for k = 0.25, 0.5, 0.8, 1, 1.25, 2, 4 and x₀ = (5°, 0). Show the minimum is at k = 1, and observe that J(kK) = J(K/k) for these k. (c) Find out what the symmetry belongs to: repeat (b) for k = 0.5 and 2 with Q replaced by diag(q₁, 0) and by diag(q₁, 10q₂), and with the initial state x₀ = (0, 0.1°/s). Say which change breaks it, and give one sentence on why a symmetric cost curve is good news for a designer who is unsure of the wheel's torque constant.

Answer form: a limit with a stability proof, a seven-row table, and a six-entry table with two sentences.
