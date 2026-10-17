# P09.15  Hints

One hint per problem. Read only the one you need.

1. det(sI − A + BK) = s² + (k₂/J)s + k₁/J. For the observer, det(sI − A + LC) = s² + l₁s + l₂ with C = [1, 0]. The peak torque of a step with a perfect estimate is at t = 0 and equals Nr.

2. (a) In the canonical form, B = e₁, so BK only changes the first row. (b) Apply Cayley–Hamilton to A_K = A − BK: φ_d(A_K) = 0. Expand A_Kᵏ = Aᵏ − Σ A^{k−1−j}BKA_Kʲ and collect the terms of the form AʲB(…). (c) A matrix and its transpose have the same eigenvalues. (d) Substitute x̂ = x − e. (e) Set ẋ = 0.

3. `ackermann(*pitch_model()[:2], poles)`, `observer_gain(A, C, 2 * poles)`, `simulate_observer_loop` with r = 5° held. The peak is at t = 0: u(0) = Nr. Look at the size of N.

4. A = [[0, 1], [a, 0]], B = (0, b)ᵀ. The closed-loop polynomial is s² + bk₂s + (bk₁ − a). For the gyro, compute [C; CA] with C = [0, 1]: CA brings in a.

5. Build the true loop matrix [[A_true, −B_true K], [LC, A − BK − LC]] and take its eigenvalues. For the steady state, solve the linear system with the reference input. For (c), use the ω₀^{3/2} law.

6. For a complex pair, L = (2σ', σ'² + ω²) with σ' = 0.3 for C = [1, 0]. In modal coordinates each component decays exactly as e^{σt}. `scipy.signal.argrelmax` finds the local maxima of ‖e‖.
