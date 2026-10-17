# P03.06  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. The fast mode lives almost entirely in the skin node, so its eigenvalue is close to the skin's own rate −k_c/C_s, corrected by the small fraction of the coupling that the structure takes up. For (c), |R| for backward Euler is 1/|1 − z|: evaluate it with z = −300-ish and z = −0.01-ish and ask what each does to the corresponding mode in one step.
2. (a) Set C h_new^q = tol·(C h^q)/E and solve; the safety factor and clamp are about what the estimate does *not* know. (b) Put y' = λy into the rule and solve for y_{n+1}/y_n; for A-stability compare |1 + z/2| with |1 − z/2| as distances in the complex plane from the points −2 and +2. (c) The characteristic polynomial is λ² − (tr A)λ + det A; the determinant of a 2×2 is a two-term difference in which the coupling term cancels.
3. The actuator eigenvalue is −1/τ; the aircraft modes are ±2πi/T. Stiffness is about the ratio of the largest to the smallest |Re λ| *or* |λ|; here the slow modes are oscillatory, so compare |λ|. The stable step is L/|λ_fast| with L = 3.307.
4. a_max ∝ |sin γ_E|: compute it before touching the code. The peak read from stored points is a lower bound on the true peak; the finer the stored grid, the closer. For (d), count how many steps are start-up (growing by 5× from 0.02 s) and how many are near the peak.
5. Copy `implicit_euler_step` and change two lines (G and G′). The first trapezoidal skin temperature after a 60 s step from 400 K is (1 + z/2)/(1 − z/2) applied to the 100 K excess with z ≈ −301; write that number down before running.
6. The coupled Jacobian is block-structured; the thermal block's fast eigenvalue is unchanged by the coupling because the trajectory does not depend on temperature. Count evaluations as steps × stages, and steps as (duration)/(stable step) for the explicit solvers.
