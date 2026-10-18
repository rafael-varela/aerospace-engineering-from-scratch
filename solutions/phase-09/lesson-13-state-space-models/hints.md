# P09.13  Hints

One hint per problem. Read only the one you need.

1. With w = y, the internal variable of the companion form is the output scaled by 1/ω_n² (C = [0, ω_n²]). At s = jω_n the real part of the denominator, ω_n² − ω², vanishes and only the damping term is left. For T, write the physical state (y, ẏ) as a linear combination of the companion state (ẇ, w).

2. (a) (sI − T⁻¹AT) = T⁻¹(sI − A)T, and det(T⁻¹MT) = det M. (b) With D = b₀, write y = b₀w⁽ⁿ⁾ + ⋯ and replace w⁽ⁿ⁾ using the ODE. (c) A scalar equals its own transpose. (d) Write adj(sI − A) as Σ N_k s^{n−k} and match the coefficient of each power of s. (e) Eliminate α and q in the Laplace domain.

3. `tf_to_controllable(np.array([1.0, 1.0]), np.polymul([1, 1], [1, 2]))`. The zero at −1 cancels a pole in the TF but not in A. The initial-condition response shows the modes the output can see: look at the decay rate.

4. Write J₁θ̈₁ = u − k(θ₁ − θ₂) − c(θ̇₁ − θ̇₂) and J₂θ̈₂ = k(θ₁ − θ₂) + c(θ̇₁ − θ̇₂). The rigid-body mode has two poles at the origin. At the antiresonance the array moves so that the hub stands still. A double eigenvalue at 0 with a single eigenvector means V is singular.

5. Look at how the pole error grew from n = 6 to n = 10 in the lesson. A TF is enough for steady-state and frequency-response questions from rest; a state model is needed for initial conditions, internal limits and multiple inputs and outputs. For the gust, find where a change in α enters rows 1 and 4.

6. The exact coefficients of Π(s − k) are integers, which Python can compute exactly. Faddeev–LeVerrier divides by k and accumulates traces: with Q = I everything stays an integer. For the roots, perturb one coefficient by a relative 10⁻¹⁰ and watch them move.
