# P01.16  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. For a 2×2 matrix the characteristic polynomial is λ² − (tr A)λ + det A, and an eigenvector of [[a, b], [b, d]] for eigenvalue λ is (b, λ − a) or (λ − d, b); normalise and take atan2.
2. For (a) and (b), the key move is **v**ᵀA**w** = (A**v**)ᵀ**w**, which needs A = Aᵀ. For (c), expand **n**ᵀ[(**r**·**r**)I − **rr**ᵀ]**n**. For (d), write out I_xx + I_yy − I_zz as a sum over masses and look at what is left.
3. S = M^{−1/2} K M^{−1/2} is symmetric because K is; M^{−1/2} is just diag(1/√m₁, 1/√m₂). The mode shapes of the original problem are M^{−1/2} times the eigenvectors of S, not the eigenvectors themselves.
4. Append the payload to the arrays returned by `satellite_point_masses()`; everything else is a function call. The angle between two unit vectors is arccos of their dot product; use the absolute value because an eigenvector's sign is arbitrary.
5. For (a), compute **H** = I **ω** with **ω** along body z and take the angle between **H** and **ω**. For (d), change the tank mass in the position/mass arrays and look at the gap between the two largest moments.
6. For a 2×2 matrix, λ = (tr A ± √(tr²A − 4 det A))/2. Power iteration on a matrix whose two eigenvalues have the *same modulus* has nothing to converge to. Jacobi needs a symmetric matrix; A is not, and its eigenvalues are not even real.
