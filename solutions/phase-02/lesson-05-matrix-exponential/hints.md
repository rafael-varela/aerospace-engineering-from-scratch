# P02.05  Hints

One hint per problem. Read only the one you need.

1. N² = 0, so the series stops after two terms. For the rotation generator, A² = −ω²I, so A^{2k} = (−1)^k ω^{2k} I and A^{2k+1} = (−1)^k ω^{2k} A: group the even powers (a cosine times I) and the odd powers (a sine times A/ω). A diagonal matrix exponentiates entry by entry. `expm_eig` needs a full set of eigenvectors.

2. Differentiate the series term by term and factor A out on the left or the right. For the group property, both X(t) = Φ(t)Φ(s) and Y(t) = Φ(t + s) satisfy Ẋ = AX with the same value at t = 0. For the determinant, (d/dt) det Φ = det Φ · tr(Φ⁻¹Φ̇) = det Φ · tr(A). For the commutation proof, the coefficient of A^j B^k in e^A e^B is 1/(j!k!) and in e^{A+B} is C(j+k, j)/(j+k)! when A and B commute; they agree. For the augmented matrix, M^k = [[A^k, A^{k−1}B], [0, 0]] for k ≥ 1.

3. A = −dI + λJ with J = [[0, 1], [−1, 0]]; the identity commutes with everything, so Problem 2(d) applies and e^{λJt} is the rotation of Problem 1(b). Amplitude decays as e^{−dt} regardless of λ; cycles = t_half × λ/(2π).

4. At the same altitude with no relative velocity the chaser is in the same orbit: look at the y₀ column of Φ. For the plan, call `rendezvous_delta_v(n, r0, v0, 0.4·P)` and `propagate` the post-burn state. To catch up the chaser must go to a lower (faster) orbit, which is x < 0. The singularity: det Φ_rv = [8(1 − cos nT) − 3nT sin nT]/n² vanishes at nT = 2π.

5. Count operations: one eigen-decomposition then 10 000 cheap diagonal exponentials versus 10 000 full exponentials; and remember which matrix is defective. For (b), compare the size of the largest Taylor term (‖At‖^k/k! near k = ‖At‖) with the size of the answer (entries of order e^{Re λ t}). For (c), think linearisation error (quadratic in separation), burn execution error (a fraction of Δv, growing over the transfer), differential drag; (d) is about what each method assumes the input does *between* samples.

6. The eigenvalues of f(A) are f(λ) for any function given by a power series, so e^{Ah} has eigenvalues e^{λh} and I + hA has 1 + λh. The stability limit comes from |1 + λh|² = (1 + h Re λ)² + (h Im λ)² < 1. For the bound, use `np.linalg.norm(A, 1)`, `math.factorial`, and compare with `np.abs(expm_taylor(A, N) - expm(A)).max()`; the squarings follow from ⌈log₂(‖A‖₁/0.5)⌉.
