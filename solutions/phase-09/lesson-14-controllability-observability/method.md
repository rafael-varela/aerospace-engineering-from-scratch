# P09.14  Method

The route to each answer, without the final numbers. Pair with `hints.md`; `worked-solutions.md` has the numbers for s = 15.

## 1. One axis

(a) B = (0, 1/J)ᵀ: AB = (1/J, 0)ᵀ, and W_c is an anti-diagonal matrix of rank 2. B = (1, 0)ᵀ: AB = 0, so W_c has rank 1 and the rate cannot be changed. (b) C = [1, 0]: CA = [0, 1], rank 2. C = [0, 1]: CA = 0, rank 1, and the null space is the angle direction (1, 0). (c) [λI − A, B] at λ = 0 is [−A, B]: rank 2 for B = (0, 1/J)ᵀ, rank 1 for B = (1, 0)ᵀ. The dual for C. (d) The four calls.

## 2. Derivations

(a) x(T) = ∫e^{A(T−τ)}Bu dτ = Σ_k A^kB ∫(T−τ)^k/k! u dτ, and by Cayley–Hamilton only k < n survive as independent terms. Conversely, the Gramian argument of (d) reaches every vector in that span. (b) w ẋ = wAx + wBu = λ wx, so wx(t) = e^{λt}wx(0) whatever the input: that component of the state is beyond the actuators' reach. (c) dF/dt = AF + FAᵀ, integrated from 0 to T: F(T) − F(0) = AW + WAᵀ, with F(0) = BBᵀ. For stable A, F(T) → 0 as T → ∞. (d) Substitute (6) into the convolution integral: ∫e^{A(T−t)}BBᵀe^{Aᵀ(T−t)}dt = W(T), so x(T) = e^{AT}x₀ + W W⁻¹d = x_f. For u + v with ∫e^{A(T−t)}Bv dt = 0, the cross term is ∫uᵀv dt = λᵀ∫e^{A(T−t)}Bv dt = 0, so the energy only adds. (e) e^{At}B = (t, 1)ᵀ/J; integrate the outer product to get W(T); then apply (6) with d = (Δθ, 0). For the bang-bang slew, half the angle comes in each half: τ = 4JΔθ/T² and the energy is τ²T.

## 3. Pyramid angle

Write W(β) and evaluate `eigvalsh(W @ W.T)` and `eigvalsh(ctrb_gramian_finite(A, B, 20.0)).min()` with and without wheel 0. Sweep β on a fine grid and take the maximiser of the smallest eigenvalue after a failure. A shallow pyramid puts most of each wheel's axis along z, which is good for z but weak in x and y. Losing a wheel then leaves a direction with very little authority.

## 4. Rendezvous

(a), (b) Call `numerical_rank(ctrb(A, B))` and `pbh_uncontrollable_modes(A, B)` for each B. For the radial-only case, find the left null vector w of W_c (the last right singular vector of W_cᵀ) and check that wA = 0 and wB = 0. Interpret it from the y equation of CW, which radial thrust does not enter: it is a conserved combination of x and ẏ that sets the along-track drift rate (the orbital period). (c) `obsv`, `pbh_unobservable_modes` and `unobservable_subspace` for each C. A pure along-track offset (x, ẋ, ẏ) = 0 with y ≠ 0 is an equilibrium of CW. A radial-only measurement never sees it.

## 5. Judgment

(a) Bisect on ε in [10⁻²⁰, 1] with `is_controllable`; evaluate λ_min W(20 s) at the two ε. Recommend the Gramian (or the smallest singular value of W_c, or the torque authority) as the health metric. (b) "Unobservable" makes a rigid design decision (cannot estimate at all); "poorly observable" means a large but finite estimation error. Compute the observability Gramian, or an estimator error covariance with the real noise levels (a Kalman filter, beyond this phase), instead of a rank. (c) Torque authority in the weakest direction (1/3 of nominal) and the momentum capacity in that direction; also the larger wheel speeds, friction and the loss of further redundancy.

## 6. Three-axis slew

(a) `min_energy_input(A, B, 0, x_f, 60)` for each `wheel_axes` configuration. (b) Evaluate u on 601 points and take the per-wheel maximum of |u|. (c) Compute J Δθ/|J Δθ| and its dot product with each wheel axis. The failed wheel's alignment with the required torque decides the cost. (d) From u(t) = Bᵀe^{Aᵀ(T−t)}λ with B = [[0], [J⁻¹W]], u = WᵀJ⁻¹ times a vector, so u lies in the row space of W. The unique solution of Wu = τ in that space is the minimum-norm one, Wᵀ(WWᵀ)⁻¹τ = W⁺τ. Check `u(0)` against `np.linalg.pinv(W) @ (W @ u(0))`.
