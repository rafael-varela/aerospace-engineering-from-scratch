# P09.14  Worked solutions

Worked for **s = 15**. Every number below is produced by `solution.py` (the functions named in each problem) or by the probe scripts that call it. All are model predictions for a representative satellite, not measurements. Units: SI.

## 1. One axis

**(a)** B = (0, 1/J)ᵀ: AB = (1/J, 0)ᵀ, **W_c = [[0, 1/J], [1/J, 0]], rank 2** (controllable). B = (1, 0)ᵀ: AB = 0, **W_c = [[1, 0], [0, 0]], rank 1**: an actuator that pushes the angle cannot change the rate.

**(b)** C = [1, 0]: CA = [0, 1], **W_o = I, rank 2**. C = [0, 1]: CA = [0, 0], **W_o = [[0, 1], [0, 0]], rank 1**, unobservable direction **(1, 0)**: the angle.

**(c)** At λ = 0, [−A, B] = [[0, −1, 0], [0, 0, 1/J]] has rank 2 (controllable); [[0, −1, 1], [0, 0, 0]] has rank 1 (not). [−A; C] = [[0, −1], [0, 0], [1, 0]] has rank 2 (observable); [[0, −1], [0, 0], [0, 1]] has rank 1 (not). Same verdicts.

**(d)** `ctrb`, `obsv`, `pbh_uncontrollable_modes` and `pbh_unobservable_modes` return these ranks, and [0] for the two failing cases.

## 2. Derivations

**(a)** From x₀ = 0, x(T) = ∫₀ᵀe^{A(T−τ)}Bu(τ)dτ = Σ_{k≥0} AᵏB ∫₀ᵀ(T − τ)ᵏ/k! u(τ)dτ. By Cayley–Hamilton, each Aᵏ with k ≥ n is a combination of I, …, A^{n−1}, so x(T) is in the column space of W_c. Conversely, if rank W_c = n then W(T) > 0 (otherwise some v ≠ 0 would have vᵀe^{At}B ≡ 0, and differentiating at t = 0 gives vᵀAᵏB = 0 for all k, i.e. vᵀW_c = 0), and (6) of the lesson reaches any state.

**(b)** If wA = λw and wB = 0, then (d/dt)(wx) = wAx + wBu = λ(wx), so **wx(t) = e^{λt}wx(0) for every input**. The state can never be driven to a target with a different wx: not controllable. Equivalently, w[λI − A, B] = 0, so the rank is below n.

**(c)** With F(t) = e^{At}BBᵀe^{Aᵀt}, dF/dt = AF + FAᵀ. Integrating, F(T) − F(0) = AW(T) + W(T)Aᵀ, i.e. **AW(T) + W(T)Aᵀ + BBᵀ = e^{AT}BBᵀe^{AᵀT}**. For stable A the right side tends to 0 and W(T) to W, giving **AW + WAᵀ + BBᵀ = 0**.

**(d)** Substituting (6): x(T) = e^{AT}x₀ + ∫e^{A(T−t)}BBᵀe^{Aᵀ(T−t)}dt · W⁻¹d = e^{AT}x₀ + W W⁻¹d = x_f. Any other input reaching x_f is u + v with ∫e^{A(T−t)}Bv dt = 0, and ∫uᵀv dt = dᵀW⁻¹∫e^{A(T−t)}Bv dt = 0, so ∫|u + v|² = ∫|u|² + ∫|v|² ≥ ∫|u|². The energy is ∫|u|² = dᵀW⁻¹WW⁻¹d = **dᵀW⁻¹d**.

**(e)** e^{At}B = (t, 1)ᵀ/J, so **W(T) = (1/J²)[[T³/3, T²/2], [T²/2, T]]**, det = T⁴/(12J⁴). For d = (Δθ, 0): W⁻¹d = (12J²/T⁴)(TΔθ, −T²Δθ/2). Then u(t) = Bᵀe^{Aᵀ(T−t)}W⁻¹d = (1/J)[(T − t), 1]W⁻¹d = **(6JΔθ/T²)(1 − 2t/T)**, with energy **12J²Δθ²/T³**. Bang-bang: each half covers Δθ/2 with acceleration a, so a(T/2)² = Δθ, a = 4Δθ/T², τ = 4JΔθ/T², and the energy is τ²T = **16J²Δθ²/T³ = (4/3) × minimum**.

## 3. Pyramid angle (s = 15: β = 47.5°)

| β (deg) | eig WWᵀ, all wheels | λ_min W(20 s), all | eig WWᵀ, wheel 0 failed | λ_min W(20 s), failed |
|---|---|---|---|---|
| 47.50 | 1.0872, 1.0872, 1.8257 | 0.2162 | 0.3095, 1.0872, 1.6034 | 0.1037 |
| 54.74 | 1.3333 (×3) | 0.2652 | 0.3333, 1.3333, 1.3333 | 0.1191 |
| 30.00 | 0.5, 0.5, 3.0 | 0.0994 | 0.1603, 0.5, 2.3397 | 0.0507 |
| 70.00 | 0.4679, 1.766, 1.766 | 0.2585 | 0.1998, 1.0342, 1.766 | 0.0950 |

**(b)** The sweep maximises the smallest eigenvalue of WWᵀ after one failure at **β = 54.75°** (value 0.3333), on the 0.05° grid of the sweep: the isotropic angle, cos²β = 1/3 (54.74°).

**(c)** At β = 30° most of each axis points along z (z authority 3.0, x and y only 0.5), so z slews are fast. After a failure the weakest direction has only 0.16 of a wheel's authority, half of what the isotropic pyramid keeps.

## 4. Rendezvous with a failed thruster

**(a)** Along-track only: **rank W_c = 4, controllable**; PBH reports no bad eigenvalue.

**(b)** Radial only: **rank 3**; PBH reports **λ = 0**. The left null vector is w = (2n, 0, 0, 1), with wA = 0 and wB = 0: the combination **ẏ + 2nx is conserved** whatever the radial thrust does. This is clear from the CW equation ÿ = −2nẋ, which the radial thrust does not enter. That combination sets the along-track drift rate, which is the difference in orbital period. A radial impulse does not change the orbital energy to first order, so a radial thruster alone cannot stop a chaser from drifting.

**(c)** Along-track measurement only: **rank W_o = 4, observable**. Radial measurement only: **rank 3**, PBH λ = 0, unobservable direction **(0, ±1, 0, 0)**: a pure along-track offset. A chaser at (0, y₀) at rest is on the target's orbit and stays at constant y (`P02.05`): it never produces a radial signal.

## 5. Judgment

**(a)** With the z wheel at fraction ε of its torque, `is_controllable` says False only for **ε below about 7.5 × 10⁻¹⁰**. The smallest eigenvalue of W(20 s) is **5.5 × 10⁻⁷** at ε = 10⁻³ and **5.5 × 10⁻¹³** at ε = 10⁻⁶, against 0.199 nominal: it falls as ε². Both cases are "controllable", but the energy to move the z axis grows by 10⁶ and 10¹². An operator should watch a continuous health number: the smallest singular value of W (the torque authority), or the smallest Gramian eigenvalue, never the rank.

**(b)** Too pessimistic for the design, but the right warning. "Unobservable" means no estimator can recover the boresight angle at all. "Poorly observable" means an estimate with a larger error. Compute the observability Gramian with the real accuracy as a weight (or an estimator's error covariance with the real noise levels) and design for the resulting error, rather than giving up the axis or trusting it as much as the others.

**(c)** It does not guarantee (1) **authority**: the weakest direction keeps 1/3 of a nominal wheel's torque (WWᵀ eigenvalue 0.3333 against 1.3333), so the torque per unit of wheel effort in that direction is √(1/3)/√(4/3) = 1/2 of nominal, and a slew there takes √2 times longer at the same wheel torque; and (2) **momentum capacity and redundancy**: the remaining wheels spin faster to store the same momentum, saturate sooner, and one more failure loses an axis (rank 4).

## 6. Three-axis slew (T = 60 s, to (5°, −3°, 8°))

**(a)**

| configuration | energy (N² m² s) |
|---|---|
| pyramid, all wheels | 1.5244 × 10⁻⁵ |
| pyramid, wheel 0 failed | 4.9940 × 10⁻⁵ |
| pyramid, wheel 2 failed | 1.5477 × 10⁻⁵ |
| orthogonal | 2.0325 × 10⁻⁵ |

**(b)** Peak torque per wheel (mN m): full pyramid **0.659, 0.035, 0.054, 0.569**; wheel 0 failed **0.694, 0.713, 1.228** (wheels 1, 2, 3); wheel 2 failed **0.713, 0.019, 0.516** (wheels 0, 1, 3); orthogonal **0.582, 0.436, 0.698**. The largest is 1.23 mN m, **well within** the 10 mN m limit.

**(c)** The torque the slew needs points along J Δθ/|J Δθ| = (0.577, −0.433, 0.693). Its dot products with the four pyramid axes are **0.871, 0.046, −0.071, 0.753**. Wheel 2 is nearly perpendicular to the job and was hardly used (0.054 mN m), so losing it changes almost nothing (energy +1.5 %). Wheel 0 carries most of it, so the other three must combine inefficiently to replace it: the energy rises by a factor 3.3.

**(d)** u(t) = Bᵀe^{Aᵀ(T−t)}λ with B = [[0], [J⁻¹W]] gives u(t) = WᵀJ⁻¹ times a vector: u lies in the row space of W at every instant. Among all u with Wu = τ, the one in the row space of W is the minimum-norm solution Wᵀ(WWᵀ)⁻¹τ = W⁺τ. Numerically, at t = 0 the full pyramid has u(0) = (6.586 × 10⁻⁴, 3.51 × 10⁻⁵, −5.40 × 10⁻⁵, 5.695 × 10⁻⁴) N m and body torque Wu(0) = (5.818 × 10⁻⁴, −4.363 × 10⁻⁴, 6.981 × 10⁻⁴) N m; `np.linalg.pinv(W)` applied to that torque returns u(0) to **10⁻¹⁹**.

**Common mistakes.** Splitting the body torque equally among wheels (not minimum energy); forgetting that the inertia weights the required torque direction (J Δθ, not Δθ).

## Check-yourself answers (lesson.md)

1. (A, B) is controllable if any state can be driven to any other in any finite time. Test: rank [B, AB, …, A^{n−1}B] = n.
2. Which eigenvalue (mode) is uncontrollable: rank [λI − A, B] < n identifies λ, and the left null vector says which combination of states is stuck.
3. W_o(A, C) = W_c(Aᵀ, Cᵀ)ᵀ, so the ranks are equal; equivalently, rank [λI − A; C] = rank [λIᵀ − Aᵀ, Cᵀ].
4. rank W_c = 4 of 6; the angle and rate about the failed wheel's axis (a double integrator, λ = 0) cannot be influenced.
5. Gyros measure rates, and the angle is the integral of the rate plus an unknown constant. Every attitude offset with zero rate gives identical gyro outputs: the angles are unobservable (rank 3 of 6).
6. AW + WAᵀ + BBᵀ = 0. It requires A stable so that ∫₀^∞ converges; the free rigid body has eigenvalues at 0, and the integral grows like T³.
7. 6JΔθ/T² = 6 × 3 × 0.17453/400 = 7.85 mN m (the torque ramps linearly from +7.85 to −7.85 mN m).
8. The rank test still says controllable (rank 6); the smallest Gramian eigenvalue (or the smallest singular value of W) falls by about 10¹² (as torque squared) and shows that the direction is practically lost.
