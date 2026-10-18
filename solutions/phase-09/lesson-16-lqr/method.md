# P09.16  Method

The route to each answer, without the final numbers. Pair with `hints.md`; `worked-solutions.md` has the numbers for s = 15.

## 1. Scalar regulator

(a) With A = a, B = b, Q = q, R = r the CARE is 2aP − (b²/r)P² + q = 0, so P = (r/b²)(a ± √(a² + b²q/r)). The closed-loop pole a − b²P/r = ∓√(a² + b²q/r) is negative for the + root: that is the stabilising P. K = bP/r. (b) r → ∞: √(a² + b²q/r) → |a|, so the pole → −|a| and K → 2a/b for a > 0 (the stable mirror image of the unstable pole: the minimum-energy stabiliser). r → 0: the pole → −b√(q/r) and K → √(q/r), both unbounded. (c) Build 1 × 1 arrays and call `care_hamiltonian`, `care_newton_kleinman` and `scipy.linalg.solve_continuous_are`.

## 2. Hamiltonian and margins

(a) Form the Hamiltonian function ℋ = xᵀQx + uᵀRu + λᵀ(Ax + Bu); ∂ℋ/∂u = 0 gives u = −½R⁻¹Bᵀλ; the costate equation λ̇ = −∂ℋ/∂x = −2Qx − Aᵀλ. Divide by 2 and stack with ẋ. (b) If (x, λ/2) = (X₁c, X₂c) then λ/2 = X₂X₁⁻¹x = Px. Differentiate λ/2 = Px along the trajectory, substitute both equations and require it for every x: the CARE. (c) J₂HJ₂⁻¹ = [[−Aᵀ, Q], [BR⁻¹Bᵀ, A]] = −Hᵀ (check the blocks), so H is similar to −Hᵀ and its eigenvalues are those of −H. (d) The algebra yields (I + L(−jω))ᵀR(I + L(jω)) = R + Gᴴ(jω)QG(jω) ⪰ R; for one input divide by R > 0. Then use the geometry of the disk |1 + L| < 1: real-axis crossings at ≤ −2, and the intersection angles of the unit circles about 0 and −1 (±120°) at |L| = 1.

## 3. Roll attitude

(a) Q = diag(1/φ_max², 1/p_max²) and R = 1/δ_max² in radians. (b) `lqr(A, B, Q, R)`; ω_n = |λ|, ζ = −Re λ/|λ| for each pole (both real here). (c) `simulate(A, B, K, (10°, 0), 5, 0.001)`, max |u|; compare K₁ × 10° (in radians, then degrees). (d) `loop_transfer` on a log grid from 10⁻² to 10³ rad/s, `phase_margin`, min |1 + L|. (e) Repeat with R/10; tabulate K, poles, peak δ, PM, ω_c and settling.

## 4. Your satellite

(a) `satellite_model(5 + s)`, `bryson_weights([5°, 0.5°/s], [0.02])`, `lqr`; compare with (9). (b) `simulate` from (5°, 0) for 300 s at 0.05 s: max |u|, the last sample outside ±2 % of θ₀ plus one step, and J × max |ω|. (c) `simulate` from (15°, 0) with and without `u_max = 0.02`; count samples at the limit; settling time as before. (d) `care_newton_kleinman` returns the count; compare with `care_hamiltonian` by the max relative difference.

## 5. Judgment

(a) Name effects that add dynamics or uncertainty inside the loop: estimators, delays, actuator dynamics, flexible modes, multiple inputs with coupling, sensor noise. (b) The LQG case: an observer in the loop breaks the return-difference inequality at the plant input (Doyle 1978: "Guaranteed margins for LQG regulators — there are none"); check margins on the actual loop, or use loop-transfer recovery. (c) Delay margin = PM (rad)/ω_c; compare the phase lost by 0.25 s at ω_c with the PM. (d) Compare the mode frequency with ω_c and evaluate |L| at the mode frequency: a lightly damped mode multiplies the loop gain near its frequency by about 1/(2ζ).

## 6. Tightness and cost curve

(a) From Problem 1, K → 2a/b as r → ∞; the loop with kK has pole a(1 − 2k): stable iff k > ½. So for this family the lower gain margin approaches ½ exactly. (b) `cost_of_gain(A, B, Q, R, k K, x0)` for the seven k. (c) Recompute K for each new Q (the gain changes with Q) and evaluate J(0.5K) and J(2K); then keep Q and change x₀.
