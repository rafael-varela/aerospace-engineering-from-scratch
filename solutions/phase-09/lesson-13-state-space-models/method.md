# P09.13  Method

The route to each answer, without the final numbers. Pair with `hints.md`; `worked-solutions.md` has the numbers for s = 15.

## 1. Second-order sensor

(a) Controllable: A = [[−2ζω_n, −ω_n²], [1, 0]], B = (1, 0)ᵀ, C = [0, ω_n²], D = 0. Observable: the transposes, A_o = Aᵀ, B_o = Cᵀ, C_o = Bᵀ. (b) G(jω_n) = ω_n²/(−ω_n² + 2jζω_n² + ω_n²) = 1/(2jζ), so the magnitude is 1/(2ζ) and the phase −90°. Check with `ss_eval(A, B, C, D, 1j*wn)` for both forms. (c) Physical x_p = (y, ẏ), with A_p = [[0, 1], [−ω_n², −2ζω_n]], B_p = (0, ω_n²)ᵀ, C_p = [1, 0]. The companion state is (ẇ, w) with y = ω_n²w, so y = ω_n² x_c,2 and ẏ = ω_n² x_c,1. Hence x_p = T⁻¹x_c with T⁻¹ = [[0, ω_n²], [ω_n², 0]]; verify with `similarity`.

## 2. Derivations

(a) C T(sI − T⁻¹AT)⁻¹T⁻¹B = C T T⁻¹(sI − A)⁻¹T T⁻¹B. The eigenvalues are the roots of det(sI − T⁻¹AT) = det(T⁻¹) det(sI − A) det(T). (b) w⁽ⁿ⁾ + Σa_k w⁽ⁿ⁻ᵏ⁾ = u and y = Σ_{k=0}^n b_k w⁽ⁿ⁻ᵏ⁾. Substitute w⁽ⁿ⁾ = u − Σa_k w⁽ⁿ⁻ᵏ⁾ into the k = 0 term: y = b₀u + Σ(b_k − b₀a_k)w⁽ⁿ⁻ᵏ⁾. (c) G = Gᵀ = (C(sI − A)⁻¹B)ᵀ + D = Bᵀ(sI − Aᵀ)⁻¹Cᵀ + D. This works only because a 1 × 1 matrix equals its transpose; for MIMO the transpose is the transfer matrix of the dual system, not the same one. (d) (sI − A)ΣN_k s^{n−k} = Σ c_k s^{n−k} I. Matching sⁿ gives N₁ = I; s^{n−k} gives N_{k+1} − AN_k = c_k I; s⁰ gives −AN_n = c_n I, which says that the characteristic polynomial annihilates A (Cayley–Hamilton). (e) As in §2.1: θ/δ = M_δ(s − Z_α)/(sΔ) is type 1 and rolls off as s⁻² (relative degree 2); h/δ = −VM_δZ_α/(s²Δ) is type 2 and rolls off as s⁻⁴.

## 3. Hidden mode

(a) From `tf_to_controllable`; the eigenvalues come from `np.linalg.eigvals`. (b) Two `simulate` calls with a unit step; take the largest absolute difference and the last value. (c) `simulate(..., np.zeros_like(t), dt, x0)`. The output decays as e^{−2t} only, because y = x₁ + x₂ and the e^{−t} mode's contribution to C cancels. (d) Use `np.linalg.matrix_rank`. This realisation is controllable but not observable (P09.14): the mode at −1 is hidden from the output.

## 4. Flexible array

(a) x = (θ₁, θ₂, ω₁, ω₂), with A rows [0 0 1 0], [0 0 0 1], [−k/J₁, k/J₁, −c/J₁, c/J₁], [k/J₂, −k/J₂, c/J₂, −c/J₂]; B = (0, 0, 1/J₁, 0)ᵀ, C = [1, 0, 0, 0]. (b) `ss_to_tf`, then `np.roots` of den and num. Compare the poles with ±j√(k(J₁ + J₂)/(J₁J₂)) (plus a small damping) and the zeros with ±j√(k/J₂). (c) `simulate` with u = 0.01 and dt = 0.05; compare with `scipy.signal.step` scaled by 0.01 and with u t²/(2(J₁ + J₂)). (d) `np.linalg.eig(A)` and `np.linalg.cond(V)`. The double zero eigenvalue (the rigid body) is defective, so there is no eigenvector basis and the modal form does not exist.

## 5. Judgment

(a) Quote the lesson's growth of the pole error from n = 6 to n = 10 (from rounding-level to order one) and the appearance of spurious complex pairs. Recommend modal coordinates (real 2 × 2 blocks) or physical or balanced coordinates for storage. (b) A TF is enough for frequency response, margins and responses from rest. A state model is needed for internal variables with limits (α), non-zero initial states (a gust upset), MIMO, and observer and state-feedback design. (c) The held elevator gives a quadratically growing climb in the model. With airspeed decaying, a real aircraft would slow down, pitch less and possibly stall. The missing state is airspeed u; the derivatives are X_u, Z_u (and M_u). (d) Rows 1 and 2 gain Z_α w_g and M_α w_g, so B_g = (Z_α, M_α, 0, 0)ᵀ. Repeat the elimination of §2.1 with the new input. You should find that α + w_g, the aerodynamic angle, carries the factor s(s − M_q)/Δ, and γ = θ − α = −Z_α(α + w_g)/s. Check the result with `ss_to_tf(A, B_g, C, D)`. The poles are the same; the numerator gains a zero at s = M_q, which also cancels one integrator.

## 6. Faddeev–LeVerrier accuracy

(a) Build Q from `np.linalg.qr(rng.standard_normal((n, n)))` with `rng = np.random.default_rng(1315)`. Compute the exact coefficients with Python integers and the relative error of both methods. (b) Repeat with Q = I and explain the exact integer arithmetic. (c) Use `np.roots` on the Faddeev–LeVerrier coefficients for the diagonal case and take the largest root error. (d) Conclude that the coefficient route is ill-conditioned at high order, and keep high-order models in state space or in factored (pole-zero) form.
