# P03.13  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and given to the figures shown.

## 1. A two-element bar by hand (s = 15: P = 25 kN)

**Assumptions.** Linear elastic, uniaxial, constant EA = 200 × 10⁹ × 200 × 10⁻⁶ = 40 MN; h = 0.5 m, EA/h = 80 MN/m.

(a) k₁ = k₂ = 80 × 10⁶ [[1, −1], [−1, 1]] N/m. Assembled: K = 80 × 10⁶ [[1, −1, 0], [−1, 2, −1], [0, −1, 1]], f = [0, 0, 25 000] N. With u₀ = 0 the free equations are 80 × 10⁶ (2u₁ − u₂) = 0 and 80 × 10⁶ (−u₁ + u₂) = 25 000, so u₂ = 2u₁ and u₁ = 25 000/80 × 10⁶ = **0.3125 mm**, u₂ = **0.625 mm**.

(b) σ₁ = σ₂ = E(0.3125 × 10⁻³)/0.5 = **125 MPa** = P/A ✓. Reaction r₀ = (Ku)₀ − f₀ = 80 × 10⁶(0 − 0.3125 × 10⁻³) = **−25 kN** (the wall pulls back). u(L) = PL/(EA) = 25 000/40 × 10⁶ = 0.625 mm ✓.

(c) Consistent load for uniform q on an element: qh/2 at each node, so f = [1250, 2500, 1250] N (the middle node collects from both elements). Solving: u₁ = **0.046875 mm**, u₂ = **0.0625 mm**; exact tip qL²/(2EA) = 5000/(80 × 10⁶) = 0.0625 mm ✓, and u₁ exact too (q(Lx − x²/2)/(EA) at 0.5 m = 0.046875 mm). The agreement is exact because a uniform bar's Green's function is piecewise linear with its kink at the source node, so the finite-element solution of any load is exact at the nodes (lesson §2.3).

**Check.** Σr = −Σf: −25 kN against 25 kN ✓; in (c) the reaction is −qL = −5 kN ✓.

**Common mistakes.** Giving the middle node qh/2 instead of qh; forgetting to strike *both* the row and the column of the fixed node.

## 2. The weak form and what follows

(a) Multiply −(EAu')' = b by v with v(0) = 0 and integrate: −∫(EAu')'v dx = ∫bv dx. Parts: −[EAu'v]₀ᴸ + ∫EAu'v' dx = ∫bv dx. At x = 0, v = 0; at x = L, EAu' = P. So ∫EAu'v' dx = ∫bv dx + Pv(L) for all admissible v. **Essential**: u(0) = 0, imposed on the trial space (and v(0) = 0 on the test space). **Natural**: EAu'(L) = P, which appears as the boundary term and is satisfied weakly by the solution.

(b) In ξ ∈ [−1, 1]: x = x_mid + (h/2)ξ, dx = (h/2)dξ, N₁ = (1 − ξ)/2, N₂ = (1 + ξ)/2, dN₁/dx = −1/h, dN₂/dx = +1/h. k_ij = EA ∫N_i'N_j' dx = EA (h/2)(2)(±1/h²) = ±EA/h: k = (EA/h)[[1, −1], [−1, 1]]. m_ij = ρA(h/2)∫N_iN_j dξ; ∫N₁² dξ = ∫(1 − ξ)²/4 dξ = 2/3, ∫N₁N₂ dξ = ∫(1 − ξ²)/4 dξ = 1/3; so m = ρA(h/2)[[2/3, 1/3], [1/3, 2/3]] = (ρAh/6)[[2, 1], [1, 2]].

(c) K_ij = ∫EAN_i'N_j' dx = K_ji. For any v, vᵀKv = Σ_ij v_iv_j∫EAN_i'N_j' dx = ∫EA(Σ_i v_iN_i')² dx = ∫EA(v_h')² dx ≥ 0. It is zero iff v_h' ≡ 0 iff v_h is constant iff v = c(1, …, 1): the ones vector spans the null space (K·1 = 0 since the hat functions sum to one and the sum has zero derivative). Fix node m: the remaining space contains no non-zero constant vector (its m-th entry would have to be zero), so vᵀKv > 0 for every non-zero admissible v: positive definite.

(d) ∂Π/∂u_i = Σ_j K_ij u_j − f_i, so ∇Π = Ku − f = 0 is (3). The Hessian is K, positive definite on the free dofs, so the stationary point is the unique minimum. At it, Π = ½uᵀKu − fᵀu = ½uᵀ(Ku) − uᵀKu = −½uᵀKu, and fᵀu = uᵀKu = 2U, i.e. U = ½fᵀu (Clapeyron).

**Common mistakes.** Dropping the boundary term at L (which loses the tip load entirely); integrating N_iN_j in x without the Jacobian h/2.

## 3. Thermal strain (E = 70 GPa, A = 600 mm², α = 23 × 10⁻⁶ K⁻¹, ΔT = 50 K, L = 2 m)

(a) The weak form of equilibrium is ∫σAv' dx = ∫bv dx + Pv(L). With σ = E(u' − αΔT): ∫EAu'v' dx − ∫EAαΔT v' dx = ∫bv dx + Pv(L). The second term moves to the right as a load; for v = N_i on an element, N_i' = ∓1/h and the integral over the element is EAαΔT(h)(∓1/h) = ∓EAαΔT: f_e = EAαΔT[−1, +1]ᵀ. K is untouched. (EAαΔT = 70 × 10⁹ × 6 × 10⁻⁴ × 23 × 10⁻⁶ × 50 = 48.3 kN.)

(b) Both ends fixed, four elements: thermal loads cancel at every interior node (−48.3 + 48.3), so f = [−48.3, 0, 0, 0, 48.3] kN and the free equations give **u = 0** at every node. Mechanical strain u' = 0; total stress σ = E(0 − αΔT) = **−80.5 MPa** in every element; reactions **+48.3 kN at x = 0 and −48.3 kN at x = L** (the walls push in: compression). Check: σ = −EαΔT = −70 × 10⁹ × 1.15 × 10⁻³ = −80.5 MPa ✓; σA = −48.3 kN = the reaction ✓.

(c) Far end free: f has +48.3 kN at the free end and nothing elsewhere free; u grows linearly to u(L) = **2.300 mm** = αΔTL = 1.15 × 10⁻³ × 2 ✓; u' = αΔT, so σ = E(αΔT − αΔT) = **0** in every element: free thermal expansion is stress-free.

(d) Equilibrium with b = 0 and both ends fixed requires σA constant, so σ is **uniform**, and ∫u' dx = 0 gives σ = −Eα × (mean ΔT) = −Eα × 25 K = −40.25 MPa; u' = α(ΔT(x) − 25 K) is linear in x, so u is quadratic (u = α(25x²/2 − 25x), −0.2875 mm at mid-length). The four-element model, with each element's thermal load evaluated from the element's mean ΔT, gets the nodal displacements exactly (nodal exactness holds for any load of a uniform bar) and therefore the element stress E(Δu/h − αΔT_mean) = −40.25 MPa exactly in every element; the run confirms both to the printed digits.

**Common mistakes.** Reporting the mechanical stress Eu' (zero) instead of the total stress; applying the thermal load as a single tip force rather than element by element (right answer for a uniform ΔT, wrong for a varying one).

## 4. Your own rotor (s = 15: Ω = 395 rpm = 41.364 rad/s, m_tip = 3.0 kg)

(a) P = m_tipΩ²L = 3 × 41.364² × 5 = **25.665 kN**; spar load ρAΩ²L²/2 = 2700 × 6 × 10⁻⁴ × 1711.0 × 25/2 = **34.648 kN**; root stress (34.648 + 25.665)/6 × 10⁻⁴ = **100.52 MPa**; tip extension (ρΩ²/2E)(2L³/3) + PL/(EA) = **5.8052 mm**; reaction **−60.313 kN**.

(b) n = 8:

| quantity | FE (8 elements) | exact | status |
|---|---|---|---|
| tip displacement | 5.8052 mm | 5.8052 mm | exact (nodal exactness) |
| root element stress | 100.22 MPa | 100.52 at x = 0; 100.30 at x = h/2 | element average: 0.30 MPa low |
| reaction | −60.313 kN | −60.313 kN | exact (equilibrium) |

(c) Element centres at 0.3125 and 0.9375 m with stresses 100.22 and 98.42 MPa; linear extrapolation to x = 0 gives **101.12 MPa**, 0.6 % *above* the exact 100.52: worse than the element value (0.3 % low). The exact stress is a downward parabola with zero slope at the root; a straight line through two points on its descending side overshoots at x = 0. Extrapolation is exact only for a linear stress (the tapered strut without body force is closer to that case).

(d) Orders: L2 **1.99, 2.00, 2.00**; energy **0.99, 1.00, 1.00**.

**Check.** Σ loads = 60.313 kN = −reaction ✓; the root stress scales with Ω² between the lesson case and this one: 88.46 × (41.364/41.89)² × (34.648 + 25.665)/(35.53 + 17.55) ... simpler: both components computed directly ✓.

**Common mistakes.** Comparing the element stress with the root value and calling the 0.3 % a bug; extrapolating to improve a quantity whose exact slope is zero.

## 5. Consistent against lumped mass (spar, not spinning; exact 254.59, 763.76, 1272.94 Hz)

(a) Ratios ω_FE/ω_exact:

| n | consistent: mode 1, 2, 3 | lumped: mode 1, 2, 3 |
|---|---|---|
| 4 | 1.0064, 1.0583, 1.1535 | 0.9936, 0.9432, 0.8469 |
| 8 | 1.0016, 1.0145, 1.0405 | 0.9984, 0.9856, 0.9603 |
| 16 | 1.0004, 1.0036, 1.0101 | 0.9996, 0.9964, 0.9900 |
| 32 | 1.0001, 1.0009, 1.0025 | 0.9999, 0.9991, 0.9975 |

(b) Every consistent ratio is above 1 and falls with n; every lumped ratio is below 1 and rises. First-mode errors: consistent 0.0064, 0.0016, 0.0004, 0.0001 → **p = 2.0**; lumped 0.0064, 0.0016, 0.0004, 0.0001 → **p = 2.0**, and to this precision equal and opposite (the average of the two is correct to a higher order, which is the basis of "averaged" mass matrices).

(c) Mode 3 has 2.5 half-waves on the bar. Consistent mass reaches 1 % at n = 16 (1.0101, just over; n ≈ 17 strictly): about **6–7 elements per half-wave**. Lumped mass reaches 1 % at n = 16 (0.9900): the same count. Rule of thumb from either: ten elements per half-wavelength for a few tenths of a per cent.

(d) Ω = 6.67 Hz; harmonics up to 10Ω reach 66.7 Hz; the first axial frequency, 254.6 Hz, is 38Ω: no resonance. But the question is empty for a real blade: its flapping and lead-lag *bending* frequencies are of the order of 1–1.5 Ω for the first modes (the blade is a slender beam, hinged or soft at the root), and those are the ones that must be placed between the harmonics, with centrifugal stiffening raising them as Ω rises. The axial stiffness EA/L is hundreds of times larger than the bending stiffness, hence the factor of forty.

**Common mistakes.** Lumping by putting the whole element mass at one node; concluding from the axial result that the blade has no resonance problem.

## 6. The tapered strut and the limits of the bar (s = 15: c = 0.45; A₀ = 800 mm², L = 1.5 m, P = 50 kN)

Exact tip stress P/(A₀(1 − c)) = 50 000/(440 × 10⁻⁶) = **113.64 MPa**; exact tip displacement 1.7793 mm.

(a)

| n | last-element stress [MPa] | error vs tip | exact at last-element centre |
|---|---|---|---|
| 4 | 103.09 | −9.28 % | 103.09 |
| 8 | 108.11 | −4.86 % | 108.11 |
| 16 | 110.80 | −2.49 % | 110.80 |
| 32 | 112.20 | −1.26 % | 112.20 |
| 64 | 112.91 | −0.64 % | 112.91 |
| 128 | 113.27 | −0.32 % | 113.27 |

Within 2 %: **n = 32**; within 0.5 %: **n = 128**. The tip displacement is already within 0.2 % at n = 4 (1.7757 mm) and 0.01 % at n = 16.

(b) With no body force the discrete equilibrium of the last node is k_e(u_n − u_{n−1}) = P with k_e = EA(x_mid)/h, so E(u_n − u_{n−1})/h = P/A(x_mid): the element stress is **exactly** the exact stress at the element's mid-point (the table confirms it to the printed digits). The tip is h/2 further out, where A is smaller by A'h/2 = −cA₀h/(2L), so the relative error is ≈ −(c/(1 − c))(h/2L) = −(0.45/0.55)(1/(2n)): −10.2 % at n = 4, −0.32 % at n = 128, matching the table to the next order. First order in h, because the error is the distance to the tip times the stress gradient.

(c) Use the natural boundary condition: σ(L) = P/A(L) = 113.64 MPa, exact with **any** mesh, because the tip force is known. It relies on the tip load being a prescribed force (a natural boundary condition) rather than a prescribed displacement; in the latter case one would extrapolate the last two element stresses (exact for a linear σ, here 0.18 % low at n = 16) or compute the reaction. (Extrapolating the last two element stresses at n = 16 gives 113.43 MPa, 0.18 % low, against the 2.5 % of the raw element value.) Codes call this stress recovery.

(d) In order: (1) **No bending**: a pinned–pinned strut under a lateral component or an eccentric load deflects laterally and carries bending stress that can exceed the axial stress; measurable as a mid-length lateral deflection and a strain that differs between the two sides of the section. (2) **The end fixity**: clevis pins concentrate the load at the lug holes (stress concentration factor 2–3); measurable as local strain at the lug far above P/A. (3) **Uniaxial stress state**: the taper and the section shape (a tube or an I) make the stress non-uniform across the section near changes; measurable as circumferential variation of axial strain. The bar model's P/A and PL/(EA) survive as nominal values; nothing in it predicts the three effects above.

**Common mistakes.** Refining the mesh to recover a tip stress that the boundary condition already gives; reading the −9 % at n = 4 as a stiffness error (the displacement is within 0.2 %).

## Check-yourself answers

1. Strong: −(EAu')' = b, u(0) = 0, EAu'(L) = P. Weak: ∫EAu'v' dx = ∫bv dx + Pv(L) for all v with v(0) = 0. Essential: u(0) = 0 (built into the space); natural: EAu'(L) = P (appears as the boundary term).
2. N₁' = −1/h, N₂' = +1/h on the element; k_ij = EA∫N_i'N_j' dx = EA h (±1/h²) = ±EA/h, giving (EA/h)[[1, −1], [−1, 1]].
3. Because a rigid translation (all nodes moving equally) stores no energy: K·1 = 0, so K is singular with the ones vector as its null vector until a displacement is fixed.
4. Uniform b₀: f_e = b₀h [½, ½]. Linear b = b₁x/h (zero at the first node, b₁ at the second): f_e = b₁h [1/6, 1/3].
5. The Green's function of the uniform bar is piecewise linear with a kink at the source, so for a source at a node it lies in the trial space; Galerkin orthogonality then forces u_h to equal u at every node for any load. With A(x) varying, the Green's function is no longer piecewise linear and the nodal values converge at order 2 instead.
6. u_h is the best approximation of u in the energy norm (slope error) among piecewise-linear functions; the best piecewise-linear interpolant has slope error O(h) and value error O(h²), so the stress (slope) converges at order 1 and the displacement at order 2.
7. Restricting the motion to the trial space constrains the structure and raises every Rayleigh quotient, so the consistent-mass eigenvalues are upper bounds that fall as the space grows. Lumping the mass over-weights the inertia of the low modes and makes the frequencies lower bounds that rise with refinement.
8. No: the element reports P/A at its mid-point, which is exact there; the tip is h/2 away where the stress is higher, a first-order recovery error, not a solution error. Use the natural boundary condition σ(L) = P/A(L), or extrapolate from element centres, or grade the mesh toward the tip.
