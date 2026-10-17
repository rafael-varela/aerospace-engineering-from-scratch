# P06.09  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for **s = 15**; numerical results are from `solution.py` (and the probe scripts that call it) and are given to the figures shown. All "simulated" values are outputs of the model, not measurements.

## 1. A satellite with a boom (s = 15: m = 13 kg)

**Assumptions.** Point mass; light boom; the centre-of-mass shift ignored (as stated).

(a) |r|² = 2.5 m², rrᵀ = [[2.25, 0.75, 0], [0.75, 0.25, 0], [0, 0, 0]], so ΔI = 13(2.5E − rrᵀ) = **[[3.25, −9.75, 0], [−9.75, 29.25, 0], [0, 0, 32.5]] kg·m²**, and

I = **[[263.65, −9.75, 0], [−9.75, 237.55, 0], [0, 0, 167.9]] kg·m²**.

Only I₁₂ appears, because r lies in the 1–2 plane (r₃ = 0): products I₁₃ and I₂₃ need a component along axis 3.

(b) Principal moments **167.90, 234.31, 266.89 kg·m²**; axes (columns of R): (0, 0, 1) for 167.90; (−0.3154, −0.9490, 0) for 234.31; (0.9490, −0.3154, 0) for 266.89.

(c) From body axis 1: **90°** (minor axis = body axis 3), **71.62°** (intermediate), **18.38°** (major). The new intermediate axis lies in the 1–2 plane, 18.4° from body axis 2: the boom has turned the principal axes in that plane.

**Check.** The trace is unchanged by the rotation: 263.65 + 237.55 + 167.9 = 669.1 = 167.90 + 234.31 + 266.89.

**Common mistakes.** Adding m|r|² to every diagonal entry and forgetting −m rrᵀ; reading the eigenvectors as rows.

## 2. The geometry and the exact solution

(a) With τ = 0: ω · Iω̇ = −ω · (ω × Iω) = 0, and since I is symmetric ω · Iω̇ = d(½ω · Iω)/dt, so **T is constant**. (Iω) · Iω̇ = −(Iω) · (ω × Iω) = 0, i.e. d(½|Iω|²)/dt = 0: **H² is constant**. In the inertial frame, dH_i/dt = Cᵀ(d(Iω)/dt + ω × Iω) = Cᵀτ = 0: **H_i is constant**.

(b) H² − 2TI₂ = Σ Iₖ²ωₖ² − I₂ Σ Iₖωₖ² = Σ Iₖ(Iₖ − I₂)ωₖ² = **I₁(I₁ − I₂)ω₁² + I₃(I₃ − I₂)ω₃²** (the k = 2 term vanishes). Setting it to zero: **I₁(I₂ − I₁)ω₁² = I₃(I₃ − I₂)ω₃²**, i.e. ω₃ = ±√(I₁(I₂ − I₁)/(I₃(I₃ − I₂))) ω₁: two planes through the ω₂ axis.

(c) ω_i · H_i = (Cᵀω) · (CᵀIω) = ω · Iω = 2T, constant, and H_i is constant, so the projection of ω_i on Ĥ is **2T/|H|**: the tip moves in a fixed plane perpendicular to H at that distance.

(d) With (6) and τ = λt, (3) becomes I₁A₁(−sn dn)λ = (I₂ − I₃)A₂A₃ sn dn, I₂A₂(cn dn)λ = (I₃ − I₁)A₃A₁ cn dn, I₃A₃(−m sn cn)λ = (I₁ − I₂)A₁A₂ sn cn. Hence

I₁A₁λ = (I₃ − I₂)A₂A₃,  I₂A₂λ = (I₃ − I₁)A₃A₁,  mI₃A₃λ = (I₂ − I₁)A₁A₂.

At τ = 0 (sn = 0, cn = dn = 1): 2T = I₁A₁² + I₃A₃², H² = I₁²A₁² + I₃²A₃². Solving the 2 × 2 system: A₁² = (2TI₃ − H²)/(I₁(I₃ − I₁)) = a/(I₁(I₃ − I₁)) and A₃² = (H² − 2TI₁)/(I₃(I₃ − I₁)) = b/(I₃(I₃ − I₁)). Multiplying the first two relations, (I₁A₁λ)(I₂A₂λ) = (I₃ − I₂)(I₃ − I₁)A₁A₂A₃², so **λ² = (I₃ − I₂)(I₃ − I₁)A₃²/(I₁I₂) = (I₃ − I₂)b/(I₁I₂I₃)**. Then A₂ = I₁A₁λ/((I₃ − I₂)A₃) gives **A₂² = a/(I₂(I₃ − I₂))**, and the third relation gives **m = (I₂ − I₁)A₁A₂/(I₃A₃λ) = (I₂ − I₁)a/((I₃ − I₂)b)**, which is (7).

(e) m → 0 means a → 0: K → π/2, so 4K/λ → 2π/λ with λ² → (I₃ − I₂)b/(I₁I₂I₃) and b → H² − 2TI₁ ≈ I₃(I₃ − I₁)ω₃² for spin ω₃: λ² → (I₃ − I₂)(I₃ − I₁)ω₃²/(I₁I₂), the square of `P04.04`'s wobble frequency. m → 1 means k′ → 0, where K(m) = ln(4/k′) + O(k′² ln k′); simultaneously the state approaches pure ω₂ spin, so b → I₂(I₂ − I₁)ω₂² and λ² → (I₃ − I₂)(I₂ − I₁)ω₂²/(I₁I₃) = σ². Hence 2K/λ → **(2/σ) ln(4/k′)**.

**Checks.** (d) for m = 0 gives the sin/cos wobble; at m = 1, sn = tanh: the rate approaches −ω₂ only as t → ∞ (separatrix).

**Common mistakes.** Using d dn/dτ = −sn cn without the factor m; mixing up which axis carries dn on the two sides of the separatrix.

## 3. Approaching a tumbling satellite (s = 15: ω = (0.02, 0.045, 0.05) rad/s)

(a) T = ½(1200 × 0.0004 + 1900 × 0.002025 + 2300 × 0.0025) = **5.039 J**; |H| = |(24, 85.5, 115)| = **145.30 N·m·s**. H² − 2TI₂ > 0: the polhode circles the **major** axis (axis 3).

(b) m = **0.4011**; period **271.2 s** (4.5 min).

(c) Body axis **3** stays nearest Ĥ; over one period its angle from H swings between **19.08° and 43.00°**. (The simulation conserves T to 2 × 10⁻¹², |H| to 1 × 10⁻¹² and H_i to 1.4 × 10⁻¹¹.)

(d) H is fixed in space for a torque-free body, so the servicer can hold station along it without chasing; the target then rotates about that line. But its axis 3 cones 19°–43° around it with the 271 s period, so a docking fixture on axis 3 would be presented at a varying angle: the servicer must match the tumble or detumble the target first.

**Check.** At t = 0 the angle is arccos(H₃/|H|) = arccos(115/145.30) = 37.7°, within the range.

**Common mistakes.** Taking the instantaneous ω as the direction to approach along (ω moves; H does not).

## 4. Your own T-handle (s = 15: 3.5 rev/s, ε = 10⁻⁴)

(a) spin = 21.991 rad/s; polhode **'minor'**; m = **0.99999998638** (1 − m = 1.362 × 10⁻⁸); λ = **12.5712 rad/s**; period **3.3226 s**; flip interval **1.6613 s**; σ = **12.5712 s⁻¹** (equal to λ to eight figures, so close is the state to the separatrix).

(b) At rtol 10⁻¹¹: flips at **0.8306, 2.4919, 4.1532, 5.8145, 7.4758, 9.1371 s**; mean interval **1.661289 s** against 2K/λ = 1.661289 s; largest deviation from the elliptic solution **6.2 × 10⁻⁷** of the spin; drifts T **3.6 × 10⁻¹¹**, |H| **1.8 × 10⁻¹¹**, H_i **7.2 × 10⁻¹¹**.

(c) At rtol 10⁻⁶: the drifts grow to **1.2 × 10⁻⁵, 6.1 × 10⁻⁶ and 1.2 × 10⁻⁵**, and the deviation from the elliptic solution to 1.9 × 10⁻⁵, but the flip times are **the same to 10⁻⁴ s** (mean interval 1.661289 s). A solver error acts like a small change of the initial state; by (9) that moves the flip time by only (1/σ) times the relative change of k′.

(d) 60/1.6613 = **36.1 flips per minute**.

**Common mistakes.** Measuring flips from |ω₂| maxima (twice per flip) instead of sign changes.

## 5. The box satellite of P04.04

(a) σ = **0.16408 s⁻¹**; flip interval **172.20 s**; first reversal **86.10 s**.

(b) The linear estimate ln(0.5/10⁻⁶)/σ = **79.97 s**, 6.1 s early. It assumes exponential growth all the way to the size of the spin, while the true turn-over slows as the disturbance becomes large; the closed form's (1/σ) ln(4/k′) has the extra factor 4 and the exact k′, and gives 86.10 s.

(c) `simulate` over 200 s with 20 001 samples: first reversal at **86.1005 s**, equal to K/λ to the digits shown.

(d) The first reversal comes ln(1000)/σ = **42.1 s** later, at **128.2 s**: a thousand-fold better initial spin buys 42 s. Intermediate-axis spin cannot be made stable by precision; it needs active control (or a different axis).

**Check.** (9) per decade of disturbance for the first reversal: ln 10/σ = 14.03 s; three decades: 42.1 s. ✓

**Common mistakes.** Using the flip interval for the first reversal (the spin starts at the top of the sn curve, so the first zero comes after a quarter period).

## 6. Judgment

(a) Two common causes: computing H_i with C instead of Cᵀ (or a quaternion convention mismatch), which makes a conserved vector look as if it turns; and a quaternion that is not renormalised before use, or q̇ written as ½(0, ω) ⊗ q with body rates. T is a body-frame scalar, so it is blind to both.

(b) The flip time depends on the logarithm of the initial disturbance, which no video measures to better than tens of per cent, so its phase is uncertain by a large fraction of a flip interval; the interval itself depends only on the spin, the moments and the logarithm of the (roughly constant) energy offset from the separatrix, and successive flips are equally spaced.

(c) The 301 axis is the intermediate one and is unstable, but σ = 0.036 × spin, barely different from the wobble frequency 0.036 × spin of the 300 axis: with moments this close, any damping, flexibility or torque decides what happens, and both 300 and 301 should be treated as marginal. The 500 axis is the clear major axis.

(d) H does not move in the inertial frame (torque-free). ω does: it moves in the body along the polhode and in space along the herpolhode on the invariable plane.

(e) The elliptic solution gives only the body rates; the attitude needs a further integration (Jacobi's closed form for it involves an elliptic integral of the third kind). The numerical quaternion is simpler and checkable.

(f) Rigidity: a spent stage carries residual propellant that sloshes and dissipates energy. Track its tumble over days (radar or light curves): a rigid body's period and polhode stay fixed; a dissipating one drifts towards spin about its major axis.

## Check-yourself answers

1. Because the inertia tensor is real and symmetric (and positive definite): such a matrix has real (positive) eigenvalues and orthogonal eigenvectors.
2. The kinetic energy and |H| (in any frame); the vector H is constant in the inertial frame.
3. H² > 2TI₂: around the major axis; H² < 2TI₂: around the minor axis.
4. The plane perpendicular to the fixed H containing the tip of the inertial ω; it lies at distance 2T/|H| from the origin.
5. The axis the polhode circles (the major axis if H² > 2TI₂, the minor if H² < 2TI₂); dn is never zero for m < 1, so that rate keeps its sign.
6. 4K(m)/λ; on the separatrix m = 1, K is infinite and the motion is not periodic.
7. By 2 ln 10/σ (σ the intermediate-axis growth rate): a fixed increment, not a factor.
8. 2π divided by the linearised wobble frequency of `P04.04`.
