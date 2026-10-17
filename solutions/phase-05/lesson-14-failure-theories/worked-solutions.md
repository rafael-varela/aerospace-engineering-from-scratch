# P05.14  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe scripts that call it, given to the figures shown. Material properties are the lesson's typical values.

## 1. A wing-skin panel (s = 15: σ_x = 195 MPa)

**Assumptions.** Plane stress; linear elastic up to the point of interest; 7075-T6 with F_ty = 503 MPa, F_tu = 572 MPa (typical).

(a) c = (195 − 40)/2 = 77.5 MPa; R = √(117.5² + 60²) = 131.93 MPa. **σ_a = 209.43 MPa, σ_b = −54.43 MPa.** With the out-of-plane zero: **σ₁ = 209.43, σ₂ = 0, σ₃ = −54.43 MPa**; `principal_stresses` agrees.

(b) **τ_max = (σ₁ − σ₃)/2 = 131.93 MPa; σ_T = 263.87 MPa; σ_vM = √(195² + 195 × 40 + 40² + 3 × 60²) = 241.30 MPa.**

(c) von Mises: **MS_yield = 503/241.30 − 1 = +1.085; MS_ult = 572/(1.5 × 241.30) − 1 = +0.580.** Tresca: **+0.906 and +0.445.** Ultimate governs (F_tu/F_ty = 1.14 < 1.5).

**Check.** σ_T/σ_vM = 1.094, inside [1, 1.155]. Here σ₁ and σ₃ have opposite signs, so the in-plane circle is the largest one and τ_max equals the in-plane R.

**Common mistakes.** Forgetting σ_y in the von Mises formula's cross term; taking σ_T = σ_a − σ_b without checking that the out-of-plane zero lies between them (it does here, but not for biaxial tension); applying 1.5 to the allowable instead of the stress.

## 2. Derivations

(a) With ε = ((1 + ν)σ − ν tr(σ)I)/E, u = ½ σ : ε = [(1 + ν) σ : σ − ν (tr σ)²]/(2E). Put σ = pI + s, p = tr(σ)/3, tr s = 0: σ : σ = 3p² + 2p tr s + s : s = 3p² + s : s. Then u = [(1 + ν)(3p² + s : s) − 9νp²]/(2E) = **(1 − 2ν)(tr σ)²/(6E)** + **(1 + ν) s : s/(2E)**, the first term using 3p² = (tr σ)²/3. In principal stresses, s : s = Σ(σ_i − p)² = ⅓[(σ₁ − σ₂)² + (σ₂ − σ₃)² + (σ₃ − σ₁)²] (expand both sides), so **u_d = ((1 + ν)/(6E))[(σ₁ − σ₂)² + (σ₂ − σ₃)² + (σ₃ − σ₁)²]**. In a tensile test at yield, the bracket is 2F_ty², so u_d,yield = (1 + ν)F_ty²/(3E). Setting u_d = u_d,yield gives ½[…] = F_ty², i.e. **σ_vM = F_ty**.

(b) Order σ₁ ≥ σ₂ ≥ σ₃ and write σ₂ = σ₃ + f D with D = σ₁ − σ₃ ≥ 0 and 0 ≤ f ≤ 1. Then σ₁ − σ₂ = (1 − f)D, σ₂ − σ₃ = fD, and σ_vM² = ½D²[(1 − f)² + f² + 1] = D²(1 − f + f²). On [0, 1], 1 − f + f² is 1 at the ends and 3/4 at f = ½. So D√(3/4) ≤ σ_vM ≤ D, i.e. **1 ≤ σ_T/σ_vM ≤ 2/√3**. The lower bound is reached when σ₂ equals σ₁ or σ₃ (uniaxial stress, equal biaxial stress, and any state plus a hydrostatic part); the upper bound when σ₂ is midway (pure shear, σ₂ = 0 with σ₁ = −σ₃).

(c) Pure shear: σ₁ = τ, σ₂ = 0, σ₃ = −τ. (7) with n = 1: τ/S_t + τ/S_c = 1, so **τ_f = S_tS_c/(S_t + S_c)**. With S_t = S_c = S, (7) reads (σ₁ − σ₃)/S = 1/n, i.e. σ₁ − σ₃ = S/n: **Tresca**.

**Common mistakes.** Losing the cross term 2p tr s (it is zero only because s is deviatoric); claiming the upper bound occurs in uniaxial tension.

## 3. A titanium lug

**Assumptions.** The FE stress is the true elastic stress at the element; Ti-6Al-4V F_ty = 880 MPa, F_tu = 950 MPa, E = 113.8 GPa, ν = 0.342 (typical).

(a) **σ₁ = 195.84, σ₂ = 95.07, σ₃ = −60.91 MPa.** **σ_vM = 224.05 MPa, σ_T = 256.74 MPa.** Margins: von Mises **MS_y = +2.928, MS_u = +1.827**; Tresca **+2.428, +1.467**. Ultimate governs (F_tu/F_ty = 1.08).

(b) With −100 MPa added to each normal stress: principal stresses **95.84, −4.93, −160.91 MPa** (each exactly 100 MPa lower); **σ_vM = 224.05 MPa and σ_T = 256.74 MPa, unchanged**, so the margins do not change. Both criteria depend only on differences of principal stresses.

(c) **u = 221.8 kJ/m³, u_v = 24.5 kJ/m³, u_d = 197.3 kJ/m³**: 89 % of the stored energy is distortion, and only that part enters the von Mises criterion. Check: (1 + ν)σ_vM²/(3E) = 1.342 × (224.05 × 10⁶)²/(3 × 113.8 × 10⁹) = 197.3 kJ/m³.

**Common mistakes.** Using only the 2 × 2 block of the tensor (τ_yz = 30 MPa matters); thinking that pressure "adds" to the equivalent stress.

## 4. Sizing the gear tube's wall (s = 15: loads × 1.15)

**Assumptions.** As the lesson; d_o = 50 mm fixed; the compression side stays critical (it does, P < 0).

(a) At t = 3 mm the scaled loads give σ_vM = 370.3 MPa and MS_y = +0.175, so a thinner wall is possible. `brentq` on MS_y(t) − 0.15: **t = 2.924 mm (von Mises), t = 3.035 mm (Tresca).**

(b) At t = 2.924 mm, **MS_ult = +0.181** by von Mises: higher than MS_yield = +0.15, so **yield governs**, as F_tu/F_ty = 1.54 > 1.5 says it must.

(c) A = π(d_o² − d_i²)/4: **3.39 kg/m (von Mises design) and 3.52 kg/m (Tresca design)**, a 3.5 % mass penalty for Tresca. In practice the wall would be rounded up to a stock gauge, and that rounding is larger than the difference between the criteria.

**Check.** At t = 3 mm the section is the lesson's and all three loads scale by the same factor, so every stress scales by it: 1.15 × 322.0 = 370.3 MPa, as computed.

**Common mistakes.** Bracketing the root badly (MS increases with t, so g(1 mm) < 0 < g(10 mm)); using the inner-surface radius for the critical point.

## 5. A brittle casting (s = 15: σ_x = ∓130 MPa)

**Assumptions.** Plane stress; typical gray cast iron S_t = 214, S_c = 752 MPa.

(a) σ_x = −130, τ = 80 MPa: **σ₁ = 38.08, σ₂ = 0, σ₃ = −168.08 MPa.** **n_CM = 1/(38.08/214 + 168.08/752) = 2.49; n_R = min(214/38.08, 752/168.08) = 4.47.**

(b) σ_x = +130, τ = 80 MPa: **σ₁ = 168.08, σ₃ = −38.08 MPa. n_CM = 1.20; n_R = 1.27.**

(c) The tensile principal stress is what opens cracks, and it is 38 MPa in (a) but 168 MPa in (b); cast iron is 3.5 times weaker in tension, so the same magnitude of loading is twice as dangerous in (b). Tresca with F = S_t would give n = 214/σ_T = 214/206.2 = 1.04 in *both* cases: it cannot tell them apart, because it treats tension and compression alike, and in (a) it would be needlessly pessimistic by a factor 2.4, in (b) pessimistic by 15 %. Rankine is unconservative in (a) because it ignores the help σ₃ gives σ₁ in opening a crack.

**Common mistakes.** Using the in-plane principal stresses without the zero, which matters when both are tensile; reporting a "margin" with the 1.5 factor for a brittle material without saying what it is applied to.

## 6. Judgment: a zero margin

A model answer, in outline. **What the margin is computed against.** The lesson's 2024-T3 F_tu = 483 MPa is a typical value. The A- or B-basis allowable for the actual product form (drawn tube, not sheet) is lower, so the true margin with a certified allowable is negative by both criteria. **Whether von Mises "wins".** It is the better yield criterion for ductile metals, but the disputed check is rupture at ultimate load, and there the 1.5σ extrapolation of an elastic stress is itself an approximation (conservative after yielding, since plastic redistribution caps the peak stress). The fact that a change of criterion flips the sign of the margin means the design sits exactly on the boundary of what the method can resolve. **What the analysis left out.** The fitting or weld at the critical section (stress concentration 2–3), fatigue under thousands of landings (2024-T3 has good fatigue strength but not at this stress level), corrosion in a gear bay, manufacturing tolerance on the wall. **Decision.** A zero margin leaves no room for any of these. Recommend keeping the 4130 tube or resizing an aluminium tube for a margin of at least +0.15 against MMPDS allowables, and ask first for the allowables of the actual product form and a detail stress analysis of the joint.

## Check-yourself answers

1. Both criteria depend only on differences between principal stresses (equivalently, on the deviatoric stress), and a hydrostatic state has none: every Mohr circle has zero radius and the distortion energy is zero.
2. Pure shear: σ₁ = 200, σ₃ = −200 MPa. **Tresca 400 MPa; von Mises √3 × 200 = 346.4 MPa.**
3. **150 MPa.** The in-plane principal stresses are both 150 MPa, but the out-of-plane principal stress is 0, so σ₁ − σ₃ = 150 MPa; the largest Mohr circle is the one through 0 and 150.
4. MS_ult = 480/(1.5 × 300) − 1 = **+0.067.**
5. F_tu/F_ty = 1.40 < 1.5, so **ultimate governs**: MS_ult < MS_yield for any stress.
6. **2/√3 = 1.155**, in pure shear (σ₂ midway between σ₁ and σ₃).
7. τ_f = 200 × 600/(200 + 600) = **150 MPa.**
8. Fatigue cracking under repeated landings; also buckling of a thin wall in compression, or fracture from a crack-like defect at a weld.
