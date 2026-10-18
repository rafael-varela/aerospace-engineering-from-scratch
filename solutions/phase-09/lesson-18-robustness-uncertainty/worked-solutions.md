# P09.18  Worked solutions

Worked for **s = 15**. Every number below is produced by `solution.py` (the functions named in each problem); the "independent check" line says how it was verified by a different route. Units: frequencies in rad/s, gains as ratios unless marked dB.

## 1. Peaks and margins

**(a)** Ms = 1.5 + 15/30 = **2.0**: gain range **(0.667, 2.0)**, phase margin ≥ **28.96°**.

**(b)** G₁ = 3/(s + 2): ‖G₁‖∞ = **1.5 at ω = 0**. G₂ with ω_n = 5, ζ = 0.1 + 0.15 = 0.25: peak 1/(2·0.25·√(1 − 0.0625)) = **2.06559 at ω = 5√(1 − 0.125) = 4.677 rad/s**.

**(c)** `tf_hinf_norm`: G₁ 1.5000000 (bisection) and 1.5 (sweep); G₂ 2.0655911 by both; the sweep locates the peak at 4.6771 rad/s.

**(d)** Not quite. Ms = 2 guarantees *at least* +6.02 dB upward and −3.52 dB downward gain margin, together with at least 29.0° of phase margin, all simultaneously. The actual gain margin may be larger (for the lesson's servo designs it is infinite on the nominal model).

**Common mistakes.** Taking the peak of the resonance at ω_n (the peak of |G| is below ω_n); quoting only the upward gain margin.

## 2. Derivations

**(a)** From y = Gu + d and u = C(r − y − n): y(1 + L) = L(r − n) + d, so **y = T(r − n) + Sd**, **r − y = S(r − d) + Tn** (since 1 − T = S), **u = KS(r − n − d)**. S + T = (1 + L)/(1 + L) = **1**.

**(b)** |1 + L(jω)| = 1/|S(jω)|, so **min_ω|1 + L| = 1/Ms**. The curve avoids the open disk of radius 1/Ms about −1, which covers the real interval (−1 − 1/Ms, −1 + 1/Ms). A crossing at −x satisfies x ≤ 1 − 1/Ms or x ≥ 1 + 1/Ms; scaling by k moves it to −kx, and the encirclement count changes only when kx = 1. For crossings with x ≥ 1 + 1/Ms that needs k ≤ Ms/(Ms + 1); for x ≤ 1 − 1/Ms it needs k ≥ Ms/(Ms − 1). So **k ∈ (Ms/(Ms + 1), Ms/(Ms − 1)) is safe**. At |L| = 1 the triangle (0, −1, L) is isosceles with sides 1, 1 and |1 + L| ≥ 1/Ms; the angle at 0 between −1 and L is 2 arcsin(|1 + L|/2) ≥ **2 arcsin(1/(2Ms)) = PM bound**.

**(c)** 1 + CG = 1 + L₀ + L₀WΔ = **(1 + L₀)(1 + TWΔ)**. If ‖WT‖∞ < 1 and ‖Δ‖∞ ≤ 1, then |TWΔ| < 1 on the whole Nyquist contour, the curve of 1 + TWΔ stays in the disk |z − 1| < 1 and cannot encircle 0, so the closed-loop poles stay those of the nominal loop's count: stable for every such Δ. Necessity: if |WT(jω*)| ≥ 1, take Δ(s) = c(a − s)/(a + s) with |c| = 1/|WT(jω*)| ≤ 1 and a chosen so that Δ(jω*) = −1/(WT)(jω*); then 1 + TWΔ = 0 at jω*, a closed-loop pole on the axis. Example: for K_c = 2.5 this Δ has size 0.627 and gives a root at 93.1671j (lesson part 4).

**(d)** Let H(γ)v = jωv with v = (p, q) and D = 0. First block row: Ap + BBᵀq/γ² = jωp, so p = (jωI − A)⁻¹B u/γ² with u = Bᵀq, and **Cp = G(jω)u/γ²**. Second block row: −CᵀCp − Aᵀq = jωq, so q = −(jωI + Aᵀ)⁻¹CᵀCp and u = Bᵀq = Bᵀ(−jωI − Aᵀ)⁻¹Cᵀ Cp = G(jω)ᴴ Cp (A, B, C real). Combining: **G(jω)ᴴG(jω)u = γ²u**, with u ≠ 0 (else p = 0 and q = 0). So γ is a singular value of G(jω). Conversely, from a right singular vector u with GᴴGu = γ²u, the same formulas build an eigenvector (p, q) with eigenvalue jω. Hence, for A Hurwitz, ‖G‖∞ ≥ γ ⇔ H(γ) has an imaginary-axis eigenvalue.

**Independent check.** `test_hamiltonian_test_switches_at_the_norm` confirms (d) numerically at 0.99 and 1.01 of the norm.

## 3. Waterbed

**Assumptions.** p = 1 + 15/10 = **2.5 s⁻¹**, ε = 0.3 for ω ≤ 2 rad/s, |S| ≈ 1 above ω₂.

**(a)** πp = 7.854. ln Ms ≥ (7.854 − 2 ln 0.3)/(30 − 2) = (7.854 + 2.408)/28 = **0.3665 ⇒ Ms ≥ 1.443**.

**(b)** ω₂ = 15 rad/s: ln Ms ≥ 10.262/13 = **0.7894 ⇒ Ms ≥ 2.202**. Halving the actuator bandwidth raises the unavoidable sensitivity peak from 1.44 to 2.20, i.e. the guaranteed phase margin from 2 arcsin(1/2.886) = 40.5° to 26.3°.

**(c)** Lesson airframe (p = 2): with leads (0.6, 1, 40), (1.0, 1.5, 60), (2.0, 2.0, 100), all stabilising, the integral is **6.283187, 6.283188, 6.283188** against 2π = 6.283185, while Ms is 1.995, 1.947 and 2.018: the controllers redistribute ln|S| but the area is fixed.

**(d)** The bound contains no controller parameter: it follows from analyticity of S and the unstable pole, for every stabilising controller. The airframe designer can reduce p (more static margin, a larger tail) or the actuator designer can raise ω₂.

## 4. Your servo (s = 15: K_c = 0.7)

**(a)** Ms = **1.16257**, Mt = **1.0000** (at ω = 0, type-1 loop), ‖KS‖∞ = **2.800 V/rad** (= K_cp/z, at high frequency); disk bounds: gain in **(0.5376, 7.151)**, PM ≥ **50.9°**.

**(b)** Sweep **1.16256651**, bisection **1.16256651**.

**(c)** ‖WT‖∞ = **0.577 < 1**: robust stability guaranteed for every actuator inside the weight.

**(d)** True loop: PM **68.3° at 23.08 rad/s**, GM **3.48**; its own Ms_true = 1.507 gives GM ≥ 2.971 and PM ≥ 38.7°: both satisfied. (The nominal loop has PM 80.6° and an infinite gain margin; the compliance mode and delay cost 12° and make the gain margin finite.)

## 5. Judgment

**(a)** The test is not wrong; it answers a different question. ‖WT‖∞ = 1.108 means that *some* plant inside the uncertainty set is destabilised (the construction of 2(c) produces it). The actual actuator is one particular member of the set, with one particular phase of Δ, and happens to be stable with 39° of margin. The test guarantees nothing about K_c = 1.6; it does not claim instability of the real servo.

**(b)** With the mode at 120 rad/s, max|Δ|/|W| = **0.997**: still covered, with 0.3 % to spare. Any lower and the weight is invalid; production data must bound the spread before the guarantee means anything.

**(c)** The test treats the gain error as any complex Δ of size up to 1/Mt. For the K_c = 1.2 design, Mt = 1.000, so it guarantees only complex relative errors below 100 %. A real gain factor k keeps the nominal loop stable for every k tested from 1.1 to 100 (and every k > 0, since the nominal phase never reaches −180°). The complex test is therefore very pessimistic about real gain errors: it cannot use the fact that a real gain does not rotate the Nyquist curve.

**(d)** Not on this evidence alone. ‖WT‖∞ = 0.974 is 2.6 % from the edge, and the true phase margin is 47° on a weight fitted to one unit. First measure the actuator's frequency response (swept sine to 300 rad/s) on several units and over temperature, and check |G/G₀ − 1| ≤ |W|.

## 6. Robust versus true limit

**Result.** (a) Largest K_c with ‖WT‖∞ < 1: **K_c = 1.2342**. (b) Largest K_c stable on the true servo: **K_c = 2.4365** (the same to five figures with Padé orders 3 to 6). (c) At K_c = 1.05 × 1.2342 = 1.2959: ‖WT‖∞ = 1.0476 at **61.62 rad/s**; the destabilising all-pass has ‖Δ‖∞ = **0.9546** and puts a closed-loop pole at **61.6210j** (real part 2 × 10⁻¹⁴).

**(d)** At 61.62 rad/s the destabilising perturbation is Δ = −1/(WT) = −0.902 − 0.311j; the true actuator's error in the same normalised units, Δ_true/W, is −0.427 − 0.146j. The phases agree (both near −161°) but the true error is only **0.45** of the size needed. So the factor-2 gap between the robust and the true gain limits comes from the *magnitude* of the weight: a first-order W that must reach 2.46 at the 150 rad/s resonance is about twice too large at 62 rad/s. A weight shaped like the resonance (Extension 1) would recover most of the gap.

**Common mistakes.** Bisecting on nominal stability (which never fails for this family); comparing the robust limit with the nominal limit.

## Check-yourself answers (lesson.md)

1. S = 1/(1 + L), T = L/(1 + L); S + T = (1 + L)/(1 + L) = 1. S should be small at low frequency: it multiplies disturbances and the reference tracking error, which are mostly slow.
2. The peak gain sup_ω|G(jω)| (largest singular value for a matrix); for 3/(s + 2) it is 1.5, at ω = 0.
3. Gain in (1.5/2.5, 1.5/0.5) = (0.6, 3) (−4.4 dB to +9.5 dB) and PM ≥ 2 arcsin(1/3) = 38.9°.
4. ‖WT‖∞ < 1 for G = G₀(1 + WΔ); robust stability means the closed loop is stable for every stable Δ with ‖Δ‖∞ ≤ 1, not only for the nominal plant.
5. Sufficient: 1 + CG = (1 + L₀)(1 + TWΔ) and |TWΔ| < 1 prevents encirclements. Necessary: if |WT| ≥ 1 at some ω*, an all-pass Δ of size ≤ 1 makes 1 + TWΔ vanish at jω*.
6. ∫₀^∞ ln|S(jω)|dω = 0 for a stable open loop (relative degree ≥ 2) and π Re p with one right-half-plane pole p.
7. The nominal model lacks the compliance mode and delay; near and above crossover they add phase lag and gain that the nominal margins do not see (K_c = 3.2 has 66° nominal and −58° true).
8. A real gain error does not rotate the Nyquist curve, but the test allows Δ of any phase, including the one that rotates L onto −1; it guards against a far larger set than the real one.
