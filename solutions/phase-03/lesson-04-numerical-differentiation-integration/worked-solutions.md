# P03.04  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and given to the figures shown.

## 1. Deceleration from a difference (s = 15: h = 47.5 km)

**Assumptions.** Default `EntryParams`: B = 280.88. SI units throughout; a in m/s².

(a) V(47.5 km) = 5317.29 m/s; a_exact = ρV²/2β = **78.745989 m/s²**. Central differences (error = numeric − exact):

| dh (m) | error (m/s²) |
|---|---|
| 100 | −6.79 × 10⁻⁶ |
| 50 | −1.69 × 10⁻⁶ |
| 25 | −4.22 × 10⁻⁷ |

Observed order **2.00**.

(b) The five-point stencil (V(h + 2δ) − 2V(h + δ) + 2V(h − δ) − V(h − 2δ))/(2δ³) with δ = 20 m gives V‴ = **−1.46 × 10⁻¹¹ m/s per m³**. Predicted error V sin γ × dh²V‴/6 = 5317 × 0.05234 × 10⁴ × (−1.46 × 10⁻¹¹)/6 = **−6.76 × 10⁻⁶ m/s²**, against −6.79 × 10⁻⁶ observed (0.5 %).

(c) h* = (3 × 2.2 × 10⁻¹⁶ × 5317/1.46 × 10⁻¹¹)^{1/3} = **0.62 m**. Sweep: dh = 10 m error −6.8 × 10⁻⁸ (still dh², ratio 100 from dh = 100 m ✓); 1 m: −7.8 × 10⁻¹⁰; 0.1 m: −7.8 × 10⁻¹⁰ (no further gain: the floor, 10⁻¹¹ relative); 0.01 m: +2.8 × 10⁻⁸ (rounding wins, and the sign is now random). The turning point is between 1 m and 0.1 m, consistent with h* = 0.6 m.

**Check.** Dimensions: V‴ in m/s/m³ times dh² in m² gives m/s/m = s⁻¹, times V sin γ in m/s gives m/s² ✓. The floor εV/dh × V sin γ at dh = 1 m is 2.2 × 10⁻¹⁶ × 5317 × 278 = 3 × 10⁻¹⁰ ✓ (same order as observed).

**Common mistakes.** Forgetting the V sin γ factor when converting the error in dV/dh to an error in a; estimating V‴ with the same dh as the main difference (the third difference is itself rounding-limited below δ ≈ 10 m).

## 2. Error laws and Allen–Eggers closed forms

(a) f(x ± h) = f ± hf′ + h²f″/2 ± h³f‴/6 + h⁴f⁗/24 ± …. Forward: [f(x + h) − f]/h = f′ + **hf″/2** + O(h²). Central: [f(x + h) − f(x − h)]/2h = f′ + **h²f‴/6** + O(h⁴). Second: [f(x + h) − 2f + f(x − h)]/h² = f″ + **h²f⁗/12** + O(h⁴). ∎

(b) On [0, h] the chord error is ½f″(ξ(t))t(t − h). With f″ continuous, the mean value theorem for integrals gives ∫₀ʰ = ½f″(ξ)∫₀ʰ(t² − ht)dt = ½f″(ξ)(h³/3 − h³/2) = **−h³f″(ξ)/12**. Summing n panels and replacing Σf″(ξᵢ) by nf″(ξ) (intermediate value theorem), with nh = b − a: E = **−(b − a)h²f″(ξ)/12**. ∎

(c) Centre the pair at 0: S = (h/3)[f(−h) + 4f(0) + f(h)]. For f = 1: 2h ✓; x: 0 ✓; x²: (h/3)(2h²) = 2h³/3 = ∫ ✓; x³: 0 = ∫ ✓. Exact for cubics. For x⁴: ∫ = 2h⁵/5, S = 2h⁵/3, S − ∫ = 4h⁵/15 = (h⁵/90)·24 = h⁵f⁗/90. Written as −(width)h⁴f⁗/180 with width 2h: −(2h)h⁴f⁗/180 has magnitude h⁵f⁗/90 ✓ (the sign convention is exact minus rule). Over n/2 pairs: **E = −(b − a)h⁴f⁗(ξ)/180**. ∎

(d) Dividing the two rates: dV/dh = ρV/(2β sin γ). With ρ = ρ₀e^{−h/H}: d ln V = [ρ₀/(2β sin γ)]e^{−h/H}dh. Integrating from h = ∞ (V = V_e): ln(V/V_e) = −[ρ₀H/(2β sin γ)]e^{−h/H} = −Be^{−h/H}, i.e. **V = V_e exp(−Be^{−h/H})** (4). Deceleration a = −dV/dt = −(dV/dh)(dh/dt) = (dV/dh)V sin γ = **ρV²/(2β)**. With u = e^{−h/H}: ln a = const − h/H − 2Bu, d/dh = −1/H + 2Bu/H = 0 ⇒ u = 1/(2B), **h_a = H ln 2B** = 7200 ln 561.8 = **45 584 m**. There V² = V_e²e⁻¹ and ρ = ρ₀/(2B) = β sin γ/H, so a_max = V_e²e⁻¹ β sin γ/(2βH) = **V_e² sin γ/(2eH)** = 7800² × 0.052336/(2 × 2.71828 × 7200) = **81.35 m/s² = 8.29 g**. Heating: ln q̇ = const + ½ln ρ + 3 ln V = const − h/(2H) − 3Bu; derivative −1/(2H) + 3Bu/H = 0 ⇒ u = 1/(6B), **h_q = H ln 6B = 53 494 m**. ∎

(e) Q = ∫q̇ dt = ∫_{0}^{h_top} q̇ dh/(V sin γ) (dh/dt = −V sin γ, and the time integral runs downward). q̇/V = k√(ρ₀/Rₙ)u^{1/2}V² = k√(ρ₀/Rₙ)V_e²u^{1/2}e^{−2Bu}. Substituting u = e^{−h/H}, dh = −H du/u: Q = [k√(ρ₀/Rₙ)V_e²H/sin γ]∫_{u_top}^{1}u^{−1/2}e^{−2Bu}du. With v = √(2Bu): u^{−1/2}du = (2/√(2B))dv and e^{−2Bu} = e^{−v²}, so ∫ = (2/√(2B))·(√π/2)[erf(√(2B)) − erf(√(2Bu_top))] = **√(π/2B)[erf(√(2B)) − erf(√(2Bu_top))]** (5). Numerically C = 1.7415 × 10⁻⁴ × 1.1068 × 7800² × 7200/0.052336 = 1.6132 × 10¹⁰ J/m²; √(π/2B) = 0.074786; erf(√561.8) = 1, erf(√(561.8 × 5.83 × 10⁻⁸)) = erf(5.72 × 10⁻³) = 6.46 × 10⁻³; Q = 1.6132 × 10¹⁰ × 0.074786 × 0.99354 = **1.1987 × 10⁸ J/m² = 119.87 MJ/m²** ✓.

**Check.** a_max independent of β and of ρ₀: a heavier or a thinner atmosphere only moves the peak in altitude ✓ (classical Allen–Eggers result). Q ∝ V_e²: doubling the entry speed quadruples the heat load, which is why lunar return (11 km/s) is a different heat-shield problem ✓.

**Common mistakes.** Integrating (4) from h = 0 instead of ∞ (gives V_e at the ground); a sign slip in d ln ρ/dh; forgetting that dt = dh/(V sin γ) brings a 1/V into the heat-load integrand (q̇ ∝ V³ becomes V² in h).

## 3. Total impulse of a thrust curve (s = 15: τ = 2.0 s)

**Assumptions.** F in N, t in s, impulse in N·s; the record is exact at the samples (no noise).

(a) I = ½ × 2000 × 0.2 + 2000τ[1 − e^{−(6 − 0.2)/τ}] = 200 + 4000(1 − e^{−2.9}) = 200 + 4000 × 0.94498 = **3979.907 N·s**.

(b) Δt = 0.1 s (n = 60, kink on sample 2): trapezoid **+0.787 N·s** (2.0 × 10⁻⁴ relative), Simpson **+1.31 × 10⁻⁴ N·s** (3.3 × 10⁻⁸). Δt = 0.15 s (n = 40, kink at t = 0.2 inside [0.15, 0.30]): trapezoid **−25.7 N·s** (6.5 × 10⁻³), Simpson **+0.0413 N·s** (1.0 × 10⁻⁵). Explanation of the Simpson pair: on the 0.1 s grid, every Simpson pair of panels lies in one smooth piece (the kink is at the boundary of the first pair), so the h⁴ law holds on each piece and the error is the Euler–Maclaurin term of the exponential tail, h⁴[F‴(6) − F‴(0.2⁺)]/180 = 10⁻⁴ × (−13.8 + 250)/180 = 1.31 × 10⁻⁴ ✓ to three figures. On the 0.15 s grid the pair [0.15, 0.45] contains the corner; a parabola through (0.15, 1500), (0.30, 1902.5), (0.45, 1766.8) misses the triangle apex at 2000 N and the error is O(h²) at the kink (300× larger), while the trapezoid rule, which also straddles the kink with the chord from 1500 to 1902.5 N, loses the whole area above the chord: −25.7 N·s.

(c) Predictions: ×¼ for the trapezoid rule, ×1/16 for Simpson. Observed at Δt = 0.05 s (n = 120): trapezoid +0.197 (ratio **4.00**), Simpson +8.20 × 10⁻⁶ (ratio **16.0**). Both laws hold because the kink sits on a node of both grids. Trapezoid prediction at 0.1 s: the first piece is linear (exact); the tail contributes +h²[F′(6) − F′(0.2⁺)]/12 with F′ = −1000e^{−(t−0.2)/τ}: 0.01 × (−55.0 + 1000)/12 = **+0.7875 N·s** against +0.7874 observed. The kink itself contributes nothing because the trapezoid rule is exact on each linear piece and the two pieces meet at a node.

**Check.** Relative Simpson error at 0.1 s, 3 × 10⁻⁸, is 1 N·s per 3 × 10⁷: for a motor test this is far below the load-cell accuracy (10⁻³), so the 0.15 s grid's 10⁻⁵ is also acceptable; the point is to know which regime you are in, not to buy more digits than the sensor has.

**Common mistakes.** Sampling with `numpy.linspace(0, 6, n)` (n points, n − 1 panels, kink off-node for n = 61); using Simpson with odd n; blaming the integrator when the record's kink falls between samples: that is a sampling choice, and a test stand that records at 0.1 s puts the ignition ramp on nodes if the ramp duration is a multiple of the period.

## 4. A heavier capsule (s = 15: β = 350 kg/m²)

(a) B = ρ₀H/(2β sin γ) = 1.225 × 7200/(2 × 350 × 0.052336) = **240.75**. h_q = H ln 6B = **52.38 km** (down 1.11 km from 53.49; H ln(300/350) = −1.110 km ✓). Peak deceleration **81.35 m/s² = 8.29 g** at h_a = H ln 2B = **44.47 km** (unchanged value, 1.11 km lower ✓). Q(120 km) = **129.537 MJ/m²**; ratio to the β = 300 value 119.871 is 1.0806 against √(350/300) = 1.0801 (the erf(√(2Bu_top)) correction accounts for the 0.05 %).

(b) Relative errors of the heat load from 120 km:

| n | trapezoid | Simpson | Gauss |
|---|---|---|---|
| 16 | 2.18 × 10⁻⁵ | 6.98 × 10⁻³ | 2.30 × 10⁻³ |
| 32 | 3.38 × 10⁻⁵ | 3.78 × 10⁻⁵ | 1.77 × 10⁻⁶ |
| 64 | 8.46 × 10⁻⁶ | 1.51 × 10⁻⁸ | 3.0 × 10⁻¹³ |
| 128 | 2.11 × 10⁻⁶ | 5.97 × 10⁻¹⁰ | 1.2 × 10⁻¹⁶ |
| 256 | 5.29 × 10⁻⁷ | 3.73 × 10⁻¹¹ | 3.5 × 10⁻¹⁶ |

Orders over 32–128: trapezoid **2.00**, Simpson **7.97** (pre-asymptotic: the n = 32 entry is still far above its h⁴ line, so the apparent slope is twice the true one; over 128–256 it is 4.0). Note the trapezoid error at n = 16 is *below* its n = 32 error: at 7.5 km panels the Euler–Maclaurin series has not settled and the error is non-monotone.

(c) Smallest even n for 10⁻⁶: trapezoid **188**, Simpson **44**, Gauss **30**. Recommendation: with the integrand as a function, Gauss (30 evaluations, and 64 buys 10⁻¹³); with samples, Simpson on 44 or more panels (≤ 2.7 km spacing), never the trapezoid rule at this tolerance.

**Check.** a_max identical to Problem 2's 81.35 m/s² although β changed ✓ (independence of β). Simpson error ratio 128 → 256 is 16.0 ✓.

**Common mistakes.** Reading the Simpson order 7.97 as evidence of a better method (it is evidence of a pre-asymptotic grid); forgetting that `heat_load_numeric` takes n as panels (Simpson needs even n).

## 5. When the laws fail (s = 15: a = 0.65)

(a) Simpson reaches 10⁻¹⁰ on ∫₀² e^{−x}cos 3x dx at **n = 380** (381 evaluations); Gauss at **n = 8** (8 evaluations). Factor 48 in cost for the same accuracy.

(b) Romberg diagonal errors to level 6: 1.1, 1.8 × 10⁻¹, 5.1 × 10⁻⁵, **1.8 × 10⁻⁴**, 1.1 × 10⁻⁶, 1.1 × 10⁻⁹. Level 3 (three extrapolations over the 1-, 2-, 4- and 8-panel trapezoid sums, error term ∝ h⁸f⁽⁸⁾) is worse than level 2 (Boole's rule on 4 panels, error ∝ h⁶f⁽⁶⁾) because each extrapolation assumes the previous column already follows its power law, and the 1- and 2-panel sums it leans on sample an integrand of period 2π/3 ≈ 2.1 once or twice over a length 2: they carry no information about the oscillation, and f⁽⁸⁾ brings a factor 3⁸ = 6561 into the constant. Richardson is only as good as the asymptotic assumption behind it; from level 4 (16 panels, 8 per period) the assumption holds and each level gains two to three digits.

(c) ∫₀¹ x^0.65 dx = 1/1.65 = 0.606061. Observed orders over n = 8 … 128: trapezoid **1.61**, Simpson **1.65**, Gauss **3.24**. The integrand has f″ ∝ x^{−1.35}, unbounded at 0, so the error laws (2)–(3), which need f ∈ C² and C⁴, do not apply; the first panel dominates and its error scales as ∫₀ʰ x^a dx ∝ h^{a+1}: order a + 1 = 1.65 for both Newton–Cotes rules (Simpson's extra accuracy buys nothing). Gauss has no node at 0 and its error falls like n^{−(2a+2)} = n^{−3.3}: algebraic, not spectral, because the singularity is on the interval rather than at a distance from it. `scipy.integrate.quad` (QUADPACK qags) bisects adaptively toward 0 and applies the ε-algorithm to the sequence of estimates, which recovers 10⁻¹⁴ in about 200 evaluations; the manual fix is the substitution x = t² (integrand 2t^{2a+1}, smooth) or a Gauss–Jacobi rule with the weight x^a built in (18.335J Lecture 36).

**Check.** Simpson n = 128 error 1.53 × 10⁻⁵ against n = 64 4.81 × 10⁻⁵: ratio 3.14 = 2^{1.65} ✓.

**Common mistakes.** Reporting the Gauss order 3.24 as "fourth order" (it is 2a + 2 and depends on a); assuming Simpson's degree-3 exactness guarantees h⁴ convergence (it needs f⁗ bounded).

## 6. A quadrature rule for a trajectory tool

(a) The trapezoid rule, Simpson and Romberg (Richardson over trapezoid sums) use equally spaced samples and run on the record directly. Gauss–Legendre needs the integrand at the irrational Legendre nodes, and adaptive Simpson needs to request new points; both require a *function*, which the samples are not unless you interpolate them (`P03.03`), which adds its own error.

(b) From the lesson's table (default capsule), Simpson's relative error is 8.4 × 10⁻³ at 16 panels (7.5 km), 2.7 × 10⁻⁵ at 32 (3.75 km), 8.8 × 10⁻⁹ at 64 (1.9 km); the 10⁻⁶ crossing is near n = 44, i.e. **2.7 km spacing**. The delivered 40–400 m spacing is 7–70 times finer: the requirement is met with margin (at 1 km spacing Simpson is at 8.3 × 10⁻¹⁰).

(c) Report |S_h − S_{2h}|/15 where S_h is Simpson on all samples and S_{2h} Simpson on every other sample. For the default capsule at 1 km spacing (n = 120 and 60): S₁₂₀ error 8.3 × 10⁻¹⁰, S₆₀ error 2.2 × 10⁻⁸, estimate |S₁₂₀ − S₆₀|/15 = **1.4 × 10⁻⁹ relative**, which bounds the true 8.3 × 10⁻¹⁰ by a factor 1.7: a conservative, cheap and honest number. (If the time samples are not uniform in altitude, integrate in time instead: Q = ∫q̇ dt, where the samples are uniform.)

(d) For the thermal group: "The numerical error of the heat load is below 10⁻⁸ at your sample spacing and we attach a Richardson estimate to every value, so 10⁻⁶ is met at no cost. The number itself is uncertain by about 10 % because of the heating correlation and the atmosphere model, and that uncertainty does not shrink with more samples. If the requirement is on the physical heat load, it should be stated as a margin on the model (e.g. +20 % on q̇), not as a digit count on the integral."

**Check.** 1.4 × 10⁻⁹ ≥ 8.3 × 10⁻¹⁰: the estimate is an upper bound here as the h⁴ law predicts (|S_h − S_{2h}| ≈ 15|E_h| when the error is in its asymptotic regime, and the estimate over-predicts when S_{2h} is still pre-asymptotic) ✓.

**Common mistakes.** Choosing Gauss for sampled data; quoting the integration error as the uncertainty of the heat load; reporting a heat load to ten figures from a model good to one.

## 7. Gauss–Legendre: exactness and positive weights

(a) The rule integrates the interpolant exactly: Σwᵢf(xᵢ) = Σf(xᵢ)∫ℓᵢ = ∫Σf(xᵢ)ℓᵢ = ∫p_{n−1}, where p_{n−1} is the interpolant of f at the nodes. If f is a polynomial of degree ≤ n − 1 its interpolant is itself (uniqueness, `P03.03` Problem 2a), so the rule is exact. ∎

(b) Divide: p = qPₙ + r, deg q, r ≤ n − 1. ∫p = ∫qPₙ + ∫r. The first integral is zero because Pₙ is orthogonal on [−1, 1] to every polynomial of lower degree (the defining property of the Legendre polynomials). The second equals Σwᵢr(xᵢ) by (a). At the nodes p(xᵢ) = q(xᵢ)Pₙ(xᵢ) + r(xᵢ) = r(xᵢ) because Pₙ(xᵢ) = 0 (the nodes are its roots). Hence ∫p = Σwᵢp(xᵢ). Orthogonality is used on the integral, the root property on the sum. ∎

(c) Let ω(x) = Π(x − xᵢ)², degree 2n. ∫ω > 0 since ω ≥ 0 and vanishes only at n points; but Σwᵢω(xᵢ) = 0. So no rule with these (or any) n nodes is exact for degree 2n: 2n − 1 is the maximum. ∎

(d) ℓᵢ² has degree 2n − 2 ≤ 2n − 1, so by (b) ∫ℓᵢ² = Σⱼwⱼℓᵢ(xⱼ)² = Σⱼwⱼδᵢⱼ = wᵢ. The left side is the integral of a non-negative, not identically zero, continuous function: wᵢ > 0. ∎

(e) `gauss_legendre_nodes(5)`: nodes ±0.9061798459, ±0.5384693101, 0; weights 0.2369268851, 0.4786286705, 0.5688888889 (×2 for the symmetric pairs), sum **2.000000000000000**. Σwᵢxᵢ⁹ = **2 × 10⁻¹⁹** (exact 0: degree 9 = 2n − 1 ✓). Σwᵢxᵢ¹⁰ = **0.178886** against the exact 2/11 = 0.181818: wrong in the second digit, as (c) predicts for degree 2n = 10.

**Check.** The degree-10 error 2.9 × 10⁻³ equals ∫Pₙ²·(leading-coefficient ratio)… more simply: by (c) the rule's error on x¹⁰ equals its error on the monic Π(x − xᵢ)² = P₅²/(63/8)², namely ∫P₅²/(63/8)² = (2/11)/(63/8)² = 2.93 × 10⁻³ ✓ (P₅ has leading coefficient 63/8 and ∫P₅² = 2/11).

**Common mistakes.** Claiming exactness to degree 2n; proving positivity from "weights are integrals of ℓᵢ", which are not positive in general (Newton–Cotes weights go negative from 9 points); using (b) without noting where each property of Pₙ enters.

## Check-yourself answers

1. Forward hf″/2 (order 1), central h²f‴/6 (order 2), second difference h²f⁗/12 (order 2).
2. Truncation h²F₃/6 against rounding εF₀/h balance at h ~ ε^{1/3} ≈ 10⁻⁵, giving an error ~ ε^{2/3} ≈ 10⁻¹¹: two thirds of the digits. The second difference has rounding 4εF₀/h², optimum h ~ ε^{1/4} ≈ 10⁻⁴ and accuracy ~ ε^{1/2} ≈ 10⁻⁸.
3. (2ᵖD(h/2) − D(h))/(2ᵖ − 1) cancels the Chᵖ term; applied to the trapezoid rule (p = 2) it gives Simpson's rule, and the first column of Romberg's table.
4. Trapezoid −(b − a)h²f″/12; Simpson −(b − a)h⁴f⁗/180. On a symmetric pair of panels the cubic term of f is odd about the midpoint, so it integrates to zero and Simpson's symmetric weights also give zero: exact, and the error starts at the quartic term.
5. When all the differences of odd derivatives at the two ends vanish (Euler–Maclaurin), i.e. for a periodic integrand over a whole period or one that decays with all derivatives at both ends; then the error decays exponentially in n. Example: ∫₀²π e^{cos x}dx to 10⁻¹⁴ with 16 panels.
6. Choosing nodes at the roots of Pₙ makes the rule exact for r (degree ≤ n − 1) by interpolation and for qPₙ by orthogonality plus Pₙ(xᵢ) = 0, so for all degrees ≤ 2n − 1; the nodes are irrational and fixed, so the integrand must be evaluated there, which sampled data does not permit.
7. Deceleration at H ln 2B = 45.6 km, heating at H ln 6B = 53.5 km: the heating peak comes first (higher), while V is still 6.6 km/s.
8. Sixteen 7.5 km panels do not resolve a 7 km pulse; Simpson's h⁴ law has not begun and its correction (4T_h − T_{2h})/3 is of the wrong size, while the trapezoid rule profits from the near-vanishing end derivatives. The test: halve h and check that the error falls by the predicted factor (4 or 16); until it does, you are not in the asymptotic regime.
