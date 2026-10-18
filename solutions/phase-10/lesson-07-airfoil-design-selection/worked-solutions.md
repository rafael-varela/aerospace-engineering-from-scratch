# P10.07  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and a probe script that calls it, given to the figures shown. All coefficients are dimensionless.

## 1. A two-term design by hand

**Assumptions.** Thin-airfoil theory; design at the ideal angle (A0 = 0).

(a) A1 = 0.5/π = **0.1592**; A2 = 0.1592 + 4(−0.08)/π = **0.0573**; α_i = A2/3 = 0.01910 rad = **1.094°**.

(b) dz/dx = α_i + A1 cos θ + A2 cos 2θ. With sin θ cos θ = ½ sin 2θ and sin θ cos 2θ = ½(sin 3θ − sin θ):
z(θ) = ½[α_i(1 − cos θ) + (A1/2) sin²θ + (A2/2)((1 − cos 3θ)/3 − (1 − cos θ))].
On a fine grid the maximum is **4.19 % at x = 0.396**. The slope is +0.236 at the leading edge and −0.083 at the trailing edge.

(c) Forward: cl = **0.5000**, cm_c/4 = **−0.0800**, A0 = −7 × 10⁻¹⁸, α_L0 = **−3.465°**.

**Check.** cl = 2π(α_i − α_L0) = 2π(0.01910 + 0.06048) = 0.500 ✓. The camber is about twice the 2412's for 1.95 times its ideal cl: the camber is linear in the load.

**Common mistakes.** Forgetting the factor ½ in dx = ½ sin θ dθ; using cm about the leading edge.

## 2. Derivations

(a) With A0 = 0, (1) gives γ/(2V∞) = Σ Aₙ sin nθ and Δcp = 2γ/V∞, so Δcp/4 = Σ Aₙ sin nθ. Multiplying by sin mθ and integrating over [0, π], orthogonality (∫ sin nθ sin mθ dθ = (π/2)δ_nm) gives **Aₙ = (2/π)∫₀^π (Δcp/4) sin nθ dθ**. The flow-tangency half of (1) is dz/dx = (α − A0) + Σ Aₙ cos nθ, and with A0 = 0 the angle is α_i: **dz/dx = α_i + Σ Aₙ cos nθ**. The coefficients are shared because the vortex sheet's induced downwash of the term sin nθ is the term cos nθ (Glauert's integral).

(b) z(1) − z(0) = ∫₀¹ dz/dx dx = ½∫₀^π (α_i + Σ Aₙ cos nθ) sin θ dθ = α_i + ½ Σ Aₙ Iₙ = 0, with Iₙ = ½∫₀^π [sin(n + 1)θ − sin(n − 1)θ] dθ = ½[(1 − cos(n + 1)π)/(n + 1) − (1 − cos(n − 1)π)/(n − 1)] = (1 + cos nπ)/(1 − n²). **Odd n: Iₙ = 0** (I₁ = ∫ cos θ sin θ dθ = 0 too). **Even n: Iₙ = −2/(n² − 1)**. So α_i = −½Σ Aₙ Iₙ = Σ_{even} Aₙ/(n² − 1); for two terms α_i = A2/3.

(c) Aₙ = (2/π)(cl_i/4)∫₀^π sin nθ dθ = (cl_i/2π)(1 − cos nπ)/n = **cl_i/(πn)** for odd n, 0 for even n; so α_i = 0. Then dz/dx = (cl_i/π)Σ_{odd} cos nθ/n = −(cl_i/2π) ln tan(θ/2). Since cos θ = 1 − 2x, tan²(θ/2) = (1 − cos θ)/(1 + cos θ) = x/(1 − x), so **dz/dx = (cl_i/4π) ln((1 − x)/x)**. Integrating, ∫ ln(1 − x) dx = −(1 − x) ln(1 − x) + (1 − x) and ∫ ln x dx = x ln x − x, so z = (cl_i/4π)[−(1 − x) ln(1 − x) − x ln x] + const, and z(0) = 0 fixes the constant at zero: **z = −(cl_i/4π)[x ln x + (1 − x) ln(1 − x)]**.

(d) The error in z is the tail E = Σ_{n>N} Aₙ Jₙ(θ). Split Jₙ = −1/(n² − 1) + ½[cos(n − 1)θ/(n − 1) − cos(n + 1)θ/(n + 1)]. The steady part gives Σ_{n>N} (cl_i/πn)(1/n²) ≈ cl_i/(2πN²). The oscillating part is a sum of terms (1/n) cos(n ± 1)θ/(n ± 1), of size 1/n² with a phase that turns with n; summation by parts bounds it by about 1/(N² sin θ), which is also O(1/N²) except right at the ends (where z → 0 anyway). So **E = O(1/N²)**, and the lesson's observed orders 1.88, 2.00, 2.00 confirm it.

**Common mistakes.** Writing Δcp = γ/V∞ (the factor 2 matters); integrating the slope with the wrong constant.

## 3. A section for a flying wing

**Assumptions.** Thin-airfoil design; the panel method for the thick section.

(a) cm = 0 gives A2 = A1 = 0.3/π = **0.0955**, α_i = **1.824°**. The camber rises to **3.08 % at x = 0.310**, then falls through zero and goes slightly negative, **−0.093 % at x = 0.940**, ending with an upward slope of +0.032 at the trailing edge: a **reflexed** mean line. The reflex unloads the rear of the section so that the load's centre sits at the quarter chord, which is what cm_c/4 = 0 means when cl > 0.

(b) 12 % thickness, panel method:

| section | α | cl | cm_c/4 |
|---|---|---|---|
| reflexed design | 1.824° | 0.344 | −0.0015 |
| reflexed design | 4.0° | 0.607 | −0.0049 |
| NACA 2412 | 0° | — | −0.0546 |
| NACA 2412 | 4° | — | −0.0601 |

The design keeps |cm| below 0.005, an order of magnitude below the 2412's.

(c) Thin-airfoil α_L0 = **−0.912°** against −2.077° for the 2412: at a given α the reflexed section gives 2π × 1.165° = **0.128 less cl**. A tailless aircraft pays for its self-trimming section with a lower cl at every angle, and so a larger wing or a higher angle.

**Check.** The thin-airfoil cm is exactly zero; the panel cm is −0.0015 at the design point because thickness moves the load slightly.

**Common mistakes.** Expecting zero cm to mean zero camber (that would also give cl_ideal = 0).

## 4. Your own roof-top (s = 15: a = 0.65)

(a) α_i = **1.171°**, cm_c/4 = **−0.0843**, maximum camber **3.64 % at x = 0.492**.

(b) Forward: cl = **0.5000**, largest load error on 2–98 % chord **0.0064**.

(c) 10 % thickness at α_i: cl = **0.551**, cm = **−0.0878**, mean load error **0.068** on 10–90 %. (b) measures the **code** (two different projections of the same theory agree to the truncation error); (c) measures the **model** (thickness, which the theory leaves out, moves cl by 10 % and the load by 0.07).

**Check.** Moving the break aft from 0.6 to 0.65 moves the load centre aft, so cm becomes more negative (−0.0792 → −0.0843) and the ideal angle smaller (1.292° → 1.171°) ✓.

**Common mistakes.** Not clipping the slope at x = 0 and 1 when thickening; reading (c) as a bug.

## 5. Is the selection robust? (s = 15)

(a) Weights [0.198, 0.230, 0.326, 0.246] (thickness, cl matching, moment, suction peak). Totals: 4415 **0.739**, 2412 **0.690**, a = 1.0 cl_i = 0.4 0.653, 4412 0.630, 0012 0.504, cl_i = 0.8 0.340, 6409 0.330. The 4415 still wins, by 0.049.

(b) Without the 0012, nominal weights:

| section | with 0012 | without |
|---|---|---|
| NACA 2412 | 0.623 | 0.497 |
| NACA 4412 | 0.654 | 0.639 |
| NACA 4415 | 0.818 | 0.817 |
| NACA 6409 | 0.403 | 0.417 |
| a = 1.0, cl_i = 0.4 | 0.629 | 0.539 |
| a = 1.0, cl_i = 0.8 | 0.494 | 0.487 |

The 0012 had the worst suction peak (5.44) and the best moment (0). Without it, the 2412 (peak 3.29) becomes the worst on the peak criterion and scores 0 there instead of 0.51, so its total falls by 0.13. Min–max normalisation makes every score relative to the set, which is why the candidate list belongs in the report.

(c) For these weights the 2412 wins once the moment weight exceeds **0.430** (and the 0012 above 0.604). Recommendation: "On the four inviscid criteria the NACA 4415 is preferred for the trainer under the nominal and the seed-15 weights. The NACA 2412 would be preferred if the pitching moment were weighted above 0.43–0.54, so the decision hinges on how much trim drag and wing torsion matter, and should be settled with drag and cl_max data rather than by the matrix."

**Common mistakes.** Normalising each candidate rather than each criterion; reporting only the winner.

## 6. What is the matrix missing?

(a) (1) **Profile drag at the cruise and climb cl**, from `P10.06` (viscous–inviscid coupling) and the measured polars of TR-824. (2) **cl_max and the type of stall** (trailing-edge, gentle; or leading-edge, abrupt), from TR-824's data and `P10.15`. (3) **Sensitivity to Reynolds number and to surface roughness** (TR-824's rough-leading-edge data, `P10.14`). A fourth: manufacturing tolerance of the trailing edge and of the camber.

(b) (1) Its viscous polar at the trainer's Reynolds number (about 4.5 × 10⁶ at 55 m/s and 1.2 m chord): run `P10.06` or XFOIL, or test it. (2) The behaviour of its log-singular camber slope near the trailing edge after thickening (a cusp-like rear), checked by the panel method at several angles and by a structural check of the trailing-edge thickness.

(c) Tuning weights to reach a predetermined answer turns the matrix into decoration: it then records the conclusion rather than testing it. The honest course is to add the criteria that make the 2412 good (drag, cl_max, stall type), keep the weights stated before the scores were seen, and report the flip points.

## Check-yourself answers

1. The load Δcp/4 is the sine series Σ Aₙ sin nθ (plus the A0 leading-edge term); the camber slope is the cosine series α − A0 + Σ Aₙ cos nθ. The same coefficients appear in both, because the sheet term sin nθ induces the downwash term cos nθ.
2. A0 multiplies the leading-edge singularity of the load (a suction spike). A0 = 0 gives smooth entry of the flow at the nose; the angle at which this happens is the ideal angle α_i (its cl is the ideal cl).
3. Zero: a symmetric load has only odd-n terms, whose closure integrals Iₙ vanish.
4. The NACA a = 1.0 mean line z = −(cl_i/4π)[x ln x + (1 − x) ln(1 − x)]; maximum (cl_i/4π) ln 2 = **0.0331 c** for cl_i = 0.6.
5. The coefficients fall as 1/n and the camber kernels Jₙ as 1/n, and the oscillation of the kernels adds a further cancellation; the tail of terms of size 1/n² sums to order 1/N² (observed 2.00).
6. Thickness effects on the surface pressures (and on cl) and viscous effects (drag, transition, separation, cl_max). XFOIL's inverse modes design thick sections with a boundary layer; the panel method plus `P10.06` checks them.
7. Each criterion is scaled by its spread over the candidate set; removing the candidate that sets a minimum or maximum changes the scale, so every other candidate's score on that criterion changes, and the weighted totals can reorder.
8. With the other weights rescaled in proportion, each total is linear in the varied weight, so the winner changes only where two straight lines cross, which can be computed exactly. They tell a reviewer how far the stated weights are from changing the decision, which the totals alone hide.
