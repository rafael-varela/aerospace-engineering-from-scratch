# P10.12  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and a probe script that calls it, given to the figures shown. Lifting line, flat sections, a0 = 2π throughout.

## 1. Taper and the stall (s = 15: AR = 9)

(a) e = **0.9286** (λ = 1), **0.9804** (0.5), **0.9855** (0.35).

(b) Stall-onset station η = **0** (root), **0.53**, **0.66**.

(c) Taper brings the *load* closer to elliptic, so e rises; but the load is divided by a chord that shrinks outboard, so the *local lift coefficient* rises towards the tips and the first-stalling station moves outboard.

**Check.** Same trend as the lesson's AR 8.33 (0.934, 0.982, 0.987; η = 0, 0.53, 0.66).

## 2. Derivations

(a) A = Bx is linear in x, so πAR Σ nAₙ² = πAR Σ n (Σ_k B_nk x_k)² = xᵀ(πAR Bᵀ diag(n) B)x = **xᵀQx**, with Q symmetric positive semi-definite. CL = πAR A₁ = πAR (Bx)₁ is linear in x.

(b) Minimise F = Σ n(αaₙ + trₙ)² subject to αa₁ + tr₁ = A₁. Lagrangian ∂/∂α: 2Σ naₙAₙ = λa₁; ∂/∂t: 2Σ nrₙAₙ = λr₁. Eliminate λ: Σ n(r₁aₙ − a₁rₙ)Aₙ = 0, i.e. Σ n qₙAₙ = 0 with qₙ = rₙ − r₁aₙ/a₁ (times −a₁). With Aₙ = pₙ + tqₙ (from the constraint): Σ nqₙpₙ + tΣ nqₙ² = 0, **t* = −Σ npₙqₙ/Σ nqₙ²**, which is (3). (q₁ = 0, so the n = 1 term drops out of both sums, as it must.)

(c) Each half-wing carries L/2, and an elliptic half-load's centroid is at 4(b/2)/(3π) = 2b/(3π) from the root. M_root = (L/2)(2b/(3π)) = Lb/(3π), and **C_B = M_root/(qSb) = CL/(3π)**.

(d) L = ρV∫Γ dy = ρV(b/2)∫₀^π Σ Gₙ sin nθ sin θ dθ = ρV(b/2)(π/2)G₁ ∝ bG₁. I = ρV∫Γy² dy = ρV(b/2)³∫ Σ Gₙ sin nθ cos²θ sin θ dθ, and sin θ cos²θ = (sin θ + sin 3θ)/4, so I ∝ b³(G₁ + G₃). D_i = (πρ/8)Σ nGₙ², independent of b. With the elliptic wing of span b₀ as unit (G₁ = 1): sG₁ = 1, s³(G₁ + G₃) = 1 ⇒ G₁ = 1/s, G₃ = 1/s³ − 1/s; D_i ratio **f(s) = 1/s² + 3(1/s³ − 1/s)²**. Positivity: Γ/G₁ = sin θ + r sin 3θ = sin θ(1 + 3r − 4r sin²θ); for r < 0 the minimum over θ is at sin θ → 0: 1 + 3r ≥ 0 ⇔ **r ≥ −1/3** (for r > 0, r ≤ 1).

**Common mistakes.** In (c), using the whole lift at b/(3π); in (d), forgetting that D_i does not depend on b.

## 3. A sailplane wing

**Assumptions.** b = 15 m, S = 11.25 m², λ = 0.4, flat sections, lifting line.

(a) Untwisted: e = **0.9649**, stall onset at **η = 0.62** (CL = 1.2).

(b) Unconstrained (5 stations): **e = 0.9974**, twist 0, +0.93, +1.55, +1.72, −2.71° (wash-in mid-span, washout at the tip), success. With a 5 % tip margin at CL = 1.2, washout only: **e = 0.9519**, twist 0, 0, 0, −2.50, −2.85°, success; the tip-margin constraint is active.

(c) Best linear twist: **−0.41°, e = 0.9653**.

| design | e at CL 0.8 | stall |
|---|---|---|
| untwisted | 0.9649 | η = 0.62 |
| linear −0.41° | 0.9653 | η = 0.55 |
| unconstrained 5-station | 0.9974 | not controlled |
| 5 % tip margin | 0.9519 | η = 0.41 |

Recommendation: the tip-margin design. It gives up 1.3 points of e at the inter-thermal CL against the untwisted wing (and 4.6 points against the unconstrained optimum), but it is the only one that guarantees an inboard stall at the thermalling CL, where the pilot flies slowly, banked and close to other gliders. Since half the flight is thermalling, a better design process would also optimise e at CL = 1.2.

**Common mistakes.** Optimising at the thermalling CL and then checking stall at cruise; reporting the unconstrained optimum without its stall behaviour.

## 4. Your own bending limit (s = 15: 0.925)

(a) e = **0.9569**; achieved C_B/C_B,ell = **0.9250**; twist 0, −0.63, −1.71, −3.03, −4.94°.

(b) Induced drag at 55 m/s: qS CL²/(πAR) = 93.0 N for e = 1, so 93.0/0.9569 − 93.0/0.99978 = **4.2 N more** than the elliptic optimum (97.2 N against 93.0 N; still below the untwisted wing's 99.6 N). The root moment is 7.5 % below elliptic and 14.6 % below the untwisted wing's (C_B 0.0325 against 0.0380): spar caps sized by the root moment could be about 7.5 % lighter than for the elliptic load.

(c) Active: the achieved ratio equals the limit to 10⁻⁷, and e is lower than the unconstrained 0.99978; an inactive constraint would leave the unconstrained optimum unchanged.

**Common mistakes.** Taking e from a different CL; forgetting that cap weight and drag are in different currencies (fuel and weight, `P03.16`).

## 5. Prandtl without the floor

(a)

| s | f(s) | G₃/G₁ |
|---|---|---|
| 1.1 | 0.9011 | −0.174 |
| 1.2247 | 0.8889 | −0.333 |
| 1.5 | 0.8560 | −0.556 |
| 2.0 | 0.6719 | −0.750 |

(b) Beyond √(3/2), G₃/G₁ < −1/3 and the load near the tips is negative (download). The proxy I = ∫L′y² dy counts download as *negative* structural weight, so the optimiser "pays" for span with tip download that reduces I; a real spar is sized by |M| and gets no lighter. Spans beyond the bound are artefacts of the proxy, not better wings.

(c) With s = 15, n_terms = 3 + 3 = 6: SLSQP returns G/G₁ = 1, −0.3333, 0, 0, 0, 0, span 1.2247, drag 0.8889. Higher harmonics enter only the drag (n Gₙ²), not L or I, so they can only add drag. A structural model whose weight depended on more than ∫Ly² dy (for example a root-bending-moment constraint, which weights the load by |y| and brings in every odd harmonic, or a shear constraint) would let them matter: that is Jones's 1950 problem (extension 3).

## 6. The trainer's planform

(a) "Taper 0.36 maximises the span efficiency of an *untwisted* wing (0.987 against 0.934, 5.4 N less induced drag at cruise), but it moves the first-stalling station to 66 % of the semi-span. Making it stall inboard with a 5 % margin at the approach CL needs about 4° of washout, which drops e at cruise to 0.766 (121 N, 22 N worse than the rectangular wing). For a trainer, the rectangular wing with about 3° of washout is better on both counts."

(b) Washout lowers the tips' angle at every CL, so at the approach CL their cl/CL falls further below the root's: the tip margin grows. Good for stall behaviour; the cost is a slightly non-optimal load at CLs far from cruise.

(c) (1) The tip sections working at a lower Reynolds number (smaller chord on a tapered wing; on a rectangular wing, the same chord but the tip vortex's local flow) or with a lower cl_max (aileron gaps, a different section). (2) Sideslip or a wing drop in a turn (one tip at a higher angle), propeller slipstream on the inboard wing delaying root stall, and manufacturing twist errors. Check with tufts in a stall flight test, and with the lattice plus a spanwise cl_max distribution.

## Check-yourself answers

1. Near λ ≈ 0.35 (0.36 at AR 8.33). It moves the stall outboard, towards the ailerons, and fixing that with washout costs more induced drag than taper saved.
2. The local lift coefficient is cl = 2Γ/(Vc): a nearly elliptic Γ divided by a chord that shrinks faster than Γ near the tips gives a higher cl there.
3. The Glauert coefficients are linear in α and the twists, and CDi is a sum of squares of the coefficients.
4. The elliptic load (e → 1), at that one design CL.
5. C_B = CL/(3π).
6. Non-negative lift along the span (with fixed lift and fixed ∫L′y² dy, span free): span × √(3/2) = 1.22, induced drag × 8/9 = 0.89.
7. Outboard, towards the tips: a tip stall takes away aileron control and can start a roll into a spin.
8. Twist fixes a basic load that does not scale with CL; at a higher CL the α-dependent (additional) load dominates and the total is no longer the optimum shape.
