# P07.14  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe script that calls it, given to the figures shown.

## 1. A sunshade tube (s = 15: h = 0.25 m, r = 0.05 m)

**Assumptions.** Diffuse surfaces; the mouth is a flat disk; each surface uniform.

(a) F₁₃ = `vf_coaxial_disks(0.05, 0.05, 0.25)` = **0.03709**.

(b) F₁₁ = 0 (flat), so F₁₂ = 1 − F₁₃ = **0.96291**. With A₁ = πr² = 7.854 × 10⁻³ m² and A₂ = 2πrh = 7.854 × 10⁻² m², reciprocity gives F₂₁ = A₁F₁₂/A₂ = **0.09629**. By symmetry F₂₃ = F₂₁ = **0.09629**. Summation: F₂₂ = 1 − 2 × 0.09629 = **0.80742**.

(c) Rows: (0, 0.96291, 0.03709), (0.09629, 0.80742, 0.09629), (0.03709, 0.96291, 0): each sums to 1. A₁F₁₂ = A₂F₂₁ = 7.563 × 10⁻³ m², and similarly for the other pairs.

(d) The side wall is concave: a point on it sees other parts of the same wall across the tube. A flat or convex surface cannot see itself. A deep, black-lined shade is effective because 96 % of what leaves the bottom (and of stray light entering) strikes the wall, where a black coating absorbs it.

**Common mistakes.** Setting F₂₂ = 0; forgetting the second end when using summation on the side.

## 2. Derivations

(a) A₁F₁₂ = ∫∫ cos θ₁ cos θ₂/(πs²) dA₁dA₂, which is symmetric in 1 ↔ 2, so **A₁F₁₂ = A₂F₂₁**. Then Σ_i Q_i = Σ_i A_iJ_i − Σ_i Σ_j A_iF_ijJ_j = Σ_i A_iJ_i − Σ_j J_j A_j Σ_i F_ji = Σ_i A_iJ_i − Σ_j A_jJ_j = **0**.

(b) Q_i = A_i(J_i − G_i), and G_i = Σ_j F_ijJ_j (by reciprocity, A_iG_i = Σ_j A_jF_jiJ_j = A_iΣ_jF_ijJ_j). Known Q: J_i − Σ_jF_ijJ_j = Q_i/A_i. Known T: J_i = ε_iσT_i⁴ + (1 − ε_i)G_i ⇒ **J_i − (1 − ε_i)Σ_jF_ijJ_j = ε_iσT_i⁴**. (Eliminating G another way gives the surface resistance form Q_i = (E_b,i − J_i)ε_iA_i/(1 − ε_i).)

(c) Q = (E_b1 − J₁)/R₁ = (J₁ − J₂)/R₁₂ = (J₂ − E_b2)/R₂ with R₁ = (1 − ε₁)/(ε₁A₁), R₁₂ = 1/(A₁F₁₂), R₂ = (1 − ε₂)/(ε₂A₂). Adding: **Q = σ(T₁⁴ − T₂⁴)/(R₁ + R₁₂ + R₂)**. (i) F = 1, A₁ = A₂ = A: q = σ(T₁⁴ − T₂⁴)/(1/ε₁ + 1/ε₂ − 1). (ii) F₁₂ = 1, A₁/A₂ → 0: Q = ε₁A₁σ(T₁⁴ − T₂⁴): a small body in a large room exchanges as if the room were black, whatever its emittance, because the room's reflections almost never return to the body.

(d) Each of the n + 1 gaps is a parallel-plate pair. The gap between plate 1 and the first shield has resistance (per unit area) 1/ε₁ + 1/ε_s − 1; gaps between shields 2/ε_s − 1; the last 1/ε_s + 1/ε₂ − 1. In series the total is (1/ε₁ + 1/ε₂ − 1) + n(2/ε_s − 1), and the same flux passes all gaps: **q = σ(T₁⁴ − T₂⁴)/[(1/ε₁ + 1/ε₂ − 1) + n(2/ε_s − 1)]**.

**Common mistakes.** Writing G with F_ji instead of F_ij without reciprocity; counting a shield's two faces as one resistance.

## 3. A thermocouple that reads low (s = 15: T_g = 575 K, T_w = 425 K)

**Assumptions.** Steady; the bead small against the duct (small-body limit); no conduction along the leads.

(a) 100(575 − T_b) = 0.8σ(T_b⁴ − 425⁴) ⇒ **T_b = 548.69 K, an error of 26.31 K** (reads low).

(b) (i) ε = 0.2: T_b = 566.98 K, **error 8.02 K**. (ii) h = 200 W/(m² K): T_b = 560.08 K, **error 14.92 K**. The error scales roughly as εσ(T_b⁴ − T_w⁴)/h.

(c) Shield: 2 × 80 × (575 − T_s) = 0.1σ(T_s⁴ − 425⁴) ⇒ **T_s = 572.35 K**. Bead seeing the shield: 100(575 − T_b) = 0.8σ(T_b⁴ − 572.35⁴) ⇒ **T_b = 574.33 K, error 0.67 K**.

(d) The shield wins by a factor of 40 over the bare bead, because the bead now sees a surface only 2.7 K below the gas. In practice the shield can foul (its emittance rises), slows the probe's response, and may block the flow if it is too closed. A low-ε coating degrades at temperature; a higher h needs a smaller bead (fragile) or a faster flow (not under your control). Combine a shield with a small bead.

**Check.** If T_w = T_g, every version reads T_g exactly: radiation error vanishes when the surroundings are at the gas temperature.

**Common mistakes.** Using T_g − T_w in a linear radiation term; forgetting that the shield convects on two faces.

## 4. A tilted solar array (s = 15: φ = 130°)

Radiator o = 0, u = (0, 1, 0), v = (0.6, 0, 0); array o = 0, u = (0, 1, 0), v = (cos 130°, 0, sin 130°).

(a)

| n | F (radiator → array) | A_r F_ra | A_a F_ar |
|---|---|---|---|
| 10 | 0.08530 | 0.051181 | 0.051181 |
| 20 | 0.08371 | 0.050229 | 0.050229 |
| 40 | 0.08286 | 0.049717 | 0.049717 |

Reciprocity holds at every n to the digits shown (the discrete kernel is symmetric).

(b) First-order Richardson: F_∞ ≈ 2F₄₀ − F₂₀ = **0.0820**. Check at 90°: n = 40 gives 0.27346, the extrapolation 0.26894, the closed form **0.26896**: the extrapolation removes 98 % of the n = 40 error.

(c) With F = 0.0820, the radiator runs at **17.46 °C**, against **35.13 °C** at 90°. Opening the array by 40° cuts the radiator's view of it by a factor of 3.3 and recovers most of the open-space performance (8.8 °C). A design rule follows: keep warm deployables out of the radiator's hemisphere, or cant them away.

**Common mistakes.** Using φ measured from the normal instead of from the radiator's plane; extrapolating with order 2.

## 5. Why a hole looks black (s = 15: ε_w = 0.4)

(a) Cube with one open face: **ε_eff = 0.7692** for ε_w = 0.4.

(b)

| cavity | ε_w = 0.1 | ε_w = 0.4 | ε_w = 0.9 |
|---|---|---|---|
| 1 × 1 × 1 cube | 0.3571 | 0.7692 | 0.9783 |
| 1 × 1 × 3, open end | 0.5903 | 0.8956 | — |

(c) Radiation entering or leaving the opening is reflected many times inside; each reflection gives the walls another chance to absorb (or emit), so the opening behaves as a surface far blacker than its walls, and more so the deeper the cavity (the opening sees mostly walls, which see mostly each other). A black-body calibration source is a deep, isothermal cavity with a small aperture for exactly this reason. A spacecraft's deep open bay or a sunshade emits like a near-black surface to space even if its walls are shiny, and also traps sunlight that enters it.

**Check.** Cube with one open face and ε_w = 1: ε_eff = 1 (the opening sees only black walls at 400 K).

## 6. Review a radiator design

(a) The linearisation about T_sink = 3 K gives 4εσT_sink³(T − T_sink) = **0.0017 W/m²** at 20 °C, against the exact **389.45 W/m²**: wrong by five orders of magnitude. Linearise about the mean of the two temperatures, or (better for a radiator facing space) use εσ(T⁴ − T_sink⁴) directly.

(b) Kirchhoff's law holds wavelength by wavelength. Sunlight is 98 % below 3 μm and the paint's absorptance there is α_s = 0.12–0.18 (Table 10.5); 0.93 is its infrared emittance. With α = 0.93 the panel would absorb 1265.7 W/m² in normal sun, more than the 389.5 W/m² it can emit at 20 °C: it would seem impossible to cool at all.

(c) At end of life with α_s = 0.30 the absorbed flux (408 W/m²) exceeds the emission at 20 °C, so no area works; the 1.079 m² panel settles at **52.6 °C** in normal sun. Size on end-of-life α_s, or (as real designs do) never let the Sun fall on the radiator.

(d) The array's view factor is 0.27–0.29, and a 330 K array face raised the 0.6 m² radiator from 8.8 to 35.1 °C; holding 20 °C needed 0.759 m² instead of 0.514 m². Neglecting it under-sizes the radiator by about a third.

(e) The ideal 20-shield stack of the lesson gives 0.536 W/m² (one bare wall would lose 376.9 W/m²), but that ideal ignores conduction between layers, seams and penetrations. Real blankets are often characterised by an effective emittance; with an illustrative ε_eff = 0.03 (not from a source in this course; untested), the loss would be **12.6 W/m²**, 23 times the ideal. Design heaters from tested blanket data, not from the ideal shield formula.

## Check-yourself answers

1. F₁₂ = (1/A₁)∫∫ cos θ₁ cos θ₂/(πs²) dA₂dA₁; multiplying by A₁ gives an integrand symmetric in 1 and 2, so A₁F₁₂ = A₂F₂₁.
2. Each row sums to one: everything leaving a face lands on some face of the closed box (a flat face sees none of itself).
3. At about 9.7 μm (λT = 2898 μm K). Sunlight is concentrated below 3 μm, so a paint can reflect visible light (low α_s) and still emit strongly in the far infrared (high ε).
4. Known T: J_i − (1 − ε_i)Σ_jF_ijJ_j = ε_iσT_i⁴. Known Q: J_i − Σ_jF_ijJ_j = Q_i/A_i.
5. With Q = 0 its radiosity equals its irradiation, J = G, whatever its emittance; its surface resistance carries no current, so ε drops out of the network.
6. 1/(1/0.8 + 1/0.6 − 1) = 1/1.9167 = 0.522.
7. The midpoint rule converges at order 2 for a smooth integrand; along a shared edge s → 0 and the kernel is singular, so the cells next to the edge carry an error that only halves with each halving of the cell size.
8. A = Q/(εσT⁴) = 200/(0.93 × 5.67 × 10⁻⁸ × 293.15⁴) ≈ 200/389.5 ≈ 0.51 m².
