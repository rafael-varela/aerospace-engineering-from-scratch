# P10.14  Worked solutions

Each solution states: assumptions, units, the derivation or method, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15; numerical results come from `solution.py` and probe scripts that call it, given to the figures shown. The bubble length rests on the stated empirical correlation (N_ℓ = 5 × 10⁴, uncertain by a factor of about 1.5). Values quoted from the Ladson file are authentic measurements (NASA TM-4074); everything else is computed.

## 1. Your own small flyer (s = 15)

(a) Reference UAV, 12 % section, all-laminar floor:

| V (m/s) | Re | cd floor |
|---|---|---|
| 8 | 1.095 × 10⁵ | 0.01005 |
| 12 | 1.643 × 10⁵ | 0.00821 |
| 18 | 2.464 × 10⁵ | 0.00670 |

The floor falls as Re^{−1/2}: a 2.25-fold increase in Re lowers it by 1.5.

(b) V = 10 m/s, c = 0.125 m: **Re = 8.56 × 10⁴**, far below the 0012's natural-transition Re of 4.2 × 10⁵ at cl = 0.445: the laminar layer will separate first and the bubble decides the performance. It is also only 1.55 times the 0012's Re_min of 5.5 × 10⁴, a slim margin given the correlation's uncertainty.

**Common mistakes.** Using the span instead of the chord; presenting the floor as the section's drag.

## 2. Derivations

(a) With Ue = Cs^m: ∫₀ˢ C⁵s^{5m} ds = C⁵s^{5m+1}/(5m + 1), so θ² = 0.45ν C⁵s^{5m+1}/((5m + 1)C⁶s^{6m}) = **0.45νs^{1−m}/(C(5m + 1))** (C = 1 in the lesson). Ue′ = mCs^{m−1}, so **λ = θ²Ue′/ν = 0.45m/(5m + 1)**, constant along the surface. λ = 0.075 at **m = 1** (the stagnation point). λ = −0.09 at 0.45m = −0.09(5m + 1), i.e. **m = −0.1**: a retarded power-law flow with m ≤ −0.1 separates at once in Thwaites' model (the exact Falkner–Skan separation value is m ≈ −0.090).

(b) θ² = ν·F(s) with F(s) = 0.45 Ue⁻⁶∫Ue⁵ds, which depends only on the shape of Ue(s). Then λ = F·Ue′ does not contain ν, so the station where λ = −0.09 does not move with ν. θ there is √(νF) ∝ √ν ∝ Re^{−1/2} at fixed V and c.

(c) The bubble closes on the airfoil if s_s + ℓ ≤ s_TE with ℓ = N_ℓν/Ue_s. In chord units (lengths divided by c, speeds by V), ν/(Vc) = 1/Re, so **ℓ/c = N_ℓ/(Re·Ue_s/V)**. Closure is lost when ℓ/c = (s_TE − s_s)/c, i.e. at **Re_min = N_ℓ/((Ue_s/V)(s_TE − s_s)/c)**.

(d) On the first interval [0, h] after a stagnation point, Ue = as, and the exact integral is ∫₀ʰa⁵s⁵ds = a⁵h⁶/6; the trapezoid rule gives h(0 + a⁵h⁵)/2 = a⁵h⁶/2, **three times too large**, and the error propagates into θ for the next few intervals, so λ near the stagnation point is wrong (the stagnation-point test caught it). For Ue linear from u₀ to u₁ over length h: ∫u⁵ds = h∫₀¹(u₀ + (u₁ − u₀)t)⁵dt = **h(u₁⁶ − u₀⁶)/(6(u₁ − u₀))**, exact for piecewise-linear Ue (use h u₀⁵ when u₁ = u₀).

**Common mistakes.** Treating λ as dimensional; forgetting that separation needs a falling Ue, not just a small one.

## 3. A rotor blade on Mars

(a) **ν = 1.1 × 10⁻⁵/0.017 = 6.47 × 10⁻⁴ m²/s** (44 times sea level). **Re = 120 × 0.08/6.47 × 10⁻⁴ = 1.48 × 10⁴**, between a bumblebee and a pigeon in the table of §2.1, despite a speed of 120 m/s.

(b) All-laminar floor for a 2 % plate: **cd = 0.0227**.

(c) At Re ≈ 1.5 × 10⁴ no bubble of the lesson's model closes on either section (the 0012's Re_min is 5.5 × 10⁴ at cl = 0.445 and the bubble would be about 2.7c long), so the flow over the aft part will be separated whatever the section. The choice is then the one of §2.5: a thin, cambered plate with a sharp leading edge, operated near its design angle, whose separation (if any) leaves little pressure recovery unmade, and whose low thickness keeps the friction and pressure drag low. The lesson's models can rank the choices but do not predict the section's performance at this Re; measured low-Re data would be needed.

**Common mistakes.** Using Earth's ν; concluding from the high speed that the flow is turbulent.

## 4. Slow flight at cl = 0.80, Re = 1.64 × 10⁵

(a)

| section | α (°) | peak Ue/V | x_s/c | ℓ/c | closes | ΔCp_rec | Re_min |
|---|---|---|---|---|---|---|---|
| 0012 | 6.63 | 2.033 | 0.029 | 0.167 | yes | 2.90 | 2.80 × 10⁴ |
| 2412 | 4.48 | 1.629 | 0.315 | 0.223 | yes | 1.43 | 5.29 × 10⁴ |
| 4402 | 3.01 | 3.131 | 0.0001 | 0.106 | yes | 7.42 | 1.73 × 10⁴ |

(b) At cl = 0.80 the 4402 is about 3° above its design angle. Its 2 % thick, sharp leading edge cannot turn the flow smoothly: the upper surface shows a narrow spike to 3.13V right at the leading edge followed by a steep deceleration, and the laminar layer separates at the leading edge itself (x/c = 0.0001). The model says the resulting short bubble closes quickly (ℓ = 0.11c), but the pressure recovery left after it, ΔCp = 7.4, is enormous; real thin plates in this condition have a leading-edge bubble that easily bursts into thin-airfoil stall. Thin plates are excellent near their design cl and intolerant away from it.

(c) Natural-transition Re (Michel reaching the separation point):

| cl | 0012 | 4402 |
|---|---|---|
| 0.2 | 3.75 × 10⁵ | 1.41 × 10⁶ |
| 0.445 | 4.20 × 10⁵ | 6.09 × 10⁵ |
| 0.8 | 5.5 × 10⁷ | none in 10³–10⁹ |

At higher cl the separation point moves to the leading-edge suction peak, where the boundary layer is too thin for any reasonable Re to transition first: a leading-edge bubble is then present at all practical Reynolds numbers, as on real airfoils near stall.

**Common mistakes.** Comparing sections at the same α instead of the same cl.

## 5. What grit does (authentic data)

(a) cl_max, free transition against wraparound grit (Ladson, NASA TM-4074): Re ≈ 2 × 10⁶: 1.408 → 1.165 (**−17.3 %**); 4 × 10⁶: 1.484 → 1.207 (**−18.7 %**); 6 × 10⁶: 1.660 → 1.217 (**−26.7 %**); 9 × 10⁶: 1.693 → 1.229 (**−27.4 %**); 12 × 10⁶: 1.617 → 1.121 (**−30.7 %**).

(b) cd near α = 0 with grit against the fully turbulent model: 0.01100/0.00987, 0.00954/0.00875, 0.00895/0.00818, 0.00862/0.00766, 0.00842/0.00730 at Re = 2.0, 4.0, 6.0, 8.95, 12.1 × 10⁶: the grit sections have 9–15 % more drag than a smooth turbulent plate model.

(c) At Re = 2 × 10⁶ the free-transition boundary layer would have transitioned on its own before separating; the grit only adds roughness drag and a thicker, older turbulent layer that separates earlier near stall, costing cl_max. At Re = 10⁵ the laminar layer separates before it can transition (§2.3) and forms a long bubble that may not close; a turbulator just ahead of the separation point makes the layer turbulent before the adverse gradient, removing the bubble. Its own drag is smaller than the bubble's penalty, so it helps.

**Common mistakes.** Generalising "roughness is bad" across Reynolds-number regimes.

## 6. How far can you trust the bubble model?

(a)

| N_ℓ | Re_min, 0012 | Re_min, 4402 |
|---|---|---|
| 3 × 10⁴ | 3.31 × 10⁴ | 1.10 × 10⁵ |
| 5 × 10⁴ | 5.51 × 10⁴ | 1.83 × 10⁵ |
| 7 × 10⁴ | 7.71 × 10⁴ | 2.57 × 10⁵ |

Over the whole range the separation points, the ranking of ΔCp_rec, and the conclusion that the UAV's 0012 bubble closes at cruise (Re = 1.64 × 10⁵ > 7.7 × 10⁴) survive. The Re at which the 0012's bubble opens is uncertain by ±40 %.

(b) (1) Free-stream turbulence and propeller wash: earlier transition, shorter bubble. (2) The bubble's displacement effect on Ue (a pressure plateau): it changes where the layer separates and usually lengthens the laminar part slightly. (3) Surface roughness and vibration of a real model wing: earlier transition, shorter bubble. (Three-dimensional and unsteady effects near the tips act in either direction.)

(c) Howarth: Thwaites separates 2.6 % late. On the 0012 at cruise s_s = 0.315 (arc from the stagnation point), so the error is about 0.008c. Re_min ∝ 1/(s_TE − s_s) with s_TE − s_s = 0.716: a 1 % change, negligible beside the ±40 % from N_ℓ.

(d) A surface-flow visualisation (oil film or tufts) or a chordwise pressure survey on the wing at cruise speed in a wind tunnel or in steady flight: it shows directly whether there is a bubble, where the layer separates and where it reattaches, testing x_s and ℓ together.

**Common mistakes.** Quoting Re_min to three figures without the correlation's uncertainty.

## Check-yourself answers

1. 1.64 × 10⁵; it lies in the band 10⁴–5 × 10⁵ where the laminar layer separates before transition and bubbles dominate.
2. θ² = (0.45ν/Ue⁶)∫₀ˢUe⁵ds; laminar separation where λ = θ²Ue′/ν = −0.09.
3. θ² is proportional to ν, so λ = θ²Ue′/ν does not contain ν: λ, and the place where it reaches −0.09, depend only on the shape of Ue(s).
4. A region where the separated laminar layer becomes turbulent in the air and reattaches; it closes if the transition length ℓ ≈ N_ℓν/Ue_s fits between the separation point and the trailing edge.
5. ℓ/c = N_ℓ/(Re·Ue_s/V), as 1/Re; Re_min = N_ℓ/((Ue_s/V)(s_TE − s_s)/c).
6. It separates late and its speed has risen and fallen gently, so the pressure recovery left after separation (ΔCp_rec = 0.36 against 1.16 for the 0012 at cl = 0.445) is small.
7. A bubble adds drag (and an open separation adds much more) on top of the laminar skin friction; the formula counts only friction.
8. The data cover 2 × 10⁶–1.9 × 10⁷, where transition precedes separation; below about 5 × 10⁵ the stall mechanism changes (bubbles, bursting), so the trend fitted at high Re does not describe it.
