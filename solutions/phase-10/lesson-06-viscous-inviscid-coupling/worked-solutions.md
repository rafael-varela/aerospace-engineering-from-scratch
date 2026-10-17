# P10.06  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Numbers come from `solution.py` and the probe scripts that call it, given to the figures shown. Measured values are Ladson's (NASA TM-4074, authentic, public domain) as transcribed in `curriculum/phase-10/data/ladson_naca0012_M015.csv`; every other number is computed.

## 1. The boundary layer on the trainer's wing

**Assumptions.** Flat plate at zero pressure gradient, 55 m/s, ISA ν; laminar from the leading edge.

(a) Re_x = 55 × 0.3/1.461 × 10⁻⁵ = **1.129 × 10⁶**. Blasius (from `falkner_skan(0)`: θ√Re_x/x = 0.66411, H = 2.5911, c_f√Re_x = 0.66411): **θ = 0.1875 mm, δ* = 0.4858 mm, c_f = 6.249 × 10⁻⁴**. Thwaites (λ = 0): **θ = 0.1894 mm, δ* = 0.4943 mm, c_f = 6.172 × 10⁻⁴**. Differences: θ +1.0 %, δ* +1.7 %, c_f −1.2 %. The whole boundary layer is about half a millimetre thick at 0.3 m.

(b) Michel's criterion on the Blasius layer is met at **Re_x = 2.02 × 10⁶**, that is **x = 0.54 m** at 55 m/s.

(c) The accelerating region near the nose stabilises the layer (Re_θ grows slowly and the critical value is far off), delaying transition; past the suction peak the deceleration destabilises it and advances transition. On the 2412 at cruise the net result is transition at 0.36c = 0.43 m on the upper surface (§4 of the lesson), ahead of the flat-plate 0.54 m.

**Common mistakes.** Using the chord Reynolds number instead of Re_x; quoting δ (the 99 % thickness) where δ* is meant.

## 2. Thwaites' method and the Squire–Young formula

(a) Multiply the momentum integral by 2u_eθ/ν: (u_e/ν) d(θ²)/ds + 2(2 + H)λ = c_f u_eθ/ν = 2ℓ, with ℓ = (θ/u_e) ∂u/∂y|_w and c_f = 2νℓ/(u_eθ). So (u_e/ν) d(θ²)/ds = 2[ℓ − (2 + H)λ] ≡ F(λ). Thwaites found that for many exact solutions F(λ) ≈ 0.45 − 6λ. Then u_e d(θ²)/ds = 0.45ν − 6θ²u_e′, which multiplied by u_e⁵ is d(θ²u_e⁶)/ds = 0.45νu_e⁵, so **θ²u_e⁶ = 0.45ν∫₀^s u_e⁵ ds**.

(b) With u_e = Cs: θ² = 0.45ν C⁵(s⁶/6)/(C⁶s⁶) = **0.075ν/C**, independent of s; λ = θ²C/ν = **0.075**. So θ = 0.27386√(ν/C), against Hiemenz's **0.29234√(ν/C)** (`falkner_skan(1)`): Thwaites is **6.3 % low** at a stagnation point (and its H = 2.358 against 2.216).

(c) With u_e,TE = V the factor is 1: **c_d = 2θ_TE/c**, the momentum deficit of two surfaces. When u_e,TE < V, the wake still accelerates to V downstream; its momentum thickness falls as it accelerates (by the momentum integral in the wake, with c_f = 0), so the far-wake θ, which is the drag, is smaller than θ_TE. The exponent (H + 5)/2 is that integral, done with an assumed relation between H and u_e in the wake.

## 3. The small UAV

(a) Re = 12 × 0.20/1.461 × 10⁻⁵ = **1.64 × 10⁵**.

(b) α = 0°: transition by **laminar separation** at x/c = 0.579 (upper) and 0.712 (lower); **c_d = 0.0118**. α = 4°: laminar separation at 0.345 (upper) and 0.874 (lower), turbulent separation predicted on the upper surface just before the trailing edge; **c_d = 0.0142**. Michel's criterion is never reached first.

(c) Not far. The method turns a laminar separation into an attached turbulent layer instantly. At Re ≈ 1.6 × 10⁵ the separated shear layer takes a sizeable fraction of the chord to go turbulent and reattach (a laminar separation bubble), adding pressure drag the model omits, and at higher α it may not reattach at all. The prediction is likely to be **low**, and its uncertainty is tens of per cent; `P10.14` treats this regime.

## 4. NACA 0012 against Ladson

(a) Table I, free transition, |α| ≤ 0.15° (measured), against the prediction with 200 panels:

| Re | measured c_d | predicted | difference |
|---|---|---|---|
| 2.00 × 10⁶ | 0.00620 | 0.00690 | +11.2 % |
| 3.94 × 10⁶ | 0.00659 | 0.00649 | −1.5 % |
| 5.97 × 10⁶ | 0.00650 | 0.00617 | −5.1 % |
| 8.86 × 10⁶ | 0.00630 | 0.00590 | −6.3 % |
| 1.19 × 10⁷ | 0.00635 | 0.00567 | −10.8 % |
| 1.89 × 10⁷ | 0.00660 | 0.00536 | −18.7 % |

(b) Table VII (wraparound grit), against a forced transition at 5 % chord:

| Re | measured c_d | predicted | difference |
|---|---|---|---|
| 2.00 × 10⁶ | 0.01096 | 0.00974 | −11.1 % |
| 4.00 × 10⁶ | 0.00953 | 0.00860 | −9.8 % |
| 6.00 × 10⁶ | 0.00895 | 0.00799 | −10.7 % |
| 8.95 × 10⁶ | 0.00855 | 0.00743 | −13.0 % |
| 1.21 × 10⁷ | 0.00835 | 0.00704 | −15.7 % |

(Measured values at repeated α ≈ 0 points are averaged.) The fixed-transition agreement is worse and uniformly low. The model's trip is a clean switch at 5 % chord; the real strip is No. 60 grit wrapping round the nose from 5 % upper to 5 % lower, so the boundary layer is turbulent from the leading edge, starts thicker, and carries the grit's own roughness drag. Both make the measured drag higher than a "transition at 5 %" model.

(c) Measured free-transition drag is nearly flat at 0.0062–0.0066 across a tenfold range of Re, while the prediction falls steadily; the grit data fall with Re, roughly parallel to the tripped prediction and 10–16 % above it.

## 5. The trainer's section polar (Re = 4.52 × 10⁶)

(a)

| α | c_l (inviscid) | c_d | x_tr upper | x_tr lower | turbulent separation |
|---|---|---|---|---|---|
| 0° | 0.260 | 0.00633 | 0.391 (Michel) | 0.222 | none |
| 2° | 0.501 | 0.00664 | 0.257 | 0.422 | none |
| 4° | 0.742 | 0.00731 | 0.132 | 0.697 | none |
| 6° | 0.981 | 0.00828 | 0.066 | 0.904 | none |
| 8° | 1.220 | 0.01029 | 0.022 (laminar separation) | 0.930 | upper surface |

(b) One displacement-thickness iteration reduces c_l by **1.7 %, 3.1 %, 4.8 % and 7.9 %** at 0°, 2°, 4° and 6°. As α grows the upper-surface boundary layer thickens (earlier transition, stronger adverse gradient) while the lower one thins, so the displacement body is effectively less cambered and thinner at the trailing edge: the viscous decambering grows with load.

(c) At cruise (α = 0.59°): c_d = 0.00637 clean and 0.00866 tripped at 5 % chord. With q = 1852.8 Pa and S = 12 m²: **D = 142 N clean, 193 N tripped**, a difference of 51 N, **2.8 kW** at 55 m/s. On a 750 kg trainer that is a noticeable share of the cruise power, which is why laminar-flow surfaces are kept clean.

## 6. Why does the prediction fall below the data at high Re?

(a) At Re = 1.19 × 10⁷: free transition (24.5 %) gives 0.00567; trips at 20 %, 10 % and 5 % give 0.00596, 0.00670 and 0.00706. The measured 0.00635 lies between the 10 % and 20 % trips: the model would need **transition near 14 % chord**, well ahead of Michel's 24.5 %.

(b) (1) **Free-stream disturbances**: tunnel turbulence and noise, even in a low-turbulence tunnel, advance transition relative to a criterion correlated on free-flight-like data; (2) **surface condition**: as Re rises, the same polished surface is relatively rougher (the roughness height is fixed while the boundary layer gets thinner), and transition moves forward. Even with the right transition point the model may be low because Head's closure and the Ludwieg–Tillmann friction law under-predict the turbulent momentum growth at high Re, and Squire–Young's extrapolation is approximate.

(c) The first improvement is the **transition model**: replace Michel's local correlation with an amplification (e^N) method that accounts for the pressure gradient history and the free-stream turbulence level (XFOIL's approach, `P10.21`). To validate it you need **measured transition locations** (surface hot films, flow visualisation) alongside the drag, at several Re.

## Check-yourself answers

1. θ = ∫(u/u_e)(1 − u/u_e) dy, the momentum thickness: the thickness of free-stream fluid whose momentum flux equals the deficit; δ* is the displacement thickness, H = δ*/θ the shape factor.
2. θ²u_e⁶ = 0.45ν∫u_e⁵ ds with λ = θ²u_e′/ν; laminar separation when λ < −0.09.
3. Re_θ ≥ 1.174(1 + 22 400/Re_x)Re_x^0.46; on a flat plate (Blasius layer) it is met at Re_x ≈ 2.0 × 10⁶.
4. Head's method integrates momentum and entrainment equations for θ and the entrainment shape factor H₁, closed by Ludwieg–Tillmann skin friction; turbulent separation is taken at H ≈ 2.4.
5. c_d = 2(θ_TE/c)(u_e,TE/V)^{(H_TE + 5)/2}, summed over both surfaces.
6. 0.6708 by Thwaites (√0.45) against 0.6641 by Blasius: 1.0 % high.
7. From 2.0 × 10⁶ to 8.86 × 10⁶ within 11 %; at 1.89 × 10⁷ the model is 19 % low because its transition point is too far aft.
