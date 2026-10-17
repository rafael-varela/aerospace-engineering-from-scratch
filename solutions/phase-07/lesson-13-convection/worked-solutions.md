# P07.13  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe script that calls it, given to the figures shown. Properties are from Lienhard Tables A.6 and A.3 at the film temperature.

## 1. A heated probe (s = 15: D = 7.5 mm, U = 35 m/s)

**Assumptions.** Steady, uniform surface temperature, air at 1 atm (at altitude the pressure is lower; Problem 3 shows how that enters), no radiation (small at these temperatures against an h of hundreds).

(a) T_f = (330 + 260)/2 = **295 K**: ν = 1.5285 × 10⁻⁵ m²/s, k = 0.02600 W/(m K), Pr = 0.7075. Re_D = 35 × 0.0075/1.5285 × 10⁻⁵ = **17 174**.

(b) Churchill–Bernstein: **Nu = 72.56**. The simplified form (7.66), which drops the high-Re bracket, gives **63.86**. At Re = 17 174 the bracket [1 + (Re/282 000)^(5/8)]^(4/5) is 1.137, not 1: the simplified form is meant for Re below about 4000. Use the full form.

(c) h = 72.56 × 0.02600/0.0075 = **251.5 W/(m² K)**; q′ = h πD (T_w − T∞) = **414.9 W/m**.

(d) At 70 m/s, Nu = 109.00: a factor **1.502**, more than √2 = 1.414. The bracket grows with Re, so the local exponent d ln Nu/d ln Re is above 1/2 here: the wake region contributes more as Re rises.

**Independent check.** The flat-plate laminar law with L = πD/2 (the run from the stagnation point to separation) gives the right order: 0.664 √(Re × π/2) Pr^(1/3) = 97, the same order as 72.6 (the rear of a cylinder transfers less than a plate).

**Common mistakes.** Using T∞ or T_w properties (small effect here); using (7.66) at Re = 17 000; forgetting the π in the perimeter.

## 2. Derivations

(a) h_x = (k/x) 0.332 (U x/ν)^(1/2) Pr^(1/3) = 0.332 k Pr^(1/3) (U/ν)^(1/2) x^(−1/2). Then h̄ = (1/L)∫₀ᴸ h_x dx = 0.332 k Pr^(1/3) (U/ν)^(1/2) × 2L^(1/2)/L, so **Nu_L = h̄L/k = 0.664 Re_L^(1/2) Pr^(1/3) = 2 Nu_x(L)**.

(b) St_x Pr^(2/3) = Nu_x/(Re_x Pr^(1/3)) = 0.332 Re_x^(−1/2) = c_f/2 ✓. Turbulent: St_x Pr^(2/3) = 0.0296 Re_x^(−0.2) Pr^(0.6 − 1/3), and c_f/2 = 0.0296 Re_x^(−0.2), so the ratio is **Pr^(0.6 − 1/3) = Pr^0.2667**: **0.9117 for air** (Pr 0.707) and **1.602 for water** (Pr 5.85). The power laws were fitted to gases; for water the analogy with Pr^(2/3) and this Nu law disagree by 60 %, a sign that (6.112) is not meant for liquids.

(c) With x_c = ν Re_c/U,
h̄L = ∫₀^(x_c) 0.332 k Pr^(1/3) (U/ν)^(1/2) x^(−1/2) dx + ∫_(x_c)^L 0.0296 k Pr^0.6 (U/ν)^0.8 x^(−0.2) dx
= 0.664 k Pr^(1/3) Re_c^(1/2) + (0.0296/0.8) k Pr^0.6 (Re_L^0.8 − Re_c^0.8).
Dividing by k: **Nu_L = 0.037 Pr^0.6 (Re_L^0.8 − Re_c^0.8) + 0.664 Re_c^(1/2) Pr^(1/3)**. At Re_L = Re_c the first term vanishes and the second is the laminar average at Re_c: continuous ✓ (the tests check it to 10⁻⁷).

(d) Nu_L = 0.037 Re_L^0.8 Pr^0.6 against Nu_x(L) = 0.0296 Re_L^0.8 Pr^0.6: ratio **1.25** = 1/0.8. Turbulent h falls only as x^(−0.2), so the leading edge contributes less extra than in laminar flow (x^(−1/2), ratio 2).

**Common mistakes.** Integrating Nu_x instead of h_x (Nu_x contains x); forgetting that Re_c is fixed while Re_L varies.

## 3. An unpressurised equipment bay

(a) ρ = p/(RT) ∝ p at fixed T; μ, k, c_p do not depend on p (ideal gas), so ν = μ/ρ ∝ 1/p and α = k/(ρc_p) ∝ 1/p, and Ra = gβΔTL³/(να) ∝ p². Laminar Nu − c ∝ Ra^(1/4) ∝ p^(1/2), so h_conv ∝ p^(1/2) roughly (the constant term weakens it). At 19 330 Pa, **(p/p₀)^(1/2) = 0.437**.

(b) With bay air and walls at 25 °C: **51.76 °C at 1 atm, 58.79 °C at 19 330 Pa**. h_conv falls from 4.807 to 2.301 W/(m² K) (a ratio of 0.479, close to the predicted 0.437; the constant terms of the correlations and the upper face's thick-boundary-layer correction soften the fall).

(c) Bay at 250 K: **9.30 °C at 1 atm, 20.91 °C at 19 330 Pa**. Radiation share: 56.3 % and 73.6 % at 25 °C; 42.7 % and 62.0 % at 250 K. At altitude, radiation does most of the work, and the paint matters even more than at sea level.

(d) That the box radiates to large surroundings at the air temperature. A real bay has neighbouring boxes at similar temperatures and structure at different ones; neighbours that are as warm as the box block radiation (lesson `P07.14`), which would raise the temperature, possibly by more than the pressure effect.

**Check.** Radiation is unchanged by pressure (h_rad 6.18 → 6.40 W/(m² K) only because the box is hotter), as it must be.

**Common mistakes.** Scaling k or μ with pressure; forgetting to restore the original property function.

## 4. Designing the cold plate (s = 15: Q = 450 W)

Assumptions as in the lesson: uniform flux, fully developed h (conservative), smooth tube, outlet wall hottest.

| D (mm) | case | ṁ (g/s) | Re | T_out (°C) | hottest wall (°C) | U (m/s) | Δp (Pa) | pumping (W) |
|---|---|---|---|---|---|---|---|---|
| 3 | minimum for 40 °C | 8.04 | 3980 | 33.39 | 40.00 | 1.141 | 10 751 | 0.087 |
| 3 | Re = 5000 | 10.48 | 5000 | 30.27 | 35.38 | 1.487 | 17 003 | 0.179 |
| 4 | minimum for 40 °C | 9.08 | 3310 | 31.85 | 40.00 | 0.725 | 3458 | 0.032 |
| 4 | Re = 5000 | 14.43 | 5000 | 27.46 | 32.52 | 1.152 | 7651 | 0.111 |
| 5 | minimum for 40 °C | 10.20 | 2930 | 30.55 | 40.00 | 0.521 | 1490 | 0.015 |
| 5 | Re = 5000 | 18.37 | 5000 | 25.86 | 30.89 | 0.938 | 4059 | 0.075 |

(d) A defensible choice is **D = 4 mm at about 14.4 g/s** (Re 5000, wall 32.5 °C, 7.5 K of margin, Δp 7.7 kPa, 0.11 W of pumping). The minimum flows all sit at Re 2900–4000, where Gnielinski is least reliable and the flow could relaminarise; the 3 mm channel needs twice the pressure for the same Re; the 5 mm channel needs more flow. Pumping power is negligible in every case compared with 450 W; the pressure a small pump can supply is the real constraint.

**Check.** First law: 0.01443 × 4180 × (27.46 − 20) ≈ 450 W ✓.

**Common mistakes.** Taking the laminar branch's root (there is none below 40 °C); using the inlet temperature for the wall; forgetting that Re depends on D at fixed ṁ as 1/D.

## 5. The analogy in a wind tunnel

(a) Re_L = 10 × 0.40/1.575 × 10⁻⁵ = **253 968**, below 5 × 10⁵: laminar if the tunnel is quiet and the leading edge sharp.

(b) C_F = 0.0383/(½ × 1.177 × 10² × 2 × 0.12) = **2.712 × 10⁻³**; laminar 1.328/√Re = **2.635 × 10⁻³** (a drag of 0.03722 N); ratio **1.029**: the measurement is 2.9 % high, inside what an edge effect or balance offset could explain.

(c) St̄ = C_F/(2 Pr^(2/3)) = 1.708 × 10⁻³, h̄ = St̄ ρ c_p U = **20.23 W/(m² K)**; the correlation gives **19.68 W/(m² K)**. The analogy turns a 2.9 % drag excess into a 2.9 % h excess.

(d) Tripped: C_F = 0.074 Re^(−0.2) = 6.14 × 10⁻³, drag **0.0867 N**, Nu_L = 633.4, **h̄ = 41.8 W/(m² K)** (against laminar Nu 298.1). The balance can tell: a tripped plate would read 2.3 times the laminar drag, far outside a few per cent.

**Common mistakes.** Using one face's area; applying the local analogy to averages without checking that the averaging is the same (for the laminar plate both average by the same factor of two, so it holds).

## 6. Review a thermal analysis

(a) Gnielinski is valid only for **2300 ≤ Re ≤ 5 × 10⁶**; at Re = 1800 the flow is laminar and fully developed Nu = **4.364** for uniform flux. The colleague's 9.5 over-predicts h by a factor of **2.2** (about 650 W/(m² K), not 1420). Use the laminar value, add an entry-length correction if the tube is short, or raise the flow above Re ≈ 5000.

(b) (8.13a) is the laminar form, recommended for Ra ≤ 10⁹; at 3 × 10⁹ it gives **121.0** while the all-Ra form (8.13b) gives **172.8**. The laminar form under-predicts h by **30 %** here. Use (8.13b).

(c) The bottom of a hot box is the **stable** case: warm air is trapped under it. With Fujii–Imura (L = width 0.12 m, Ra = 3.42 × 10⁶, inside its 10⁶–10¹¹ range), h = **2.66 W/(m² K)**, less than half the top's **5.84 W/(m² K)**.

(d) For the black box at 60 W, radiation carries 56 % of the heat: neglecting it raises the predicted surface temperature from **51.8 to 77.5 °C**. Near room temperature h_rad = 4εσT³ ≈ 6 W/(m² K) is comparable with natural convection; it is never "small" for a high-emittance surface in still air.

(e) Natural-convection correlations reproduce their data to several per cent at best and ±10–20 % in practice. With h_conv × 0.8 and × 1.2 the box runs at **53.9 and 49.9 °C**: report about ±2 K, not ±0.1 K.

## Check-yourself answers

1. Nu = hL/k (wall gradient in units of ΔT/L), Re = UL/ν (inertia/viscosity), Pr = ν/α (momentum/thermal diffusivity), Ra = gβΔTL³/(να) (buoyancy/diffusion), St = h/(ρc_pU) = Nu/(RePr). Pr is a property of the fluid alone.
2. Because the local h falls as x^(−1/2), whose integral from 0 to L is 2L^(1/2): the average of x^(−1/2) over [0, L] is twice its value at L.
3. St Pr^(2/3) = c_f/2. Exact for the laminar plate at uniform wall temperature, because Nu_x and c_f come from the same similarity solution; with uniform heat flux the thermal boundary layer grows differently and the constant becomes 0.4587, not 0.332.
4. At T_f = (50 + 20)/2 = 35 °C. Air's properties change by a few per cent over 30 K, and h depends on them weakly (powers of 1/4 to 1/2), so the box temperature moves by 0.2 K between ambient, film and wall properties.
5. In the turbulent range Nu ∝ Ra^(1/3) ∝ L, so h = Nu k/L does not contain L: the heat transfer is governed by local plumes, not by the plate's size.
6. Gnielinski is valid only above Re = 2300; at 1500 the flow is laminar, and fully developed Nu = 4.364 for uniform heat flux (3.657 for uniform wall temperature).
7. A black surface radiates with h_rad ≈ 6 W/(m² K) near 40 °C, more than the natural-convection h of about 5; bare aluminium (ε ≈ 0.09) radiates almost nothing, so the convection alone must carry the heat at a larger ΔT (77.5 against 51.8 °C for 60 W).
8. At 8.8 g/s the flow is at Re 3171, in the transitional range where Gnielinski is least reliable and the flow can relaminarise, so h could be half the predicted value; a design at Re ≥ 5000 (15 g/s, wall 31.6 °C) has margin on both h and temperature.
