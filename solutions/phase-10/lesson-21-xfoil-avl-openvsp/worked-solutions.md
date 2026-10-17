# P10.21  Worked solutions

Each solution states assumptions, units, the derivation, the result, an independent check, and the common mistakes. Seeded problems are worked for s = 15. Numbers come from `solution.py` and a probe script that calls it. No XFOIL, AVL or OpenVSP output appears anywhere in these solutions: the tools were not available when the lesson was written. "Own code" means the compact panel method and one-way boundary-layer estimate given in the lesson; "published" means Ladson, NASA TM-4074 (`curriculum/phase-10/data/ladson_naca0012_M015.csv`).

## 1. What a polar file must say

(a) `read_tool_polar("naca0012.csv")` raises

```
ValueError: naca0012.csv: unlabelled polar, header lacks ['tool', 'version', 'settings']; record tool, version and settings
```

A header that passes:

```
# tool: XFOIL
# version: 6.99
# settings: re=6e6; mach=0.15; ncrit=9; panels=160; iter=100; xtr_top=1.0; xtr_bot=1.0
# airfoil: NACA 0012
# produced_by: <name>, <date>, <machine/OS>
# note: free transition; polar accumulated with PACC; unconverged points removed (list them)
```

(b) **Re**: sets the boundary-layer thickness and the transition location, so cd and cl_max change strongly. **Mach**: XFOIL applies a compressibility correction to the inviscid solution; 0.15 raises the lift slope by about 1 % over M = 0. **Ncrit**: the amplification exponent at which transition is predicted; a higher value delays transition and lowers drag. **Panel count**: changes the resolution of the leading-edge suction peak and the boundary-layer march; too few panels change cd by several per cent (Problem 5 shows 6 % for the own code). **Iteration limit**: too few iterations leave points unconverged; XFOIL reports them, and a polar that silently includes or drops them is not reproducible.

(c) XFOIL's polar file header records the version, the airfoil name, xtrf (forced transition, top and bottom), Mach, Re and Ncrit (`parse_xfoil_polar` reads exactly those). The panel count, the panel-spacing parameters and the iteration limit are not in it, so `xfoil_polar_to_labelled_csv` requires them as arguments, with who ran it.

**Common mistakes.** Recording Re but not Ncrit; editing a polar file by hand to remove "bad" points without saying so.

## 2. Comparison metrics

(a) α_L0 = −c/m. To first order, u²(α_L0) = (1/m²)[u²(c) + (c/m)² u²(m) + 2(c/m) Cov(m, c)] (signs: ∂α_L0/∂c = −1/m, ∂α_L0/∂m = c/m²). If the fitted α range lies far from α_L0, the intercept is an extrapolation: u(c) is large and strongly correlated with u(m), and α_L0 is poorly determined. A symmetric range around α_L0 decorrelates them.

(b) Comparing cl at different α mixes slope and offset errors with the sampling. Linear interpolation between published points has an error of order (Δα²/8)|cl″|, negligible where cl is linear in α, but near stall cl″ is large and the stall angle itself may fall between samples, so an interpolated difference there measures the interpolation, not the codes.

(c) If cl₂ = (1 + ε)mα and cl₁ = mα, the difference is εmα and RMS = εm√(⟨α²⟩). For a range symmetric about zero with half-width A and uniform sampling, √⟨α²⟩ = A/√3: the RMS grows linearly with A. An RMS difference without its α range is therefore meaningless; the lesson always quotes the range (−4.5° to 8.5°) and the slope difference separately.

## 3. A suspicious drag value

(a) Published at α = 0: free transition (Table I, Re = 5.97 × 10⁶) **cd = 0.0065**; wraparound grit (Table VII, Re = 6.00 × 10⁶) **cd = 0.00895**. 0.0045 is 31 % below the free-transition value of a polished model in a low-turbulence tunnel; for a real wing with rivets, joints and insect debris the flow transitions earlier and the grit-strip value is closer. 0.0045 is not plausible for either.

(b) A high Ncrit (laminar flow assumed to persist longer than in the tunnel, let alone in flight); no forced transition where the real wing has roughness; a different Re than stated; an unconverged point; too few panels; reading cdp (pressure drag only) instead of cd.

(c) The header of Problem 1 with the actual values, and: "Section drag from XFOIL with free transition at Ncrit = …; the published free-transition value at this Re is 0.0065 (Ladson, TM-4074), and roughness on a production wing raises it toward the fixed-transition value 0.009; the sizing uses 0.009."

**Common mistakes.** Treating a single tool run as data; quoting four significant figures without the settings.

## 4. Your own Reynolds number (s = 15: index 3, Re = 8.86 × 10⁶)

(a)

| | published (Table I) | own code | difference |
|---|---|---|---|
| cl_α [/deg] | 0.10866 | 0.12050 | +10.9 % |
| α_L0 [deg] | −0.024 | +0.001 | +0.025° |
| cd_min | 0.0063 | 0.00587 | −6.8 % |
| RMS cl difference, −4.5°…8.5° | | 0.0524 (max 0.090, 7 points) | |

(b) cd(own)/cd(published):

| α [deg] | −4.00 | −1.98 | −0.10 | 2.13 | 4.13 | 6.08 | 8.02 |
|---|---|---|---|---|---|---|---|
| published | 0.0068 | 0.0063 | 0.0063 | 0.0065 | 0.0070 | 0.0072 | 0.0085 |
| own | 0.00646 | 0.00610 | 0.00587 | 0.00608 | 0.00646 | 0.00703 | 0.00837 |
| ratio | 0.950 | 0.969 | 0.932 | 0.935 | 0.923 | 0.977 | 0.985 |

The own code is low by 2–8 % across the range, worst near 0–4°, where Michel's criterion puts transition at 0.25–0.08c on the upper surface and the estimate of the turbulent run (fixed H = 1.4) is crudest.

(c) `p = own_polar("0012", alphas, 8.86e6)`; `write_tool_polar("own_naca0012_Re8.86e6.csv", p["meta"], ...)`; the file reads back with tool "own code (P10.05 panel + P10.06 boundary layer)", version "lesson P10.21" and settings panels = 200, re = 8.86e6, transition = Michel, squire_young_x = 0.99. Your code is one more model; next month you will not remember which panel count or transition setting produced a polar, and a reader must be able to tell a computed number from a measured one.

**Common mistakes.** Comparing at the default 5.97 × 10⁶ instead of your Re; calling the own-code polar "measured".

## 5. Which setting matters?

(a)

| panels | cl at 4° | cd at 0° | cd at 4° |
|---|---|---|---|
| 80 | 0.48237 | 0.00596 | 0.00661 |
| 160 | 0.48257 | 0.00622 | 0.00683 |
| 320 | 0.48262 | 0.00633 | 0.00691 |

cl changes in the fourth figure (it converges at second order); cd changes by 6 % from 80 to 320 panels and is not converged.

(b) Squire–Young at x/c = 0.98, 0.99, 0.995: cd = 0.00618, 0.00627, 0.00635 at 0°; 0.00678, 0.00688, 0.00697 at 4°: about ±1.4 %.

(c) Forced at 0.05c: own cd = **0.00822** at 0°, against Table VII **0.00895** at 6.00 × 10⁶. The wraparound No. 60 grit strip is roughness, not just a trip: it adds its own drag and thickens the turbulent boundary layer from the leading edge, so a pure transition model underestimates it (−8 %).

(d) Ranking for cd: transition setting (+31 % between free and forced at 0.05c) ≫ panel count (6 %) > Squire–Young station (±1.4 %). All three belong in the header; for a tool, the analogues are xtr/Ncrit, the panel count, and the convergence settings.

## 6. Tool, code and tunnel disagree

(a) The wind-tunnel data are the reference: they measured the flow; both codes model it. Ladson's coefficients already include the standard two-dimensional wall corrections of Allen and Vincenti, which the report says amount to about 2 % of the measured values; the comparison must not apply them twice.

(b) The boundary layer's displacement thickness effectively decambers and thickens the section, more on the suction side, and reduces circulation; one-way coupling computes the boundary layer from the inviscid solution but never feeds it back. XFOIL's "Viscous Formulation" (its user primer) couples a two-equation integral boundary layer and wake to the inviscid solution simultaneously (strong coupling), so it captures this effect.

(c) Ncrit is the logarithm of the amplification factor of the most unstable disturbance at which transition is assumed (the e^N method); it stands in for the free-stream turbulence and disturbance environment. A low-turbulence tunnel and calm flight have different environments, so the Ncrit that matches one need not match the other; the value used must be recorded.

(d) "XFOIL, AVL and OpenVSP were not available when this work was done, so no tool output is included. The codes were instead validated against published wind-tunnel data (Ladson, NASA TM-4074), with the differences quantified above."

## Check-yourself answers

1. Because a polar without its tool, version and settings cannot be reproduced or checked; the same airfoil gives different polars at different Re, Ncrit, panel counts and convergence settings.
2. Reynolds number, Mach number, Ncrit (or forced transition locations), the panel count and the iteration limit; and the airfoil coordinates used.
3. No XFOIL, AVL or OpenVSP output ships with this lesson; the tools were not available when it was written. The tests use small synthetic fixtures written by the tests themselves.
4. The lift-curve slope and zero-lift angle from a line fit over a stated linear range, the RMS cl difference at common angles in that range, and the minimum drag.
5. About +11 to +16 % in cl_α (0.1205 per degree against 0.104–0.109); the missing effect is the boundary layer's displacement, which one-way coupling leaves out.
6. Within about ±7 % of the free-transition minimum drag between 3.94 and 8.86 × 10⁶ (+18 % at 2 × 10⁶); transition location is the dominant uncertainty.
7. XFOIL and AVL: GNU GPL; OpenVSP: NASA Open Source Agreement 1.3.
