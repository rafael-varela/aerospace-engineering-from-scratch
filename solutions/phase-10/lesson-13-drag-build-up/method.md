# P10.13  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. V = 40 + s; Re = 1.2V/ν. 2. `cf_mixed(Re, x)` for the five x. 3. Saving = 1 − C_f(0.5)/C_f(0).

## 2.
(a) Drag per unit width ∝ C_f L: turbulent over L, minus turbulent over x_tr (C_f,t(Re_tr) x_tr), plus laminar over x_tr (C_f,l(Re_tr) x_tr); divide by L. (b) d/dC_L of C_L/(C_D0 + kC_L²) = 0 ⇒ C_D0 = kC_L². (c) V* = √(2(W/S)/(ρC_L*)), C_L* ∝ C_D0^{1/2}.

## 3.
1. c = S/b, AR = b²/S. 2. Wing: `component("wing", c, 2·perimeter_ratio("2412")·(S − 0.6c), FF(0.14), 1.0, 30, 0.5)`. 3. Fuselage: profile on a grid, `body_wetted_area`, `form_factor_body(6.5/0.6)`, turbulent. 4. Tails: `component("tails", 0.6, 2·perimeter_ratio("0010")·2.2, FF(0.10), 1.05, 30)`. 5. Sum, × 12/11, × 1.05. 6. e from (4) with `lifting_line_e(AR, 0.4)`; (5); sink = V*/(L/D).

## 4.
1. `cd0_build_up(trainer_components(), gear_dq=g)`, `oswald_e`, `best_ld` for g = 0 and 0.25. 2. P₀ = `cruise_drag(55, CD0, e)["P"]` of the baseline; `brentq(lambda V: cruise_drag(V, ...)["P"] − P₀, 30, 120)`.

## 5.
1. Group Table I rows by Re; nearest to 4.2° and to 0°. 2. `effective_transition(Re, cd)` for both. 3. Compare the measured cd with `section_cd(Re, 0.12, 0.0)`. 4. Relate to the trainer's cruise C_L of 0.33.

## 6.
(a) Use the lesson's Experiments 3 and 4 and your Problem 4; add the interference factors (vary Q by ±0.1) and the misc allowance. (b) Combine independent contributions in quadrature. (c) Fraction of C_D0 that is wing friction; the build-up's x_tr choice. (d) Power-off glides at several speeds, or level flight with measured power and propeller efficiency; plot C_D against C_L².
