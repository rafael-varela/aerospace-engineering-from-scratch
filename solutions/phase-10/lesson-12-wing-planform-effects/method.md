# P10.12  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. AR, S = 100/AR. 2. `efficiency_vs_taper(AR, [1, 0.5, 0.35])`. 3. `stall_onset(10, S, tapered_chord(10, S, λ), zero, 1.0)`. 4. Relate cl = 2Γ/(Vc) to the chord at the tips.

## 2.
(a) Substitute A = Bx into πAR ΣnAₙ²; CL = πAR(Bx)₁. (b) L = Σn(αaₙ + trₙ)² − λ(αa₁ + tr₁ − A₁); two stationarity equations and the constraint; eliminate λ. (c) Half-lift L/2 at 4(b/2)/(3π) from the root; divide by qSb. (d) Fourier expressions for L and I; substitute G₁ = 1/s; positivity near θ = 0.

## 3.
1. `tapered_chord(15, 11.25, 0.4)`; `ll_solve`; `wing_from_coefficients`; `stall_onset(..., 1.2)`. 2. `optimize_twist(..., 0.8, eta5)` and with `tip_margin=0.05, CL_high=1.2, washout_only=True`. 3. `linear_twist_optimum`; tabulate e and stall stations; judge.

## 4.
1. `bending_trade(10, 12, rect, 0.3309, eta5, [ratio])`. 2. ΔD_i = qS CL²/(πAR)(1/e − 1/e_opt) and the cap-area change ∝ C_B. 3. Compare achieved C_B with the limit.

## 5.
1. Evaluate (4). 2. Tip load sign. 3. `prandtl_bell(n_terms)`; the structural model's dependence on higher harmonics (only G₁ and G₃ enter I).

## 6.
(a) Induced drag gain vs washout penalty, both in newtons, and the stall station. (b) Washout lowers the tip cl at every CL. (c) Tip section Reynolds number, manufacturing twist, ailerons, sideslip; tuft tests.
