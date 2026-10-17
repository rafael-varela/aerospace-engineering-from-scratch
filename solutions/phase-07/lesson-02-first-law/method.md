# P07.02  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
(a) v₁ = R T₁/p₁; v₂ = 3v₁; p₂ = p₁(1/3)ⁿ; T₂ = p₂v₂/R; `polytropic_work(p1, v1, p2, v2, n)`; Δu = cv(T₂ − T₁); q = Δu + w. Compare across paths: what depends on the end states only? (b) m = p₀V/(R T₀); p₂ = p₀T₂/T₀ (rigid, closed); Q = m cv ΔT; substitute m and cv = R/(γ − 1).

## 2.
(a) E_cv = m u; d(mu)/dt = −ṁh + Q̇; m du/dt = −ṁ(h − u) + Q̇; du = cv dT; h − u = R T. (b) T₀/T = 1 + (γ − 1)M²/2; ρ/ρ₀ = (T/T₀)^(1/(γ−1)); ṁ/A = ρ M √(γRT); set M = 1; critical ratio p*/p₀ from T*/T₀ = 2/(γ + 1). (c) Entropy of the remaining gas is constant (no heat, no dissipation inside), so p ∝ ρ^γ; separate dr/dt = −r^((γ+1)/2)/τ; integrate from r = 1; τ = Vρ₀/ṁ₀. (d) Show that d(U + H_out − Q_in)/dt = 0 and that an RK step is a linear combination of right-hand-side values.

## 3.
(a) Energy balance for the filling CV with m₁ = 0: m₂u₂ = m₂h_line; cv T₂ = cp T_line. (b) m = p V/(R T₂); after cooling at constant V, p ∝ T; isothermal mass p V/(R T_line). (c) Write cv T₂ m₂ − cv T₁ m₁ = cp T_line (m₂ − m₁) with m₂ = p₂V/(R T₂); solve for T₂ (linear after multiplying out, or with `brentq`). (d) Heat transfer, overshoot of temperature limits, top-up after cooling.

## 4.
(a) A = πd²/4; `mass_flow(30e6, 300, 0.5e6, 0.8, A, R, 1.4)`; `blowdown_time_constant`; `adiabatic_time_to_pressure(3e6, ...)`; τ ln 10. (b) Two `simulate_blowdown` calls with `p_stop=3e6`; read `t_stop`, `T[-1]`, `m[-1]`; `energy_residual`. (c) `recovery_pressure(3e6, T_end, 300)`. (d) (2/d)².

## 5.
(a) Six `simulate_blowdown(..., mode="wall", hA=...)` calls (hA = 0 is adiabatic). (b) m₀cv/hA; the crossover hA = m₀cv/τ. (c) Compare the two crossover hA values and the τ values. (d) Think about tank mass per unit of pressurant mass, the boiling point, and the ullage of a liquid-oxygen or liquid-hydrogen tank.

## 6.
For each case, identify whether the issue is the equation of state, stratification, the wall model, a phase change, or the direction of the error; name the better model or the test.
