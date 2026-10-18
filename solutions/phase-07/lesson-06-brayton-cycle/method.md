# P07.06  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. π = 6 + s/5; T₂ = T₁π^m, T₄ = T₃/π^m; p₂ = p₃ = πp₁; v = RT/p. 2. Works and heats with cp = 1004.7 J/(kg·K); η = w_net/q_in; bwr = w_c/w_t. 3. Carnot 1 − 288.15/1200; ratio. 4. `brayton_performance(AIR, 288.15, 101325, pi, 1200)`; compare T₂ and T₄ through `brayton_states`.

## 2.
(a) q_out/q_in = (T₄ − T₁)/(T₃ − T₂) and T₄/T₁ = T₃/T₂. Δs of 2→3 equals Δs of 1→4; divide q by Δs. T̄_in lies between T₂ and T₃; T̄_out between T₁ and T₄. (b) Substitute T₂ = xT₁, T₄ = T₃/x into w_t − w_c; differentiate with respect to x. (c) Write w_net = h₃ − h₄ − h₂ + h₁; differentiate h₂ and h₄ with respect to ln π using dh = v dp on each isentrope (for h₄, the inlet pressure of the turbine moves). (d) ∮d(pv) = 0; product rule; shaft work of each device is −∫v dp for reversible adiabatic flow and zero for the isobars.

## 3.
1. R = 8.314/0.040; `constant_cp_gas(5/3, R)`. 2. π_opt from the closed form; `brayton_performance` there. 3. Same with γ = 1.4 for air. 4. η at x = √τ. 5. For T₁ = 400 and 500 K: π_opt, q_out, w_net, (q_out/w_net)/T₁⁴; ratio of the two.

## 4.
1. T₃ = 1300 + 10s. 2. `optimum_pressure_ratio(gas, 288.15, 101325, T3)` and `brayton_performance` there; closed form. 3. `brayton_states` at the NASA optimum: T₂ − T₄. 4. Loop n over 20–160: `cycle_path`, `enclosed_area(path["v"], path["p"])` and `enclosed_area(path["s"], path["T"])`; `observed_order`.

## 5.
(a) Compute T̄_in and T̄_out (q/Δs) and compare with T₃ and T₁; note that the ideal cycle has no friction. (b) Component losses; pressure losses and cooling air; real combustion and the fuel mass; for each, the direction and the lesson. (c) Work per kilogram near π_max, and delivery temperature and stage count.
