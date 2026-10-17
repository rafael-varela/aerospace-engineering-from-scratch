# P07.07  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. π = 8 + s/5; x = π^0.2857; τ = 1400/288.15; cp = 1004.7 J/(kg·K). 2. T₂ = T₁[1 + (x − 1)/η_c]; T₇ = T₃[1 − η_t(1 − 1/x)]. 3. w_c, w_t, w_net, q_in = cp(T₃ − T₂), η. 4. Ideal: η_c = η_t = 1. 5. `cycle_performance(gas, CycleSpec(pi=..., T_tit=1400, eta_c=0.85, eta_t=0.88))` for both gases. 6. Set w_net = 0 in (1): solve for x.

## 2.
(a) T₅ = T₇ = T₃/x; q_in = cp(T₃ − T₃/x) = w_t; η = 1 − w_c/w_t; equate with 1 − 1/x. (b) w_c ∝ a^m + (π/a)^m − 2; differentiate in a; check the second derivative. (c) w_c = 2cpT₁(x^½ − 1); q_in = w_t; η = 1 − w_c/w_t. (d) Write s_gen for each component with its reservoir terms; sum around the loop (stream entropies cancel); multiply by T₀; substitute q_ic + q_ex = q_in − w_net.

## 3.
1. `CycleSpec(pi=4, T_tit=1250, **REAL)` and with `**REC` at each ε (override `eps`). 2. `cycle_performance`: η, w_net, q_in. 3. Fuel per unit work: q_in/w_net; ratio recuperated to simple. 4. Think about mass, volume and what an aircraft carries per kilogram.

## 4.
1. T_tit = 1300 + 5s. 2. `pis = np.linspace(2, 40, 77)`; for each entry of `CONFIGS`, `sweep(AIR, CycleSpec(T_tit=..., **REAL, **extra), pis)`. 3. argmax of η and of w; η at π = 12. 4. Weigh specific work (size) against efficiency (fuel) for each application.

## 5.
1. `entropy_generation(AIR, CycleSpec(pi=12, T_tit=1400, **REAL, **IC, **REC))`; multiply by 288.15. 2. `cycle_performance` for q_in and w_net; evaluate the right side. 3. Compare the temperature differences across which heat flows in each component, and which of them a design parameter (ε, intercooler approach, π) controls.

## 6.
(a) `cycle_performance` at π = 25, 30, 40 for the three configurations; `cycle_states` to compare T₄ and T₇. (b) Override the four loss fractions with doubled values. (c) Weight, volume and frontal area, extra pressure losses and the intercooler's cooling-air drag; reliability.
