# P07.05  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. π = 1.2 + s/100. 2. `isentropic_temperature(gas, 288.15, pi)` and `compressor_exit_temperature(gas, 288.15, pi, 0.90)` for `constant_cp_gas(1.4)` and `AIR`. 3. `specific_work`, `entropy_generation`, then multiply s_gen by 288.15 K. 4. `polytropic_efficiency(gas, T1, T2, pi)`; expand (8) of the lesson for small ln π to see why η_p ≈ η_s.

## 2.
(a) ds = cp dT/T − R dp/p; integrate with ds = 0; define s°(T) = ∫cp dT/T; for constant cp take exponentials and use R/cp = (γ − 1)/γ. (b) cp dT = RT dp/(η_p p); divide by T; integrate; for constant cp solve for T₂/T₁ and substitute into (h₂ₛ − h₁)/(h₂ − h₁) = (T₂ₛ − T₁)/(T₂ − T₁). Turbine: same steps with η_p multiplying. (c) T ds = dh − v dp with v dp = η_p dh; integrate. s₂ − s₁ = s°(T₂) − s°(T₁) − R ln π, then use (b). (d) With x = π^m and k = 1/η_p, η_c = (x − 1)/(x^k − 1); convexity gives x^k − 1 > k(x − 1) for x > 1, so η_c < 1/k = η_p.

## 3.
1. π = 75 000/22 632. 2. `compressor_exit_temperature(AIR, 216.65, pi, 0.80)`; `specific_work`; times 1.2 kg/s. 3. `entropy_generation`; times T₀ = 216.65 K. 4. `polytropic_efficiency`. 5. Compare T₂ with 290 K and say where the heat goes.

## 4.
1. n = 8 + (s mod 7), π = 20 + s. 2. `T = multistage_compressor(AIR, 288.15, pi, n, 0.88)`; `isentropic_efficiency` and `polytropic_efficiency` on (T[0], T[-1]). 3. `polytropic_exit_temperature(AIR, 288.15, pi, 0.88)`; the difference. 4. Loop over n until |T_exit − limit| < 1 K. 5. `specific_work(AIR, T[i], T[i+1])` for i = 0 and n − 1.

## 5.
1. `polytropic_exit_temperature(gas, 1600, 0.25, 0.90)`, `isentropic_efficiency(gas, 1600, T2, 0.25)`, `specific_work` for `AIR`, `constant_cp_gas(1.4)` and `constant_cp_gas(1.33)`. 2. γ = cp/(cp − R) from `AIR.cp(1600)`. 3. Reheat: the argument of §2.5 of the lesson.

## 6.
(a) Loop over η_p with `polytropic_exit_temperature`, `specific_work`, `isentropic_efficiency`; compare with the constant-cp row of the lesson. (b) Read the lesson's sweep: where is the work error below 1 %, and how large is the temperature error there? (c) Think about what each effect does to the delivery temperature and to the shaft work separately.
