# P07.06  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact, which are closed forms and which are numerical. Sketch the cycle on T–s axes before writing any equation.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

Gas data: constant cp with γ = 1.4 and R = 287.05 J/(kg·K), or the NASA dry-air model (`AIR`), as each problem says. Engine numbers are typical, illustrative values.

## 1. (Warm-up: an ideal cycle by hand)

An ideal Brayton cycle takes air at 288.15 K and 101 325 Pa through a pressure ratio π = 6 + s/5 to a turbine inlet temperature of 1200 K.

(a) With constant cp, compute T, p and v at the four states, and w_c, w_t, w_net, q_in, q_out, η and the back-work ratio. Check η against 1 − π^(−(γ−1)/γ).
(b) Compute the Carnot efficiency between 288.15 K and 1200 K and the ratio η/η_Carnot.
(c) Repeat (a) with `brayton_performance(AIR, ...)` and explain the changes in η and w_net.

Answer form: a state table and seven numbers; two numbers; two numbers and a sentence.

## 2. (Derivation: efficiency, work, optimum and the areas)

(a) For constant cp, derive η = 1 − π^(−(γ−1)/γ) from q_in and q_out. Then show that η = 1 − T̄_out/T̄_in for any ideal gas, where T̄ = q/Δs is the entropic mean temperature of each heat transfer, and use it to prove η < 1 − T₁/T₃.
(b) Derive w_net = cpT₁[τ(1 − 1/x) − (x − 1)], its maximum at x = √τ, and T₂ = T₄ = √(T₁T₃) there.
(c) For any ideal gas (variable cp), show that dw_net/d ln π = R(T₄ − T₂) at fixed T₁ and T₃.
(d) Show that ∮p dv = −∮v dp around any closed path, and explain why both equal the net shaft work of the steady-flow cycle.

Answer form: four derivations.

## 3. (Unfamiliar scenario: closed-cycle Brayton power in space)

A space nuclear power system runs a closed Brayton cycle on a helium–xenon mixture (molar mass 40 g/mol, monatomic, γ = 5/3; typical, illustrative). The reactor heats the gas to T₃ = 1150 K and a radiator cools it to T₁ = 400 K. Take the components as ideal.

(a) Find the pressure ratio for maximum specific work, and at it w_net, η and the back-work ratio.
(b) Show that the efficiency at the maximum-work pressure ratio is 1 − √(T₁/T₃) for any γ. Compare the He–Xe cycle with an air cycle between the same temperatures: which differs, the efficiency or the pressure ratio and the work?
(c) In space all heat is rejected by radiation, and a radiator's area scales as the heat it rejects divided by T⁴. Using T₁ as the radiating temperature (a conservative simplification), compare the radiator area per watt of net power for T₁ = 400 K and T₁ = 500 K (each at its own optimum π). Which would you choose for a spacecraft, and why does the answer differ from a power station's?

Answer form: four numbers; a proof and a sentence; two numbers, a ratio and a short argument.

## 4. (Numerical, uses your code: optimum and convergence)

Take T₃ = 1300 + 10s K and π = 12.

(a) Find the pressure ratio of maximum work numerically with `optimum_pressure_ratio`, for constant cp and for NASA air, with w_net and η at each. Compare with (T₃/T₁)^(γ/(2(γ−1))) and check T₂ = T₄ at the NASA optimum.
(b) For NASA air at π = 12, compute |∮p dv − w_net| and |∮T ds − w_net| with `cycle_path` at n = 20, 40, 80 and 160 points per process, and the observed orders.
(c) Why is the T–s integral more accurate than the p–v integral at the same n?

Answer form: two triples of numbers, a comparison and a check; a table and two orders; a sentence.

## 5. (Judgment: what the ideal cycle can and cannot say)

A colleague writes: "The ideal Brayton efficiency at π = 12 and 1400 K is only 64 % of Carnot, so the gas turbine wastes a third of what is possible through friction and should be replaced."

(a) What is wrong with the reasoning? Use the entropic mean temperatures of the lesson's example.
(b) List three things the air-standard ideal cycle leaves out of a real engine. For each, say whether including it raises or lowers the efficiency, and which later lesson quantifies it.
(c) The ideal cycle's efficiency rises with π without limit up to π_max. Give two reasons a real design does not use π close to π_max.

Answer form: a paragraph with two numbers; a three-row table; two reasons.
