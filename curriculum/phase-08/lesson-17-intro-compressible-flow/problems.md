# P08.17  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact for a perfect gas (γ = 1.4, R = 287.05287 J/(kg·K)), which assume isentropic flow, and which assume a normal shock.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

Airspeed in knots appears only where aviation practice uses it (1 kt = 0.514444 m/s). Flight level 350 is a pressure altitude of 35 000 ft = 10 668 m; treat it here as an ISA geopotential altitude of 10 668 m.

## 1. (Warm-up: an airliner in cruise)

An airliner cruises at FL350 at a true airspeed V = 230 + s m/s.

(a) Find T, p, ρ and a from the ISA, and the Mach number.
(b) Find the total temperature, the total pressure and the impact pressure q_c = p₀ − p, and compare q_c with ½ρV².
(c) What airspeed would the incompressible Bernoulli formula of `P08.04` give from q_c and the true static density, and what is its error? By how much does the density rise when the air is brought to rest isentropically?

Answer form: four properties and a Mach number; three pressures, a temperature and a comparison; a speed, a percentage and a ratio.

## 2. (Derivation: total conditions and the shock)

(a) From the steady adiabatic energy equation h + V²/2 = h₀ for a calorically perfect gas, derive T₀/T = 1 + (γ − 1)M²/2. Then use the isentropic relation p ∝ ρ^γ ∝ T^{γ/(γ−1)} to derive p₀/p and ρ₀/ρ. Which of the three relations needs the flow to be isentropic, and which only adiabatic?
(b) Expand p₀/p for small M with the binomial series to show that p₀ − p = ½ρV²(1 + M²/4 + (2 − γ)M⁴/24 + …). Deduce that the incompressible airspeed is too high by about M²/8.
(c) For a normal shock, write the mass, momentum and energy balances across a thin control volume (`P08.03`) and derive p₂/p₁ = 1 + 2γ(M₁² − 1)/(γ + 1). Then derive Rayleigh's pitot formula by combining the shock with an isentropic compression from M₂ to rest.

Answer form: three derivations with the assumptions named.

## 3. (Unfamiliar scenario: a high-speed wind tunnel)

A blowdown wind tunnel's settling chamber holds air at p₀ = 200 kPa, T₀ = 290 K, and the test section runs at M = 0.8 (the flow from the chamber to the test section is isentropic).

(a) Find the test-section static pressure, temperature, density and speed, and the mass flux ρV per unit area.
(b) A technician estimates the speed from Bernoulli with the measured p₀ − p, once with the static density and once with the chamber density. What does each give, and how large are the errors?
(c) Which density would you have to use for the incompressible formula to be right, and why is there no single correct one?

Answer form: five numbers; two speeds and two percentages; a sentence.

## 4. (Numerical, uses your code: pitot readings on a supersonic rocket)

A sounding rocket's nose pitot probe and a static port read during the climb through the altitude h = 5 + s/3 km.

(a) For M = 1.5, 2.0, 2.5 and 3.0 compute the ratio of pitot to static pressure with `rayleigh_pitot`, the pitot reading in kPa, and the Mach number recovered by `mach_from_pitot`.
(b) For each, what Mach number would the subsonic (isentropic) inversion give from the same reading, and what speed would incompressible Bernoulli give? Tabulate the errors.
(c) Give the total-pressure ratio across the shock for each Mach number, and explain why the subsonic inversion under-reads.
(d) Compute the stagnation temperature at the probe tip at M = 3. What does it imply for the probe's materials and for the pressure transducer's location?

Answer form: a table of four rows; two columns added to it; a column and a sentence; a temperature and a sentence.

## 5. (Numerical, extends your code: calibrated and equivalent airspeed)

An airspeed indicator is calibrated to show the airspeed a pitot-static system would read at sea level: calibrated airspeed CAS = a_SL √(5[(q_c/p_SL + 1)^{2/7} − 1]). Equivalent airspeed EAS = V√(ρ/ρ_SL).

(a) Write a function `cas(qc)` and compute CAS and EAS (in m/s and knots) for the light aeroplane at 3000 m and 60 m/s, and for the airliner of Problem 1 at FL350.
(b) Compute the difference CAS − EAS for each (the "compressibility correction" pilots apply). Why is it negligible for one and not for the other?
(c) Compare CAS with √(2q_c/ρ_SL), the incompressible formula with the sea-level density. Which of CAS, EAS and TAS does a stall speed or a structural limit naturally belong to, and why?

Answer form: four speeds in two units; two differences and a sentence; a comparison and two sentences.

## 6. (Judgment: is the flow compressible?)

(a) The rule of thumb says flow is incompressible below M ≈ 0.3. Compute the stagnation density change ρ₀/ρ − 1 and the Bernoulli airspeed error at M = 0.1, 0.3 and 0.5. What does "incompressible" mean for each quantity?
(b) For each case, say whether you would use the incompressible formula and why: (i) the light aeroplane's airspeed indicator (accuracy wanted about 1 %); (ii) a cabin-air duct at 10 m/s; (iii) the airliner's air-data computer; (iv) the wind-tunnel speed of Problem 3.
(c) A pitot probe's isentropic stagnation assumes the flow is brought to rest without losses. Name two situations in which a real probe at subsonic speed does not satisfy that, and what each does to the reading.

Answer form: six numbers and two sentences; four verdicts with reasons; two situations.
