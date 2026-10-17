# P07.09  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer, say which inputs are illustrative, and label any dew point below 0.01 °C as an estimate. Draw the control volume and write the dry-air, water and energy balances before any other equation.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

Data: the module's dry-air model and constants; saturation pressures from `p_sat` (the steam table of `P07.03`); ISA pressures from `isa_state`. Moisture rates and flows below are typical, illustrative values, not data.

## 1. (Warm-up: the cabin's air at altitude)

A cabin is pressurised to the ISA pressure of a cabin altitude of (6000 + 100 s) ft (1 ft = 0.3048 m).

(a) Find the cabin pressure and the partial pressure of oxygen, and compare with sea level. What fraction of the sea-level oxygen partial pressure does a passenger breathe?
(b) Compute the mass fractions of dry air's four components and the mole-specific and mass-specific entropy of mixing them from pure gases.
(c) Cabin air at 22 °C and 50 % RH: find ω, the dew point and h at the cabin pressure and at sea level. Which of the three change with pressure, and why?

Answer form: two pressures and a ratio; four fractions and two entropies; a table of six numbers and a sentence.

## 2. (Derivation: the three results the lesson rests on)

(a) Starting from Dalton's law and the ideal-gas law for each component, derive ω = ε p_v/(p − p_v) and its inverse, and evaluate ε for the phase's dry air.
(b) Prove that the outlet state of adiabatic, work-free mixing of two moist streams lies on the straight line between the inlet states in the (ω, h) plane and divides it in the ratio ṁ₂/ṁ₁. Then explain, using h = h_a(T) + ω h_v(T), why the outlet *temperature* is not exactly the dry-air-mass-weighted mean of the inlet temperatures.
(c) For the cabin and mixing manifold of §2.6, derive ω_cabin = ω_out + ṁ_w/ṁ_fresh and show that the pack outlet temperature does not depend on the recirculated flow.

Answer form: three derivations.

## 3. (Unfamiliar scenario: a spacecraft cabin with no fresh air)

A crewed spacecraft cabin is held at 101.325 kPa and 22 °C. It has no outside air: four crew each release 1.8 kg/day of water vapour (typical, illustrative), and a condensing heat exchanger (CHX) cools part of the cabin air to (4 + s/10) °C, condenses water, and returns the air saturated at that temperature.

(a) Write the water balance of the cabin with the CHX and find the dry-air flow through the CHX that holds the cabin dew point at 10 °C.
(b) What is the cabin relative humidity?
(c) Using `cool_dehumidify`, find the CHX heat load and how much of it is sensible (cooling the dry air) rather than latent and vapour cooling.
(d) What happens to the cabin humidity if the CHX flow halves? Answer from your balance without running the code again.

Answer form: a mass flow; a percentage; two powers; a sentence with a number.

## 4. (Numerical, uses your code: your own airliner)

An aircraft carries n = 150 + 2s occupants, supplies (0.2 + s/100) kg/min of fresh air per occupant, recirculates 40 % of the supply, and otherwise uses the defaults of `ecs_cruise`.

(a) Run `ecs_cruise` and report the fresh, recirculated and supply flows, ω_cabin, the cabin RH and dew point, the supply and pack outlet temperatures and the pack duty.
(b) Reproduce ω_cabin and the RH by hand.
(c) Find, with `brentq`, the fresh air per occupant that would give 20 % RH. What pack outlet temperature does that need, and is it a practical design point?

Answer form: nine numbers; a hand calculation; a flow, a temperature and a sentence.

## 5. (Numerical, uses your code: how ideal is water vapour?)

(a) Compare `p_sat` and `p_sat_gibbs` at 5, 20, 40 and 60 °C.
(b) Using v_g from Yan Appendix A Table A1 (147.011, 57.7567, 19.5151 and 7.66723 m³/kg at those temperatures), compute the compressibility factor Z_g = p_sat v_g/(R_w T) of saturated vapour at each, with R_w = R/M_w.
(c) Estimate h_fg from the Clausius–Clapeyron slope between 20 and 40 °C, once with each route, and compare with the table's h_fg at 30 °C (2555.55 − 125.73 kJ/kg).
(d) Explain the pattern in (a) using (b). Which deviation at 5 °C is not real-gas behaviour?

Answer form: two tables; two estimates and a reference value; a paragraph.

## 6. (Judgment: should the airline humidify the cabin?)

A passenger-comfort study proposes raising the cruise cabin of the lesson's aircraft to 30 % RH at 22 °C.

(a) Using the lesson's ECS (186 occupants, 0.25 kg/min fresh air each), how much water per hour must be added, beyond what the occupants supply, to hold 20 % and 30 %? How much is that over a 10-hour flight?
(b) The inner face of the fuselage skin behind the insulation runs at about −30 °C in cruise (typical, illustrative). Compare the cabin dew point at 18.4 %, 20 % and 30 % with it. What happens to the water that reaches the skin, and where does it go on descent?
(c) Give two further reasons for and two against humidifying, and say which number in this lesson you would most want measured before deciding.

Answer form: two flows and a mass; three temperatures and two sentences; a short list and a sentence.
