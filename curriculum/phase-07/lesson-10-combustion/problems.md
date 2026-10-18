# P07.10  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer, and say which numbers are complete-combustion values. **Everything here is computed: there is no combustion, flame, burning, fuel handling or fuel of any kind.**

Several problems use your personal seed **s**: the day of the month you were born (1–31).

Data: the module's NASA polynomials and dry-air model; the steam-table h_f and h_g at 25 °C; your `P07.09` module for dew points. Combustor temperatures, pressures and flows are typical, illustrative values.

## 1. (Warm-up: Jet-A stoichiometry and the exhaust's dew point)

(a) By hand, balance the complete combustion of C₁₂H₂₃ in the phase's dry air at φ = 1 and find f_st and the air–fuel ratio.
(b) At φ = 0.3 + s/50, find f, the reactant and product amounts per mole of fuel, and the products' mole fractions.
(c) At 1 atm, what is the partial pressure of water in those products and their dew point (use your `P07.09` code)? What does this say about contrails?

Answer form: a balanced equation and two numbers; a number, two tables; a pressure, a temperature and a sentence.

## 2. (Derivation)

(a) Derive a = x + y/4 − z/2 for C_xH_yO_z.
(b) Show that HHV − LHV = (y/2)(M_w/M_fuel) h_fg(298.15 K) and evaluate it per kilogram of methane by hand.
(c) Split the adiabatic-combustor first law into two steps and show that T_ad − T₀ = m_fuel LHV/(m_products c̄_p) for reactants entering at T₀ = 298.15 K, with c̄_p the products' mean cp over [T₀, T_ad].
(d) Show that the complete-combustion T_ad of ideal gases does not depend on pressure. Which assumption makes it so?

Answer form: four derivations.

## 3. (Unfamiliar scenario: a hydrogen-fuelled turbofan)

A turbofan is converted to burn gaseous hydrogen. Compressor delivery is 800 K and the turbine inlet must stay at 1700 K (typical, illustrative).

(a) Find φ and f for hydrogen and for Jet-A, and the fuel mass per 100 kg of air for each.
(b) Find f × LHV per kilogram of air for each. Why are they nearly, but not exactly, equal?
(c) Find the mass of water in the products per kilogram of air for each. Why does this matter for contrails?
(d) Compare the stoichiometric T_ad of each at the same inlet. Which fuel puts the hotter primary zone into the combustor?

Answer form: four numbers per fuel and a sentence each for (b)–(d).

## 4. (Numerical, uses your code: your own combustor)

Compressor delivery T3 = 700 + 10 s K; turbine inlet T4 = 1600 + 10 s K; Jet-A vapour at 298.15 K.

(a) Find φ, f and the fuel flow for a 50 kg/s core air flow.
(b) Find the stoichiometric T_ad at the same inlet.
(c) Estimate dT4/dφ by a central difference and the change in T4 from 1 % more fuel.
(d) Check your (a) with the energy balance residual of (5).

Answer form: three numbers; one number; two numbers; a residual.

## 5. (Numerical, uses your code: the heat of reaction is not constant)

(a) For CH₄, H₂ and Jet-A compute LHV and HHV and check HHV − LHV against the steam table.
(b) Compute the heat released per kilogram of fuel when stoichiometric reactants and complete products are both at 1000 K (use `mixture_enthalpy` and `h_species`), and compare with the 298.15 K LHV.
(c) Explain the sign of each change from the heat capacities of reactants and products (Kirchhoff's law: d(ΔH_r)/dT = Δc_p).

Answer form: a table; three numbers and three differences; a paragraph.

## 6. (Judgment: where complete combustion stops)

(a) Which of these complete-combustion results would you use in a design review, and which not: the turbine-inlet φ of Problem 4, the stoichiometric Jet-A T_ad at 800 K, the CH₄/O₂ stoichiometric T_ad of 5166 K?
(b) A methane-fuelled rocket's fuel-rich preburner runs at φ = 3 with oxygen. What products and temperature does this lesson's model give? Name two things that are physically wrong with them.
(c) The Jet-A in this lesson is a vapour at 298.15 K. In which direction does a liquid fuel move the turbine-inlet φ, and why?

Answer form: three verdicts with reasons; a composition, a temperature and two reasons; a direction and a reason.
