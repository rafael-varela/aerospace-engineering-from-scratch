# P07.10  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. a = 12 + 23/4; inerts per O₂ from AIR_X; f_st = M_fuel/(a/0.2095 × M_air). 2. φ = 0.3 + s/50; `reactants`, `complete_products`; normalise to mole fractions. 3. p_v = x_H2O × 101 325; `humidity_ratio` and `dew_point` from P07.09 at 101 325 Pa.

## 2.
(a) Atom balances for C, H, O. (b) Subtract the two Hess's-law expressions; Hf(g) − Hf(l) per mole = M_w h_fg. (c) H_r(T₀) = H_p(T₀) + m_fuel LHV (definition of LHV, lean products contain the spare O₂ and inerts on both sides); H_p(T_ad) − H_p(T₀) = m_p c̄_p (T_ad − T₀). (d) h = h(T) for ideal gases; complete products do not depend on p.

## 3.
1. `phi_for_temperature(fuel, 1700, "air", 298.15, 800)` for H₂ and Jet-A; `fuel_oxidizer_ratio`. 2. f × `heating_values(fuel)[0]`. 3. n_H2O × M_w × f/M_fuel. 4. `adiabatic_flame_temperature(fuel, 1, "air", 298.15, 800)`.

## 4.
1. T3, T4 from s; solve φ; f; 50 f. 2. T_ad at φ = 1. 3. Central difference; ΔT ≈ (dT/dφ)(0.01 φ). 4. `mixture_enthalpy(products, T4) − reactant_enthalpy(...)`.

## 5.
1. `heating_values`; (y/2) M_w/M_fuel × (h_g − h_f at 25 °C). 2. ΔH_r(1000 K) per kg. 3. Σ cp of products against reactants over 298–1000 K.

## 6.
(a) Compare T with 2000 K and with the `P07.11` gaps. (b) `complete_products("CH4", 3, "O2")`, `adiabatic_flame_temperature`; think soot, CH₄ survival, equilibrium water-gas split. (c) The liquid enters with less enthalpy.
