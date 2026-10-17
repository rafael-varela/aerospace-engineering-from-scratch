# P04.14  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. `transition_wavelength(n_u, 2)`, `photon_energy_ev`, `photons_per_pulse`, `lidar_range`; E₁ scales as Z²; dose = ΦS in MeV/g, and 1 MeV/g = 1.602 × 10⁻¹⁰ Gy.
2. Multiply the ODE by e^{λ₂t}: the left side becomes (N₂e^{λ₂t})′. For the peak, differentiate e^{−λ₁t} − e^{−λ₂t}. For the ratio, divide λ₂N₂ by λ₁N₁ and let t → ∞.
3. `q_value_mev(241.056829, 237.048173, ATOMIC_MASS_U["He-4"])` and `decay_heat_w_per_g`. The inventory is 110/(η · fraction remaining · specific power). Equal inventories: a ratio of specific powers against a ratio of decay factors.
4. `fit_attenuation(0.27, 20, 0.81, D2)`; half-value thickness = L ln 2 / 0.27 mm; `areal_density_for_dose(rating/20, D0, L)`; mass = areal density × 600 cm².
5. `seu_rate_per_day`; P(none) = e^{−rt}; per-word rate r_w = σφ × 64; N_words = bits/64. Solve N_words (r_w T)²/2 = 10⁻⁶ for T.
6. Ask for each case: is it cumulative or sudden? Does it depend on the dose, on displaced atoms, or on one particle? Does shielding stop the particles that cause it?
