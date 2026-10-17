# P04.14  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
(a) n_u = 4 + (s mod 4); λ from `transition_wavelength(n_u, 2)`; E from `photon_energy_ev(λ)`. (b) `photons_per_pulse((5 + s)e-3, 532e-9)`; `lidar_range(s × 1e-6)`. (c) −Z² × 13.6057 eV. (d) `dose_from_fluence_gy(1e10, 34.6)`; divide by 10 for krad.

## 2.
(a) (N₂e^{λ₂t})′ = λ₁N₁(0)e^{(λ₂−λ₁)t}; integrate from 0 to t; multiply back by e^{−λ₂t}. (b) Set d/dt(e^{−λ₁t} − e^{−λ₂t}) = 0; solve for t. (c) λ₂N₂/(λ₁N₁) = λ₂/(λ₂ − λ₁) × (1 − e^{−(λ₂−λ₁)t}); limit; then Sr-90/Y-90 with consistent units (hours). (d) Atoms per gram N_A/M; decays per second λN_A/M; times Q in joules.

## 3.
(a) Q from masses; `decay_heat_w_per_g(432.6 y in s, 241.0568, Q)`. (b) `rtg_fuel_mass_g` for Pu-238; the same formula with the Am-241 half-life and specific power. (c) Write m_Pu(t) = m_Am(t) as w_Am/w_Pu = 2^{−t/T_Pu}/2^{−t/T_Am}, take logarithms, solve for t (or use `brentq`); compare with realistic mission lengths.

## 4.
(a) `fit_attenuation`; L ln 2 in g/cm², divide by 0.27. (b) Allowed rate = rating/(10 × 2); `areal_density_for_dose`; × 600 cm² → grams. (c) Refit with the perturbed point; compare relative changes; think about the logarithm. (d) What the single exponential leaves out at larger depth.

## 5.
(a) r = `seu_rate_per_day(σ, 1, 4·2³⁰)`; P = exp(−r × 95/1440); `default_rng(s).poisson(mean, 100000) == 0`. (b) r_w = σ × 64 per second; N_w = 4·2³⁰/64; evaluate. (c) T = √(2 × 10⁻⁶/N_w)/r_w; check r_w T ≪ 1.

## 6.
For each case decide cumulative (TID, displacement) vs single-particle (SEU, SEL) vs environmental charging; then match: TID → shielding or harder part with margin; displacement → design margin (end-of-life sizing); SEU → EDAC and scrubbing; SEL → current limiting and power cycling; charging → grounding and shielding per the handbook. For (e) ask whether the model is valid at 12 mm.
