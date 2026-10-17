# P04.07  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. T = 200 + 5s. 2. `most_probable_speed`, `mean_speed`, `rms_speed` with M_N2. 3. `number_density(101325, T)`; `pressure_from_momentum_transfer(n, M/N_A, v_rms)`. 4. (3/2)kT; divide by the electron charge. 5. `fraction_above(2 v_rms, T, M)`. 6. `speed_of_sound(T)` for air; `sound_to_rms_ratio(1.4)`.

## 2.
(a) Force of one molecule on the x-wall mv_x²/L; sum; ⟨v_x²⟩ = ⟨v²⟩/3 by isotropy; divide by L². (b) Molecules with speed in dv and direction in dΩ cross unit area in unit time: n f(v) dv (dΩ/4π) v cos θ for θ < π/2; ∫cos θ dΩ over the hemisphere = π; so the flux is n⟨v⟩π/(4π) = n⟨v⟩/4. (c) d/dv(v² e^{−av²}) = 0 → v_p = 1/√a; ⟨v⟩ = ∫ v f dv, ⟨v²⟩ = ∫ v² f dv with the given integral; a = m/2kT. (d) p/ρ = kT/m = v_rms²/3; a² = γ v_rms²/3; γ ≤ 5/3 < 3.

## 3.
(a) ρ = pM/(RT); `speed_of_sound(210, 1.29, 0.04401)`; `rms_speed(210, 0.04401)`. (b) f = 2/(γ − 1); interpret against 5 (linear molecule, no vibration) and 7 (one vibrational mode). (c) `mean_free_path(210, 610, 3.9e-10)`; divide by the lengths. (d) `fraction_above(5027, 200, 1.008e-3)` and with M_H2.

## 4.
(a) `isa_temperature`, `isa_pressure`, `isa_density`, `isa_speed_of_sound` at h; σ = ρ/1.225. (b) V = 0.85 a; q = ½ρV². (c) `hydrostatic_residual(np.array([h]))`. (d) T and p at 80 km; `mean_free_path`; Kn = λ/4; scan h to 86 km for Kn > 0.01 for ℓ = 4 m.

## 5.
(a) rng = default_rng(s); random unit directions times v_rms(300); for each target per-molecule count c, call `relax_to_maxwellian(v, int(c N) − done, rng)`; record the statistics. (b) Read off the first stage within 1 %. (c) Another seeded generator; 50 × `sample_speeds(20000, 300, M_N2, rng2)`; mean and std of the kurtosis; compare. (d) The collision count from (b) divided by 7.1 × 10⁹ s⁻¹ (each pair collision involves two molecules: decide how to convert).

## 6.
Classify each by Kn, temperature against the vibrational and dissociation scales, pressure against the ideal-gas range, and the difference between the real day and the ISA.
