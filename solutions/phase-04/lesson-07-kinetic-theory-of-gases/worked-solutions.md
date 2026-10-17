# P04.07  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for $s = 15$; numbers come from `solution.py` and probe calculations that call it. Output of the two simulations is labelled SIMULATION.

## 1. Nitrogen at 275 K ($s = 15$)

**Assumptions.** Ideal gas, equilibrium (Maxwell–Boltzmann), M(N₂) = 28.0134 g/mol.

(a) v_p = √(2RT/M) = **404.0 m/s**; ⟨v⟩ = √(8RT/πM) = **455.9 m/s**; v_rms = √(3RT/M) = **494.8 m/s**.

(b) n = p/kT = **2.669 × 10²⁵ m⁻³**; (1/3)nm v_rms² = 101 325.0 Pa, equal to p (it is the same identity written twice: (1/3)m v_rms² = kT).

(c) (3/2)kT = **5.695 × 10⁻²¹ J = 0.0355 eV**: independent of the gas.

(d) **7.38 × 10⁻³** of the molecules (0.74 %) exceed 2v_rms = 989.7 m/s.

(e) Air: a = √(1.4 × 287.05 × 275) = **332.4 m/s**. For N₂, a/v_rms = √(1.4/3) = **0.6831**.

**Check.** Ratios v_p : ⟨v⟩ : v_rms = 1 : 1.128 : 1.225 (404.0 × 1.128 = 455.9; × 1.225 = 494.8). The fraction above 2v_rms does not depend on T or M (it is a property of the shape of the distribution): the same 0.74 % as for air at any temperature.

**Common mistakes.** Using the molar mass in g/mol inside √(RT/M); using the mass of a mole instead of a molecule in (b); computing the speed of sound with γ = 1.

## 2. Derivations

(a) A molecule with x-velocity v_x rebounds elastically from the wall at x = L (momentum transfer 2m|v_x|) and returns after 2L/|v_x|, so its average force is mv_x²/L. Summing over N molecules: F = Nm⟨v_x²⟩/L. **Isotropy** gives ⟨v_x²⟩ = ⟨v²⟩/3. Pressure is F/L² = (N/L³)m⟨v²⟩/3 = **(1/3)nm⟨v²⟩**. Collisions between molecules exchange velocities but do not change ⟨v_x²⟩ in equilibrium, so they do not change the result.

(b) Molecules with speed in [v, v + dv] moving into a solid angle dΩ at angle θ to the wall's normal have density n f(v) dv dΩ/(4π) and cross unit area at rate (that density) × v cos θ. Integrate over the hemisphere facing the wall: ∫ cos θ dΩ = ∫₀^{2π}∫₀^{π/2} cos θ sin θ dθ dφ = π. So the rate is n (π/4π) ∫ v f dv = **n⟨v⟩/4**.

(c) Write a = m/2kT. f ∝ v² e^{−av²}; d/dv = (2v − 2av³)e^{−av²} = 0 → **v_p = 1/√a = √(2kT/m)**. With x = √a v and the normalisation 4π(a/π)^{3/2}: ⟨v⟩ = 4π(a/π)^{3/2} ∫ v³ e^{−av²} dv = 4π(a/π)^{3/2} × 1/(2a²) = **√(8kT/πm)**; ⟨v²⟩ = 4π(a/π)^{3/2} × 3√π/(8a^{5/2}) = 3/(2a) = **3kT/m** (the even moment uses ∫ x⁴ e^{−x²} dx = 3√π/8).

(d) p/ρ = nkT/(nm) = kT/m = v_rms²/3. So a² = γp/ρ = γ v_rms²/3 and **a/v_rms = √(γ/3)**. Equipartition gives γ = (f + 2)/f ≤ 5/3 (f ≥ 3), so a/v_rms ≤ √(5/9) = 0.745 < 1. A pressure pulse is carried by the molecules themselves and cannot outrun them.

**Common mistakes.** In (a), using momentum transfer m v_x instead of 2m v_x and then "fixing" the factor 2 elsewhere; in (b), forgetting the cos θ (that gives n⟨v⟩/2); in (c), using the one-dimensional normal moments instead of the speed density.

## 3. Mars

**Assumptions.** Ideal gas, pure CO₂, γ as measured, hard-sphere collisions.

(a) ρ = pM/(RT) = **0.01538 kg/m³** (1.3 % of sea-level air). a = √(1.29 × 8.314 × 210/0.04401) = **226.2 m/s**, 66 % of the Earth's 340.3 m/s; v_rms = **345.0 m/s**. A rotor tip at Mach 0.7 on Mars moves at only 158 m/s.

(b) f = 2/(γ − 1) = **6.9**. A linear molecule has 3 translational and 2 rotational degrees of freedom (f = 5, γ = 1.4); the extra ~2 come from CO₂'s low-frequency bending vibration, which is already partly active at 210 K. CO₂ is "softer" than N₂ in this sense, and its γ is lower.

(c) λ = **7.0 μm**. Rotor: Kn = **5.9 × 10⁻⁶**, continuum. Dust grain: Kn = **0.070**, transitional: its drag needs a slip correction and settling speeds computed from Stokes' law would be wrong.

(d) Above 5.03 km/s at 200 K: H atoms **1.56 × 10⁻³**, H₂ **1.0 × 10⁻⁶**. One H atom in 640 at the exobase is above escape speed, about 500 times the Earth's fraction for H₂ at 1000 K: Mars loses hydrogen efficiently, and hydrogen escaping from photo-dissociated water is a route by which Mars lost much of its water.

**Check.** a/v_rms = 0.656 = √(1.29/3): the same ratio formula.

**Common mistakes.** Using γ = 1.4 for CO₂; using the Earth's escape speed; using the Martian surface temperature instead of the exobase temperature for escape.

## 4. An airliner at 16 km ($s = 15$)

(a) T = **216.65 K** (isothermal layer), p = **10 287.6 Pa**, ρ = **0.16542 kg/m³**, a = **295.07 m/s**, σ = **0.1350**.

(b) V = 0.85 × 295.07 = **250.8 m/s**; q = ½ρV² = **5203 Pa**.

(c) Residual **−4.1 × 10⁻⁹**: the isothermal-layer formula satisfies dp/dh = −ρg₀ to the finite-difference error.

(d) At 80 km, λ = **5.04 mm**, Kn = **1.26 × 10⁻³** for the 4 m capsule: below 0.01, so the continuum model is still adequate for its overall forces (though thin regions such as the shock and the boundary layer are already rarefied). For ℓ = 4 m, Kn stays below 0.01 everywhere up to 86 km, so the ISA cannot answer where it crosses: the crossing is above the model's range, and a thermospheric model is needed.

**Check.** 11–20 km is isothermal: p(16 km)/p(11 km) = exp(−5000 × 9.80665/(287.05 × 216.65)) = 0.4545, and 22 632.3 × 0.4545 = 10 287.6 Pa.

**Common mistakes.** Using km in the ISA functions; using the sea-level speed of sound to convert Mach to true airspeed (that gives 289 m/s, 15 % too high).

## 5. How fast does a gas forget? (SIMULATION, seed 15)

(a)

| pair collisions per molecule | ⟨v⁴⟩/⟨v²⟩² | ⟨v⟩/v_rms | T_kin (K) |
|---|---|---|---|
| 0 (start) | 1.000 | 1.000 | 300.00 |
| 0.25 | 1.1029 | 0.9845 | 300.00 |
| 0.5 | 1.1844 | 0.9734 | 300.00 |
| 1 | 1.3163 | 0.9566 | 300.00 |
| 2 | 1.4888 | 0.9378 | 300.00 |
| 4 | 1.6235 | 0.9244 | 300.00 |
| 8 | 1.6619 | 0.9211 | 300.00 |

(b) Within 1 % of 5/3 at **8** pair collisions per molecule (−0.29 %); at 4 it is still −2.6 %.

(c) Exact Maxwell–Boltzmann samples of 20 000: kurtosis **1.6647 ± 0.0086** (mean ± standard deviation over 50 samples). The relaxed gas at 1.6619 is 0.3σ from that mean: **indistinguishable** at this sample size. (At 4 collisions per molecule, 1.6235 is 4.8σ low: distinguishable.)

(d) Each pair collision involves two molecules, so 8 pair collisions per molecule means about 16 collisions experienced by each molecule. At 7.1 × 10⁹ collisions per second: **about 2 ns**. Real air re-establishes its Maxwell–Boltzmann distribution within nanoseconds; that is why a single temperature describes it even in a sound wave or a shock layer many mean free paths thick.

**Common mistakes.** Calling the relaxation "proof" that Maxwell–Boltzmann is right (it shows that this collision model has it as its equilibrium); judging convergence without the sampling scatter; mixing pair collisions and collisions per molecule.

## 6. Judgment

(a) **The standard day.** At 11 km the real temperature may differ from ISA by 15–20 K; the air-data computer uses measured total temperature to get the true speed of sound, and ISA only to define pressure altitude. Kinetic theory itself is fine here.

(b) **Constant γ and a single species.** Behind the shock the gas reaches several thousand kelvin: vibration is excited (γ is already 5 % low at 1191 K in the Einstein model), O₂ and N₂ dissociate and some ionise. Use a high-temperature gas model with chemistry; the ideal-gas law per species still holds.

(c) **The continuum.** At 400 km the mean free path is far larger than the satellite (Kn ≫ 10): free-molecular flow, with drag from molecules hitting the surface directly (a drag coefficient from gas–surface interaction models) and a thermospheric density model, not ISA.

(d) **The ideal gas.** At 30 MPa and 90 K helium's molecules are close enough that non-ideal behaviour matters for the mass stored in the tank; use a real-gas equation of state (a compressibility factor) when sizing the tank.

(e) **Mostly fine, with a measurement caveat.** At 30 km the ISA is a good first model and the continuum holds for the balloon (Kn = 4.4 × 10⁻⁶ for a 1 m body), but the thin air couples weakly to the thermistor: sunlight heating it can bias the reading noticeably. The issue is heat transfer to the sensor, not kinetic theory.

## Check-yourself answers

1. p = (1/3)nm v_rms². The 1/3 is isotropy: only the velocity component normal to the wall carries momentum to it, and on average it holds a third of ⟨v²⟩.
2. (3/2)kT.
3. v_p < ⟨v⟩ < v_rms, in the ratios 1 : 1.128 : 1.225 (√2 : √(8/π) : √3).
4. At the exobase (≈ 1000 K) about 3 × 10⁻⁶ of H₂ molecules exceed the escape speed and are continuously replaced from the distribution's tail; for N₂ the fraction is 10⁻⁸⁵. Over billions of years hydrogen leaks away and nitrogen does not.
5. 5/3 for a monatomic gas (f = 3); 7/5 for air (f = 5: three translational and two rotational degrees of freedom; vibration frozen). More ways to store energy mean a larger C_v and a smaller γ.
6. 295.1 m/s at 216.65 K; a ∝ √T, and the air at 11 km is 71.5 K colder.
7. Kn = λ/ℓ, mean free path over body size; about 0.01.
8. Momentum and kinetic energy of the pair are conserved exactly by each collision; the speed distribution (and so ⟨v⁴⟩/⟨v²⟩²) is not, and relaxes to Maxwell–Boltzmann with kurtosis 5/3.
