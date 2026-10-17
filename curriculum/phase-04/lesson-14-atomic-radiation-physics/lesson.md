# P04.14  Atomic & Radiation Physics

> The same physics that makes a hydrogen lamp glow red, a plutonium pellet glow hot and a memory chip flip a bit is one idea: energy comes in packets, and a packet deposited in the wrong place does damage. A radiation budget is bookkeeping of those packets, done before launch.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the Bohr energy levels of hydrogen and the Rydberg formula, compute line wavelengths and photon energies, and use photon counting to size a lidar pulse. 2. Compute nuclear binding energies and decay Q-values from atomic masses, compare them with the semi-empirical mass formula, and explain from the binding-energy curve why heavy nuclei alpha-decay. 3. Derive the half-life law, activity, specific activity and decay heat, solve a two-member decay chain (Bateman) and check it against `solve_ivp`, and size the Pu-238 inventory of an RTG for an end-of-mission power with its half-life decline. 4. Distinguish absorbed dose (Gy) from equivalent dose (Sv), and total ionizing dose, displacement damage and single-event effects in electronics; compute a TID with a design margin, an upset rate and a Poisson probability. 5. Fit a first-order exponential shielding model to two dose-depth points, compute the aluminium thickness for a dose target, and state where the model fails. |
| **Prerequisites** | `P04.13` (relativistic energy; E = mc² for the mass defect), `P04.11` (photons and the electromagnetic spectrum), `P04.02` (energy bookkeeping, power), `P02.01` (first-order linear ODEs: the decay law), `P02.04` (systems of ODEs: the chain), `P02.16` (the Poisson distribution), `P03.05` (what `solve_ivp` does). Skill check: you can solve dN/dt = −λN and say what e^{−λt} is after one half-life. |
| **Engineering question** | A deep-space probe needs 110 W of electrical power after 14 years and carries electronics rated for 30 krad(Si): how much Pu-238 must the RTG carry, how does its power decline over the mission, and how much aluminium must surround the electronics if a dose-depth curve gives 20 krad/yr behind 1 mm and 4 krad/yr behind 3 mm and the design margin is 2? |
| **Assumptions and limits** | Hydrogen: the Bohr model with the reduced-mass correction (exact for the energy levels of a one-electron atom without fine structure; wrong for multi-electron atoms). Nuclei: atomic masses from the published mass evaluation; the semi-empirical formula is a model good to about 1 % for A > 30 and poor for light nuclei. Decay: each nucleus decays independently at a constant rate (exact for radioactive decay). RTG: all alpha energy is deposited in the fuel, the conversion efficiency is an assumed constant (6 %) with an optional assumed degradation rate; real thermoelectric performance also depends on the hot- and cold-side temperatures. Dose: thin-target deposition for the fluence formula. Shielding: **a model** — the dose behind an aluminium slab is taken to fall exponentially with areal density, with an attenuation length fitted to two points of a dose-depth curve; this is exact only for a narrow beam of monoenergetic photons, holds over the fitted range for the trapped-electron part of a spectrum, and fails for high-energy protons and galactic cosmic rays, which barely attenuate and produce secondaries. The two dose-depth points are illustrative, not the output of an environment model. |
| **Analytical work** | Bohr levels from the quantized angular momentum and the Coulomb force; the Rydberg formula with the reduced mass; photon energy and photon count; mass defect and binding energy; the Q-value of Pu-238 → U-234 + α and its split between the α and the recoil; the decay law, activity, specific activity, decay heat; the Bateman solution for a two-member chain and its secular-equilibrium limit; dose from fluence and stopping power; TID with margin; the Poisson probability of at least one upset; the RTG inventory; the exponential model, its two-point fit and the thickness for a target. |
| **Computational deliverable** | `starter/main.py`: `energy_level`, `photon_energy_ev`, `transition_wavelength`, `photons_per_pulse`, `lidar_range`, `binding_energy_mev`, `semf_binding_energy_mev`, `q_value_mev`, `decay_constant`, `remaining_fraction`, `activity_bq`, `specific_activity_bq_per_g`, `decay_heat_w_per_g`, `bateman_two`, `chain_numeric`, `absorbed_dose_gy`, `equivalent_dose_sv`, `dose_from_fluence_gy`, `gy_to_krad`, `total_ionizing_dose`, `seu_rate_per_day`, `prob_at_least_one`, `rtg_electrical_power`, `rtg_fuel_mass_g`, `areal_density_g_cm2`, `shielded_dose`, `fit_attenuation`, `areal_density_for_dose`. Tests prove: the 1/n² law, −13.6057 and −13.5983 eV, the Rydberg formula to 10⁻⁸ and the Lyman-α and H-α wavelengths; hc = 1239.84 eV nm and the photon count of a 10 mJ pulse; the binding energy of He-4 (28.30 MeV) and of Fe-56 per nucleon (8.790 MeV), the peak of the curve at Fe/Ni; the liquid-drop model within 1 % for heavy nuclei and not for He-4; Q(Pu-238) = 5.593 MeV, equal to the binding-energy gain; the half-life law (½, ⅛, memoryless); 6.336 × 10¹¹ Bq/g and 0.568 W/g for Pu-238; the Bateman solution equal to `solve_ivp` to 10⁻⁷ and its equilibrium ratio λ₂/(λ₂ − λ₁); dose units and weighting; TID and SEU scaling; the Poisson probability against 20 000 seeded draws; the RTG halving every 87.7 years and the fuel-mass round trip; zero-thickness shielding returning the model's D₀, monotone decay, the fit round trip and the 4.22 mm answer; the photon half-value layer; the figure. |
| **Practical activity** | **computational experiment** (no sources, no detectors at home): an RTG sizing study over mission length and assumed degradation; a shielding trade over the dose target and the sensitivity of the answer to the fitted points; an upset-rate and scrubbing study; a comparison of Pu-238 with Am-241 as a heat source. |
| **Equipment** | None. A computer with the environment from `P00.01`. No radioactive source, detector, smoke-detector disassembly or "radiation experiment" of any kind is part of this lesson. |
| **Safety** | None beyond normal computer use. Do not open smoke detectors, buy sources or build detectors for this lesson. |
| **Evidence** | Completed `main.py`, passing tests, `radiation_dose.png`, the hand calculations (H-α, Q-value, specific power, RTG inventory, shielding thickness), the four study tables, and a half-page design note giving the RTG inventory and the shielding thickness with the margin, the assumptions, and the single input that most changes each answer. |
| **Assessment** | `assessments/phase-04/grading-rubric.md`; competencies 1–4. Pass: the Bateman solution and the decay heat are derived; all public tests pass; the note states that the shielding law is a model, gives its range of validity, and shows the sensitivity of the thickness to the fitted points. |
| **References** | `@openstax-university-physics-vol3` Ch. 8 "Atomic Structure": 8.1 The Hydrogen Atom; 8.4 The Exclusion Principle and the Periodic Table; 8.5 Atomic Spectra and X-rays; 8.6 Lasers. Ch. 10 "Nuclear Physics": 10.1 Properties of Nuclei; 10.2 Nuclear Binding Energy; 10.3 Radioactive Decay; 10.4 Nuclear Reactions; 10.7 Medical Applications and Biological Effects of Nuclear Radiation (required-free): the levels, the spectra, binding energy, the decay law, gray and sievert. `@nasa-why-space-radiation-matters` "What is Space Radiation?", "Where Does Radiation Come From?", "What is Non-Ionizing versus Ionizing Radiation?", "What Is Galactic Cosmic Radiation?", "What Factors Determine the Amount of Radiation Astronauts Receive?" (required-free): the environment. `@ntrs-pellish-radiation-effects-tutorial` "What are radiation effects?", "What is total ionizing dose? (TID)", "What is displacement damage? (DDD)", "What is NIEL?", "What are single-event effects?", "Electronics for Space Use" (required-free): effects on electronics. `@doe-what-is-radioisotope-power-system` "Types of Radioisotope Systems", "How Does an RTG Work?" and `@nasa-basics-of-space-flight` Chapter 11 "Onboard Systems" page 3, "Radioisotope Power Systems" (required-free): RTGs in context. `@nasa-hdbk-4002b-in-space-charging` (supplementary-free): surface and internal charging, the other radiation hazard to electronics, not treated here. |
| **Extensions** | ★ Extend `bateman_two` to an n-member chain and apply it to the first three steps of the U-238 series. ★ Add a buildup factor B(x) = 1 + x/L to the photon attenuation law and see how much the thickness changes. ★ Replace the two-point fit by a least-squares fit to five points of a dose-depth curve you construct from two exponentials (a soft electron part and a hard proton part) and show the one-exponential model's error outside the fitted range. |
| **Time** | unmeasured |

## 1. The engineering problem

Beyond Mars, sunlight is too weak for solar panels of reasonable size, and the probes that went there (and to the poles of the Moon, and to the surface of Mars through dust storms) took their power from a lump of plutonium-238. Pu-238 emits alpha particles; they stop within micrometres in the fuel, and their energy becomes heat; thermocouples across the temperature difference between the hot fuel and the cold radiator convert a few percent of it to electricity. The heat falls with the half-life of 87.7 years, so a probe that needs 110 W after 14 years must launch with more. How much more, and how many kilograms of a scarce, specially produced isotope that means, is a calculation from the decay law.

The same probe flies through radiation: electrons and protons trapped in planetary magnetic fields, protons from solar storms, and galactic cosmic rays. In its electronics this does three different kinds of damage (Pellish's tutorial is the map): a slow accumulation of trapped charge in insulators (**total ionizing dose**, TID), a slow accumulation of displaced atoms in the crystal (**displacement damage**), and sudden events where one particle deposits enough charge in one transistor to flip a bit or, worse, switch on a parasitic short circuit (**single-event effects**: an upset, SEU, or a latch-up, SEL). The first is fought with shielding and with parts rated for a dose; the third with error-correcting memory and current limiters, because no practical shielding stops a cosmic-ray iron nucleus.

This lesson builds both budgets from the physics up: atoms and photons (lidar is a laser rangefinder and its link budget counts photons), nuclei and their binding energy (why Pu-238 decays at all, and how much energy each decay releases), the decay law (power over the mission), dose and its units, and a first-order shielding model. Success criteria: the hydrogen lines to four figures; Q(Pu-238) from masses to the keV; the specific power of Pu-238 to three figures; the RTG inventory for the requirement; the shielding thickness with its sensitivity; and a clear statement, in the note, of what the shielding model cannot do.

## 2. Understand and derive

### 2.1 The hydrogen atom

Bohr's model (OpenStax 8.1 gives the full quantum treatment; the energies agree) puts the electron on a circular orbit where the Coulomb force supplies the centripetal force, mv²/r = ke²/r², and quantizes the angular momentum, mvr = nħ. Eliminating v gives r_n = n²a₀ and

**E_n = −13.6057 eV / n²**,   (1)

the Rydberg energy R_y = mₑe⁴k²/(2ħ²) over n². The nucleus is not infinitely heavy; both particles orbit their centre of mass, and mₑ is replaced by the reduced mass μ = mₑ/(1 + mₑ/m_p), which lowers R_y by 1 part in 1836: the hydrogen ground state is **−13.5983 eV**. A transition from n_u to n_l emits a photon of energy E_u − E_l = hc/λ, so 1/λ = R_H(1/n_l² − 1/n_u²), the Rydberg formula (8.5), with R_H = 1.096 776 × 10⁷ m⁻¹.

| line | transition | λ (vacuum) | photon energy |
|---|---|---|---|
| Lyman α | 2 → 1 | 121.57 nm | 10.20 eV |
| H-α | 3 → 2 | 656.47 nm | 1.889 eV |
| H-β | 4 → 2 | 486.27 nm | 2.550 eV |
| Balmer limit | ∞ → 2 | 364.71 nm | 3.400 eV |
| Lyman limit | ∞ → 1 | 91.18 nm | 13.598 eV |

(Tables of the red hydrogen line often say 656.28 nm: that is the wavelength in air, whose refractive index is 1.00028.) Lyman α is the line that ultraviolet space instruments map the hydrogen of planetary exospheres with; H-α is what a solar telescope filters on to watch flares coming.

A **laser** (8.6) is stimulated emission in a medium with a population inversion: each photon triggers the emission of an identical one, so the output is coherent, nearly monochromatic and collimated. For **lidar** (a laser rangefinder or altimeter) two numbers matter. The range is R = ct/2 from the round-trip time t: 6.67 µs is 1000 m. And the photon budget: a 1064 nm Nd:YAG photon carries hc/λ = **1.165 eV**, so a 10 mJ pulse is **5.356 × 10¹⁶ photons**; after beam spreading, a diffuse reflection from the ground and a small receiving telescope, only a tiny fraction of them comes back (the inverse-square law of `P04.11` applies to the return), which is why sensitive lidar receivers are built to detect single photons.

### 2.2 Nuclei and binding energy

A nucleus of Z protons and N = A − Z neutrons weighs less than its parts; the difference, times c² (`P04.13`), is the **binding energy** (10.2):

**B = (Z m_H + N m_n − M_atom) c²**,   (2)

using the hydrogen-atom mass m_H so that the Z electron masses cancel against those in the atomic mass M_atom; 1 u = 931.494 MeV. From published masses: He-4 28.30 MeV (7.07 MeV per nucleon, unusually tight for its size), C-12 7.68, O-16 7.98, **Fe-56 8.790**, Ni-62 8.795, U-238 7.570 MeV per nucleon. The curve rises steeply, peaks at iron and nickel, and falls slowly: light nuclei release energy by fusing, heavy ones by splitting or by shedding a tightly bound He-4.

The **semi-empirical mass formula** is the liquid-drop model of that curve: a volume term a_V A, minus a surface term a_S A^{2/3}, minus the Coulomb repulsion a_C Z(Z − 1)/A^{1/3}, minus an asymmetry term a_A(A − 2Z)²/A, plus a pairing term. With a common textbook set of coefficients (15.75, 17.8, 0.711, 23.7, 11.18 MeV) it gives 8.846 MeV per nucleon for Fe-56 (+0.6 %) and 7.625 for U-238 (+0.7 %), but 5.71 for He-4 against 7.07: a model of a drop is a poor model of four particles. It is the reason the code computes B from masses and uses the formula only as a check.

**Alpha decay of Pu-238.** Pu-238 → U-234 + He-4. The **Q-value** is the mass difference times c²:

Q = (238.049 560 − 234.040 952 − 4.002 603) u × 931.494 MeV/u = **5.593 MeV**.

It equals the gain in total binding energy (the tests compute it both ways). Momentum conservation splits it: the alpha carries Q × 234/238 = 5.50 MeV and the uranium recoil 0.094 MeV. Both stop within micrometres in an oxide pellet, so all 5.593 MeV becomes heat.

### 2.3 The decay law

Each nucleus decays independently with probability λ dt in time dt (10.3), so dN/dt = −λN and

**N(t) = N₀ e^{−λt} = N₀ 2^{−t/T½},  λ = ln 2 / T½**.   (3)

The law is memoryless: the fraction surviving t₁ + t₂ is the product of the fractions surviving each. The **activity** A = λN counts decays per second (becquerel, Bq; the old unit curie is 3.7 × 10¹⁰ Bq). One gram of a pure nuclide of molar mass M contains N_A/M atoms, so its **specific activity** is λN_A/M; for Pu-238 (T½ = 87.7 y = 2.768 × 10⁹ s, λ = 2.505 × 10⁻¹⁰ s⁻¹) it is **6.336 × 10¹¹ Bq/g = 17.1 Ci/g**. Multiplying by Q:

**decay heat = λ (N_A/M) Q = 0.568 W/g**,   (4)

about 568 W per kilogram. The same formula explains the choice of isotope: U-238 (T½ = 4.47 × 10⁹ y) has 5 × 10⁷ times less power per gram; Po-210 (T½ = 138 d) has far more but would be gone within a year; Am-241 (T½ = 432.6 y, Q = 5.638 MeV) gives 0.115 W/g, five times less than Pu-238 but declining five times more slowly. A half-life of decades is the sweet spot for a mission of decades.

### 2.4 Decay chains

When the daughter is itself radioactive (10.3, 10.4), N₁′ = −λ₁N₁ and N₂′ = λ₁N₁ − λ₂N₂. With N₂(0) = 0 the solution (the **Bateman solution**) is

**N₂(t) = N₁(0) λ₁/(λ₂ − λ₁) (e^{−λ₁t} − e^{−λ₂t})**,   (5)

found by the integrating factor of `P02.01` (Problem 2). The daughter activity λ₂N₂ rises, peaks at t_max = ln(λ₂/λ₁)/(λ₂ − λ₁), and then, if the parent lives longer (λ₁ < λ₂), settles into a constant ratio λ₂N₂/(λ₁N₁) → λ₂/(λ₂ − λ₁). For a very long-lived parent this ratio → 1: **secular equilibrium**, the daughter decaying as fast as it is produced. For half-lives 2 and 0.3 (any unit), t_max = 0.966. U-234 from Pu-238 has a half-life of 245 000 years, so it is effectively stable over a mission and the RTG's heat is Pu-238's alone.

### 2.5 Dose and its units

**Absorbed dose** is energy deposited per unit mass: 1 gray = 1 J/kg (10.7). Electronics people still use the rad (0.01 Gy), so 1 krad = 10 Gy, and they specify the material because the deposit depends on it: krad(Si). For a thin layer crossed by a fluence Φ (particles per cm²) of particles each losing S MeV cm²/g (the mass stopping power), the dose is ΦS in MeV/g, and 1 MeV/g = 1.602 × 10⁻¹⁰ Gy. For living tissue, the **equivalent dose** in sieverts weights each radiation by its biological effect: w_R = 1 for photons and electrons, 2 for protons, 20 for alphas. An alpha dose of 1 mGy is 20 mSv. (The sievert is for people; electronics budgets are always in gray or rad.)

The space environment (NASA, "Why Space Radiation Matters") supplies three sources: particles trapped in the magnetic field (electrons and protons of the radiation belts), solar particle events, and galactic cosmic rays (protons and heavy nuclei from outside the solar system, very energetic and very penetrating). The Earth's atmosphere and field protect the ground from most of it; a spacecraft has only its own walls.

### 2.6 Three ways radiation breaks electronics

From Pellish's tutorial, in this lesson's words:

- **Total ionizing dose.** Ionization in the oxide layers of transistors leaves trapped positive charge, which shifts threshold voltages and increases leakage. The effect is cumulative and gradual; a part is rated by the dose at which it stops meeting specification (commercial parts often a few krad to tens of krad, radiation-hardened parts hundreds). The design rule: the dose the part will receive behind its shielding, times a **radiation design margin** (commonly 2), must be below its rating: TID_design = dose rate × years × margin.
- **Displacement damage.** A proton or neutron knocks silicon atoms out of the lattice; the defects shorten carrier lifetimes. Solar cells, optical sensors and bipolar transistors suffer most. It scales with the non-ionizing energy loss (NIEL), not with the ionizing dose, so it needs its own budget.
- **Single-event effects.** One energetic particle leaves a dense trail of charge; if a transistor node collects enough, a memory cell flips (SEU, a "soft" error, harmless once rewritten) or a parasitic thyristor turns on (SEL, a "hard" error that draws current until the power is cycled and can burn the part). The rate is a cross-section σ (cm² per bit) times a flux of upsetting particles: upsets per day = σ φ N_bits × 86 400. Because upsets are independent random events, the number in time t is Poisson with mean rt (`P02.16`), and the probability of at least one is 1 − e^{−rt}. With illustrative values σ = 10⁻¹⁴ cm²/bit, φ = 1 cm⁻² s⁻¹ and 1 GiB (8.59 × 10⁹ bits), r = 7.42 per day and the chance of at least one upset in an hour is 0.266. The defence is error-correcting code plus **scrubbing** (periodically reading and rewriting memory) so that two upsets never accumulate in one word.

Shielding helps a great deal against TID from trapped electrons, somewhat against protons, and almost not at all against the heavy ions that cause most single-event effects.

### 2.7 The RTG budget

The thermal power is P_th(t) = P_th(0) 2^{−t/T½}; the electrical power is P_e = η P_th with η the thermoelectric conversion efficiency. η here is an **assumed** 6 % (thermoelectric generators convert a few percent; the DOE article describes the principle, not this number); the thermocouples also age, which the model represents with an assumed fractional loss d per year: P_e(t) = η P_th(0) 2^{−t/T½}(1 − d)^t. Requiring P_e(14 y) = 110 W with d = 0:

P_th(0) = 110 / (0.06 × 2^{−14/87.7}) = 110 / (0.06 × 0.8953) = **2048 W**, mass = 2048 / 0.568 = **3.61 kg of Pu-238**,

which delivers 122.9 W at launch. With d = 1.5 %/year (assumed) the requirement grows to 4.46 kg: over 14 years the ageing costs more than the decay. Real fuel is plutonium dioxide containing other plutonium isotopes, so the mass of fuel is larger than the mass of Pu-238; the model counts only the Pu-238.

| years | 0 | 5 | 10 | 14 | 20 | 30 |
|---|---|---|---|---|---|---|
| P_e, decay only (W) | 122.9 | 118.1 | 113.5 | 110.0 | 104.9 | 96.9 |
| P_e, + 1.5 %/y (W) | 122.9 | 109.5 | 97.6 | 89.0 | 77.5 | 61.6 |

### 2.8 Shielding: a first-order model

Real shielding calculations transport a particle spectrum through the geometry with a Monte Carlo code and produce a **dose-depth curve**: dose behind a spherical shell of aluminium against its areal density ρx (g/cm²; 1 mm of aluminium is 0.27 g/cm²). This lesson does not do transport. It uses a model:

**D(x) = D₀ e^{−x/L}**,   (6)

with x the areal density and L an effective attenuation length. For a narrow beam of monoenergetic photons (6) is exact with L = 1/(μ/ρ): at 1 MeV in aluminium μ/ρ ≈ 0.0615 cm²/g, so the half-value layer is ln 2/0.0615 = 11.3 g/cm² (4.2 cm). For the electron-dominated part of a trapped-particle spectrum, dose-depth curves fall roughly exponentially over the first few millimetres, so (6) fitted to two points is a usable interpolation there and nowhere else. Fitting through (x₁, D₁) and (x₂, D₂): L = (x₂ − x₁)/ln(D₁/D₂), D₀ = D₁e^{x₁/L}. The thickness for a target D_t is x = L ln(D₀/D_t), and zero if D_t ≥ D₀.

**Worked design.** Illustrative points (not from an environment model): 20 krad/yr behind 1 mm (0.27 g/cm²) and 4 krad/yr behind 3 mm (0.81 g/cm²). Then L = 0.54/ln 5 = **0.3355 g/cm²** (a half-value layer of 0.86 mm) and D₀ = **44.7 krad/yr**. The part is rated 30 krad, the mission is 10 years, the margin 2: the allowed rate is 30/(10 × 2) = 1.5 krad/yr, so x = 0.3355 ln(44.7/1.5) = **1.139 g/cm² = 4.22 mm of aluminium**. A 10 cm cube of electronics has 600 cm² of wall, so that shield weighs 0.68 kg. If the 3 mm point were 20 % higher (4.8 krad/yr), L would be 0.378 g/cm² and the answer 4.63 mm; 20 % lower, 3.83 mm. The thickness moves by about half the perturbation, because it depends on the dose only through a logarithm and the target lies just beyond the fitted range; extrapolating much further would amplify the uncertainty of L, which two points barely constrain.

Where (6) fails: D₀ is not the unshielded dose (the softest electrons stop in the first fraction of a millimetre, so the true curve at x → 0 is far above D₀); beyond the fitted range the dose stops falling exponentially because the hard proton component and bremsstrahlung photons dominate, so (6) **underestimates** the dose at large thickness; heavy ions and secondaries are not in it at all. A real design uses the dose-depth curve of a validated environment model for the actual orbit and mission dates.

## 3. Build the computational model

Open `starter/main.py`. The constants and the mass table are given.

### Step 1: atoms and photons

`energy_level(n, Z, reduced_mass)` from (1) with arrays allowed and `ValueError` for n < 1; `photon_energy_ev`; `transition_wavelength(n_upper, n_lower)`; `photons_per_pulse`; `lidar_range`.

### Step 2: nuclei

`binding_energy_mev(Z, A, mass_u)` from (2); `semf_binding_energy_mev(Z, A)` with the pairing term; `q_value_mev(parent, *products)`.

### Step 3: decay

`decay_constant`, `remaining_fraction`, `activity_bq`, `specific_activity_bq_per_g`, `decay_heat_w_per_g` from (3) and (4); `bateman_two` from (5) with a `ValueError` for equal constants; `chain_numeric` with `solve_ivp(method="Radau", rtol=1e-10)` as the oracle.

### Step 4: dose, TID, SEU

`absorbed_dose_gy`, `equivalent_dose_sv(dict)`, `dose_from_fluence_gy`, `gy_to_krad`, `total_ionizing_dose`, `seu_rate_per_day`, `prob_at_least_one` with `expm1` so that tiny rates do not cancel.

### Step 5: RTG and shielding

`rtg_electrical_power(t, p_th0, efficiency, half_life, degradation)`, `rtg_fuel_mass_g`; `areal_density_g_cm2`, `shielded_dose` (6), `fit_attenuation`, `areal_density_for_dose`.

Running `main()` prints

```
hydrogen ground state: -13.5983 eV (reduced mass), -13.6057 eV (infinite mass)
  Lyman alpha    121.57 nm
  H-alpha        656.47 nm
  H-beta         486.27 nm
Nd:YAG lidar 1064 nm: photon 1.165 eV, 10 mJ pulse = 5.356e+16 photons; 6.67 us round trip = 1000 m
  He-4    B =    28.30 MeV, B/A = 7.074 MeV; SEMF B/A = 5.710
  Fe-56   B =   492.26 MeV, B/A = 8.790 MeV; SEMF B/A = 8.846
  U-238   B =  1801.69 MeV, B/A = 7.570 MeV; SEMF B/A = 7.625
  Pu-238  B =  1801.27 MeV, B/A = 7.568 MeV; SEMF B/A = 7.621
Pu-238 alpha decay Q = 5.594 MeV; specific activity 6.336e+11 Bq/g = 17.12 Ci/g; decay heat 0.5678 W/g
  after  0 y: 100.00 % of the heat
  after 14 y:  89.53 % of the heat
  after 17 y:  87.43 % of the heat
  after 50 y:  67.36 % of the heat
Pu-238 for 110 We at end of 14 y, 6 % efficiency: 3.61 kg; with 1.5 %/y degradation: 4.46 kg; thermal power at launch 2048 W, electrical 122.9 W
shielding MODEL (illustrative dose-depth points 20 and 4 krad/yr behind 1 and 3 mm Al): D0 = 44.7 krad/yr, L = 0.3355 g/cm^2; 30 krad part, 10 y, margin 2 -> 1.139 g/cm^2 = 4.22 mm Al
SEU (illustrative): sigma 1e-14 cm^2/bit, 1 /cm^2/s, 1 GiB -> 7.42 upsets/day; P(>=1 in 1 h) = 0.266
```

(the Q-value prints as 5.594 with three decimals of the unrounded 5.5935; the text's 5.593 is the same number) and writes `radiation_dose.png`: RTG power over 30 years with and without the assumed degradation; the shielding model on a log scale with its two fitted points; the activities of a two-member chain.

### Verify

Run `python -m pytest curriculum/phase-04/lesson-14-atomic-radiation-physics`. The tests check, in order: the Bohr levels (1/n², the two ground-state values, Z² scaling, the n = 0 error); the Rydberg formula to 10⁻⁸ for four lines and the Lyman-α and H-α values; hc = 1239.84 eV nm, the photon count and its energy round trip, the lidar range; B = 0 for a lone proton, 28.30 MeV for He-4, 8.790 MeV/nucleon for Fe-56, the maximum at Fe or Ni and the ordering U < Fe, He < C; the liquid-drop model within 1 % for Fe-56, U-238 and Pu-238 and more than 10 % off for He-4; Q(Pu-238) = 5.593 MeV, equal to the binding-energy gain, and negative for an impossible alpha decay of Fe-56; the half-life law; specific activity and decay heat; the Bateman chain against `solve_ivp` to 10⁻⁷, its equilibrium ratio and secular limit, and the equal-constant error; dose units, weighting, fluence linearity; TID and SEU scaling; the Poisson probability, its small-rate limit to 10⁻⁶ relative, and 20 000 seeded draws (tolerance five standard errors); the RTG halving at 87.7 and 175.4 years, the 0.8953 factor at 14 years and the fuel round trip; the shielding model at zero thickness, its monotone decay and constant ratio per step, the fit round trip, the 4.22 mm design, the zero-thickness answer for an easy target and the photon half-value layer; and the figure.

### Compare

Three independent checks are built in. Q from masses equals Q from binding energies: two routes through the same table that would disagree if (2) had the wrong electron bookkeeping. The Bateman formula against `solve_ivp`: an analytic solution and a numerical one with no shared code. The liquid-drop formula against the measured binding energies: a model against data, which agrees to 1 % where the model should hold and fails by 20 % where it should not. For the decay heat, the number 0.568 W/g is also the one quoted in RTG literature for pure Pu-238, which is a check on the half-life and the Q-value together.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** There is no physical activity in this lesson: no sources, no detectors, no opened smoke detectors. Write each prediction first.

1. **RTG sizing over mission length.** For missions of 5, 10, 14, 20, 30 years, compute the Pu-238 mass for 110 W at end of mission with d = 0, 1 % and 2 % per year. Predict that with d = 0 the mass grows only about 20 % from 5 to 30 years (the half-life is long), and that degradation, not decay, dominates the growth for any d above about 0.8 %/year (ln 2/87.7 = 0.79 %/year is the decay's own fractional rate).
2. **Shielding trade.** With the lesson's two points, compute the thickness for part ratings of 10, 30, 100, 300 krad (10 years, margin 2). Predict that every factor of 10 in rating saves the same thickness, L ln 10 = 0.77 g/cm² = 2.86 mm, and that a 400 krad part needs only the 1 mm of the first fitted point (20 krad/yr × 10 y × 2). The model would give x = 0 at 894 krad (D₀ × 20), but that lies outside the fitted range and means nothing.
3. **Sensitivity.** Perturb the 3 mm point by ±10 % and ±20 % and recompute L and the 30 krad thickness. Predict whether the thickness changes by more or less than the perturbation, given that it depends on the dose through a logarithm and the target lies just beyond 3 mm.
4. **Upsets and scrubbing.** For a memory of 4 Gbit with σ = 10⁻¹⁴ cm²/bit and φ = 1 cm⁻² s⁻¹ (illustrative), compute the daily upset count and, for 64-bit words protected by single-error-correcting code, the probability that some word receives two upsets within a scrub interval of 1 s, 1 min, 1 h, 1 day (approximately N_words (r_w T)²/2). Predict that it scales as T².
5. **Pu-238 against Am-241.** Repeat item 1 with Am-241 (molar mass 241.0568 g/mol, Q from masses with Np-237 at 237.048 173 u). Predict a fuel mass about 4.5 times larger at 14 years, with the ratio shrinking for longer missions.

Raw data format: CSV tables `mission_y, degradation_per_y, pu238_kg` (1), `rating_krad, areal_density_g_cm2, thickness_mm` (2), `perturbation, L_g_cm2, thickness_mm` (3), `scrub_s, p_double` (4), `isotope, mission_y, mass_kg` (5); plus `radiation_dose.png`.

What you should find: 3.36, 3.49, 3.61, 3.78, 4.09 kg for 5–30 years with d = 0 (+22 %), and 7.50 kg at 30 years with 2 %/year, nearly double; 5.58, 4.22, 2.72, 1.36 mm for 10, 30, 100, 300 krad, 2.86 mm per decade of rating; a ±20 % change in the 3 mm point moves the 30 krad thickness from 4.22 mm to 4.63 mm (+10 %) and to 3.83 mm (−9 %), ±10 % to 4.42 and 4.02 mm: about half the perturbation; the double-upset probability rising by 3600² from 1 s to 1 h; Am-241 at 16.4 kg against Pu-238 at 3.61 kg for 14 years (a factor 4.5), and 16.7 against 3.93 kg for 25 years (4.2).

## 5. Validate and improve

| quantity | hand (§2) | code | reference value | difference |
|---|---|---|---|---|
| H ground state | −13.598 eV | −13.5983 eV | ionization energy of H, 13.598 eV | < 10⁻⁴ |
| H-α (vacuum) | 656.47 nm | 656.47 nm | 656.28 nm in air × 1.00028 | 0 |
| B(He-4) | 28.30 MeV | 28.30 MeV | — | 0 |
| B/A(Fe-56) | 8.790 MeV | 8.790 MeV | SEMF 8.846 | 0.6 % (model) |
| Q(Pu-238 → U-234 + α) | 5.593 MeV | 5.593 MeV | binding-energy route: same | 10⁻⁹ |
| Pu-238 specific power | 0.568 W/g | 0.5678 W/g | commonly quoted ≈ 0.57 W/g | 0 |
| Bateman N₂(t) | (5) | `solve_ivp` | — | < 10⁻⁷ |
| RTG inventory, 14 y | 3.61 kg | 3.61 kg | — | 0 |
| shielding thickness (model) | 4.22 mm | 4.22 mm | none: illustrative inputs | — |

Everything with a definite answer is reproduced to the figures shown. The shielding row has no reference value, and it should not: it is a model through two illustrative points, and its honest uncertainty is the sensitivity found in experiment 3 plus the systematic error of extrapolating an exponential. Report the thickness as "4.2 mm from a two-point exponential fit, ±0.4 mm for ±20 % in the 3 mm dose, model error not included", never as a bare number.

What to improve: for the RTG, the efficiency and its ageing are the weak inputs (the decay is known to four figures); bound them, don't refine the half-life. For shielding, the model is the weak link: replace it with a dose-depth curve from a validated environment model, then account for the actual shielding geometry (the spacecraft structure around a box adds areal density in most directions). For upsets, the cross-section must come from a heavy-ion test of the actual part; the illustrative value here only shows how the numbers combine.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the hand derivations of (1), (2), the Q-value, (4), (5) and the RTG inventory; `radiation_dose.png`; the five study tables; and a half-page design note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6, 7) giving the Pu-238 inventory and the shielding thickness, each with its assumptions and the single input that most changes it.

### Adversarial self-checks

- **What if the mission doubles from 14 to 28 years?** With d = 0, the decay factor goes from 0.895 to 0.801, so the inventory rises 12 %; with d = 1.5 %/year it rises by another factor 1.24 from ageing. The shielding thickness for the same part rises by L ln 2 = 0.23 g/cm² (0.86 mm), because the allowed rate halves.
- **What if the shield thickness goes to zero?** The model gives D₀ = 44.7 krad/yr, but the real unshielded dose is much larger because the softest electrons are stopped in the first fraction of a millimetre: the model is not valid below the first fitted point. A zero-thickness test checks only the code, not the physics.
- **Which assumption fails first in the real world?** For shielding, the single exponential: beyond the fitted range the hard proton component makes the dose fall more slowly than (6) predicts, so adding thickness buys less than the model promises. For the RTG, the constant efficiency.
- **What measurement could falsify the model?** A dosimeter inside the box reading more than the model at the design thickness; or an RTG whose power falls faster than 0.79 %/year, which the half-life cannot explain and degradation must.
- **Does the result survive a change in input?** The RTG answer changes by 0.1 % per 0.1 year of half-life and inversely with efficiency; the shielding answer changes by +10 %/−9 % for ±20 % in one fitted point (+0.41/−0.39 mm).
- **Which decision changes if the dose uncertainty were twice as large?** At +40 % in the 3 mm point the 30 krad part needs 5.07 mm (0.82 kg instead of 0.68 kg for a 10 cm box): the logarithm makes thickness forgiving of dose uncertainty. What it does not forgive is the model itself; if the real curve flattens beyond 3 mm, no thickness of aluminium reaches the target, and a 100 krad part (2.72 mm in the model) is the robust design: buy hardness, not aluminium.

## Check yourself

1. Write the Bohr energy of level n and the wavelength of the 3 → 2 line of hydrogen. Why does the reduced mass matter?
2. How many photons are in a 10 mJ pulse at 1064 nm, and what range corresponds to a 6.67 µs round trip?
3. Define binding energy and give its value per nucleon for Fe-56. Why does Pu-238 alpha-decay?
4. Write the decay law and the half-life in terms of λ. What fraction of an RTG's heat remains after 14 years of Pu-238 decay?
5. Compute the specific power of Pu-238 from its half-life, molar mass and Q.
6. What is the difference between a gray and a sievert? What is 1 krad in gray?
7. Name the three classes of radiation effects in electronics and say which one shielding helps least.
8. In the exponential shielding model, what thickness halves the dose, and why is the model's value at zero thickness not the unshielded dose?
