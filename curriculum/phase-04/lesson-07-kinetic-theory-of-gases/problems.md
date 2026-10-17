# P04.07  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact, which numerical and which from a SIMULATION (label simulation output as such).

Several problems use your personal seed **s**: the day of the month you were born (1–31). Constants: k = 1.380649 × 10⁻²³ J/K, N_A = 6.02214076 × 10²³ mol⁻¹, M(N₂) = 28.0134 g/mol, M(air) = 28.9647 g/mol.

## 1. (Warm-up: nitrogen at your temperature)

Take N₂ at T = 200 + 5s K and p = 101 325 Pa.

(a) Compute v_p, ⟨v⟩ and v_rms.
(b) Compute n, and check that (1/3)nm v_rms² returns p.
(c) Compute the mean translational kinetic energy per molecule in joules and in electronvolts (1 eV = 1.602176634 × 10⁻¹⁹ J).
(d) What fraction of molecules move faster than 2v_rms?
(e) Compute the speed of sound in air at this temperature and the ratio a/v_rms for N₂ with γ = 1.4.

Answer form: eight numbers with units.

## 2. (Derivation: pressure, collision rate, the three speeds, the speed-of-sound ratio)

(a) Derive p = (1/3)nm⟨v²⟩, stating where you use isotropy and elastic wall collisions.
(b) Show that the number of molecules striking unit area of wall per unit time is n⟨v⟩/4. (Hint: the flux through the wall from molecules with speed v at angle θ to its normal is n f(v) v cos θ, averaged over the hemisphere of directions.)
(c) From the Maxwell–Boltzmann density, derive v_p, ⟨v⟩ and ⟨v²⟩. (You may use ∫₀^∞ x^{2k+1} e^{−x²} dx = k!/2.)
(d) Starting from a = √(γp/ρ), show that a/v_rms = √(γ/3), and explain why a must be less than v_rms for any γ allowed by equipartition.

Answer form: four derivations.

## 3. (Unfamiliar scenario: the air of Mars)

Near the Martian surface the atmosphere is mostly CO₂ (M = 44.01 g/mol) at about 210 K and 610 Pa; measured γ ≈ 1.29. Take a hard-sphere diameter of 3.9 × 10⁻¹⁰ m for CO₂.

(a) Compute the density, the speed of sound and v_rms. Compare the speed of sound with the Earth's at 288.15 K.
(b) How many active degrees of freedom does γ = 1.29 imply? What does that say about CO₂ at 210 K?
(c) Compute the mean free path and the Knudsen number for a 1.2 m helicopter rotor and for a 0.1 mm dust grain. Which is in the continuum regime?
(d) Mars' escape speed from the surface is 5.03 km/s (`P04.05`, Problem 3). Take the exobase temperature as 200 K and compute the fraction of H atoms (M = 1.008 g/mol) and of H₂ molecules above that speed. What does this suggest about the history of water on Mars?

Answer form: three numbers and a comparison, one number and a sentence, two Knudsen numbers, two fractions and a sentence.

## 4. (Numerical, uses your code: an airliner at altitude)

At geopotential altitude h = (1 + s) km in the ISA:

(a) Compute T, p, ρ and a with your ISA functions, and the density ratio σ = ρ/ρ₀.
(b) An airliner flies at Mach 0.85. Compute its true airspeed and the dynamic pressure ½ρV².
(c) Check the hydrostatic equation at this altitude with `hydrostatic_residual`.
(d) A 4 m re-entry capsule passes 80 km. Compute λ and Kn there. Is the continuum model adequate at 80 km for the capsule? At what altitude would Kn reach 0.01 for it, and can the ISA tell you?

Answer form: five numbers, two numbers, one residual, two numbers and two sentences.

## 5. (Simulation, uses your code: how fast does a gas forget?)

Use `numpy.random.default_rng(s)`. Start 20 000 N₂ molecules at the single speed v_rms(300 K) in random directions. This is a SIMULATION.

(a) Apply `relax_to_maxwellian` in stages so that the cumulative number of pair collisions per molecule is 0.25, 0.5, 1, 2, 4 and 8. After each, record ⟨v⁴⟩/⟨v²⟩², ⟨v⟩/v_rms and the kinetic temperature.
(b) How many pair collisions per molecule are needed to bring the kurtosis within 1 % of 5/3?
(c) Draw 50 exact Maxwell–Boltzmann samples of 20 000 speeds with `sample_speeds` (another seeded generator) and compute the mean and standard deviation of the kurtosis. Is your relaxed gas distinguishable from Maxwell–Boltzmann at this sample size?
(d) At sea level each molecule collides 7.1 × 10⁹ times per second. How long does real air take to "forget" a disturbance of its speed distribution?

Answer form: a 6-row table, one number, two numbers and a sentence, one time.

## 6. (Judgment: is the ideal-gas, ISA, continuum model good enough?)

For each case say which assumption fails first (if any), what the consequence is, and what you would use instead. Two sentences each.

(a) The air-data computer of a subsonic airliner at 11 km on a hot day.
(b) The gas behind the bow shock of a capsule re-entering at 7.5 km/s.
(c) Drag on a CubeSat at 400 km altitude.
(d) Helium pressurant in a 30 MPa composite tank at 90 K.
(e) A weather balloon at 30 km measuring temperature with a thermistor.

Answer form: five short paragraphs.
