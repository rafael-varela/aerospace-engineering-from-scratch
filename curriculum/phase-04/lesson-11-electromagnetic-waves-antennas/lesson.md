# P04.11  Electromagnetic Waves & Antennas

> A watt radiated from orbit spreads over a sphere two thousand kilometres across; the ground station catches a few hundred attowatts of it, and that is enough. The link budget is the bookkeeping of where the rest went.

| | |
|---|---|
| **Learning outcomes** | 1. Show that Maxwell's equations in vacuum give waves at c = 1/√(μ₀ε₀) with B = E/c, compute the Poynting intensity of a wave and the inverse-square fall of a point source, and compute the radiation pressure on a solar sail. 2. State the far-field patterns of the short and the half-wave dipole, compute their directivities (1.5 and 1.64) and beamwidths by integrating over the sphere with a `scipy.integrate.quad` oracle, and show numerically that the half-wave pattern is the sum of short-dipole elements. 3. Relate gain, efficiency, effective aperture and polarisation, and derive the Friis transmission equation from the inverse-square law and the effective aperture. 4. Build and defend a link budget for a 437 MHz CubeSat downlink from 500 km: slant range, path loss, received power, noise, Eb/N0 and margin as functions of elevation, and decide the highest data rate the link supports. |
| **Prerequisites** | `P04.06` (waves: wavelength, frequency, the wave equation, energy carried by a wave), `P04.09` (Biot–Savart, Faraday's law), `P04.10` (phasors, impedance, decibels, the matched load), `P02.11` (the wave equation as a PDE), `P02.07` (phasors), `P01.07` (vectors, cross product), `P03.04` (quadrature). Skill check: you can convert a power ratio to decibels and back, and integrate a function of θ over a sphere with the sin θ weight. |
| **Engineering question** | A CubeSat at 500 km transmits 1 W at 437 MHz from a half-wave dipole to a ground station with a 14 dBi Yagi antenna. At 10° elevation (slant range 1695 km), what power arrives, what is the signal-to-noise ratio per bit at 9600 bit/s, and what data rate could the pass support with a 3 dB margin? Answer: −109.7 dBm, Eb/N0 = 22.1 dB, a 12.1 dB margin at 9600 bit/s, and about 78 kbit/s at 3 dB, under the stated assumptions on losses and noise temperature. |
| **Assumptions and limits** | Vacuum propagation (no ionospheric absorption or scintillation, which matter at low elevation and below about 1 GHz); far field (distance ≫ λ and ≫ D²/λ); thin, centre-fed dipoles with sinusoidal current, in free space (no spacecraft body, no ground reflection); ideal matched receivers; losses (polarisation, cable, pointing, atmosphere) and the system noise temperature are budget assumptions stated in the code, not computed; Eb/N0 required for the modulation is an assumption. The patterns are exact for the model; a dipole on a CubeSat has a pattern distorted by the body and its gain toward the ground station may be well below 2.15 dBi. |
| **Analytical work** | The wave equation and c from Maxwell's equations; E/B = c; the Poynting vector and I = E₀²/(2η₀); inverse square from energy conservation; radiation pressure I/c and 2I/c; the Hertzian far field and its sin²θ pattern; directivity of the short dipole (1.5) by hand; the half-wave pattern and its directivity by numerical integration; radiation resistance; A_e = Gλ²/(4π); the Friis equation; polarisation loss factors; slant range from the elevation angle; Eb/N0 and margin in decibels. |
| **Computational deliverable** | `starter/main.py`: `wave_speed`, `b_from_e`, `plane_wave_intensity`, `e_from_intensity`, `isotropic_intensity`, `radiation_pressure`, `sail_acceleration`, `hertzian_pattern`, `dipole_pattern`, `halfwave_pattern`, `dipole_pattern_numeric`, `directivity`, `directivity_quad`, `half_power_beamwidth`, `radiation_resistance_short`, `radiation_resistance_dipole`, `to_db`, `from_db`, `gain`, `effective_aperture`, `dish_gain`, `polarization_loss`, `free_space_path_loss`, `friis_received_power`, `slant_range`, `hertzian_far_field`, `hertzian_radiated_power`, `power_through_sphere` (the `Link` budget class and `main` are given). Tests prove: 1/√(μ₀ε₀) = c and c/n in a medium; B = E/c and I = E₀B₀/(2μ₀); the inverse square; the power through spheres of 5 m to 1000 km around a Hertzian dipole equals its radiated power to 10⁻⁶; pressure doubles on a mirror and the sail scales as 1/r², cos² and A/m; patterns symmetric, null on the axis, maximum broadside; the short-dipole limit of any dipole is sin²θ; directivities 1.5 and 1.641 against an independent `quad`; beamwidths 90° and 78.1°; the numerical element sum equals the closed form and converges at order 2; R_rad = 73.08 Ω and D = η₀/(πR); A_e = ηπD²/4 for a dish; polarisation factors; Friis reciprocity, 20 dB per decade, and Friis = intensity × aperture; slant range at zenith and horizon and by the law of cosines; the link budget reproduces Friis and kTB and loses 3 dB per doubling of bit rate. |
| **Practical activity** | **computational experiment**: dipole patterns and directivities for lengths from 0.1λ to 2λ, a convergence study of the element sum, the inverse-square check by integrating the Poynting flux, a solar-sail acceleration study, and a link budget swept over elevation, frequency band, antenna and data rate. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. Nothing in this lesson involves transmitting; transmitting on any frequency requires the appropriate licence. |
| **Evidence** | Completed `main.py`, passing tests, `dipole_pattern.png`, the derivations of c, the Poynting intensity, the short-dipole directivity and the Friis equation, the four experiment tables, and a one-page link budget for a stated mission with every assumption listed and the margin at the lowest usable elevation. |
| **Assessment** | `assessments/phase-04/grading-rubric.md`; competencies 1–4. Pass: the Friis equation is derived from the inverse square and the effective aperture; all public tests pass; the link budget lists each term with its source (computed, assumed or datasheet) and states the data rate decision with its margin. |
| **References** | `@openstax-university-physics-vol2` Ch. 16 "Electromagnetic Waves": 16.1 Maxwell's Equations and Electromagnetic Waves, 16.2 Plane Electromagnetic Waves, 16.3 Energy Carried by Electromagnetic Waves, 16.4 Momentum and Radiation Pressure (required-free); 16.5 The Electromagnetic Spectrum (supplementary-free: band names). `@ellingson-electromagnetics-vol2` Ch. 9 "Radiation": 9.4 Radiation from a Hertzian Dipole, 9.5 Radiation from an Electrically-Short Dipole, 9.7 Far-Field Radiation from a Half-Wave Dipole; Ch. 10 "Antennas": 10.1 How Antennas Radiate, 10.2 Power Radiated by an Electrically-Short Dipole, 10.7 Directivity and Gain, 10.8 Radiation Pattern, 10.13 Effective Aperture, 10.14 Friis Transmission Equation (required-free: use the results; the vector-potential derivations are beyond this lesson); 10.5 Equivalent Circuit Model for Transmission; Radiation Efficiency, 10.10 Reciprocity (supplementary-free). `@nasa-basics-of-space-flight` Chapter 10 "Telecom": "Transmitters and Receivers", "Signal Power", "Uplink and Downlink"; Chapter 11 "Onboard Systems", "Telecommunications Subsystems" ("High-Gain Antennas (HGA)", "Low-Gain Antennas") (supplementary-free: how deep-space missions trade antenna gain against pointing). |
| **Extensions** | ★ Add a ground-reflected ray (a perfectly conducting ground under the Yagi at height h) and plot the received power against elevation: find the nulls. ★ Compute the pattern of a two-element array of half-wave dipoles spaced λ/4 and fed 90° apart (an endfire array) by adding the element sum of `dipole_pattern_numeric` with phase shifts, and find its front-to-back ratio. ★ Redo the link budget for a lunar distance (384 000 km) at 8.4 GHz and size the ground dish for 1 Mbit/s. |
| **Time** | unmeasured |

## 1. The engineering problem

A satellite is useful only if its data reach the ground. For a CubeSat the radio is a one-watt transmitter and an antenna that is often no more than two strips of tape-measure steel unrolled after deployment: a dipole. The ground station is a directional antenna on a rotator. Between them are hundreds to thousands of kilometres of vacuum, and the question every mission must answer before launch is whether the signal that arrives is far enough above the receiver's noise to be decoded, at what data rate, and over which part of each pass. That calculation, the **link budget**, is a sum of decibels, but each term is a piece of physics: the wave that carries the energy, the inverse-square law that dilutes it, the pattern of the transmitting antenna that concentrates it, the effective area of the receiving antenna that collects it, and the thermal noise that competes with it.

This lesson builds those terms from Maxwell's equations. Part 2 derives the electromagnetic wave and its intensity, uses the momentum it carries to size a solar sail, then derives the radiation pattern of a dipole, the directivity that measures how much it concentrates power, the effective aperture of a receiving antenna, and the Friis equation that combines them. Part 3 computes patterns and directivities numerically and checks them against closed forms and `scipy.integrate.quad`, then evaluates the link. Part 4 sweeps the design: elevation, frequency, antennas, data rate.

Success criteria, written first: the directivities 1.5 (short dipole) and 1.641 (half-wave) are reproduced to 10⁻⁶ by two independent integrations; the half-wave pattern computed as a sum of short-dipole elements agrees with the closed form and converges at second order; the Poynting flux through every sphere around a dipole equals its radiated power; the link budget reproduces the Friis equation and kTB exactly and gives a positive margin at 10° elevation for 9600 bit/s.

## 2. Understand and derive

### 2.1 From Maxwell's equations to waves

In vacuum, with no charges or currents, Maxwell's equations are ∇·E = 0, ∇·B = 0, ∇ × E = −∂B/∂t (Faraday, `P04.09`) and ∇ × B = μ₀ε₀ ∂E/∂t (Ampère with Maxwell's displacement current). Take the curl of the third: ∇ × (∇ × E) = ∇(∇·E) − ∇²E = −∇²E, and the right side becomes −∂(∇ × B)/∂t = −μ₀ε₀ ∂²E/∂t². So

**∇²E = μ₀ε₀ ∂²E/∂t²,**   (1)

the wave equation of `P02.11` and `P04.06` with speed **c = 1/√(μ₀ε₀)** = 2.997925 × 10⁸ m/s (OpenStax 16.1). Two constants measured with charges and currents on a bench predicted the speed of light; that is how light was recognised as an electromagnetic wave. In a medium of permittivity ε = n²ε₀ the speed is c/n.

A **plane wave** travelling along z has E along x and B along y, perpendicular to each other and to the direction of travel, in phase, with **B = E/c** (16.2): Faraday's law applied to E_x = E₀ cos(kz − ωt) gives ∂B_y/∂t = −∂E_x/∂z, so B_y = (k/ω)E_x = E_x/c. The direction of E is the wave's **polarisation**.

### 2.2 Energy, intensity, and the inverse-square law

The energy flux of the field is the **Poynting vector** S = E × B/μ₀ (W/m²), along the direction of travel. Time-averaged over a cycle of a plane wave,

**I = ⟨S⟩ = E₀B₀/(2μ₀) = E₀²/(2η₀),  η₀ = μ₀c = 376.73 Ω,**   (2)

the **impedance of free space** (16.3). Sunlight at 1 au, I = 1361 W/m², has E₀ = 1013 V/m and B₀ = 3.38 μT. A source radiating P equally in all directions spreads it over spheres of area 4πr²; energy conservation then requires

**I(r) = P/(4πr²),**   (3)

the **inverse-square law**: 1 W at 500 km gives 3.18 × 10⁻¹³ W/m², at 1000 km a quarter of that. The code checks (3) for a real (non-isotropic) source by integrating the Poynting flux of a short dipole's field over spheres of 10 m, 1 km and 100 km: 0.986278 W through each, equal to the dipole's radiated power.

### 2.3 Momentum: radiation pressure and solar sails

A wave carrying energy U also carries momentum U/c (16.4). Absorbed by a surface facing it, light of intensity I therefore exerts a pressure I/c; reflected straight back, it delivers twice the momentum, 2I/c. At 1 au this is 4.54 μPa absorbed and 9.08 μPa reflected. A flat mirror sail of area A at angle α between its normal and the sunline catches a beam of area A cos α and turns it about the normal, so the force is along the normal with magnitude **F = 2IA cos²α/c**. A 32 m² sail on a 5 kg spacecraft facing the Sun at 1 au: 291 μN, a = 5.81 × 10⁻⁵ m/s², 5.0 m/s of Δv per day without propellant; tilted 35° to push sideways (the angle that maximises the transverse component), 3.90 × 10⁻⁵ m/s²; at Mars's distance (1.52 au) 2.52 × 10⁻⁵ m/s², because I falls as 1/r². Tiny, and unlimited.

### 2.4 Radiation from a dipole

Charges radiate when they accelerate. The simplest antenna is a short wire of length Δl ≪ λ carrying an oscillating current I₀ cos ωt, the **Hertzian dipole**. Far from it (r ≫ λ), the field is a spherical wave whose electric field lies along the θ direction (θ measured from the wire's axis), with peak amplitude

**|E_θ| = η₀kI₀Δl sin θ/(4πr),**   (4)

k = 2π/λ (Ellingson 9.4; the derivation goes through the vector potential and is not repeated here). Three things are visible: the field falls as 1/r, so the intensity falls as 1/r² (2.2); it is proportional to the current-length product; and it has a **pattern**, sin θ, zero along the wire and largest broadside. The radiation intensity (power per steradian) is U(θ) ∝ sin²θ, a doughnut. Integrating |E|²/(2η₀) over a sphere gives the radiated power P = ½I₀²R_rad with **radiation resistance** R_rad = (2π/3)η₀(Δl/λ)² ≈ 80π²(Δl/λ)² Ω (10.2): a uniform-current dipole 0.05λ long has 1.97 Ω, which is why short antennas are hard to feed efficiently (a few ohms of radiation resistance against a fraction of an ohm of wire loss).

A real antenna is a longer wire with a current that varies along it. For a centre-fed thin dipole of total length L the current is close to sinusoidal, I(z) = I₀ sin(k(L/2 − |z|)), zero at the ends. Each element dz radiates like (4) with a phase e^{jkz cos θ} set by its position, and the far field is the sum. Doing the integral,

**U(θ) ∝ F(θ) = [(cos(kL/2 · cos θ) − cos(kL/2))/sin θ]²,**   (5)

and for the **half-wave dipole**, kL/2 = π/2, **F = [cos((π/2) cos θ)/sin θ]²** (9.7). The code computes the same pattern by summing 200 short elements numerically (`dipole_pattern_numeric`): the maximum difference from (5) is 2.1 × 10⁻⁵, falling as the square of the element length. Its radiation resistance referred to the feed current is R = (η₀/2π)∫F sin θ dθ = **73.08 Ω**, close to the 50 or 75 Ω of standard cables, which is one reason the half-wave dipole is everywhere.

### 2.5 Directivity, gain and beamwidth

The **directivity** compares an antenna's peak intensity with the intensity an isotropic antenna radiating the same total power would give:

**D = 4πU_max/P_rad,  P_rad = 2π∫₀^π U(θ) sin θ dθ**   (6)

for a pattern symmetric about the axis (Ellingson 10.7). For the short dipole by hand: ∫sin³θ dθ = 4/3, so P = 8π/3 U_max and **D = 1.5 (1.76 dBi)**. For the half-wave dipole the integral has no elementary form; numerically **D = 1.641 (2.15 dBi)**, by the trapezoid rule on 2001 points and by `quad`, agreeing to 10⁻⁹. The **half-power beamwidth** is the angle between the directions where U falls to half its maximum: 90.0° for the short dipole and 78.1° for the half-wave. Longer dipoles concentrate more: 1λ gives D = 2.41 (3.82 dBi), 1.25λ gives 3.28 (5.16 dBi, the most a single dipole achieves); at 1.5λ the broadside lobe splits and the maximum moves to 42.6° from the axis.

**Gain** is directivity times radiation efficiency, G = e_rad D: the fraction of the power delivered to the antenna that is radiated rather than dissipated in its own resistance (10.5). A half-wave dipole of good conductor has e_rad > 0.95; an electrically short antenna on a small spacecraft can have e_rad of 0.5 or less. Gains are quoted in **dBi**, decibels relative to an isotropic radiator.

A **dish** of diameter D at wavelength λ has gain η(πD/λ)², with an aperture efficiency η ≈ 0.55: a 3 m dish at 2.25 GHz (S band) has 34.4 dBi, and a beamwidth of about λ/D radians, a few degrees, which must be pointed.

### 2.6 Receiving: effective aperture and the Friis equation

A receiving antenna in a wave of intensity I delivers to a matched load the power P = IA_e, with A_e its **effective aperture**. By reciprocity (an antenna's pattern is the same transmitting and receiving, 10.10) and a thermodynamic argument (10.13),

**A_e = Gλ²/(4π).**   (7)

An isotropic antenna collects λ²/(4π) whatever its physical size: 0.0375 m² at 437 MHz, 0.0014 m² at 2.25 GHz. A half-wave dipole at 437 MHz has 0.0615 m²; a 3 m dish at 2.25 GHz has 3.89 m². High frequencies need large dishes to collect the same power, which is why the "path loss" below grows with frequency even in vacuum: the loss is a property of the receiving antenna, not of space.

Combine (3) with the transmitter's gain (its intensity in the receiver's direction is G_t times the isotropic value) and (7):

**P_r = P_t G_t G_r (λ/(4πd))²**   (8)

the **Friis transmission equation** (10.14). In decibels, P_r = P_t + G_t + G_r − FSPL, with the free-space path loss FSPL = 20 log₁₀(4πd/λ): 20 dB more for every factor 10 in distance or frequency. It assumes the far field, matched impedances and **matched polarisation**.

### 2.7 Polarisation

The receiving antenna collects only the part of the wave's field along its own polarisation. Two linear antennas at angle ψ couple a fraction **cos²ψ** of the power; crossed (ψ = 90°) they couple none. A tumbling or rotating spacecraft with a linear dipole presents an unknown, changing ψ, so ground stations use **circular polarisation** (a field that rotates once per cycle, right- or left-handed): a linear wave received on a circular antenna loses exactly half the power, **3.01 dB**, at any orientation, a fixed and budgetable loss instead of fades to zero. Two circular antennas of the same hand couple fully; of opposite hands, not at all. The budget below carries 3 dB for the linear CubeSat dipole into a circular ground antenna.

### 2.8 The link budget

Geometry first. A satellite at altitude h seen at elevation e over a spherical Earth of radius R is at slant range

**d = √((R + h)² − R² cos²e) − R sin e**   (9)

(the law of cosines in the triangle of Earth's centre, station and satellite): 500 km at the zenith, 1695 km at 10°, 2077 km at 5°.

Then noise. A receiver at system noise temperature T_sys sees thermal noise of spectral density N₀ = kT_sys; at UHF, with sky, ground and receiver contributions, T_sys ≈ 500 K is a reasonable assumption, N₀ = −171.6 dBm/Hz. At bit rate R_b the energy per bit is E_b = P_r/R_b, and the figure of merit is **E_b/N₀ = P_r/(kT_sys R_b)**. A modulation and coding scheme needs a certain E_b/N₀ (assumed here 10 dB) for an acceptable bit error rate; the excess is the **margin**.

The budget at 10° elevation (`Link`, all assumptions in its fields):

| term | value | source |
|---|---|---|
| transmitter power | 1 W = 0 dBW | assumption (typical CubeSat UHF radio) |
| transmit antenna | +2.15 dBi | half-wave dipole, peak (2.5) |
| receive antenna | +14.0 dBi | ground Yagi, datasheet-class assumption |
| slant range | 1695 km | (9) |
| free-space path loss | −149.8 dB | (8) at 437 MHz |
| polarisation, cable, pointing and atmosphere | −3.0 − 1.0 − 2.0 dB | (2.7) and assumptions |
| received power | −109.7 dBm | sum |
| N₀ at 500 K | −171.6 dBm/Hz | kT |
| E_b/N₀ at 9600 bit/s | 22.1 dB | −109.7 + 171.6 − 39.8 |
| required E_b/N₀ | 10.0 dB | assumption for the modulation |
| **margin** | **12.1 dB** | |

The margin is 10.3 dB at 5°, 17.5 dB at 30° and 22.7 dB at the zenith. With 3 dB kept in reserve, the link at 10° supports 9600 × 10^{(12.1 − 3)/10} ≈ **78 kbit/s**. NASA's *Basics of Space Flight* (Chapter 10, "Signal Power") describes the same budget for deep-space links, where the distances are a million times larger and the 60 dB more path loss is paid for by 70 m dishes and very low data rates.

### 2.9 Where the models stop

The dipole's gain is 2.15 dBi only broadside: toward its axis it has a null, and a CubeSat with a dipole pointing at the ground station delivers nothing; that is why CubeSats use crossed dipoles or turnstiles, and why the "pointing" loss above is a placeholder, not a calculation. The far-field formulas need r ≫ λ and r ≫ 2D²/λ (for a 3 m dish at S band, about 140 m: irrelevant for orbit, relevant on an antenna test range). Below about 1 GHz the ionosphere rotates the polarisation (Faraday rotation) and scintillates at low elevation; at low elevation, ground reflections add and cancel. T_sys and the required E_b/N₀ are the least certain inputs and must come from the receiver and the modem.

## 3. Build the computational model

Open `starter/main.py`. The `Link` dataclass and `main` are given; `Link.budget` calls your functions.

### Step 1: waves and pressure

`wave_speed`, `b_from_e`, `plane_wave_intensity` (2), `e_from_intensity`, `isotropic_intensity` (3), `radiation_pressure`, `sail_acceleration` (2.3).

### Step 2: patterns

`hertzian_pattern` is sin²θ; `dipole_pattern(theta, length_wl)` is (5), returning 0 where sin θ = 0 (use `np.where`); `halfwave_pattern` calls it with 0.5. `dipole_pattern_numeric` sums the elements:

```python
z = -L/2 + dz * (np.arange(n_seg) + 0.5)
I = np.sin(k * (L/2 - np.abs(z)))
af = np.sum(I * np.exp(1j * k * np.outer(np.cos(theta), z)), axis=1) * dz
return np.abs(np.sin(theta) * af * k / 2) ** 2
```

with k = 2π (lengths in wavelengths); the factor k/2 matches the normalisation of (5).

### Step 3: integrals over the sphere

`directivity(pattern)` applies (6) with `np.trapezoid` on θ ∈ [0, π]; `directivity_quad` with `quad`; `half_power_beamwidth` finds the two half-power angles on a fine grid; `radiation_resistance_short` and `radiation_resistance_dipole` (by `quad`).

### Step 4: gains and the link

`to_db`, `from_db`, `gain`, `effective_aperture` (7), `dish_gain`, `polarization_loss` (2.7), `free_space_path_loss`, `friis_received_power` (8), `slant_range` (9). Then `hertzian_far_field` (4), `hertzian_radiated_power` and `power_through_sphere`, the inverse-square check.

Running `main()` prints

```
1/sqrt(mu0 eps0) = 2.997925e+08 m/s (c = 2.997925e+08); eta0 = 376.73 ohm
sunlight at 1 au, 1361 W/m^2: E0 = 1013 V/m, B0 = 3.38 uT; pressure absorbed 4.54 uPa, reflected 9.08 uPa
sail 32 m^2, 5 kg, perfect mirror at 1 au: force 291 uN, a = 5.81e-05 m/s^2, dv per day 5.02 m/s; at 35 deg cone 3.90e-05 m/s^2; at 1.52 au 2.52e-05 m/s^2
1 W isotropic at 500 km: 3.18e-13 W/m^2, at 1000 km: 7.96e-14 W/m^2

short dipole: D = 1.50000 (quad 1.50000, exact 1.5) = 1.76 dBi, HPBW 90.00 deg
half-wave dipole: D = 1.64092 (quad 1.64092) = 2.15 dBi, HPBW 78.08 deg, R_rad = 73.08 ohm
half-wave pattern by summing 200 Hertzian elements: max |difference| 2.1e-05
  L = 0.10 lambda: D = 1.505 (1.78 dBi), max at 90.0 deg
  L = 0.50 lambda: D = 1.641 (2.15 dBi), max at 90.0 deg
  L = 1.00 lambda: D = 2.411 (3.82 dBi), max at 90.0 deg
  L = 1.25 lambda: D = 3.282 (5.16 dBi), max at 90.0 deg
  L = 1.50 lambda: D = 2.226 (3.48 dBi), max at 42.6 deg
```

then the sphere-flux check, the apertures, and the link:

```
CubeSat downlink 437 MHz, 1.0 W, 2.15 dBi -> 14.0 dBi, 500 km at 10 deg elevation: slant range 1695 km
  FSPL 149.8 dB, other losses 6.0 dB, P_r = -109.7 dBm, N0 = -171.6 dBm/Hz, Eb/N0 = 22.1 dB at 9600 bit/s, margin 12.1 dB
  elevation    5 deg: range   2077 km, FSPL 151.6 dB, margin  10.3 dB
  elevation   10 deg: range   1695 km, FSPL 149.8 dB, margin  12.1 dB
  elevation   30 deg: range    909 km, FSPL 144.4 dB, margin  17.5 dB
  elevation   90 deg: range    500 km, FSPL 139.2 dB, margin  22.7 dB
  polarisation: linear-linear at 90 deg 0.00, linear-circular -3.01 dB, RHCP-LHCP 0
  highest bit rate with margin >= 3 dB at 10 deg: 77999 bit/s
```

and writes `dipole_pattern.png` (polar patterns of the short, half-wave and 1.25λ dipoles; the half-wave pattern in dB, closed form against the element sum; the margin against elevation).

### Verify

Run `python -m pytest curriculum/phase-04/lesson-11-electromagnetic-waves-antennas`. The tests:

1. **c**: 1/√(μ₀ε₀) = c to 10⁻⁹; c/2 for ε = 4ε₀.
2. **Plane wave**: B = E/c; I = E₀B₀/(2μ₀) = cε₀E₀²/2; inversion; I ∝ E².
3. **Inverse square**: I·4πr² = P; ratios 4 and 100.
4. **Sphere flux**: the Poynting flux of (4) through spheres of 5 m, 1 km and 1000 km equals ½I₀²R_rad to 10⁻⁶; |E| ∝ 1/r and ∝ sin θ; P ∝ I₀².
5. **Pressure and sail**: I/c and 2I/c; a ∝ 1/r², cos²α, A/m.
6. **Patterns**: symmetric about 90°, null on the axis, unity broadside, monotone toward broadside; half-wave narrower; a 0.01λ dipole normalised is sin²θ to 10⁻³.
7. **Directivity**: 1.5 to 10⁻⁶; the half-wave value against an independent `quad` integral to 10⁻⁶ (trapezoid) and 10⁻⁹ (`directivity_quad`); 2.15 dBi.
8. **Beamwidth**: 90.0° and 78.1°.
9. **Element sum**: equals (5) for L = 0.5, 1 and 1.25λ to 10⁻⁴ of the peak; observed order 2.0 ± 0.2 in the element length.
10. **R_rad**: 73.08 Ω; the short-dipole formula and its L² scaling; D = η₀/(πR_rad) for the half-wave dipole.
11. **Aperture**: λ²/(4π) isotropic; A_e of a dish equals η × area; dish gain ∝ 1/λ²; dB round trip.
12. **Polarisation**: cos²ψ; 0.5 linear to circular at any angle; same hand 1, opposite 0; invalid input raises.
13. **Friis**: reciprocity in G_t ↔ G_r; −20 dB per decade of distance; P_r = I_iso G_t × A_e; FSPL in dB.
14. **Slant range**: h at zenith, the horizon distance at 0°, decreasing with elevation, the law of cosines.
15. **Link**: the budget equals Friis with losses and kTB to 10⁻⁹ dB; −3.01 dB per doubling of bit rate; zenith gain 20 log₁₀(d/h); a margin between 0 and 30 dB.
16. **main** writes the figure.

### Compare

`scipy.integrate.quad` is the oracle for every integral over the sphere: it reproduces the trapezoid-rule directivities to 10⁻⁹ and gives R_rad = 73.08 Ω. The deeper check is physical: the half-wave pattern, a closed form from a textbook integral, equals the numerical superposition of 200 short dipoles of (4) to 2 × 10⁻⁵, with the second-order convergence of the midpoint rule; and the Poynting flux of (4) is the same through spheres four decades apart in radius, which is the inverse-square law as a conservation law rather than an assumption. Professional antenna codes (method-of-moments solvers) compute the current on the wire instead of assuming it is sinusoidal; for a thin half-wave dipole the assumption is good to a few per cent in pattern and directivity, which is far below the uncertainty of the budget's loss terms.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction first.

1. **Dipole length.** Directivity, beamwidth, direction of maximum, and the resistance seen at the feed, R_feed = R_rad/sin²(kL/2) with R_rad from `radiation_resistance_dipole` (referred to the current amplitude I₀), for L = 0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5 and 2λ. Predict D rising from 1.5 to a maximum near 1.25λ and the pattern splitting beyond it.
2. **Element-sum convergence.** Maximum error of `dipole_pattern_numeric` against (5) for 25, 50, 100, 200, 400 elements at L = 0.5λ and 1.25λ; predict order 2.
3. **Solar sail.** Acceleration of a 32 m², 5 kg sail against heliocentric distance (0.5–5 au) and cone angle (0–80°), and the ratio to the Sun's gravity at each distance. Predict that the ratio is independent of distance (both fall as 1/r²).
4. **Link sweep.** Margin against elevation for (a) the baseline, (b) a crossed-dipole CubeSat antenna that is circularly polarised (remove the 3 dB polarisation loss, keep the other 3 dB, assume 0 dBi average gain), (c) 2.25 GHz with a 2 W transmitter, a 6 dBi circularly polarised patch and a 1 m dish (55 % efficiency) at T_sys = 300 K with 4 dB of other losses; for each, the highest bit rate with 3 dB margin at 10°. Predict that (c) supports more data despite 14 dB more path loss.

Raw data format: CSV tables `L_lambda, D, D_dBi, HPBW_deg, theta_max_deg, R_rad_ohm` (1); `n_seg, L_lambda, max_err` (2); `r_au, cone_deg, a_mps2, a_over_gsun` (3); `case, elevation_deg, slant_km, fspl_dB, Pr_dBm, EbN0_dB, margin_dB, max_rate_bps` (4); plus `dipole_pattern.png`.

What you should find: D = 1.505, 1.532, 1.641, 1.882, 2.411, 3.282, 2.226, 2.529 for the eight lengths, the maximum leaving broadside at 1.5λ (42.6°) and 2λ (57.4°); feed resistances 2.0, 13.4, 73.1, 371, ∞, 213, 105 and ∞ Ω (infinite where the centre feed sits at a current node, L = λ and 2λ: such a dipole must be fed elsewhere); the element-sum error falling at order 2.00 for 0.5λ from 25 elements on (1.3 × 10⁻³ to 5.1 × 10⁻⁶), and at 1.25λ only from 100 elements on; a sail-to-gravity ratio of 0.98 % at every distance for the 32 m², 5 kg sail facing the Sun; for (b) a margin of 12.9 dB at 10° and about 95 kbit/s, better than the baseline although the antenna has less peak gain, because it has no null and no polarisation loss; for (c) a dish of 24.9 dBi and a margin of −0.4 dB at 1 Mbit/s, so the highest rate with 3 dB margin is about 460 kbit/s, six times the UHF link's 78 kbit/s.

## 5. Validate and improve

| quantity | theory (§2) | code | difference | cause |
|---|---|---|---|---|
| c from μ₀, ε₀ | 2.997925 × 10⁸ m/s | same | < 10⁻⁹ | — |
| E₀ of sunlight | √(2η₀I) = 1013 V/m | 1013 V/m | 0 | — |
| short-dipole D | 1.5 | 1.50000 | < 10⁻⁶ | — |
| half-wave D | `quad`: 1.64092 | trapezoid 1.64092 | < 10⁻⁹ | — |
| half-wave pattern | (5) | element sum | 2.1 × 10⁻⁵ | 200 elements, O(h²) |
| R_rad, half-wave | 73.08 Ω | 73.08 Ω | — | — |
| D = η₀/(πR) | 1.641 | 1.641 | < 10⁻⁴ | — |
| sphere flux, r = 10 m to 100 km | ½I₀²R_rad = 0.986278 W | 0.986278 W | < 10⁻⁶ | — |
| Friis vs budget | equal | equal to 10⁻⁹ dB | — | — |
| slant range at zenith | h = 500 km | 500 km | 0 | — |

Every computed quantity with an independent value agrees with it. That validates the physics and the code; it does not validate the budget, whose largest uncertainties are the assumed terms (T_sys ± 2 dB, the satellite antenna's gain toward the station, which can be anywhere from +2 dB to a null, pointing and polarisation). The engineering answer is a **12 dB margin at 9600 bit/s from 10° elevation**, robust to any single assumption being wrong by a few decibels, and a usable rate of about **78 kbit/s with 3 dB in reserve**; a mission that needs more data per pass should move to S band with a dish (Part 4, item 4c), not to more transmit power, since 1 W more buys 3 dB and a 1 m dish buys 11.

What to improve: a measured or simulated pattern of the satellite's antenna on its body, replacing the 2.15 dBi peak; a pass simulation that integrates the data volume over the elevation profile; a T_sys from the actual receiver. What not to improve: the directivities beyond three digits.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (1), (2), the short-dipole directivity and (8); `dipole_pattern.png`; the four experiment tables; and a one-page link budget (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) for a mission of your choice, every term marked computed, assumed or from a datasheet, the margin at the lowest elevation you intend to use, and the data-rate decision.

### Adversarial self-checks

- **What if the distance doubles?** P_r falls by 6.02 dB; at 10° the margin becomes 6.1 dB at 9600 bit/s, still closing; the usable rate halves twice (÷ 4). From 500 km to geostationary distance (about 36 000 km slant) the loss grows by 26.5 dB more than at 1695 km: a CubeSat UHF dipole link from GEO does not close.
- **What if the satellite's dipole points at the station?** Its gain in that direction is zero (the null of sin²θ); the link fails for as long as the geometry lasts. A crossed-dipole or turnstile antenna, circularly polarised, fills the null and removes the 3 dB polarisation loss.
- **Which assumption fails first in the real world?** The satellite antenna gain toward the ground (body blockage, attitude), then T_sys (urban noise at UHF can add several decibels), then the low-elevation atmosphere and multipath.
- **What measurement could falsify the budget?** The received signal strength logged by the ground station over a pass, against elevation: the model predicts a 12.4 dB rise from 5° to the zenith (FSPL 151.6 → 139.2 dB, the rest fixed). A flat or erratic record means attitude or pattern effects dominate.
- **Does the result survive a change of method?** The trapezoid rule and `quad` agree on D to 10⁻⁹; the closed-form pattern and the element sum to 10⁻⁵; Friis and the dB budget to 10⁻⁹ dB. The physics is not in doubt; the assumptions are.
- **Which decision changes if the uncertainty doubles?** If the combined uncertainty of the assumed terms grows from ±3 to ±6 dB, the 3 dB reserve is no longer enough to promise 78 kbit/s; a conservative mission would plan about 19 kbit/s (6 dB below the 78 kbit/s point) and adapt the rate during the pass.

## Check yourself

1. Which two constants give the speed of light, and what is E/B in a plane wave?
2. The Sun delivers 1361 W/m² at 1 au. What does it deliver at 2 au, and what pressure does it exert on a perfect mirror there?
3. Why does a short dipole radiate nothing along its axis, and what is its directivity in dBi?
4. A half-wave dipole has D = 1.64. What does that number mean, and what is its radiation resistance?
5. What is the effective aperture of an isotropic antenna at 437 MHz, and why is it smaller at 2.25 GHz?
6. Write the Friis equation in decibels. How many decibels does doubling the distance cost?
7. A linearly polarised signal is received by a circularly polarised antenna. What fraction of the power is lost, and why do ground stations accept this?
8. A link has a 12 dB margin at 9600 bit/s. Roughly what bit rate gives a 3 dB margin, and why?
