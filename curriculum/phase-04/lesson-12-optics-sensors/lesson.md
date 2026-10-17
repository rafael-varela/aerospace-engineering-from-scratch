# P04.12  Optics & Sensors

> A star is a point; a good star tracker deliberately blurs it over a few pixels, because a point on one pixel says only which pixel, while a spot over nine says where inside it to a twentieth.

| | |
|---|---|
| **Learning outcomes** | 1. Apply the laws of reflection and refraction, find the critical angle for total internal reflection and the normal-incidence reflection loss, and use the thin-lens equation and the f-number to describe a camera. 2. Size a star tracker: pixel scale, field of view, diffraction-limited spot against pixel size, expected number of stars in the field and photo-electrons per star, from first principles and stated assumptions. 3. Estimate the centroiding accuracy of a star image with a seeded Monte Carlo on SYNTHETIC images, show the photon-limited 1/√N law, the cost of read noise and window size, and the pixel-locking bias of an over-sharp spot, and turn the result into an attitude accuracy. 4. Model four sensor archetypes (the pitot-static probe, the strain gauge in a Wheatstone bridge, the NTC thermistor, the spring-mass accelerometer) and express an ADC's resolution in the sensor's physical unit. |
| **Prerequisites** | `P04.10` (the ADC, LSB, op-amp gain, dividers), `P04.11` (light as an electromagnetic wave, intensity), `P04.08` (resistance, dividers, Kirchhoff), `P04.06` (waves; diffraction is wave interference), `P02.16` (random variables, variance of a mean), `P00.14` (uncertainty), `P01.07` (vectors and angles). Skill check: you can propagate a variance through a sum and convert between radians and arcseconds. |
| **Engineering question** | A star tracker with a 25 mm f/1.2 lens and a 1024 × 1024 detector of 5.5 μm pixels must report the spacecraft's pointing to a few arcseconds. What are its pixel scale and field of view, how many stars of magnitude 6 or brighter does it see, should the stars be in sharp focus, how many photo-electrons does a magnitude-6 star give in 100 ms, and what centroid and attitude accuracy follow? Answer: 45.4″ per pixel, a 12.85° field with about 17 stars, a focused spot of 0.16 pixel that must be defocused to about one pixel, 543 electrons, a centroid of 0.068 pixel (3.1″) per star with a 5 × 5 window and 3 e⁻ read noise, and 0.75″ across the boresight from 17 stars, under the stated photometric assumptions. |
| **Assumptions and limits** | Geometric optics for imaging (paraxial thin lens, no aberrations) and scalar diffraction for the spot size; a Gaussian spot in place of the defocused, aberrated real spot; a uniform star density over the sky from an empirical whole-sky count (the galactic plane has several times the average, the poles fewer); an assumed photometric zero point for a broad silicon band, good to a factor of two; Poisson photon noise and Gaussian read noise only (no dark current, stray light, smear from rotation or detector non-uniformity). The sensor models are the textbook ideal versions: incompressible and subsonic isentropic pitot, ideal bridge, beta-model thermistor, undamped static accelerometer. |
| **Analytical work** | Snell's law, the critical angle and fibre acceptance; normal-incidence reflectance; the thin-lens equation and magnification; pixel scale p/f and field of view 2 arctan(np/2f); the Airy radius 1.22λN; the solid angle of a rectangular field; magnitudes and photon counts; the variance of a centroid σ²/N plus the pixel's 1/12; the attitude error over n stars; Bernoulli for the pitot and its compressible correction; the quarter-bridge output and its linearisation; the beta model and its inverse; x = a/ω_n²; ADC resolution and bits needed. |
| **Computational deliverable** | `starter/main.py`: `refraction_angle`, `critical_angle`, `fresnel_reflectance_normal`, `thin_lens_image_distance`, `magnification`, `f_number`, `pixel_scale`, `field_of_view`, `airy_radius`, `airy_gaussian_sigma`, `stars_brighter_than`, `fov_solid_angle`, `mean_stars_in_fov`, `star_electrons`, `gaussian_spot`, `centroid`, `simulate_centroid_errors`, `centroid_sigma_photon_limit`, `pixel_locking_bias`, `attitude_error_cross_boresight`, `pitot_airspeed`, `pitot_delta_p_compressible`, `calibrated_airspeed`, `bridge_output`, `quarter_bridge`, `full_bridge`, `thermistor_resistance`, `thermistor_temperature`, `divider_voltage`, `accelerometer_deflection`, `adc_resolution`, `bits_needed` (the `StarTracker` dataclass and `main` are given). Tests prove: Snell's law, reversibility, bending toward the normal, TIR just past 41.81°, 4 % reflectance; the thin-lens equation with its infinity, 2f and virtual-image limits; p/f, the arctan FOV and its small-angle limit, the solid angle's two limits; the Airy radius ∝ λN and below one pixel; star and photon counts scaling with magnitude, area, time and solid angle; the pixel-integrated spot conserving flux with an exact centroid; SYNTHETIC centroid errors falling with slope −0.50 ± 0.07 and matching √(σ² + 1/12)/√N within 12 %, unbiased and reproducible; read noise punishing large windows; pixel locking above 0.1 px at σ = 0.2 px and below 10⁻³ px at 1 px; the compressible pitot reducing to Bernoulli at low Mach and giving CAS = TAS at sea level; the bridge's balance, exact quarter-bridge output, linear limit and full-bridge factor 4; the thermistor round trip; x ∝ 1/f_n²; and bits for a stated resolution. |
| **Practical activity** | **computational experiment**: a star-tracker optics budget with a design sweep over focal length, f-number and pixel size, a SYNTHETIC-image Monte Carlo of centroiding against signal, read noise, window size and defocus, and a sensor-chain study expressing each archetype's ADC resolution in its own unit. |
| **Equipment** | None. A computer with the environment from `P00.01`. (The thermistor and divider of this lesson are the circuit measured in `P04.15`.) |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `star_tracker.png`, the derivations of the pixel scale and FOV, the centroid variance and the quarter-bridge output, the four experiment tables (each SYNTHETIC one labelled), and a one-page star-tracker budget for a stated requirement with every assumption listed. |
| **Assessment** | `assessments/phase-04/grading-rubric.md`; competencies 1–4. Pass: the optics budget is derived and reproduced in code; all public tests pass; the Monte Carlo is seeded, labelled SYNTHETIC and compared with the √(σ² + 1/12)/√N law; the budget states the defocus decision and the attitude accuracy with its assumptions. |
| **References** | `@openstax-university-physics-vol3` Ch. 1 "The Nature of Light": 1.1 The Propagation of Light, 1.2 The Law of Reflection, 1.3 Refraction, 1.4 Total Internal Reflection (required-free); 1.7 Polarization (supplementary-free). Ch. 2 "Geometric Optics and Image Formation": 2.4 Thin Lenses, 2.6 The Camera (required-free); 2.8 Microscopes and Telescopes (supplementary-free). `@openstax-university-physics-vol2` Ch. 9 "Current and Resistance": 9.3 Resistivity and Resistance (required-free: temperature dependence of resistance, the physics of the thermistor). `@leishman-intro-aerospace-flight-vehicles` "Airspeed Definitions & Measurement", https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/determination-of-airspeed/ (required-free: pitot-static measurement and position error). `@sparkfun-learn-tutorials` "Getting Started with Load Cells": "Strain Gauge Basics" (Gauge Factor; Small Changes in Strain; Amplifiers and Wheatstone Bridge); "Accelerometer Basics": "How an Accelerometer Works", "How to Select an Accelerometer" (Range); "Analog to Digital Conversion": "Relating ADC Value to Voltage" (required-free). `@fiore-opamps-linear-ics` Ch. 12 "Analog-to-Digital-to-Analog Conversion": 12.2 The Sampling Theorem, 12.3 Resolution and Sampling Rate, 12.5 Analog-to-Digital Conversion (required-free). `@nasa-basics-of-space-flight` Chapter 11 "Onboard Systems", "Attitude and Articulation Control Subsystems (AACS)" ("Celestial Reference", "Inertial Reference") (supplementary-free: where star trackers sit in an attitude system). |
| **Extensions** | ★ Replace the Gaussian spot by a defocused Airy pattern (a uniform disc convolved with the Airy core) and repeat the pixel-locking study. ★ Add rotation smear: the spacecraft turns at 0.5°/s during the 100 ms exposure; compute the streak length in pixels and the centroid error along and across it. ★ Estimate the roll accuracy about the boresight from the star positions at the edge of the field (it is worse than the cross-boresight accuracy by roughly the ratio of the field radius to the distance between stars, a factor of order ten). |
| **Time** | unmeasured |

## 1. The engineering problem

Every spacecraft must know where it is pointing. Gyroscopes measure rotation but drift; sun sensors and magnetometers give a direction to a few tenths of a degree; for anything better, an Earth-observation camera, a telescope, a laser link, a planetary approach, the reference is the stars. A **star tracker** is a small camera that photographs a patch of sky, finds the star images, measures their positions to a fraction of a pixel, identifies the pattern against a catalogue and computes the attitude. Commercial units give a few arcseconds (one arcsecond is 4.8 μrad: a coin at 4 km). How is a lens and a detector turned into an arcsecond instrument, when each pixel spans 45 arcseconds?

The answer is centroiding: if the star's light is spread over several pixels, the brightness-weighted mean position can be located to a small fraction of a pixel, limited by the number of photons and the detector noise. This lesson builds the star tracker's optics budget from the physics of light: reflection and refraction, the thin lens, the f-number, the field of view and the pixel scale, diffraction (which sets the smallest possible spot), and the photon counting that sets the noise. A SYNTHETIC-image Monte Carlo then measures centroid accuracy against the 1/√N law. The second half meets the sensors that fill the rest of an aircraft or spacecraft (the pitot-static probe, the strain gauge, the thermistor, the accelerometer) and ends where `P04.10` left off: what the ADC's least significant bit is worth in pascals, microstrain, kelvin or g.

Success criteria, written first: the thin-lens, Snell and critical-angle relations hold to rounding; the code reproduces 45.4″ per pixel and a 12.85° field for the baseline design; the SYNTHETIC centroid error falls as N^{−1/2} with the level √(σ² + 1/12)/√N within 12 %; pixel locking appears for a spot narrower than half a pixel and vanishes at one pixel; the compressible pitot formula reduces to Bernoulli at low Mach.

## 2. Understand and derive

### 2.1 Reflection, refraction and total internal reflection

Light travels in straight lines in a uniform medium at c/n, n the refractive index (OpenStax Vol. 3, 1.1). At a smooth surface it reflects at the angle of incidence (1.2) and the transmitted part refracts by **Snell's law**

**n₁ sin θ₁ = n₂ sin θ₂**   (1)

(1.3), bending toward the normal when it enters a denser medium: 30° in air becomes 19.47° in glass of n = 1.5. Going the other way, from n₁ > n₂, sin θ₂ reaches 1 at the **critical angle** θ_c = arcsin(n₂/n₁), and beyond it no light is transmitted: **total internal reflection** (1.4). For glass to air θ_c = 41.81°; for an optical fibre with a core of 1.48 and a cladding of 1.46 it is 80.57°, so rays within about 10° of the axis are trapped and guided round bends. Even below θ_c some light is reflected: at normal incidence the fraction is ((n₁ − n₂)/(n₁ + n₂))², **4 % per uncoated glass surface**, which is why every lens in a star tracker is anti-reflection coated and why stray reflections of the Sun ("ghosts") are a design driver.

### 2.2 The thin lens and the camera

A thin lens of focal length f images an object at distance d_o to a distance d_i with

**1/d_o + 1/d_i = 1/f,  m = −d_i/d_o**   (2)

(2.4). A star is at infinity, so its image forms exactly at d_i = f: a star tracker is focused at infinity and the detector sits in the focal plane. An object at 1 m images at 25.641 mm behind a 25 mm lens; an object at 2f images at 2f, inverted, life-size. The light-gathering power is set by the aperture diameter D, and the brightness of an extended image by the **f-number N = f/D** (2.6): image irradiance ∝ 1/N². Star trackers use fast lenses, f/1 to f/2, to collect photons in short exposures.

### 2.3 Pixel scale and field of view

A star at angle α from the optical axis images at x = f tan α in the focal plane. A detector of n pixels of pitch p therefore covers a **field of view**

**FOV = 2 arctan(np/(2f)),**   (3)

and near the axis one pixel subtends the **pixel scale** p/f radians. The baseline tracker: f = 25 mm, p = 5.5 μm, n = 1024: p/f = 220 μrad = **45.4″ per pixel**, **FOV = 12.85°** square, solid angle 4 arcsin(sin²(FOV/2)) = 0.0501 sr. Focal length trades scale against field: a longer lens gives finer pixels and fewer stars.

### 2.4 Diffraction: the smallest spot

Light is a wave (`P04.11`), and a circular aperture cannot focus a plane wave to a point. The image of a star is the **Airy pattern**, a bright core surrounded by faint rings, with the first dark ring at radius

**r_Airy = 1.22 λN**   (4)

in the focal plane (the angular form, 1.22λ/D, is the resolution limit of a telescope, 2.8). At λ = 600 nm and f/1.2, r_Airy = 0.88 μm = **0.16 pixel**; the best-fit Gaussian to the core has σ ≈ 0.42λN = 0.055 pixel. A perfectly focused star lands on one pixel, and then its position is known only to within that pixel. The cure, counter-intuitive and universal in star trackers, is to **defocus** the lens until the spot spans about 3 × 3 pixels (σ ≈ 1 pixel); 2.7 shows why.

### 2.5 How many stars, how many photons

Star brightness is measured in **magnitudes**: five magnitudes are a factor 100 in flux, so flux ∝ 10^{−0.4m}, and fainter stars have larger m. The number of stars over the whole sky brighter than m grows roughly exponentially; this lesson uses the empirical fit N(m) ≈ 6.57 e^{1.08m}, an assumption good to tens of per cent for 3 < m < 7 (it gives 4283 stars brighter than 6; catalogues list about 5000). With stars spread uniformly, the tracker's 0.0501 sr field holds on average 4283 × 0.0501/(4π) = **17.1 stars** brighter than magnitude 6: enough to identify the pattern (three or four are needed) and to average.

For photons, assume a broad silicon band receives about 10¹⁰ photons per second per square metre from a magnitude-0 star (an order-of-magnitude zero point, stated as an assumption). With D = 20.8 mm, a 100 ms exposure, 50 % quantum efficiency and 80 % optical transmission, a magnitude-6 star gives **543 photo-electrons**, magnitude 5 gives 1364, magnitude 3 gives 8603. Photon arrivals are random (Poisson), so the signal itself has a standard deviation √543 = 23 electrons.

### 2.6 Centroiding and its error

The **centroid** of a spot is its brightness-weighted mean position, x̂ = Σx_iS_i/ΣS_i over a window of pixels. For a spot made of N photons each landing at a random position with spread σ, x̂ is the mean of N independent samples, and its variance is σ²/N (`P02.16`). Each photon's position is further known only to its pixel, which adds a uniform error of variance 1/12 pixel². The **photon-limited centroid error** is therefore

**σ_c ≈ √(σ² + 1/12)/√N  pixels.**   (5)

For σ = 1 pixel and N = 543: 0.0447 pixel. The SYNTHETIC Monte Carlo (2000 images, true centres uniformly spread within a pixel, Poisson photons, seeded) gives 0.0454 pixel with a 9 × 9 window and no read noise. Read noise changes the picture: every pixel in the window adds noise of σ_r electrons whether or not it contains starlight, and the centroid multiplies it by the pixel's distance from the centre, so large windows are punished. SYNTHETIC results for magnitude 6 (rms error in pixels):

| read noise | 3 × 3 | 5 × 5 | 7 × 7 | 9 × 9 |
|---|---|---|---|---|
| 0 e⁻ | 0.1435 | 0.0522 | 0.0446 | 0.0454 |
| 3 e⁻ | 0.1450 | 0.0680 | 0.0930 | 0.1420 |
| 10 e⁻ | 0.1572 | 0.1532 | 0.2791 | 0.4770 |

The 3 × 3 window is too small for a σ = 1 pixel spot: it truncates the wings asymmetrically, a systematic error of 0.14 pixel. The best window balances truncation against read noise; at 3 e⁻ it is 5 × 5, giving **0.068 pixel = 3.1″ per star**. Real trackers also threshold the window to reject pixels below the noise, which this simple estimator does not do.

### 2.7 Why defocus: pixel locking

Without noise, a centroid should return the true position exactly. It does for a wide spot, but not for a narrow one: if the spot is much smaller than a pixel, almost all the light falls in one pixel whatever its position within it, and the centroid "locks" to the pixel centre. The bias, centroid minus truth, is an S-shaped function of the sub-pixel position: noise-free maximum |bias| **0.145 pixel at σ = 0.2 pixel, 0.0023 at 0.5, 0.0001 at 1.0**. A focused star tracker (σ = 0.055 pixel) would have a systematic error of up to 0.37 pixel, 17″, no matter how many photons it collected. Defocusing to σ ≈ 1 pixel removes the bias at the cost of spreading the light over more pixels (more read noise), and the optimum is the 3 × 3-to-5 × 5 spot used here.

### 2.8 From centroids to attitude

Each star's centroid gives a direction with error σ_c × (p/f). With n stars spread over the field, the errors average for the two axes across the boresight: **σ_att ≈ σ_c(p/f)/√n**, here 3.1″/√17 = **0.75″**. The roll about the boresight is determined by the stars' positions around the axis, with a lever arm of half the field rather than the full focal length, and is typically an order of magnitude worse (Extension 3). The budget omits calibration of the focal length and lens distortion, which in practice set a floor of about an arcsecond, and the motion of the spacecraft during the exposure.

### 2.9 Four sensor archetypes

**The pitot-static probe.** A forward-facing tube brings the air to rest and measures the total pressure p_t; side holes measure the static pressure p. For incompressible flow, Bernoulli gives **V = √(2(p_t − p)/ρ)**: 0.6 kPa at sea level is 31.30 m/s. At higher speed the air compresses as it stops, and the subsonic isentropic relation q_c = p[(1 + (γ − 1)M²/2)^{γ/(γ−1)} − 1] replaces ½ρV²; the airspeed indicator inverts it with sea-level constants to give the **calibrated airspeed**. At Mach 0.6 at sea level q_c = 27.92 kPa; the compressible formula returns 204.18 m/s, correct, and Bernoulli would read 213.49 m/s, 4.6 % high (Leishman, "Airspeed Definitions & Measurement", treats the airspeed definitions and the position error of the static port).

**The strain gauge and the Wheatstone bridge.** A strain gauge is a foil resistor bonded to a structure; stretched by strain ε its resistance changes by ΔR/R = GF·ε, with gauge factor GF ≈ 2 (SparkFun, "Strain Gauge Basics"). At 1000 microstrain the change is 0.2 %, far too small to measure as an absolute resistance, so the gauge is placed in a **Wheatstone bridge**, two dividers whose outputs are subtracted. With one active gauge and three fixed resistors (a quarter bridge), x = GF·ε,

**V_out = V_ex x/(2(2 + x)) ≈ V_ex GF ε/4:**   (6)

at 3.3 V excitation and 1000 με, 1.6484 mV (the linear formula gives 1.6500 mV). Four active gauges, two in tension and two in compression on a bending beam, give V_ex GF ε exactly, four times larger and linear: 6.600 mV. Millivolts need the amplifier of `P04.10` before the ADC.

**The thermistor.** An NTC thermistor's resistance falls steeply with temperature (OpenStax Vol. 2, 9.3 gives the physics of temperature-dependent resistance; a semiconductor's falls because more carriers are excited). The beta model R(T) = R₀ exp(B(1/T − 1/T₀)) with R₀ = 10 kΩ at 25 °C and B = 3950 K gives 105 385 Ω at −20 °C and 2486 Ω at +60 °C. In a divider with a 10 kΩ top resistor on 3.3 V, the output is 3.014, 1.650 and 0.657 V at −20, 25 and 60 °C: one ADC pin reads a battery's temperature (`P04.15` measures exactly this circuit). Thermocouples, the other common temperature sensor, generate tens of microvolts per kelvin from the junction of two metals and need amplification and a reference-junction temperature.

**The accelerometer.** A proof mass on a spring deflects by x = a/ω_n² under a steady acceleration a (SparkFun, "How an Accelerometer Works"); a MEMS device with a 5 kHz resonance deflects 9.94 nm under 1 g, measured capacitively. A high natural frequency gives bandwidth and a small deflection, which is why MEMS accelerometers are noisy. The range is a design choice: a ±16 g part read with 12 bits resolves 7.81 mg per count, and 1 mg over 32 g needs 15 bits (Fiore 12.3; SparkFun, "How to Select an Accelerometer").

### 2.10 The chain, end to end

Each archetype fits the chain of `P04.10`: transducer → (bridge or divider) → amplifier with gain chosen to fill the ADC range → anti-aliasing filter → ADC → counts → engineering unit by the inverse of the transducer law. The **resolution** in the physical unit is the ADC's LSB divided by the chain's sensitivity (V per unit), and the **accuracy** is something else entirely: calibration, temperature drift and nonlinearity. A 12-bit converter on a ±1500 με full bridge resolves 0.76 με (Problem 5); whether the reading is right to 10 με depends on the gauge's bonding and temperature compensation, not on the ADC.

## 3. Build the computational model

Open `starter/main.py`. The `StarTracker` dataclass and `main` are given.

### Step 1: geometric optics

`refraction_angle` returns `nan` where |n₁ sin θ₁/n₂| > 1; `critical_angle` raises for n₂ ≥ n₁; `fresnel_reflectance_normal`; `thin_lens_image_distance` (accept `np.inf`); `magnification`; `f_number`; `pixel_scale`; `field_of_view` (3); `airy_radius` (4) and `airy_gaussian_sigma`.

### Step 2: stars and photons

`stars_brighter_than` uses the module constants `STAR_COUNT_A`, `STAR_COUNT_B`; `fov_solid_angle` is 4 arcsin(sin(a/2) sin(b/2)); `mean_stars_in_fov`; `star_electrons` with `PHOTON_FLUX_MAG0`.

### Step 3: SYNTHETIC images and centroids

`gaussian_spot(x0, y0, sigma, flux, size)` integrates the Gaussian over each pixel exactly with the error function, so the image is right even when the spot is smaller than a pixel:

```python
edges = np.arange(size + 1) - 0.5 - size // 2
fx = 0.5 * np.diff(erf((edges - x0) / (np.sqrt(2) * sigma)))
fy = 0.5 * np.diff(erf((edges - y0) / (np.sqrt(2) * sigma)))
return flux * np.outer(fy, fx)
```

`centroid` subtracts a background and returns the weighted mean (x, y) relative to the central pixel. `simulate_centroid_errors` draws a true position in the central pixel, `rng.poisson` for photons and `rng.normal` for read noise, with `numpy.random.default_rng(seed)`. `pixel_locking_bias` scans x₀ without noise. `centroid_sigma_photon_limit` is σ/√N; `attitude_error_cross_boresight` divides by √n.

### Step 4: the sensors

`pitot_airspeed`, `pitot_delta_p_compressible`, `calibrated_airspeed` (2.9); `bridge_output` for any four resistors, `quarter_bridge` (6) and `full_bridge`; `thermistor_resistance`, `thermistor_temperature`, `divider_voltage`; `accelerometer_deflection`; `adc_resolution`, `bits_needed`.

Running `main()` prints

```
glass n = 1.5 to air: critical angle 41.81 deg; normal reflectance 4.00 %; 30 deg in air -> 19.47 deg in glass
fibre core 1.48 / cladding 1.46: critical angle 80.57 deg
thin lens f = 25 mm: object at 1 m -> image at 25.641 mm; at infinity -> 25.000 mm; at 50 mm -> 50.0 mm, m = -1.00

star tracker: f = 25 mm, f/1.2 (D = 20.8 mm), 1024 x 1024 px of 5.5 um
  pixel scale 220 urad = 45.4 arcsec/px; FOV 12.85 deg square, 0.0501 sr
  Airy radius at 600 nm: 0.88 um = 0.16 px (Gaussian sigma 0.055 px): the focused star is smaller than a pixel
  stars brighter than m = 6 over the sky (assumed fit): 4283; expected in FOV: 17.1
  m = 3 star, 100 ms: 8603 electrons
  m = 5 star, 100 ms: 1364 electrons
  m = 6 star, 100 ms: 543 electrons
  SYNTHETIC centroiding, m = 6, defocused sigma = 1 px, 9 x 9 window, no read noise: rms 0.0454 px (sigma/sqrt(N) = 0.0429; with the pixel's own 1/12 px^2: 0.0447)
```

then the window-and-read-noise table of 2.6, the pixel-locking line of 2.7,

```
  baseline 5 x 5 window, 3 e- read noise: 0.0680 px = 3.09 arcsec per star; cross-boresight attitude error with 17 stars: 0.75 arcsec
```

and the sensor lines of 2.9, and writes `star_tracker.png` (a SYNTHETIC magnitude-6 star image with its true and measured centres; the pixel-locking S-curves; the centroid error against electrons with the σ/√N line).

### Verify

Run `python -m pytest curriculum/phase-04/lesson-12-optics-sensors`. The tests:

1. **Snell**: n₁ sin θ₁ = n₂ sin θ₂; normal incidence undeviated; the reverse path returns the original angle; bending toward the normal.
2. **TIR**: sin θ_c = n₂/n₁, 41.81°; finite just below, `nan` just above, grazing at θ_c; n₂ ≥ n₁ raises; 4 % reflectance, symmetric, zero for equal indices.
3. **Thin lens**: (2) on four objects; infinity → f; 2f → 2f with m = −1; virtual image inside f.
4. **Camera**: N = f/D; p/f; the arctan FOV below n·p/f and equal to it for 10 pixels to 10⁻⁶; solid angle → FOV² for small fields and 2π for a hemisphere.
5. **Airy**: 1.22λN, doubling with N and with λ; Gaussian σ smaller; below one pixel for the baseline.
6. **Counts**: e^{B} per magnitude; stars ∝ solid angle and all of them over 4π; electrons ×100 for 5 magnitudes, ∝ D² and ∝ t.
7. **Spot and centroid**: flux conserved to 10⁻⁶ on an 11 × 11 window; symmetric; exact centroid when centred and within 10⁻⁴ off centre for σ = 1; background subtraction.
8. **1/√N**: SYNTHETIC rms over 200–3200 electrons with log–log slope −0.50 ± 0.07 and level √(σ² + 1/12)/√N within 12 %; mean within four standard errors of zero; identical with the same seed.
9. **Read noise**: a 9 × 9 window more than 1.5 times worse than 5 × 5 at 10 e⁻.
10. **Pixel locking**: above 0.1 px for σ = 0.2, below 10⁻³ for σ = 1, odd in x₀, zero at the centre; attitude error ∝ 1/√n.
11. **Pitot**: Bernoulli inverts ½ρV²; the compressible chain returns M·a₀ at sea level for M = 0.05–0.9; Bernoulli agrees at M = 0.05 to 10⁻³ and overreads by more than 3 % at 0.6.
12. **Bridge**: balance for equal resistors and for R₁R₃ = R₂R₄; `quarter_bridge` equals `bridge_output` exactly; linear within 2 × 10⁻⁴ at ±100 με, visibly not at 2000 με; full bridge = 4 × the linear quarter bridge.
13. **Thermistor**: 10 kΩ at 25 °C, decreasing, round trip to 10⁻¹², divider 1.65 V at 25 °C and decreasing.
14. **Accelerometer and ADC**: x = a/ω_n², ¼ for twice f_n; 32/4096; 15 bits for 1 mg over 32 g, 14 not enough.
15. **main** writes the figure.

### Compare

There is no library star tracker to compare with, so the checks are internal and analytic: the pixel-integrated Gaussian conserves flux and has an exact centroid; the Monte Carlo matches the variance law (5), including the 1/12 term that a naïve σ/√N misses (0.0454 measured, 0.0447 with the term, 0.0429 without); the noise-free centroid of a wide spot is unbiased to 10⁻⁴ pixel. The sensors are compared with their own limits: Bernoulli as the low-Mach limit of the isentropic pitot law, the linear bridge as the small-strain limit of the exact one. The photometric zero point and the star-count fit are the weakest inputs; a design review would replace them with a catalogue query (for example, the stars of the Hipparcos catalogue in each field) and the detector's measured quantum efficiency.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction first.

1. **Design sweep.** For f = 16, 25, 35, 50 mm at f/1.2 with the baseline detector, tabulate the pixel scale, FOV, expected stars to magnitude 6, electrons for a magnitude-6 star in 100 ms (D grows with f), the SYNTHETIC centroid error (5 × 5, 3 e⁻, σ = 1 pixel) and the cross-boresight attitude error. Predict the attitude error to fall roughly as 1/f² (finer pixels and more photons per star) until too few stars remain.
2. **Signal and noise.** For electrons N = 100 … 10⁴ and read noise 0, 3, 10 e⁻, the SYNTHETIC centroid error with the best window at each point. Predict the photon-limited slope −½ at high N and a steeper read-noise-limited slope −1 at low N.
3. **Defocus.** For σ = 0.2 … 2 pixels at N = 543 and 3 e⁻, the total error (pixel-locking bias plus noise). Predict an optimum near σ = 0.7–1 pixel.
4. **Sensor resolution.** For each archetype, the physical resolution of a 12-bit, 3.3 V ADC behind it: the pitot with a ±2 kPa differential sensor of 1 mV/Pa (in m/s at 30 m/s and 100 m/s), the full bridge at 3.3 V excitation and a gain of 160 (in με), the thermistor divider (in K at −20, 25 and 60 °C), a ±16 g accelerometer (in mg). Predict which are limited by the ADC and which by the sensor.

Raw data format: CSV tables `f_mm, scale_arcsec, fov_deg, n_stars, electrons_m6, centroid_px_SYNTHETIC, attitude_arcsec` (1); `electrons, read_noise_e, window, rms_px_SYNTHETIC` (2); `sigma_px, bias_px, noise_px_SYNTHETIC, total_px` (3); `sensor, operating_point, resolution, unit` (4); plus `star_tracker.png`. Every SYNTHETIC column says so in its header and the seed is recorded.

What you should find (for the baseline and the checks in `main()`): at f = 25 mm, 45.4″ per pixel, 17 stars, 543 electrons, 0.068 pixel and 0.75″; the slope −0.5 at high N; pixel-locking bias negligible above σ = 0.5 pixel; for the thermistor about 0.02–0.05 K per count across −20 to 60 °C with the top resistor matched at 20 °C (Problem 3); for the accelerometer 7.81 mg per count. The remaining entries are yours.

## 5. Validate and improve

| quantity | theory (§2) | code | difference | cause |
|---|---|---|---|---|
| critical angle, glass | arcsin(1/1.5) = 41.81° | 41.81° | 0 | — |
| image of a 1 m object, f = 25 mm | 25.641 mm | 25.641 mm | 0 | — |
| pixel scale | p/f = 45.4″ | 45.4″ | 0 | — |
| FOV | 2 arctan(5.632/50) = 12.85° | 12.85° | 0 | — |
| Airy radius | 1.22λN = 0.88 μm | 0.88 μm | 0 | — |
| centroid, photon-limited | √(1 + 1/12)/√543 = 0.0447 px | SYNTHETIC 0.0454 px | 1.6 % | Monte Carlo (2000 images: ±1.6 %) |
| pixel locking, σ = 1 px | 0 | 10⁻⁴ px | — | Gaussian tails beyond 9 × 9 |
| pitot, M = 0.6 | M a₀ = 204.18 m/s | 204.18 m/s | 0 | — |
| quarter bridge, 1000 με | linear 1.6500 mV | exact 1.6484 mV | 0.1 % | x/(2 + x) nonlinearity |
| thermistor at 25 °C | 10 kΩ, 1.650 V | same | 0 | — |

The analytic parts are exact; the one statistical comparison, the Monte Carlo against (5), agrees within its own sampling error. That validates the code and the variance law, not the star tracker: the attitude accuracy of 0.75″ rests on the photometric zero point (a factor of 2 fewer photons is a factor of about 1.6 in centroid error at 3 e⁻ read noise), the uniform star density (a field near a galactic pole may hold half as many stars), and the absence of calibration errors, which in real instruments set a floor near one arcsecond. The design decision that does not depend on those assumptions is the **defocus**: at f/1.2 a focused star is 0.16 pixel and pixel locking would cost up to 0.37 pixel; spreading it to σ ≈ 1 pixel costs a little read noise and removes the bias.

What to improve: a catalogue-based star count for the actual pointing directions; a defocused Airy spot (Extension 1); thresholded windows; rotation smear (Extension 2). What not to improve: the Monte Carlo beyond 2000 images, since its 1.6 % uncertainty is far below that of the photometric assumptions.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (3), (5) and (6), and of the compressible pitot law's low-Mach limit; `star_tracker.png`; the four experiment tables, SYNTHETIC ones labelled with their seeds; and a one-page star-tracker budget (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) for a stated attitude requirement: optics, detector, defocus, exposure, expected stars, centroid and attitude accuracy, and every assumption.

### Adversarial self-checks

- **What if the exposure halves?** Electrons halve to 272 for a magnitude-6 star; the photon-limited error would grow by √2, and with read noise it grows more, from 0.067 to 0.108 pixel (SYNTHETIC, 5 × 5, 3 e⁻); the attitude error rises from about 0.75″ to about 1.2″. In exchange, rotation smear halves, which matters more on an agile spacecraft.
- **What if the spot is perfectly focused?** σ = 0.055 pixel: pixel locking dominates with biases of up to 0.37 pixel (17″ per star) that do not average down with more photons; averaging over 17 stars at random sub-pixel positions reduces it only as 1/√17. Defocus is not optional.
- **Which assumption fails first in the real world?** Stray light: the Sun, the Earth or the Moon near the field raises the background by orders of magnitude, and a baffle with an exclusion angle of 30–45° is part of every star tracker; then the uniform-star assumption near the galactic poles; then the photometric zero point.
- **What measurement could falsify the budget?** Night-sky images through the actual lens and detector: count the detected stars against the prediction for the field, and measure the scatter of centroids of the same star across consecutive frames; a scatter well above the predicted 0.07 pixel means more noise (or atmospheric seeing, on the ground) than the model has.
- **Does the result survive a change of seed or estimator?** Different seeds change the Monte Carlo rms by about 1.6 %; a thresholded estimator improves the 7 × 7 and 9 × 9 cases substantially, so the 5 × 5 optimum is a property of the simple estimator, not of the physics.
- **Which decision changes if the photon count is a factor two lower than assumed?** The centroid error rises by a factor of about 1.6 (read noise weighs more at low signal) and the attitude error to about 1.2″; if the requirement is 1″, the design then needs a longer exposure or a larger aperture. The defocus and the window size do not change.

## Check yourself

1. Light passes from glass (n = 1.5) into air. Above what angle of incidence is it totally reflected, and what fraction is reflected at normal incidence from an uncoated surface?
2. Where does a thin lens of focal length f image a star, and why does a star tracker's detector sit there?
3. A 1024-pixel detector of 5.5 μm pitch sits behind a 25 mm lens. What are the pixel scale in arcseconds and the field of view in degrees?
4. Why does a star tracker defocus its stars, and what goes wrong if it does not?
5. A star gives 543 photo-electrons spread with σ = 1 pixel. Estimate the centroid error, and say why it is slightly more than σ/√N.
6. Why does a large centroiding window hurt when there is read noise?
7. A quarter-bridge strain gauge (GF = 2) at 3.3 V sees 1000 με. Roughly what is the output, and why does a full bridge help?
8. A 12-bit ADC reads a ±16 g accelerometer. What is one count in milli-g, and how many bits would 1 mg need?
