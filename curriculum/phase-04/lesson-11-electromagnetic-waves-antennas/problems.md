# P04.11  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical. Gains in dBi, powers in W or dBm.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: how weak is a navigation signal?)

A navigation satellite at 20 200 km altitude radiates an assumed 50 W through an antenna with 13 dBi gain toward the Earth. A receiver at the sub-satellite point uses an isotropic antenna at 1575.42 MHz.

(a) Compute the effective isotropic radiated power (EIRP), the intensity at the receiver and the peak electric field.
(b) Compute the effective aperture of the receiving antenna and the received power in W and dBm.
(c) Compare with the thermal noise power kTB in a 2 MHz bandwidth at 290 K. Is the signal above or below the noise? (Navigation receivers recover it by correlating over many code periods; you need not explain how.)

Answer form: three numbers with units, two numbers, one comparison in dB.

## 2. (Derivation: waves, directivity, Friis and geometry)

(a) From Maxwell's equations in vacuum, derive the wave equation for E and the wave speed. Show that a plane wave E_x = E₀ cos(kz − ωt) has B_y = E_x/c.
(b) Show by hand that the short dipole's directivity is exactly 3/2. Then show that any pattern U(θ) = sinⁿθ has D = 4π/(2π∫sinⁿ⁺¹θ dθ) and evaluate it for n = 4.
(c) Derive the Friis equation (8) from the inverse-square law (3), the definition of gain, and the effective aperture (7). State every assumption it carries.
(d) Derive the slant range (9) from the triangle formed by Earth's centre, the ground station and the satellite.

Answer form: four derivations, one number.

## 3. (Unfamiliar scenario: a solar-sail cargo tug)

A flat, perfectly reflecting sail of area A = (10 + s)² m² carries a total mass of 10 + s/3 kg at 1 au.

(a) Compute the force and acceleration facing the Sun, and the ratio to the Sun's gravitational acceleration at 1 au (GM_sun = 1.32712 × 10²⁰ m³/s²). Explain why the ratio does not depend on the distance from the Sun.
(b) To raise the orbit, the sail must push along the orbital velocity, which is perpendicular to the sunline. Show that the transverse component of the force is maximised at a cone angle α with tan α = 1/√2, and compute the transverse acceleration there.
(c) Starting from the transverse acceleration, estimate the Δv accumulated in one year. Why is a sail a propulsion system with "infinite specific impulse" and what limits it in practice?

Answer form: three numbers and a sentence, a derivation and a number, a number and two sentences.

## 4. (Numerical, uses your code: an S-band downlink)

A CubeSat moves from UHF to S band (2.25 GHz): 2 W transmitter, 6 dBi circularly polarised patch antenna, ground dish of diameter D = 0.5 + s/30 m with 55 % aperture efficiency, T_sys = 300 K, 4 dB of other losses, required E_b/N₀ = 10 dB, altitude 500 km.

(a) Compute the dish gain, its effective aperture and its approximate beamwidth λ/D in degrees. Why must the dish track the satellite?
(b) Using `Link`, compute the margin at 10° elevation for 1 Mbit/s and the highest bit rate with a 3 dB margin.
(c) Compare with the lesson's UHF link at the same 1 Mbit/s and 2 W. Explain the difference in terms of the Friis equation: what does the higher frequency cost, and what does the dish buy?

Answer form: three numbers, two numbers, one comparison and two sentences.

## 5. (Numerical, uses your code: a dipole of your own length)

A dipole of total length L = 0.25 + s/40 wavelengths is centre-fed.

(a) With `dipole_pattern`, `directivity`, `directivity_quad` and `half_power_beamwidth`, find D (linear and dBi), the beamwidth and the direction of maximum.
(b) Find the radiation resistance referred to the current amplitude and the resistance at the feed, R_rad/sin²(kL/2). Would you feed it from a 50 Ω cable directly?
(c) Run `dipole_pattern_numeric` with 25, 50, 100, 200, 400 elements and tabulate the maximum error against (5) and the observed order. Explain any pre-asymptotic behaviour.

Answer form: four numbers, two resistances and a sentence, a 5-row table.

## 6. (Judgment: antennas on spacecraft)

Answer each in two or three sentences, with a number where natural.

(a) A deep-space probe carries both a 3 m high-gain dish and a low-gain antenna (see Basics of Space Flight, Chapter 11, "Telecommunications Subsystems"). Why both? Estimate the gain difference at 8.4 GHz between the dish (55 %) and a 2 dBi antenna.
(b) A CubeSat with a single dipole tumbles slowly. Describe what the ground station's received power does during a pass and what antenna change fixes it.
(c) A ground station designer wants to double the data rate. List three changes ranked by cost per decibel.
(d) Your link budget shows a 1 dB margin. A colleague says "it closes, ship it". What do you answer?

Answer form: four short paragraphs, one number.
