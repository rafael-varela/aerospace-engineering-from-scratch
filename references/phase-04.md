# Reference map: Phase 04: Physics Foundations

**Update 2026-10-16:** the GAP lessons 04-10 through 04-15 were re-mapped by loading the live pages of LibreTexts (Fiore op-amp text, Kahn digital-circuit text, Ellingson *Electromagnetics* Vol. 2), SparkFun Learn, the MicroPython docs, NASA *Basics of Space Flight*, NASA "Why Space Radiation Matters", NTRS 20170004674, DOE "What is a Radioisotope Power System?", the NASA-HDBK-4002B record page, MIT OCW 8.033 (session 23), and the PubMed Central copy of Ashby 2003; eleven ids were added to the manifest and `@mit-ocw-8-033-relativity` was extended. Not added: Arduino docs (docs.arduino.cc pages are JavaScript-rendered and show no titles or license to the checker; the repository README states CC BY-SA 4.0 but no page content could be verified), Adafruit Learn (no license line on the guide page), MIT OCW 6.013 textbook page (server error on two attempts), the Springer page of Ashby 2003 (rate-limited). Licenses added: CC BY 4.0 (Kahn), CC BY-SA 4.0 (Ellingson, SparkFun), CC BY-NC-SA 4.0 (Fiore, MIT OCW), MIT (MicroPython), US Government public domain (NASA, NTRS, DOE).

Verified 2026-10-11 by loading the live pages; section numbers and titles are copied from the publishers' own chapter outlines or section pages. Every source here is **required-free** unless marked otherwise; sources are identified by an `@`-prefixed identifier and recorded in [sources.yaml](sources.yaml) with license, version, and caveats. **GAP** marks a lesson that no listed free source covers adequately: the lesson must carry original instruction there, and any "candidate" named is unverified until it is added to the manifest.

Conventions: OpenStax section URLs follow `https://openstax.org/books/<book-slug>/pages/<n>-<m>-<title-slug>`; a URL written as `.../pages/...` continues the last full URL above it. MIT OCW video URLs follow `https://ocw.mit.edu/courses/<course-slug>/resources/lecture-<n>-<slug>/`. OpenStax and MIT OCW are **CC BY-NC-SA 4.0**, Leishman is **CC BY-NC-ND 4.0**: the curriculum links to them and writes its own explanations; nothing is copied or adapted.


## 04-01 Kinematics & Newton's Laws (Vol. 1)

Sources: `@openstax-university-physics-vol1`

- Ch. 3 "Motion Along a Straight Line": 3.1 Position, Displacement, and Average Velocity; 3.2 Instantaneous Velocity and Speed; 3.3 Average and Instantaneous Acceleration; 3.4 Motion with Constant Acceleration; 3.5 Free Fall
  - https://openstax.org/books/university-physics-volume-1/pages/3-1-position-displacement-and-average-velocity (loaded)
- Ch. 4 "Motion in Two and Three Dimensions": 4.3 Projectile Motion; 4.4 Uniform and Nonuniform Circular Motion; 4.5 Relative Motion in One and Two Dimensions
- Ch. 5 "Newton's Laws of Motion": 5.1 Forces; 5.2 Newton's First Law; 5.3 Newton's Second Law; 5.5 Newton's Third Law; 5.7 Drawing Free-Body Diagrams
- Ch. 6 "Applications of Newton's Laws": 6.1 Solving Problems with Newton's Laws; 6.2 Friction; 6.3 Centripetal Force; 6.4 Drag Force and Terminal Speed
  - https://openstax.org/books/university-physics-volume-1/pages/6-4-drag-force-and-terminal-speed (pattern)
- Extract: free-body diagrams as the core skill; relative motion (wind triangle); quadratic drag and terminal speed as first aerodynamics contact.

## 04-02 Work, Energy & Power (Vol. 1)

Sources: `@openstax-university-physics-vol1`

- Ch. 7 "Work and Kinetic Energy": 7.1 Work; 7.2 Kinetic Energy; 7.3 Work-Energy Theorem; 7.4 Power
- Ch. 8 "Potential Energy and Conservation of Energy": 8.1 Potential Energy of a System; 8.2 Conservative and Non-Conservative Forces; 8.3 Conservation of Energy; 8.4 Potential Energy Diagrams and Stability
  - https://openstax.org/books/university-physics-volume-1/pages/8-3-conservation-of-energy (pattern)
- Extract: energy bookkeeping with non-conservative losses; power = F·v (climb rate from excess power); reading a potential-energy diagram for stability.

## 04-03 Momentum & Collisions (Vol. 1)

Sources: `@openstax-university-physics-vol1`

- Ch. 9 "Linear Momentum and Collisions": 9.1 Linear Momentum; 9.2 Impulse and Collisions; 9.3 Conservation of Linear Momentum; 9.4 Types of Collisions; 9.6 Center of Mass; 9.7 Rocket Propulsion
  - https://openstax.org/books/university-physics-volume-1/pages/9-7-rocket-propulsion (pattern)
- Extract: impulse-momentum; elastic vs. inelastic; 9.7 derives the rocket equation (Tsiolkovsky) from momentum conservation — the key aerospace payoff.

## 04-04 Rotation & Angular Momentum (Vol. 1)

Sources: `@openstax-university-physics-vol1`

- Ch. 10 "Fixed-Axis Rotation": 10.1 Rotational Variables; 10.3 Relating Angular and Translational Quantities; 10.4 Moment of Inertia and Rotational Kinetic Energy; 10.5 Calculating Moments of Inertia; 10.6 Torque; 10.7 Newton's Second Law for Rotation
- Ch. 11 "Angular Momentum": 11.2 Angular Momentum; 11.3 Conservation of Angular Momentum; 11.4 Precession of a Gyroscope
  - https://openstax.org/books/university-physics-volume-1/pages/11-4-precession-of-a-gyroscope (pattern)
- Extract: I, torque, and L for a rigid body; conservation of L (reaction wheels, spin stabilization); gyroscopic precession (attitude sensors).

## 04-05 Gravitation & Kepler's Laws (Vol. 1)

Sources: `@openstax-university-physics-vol1`, `@leishman-intro-aerospace-flight-vehicles`

- Ch. 13 "Gravitation": 13.1 Newton's Law of Universal Gravitation; 13.2 Gravitation Near Earth's Surface; 13.3 Gravitational Potential Energy and Total Energy; 13.4 Satellite Orbits and Energy; 13.5 Kepler's Laws of Planetary Motion
  - https://openstax.org/books/university-physics-volume-1/pages/13-4-satellite-orbits-and-energy (pattern)
- Supplement (Leishman): chapter "Astronautics, Space & Astrodynamics" (title appears in site menu; URL not confirmed, see caveat below).
- Extract: circular orbit speed and period, specific orbital energy, escape velocity, Kepler's three laws with the vis-viva relation.

## 04-06 Oscillations & Waves (Vol. 1)

Sources: `@openstax-university-physics-vol1`

- Ch. 15 "Oscillations": 15.1 Simple Harmonic Motion; 15.2 Energy in Simple Harmonic Motion; 15.4 Pendulums; 15.5 Damped Oscillations; 15.6 Forced Oscillations
- Ch. 16 "Waves": 16.1 Traveling Waves; 16.2 Mathematics of Waves; 16.4 Energy and Power of a Wave; 16.5 Interference of Waves; 16.6 Standing Waves and Resonance
  - https://openstax.org/books/university-physics-volume-1/pages/15-6-forced-oscillations (pattern)
- Extract: natural frequency, damping ratio, resonance curve (structural vibration, flutter preview); wave equation and standing waves.

## 04-07 Kinetic Theory of Gases (Vol. 2)

Sources: `@openstax-university-physics-vol2`, `@leishman-intro-aerospace-flight-vehicles`

- Ch. 2 "The Kinetic Theory of Gases": 2.1 Molecular Model of an Ideal Gas; 2.2 Pressure, Temperature, and RMS Speed; 2.3 Heat Capacity and Equipartition of Energy; 2.4 Distribution of Molecular Speeds
  - https://openstax.org/books/university-physics-volume-2/pages/2-2-pressure-temperature-and-rms-speed (pattern)
- Supplement (Leishman): "Atmospheric Properties & the ISA" — https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/international-standard-atmosphere-isa/ (loaded; heading "13 Atmospheric Properties & the ISA" on 2026-10-11)
- Extract: pressure from molecular momentum transfer; RMS speed vs. speed of sound; gamma = Cp/Cv from equipartition; then apply to the ISA model.

## 04-08 Electrostatics & Circuits (Vol. 2)

Sources: `@openstax-university-physics-vol2`

- Ch. 5 "Electric Charges and Fields": 5.1 Electric Charge; 5.3 Coulomb's Law; 5.4 Electric Field
- Ch. 7 "Electric Potential": 7.2 Electric Potential and Potential Difference; 7.5 Equipotential Surfaces and Conductors
- Ch. 9 "Current and Resistance": 9.1 Electrical Current; 9.3 Resistivity and Resistance; 9.4 Ohm's Law; 9.5 Electrical Energy and Power
- Ch. 10 "Direct-Current Circuits": 10.1 Electromotive Force; 10.2 Resistors in Series and Parallel; 10.3 Kirchhoff's Rules; 10.4 Electrical Measuring Instruments; 10.5 RC Circuits
  - https://openstax.org/books/university-physics-volume-2/pages/10-3-kirchhoffs-rules (pattern)
- Extract: field/potential vocabulary; Ohm's law and power; Kirchhoff's rules on a two-loop circuit; RC time constant; how a voltmeter/ammeter loads a circuit.

## 04-09 Magnetism & Induction (Vol. 2)

Sources: `@openstax-university-physics-vol2`

- Ch. 11 "Magnetic Forces and Fields": 11.2 Magnetic Fields and Lines; 11.3 Motion of a Charged Particle in a Magnetic Field; 11.4 Magnetic Force on a Current-Carrying Conductor; 11.5 Force and Torque on a Current Loop; 11.6 The Hall Effect
- Ch. 12 "Sources of Magnetic Fields": 12.2 Magnetic Field Due to a Thin Straight Wire; 12.5 Ampère's Law; 12.6 Solenoids and Toroids; 12.7 Magnetism in Matter
- Ch. 13 "Electromagnetic Induction": 13.1 Faraday's Law; 13.2 Lenz's Law; 13.3 Motional Emf; 13.5 Eddy Currents; 13.6 Electric Generators and Back Emf
- Ch. 14 "Inductance": 14.2 Self-Inductance and Inductors; 14.4 RL Circuits
  - https://openstax.org/books/university-physics-volume-2/pages/13-1-faradays-law (pattern)
- Extract: Lorentz force (magnetorquers, Hall sensors); Faraday/Lenz (motors, generators, back-emf in brushless motors); inductor time constant.

## 04-10 AC Circuits & Electronics Basics (Vol. 2 + Vol. 3)

Sources: `@openstax-university-physics-vol2`, `@openstax-university-physics-vol3`, `@fiore-opamps-linear-ics`, `@kahn-digital-circuit-projects`, `@sparkfun-learn-tutorials`, `@micropython-docs`

- Vol. 2 Ch. 14 "Inductance": 14.5 Oscillations in an LC Circuit; 14.6 RLC Series Circuits
- Vol. 2 Ch. 15 "Alternating-Current Circuits": 15.1 AC Sources; 15.2 Simple AC Circuits; 15.3 RLC Series Circuits with AC; 15.4 Power in an AC Circuit; 15.5 Resonance in an AC Circuit; 15.6 Transformers
  - https://openstax.org/books/university-physics-volume-2/pages/15-3-rlc-series-circuits-with-ac (loaded)
- Vol. 3 Ch. 9 "Condensed Matter Physics": 9.5 Band Theory of Solids; 9.6 Semiconductors and Doping; 9.7 Semiconductor Devices
  - https://openstax.org/books/university-physics-volume-3/pages/9-7-semiconductor-devices (pattern)
- Op-amps (Fiore): Ch. 4 "Basic Op Amp Circuits": 4.1 Introduction; 4.2 Inverting and Noninverting Amplifiers; 4.3 Single Supply Biasing; 4.4 Current Boosting; 4.5 Summary; 4.6 Problems (loaded 2026-10-16). Book chapters 2 "Operational Amplifier Internals" and 3 "Negative Feedback" exist for the learner who wants the model behind the ideal-op-amp rules (titles from the book TOC; not loaded).
  - https://eng.libretexts.org/Bookshelves/Electrical_Engineering/Electronics/Operational_Amplifiers_and_Linear_Integrated_Circuits_-_Theory_and_Application_(Fiore)/04%3A_Basic_Op_Amp_Circuits (loaded)
- Logic gates (Kahn): Ch. 4 "Gates"; Ch. 5 "Associative Boolean Operators"; Ch. 6 "Adders"; Ch. 9 "Memory Basics - Flip-Flops and Latches" (chapter titles from the book TOC, loaded 2026-10-16). SparkFun "Digital Logic": sections "Combinational Logic" (Truth Tables; Written Boolean Logic), "Sequential Logic" (D-type Flip-Flop; Setup, Hold, and Propagation Times), "Boolean Logic in Programming" (Bitwise Logic; Bit Shifting).
  - https://learn.sparkfun.com/tutorials/digital-logic/all (loaded)
- Microcontroller basics: SparkFun "What is an Arduino?": "What Does it Do?"; "What's on the board?" (Power; Pins (5V, 3.3V, GND, Analog, Digital, PWM, AREF); Main IC; Voltage Regulator). MicroPython "Quick reference for the RP2": "General board control"; "Pins and GPIO"; "PWM (pulse width modulation)"; "ADC (analog to digital conversion)"; "Hardware I2C bus"; "Hardware SPI bus".
  - https://learn.sparkfun.com/tutorials/what-is-an-arduino/all (loaded); https://docs.micropython.org/en/latest/rp2/quickref.html (loaded)
- Extract: impedance and phasors, RLC resonance and filters, RMS power; p-n junction, diode and transistor at the conceptual level; the two ideal-op-amp rules applied to the inverting and non-inverting amplifier (Fiore 4.2); truth tables and a D flip-flop as the memory primitive; a microcontroller as GPIO + timers + ADC + serial buses around a CPU.
- Remaining partial GAP: OpenStax has no op-amp, logic gate, or microcontroller content; the sources above fill it. No single free source gives a sensor-amplifier-ADC-microcontroller chain end to end; the lesson writes that connecting narrative.

## 04-11 Electromagnetic Waves & Antennas (Vol. 2)

Sources: `@openstax-university-physics-vol2`, `@ellingson-electromagnetics-vol2`, `@nasa-basics-of-space-flight`

- Ch. 16 "Electromagnetic Waves": 16.1 Maxwell's Equations and Electromagnetic Waves; 16.2 Plane Electromagnetic Waves; 16.3 Energy Carried by Electromagnetic Waves; 16.4 Momentum and Radiation Pressure; 16.5 The Electromagnetic Spectrum
  - https://openstax.org/books/university-physics-volume-2/pages/16-5-the-electromagnetic-spectrum (pattern)
- Antennas (Ellingson Vol. 2): Ch. 9 "Radiation": 9.4 Radiation from a Hertzian Dipole; 9.5 Radiation from an Electrically-Short Dipole; 9.7 Far-Field Radiation from a Half-Wave Dipole. Ch. 10 "Antennas": 10.1 How Antennas Radiate; 10.2 Power Radiated by an Electrically-Short Dipole; 10.5 Equivalent Circuit Model for Transmission; Radiation Efficiency; 10.7 Directivity and Gain; 10.8 Radiation Pattern; 10.10 Reciprocity; 10.13 Effective Aperture; 10.14 Friis Transmission Equation (section titles loaded 2026-10-16 from the chapter pages).
  - https://eng.libretexts.org/Bookshelves/Electrical_Engineering/Electro-Optics/Book%3A_Electromagnetics_II_(Ellingson)/10%3A_Antennas (loaded); .../09%3A_Radiation (loaded)
- Spacecraft context (NASA Basics of Space Flight): Chapter 10 "Telecom": "Transmitters and Receivers"; "Signal Power"; "Uplink and Downlink"; "Modulation and Demodulation, Carrier and Subcarrier"; "Symbols and Bits and Coding". Chapter 11 "Onboard Systems" page 2: "Telecommunications Subsystems" with "High-Gain Antennas (HGA)", "Low-Gain Antennas", "Medium-gain Antennas", "Spacecraft Transmitters", "Spacecraft Receivers".
  - https://science.nasa.gov/learn/basics-of-space-flight/chapter10-1/ (loaded); https://science.nasa.gov/learn/basics-of-space-flight/chapter11-2/ (loaded)
- Extract: c from Maxwell, E/B relation, Poynting intensity and inverse-square falloff, radiation pressure (solar sails), band names; dipole pattern and directivity, gain = efficiency x directivity, effective aperture, and the Friis equation as the link budget; HGA vs. LGA trade on a real spacecraft.
- Remaining partial GAP: polarization of antennas (linear vs. circular, polarization-mismatch loss) is not a titled section in Ellingson Ch. 10 or in Basics of Space Flight; the lesson supplies it, building on OpenStax 16.2 for wave polarization.

## 04-12 Optics & Sensors (Vol. 3)

Sources: `@openstax-university-physics-vol3`, `@openstax-university-physics-vol2`, `@leishman-intro-aerospace-flight-vehicles`, `@sparkfun-learn-tutorials`, `@fiore-opamps-linear-ics`, `@nasa-basics-of-space-flight`

- Ch. 1 "The Nature of Light": 1.1 The Propagation of Light; 1.2 The Law of Reflection; 1.3 Refraction; 1.4 Total Internal Reflection; 1.7 Polarization
- Ch. 2 "Geometric Optics and Image Formation": 2.4 Thin Lenses; 2.6 The Camera; 2.8 Microscopes and Telescopes
  - https://openstax.org/books/university-physics-volume-3/pages/2-6-the-camera (pattern)
- Supplement (Leishman): "Airspeed Definitions & Measurement" — https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/determination-of-airspeed/ (loaded; heading "41 Airspeed Definitions & Measurement"; covers pitot-static measurement and position error)
- Strain gauges (SparkFun "Getting Started with Load Cells"): "Load Cell Basics" (Types of Load Cells: Strain Gauge Load Cells); "Strain Gauge Basics" (Gauge Factor; Small Changes in Strain; Amplifiers and Wheatstone Bridge).
  - https://learn.sparkfun.com/tutorials/getting-started-with-load-cells/all (loaded 2026-10-16)
- Accelerometers (SparkFun "Accelerometer Basics"): "What is an Accelerometer?"; "How an Accelerometer Works"; "How to Connect to an Accelerometer" (Communication Interface; Power); "How to Select an Accelerometer" (Range).
  - https://learn.sparkfun.com/tutorials/accelerometer-basics/all (loaded)
- Thermistors / resistive temperature sensing: OpenStax Vol. 2 9.3 Resistivity and Resistance (temperature dependence of resistance; already mapped in 04-08) is the physics; no free source with a thermistor-specific section was verified.
- ADC sampling (Fiore): Ch. 12 "Analog-to-Digital-to-Analog Conversion": 12.1 Introduction to Analog to Digital to Analog Conversion; 12.2 The Sampling Theorem; 12.3 Resolution and Sampling Rate; 12.5 Analog-to-Digital Conversion. SparkFun "Analog to Digital Conversion": "What is the ADC?"; "Relating ADC Value to Voltage".
  - https://eng.libretexts.org/Bookshelves/Electrical_Engineering/Electronics/Operational_Amplifiers_and_Linear_Integrated_Circuits_-_Theory_and_Application_(Fiore)/12%3A_Analog-to-Digital-to-Analog_Conversion (loaded); https://learn.sparkfun.com/tutorials/analog-to-digital-conversion/all (loaded)
- Star trackers in context (NASA Basics of Space Flight): Chapter 11 "Onboard Systems" page 2, "Attitude and Articulation Control Subsystems (AACS)" (run-in labels Spin, 3-Axis, Articulation, Celestial Reference, Inertial Reference). Qualitative only.
  - https://science.nasa.gov/learn/basics-of-space-flight/chapter11-2/ (loaded)
- Extract: thin-lens equation and f-number for a camera/star tracker; total internal reflection (fiber); pitot-static as the archetypal aerospace sensor; gauge factor and the Wheatstone bridge; accelerometer range and interface; Nyquist rate, resolution (LSB) and quantization for an ADC.
- Remaining partial GAP: thermocouples and a quantitative star-tracker optical design (plate scale, centroiding) have no verified free source; the lesson supplies them.

## 04-13 Special Relativity for GPS (Vol. 3)

Sources: `@openstax-university-physics-vol1`, `@openstax-university-physics-vol3`, `@mit-ocw-8-033-relativity`, `@ashby-relativity-in-gps`

- Ch. 5 "Relativity": 5.1 Invariance of Physical Laws; 5.3 Time Dilation; 5.4 Length Contraction; 5.5 The Lorentz Transformation; 5.7 Doppler Effect for Light; 5.9 Relativistic Energy
  - https://openstax.org/books/university-physics-volume-3/pages/5-3-time-dilation (loaded; note: the section does not mention GPS — it uses muon and Hafele-Keating examples)
- Vol. 1 Ch. 13 "Gravitation": 13.7 Einstein's Theory of Gravity (qualitative gravitational time dilation)
- GPS worked example (MIT OCW 8.033): Session 23 "General Relativity and Black Holes", topics "Gravitational Redshift" and "Application to the GPS System" (slides; the file is served as resource lecture21_bh1, pdf 281 kB, despite the session number).
  - https://ocw.mit.edu/courses/8-033-relativity-fall-2006/pages/lecture-notes/ (loaded 2026-10-16); https://ocw.mit.edu/courses/8-033-relativity-fall-2006/resources/lecture21_bh1/ (loaded)
- Supplement (Ashby 2003, free to read on PubMed Central; see license caveat in the manifest): section "Relativistic Effects on Satellite Clocks" (run-in headings "Satellite orbits", "The eccentricity correction") derives the pre-launch factory frequency offset that combines the velocity and gravitational terms; also "Reference Frames and the Sagnac Effect" and "GPS Coordinate Time and TAI".
  - https://pmc.ncbi.nlm.nih.gov/articles/PMC5253894/ (loaded)
- Extract: time dilation factor for ~3.9 km/s satellite speed; the net GPS clock offset (~38 us/day) assembled from OpenStax 5.3 (velocity term) and 13.7 (gravitational term), checked against the 8.033 slides and Ashby; the lesson still writes the worked example in its own words, since OpenStax does not do it and the two free sources present it as slides and as a research review.

## 04-14 Atomic & Radiation Physics (Vol. 3)

Sources: `@openstax-university-physics-vol3`, `@nasa-why-space-radiation-matters`, `@ntrs-pellish-radiation-effects-tutorial`, `@doe-what-is-radioisotope-power-system`, `@nasa-basics-of-space-flight`, `@nasa-hdbk-4002b-in-space-charging`

- Ch. 8 "Atomic Structure": 8.1 The Hydrogen Atom; 8.4 The Exclusion Principle and the Periodic Table; 8.5 Atomic Spectra and X-rays; 8.6 Lasers
- Ch. 10 "Nuclear Physics": 10.1 Properties of Nuclei; 10.2 Nuclear Binding Energy; 10.3 Radioactive Decay; 10.4 Nuclear Reactions; 10.7 Medical Applications and Biological Effects of Nuclear Radiation
  - https://openstax.org/books/university-physics-volume-3/pages/10-3-radioactive-decay (pattern)
- Space radiation environment (NASA "Why Space Radiation Matters"): "What is Space Radiation?"; "How much Space Radiation are Astronauts Exposed to?"; "Where Does Radiation Come From?"; "What is Non-Ionizing versus Ionizing Radiation?"; "What Is Galactic Cosmic Radiation?"; "Are We Protected from Space Radiation on Earth?"; "What Factors Determine the Amount of Radiation Astronauts Receive?".
  - https://www.nasa.gov/missions/analog-field-testing/why-space-radiation-matters/ (loaded 2026-10-16)
- Effects on electronics (Pellish, NTRS 20170004674, 22 slides): "What are radiation effects?"; "What is total ionizing dose? (TID)"; "What is displacement damage? (DDD)"; "What is NIEL?"; "What are single-event effects?" (SEE taxonomy incl. SEU and SEL); "Single-event effects processes"; "Short history of single-event effects"; "Electronics for Space Use".
  - https://ntrs.nasa.gov/citations/20170004674 (loaded; PDF at /api/citations/20170004674/downloads/20170004674.pdf loaded)
- RTGs (DOE): "Types of Radioisotope Systems"; "How Does an RTG Work?"; "Developing RTGs"; "History of RTGs". NASA Basics of Space Flight Chapter 11 "Onboard Systems" page 3: "Electrical Power Supply and Distribution Subsystems"; "Photovoltaics"; "Radioisotope Power Systems"; "Electrical Power Storage".
  - https://www.energy.gov/ne/articles/what-radioisotope-power-system (loaded); https://science.nasa.gov/learn/basics-of-space-flight/chapter11-3/ (loaded)
- Supplement: NASA-HDBK-4002B "Mitigating In-Space Charging Effects - A Guideline" (2022-06-07) for surface/internal charging; record page only, PDF not loaded.
  - https://standards.nasa.gov/standard/NASA/NASA-HDBK-4002 (loaded)
- Extract: spectra and lasers (lidar, spectroscopy); half-life and activity; dose units from 10.7; trapped-belt, solar-particle and galactic-cosmic-ray sources; TID vs. displacement damage vs. single-event effects and why an SEU is a bit flip while SEL is a current runaway; an RTG as Pu-238 decay heat across a thermoelectric junction, with the half-life budget from 10.3.
- Remaining partial GAP: quantitative shielding (dose-depth curves, aluminium-equivalent thickness) has no verified free source with a titled section; the lesson supplies a first-order treatment.

## 04-15 Electronics Bench Lab (breadboard, multimeter, microcontroller)

Sources: `@openstax-university-physics-vol2`, `@sparkfun-learn-tutorials`, `@micropython-docs`

- Vol. 2 10.4 Electrical Measuring Instruments (how voltmeters/ammeters/ohmmeters work) and 10.6 Household Wiring and Electrical Safety are the only relevant textbook sections.
  - https://openstax.org/books/university-physics-volume-2/pages/10-4-electrical-measuring-instruments (pattern)
- Breadboarding (SparkFun "How to Use a Breadboard"): "Why Use Breadboards?"; "Anatomy of a Breadboard" (Terminal Strips; Power Rails; DIP Support; Rows and Columns); "Providing Power to a Breadboard"; "Building Your First Breadboard Circuit" (What You'll Need; Build the Circuit; Circuit Schematics).
  - https://learn.sparkfun.com/tutorials/how-to-use-a-breadboard/all (loaded 2026-10-16)
- Multimeter practice (SparkFun "How to Use a Multimeter"): "Parts of a Multimeter"; "Probe Types"; "Measuring Voltage" (Overload; Selection Knob); "Measuring Resistance"; "Measuring Current"; "Continuity"; "Changing the Fuse".
  - https://learn.sparkfun.com/tutorials/how-to-use-a-multimeter/all (loaded)
- Microcontroller ADC reading (MicroPython): "Quick reference for the RP2" sections "Installing MicroPython", "Pins and GPIO", "ADC (analog to digital conversion)"; "class ADC - analog to digital conversion": Constructors (machine.ADC(id, *, sample_ns, atten)); Methods (ADC.read_u16(); ADC.read_uv()). SparkFun "Analog to Digital Conversion": "Relating ADC Value to Voltage"; "Arduino ADC Example".
  - https://docs.micropython.org/en/latest/rp2/quickref.html (loaded); https://docs.micropython.org/en/latest/library/machine.ADC.html (loaded); https://learn.sparkfun.com/tutorials/analog-to-digital-conversion/all (loaded)
- Extract: breadboard row/rail topology; voltage in parallel, current in series, resistance with the circuit unpowered; continuity for debugging; reading a potentiometer or sensor with read_u16() and converting counts to volts; the vendor purchase sections of the SparkFun pages are not part of the reading.
- Remaining note: the Arduino documentation (CC BY-SA 4.0 per its repository README) could not be verified page-by-page on 2026-10-16 and is not cited; the lesson is written so that either an Arduino-compatible board or a MicroPython board works.


## Appendix: Leishman, Introduction to Aerospace Flight Vehicles — relevant chapters

Site menu (loaded 2026-10-11) lists the book in parts; chapter titles relevant to the lessons above:
"Mathematics for Engineering"; "Units & Conversion Factors"; "Fundamental Properties of Fluids";
"Fluid Statics & the Hydrostatic Equation"; "Atmospheric Properties & the ISA"; "Momentum Equation";
"Energy Equation & Bernoulli's Equation"; "Airfoil Shapes"; "Aerodynamics of Airfoil Sections";
"Airplane Equations of Motion"; "Fundamentals of Propulsion Systems"; "Airspeed Definitions & Measurement";
"Altitude Definitions & Measurement"; "Rocket Engines"; "Rockets & Launch Vehicles";
"Astronautics, Space & Astrodynamics"; "Electrically-Powered Aircraft".

Chapter pages confirmed to load, with the number shown in the live heading on 2026-10-11:
- 1 What is Aerospace Engineering? — .../chapter/what-is-aerospace-engineering/
- 5 Units & Conversion Factors — .../chapter/units-conversion-factors/
- 13 Fundamental Properties of Fluids — .../chapter/fundamental-properties-of-fluids/
- 13 Atmospheric Properties & the ISA — .../chapter/international-standard-atmosphere-isa/
- 32 Airfoil Shapes — .../chapter/airfoil-geometries/ (cites Abbott & von Doenhoff and the NACA families — pairs with the NTRS report below)
- 37 Airplane Equations of Motion — .../chapter/airplane-equations-of-motion/
- 41 Airspeed Definitions & Measurement — .../chapter/determination-of-airspeed/
- 47 Rocket Engines — .../chapter/rocket-engines/
(prefix: https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/)

Caveat: chapter numbers on this site are not stable. Two menu snapshots fetched the same day gave
different numbering (e.g. "Atmospheric Properties" as 12 or 13, "Rocket Engines" as 45 or 47), and two
chapters both showed "13" in their headings, so cite Leishman chapters by title and URL, never by number.
URL slugs also differ from titles (airfoil-geometries, determination-of-airspeed), so guessing slugs fails.

License as shown on the site (front matter, verbatim): "Introduction to Aerospace Flight Vehicles
Copyright © 2022–2025 by J. Gordon Leishman is licensed under a Creative Commons
Attribution-NonCommercial-NoDerivatives 4.0 International License, except where otherwise noted."
Publisher: Embry-Riddle Aeronautical University (Eagle Pubs), Daytona Beach, FL. Ebook ISBN 979-8-9852614-0-0.
DOI https://doi.org/10.15394/eaglepub.2022.1066. Published 2023-01-01; "last updated" shown as 2025-01-01.
Because it is NoDerivatives, link to it; do not excerpt or adapt figures into the curriculum.

---

### 00-11 — "Reading a Technical Paper": NTRS pick

- Title: Summary of Airfoil Data
- Authors: Ira H. Abbott, Albert E. von Doenhoff, Louis S. Stivers, Jr. (NTRS citation metadata lists only Stivers; the report cover lists all three)
- Report number: NACA-TR-824 (NACA Report No. 824)
- Year: 1945
- NTRS URL: https://ntrs.nasa.gov/citations/19930090976 (loaded 2026-10-11; abstract, authors, and PDF download link present; the inline preview shows "No Preview Available", so the learner downloads the PDF)
- Why this one: public-domain US Government work; classic, heavily cited; the first ~30 pages are readable prose on how NACA airfoil families are designated and tested, followed by data plots the learner can practice reading (Cl vs. alpha, drag polars). It pairs directly with Leishman's "Airfoil Shapes" chapter, which cites it.
- Annotation targets: abstract vs. summary; what "Reynolds number of 3, 6, 9 million" means in the plots; how the 4-digit designation encodes camber and thickness; locate one figure and restate its axes in words.

---
