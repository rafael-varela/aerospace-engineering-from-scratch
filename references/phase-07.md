# Reference map: Phase 07: Thermodynamics

**Note 2026-10-17:** first mapping of Phase 07. Pages were loaded with curl through the session proxy because WebFetch could not resolve hosts in this environment. Every URL marked (loaded) below was downloaded again on 2026-10-17 into a dedicated folder and returned HTTP 200 with the expected page title, and the license statements recorded in sources.yaml were re-read from the downloaded pages. Verified 2026-10-17 by loading the live pages of NTRS (records and PDFs of NASA/TP-2002-211556 and NASA RP-1311 Part I; the NTRS record of RP-1311 Part II), the NASA CEA repository on GitHub (README, LICENSE.txt and `data/thermo.inp` via raw.githubusercontent.com), LibreTexts (Yan *Introduction to Engineering Thermodynamics*: book, front-matter, chapter, section and appendix pages; Ellgen *Thermodynamics and Chemical Equilibrium*: book and chapter pages), ahtt.mit.edu and the full PDF of Lienhard & Lienhard *A Heat Transfer Textbook* v6.00, web.mit.edu (16.Unified Thermodynamics and Propulsion notes: contents page, front-matter nodes and the section nodes cited below; the Fall thermodynamics and Spring propulsion home pages), MIT OCW 16.050 (lecture-notes page, resource pages, and the contents, Part 1C, Part 2A and Part 3 PDFs), and OpenStax University Physics Vol. 2 Ch. 1–4 and Chemistry 2e §§4.1–4.4, 5.1–5.3, 9.3, 13.1–13.4, 16.1–16.4 section pages. Eight ids were added to the manifest; `@openstax-university-physics-vol2` and `@openstax-chemistry-2e` were extended. Not cited: the NIST Chemistry WebBook (webbook.nist.gov unreachable from the verification environment; it is named in prose only, as the original source of the Yan appendix tables) and IAPWS (iapws.org unreachable); the BCcampus Pressbooks original of Yan (pressbooks.bccampus.ca unreachable; the LibreTexts copy is cited); MIT 2.005 and 2.006 Thermal-Fluids Engineering I/II (no OCW course exists: the guessed URLs return 404 and neither course appears in the OCW sitemap); MIT OCW 10.40 Chemical Engineering Thermodynamics, Fall 2003 (loaded: its materials are outline "Supplementary Notes" keyed to the paid Tester & Modell text, partly scanned, with no Gibbs-minimization or exergy derivation a learner could follow alone); MIT OCW 2.051 Introduction to Heat Transfer, Fall 2015 (loaded: readings are only from a paid textbook); MIT OCW 16.01 Unified Engineering Fall 2005 "Thermodynamics and Propulsion" page (loaded: problem sets, quizzes and learning objectives only); NASA spacecraft thermal-control handbooks (none located and loaded, so none cited).

Verified 2026-10-17 by loading the live pages; section numbers and titles are copied from the publishers' own tables of contents, section pages or PDF outlines. Every source here is **required-free** unless marked otherwise; sources are identified by an `@`-prefixed identifier and recorded in [sources.yaml](sources.yaml) with license, version, and caveats. **GAP** marks a lesson that no listed free source covers adequately: the lesson must carry original instruction there, and any "candidate" named is unverified until it is added to the manifest.

Conventions: Yan section URLs follow `https://eng.libretexts.org/Bookshelves/Mechanical_Engineering/Introduction_to_Engineering_Thermodynamics_(Yan)/<NN>%3A_<Chapter_Title>/<N>.<MM>%3A_<Section_Title>`; a URL written as `Yan .../<NN>%3A_...` continues that prefix. Yan's subsection numbers lag their section numbers by one (section 6.8 holds subsections 6.7.1–6.7.4), as LibreTexts shows them; cite by section number and title. The 16.Unified notes are HTML pages `https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/nodeNN.html` (the Spring propulsion URL serves the same notes); titles are from the notes' contents page, and where a page heading differs it is noted. MIT OCW 16.050 resources follow `https://ocw.mit.edu/courses/16-050-thermal-energy-fall-2002/resources/<code>/` (codes `02_part0`, `03_part1a_to_cms`, `04_part1b`, `05_part1c`, `06_part1d`, `07_part2a`, `08_part2b`, `09_part2c`, `10_part3`), each with one PDF. Lienhard's subsections are unnumbered: cite as "§N.M, subsection title"; the book is one PDF, https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf (loaded). OpenStax pages follow `https://openstax.org/books/<book>/pages/<section-slug>`. Licenses: LibreTexts Yan (book and Chapters 1–6), MIT OCW 16.050 and OpenStax are **CC BY-NC-SA 4.0**; the Yan Chapter 7 appendix pages say "shared under a not declared license" (their values are NIST SRD 69 data, used as data with a citation); Ellgen is **CC BY-SA 4.0**; NASA reports are US Government works ("Public Use Permitted"); the NASA CEA repository is **Apache-2.0**; Lienhard is copyrighted, free to download for personal or non-profit instructional use; the 16.Unified notes are all rights reserved ("No portion of this web site may be copied, borrowed or used without written consent from the authors") and are linked by section title only. The curriculum links to these sources and writes its own explanations; nothing is copied or adapted, and property values are used as data with a citation.


## 07-01 Systems, Properties & State

Sources: `@nasa-tp-2002-211556-glenn-coefficients`, `@nasa-cea-github`, `@libretexts-yan-engineering-thermodynamics`, `@mit-16-unified-thermodynamics-propulsion-notes`, `@openstax-university-physics-vol2`, `@lienhard-heat-transfer-textbook`

- NASA/TP-2002-211556: Introduction — "Standard States, Reference States, and Fundamental Constants"; "Empirical Equations for Fitting Thermodynamic Functions" (eqs. for C°p/R, H°/RT, S°/R); "Assigned Enthalpy Values"; "Heats of Formation and Equilibrium Constants"; "Thermodynamic Coefficients Data" → "Temperature Ranges"; "Using the Coefficients to Obtain Thermodynamic Functions"; Appendix B "Thermodynamic Properties at 0 and 298.15 K" (Table B1.—Thermodynamic Properties for Gases at 0 and 298.15 K); Appendix C "Format for Thermodynamic Data Coefficients" (Table C1.—FORTRAN Format Used for Data in Appendix D); Appendix D "Listing of NASA Glenn Thermodynamic Data Coefficients"
  - https://ntrs.nasa.gov/citations/20020085330 (loaded; PDF loaded)
- NASA CEA repository: `data/thermo.inp`, the machine-readable coefficient file the phase's constants were extracted from
  - https://raw.githubusercontent.com/nasa/cea/main/data/thermo.inp (loaded); https://github.com/nasa/cea (README and LICENSE.txt loaded)
- Yan Ch. 1 "Basic Concepts and Definitions": 1.2 What is thermodynamics about?; 1.3 System and surroundings; 1.4 Extensive and intensive properties; 1.5 State, process, and cycle
  - Yan .../01%3A_Basic_Concepts_and_Definitions (loaded); .../1.03%3A_System_and_surroundings (loaded)
- Yan Ch. 3 "Ideal and Real Gasses": 3.2 Ideal gas and ideal gas equation of state; Yan 4.2 Internal energy in a system (4.1.2 Constant-volume specific heat; 4.1.3 Using Cv to calculate Δu for ideal gases); Yan 5.2 Enthalpy (5.1.2 Constant-pressure specific heat; 5.1.3 Using Cp to calculate Δh for ideal gases); 7.7 Appendix G (Table G1: Properties of various ideal gases (at 300K))
  - Yan .../03%3A_Ideal_and_Real_Gasses/3.02%3A_Ideal_gas_and_ideal_gas_equation_of_state (loaded); .../05%3A_The_First_Law_of_Thermodynamics_for_a_Control_Volume/5.02%3A_Enthalpy (loaded)
- 16.Unified notes 1.2 Definitions and Fundamental Ideas of Thermodynamics (1.2.1 The Continuum Model; 1.2.2 The Concept of a "System"; 1.2.3 The Concept of a "State"; 1.2.4 The Concept of "Equilibrium"; 1.2.5 The Concept of a "Process"; 1.2.6 Quasi-Equilibrium Processes; 1.2.7 Equations of state); 2.4 Specific Heats (2.4.1 Specific Heats of an Ideal Gas)
  - https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node11.html (loaded); .../node18.html (loaded)
- OpenStax Vol. 2 Ch. 2: 2.1 Molecular Model of an Ideal Gas; 2.3 Heat Capacity and Equipartition of Energy (the P04.07 link from molecules to cv)
  - https://openstax.org/books/university-physics-volume-2/pages/2-1-molecular-model-of-an-ideal-gas (loaded); .../2-3-heat-capacity-and-equipartition-of-energy (loaded)
- Air check data: Lienhard Appendix A, Table A.6 "Thermophysical properties of gases at atmospheric pressure (101325 Pa)" (contents: "Properties of gases at 1 atm")
- Extract: system, control volume, extensive/intensive and specific properties; the state postulate (two independent intensive properties fix a simple compressible substance); cp(T), h(T) and s°(T) from the nine-coefficient form, with h including the enthalpy of formation; dry air as a mole-fraction-weighted mixture of species polynomials.
- Partial GAP: the phase's fixed dry-air composition (N2 0.7808, O2 0.2095, Ar 0.0093, CO2 0.0004) is a curriculum convention that no listed source states; say so in the lesson.

## 07-02 The First Law

Sources: `@libretexts-yan-engineering-thermodynamics`, `@mit-16-unified-thermodynamics-propulsion-notes`, `@openstax-university-physics-vol2`, `@lienhard-heat-transfer-textbook`

- Yan Ch. 4 "The First Law of Thermodynamics for Closed Systems": 4.2 Internal energy in a system; 4.3 Heat transfer across a boundary; 4.4 Work (4.3.1 Boundary work; 4.3.2 Polytropic process and its boundary work; 4.3.3 Spring work); 4.5 The first law of thermodynamics for closed systems
  - Yan .../04%3A_The_First_Law_of_Thermodynamics_for_Closed_Systems (loaded); .../4.04%3A_Work (loaded)
- Yan Ch. 5 "The First Law of Thermodynamics for a Control Volume": 5.2 Enthalpy; 5.3 Mass and energy conservation equations in a control volume (5.2.1 Steady flow and transient flow; 5.2.2 Mass conservation equation; 5.2.3 Energy conservation equation); 5.4 Applications of the mass and energy conservation equations in steady flow devices (5.3.2 Throttling valves; 5.3.3 Nozzles and diffusers)
  - Yan .../05%3A_The_First_Law_of_Thermodynamics_for_a_Control_Volume/5.03%3A_Mass_and_energy_conservation_equations_in_a_control_volume (loaded)
- 16.Unified notes Ch. 2 "The First Law of Thermodynamics": 2.1 First Law of Thermodynamics; 2.2 Corollaries of the First Law; 2.3 Example Applications of the First Law; enthalpy (page heading: "...to motivate the use of a property called 'enthalpy'"): 2.3.1 Adiabatic, steady, throttling of a gas; 2.3.2 Quasi-Static Expansion of a Gas; 2.3.3 Transient filling of a tank; 2.3.4 The First Law in Terms of Enthalpy; 2.5 Control volume form of the conservation laws (2.5.1 Conservation of mass; 2.5.2 Conservation of energy; 2.5.3 Stagnation Temperature and Stagnation Enthalpy)
  - https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node17.html (loaded); .../node19.html (loaded)
- OpenStax Vol. 2 Ch. 3: 3.1 Thermodynamic Systems; 3.2 Work, Heat, and Internal Energy; 3.3 First Law of Thermodynamics; 3.4 Thermodynamic Processes; 3.6 Adiabatic Processes for an Ideal Gas
  - https://openstax.org/books/university-physics-volume-2/pages/3-3-first-law-of-thermodynamics (loaded); .../3-6-adiabatic-processes-for-an-ideal-gas (loaded)
- Exit flow: 16.Unified notes 14.3 Rocket Nozzles: Connection of Flow to Geometry (14.3.1 Quasi-one-dimensional compressible flow in a variable area duct — mass flow in terms of area, Mach number and stagnation conditions, evaluated at the throat where M = 1; 14.3.2 Thrust in terms of nozzle geometry)
  - https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node104.html (loaded)
- Condensation check data: Lienhard Appendix A, Table A.5 "Thermophysical properties of saturated vapors (p ≠ 1 atm)", nitrogen rows (saturation pressure 0.3605 MPa at 90 K, 0.7783 MPa at 100 K)
  - https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf (loaded)
- Extract: the unsteady control-volume energy balance d(mu)/dt = −ṁh + Q̇ for a tank with one exit; 16.Unified 2.3.3 (filling) is the mirror image of blowdown; adiabatic blowdown follows p ∝ ρ^γ for the gas left in the tank; the choked mass flow from 14.3.1 with the tank as the stagnation state; the isothermal limit needs heat transfer from the wall.
- Partial GAP: no verified free source couples the choked exit flow to the tank balance or derives the closed-form adiabatic blowdown p(t) and the blowdown time; that derivation is original. Pressure vessels are never built, filled or vented at home in this course.

## 07-03 Properties of Pure Substances

Sources: `@libretexts-yan-engineering-thermodynamics`, `@libretexts-ellgen-thermodynamics-chemical-equilibrium`, `@mit-16-unified-thermodynamics-propulsion-notes`, `@openstax-university-physics-vol2`

- Yan Ch. 2 "Thermodynamic Properties of a Pure Substance": 2.2 Pure substance; 2.3 Thermodynamic properties; 2.4 Phase diagrams; 2.5 Thermodynamic tables
  - Yan .../02%3A_Thermodynamic_Properties_of_a_Pure_Substance (loaded); .../2.05%3A_Thermodynamic_tables (loaded)
- Yan 3.3 Real gas and compressibility factor; 4.2 (4.1.1 Using thermodynamic tables to determine specific internal energy u); 5.2 (5.1.1 Using thermodynamic tables to determine specific enthalpy h); 6.8 (6.7.1 Determining the specific entropy of pure substances by using thermodynamic tables)
- Yan Ch. 7 "Appendices": 7.1 Appendix A - Thermodynamic Properties of Water (Table A1: Saturated Water; Table A2: Superheated Vapor, Water; Table A3: Compressed Liquid Water); 7.2 Appendix B - Thermodynamic Properties of Ammonia (Table B1: Saturated Ammonia; Table B2: Superheated Ammonia); 7.4 Appendix D - Thermodynamic Properties of Carbon Dioxide (Table D1: Saturated CO2; Table D2: Superheated CO2); 7.5 Appendix E - Critical Properties of Selected Fluids (Table E1: Critical Constants of Selected Fluids); 7.6 Appendix F - Triple Point of Selected Substances (Table F1: Triple Point of Selected Fluids). Each table states its source: Lemmon, McLinden and Friend, "Thermophysical Properties of Fluid Systems", NIST Chemistry WebBook, NIST SRD 69.
  - Yan .../07%3A_Appendices/7.01%3A_Appendix_A_-_Thermodynamic_Properties_of_Water (loaded); .../7.02%3A_Appendix_B_-_Thermodynamic_Properties_of_Ammonia (loaded); .../7.04%3A_Appendix_D_-_Thermodynamic_Properties_of_Carbon_Dioxide (loaded); .../7.05%3A_Appendix_E_-_Critical_Properties_of_Selected_Fluids (loaded)
- OpenStax Vol. 2 1.5 Phase Changes (Phase Diagrams; Equilibrium; Phase Change and Latent Heat)
  - https://openstax.org/books/university-physics-volume-2/pages/1-5-phase-changes (loaded)
- 16.Unified notes 8.1 Behavior of Two-Phase Systems; 8.4 The Clausius-Clapeyron Equation (page heading adds "(application of 1st and 2nd laws of thermodynamics)")
  - https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node64.html (loaded)
- Ellgen (supplementary-free) Ch. 12 "Applications of the Thermodynamic Criteria for Change": 12.6 Phase Equilibria - Temperature Dependence of the Boiling Point; 12.8 The Clapeyron Equation; 12.9 The Clausius-Clapeyron Equation
  - https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Thermodynamics_and_Chemical_Equilibrium_(Ellgen)/12%3A_Applications_of_the_Thermodynamic_Criteria_for_Change (loaded)
- Extract: saturation dome, quality x = (v − vf)/vfg and property mixing rules; phase identification from (T, p) or (p, v); linear and bilinear interpolation (P03.03); dp/dT = hfg/(T vfg) as a table consistency check; Z → 1 at low pressure.

## 07-04 The Second Law & Entropy

Sources: `@libretexts-yan-engineering-thermodynamics`, `@mit-16-unified-thermodynamics-propulsion-notes`, `@mit-ocw-16-050-thermal-energy`, `@openstax-university-physics-vol2`

- Yan Ch. 6 "Entropy and the Second Law of Thermodynamics": 6.2 Heat engine; 6.3 Refrigerator and heat pump; 6.4 The second law of thermodynamics- Kelvin-Planck and Clausius statements; 6.5 Carnot cycles; 6.6 Entropy and entropy generation (6.5.1 The inequality of Clausius; 6.5.2 Definition of entropy; 6.5.3 Entropy generation, Sgen; 6.5.4 The second law of thermodynamics expressed in terms of entropy generation, Sgen); 6.7 The second law of thermodynamics for closed systems; 6.9 Applications of the second law of thermodynamics in closed systems; 6.10 The second law of thermodynamics for open systems; 6.11 Applications of the second law of thermodynamics in open systems
  - Yan .../06%3A_Entropy_and_the_Second_Law_of_Thermodynamics (loaded); .../6.06%3A_Entropy_and_entropy_generation (loaded)
- 16.Unified notes Ch. 4 "Background to the Second Law of Thermodynamics": 4.1 Reversibility and Irreversibility in Natural Processes; 4.2 Difference between Free and Isothermal Expansions; 4.3 Features of reversible processes. Ch. 5 "The Second Law of Thermodynamics": 5.1 Concept and Statements of the Second Law; 5.3 Combined First and Second Law Expressions; 5.4 Entropy Changes in an Ideal Gas; 5.5 Calculation of Entropy Change in Some Basic Processes (worked cases: heat transfer from, or to, a heat reservoir; heat transfer between two heat reservoirs; possibility of obtaining work from a single heat reservoir; free versus reversible isothermal expansion; the "hot brick problem")
  - https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node39.html (loaded); .../node41.html (loaded)
- MIT OCW 16.050 Part 1C "Applications of the Second Law": 1.C.7 Examples of Lost Work in Engineering Processes, items c) Isothermal Expansion with Friction and d) Entropy Generation, Irreversibility, and Cycle Efficiency
  - https://ocw.mit.edu/courses/16-050-thermal-energy-fall-2002/resources/05_part1c/ (loaded; PDF loaded)
- OpenStax Vol. 2 Ch. 4: 4.1 Reversible and Irreversible Processes; 4.2 Heat Engines; 4.3 Refrigerators and Heat Pumps; 4.4 Statements of the Second Law of Thermodynamics; 4.5 The Carnot Cycle; 4.6 Entropy
  - https://openstax.org/books/university-physics-volume-2/pages/4-6-entropy (loaded); .../4-5-the-carnot-cycle (loaded)
- Extract: entropy as a state property (path independence via T ds = du + p dv); S_gen ≥ 0 as the bookkeeping form of the second law; the four generators the build tracks — heat across finite ΔT, throttling, mixing, friction; the Clausius integral ∮δQ/T ≤ 0 and the Carnot bound.

## 07-05 Isentropic Processes & Efficiencies

Sources: `@libretexts-yan-engineering-thermodynamics`, `@mit-16-unified-thermodynamics-propulsion-notes`, `@mit-ocw-16-050-thermal-energy`, `@nasa-tp-2002-211556-glenn-coefficients`, `@openstax-university-physics-vol2`, `@lienhard-heat-transfer-textbook`

- Yan 6.8 Specific entropy of a state (6.7.3 Determining the specific entropy of ideal gases; 6.7.4 Isentropic relations for an ideal gas); 4.4 Work (4.3.2 Polytropic process and its boundary work); 5.4 Applications of the mass and energy conservation equations in steady flow devices (its opening subsection on turbines and compressors, whose heading is blank on the LibreTexts page); 6.11 Applications of the second law of thermodynamics in open systems (worked two-stage compressor with intercooler)
  - Yan .../06%3A_Entropy_and_the_Second_Law_of_Thermodynamics/6.08%3A_Specific_entropy_of_a_state (loaded); .../05%3A_The_First_Law_of_Thermodynamics_for_a_Control_Volume/5.04%3A_Applications_of_the_mass_and_energy_conservation_equations_in_steady_flow_devices (loaded)
- 16.Unified notes 2.4 Specific Heats (2.4.2 Reversible adiabatic processes for an ideal gas); 6.7 Examples of Lost Work in Engineering Processes, subsection "Adiabatic Efficiency of a Propulsion System Component (Turbine)"; 12.4 Multistage Axial Compressors (context)
  - https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node50.html (loaded); .../node92.html (loaded)
- MIT OCW 16.050 1.C.7 b) Adiabatic Efficiency of a Propulsion System Component (Turbine)
  - https://ocw.mit.edu/courses/16-050-thermal-energy-fall-2002/resources/05_part1c/ (loaded; PDF loaded)
- OpenStax Vol. 2 3.6 Adiabatic Processes for an Ideal Gas
- Variable cp: NASA/TP-2002-211556 "Empirical Equations for Fitting Thermodynamic Functions" (s°(T)/R, so that s°(T2) − s°(T1) = R ln(p2/p1) defines the isentropic T2)
- Air check data: Lienhard Appendix A, Table A.6 "Thermophysical properties of gases at atmospheric pressure (101325 Pa)" (cp of air at 300–1400 K, against the NASA dry-air model)
  - https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf (loaded)
- Extract: T2s = T1 π^((γ−1)/γ) for constant cp and its variable-cp generalisation through s°(T); isentropic efficiency of a compressor as work_ideal/work_actual; the entropy generated by a non-ideal stage; the constant-cp limit as a test.
- Partial GAP: no verified free source defines polytropic (small-stage) efficiency or derives its relation to the overall isentropic efficiency of a multi-stage compressor (the "preheat" effect); the 16.Unified and 16.050 efficiency subsections treat the turbine only. Polytropic efficiency and the compressor form are original.

## 07-06 Power Cycles: Brayton

Sources: `@mit-16-unified-thermodynamics-propulsion-notes`, `@mit-ocw-16-050-thermal-energy`, `@openstax-university-physics-vol2`, `@nasa-tp-2002-211556-glenn-coefficients`

- 16.Unified notes Ch. 3 "The First Law Applied to Engineering Cycles": 3.1 Some Properties of Engineering Cycles; Work and Efficiency; 3.2 Generalized Representation of Thermodynamic Cycles; 3.3 The Carnot Cycle; 3.7 Brayton Cycle (3.7.1 Work and Efficiency; 3.7.2 Gas Turbine Technology and Thermodynamics — the compressor pressure ratio for maximum work per unit mass; 3.7.3 Brayton Cycle for Jet Propulsion: the Ideal Ramjet)
  - https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node28.html (loaded); .../node22.html (loaded)
- 16.Unified notes Ch. 6 "Applications of the Second Law": 6.1 Limitations on the Work that Can be Supplied by a Heat Engine; 6.3 Representation of Thermodynamic Processes in T-s coordinates; 6.4 Brayton Cycle in T-s Coordinates (6.4.1 Net work per unit mass flow in a Brayton cycle)
  - https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node47.html (loaded)
- MIT OCW 16.050 Part 1A "Background to the Second Law of Thermodynamics": 1.A.3 Brayton Cycles (or Joule Cycles): The Power Cycle for a Gas Turbine Jet Engine; 1.A.4 Gas Turbine Technology and Thermodynamics; Part 1C: 1.C.4 Brayton Cycle in T-s coordinates; Part 2A "Gas Power and Propulsion Cycles": 2.A.3 Brayton Cycle
  - https://ocw.mit.edu/courses/16-050-thermal-energy-fall-2002/resources/03_part1a_to_cms/ (loaded; contents PDF loaded); .../07_part2a/ (loaded; PDF loaded)
- OpenStax Vol. 2 4.2 Heat Engines; 4.5 The Carnot Cycle
- Variable-cp extension: NASA/TP-2002-211556 coefficient functions
- Extract: the four-process state table; η = 1 − π^(−(γ−1)/γ); specific work as the enclosed area in p–v and T–s; the optimum pressure ratio for maximum work at fixed turbine inlet temperature (16.Unified 3.7.2); why η < η_Carnot between the same extreme temperatures (16.Unified 6.4 argues it with small Carnot cycles).

## 07-07 Real Brayton Cycles

Sources: `@mit-ocw-16-050-thermal-energy`, `@mit-16-unified-thermodynamics-propulsion-notes`, `@lienhard-heat-transfer-textbook`, `@libretexts-yan-engineering-thermodynamics`

- MIT OCW 16.050 2.A.7 Effect of Departures from Ideal Behavior - Real Cycle behavior (compressor and turbine losses, combustor stagnation-pressure drop, heat transfer)
  - https://ocw.mit.edu/courses/16-050-thermal-energy-fall-2002/resources/07_part2a/ (loaded; PDF loaded)
- 16.Unified notes 11.6 Performance of Jet Engines (11.6.5 Effect of Departures from Ideal Behavior); 3.7.4 MIT Cogenerator; 18.5 Heat Exchangers (18.5.1 Simplified Counterflow Heat Exchanger (With Uniform Wall Temperature); 18.5.2 General Counterflow Heat Exchanger; 18.5.3 Efficiency of a Counterflow Heat Exchanger)
  - https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node131.html (loaded)
- Recuperator: Lienhard Ch. 3 "Heat exchanger design": 3.1 Function and configuration of heat exchangers; 3.3 Heat exchanger effectiveness
- Intercooling: Yan 6.11 Applications of the second law of thermodynamics in open systems (two-stage compressor with intercooler example)
- Extract: component isentropic efficiencies and pressure-loss factors in the state table; recuperator effectiveness ε = (T_cold,out − T_cold,in)/(T_hot,in − T_cold,in); why recuperation helps only below the pressure ratio where turbine exit temperature equals compressor exit temperature; intercooling cuts compressor work but adds combustor heat.
- Partial GAP: no verified free source analyses the recuperated or intercooled gas-turbine cycle as a cycle (the ideal recuperated-cycle efficiency 1 − (T1/T3)π^((γ−1)/γ), the crossover pressure ratio, or turboprop/APU-scale numbers); that analysis is original, with Lienhard §3.3 for the effectiveness concept.

## 07-08 Exergy

Sources: `@mit-16-unified-thermodynamics-propulsion-notes`, `@mit-ocw-16-050-thermal-energy`, `@nasa-tp-2002-211556-glenn-coefficients`

- 16.Unified notes 6.5 Irreversibility, Entropy Changes, and "Lost Work"; 6.6 Entropy and Unavailable Energy (page heading adds "(Lost Work by Another Name)"); 6.7 Examples of Lost Work in Engineering Processes (subsections "Lost work in Adiabatic Throttling: Entropy and Stagnation Pressure Changes"; "Adiabatic Efficiency of a Propulsion System Component (Turbine)"; "Isothermal Expansion with Friction"; "Propulsive Power and Entropy Flux")
  - https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node48.html (loaded); .../node49.html (loaded); .../node50.html (loaded)
- MIT OCW 16.050 Part 1C: 1.C.5 Irreversibility, Entropy Changes, and Lost Work; 1.C.6 Entropy and Unavailable Energy (Lost Work by Another Name); 1.C.7 Examples of Lost Work in Engineering Processes (a) Lost work in Adiabatic Throttling: Entropy and Stagnation Pressure Changes; b) Adiabatic Efficiency of a Propulsion System Component (Turbine); c) Isothermal Expansion with Friction; d) Entropy Generation, Irreversibility, and Cycle Efficiency; e) Propulsive Power and Entropy Flux)
  - https://ocw.mit.edu/courses/16-050-thermal-energy-fall-2002/resources/05_part1c/ (loaded; PDF loaded)
- Property data: NASA/TP-2002-211556 h(T), s°(T) for the dead state and each station
- Extract: maximum work from heat Q at T against surroundings at T0 is Q(1 − T0/T); lost work = T0·ΔS_total (the Gouy–Stodola result, in 16.050 1.C.6 as (W_rev − W) = T0 ΔS_total); the turbojet-core audit sums T0·S_gen per component.
- Partial GAP: no verified free source defines flow exergy ψ = (h − h0) − T0(s − s0), writes the full control-volume exergy balance, defines fuel (chemical) exergy, or defines a second-law efficiency for an engine; these are original (Ch. 14 of the paid Tester & Modell text, which MIT OCW 10.40 supplements, was not used).

## 07-09 Gas Mixtures & Psychrometrics

Sources: `@openstax-chemistry-2e`, `@mit-16-unified-thermodynamics-propulsion-notes`, `@lienhard-heat-transfer-textbook`, `@libretexts-yan-engineering-thermodynamics`, `@nasa-tp-2002-211556-glenn-coefficients`

- OpenStax Chemistry 2e 9.3 Stoichiometry of Gaseous Substances, Mixtures, and Reactions (The Pressure of a Mixture of Gases: Dalton's Law; Collection of Gases over Water)
  - https://openstax.org/books/chemistry-2e/pages/9-3-stoichiometry-of-gaseous-substances-mixtures-and-reactions (loaded)
- 16.Unified notes 7.1 Entropy Change in Mixing of Two Ideal Gases
  - https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node54.html (loaded)
- Lienhard §11.2 Mixture compositions and species fluxes, subsection "Ideal gases" (Dalton's law of partial pressures, eqn. 11.14); §11.7 Simultaneous heat and mass transfer at low rates, subsections "The wet bulb temperature" and "Dew-point temperature and dehumidification" (relative humidity is defined in the section's evaporative-cooling worked example, a cup of coffee in a 50 % RH room)
- Saturation pressure: Yan 7.1 Appendix A (Table A1: Saturated Water); single-component stream mixing: Yan 5.4 (5.3.5 Mixing chambers)
- Water vapour as an ideal gas: NASA/TP-2002-211556 coefficient functions (H2O)
- Extract: mole fraction = pressure fraction for ideal gases; relative humidity φ = p_v/p_sat(T) and dew point; wet-bulb temperature as an energy balance; entropy of mixing as a source of S_gen.
- Partial GAP: no verified free source defines the humidity ratio ω, moist-air enthalpy per kilogram of dry air, Amagat's law, adiabatic mixing of moist streams as a straight line in (ω, h), or the psychrometric chart, and none describes an airliner environmental control system (bleed air, air-cycle pack, recirculation) at cruise; all are original.

## 07-10 Chemical Reactions & Combustion

Sources: `@mit-16-unified-thermodynamics-propulsion-notes`, `@mit-ocw-16-050-thermal-energy`, `@openstax-chemistry-2e`, `@nasa-tp-2002-211556-glenn-coefficients`, `@nasa-rp-1311-cea-analysis`, `@nasa-cea-github`

- 16.Unified notes Ch. 15 "Generating Heat: Thermochemistry": 15.1 Fuels; 15.2 Fuel-Air Ratio; 15.3 Enthalpy of formation; 15.4 First Law Analysis of Reacting Systems; 15.5 Adiabatic Flame Temperature (15.5.1 Approximate solution using "average" values of specific heat; 15.5.2 Solution for adiabatic flame temperature using evolutions of specific heats with temperature; 15.5.3 Solution for adiabatic flame temperature using tabulated values for gas enthalpy)
  - https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node111.html (loaded); .../node108.html (loaded)
- MIT OCW 16.050 Part 2C "Introduction to Thermochemistry": 2.C.1 Fuels; 2.C.2 Fuel-Air Ratio; 2.C.3 Enthalpy of Formation; 2.C.4 First Law analysis of Reacting systems; 2.C.5 Adiabatic Flame Temperature
  - https://ocw.mit.edu/courses/16-050-thermal-energy-fall-2002/resources/09_part2c/ (loaded; titles from the loaded contents PDF)
- OpenStax Chemistry 2e 4.1 Writing and Balancing Chemical Equations; 4.3 Reaction Stoichiometry; 4.4 Reaction Yields (Limiting Reactant); 5.3 Enthalpy (Standard Enthalpy of Combustion; Standard Enthalpy of Formation; Hess's Law)
  - https://openstax.org/books/chemistry-2e/pages/5-3-enthalpy (loaded); .../4-1-writing-and-balancing-chemical-equations (loaded)
- NASA/TP-2002-211556 "Assigned Enthalpy Values"; Table B1 (ΔfH°(298.15) of gases); Table B3.—Assigned Enthalpies for Species That Are Used in CEA as Reactants Only
- NASA RP-1311 Part I Ch. 1 Introduction (the hp problem gives constant-pressure combustion); Ch. 4 "Thermodynamic Data": 4.1 Assigned Enthalpies; 4.2 Least-Squares Coefficients; Ch. 9 "Input Calculations" (equivalence ratios)
  - https://ntrs.nasa.gov/citations/19950013764 (loaded; PDF loaded)
- Extract: element balances for CxHy + air at any equivalence ratio; h_reactants(T_in) = h_products(T_ad) solved by root finding (P03.02) with NASA polynomials; why T_ad peaks near stoichiometric; 16.Unified 15.5 notes that dissociation lowers the real flame temperature (lesson 07-11).
- Partial GAP: no verified free source states the higher/lower heating value distinction (product water liquid vs. vapour) with numbers, or a Jet-A (C12H23) surrogate's properties beyond the NASA database entry; both are original, computed from the NASA enthalpies. Computational only — no combustion of any kind at home.

## 07-11 Chemical Equilibrium

Sources: `@nasa-rp-1311-cea-analysis`, `@libretexts-ellgen-thermodynamics-chemical-equilibrium`, `@openstax-chemistry-2e`, `@nasa-tp-2002-211556-glenn-coefficients`, `@nasa-cea-github`

- NASA RP-1311 Part I Ch. 2 "Equations Describing Chemical Equilibrium": 2.1 Units; 2.2 Equation of State; 2.3 Minimization of Gibbs Energy (2.3.1 Gibbs Iteration Equations; 2.3.2 Reduced Gibbs Iteration Equations); Ch. 3 "Procedure for Obtaining Equilibrium Compositions": 3.1 Initial Estimates; 3.2 Magnitude of Species Used During Iteration; 3.3 Convergence; 3.6 Singularities; Ch. 6 "Theoretical Rocket Performance": 6.3.1 Combustion Conditions (chamber equilibrium context)
  - https://ntrs.nasa.gov/citations/19950013764 (loaded; PDF loaded)
- Ellgen (supplementary-free) Ch. 13 "Equilibria in Reactions of Ideal Gases": 13.1 The Gibbs Free Energy of an Ideal Gas; 13.2 The Gibbs Free Energy Change for A Reaction of Ideal Gases; 13.3 The Thermodynamics of Mixing Ideal Gases; 13.5 ∆rG is the rate at which the Gibbs Free Energy Changes with The Extent of Reaction; 13.6 The Standard Gibbs Free Energy Change and Equilibrium in Ideal Gas Reactions; 13.7 The Gibbs Free Energy of Formation and Equilibrium in Ideal Gas Reactions
  - https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Thermodynamics_and_Chemical_Equilibrium_(Ellgen)/13%3A_Equilibria_in_Reactions_of_Ideal_Gases (loaded)
- OpenStax Chemistry 2e 13.1 Chemical Equilibria; 13.2 Equilibrium Constants (Homogeneous Equilibria; Heterogeneous Equilibria; Coupled Equilibria); 13.3 Shifting Equilibria: Le Châtelier's Principle; 13.4 Equilibrium Calculations; 16.4 Free Energy (Free Energy and Equilibrium; Temperature Dependence of Spontaneity)
  - https://openstax.org/books/chemistry-2e/pages/13-2-equilibrium-constants (loaded); .../16-4-free-energy (loaded)
- g°(T) = h(T) − T s°(T) from NASA/TP-2002-211556; its "Heats of Formation and Equilibrium Constants" explains how the database relates to K
- Reference implementation (not to be called from tested code): NASA CEA repository (README "samples/rp1311_examples")
- Extract: minimise G = Σ n_j μ_j subject to element conservation; Lagrange multipliers as element potentials (RP-1311 §2.3); K_p = exp(−ΔG°/RT) as the independent check (Ellgen 13.6); why high-temperature equilibrium (dissociation to CO, H2, OH, H, O) gives a lower T_ad than complete combustion.
- Partial GAP: RP-1311 is a research-level derivation and Ellgen works one reaction at a time; no verified free source builds the element-potential method step by step for a learner or compares it with a general-purpose constrained minimizer (P03.15). That bridge is original. Computational only.

## 07-12 Heat Transfer: Conduction

Sources: `@lienhard-heat-transfer-textbook`, `@mit-16-unified-thermodynamics-propulsion-notes`, `@mit-ocw-16-050-thermal-energy`, `@openstax-university-physics-vol2`

- Lienhard Ch. 2 "Heat conduction concepts, thermal resistance, and the overall heat transfer coefficient": 2.1 The heat conduction equation (subsection "Fourier's law"); 2.2 Steady heat conduction in a slab: method; 2.3 Thermal resistance and the electrical analogy (Resistances for cylinders and for convection; Resistance for thermal radiation; Contact resistance)
- Lienhard Ch. 4 "Conduction analysis, dimensional analysis, and fin design": 4.1 The well-posed problem; 4.2 General solution of the heat conduction equation (One-dimensional steady heat conduction; Separation of variables: A general solution of multidimensional problems — Example 4.1 is the slab with a sinusoidal temperature on one edge, the build's Fourier-series check); 4.3 Dimensional analysis; 4.5 Fin design (Analysis of a one-dimensional fin; Fin performance and design; Fin thermal resistance; Fin Arrays; Fins of variable cross section; Wall thermal resistance and fin root temperature)
- Lienhard Ch. 5 "Transient and multidimensional heat conduction": 5.6 Transient heat conduction to a semi-infinite region (Table 5.3 Error function and complementary error function); 5.7 Steady multidimensional heat conduction (The flux plot; The shape factor; Table 5.4 Conduction shape factors); 5.8 Transient multidimensional heat conduction—the tactic of superposition
- Properties: Lienhard Table A.1 Properties of metallic solids; Table A.2 Properties of nonmetallic solids
- 16.Unified notes Ch. 16 "Conductive Heat Transfer": 16.1 Heat Transfer Modes; 16.2 Introduction to Conduction; 16.3 Steady-State One-Dimensional Conduction; 16.4 Thermal Resistance Circuits; 16.5 Quasi-One-Dimensional Heat Flow; Ch. 18: 18.1 Temperature Distributions in the Presence of Heat Sources; 18.2 Heat Transfer From a Fin
  - https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node115.html (loaded); .../node128.html (loaded)
- MIT OCW 16.050 Part 3 "Introduction to Engineering Heat Transfer": 2.0 Conduction Heat Transfer; 2.1 Steady-State One-Dimensional Conduction; 2.2 Thermal Resistance Circuits; 2.3 Steady Quasi-One-Dimensional Heat Flow in Non-Planar Geometry; 4.0 Temperature Distributions in the Presence of Heat Sources; 5.0 Heat Transfer From a Fin
  - https://ocw.mit.edu/courses/16-050-thermal-energy-fall-2002/resources/10_part3/ (loaded; PDF loaded)
- OpenStax Vol. 2 1.6 Mechanisms of Heat Transfer (Conduction)
  - https://openstax.org/books/university-physics-volume-2/pages/1-6-mechanisms-of-heat-transfer (loaded)
- Extract: Fourier's law and the heat equation; thermal-resistance networks for a nozzle wall or tile stack; fin efficiency; the separated-variables rectangle solution and the erfc semi-infinite solution as checks on the numerical solver.
- Partial GAP: none of these sources discretises the heat equation (Lienhard's 6th edition treats 2-D steady conduction analytically and by flux plots); the finite-difference Laplace/Poisson solver and the explicit/implicit transient schemes build on P03.11 and are original, as is the nozzle-wall/heat-shield context.

## 07-13 Heat Transfer: Convection

Sources: `@lienhard-heat-transfer-textbook`, `@mit-16-unified-thermodynamics-propulsion-notes`, `@mit-ocw-16-050-thermal-energy`, `@openstax-university-physics-vol2`

- Lienhard Ch. 6 "Laminar and turbulent boundary layers": 6.1 Some introductory ideas; 6.2 Laminar incompressible boundary layer on a flat surface; 6.3 The energy equation; 6.4 The Prandtl number and the boundary layer thicknesses; 6.5 Heat transfer coefficient for laminar, incompressible flow over a flat surface (Predicting the temperature distribution in the laminar thermal boundary layer; Predicting the heat flux in the laminar boundary layer; The problem of uniform wall heat flux); 6.6 The Reynolds-Colburn analogy; 6.7 Turbulent boundary layers; 6.8 Heat transfer in turbulent boundary layers (The Reynolds-Colburn analogy for turbulent flow; Power-law correlations for local turbulent heat transfer); 6.9 Turbulent transition and overall heat transfer (A correlation for laminar, transitional, and turbulent flow; Average Nusselt number for the entire plate)
- Lienhard Ch. 7 "Forced convection in a variety of configurations": 7.2 Heat transfer to or from laminar flows in pipes; 7.3 Turbulent pipe flow; 7.5 Heat transfer coefficients for noncircular ducts; 7.6 Heat transfer during cross flow over cylinders; 7.7 Finding and assessing correlations for other configurations
- Lienhard Ch. 8 "Natural convection in single-phase fluids and during film condensation": 8.3 Laminar natural convection on a vertical isothermal surface; 8.4 Natural convection in other situations (Natural convection from horizontal isothermal cylinders; Natural convection from vertical cylinders; Heat transfer from general submerged bodies; Laminar heat transfer from inclined and horizontal plates; Natural convection with uniform heat flux)
- Properties: Lienhard Table A.6 Thermophysical properties of gases at atmospheric pressure (101325 Pa); Table A.3 Thermophysical properties of saturated liquids
- 16.Unified notes Ch. 17 "Convective Heat Transfer": 17.1 The Reynolds Analogy; 17.2 Combined Conduction and Convection; 17.3 Dimensionless Numbers and Analysis of Results
  - https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node122.html (loaded); .../node124.html (loaded)
- MIT OCW 16.050 Part 3: 3.0 Convective Heat Transfer; 3.1 The Reynolds Analogy; 3.2 Combined Conduction and Convection; 3.3 Dimensionless Numbers and Analysis of Results
- OpenStax Vol. 2 1.6 Mechanisms of Heat Transfer (Convection)
- Extract: Re, Pr, Nu, Gr, Ra and St; local vs. average Nusselt number on a flat plate (average = 2 × local at the trailing edge for laminar flow); St·Pr^(2/3) = c_f/2; each correlation's stated range of validity; property evaluation at the film temperature by interpolating Table A.6.
- Avionics-box, cold-plate and regenerative-cooling-channel examples are original applications of the cited correlations.

## 07-14 Heat Transfer: Radiation

Sources: `@lienhard-heat-transfer-textbook`, `@mit-16-unified-thermodynamics-propulsion-notes`, `@mit-ocw-16-050-thermal-energy`, `@openstax-university-physics-vol2`

- Lienhard §1.3 Modes of heat transfer, subsection "Radiation" (radiant heat exchange, a first view factor, radiation shielding)
- Lienhard Ch. 10 "Radiative heat transfer": 10.1 The problem of radiative exchange (Some definitions; Black body; Table 10.1 Total emittances for a variety of surfaces); 10.2 Kirchhoff's law (Total absorptance during radiant exchange; The gray body approximation); 10.3 Radiant heat exchange between two finite black bodies (View factor reciprocity; Calculation of the black-body view factor, F1–2; Table 10.2 View factors for a variety of two-dimensional configurations (infinite in extent normal to the paper); Table 10.3 View factors for some three-dimensional configurations); 10.4 Heat transfer among gray bodies (Additional two-body exchange problems; The electrical circuit analogy when more than two gray bodies are involved in heat exchange; Algebraic solution of multisurface enclosure problems; Enclosures with nonisothermal, nongray, or nondiffuse surfaces); 10.6 Solar energy (The distribution of the Sun's energy and atmospheric irradiation; Selective emitters, absorbers, and transmitters; Table 10.5 Solar absorptance and infrared emittance for several surfaces near 300 K)
- Lienhard 2.3, subsection "Resistance for thermal radiation" (linearised radiation coefficient)
- 16.Unified notes Ch. 19 "Radiation Heat Transfer": 19.1 Ideal Radiators; 19.2 Kirchhoff's Law and "Real Bodies"; 19.3 Radiation Heat Transfer Between Planar Surfaces (19.3.1 Example 1: Use of a thermos bottle to reduce heat transfer; 19.3.2 Example 2: Temperature measurement error due to radiation heat transfer); 19.4 Radiation Heat Transfer Between Arbitrary Surfaces (page heading: "...Between Black Surfaces of Arbitrary Geometry"; 19.4.1 Example: Concentric cylinders or concentric spheres)
  - https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node136.html (loaded); .../node137.html (loaded)
- MIT OCW 16.050 Part 3: 9.0 Radiation Heat Transfer (Heat transfer by thermal radiation); 9.1 Ideal Radiators; 9.2 Kirchhoff's Law and "Real Bodies"; 9.3 Radiation Heat Transfer Between Planar Surfaces; 9.4 Radiation Heat Transfer Between Black Surfaces of Arbitrary Geometry (its view-factor charts are third-party figures from a commercial text)
- OpenStax Vol. 2 1.6 Mechanisms of Heat Transfer (Radiation)
- Extract: the view-factor double-area integral; reciprocity A_i F_ij = A_j F_ji and summation Σ_j F_ij = 1; the radiosity network and its linear system (P01.13); solar absorptance vs. infrared emittance for a radiator coating (Table 10.5).
- Partial GAP: no verified free source treats numerical evaluation of the view-factor double integral or a spacecraft radiator's equilibrium temperature with solar, albedo and planetary loads; both are original (Lienhard 10.6 and Table 10.5 supply the surface properties and solar input).

## 07-15 Home Lab: Cooling Curve

Sources: `@lienhard-heat-transfer-textbook`, `@mit-16-unified-thermodynamics-propulsion-notes`, `@mit-ocw-16-050-thermal-energy`, `@openstax-university-physics-vol2`

- Lumped model: Lienhard §1.3, subsection "Heat Convection" (run-in "Lumped-capacity solution" and the Biot-number condition for using it); 5.2 Lumped-capacity solutions (Electrical and mechanical analogies to the lumped-thermal-capacity problem; Lumped-capacity solution with a variable ambient temperature; Second-order lumped-capacity systems)
- Distributed model: Lienhard 5.3 Transient conduction in a one-dimensional slab (Table 5.1 Terms of series solutions for slabs, cylinders, and spheres); 5.4 Temperature-response charts; 5.5 One-term solutions (Table 5.2 One-term coefficients for transient conduction with surface convection)
- Prediction of h: Lienhard 8.3 Laminar natural convection on a vertical isothermal surface; 8.4 Natural convection in other situations (Natural convection from horizontal isothermal cylinders; Laminar heat transfer from inclined and horizontal plates); 2.3 "Resistance for thermal radiation"; Table 10.1 Total emittances for a variety of surfaces; properties from Table A.1 Properties of metallic solids and Table A.6 Thermophysical properties of gases at atmospheric pressure (101325 Pa)
- 16.Unified notes 18.3 Transient Heat Transfer (Convective Cooling or Heating); 17.3 Dimensionless Numbers and Analysis of Results (Biot number)
  - https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node129.html (loaded)
- MIT OCW 16.050 Part 3: 6.0 Transient Heat Transfer (Convective Cooling or Heating)
- OpenStax Vol. 2 1.4 Heat Transfer, Specific Heat, and Calorimetry; 1.6 Mechanisms of Heat Transfer
  - https://openstax.org/books/university-physics-volume-2/pages/1-4-heat-transfer-specific-heat-and-calorimetry (loaded)
- Extract: τ = ρcV/(hA) and the linear fit of ln θ vs. t; Bi = hL/k and when the lumped model holds; the one-term series solution as the distributed comparison; h_total = h_conv + h_rad with h_rad ≈ 4εσT_m³; uncertainty propagation from P00.14.
- GAP: the home procedure (60 °C tap-hot water limit, tongs, insulating stand, kitchen probe or the P04.15 thermistor, safety box, data logging) and the synthetic dataset are original; no free source describes this experiment.


## Optional textbooks for deeper study

Paid; never required. Not recorded in `sources.yaml`.

- Y. A. Çengel and M. A. Boles, *Thermodynamics: An Engineering Approach* (McGraw-Hill) — psychrometrics, exergy and gas-turbine cycles with regeneration and intercooling, with extensive property tables (07-01 … 07-10).
- M. J. Moran, H. N. Shapiro, D. D. Boettner and M. B. Bailey, *Fundamentals of Engineering Thermodynamics* (Wiley) — the standard treatment of control-volume exergy balances and second-law efficiency (07-08) and of psychrometric processes (07-09).
- T. L. Bergman, A. S. Lavine, F. P. Incropera and D. P. DeWitt, *Fundamentals of Heat and Mass Transfer* (Wiley) — finite-difference conduction, correlations and view-factor charts (07-12 … 07-15).
- S. R. Turns, *An Introduction to Combustion: Concepts and Applications* (McGraw-Hill) — heating values, adiabatic flame temperature and equilibrium with dissociation (07-10, 07-11).
- H. I. H. Saravanamuttoo, G. F. C. Rogers, H. Cohen, P. V. Straznicky and A. C. Nix, *Gas Turbine Theory* (Pearson) — real cycles, polytropic efficiency and component losses (07-05 … 07-07).
- J. L. Kerrebrock, *Aircraft Engines and Gas Turbines* (MIT Press) — cycle analysis and turbomachinery in the 16.Unified tradition (07-05 … 07-08).
- D. G. Gilmore (ed.), *Spacecraft Thermal Control Handbook, Vol. I: Fundamental Technologies* (The Aerospace Press/AIAA) — radiator sizing, coatings and orbital heat loads (07-14).
