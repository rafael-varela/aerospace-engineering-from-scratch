# Reference map: Phase 08: Fluid Mechanics

**Note 2026-10-17:** first mapping of Phase 08. The WebFetch tool failed DNS resolution for every host except github.com in this session, so all pages except the CFD Python README were loaded live with curl through the session's HTTPS proxy and read as text (PDFs with pdftotext, and page images for scanned PDFs). Verified 2026-10-17 by loading: MIT OCW 2.25 (home page, syllabus, all eleven topic pages, five class-note PDFs), 1.060 (lecture-notes page, five lecture PDFs), 16.100 (lecture-notes page, nine lecture PDFs), 16.01–16.04 Unified Engineering (fluid-mechanics page, F8 and F10 PDFs), 2.27 (lecture-notes page, four outline PDFs) and 2.29 (lecture page, lectures 17–25 PDFs); the NASA Glenn *Beginner's Guide to Aeronautics* (legacy index and pages, new-site pages); LibreTexts (Bar-Meir, Liburdy, Smyth, Ahmari & Kabir, Marchman: chapter pages, and section lists read by following each book's next-page links because the contents panel is script-rendered); Leishman's contents page and nine chapter pages; OpenStax University Physics Vol. 1 Ch. 14 (14.2, 14.6, 14.7); barbagroup/CFDPython (README, LICENSE, Step 10–12 notebooks); NTRS API records and files for the 1976 U.S. Standard Atmosphere, NACA TR 1235, NACA TR 1135 and Malik (1987); JHTDB (home, datasets, Citing); the Lee & Moser channel-data README and data pages; Schlatter's KTH data page; the NPARC verification tutorial; Crossref and doi.org for Ghia, Ghia & Shin (1982). Twenty-one ids were added; `@leishman-intro-aerospace-flight-vehicles`, `@openstax-university-physics-vol1` and `@nasa-nparc-cfd-vv-tutorial` were extended. Leishman's site now shows a 2027 edition and chapters are renumbered; it is cited by chapter title only.

Not cited: the MIT NCFMF film series ("Illustrated Experiments in Fluid Mechanics") — web.mit.edu/hml/ncfmf.html, web.mit.edu/hml/ and the film-notes PDFs linked from 2.25 (e.g. .../ncfmf/09VOR.pdf) all return MIT's 404 page, and the YouTube copies linked from 2.25 are uploads by a private channel ("Barry Belmont") with no usage terms, so no license could be confirmed; MIT OCW 2.06 Fluid Dynamics (Spring 2013) — loaded, but its readings are chapters of White's paid textbook and it has no notes of its own; the Ghia, Ghia & Shin paper as a free source — ScienceDirect returned 403 and Crossref shows only Elsevier's text-mining license, so it is recorded as **optional-paid** (`@ghia-ghia-shin-1982-cavity`) and the P08.18 lesson embeds the published Re = 100 centreline values as data, cited as "published data from Ghia et al. 1982, Table I/II" (no free, licensed reproduction of the tables was found that could be verified); Zang, Chang & Ng (1992, NTRS 19930018240), a status survey of transition tools that adds nothing over Malik (1987) for this level.

Verified 2026-10-17 by loading the live pages; section numbers and titles are copied from the publishers' own tables of contents, section pages or PDF outlines. Every source here is **required-free** unless marked otherwise; sources are identified by an `@`-prefixed identifier and recorded in [sources.yaml](sources.yaml) with license, version, and caveats. **GAP** marks a lesson that no listed free source covers adequately: the lesson must carry original instruction there, and any "candidate" named is unverified until it is added to the manifest.

Conventions: Leishman chapters are at `https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/<slug>/` and are cited **by title only** (numbering changes between builds; several slugs differ from titles). MIT OCW 2.25 pages follow `https://ocw.mit.edu/courses/2-25-advanced-fluid-mechanics-fall-2013/pages/<topic>/` and class notes `.../resources/mit2_25f13_<code>/`; other OCW courses link from the lecture-notes page given in each entry. NASA Glenn pages are either new-site `https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/<slug>/` or legacy `https://www.grc.nasa.gov/www/k-12/airplane/<name>.html`; cite whichever is marked loaded (several legacy pages are now "moved" stubs and several topics are not yet on the new site). LibreTexts URLs begin `https://eng.libretexts.org/Bookshelves/`; a URL written as `.../<NN>%3A_...` continues the book prefix given at its first use. Licenses: MIT OCW, Liburdy and OpenStax are **CC BY-NC-SA 4.0**; Smyth is CC BY-NC 4.0; Ahmari & Kabir, Marchman and CFD Python content are CC BY 4.0 (CFD Python code BSD-3-Clause); Bar-Meir is GNU FDL 1.3 on LibreTexts (author's notice: GFDL 1.2 or later, or Potto license); Leishman is **CC BY-NC-ND 4.0** (link only, never adapt); NASA/NACA material is US Government work ("Public Use Permitted" on NTRS); JHTDB data are ODC-By; the Lee & Moser and KTH datasets carry citation terms but no license, so they are linked for learners to download and are not redistributed — lesson-bundled turbulence data are synthetic. The curriculum links to these sources and writes its own explanations; nothing is copied or adapted.


## 08-01 Fluid Properties & Hydrostatics

Sources: `@leishman-intro-aerospace-flight-vehicles`, `@ntrs-us-standard-atmosphere-1976`, `@libretexts-bar-meir-fluid-mechanics`, `@libretexts-marchman-aerodynamics-aircraft-performance`, `@nasa-glenn-beginners-guide-aeronautics`, `@openstax-university-physics-vol1`, `@mit-ocw-1-060-engineering-mechanics-ii`, `@mit-ocw-16-01-unified-engineering-fluids`, `@ntrs-naca-tr-1235-standard-atmosphere` (supplementary)

- Leishman, "Fundamental Properties of Fluids": What is a Fluid?; Concept of a Continuum; Fluid Pressure; Fluid Density; Equation of State; Bulk Modulus; Viscosity of Fluids; Speed of Sound in Fluid Media; Summary of Properties of Air at MSL Standard
  - .../chapter/fundamental-properties-of-fluids/ (listed on the loaded contents page)
- Leishman, "Fluid Statics & the Hydrostatic Equation": Pascal's Law; Aerostatic/Hydrostatic Pressure Field; Effects of Gravity Alone; Solutions to the Aero-Hydrostatic Equation; Archimedes' Principle; Pressure Forces on Submerged Surfaces
  - .../chapter/fluid-statics/ (loaded)
- Leishman, "Atmospheric Properties": General Properties of the Atmosphere; Reasons for the ISA model; Geopotential Altitude vs Geometric Altitude; Equations of the ISA in the Troposphere; Pressure Altitude and Density Altitude; Tables of ISA Properties; Interpolation of ISA Table Values
  - .../chapter/international-standard-atmosphere-isa/ (loaded)
- U.S. Standard Atmosphere, 1976, Part 1 "Defining constants and equations": 1.2 Basic assumptions and formulas (1.2.1 Adopted constants; 1.2.3 Gravity and geopotential altitude; 1.2.5 Molecular-scale temperature vs. geopotential altitude 0.0 to 84.8520 km); 1.3 Computational equations (1.3.1 Pressure; 1.3.4 Mass density; 1.3.10 Speed of sound; 1.3.11 Dynamic viscosity); 1.4 Selected tabular values (1.4.1 Sea-level values). Part 4 "Main tables": Table I Temperature, pressure, and density for geopotential and geometric altitudes in metric units
  - https://ntrs.nasa.gov/citations/19770009539 (NTRS record loaded via API; PDF loaded, contents pages read)
- Bar-Meir Ch. 1 "Introduction to Fluids": 1.4 Shear Stress; 1.5 Viscosity (1.5.1 General Discussion; 1.5.3 Kinematic Viscosity; 1.5.4 Estimation of the Viscosity); 1.6 Fluid Properties (1.6.1 Fluid Density; 1.6.2 Bulk Modulus). Ch. 4 "Fluids Statics": 4.2 The Hydrostatic Equation; 4.3 Pressure and Density in a Gravitational Field (4.3.1 Constant Density in Gravitational Field; 4.3.2.2 Pressure Measurement; 4.3.3.1 Gas Phase under Hydrostatic Pressure; 4.3.4 The Pressure Effects Due To Temperature Variations); 4.5 Fluid Forces on Surfaces; 4.6 Buoyancy and Stability
  - https://eng.libretexts.org/Bookshelves/Civil_Engineering/Fluid_Mechanics_(Bar-Meir)/00%3A_Introduction (loaded); .../04%3A_Fluids_Statics (loaded); .../4.3.3.1%3A_Gas_Phase_under_Hydrostatic_Pressure (loaded)
- Marchman Ch. 1 "Introduction to Aerodynamics": 1.2 Air, Our Flight Environment; 1.4 The Standard Atmosphere; 1.5 The Troposphere; 1.6 The Stratosphere; 1.7 Altitude Measurement
  - https://eng.libretexts.org/Bookshelves/Aerospace_Engineering/Aerodynamics_and_Aircraft_Performance_3e_(Marchman)/01%3A_Introduction_to_Aerodynamics (loaded)
- NASA Glenn "Earth Atmosphere Model - Metric Units"
  - https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html (loaded; not yet on the new site)
- OpenStax Vol. 1 Ch. 14: 14.1 Fluids, Density, and Pressure; 14.2 Measuring Pressure; 14.4 Archimedes' Principle and Buoyancy
  - https://openstax.org/books/university-physics-volume-1/pages/14-2-measuring-pressure (loaded); .../14-1-fluids-density-and-pressure (pattern)
- MIT 1.060 lecture notes: L1 "Subject Outline, Fluid Concepts, The Continuum Hypothesis"; L3 "Hydrostatics, Forces on Plane and Curved Surfaces"; L4 "Manometry, Stability of Dams and Gates" (handwritten scans). MIT 16.01 fluids F2 "Hydrostatics"
- Supplementary: NACA TR 1235 *Standard Atmosphere — Tables and Data for Altitudes to 65,800 Feet* (1955), the ICAO-standard predecessor, for the historical table format — https://ntrs.nasa.gov/citations/19930090991 (NTRS record and text loaded)
- Extract: dp/dz = −ρg with p = ρRT; isothermal and constant-lapse-rate layers integrated in closed form; geopotential vs. geometric altitude; ISA anchor values (288.15 K, 101 325 Pa, 1.225 kg/m³, −6.5 K/km to 11 km, 216.65 K) taken from USSA-1976 1.4.1 and Table I and used as the build's unit tests; Sutherland viscosity (USSA 1.3.11).

## 08-02 Fluid Kinematics

Sources: `@mit-ocw-16-01-unified-engineering-fluids`, `@mit-ocw-2-25-advanced-fluid-mechanics`, `@libretexts-smyth-all-things-flow` (supplementary), `@nasa-glenn-beginners-guide-aeronautics`, `@leishman-intro-aerospace-flight-vehicles`

- MIT 16.01 fluids notes: F8 "Streamlines, Pathlines, Streaklines" (Fluids – Lecture 8 Notes: 1. Streamlines 2. Pathlines 3. Streaklines); F10 "Substantial Derivative" (1. Substantial Derivative 2. Recast Governing Equations)
  - https://ocw.mit.edu/courses/16-01-unified-engineering-i-ii-iii-iv-fall-2005-spring-2006/resources/f08_new_fall/ (loaded; PDF loaded); .../resources/f10_fall/ (loaded; PDF loaded)
- MIT 2.25 "The Continuum Viewpoint and the Equations of Motion" class notes: "Overview of Lagrangian and Eulerian Descriptions"; "Particle Kinematics"; "Geometric Interpretation of Fluid Kinematics In Steady Shear Flow"
  - https://ocw.mit.edu/courses/2-25-advanced-fluid-mechanics-fall-2013/pages/the-continuum-viewpoint-and-the-equations-of-motion/ (loaded)
- Smyth Ch. 5 "Fluid Kinematics": 5.1 Lagrangian and Eulerian descriptions (5.1.1 The material derivative); 5.2 The Streamfunction; 5.3 Rotation and strain- the relative motion of two nearby particles
  - https://eng.libretexts.org/Bookshelves/Civil_Engineering/All_Things_Flow_-_Fluid_Mechanics_for_the_Natural_Sciences_(Smyth)/05%3A_Fluid_Kinematics/5.01%3A_Lagrangian_and_Eulerian_descriptions (loaded)
- NASA Glenn "Definition of Streamlines" — https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/definition-of-streamlines/ (loaded)
- Leishman, "Fundamental Properties of Fluids": Flow Velocity; Streamtubes. "Types of Fluid Flows": Steady Versus Unsteady Flows; Single Versus Multi-Dimensional Flow
- Extract: D/Dt = ∂/∂t + u·∇ as local plus convective change; integrating dx/dt = u(x,t) for pathlines vs. tracing streaklines from a fixed release point; the three coincide only in steady flow — the visualizer's test case is an unsteady (oscillating) field where they differ.

## 08-03 Control Volume: Mass & Momentum

Sources: `@mit-ocw-2-25-advanced-fluid-mechanics`, `@leishman-intro-aerospace-flight-vehicles`, `@libretexts-bar-meir-fluid-mechanics`, `@mit-ocw-1-060-engineering-mechanics-ii`, `@mit-ocw-16-01-unified-engineering-fluids`, `@nasa-glenn-beginners-guide-aeronautics`

- Sonin, "Fundamental Laws of Motion for Particles, Material Volumes, and Control Volumes" (MIT 2.25, March 2003): 1 Basic laws for material volumes; 2 The transformation to control volumes (2.4 Reynolds' material-volume to control-volume transformation); 3 Basic laws for control volumes (3.1 Mass conservation; 3.2 Linear momentum theorem; 3.3 Angular momentum theorem). Also "On Choosing and Using Control Volumes: Six Ways of Applying the Integral Mass Conservation Theorem to a Simple Problem"
  - https://ocw.mit.edu/courses/2-25-advanced-fluid-mechanics-fall-2013/resources/mit2_25f13_fundam_law-son/ (loaded; PDF loaded); topic page .../pages/control-volume-theorems-and-applications/ (loaded)
- Leishman, "Equations of Fluid Motion": Finite Control Volume Approach; Mass Flow and Mass Flux; Momentum & Energy Flow Rates; Extensive & Intensive Properties; Reynolds Transport Theorem. "Continuity Equation": Continuity Equation – Integral Form; Flow Through a Converging/Diverging Duct. "Momentum Equations": Momentum Equation – Integral Approach; Simplifications to the Momentum Equation. "Applications of the Conservation Laws": Force on a Pipe Bend; Jet Engine Performance
- Bar-Meir Ch. 5 "The Control Volume and Mass Conservation": 5.2 Control Volume; 5.3 Continuity Equation; 5.4 Reynolds Transport Theorem. Ch. 6 "Momentum Conservation for Control Volume": 6.1 Momentum Governing Equation; 6.2 Momentum Equation Application
  - .../Fluid_Mechanics_(Bar-Meir)/05%3A_The_Control_Volume_and_Mass_Conservation/5.4%3A_Reynolds_Transport_Theorem (loaded); .../06%3A_Momentum_Conservation_for_Control_Volume (loaded)
- MIT 1.060: L5 "Mass and Volume Conservation"; L9 "Conservation Laws for Finite Control Volumes, Reynolds Transport Theorem"; L10 "Conservation of Momentum for Finite Control Volumes"; L11 "Applications of Momentum"
- MIT 16.01 fluids: F6 "Control Volume, Mass Conservation"; F7 "Momentum Conservation"; F9 "Momentum Integral Applications"
- NASA Glenn "Conservation of Mass" — https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/conservation-of-mass/ (loaded)
- Extract: RTT as the bridge from system to control volume; choosing the CV so unknown forces appear once; jet thrust F = ṁ(Ve − V0) + (pe − p0)Ae from a CV around an engine (Leishman "Jet Engine Performance"), feeding the Phase 12 thrust equation.

## 08-04 Bernoulli & Energy Equation

Sources: `@leishman-intro-aerospace-flight-vehicles`, `@nasa-glenn-beginners-guide-aeronautics`, `@libretexts-marchman-aerodynamics-aircraft-performance`, `@openstax-university-physics-vol1`, `@mit-ocw-2-25-advanced-fluid-mechanics`, `@libretexts-bar-meir-fluid-mechanics`, `@libretexts-liburdy-intermediate-fluid-mechanics`, `@mit-ocw-1-060-engineering-mechanics-ii`, `@mit-ocw-16-01-unified-engineering-fluids`, `@libretexts-ahmari-kabir-fluid-mechanics-lab-manual`

- Leishman, "Energy Equation & Bernoulli's Equation": Bernoulli's Equation; Bernoulli's Equation from Euler's Equation; Other Forms of the Bernoulli Equation. "Applications of the Conservation Laws": Flow in a Venturi; Total, Static, & Dynamic Pressures; Pitot Tube & Pitot-Static Tube; Flow Speed in a Wind Tunnel (Take 2). "Airspeed Definitions & Measurement": Definitions of Airspeeds; Pitot-Static System
  - .../chapter/energy-equation/ and .../chapter/applications-of-the-conservation-laws/ (listed on the loaded contents page); .../chapter/determination-of-airspeed/ (loaded)
- NASA Glenn: "Bernoulli's Equation"; "Pitot – Static Tube – Speedometer"; "Dynamic Pressure"; "Venturi Theory"
  - https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/bernoullis-equation/ (loaded); .../pitot-static-tube-speedometer/ (loaded); .../dynamic-pressure/ (loaded); .../venturi-theory/ (pattern)
- Marchman Ch. 1: 1.8 Bernoulli's Equation; 1.9 Airspeed Measurement; 1.10 Bernoulli's Equation for Compressible Flow
- OpenStax Vol. 1: 14.5 Fluid Dynamics; 14.6 Bernoulli's Equation (Energy Conservation and Bernoulli's Equation; Applications of Bernoulli's Principle — Velocity measurement)
  - https://openstax.org/books/university-physics-volume-1/pages/14-6-bernoullis-equation (loaded)
- MIT 2.25 "Inviscid Flow and Bernoulli" class notes: "Equation of Motion in Streamline Coordinates"; "Unsteady Bernoulli Equation"
  - https://ocw.mit.edu/courses/2-25-advanced-fluid-mechanics-fall-2013/pages/inviscid-flow-and-bernoulli/ (loaded)
- Bar-Meir Ch. 7 "Energy Conservation": 7.1 The First Law of Thermodynamics; 7.3 Approximation of Energy Equation (7.3.2 Energy Equation in Frictionless Flow and Steady State); 7.4.5 Energy Losses in Incompressible Flow. Liburdy Ch. 3 "Bernoulli Equation"
  - .../Fluid_Mechanics_(Bar-Meir)/07%3A_Energy_Conservation (loaded); https://eng.libretexts.org/Bookshelves/Civil_Engineering/Intermediate_Fluid_Mechanics_(Liburdy)/03%3A_Bernoulli_Equation (loaded)
- MIT 1.060: L7 "The Bernoulli Equation for Ideal Fluids"; L12 "The Bernoulli Equation from Energy Considerations"; L13 "Special Forms and Interpretations of the Bernoulli Equation". MIT 16.01 fluids F13 "Bernoulli's Equation". Bench version: Ahmari & Kabir 1.2 "Experiment #2: Bernoulli's Theorem Demonstration"
- Extract: p + ½ρV² + ρgz = const along a streamline (steady, inviscid, incompressible); V = √(2(p0 − p)/ρ) for a pitot-static probe; IAS/EAS/TAS distinctions (feeds P15.01); energy equation with head loss as the bridge to 08-07.

## 08-05 Dimensional Analysis & Similarity

Sources: `@leishman-intro-aerospace-flight-vehicles`, `@mit-ocw-2-25-advanced-fluid-mechanics`, `@libretexts-bar-meir-fluid-mechanics`, `@mit-ocw-16-100-aerodynamics`, `@libretexts-marchman-aerodynamics-aircraft-performance`, `@mit-ocw-16-01-unified-engineering-fluids`, `@nasa-glenn-beginners-guide-aeronautics`

- Leishman, "Dimensional Analysis": Basis of Dimensional Analysis; The Buckingham Π Method; Π Method Applied to the Drag on a Sphere; Π Method Applied to a Wing; When Things Go Wrong in the Buckingham Π Method. "Dynamic Similarity": Similarity Criteria; Aerodynamic Similarity; Helpful Scaling Laws; Partial Similarity; Summary of Similarity Parameters. "Mach Number & Reynolds Number": Definition of Reynolds Number; Physical Interpretation of Reynolds Number. "Wind Tunnel Testing": Model Scaling & Similarity Parameters
- MIT 2.25 "Dimensional Analysis" class notes: "The Physical Basis of Dimensional Analysis"; "Similarity"; "The Buckingham Pi Theorem in Dimensional Analysis" (Historical Note; The Theorem); Rayleigh, "The Principle of Similitude" (1915) linked as a reading
  - https://ocw.mit.edu/courses/2-25-advanced-fluid-mechanics-fall-2013/resources/mit2_25f13_the_buckingham/ (loaded; PDF loaded); .../pages/dimensional-analysis/ (loaded)
- Bar-Meir Ch. 9 "Dimensional Analysis": 9.1 Introductory Remarks; 9.2 Buckingham–π–Theorem (9.2.1 Construction of the Dimensionless Parameters; 9.2.3 Implementation of Construction of Dimensionless Parameters; 9.2.4 Similarity and Similitude); 9.4 Summary of Dimensionless Numbers
  - .../Fluid_Mechanics_(Bar-Meir)/09%3A_Dimensional_Analysis/9.2%3A_Buckingham%E2%80%93%CF%80%E2%80%93Theorem (loaded)
- MIT 16.100 Lecture 31 "Similarity in Wind Tunnel Testing" (why M∞ and Re cannot both be matched)
  - https://ocw.mit.edu/courses/16-100-aerodynamics-fall-2005/resources/16100lectre31_cj/ (loaded; PDF loaded)
- Marchman 1.12 Force Coefficients; 1.13 "Similarity Parameters". MIT 16.01 fluids F4 "Dimensional Analysis and Similarity"; F5 "Reynolds and Mach Number". NASA Glenn "Similarity parameters" — https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/similarity-parameters-2/ (loaded)
- Extract: Π groups as the null space of the dimension matrix (ties to P01.14 rank/null space — the build computes them that way); Re and M as the two aero similarity parameters; the model-scaling conflict (Re matching needs higher speed or pressure, cryogenic tunnels) and partial similarity.

## 08-06 Navier-Stokes Equations

Sources: `@mit-ocw-2-25-advanced-fluid-mechanics`, `@leishman-intro-aerospace-flight-vehicles`, `@libretexts-liburdy-intermediate-fluid-mechanics`, `@libretexts-bar-meir-fluid-mechanics`, `@libretexts-smyth-all-things-flow` (supplementary)

- MIT 2.25 "Equations of Viscous Flow" class notes: "Equation of Motion for Viscous Flow"; "Couette & Poiseuille Flows" (October 2008: plane wall-driven Couette flow, pressure-driven Poiseuille flow, flow rates and wall stresses); "Criteria for Locally Fully Developed Viscous Flow". "More Complex Viscous-Dominated Flows": "Flow Inside a Cylinder Which is Suddenly Rotated"
  - https://ocw.mit.edu/courses/2-25-advanced-fluid-mechanics-fall-2013/resources/mit2_25f13_couet_and_pois/ (loaded; PDF loaded); .../pages/equations-of-viscous-flow/ (loaded); .../pages/more-complex-viscous-dominated-flows/ (loaded)
- Leishman, "Momentum Equations": Momentum Equation – Differential Form; Navier-Stokes Equations. "Viscous-Dominated Flows": Stokes Flow Between Parallel Plates; Poiseuille Flow in a Rectangular Channel; Stokes' First and Second Problems; Taylor-Couette Flow
- Liburdy Ch. 8 "Introduction to Viscous Flows": Basics of Viscous Forces; Viscous Stress; Navier-Stokes Equation; Examples of Exact Solutions to Navier-Stokes Equation (includes Stokes' first problem)
  - https://eng.libretexts.org/Bookshelves/Civil_Engineering/Intermediate_Fluid_Mechanics_(Liburdy)/08%3A_Introduction_to_Viscous_Flows (loaded)
- Bar-Meir Ch. 8 "Differential Analysis": 8.2 Mass Conservation; 8.4 Momentum Conservation; 8.5 Derivations of the Momentum Equation; 8.6 Boundary Conditions and Driving Forces; 8.7 Examples for Differential Equation (Navier-Stokes)
  - .../Fluid_Mechanics_(Bar-Meir)/08%3A_Differential_Analysis/8.5%3A_Derivations_of_the_Momentum_Equation (loaded)
- Smyth Ch. 6 "Fluid Dynamics": 6.2 Mass conservation; 6.3 Momentum conservation; 6.8 Summary: the equations of motion; 6.9 Boundary conditions
- Extract: Cauchy momentum equation plus Newtonian constitutive law → N-S; no-slip; exact solutions as the solver's analytic test cases — Couette (linear), Poiseuille (parabolic), Stokes' first problem (erf similarity solution) — checked against the P03 PDE solver.

## 08-07 Pipe Flow & Losses

Sources: `@leishman-intro-aerospace-flight-vehicles`, `@mit-ocw-1-060-engineering-mechanics-ii`, `@libretexts-ahmari-kabir-fluid-mechanics-lab-manual`, `@openstax-university-physics-vol1`

- Leishman, "Internal Flows": Entrance Length; Hydraulic Diameter; Laminar Flow in Pipes; Turbulent Flows in Pipes; Colebrook-White Equation; Moody Chart; Pressure Drops Over Bends & Corners; Pumps
  - .../chapter/internal-flows/ (loaded)
- MIT 1.060: L16 "Pipe Friction Losses, Moody Diagram, Minor Losses" (handwritten: frictional head loss in circular pipes, Darcy–Weisbach f); L17 "Pipe Flow Analysis"; L18 "Pump Specification and Requirements"; L15 "Velocity Distribution in Laminar Flows, Turbulence, The Reynolds Number"
  - https://ocw.mit.edu/courses/1-060-engineering-mechanics-ii-spring-2006/resources/lecture16/ (loaded; PDF loaded, image-only)
- Ahmari & Kabir: 1.3 "Experiment #3: Energy Loss in Pipe Fittings"; 1.4 "Experiment #4: Energy Loss in Pipes"
  - https://eng.libretexts.org/Bookshelves/Mechanical_Engineering/Applied_Fluid_Mechanics_Lab_Manual_(Ahmari_and_Kabir)/01%3A_Lab_Manual/01.4%3A_Experiment_%234%3A_Energy_Loss_in_Pipes (loaded during contents crawl)
- OpenStax Vol. 1 14.7 Viscosity and Turbulence: Laminar Flow Confined to Tubes: Poiseuille's Law; Flow and Resistance as Causes of Pressure Drops
  - https://openstax.org/books/university-physics-volume-1/pages/14-7-viscosity-and-turbulence (loaded)
- Extract: h_f = f(L/D)V²/2g; f = 64/Re laminar; Colebrook solved by Newton iteration (P03.03) and checked against an explicit fit; minor losses K·V²/2g summed along the feed line; Δp budget from tank to injector for P12.06.
- Partial GAP: no verified free source sizes a rocket propellant feed line (tank pressure, line losses, injector Δp, cavitation margin); that worked case is original.

## 08-08 Potential Flow

Sources: `@leishman-intro-aerospace-flight-vehicles`, `@libretexts-liburdy-intermediate-fluid-mechanics`, `@mit-ocw-2-25-advanced-fluid-mechanics`, `@libretexts-bar-meir-fluid-mechanics`, `@mit-ocw-16-01-unified-engineering-fluids`

- Leishman, "Potential Flow Theory": Definition of a Potential Flow; Potential Flow Equations; Velocity Potential & Stream Function; Elementary Flows; Understanding Circulation and Vorticity; Superposition of Elementary Flows (includes the Rankine oval); Method of Images; Complex Variable Methods in Potential Flow
  - .../chapter/potential-flows/ (loaded)
- Liburdy Ch. 4 "Potential Flow Basics" (Continuity Equation; Streamfunction) and Ch. 6 "Potential Flows": Basic Flows (Uniform Flow; Source/Sink Flow; Vortex Flow); Superposition; DOUBLET; Uniform Flow over a Cylinder; Rankine Half-Body; Enclosed Bodies; Blasius Theorem and Lift Force for an Arbitrary Body
  - https://eng.libretexts.org/Bookshelves/Civil_Engineering/Intermediate_Fluid_Mechanics_(Liburdy)/06%3A_Potential_Flows (loaded)
- MIT 2.25 "Potential Flow Theory": class notes "Equation of Motion in Streamline Coordinates", "Inviscid Flow Equation Sheet"; problem "The Magnus Effect" with solution
  - https://ocw.mit.edu/courses/2-25-advanced-fluid-mechanics-fall-2013/pages/potential-flow-theory/ (loaded)
- Bar-Meir Ch. 10 "Inviscid Flow or Potential Flow": 10.2.1 Streamline and Stream function; 10.2.3 The Connection Between the Stream Function and the Potential Function; 10.3 Potential Flow Functions Inventory (10.3.1 Flow Around a Circular Cylinder; 10.3.1.1 Adding Circulation to a Cylinder); 10.4.1 Complex Potential and Complex Velocity
  - .../Fluid_Mechanics_(Bar-Meir)/10%3A_Inviscid_Flow_or_Potential_Flow/10.3_Potential_Flow_Functions_Inventory (loaded)
- MIT 16.01 fluids: F12 "Stream Function and Velocity Potential"; F14 "Incompressible and Irrotational Flows"; F15 "Sources, Sinks and Doublets"; F16 "Vortex, Lifting Flow"; F18 "Lift Prediction, Vortex Sheets"
- Extract: ψ and φ for uniform flow, source/sink, doublet, vortex; superposition to the Rankine oval and the cylinder with circulation; stagnation-point location vs. Γ; Cp = 1 − 4 sin²θ; L′ = ρ∞V∞Γ (Kutta–Joukowski; Leishman/Liburdy "Blasius Theorem") — prerequisite for P10.03.

## 08-09 Vorticity & Circulation

Sources: `@leishman-intro-aerospace-flight-vehicles`, `@mit-ocw-2-25-advanced-fluid-mechanics`, `@libretexts-smyth-all-things-flow` (supplementary), `@mit-ocw-16-01-unified-engineering-fluids`

- Leishman, "Potential Flow Theory": Understanding Circulation and Vorticity. "Lifting Line & Finite Wing Theory": Helmholtz's Vortex Theorems; Theory of Vortex Lines (the chapter states the Biot–Savart law)
  - .../chapter/lifting-line-theory/ (loaded)
- MIT 2.25 "Vorticity and Circulation" (lectures 19–20; readings in Kundu & Cohen Ch. 5 — 5.2 Kelvin's Circulation Theorem, 5.5 the law of Biot and Savart — are paid). Class note "Viscous Decay of a Free Vortex" (a problem set as worked notes: a spinning rod stopped in a fluid, solved by a similarity variable built from r and νt, with a profile of the form 1 − exp(−η) — the Lamb–Oseen vortex, not named as such)
  - https://ocw.mit.edu/courses/2-25-advanced-fluid-mechanics-fall-2013/pages/vorticity-and-circulation/ (loaded); .../resources/mit2_25f13_viscousdecay/ (loaded; PDF loaded)
- Smyth Ch. 4 "Tensor Calculus": 4.3 Vorticity and Circulation. Ch. 7 "Vortices": 7.2 Vortex dynamics in a homogeneous, inviscid fluid (7.2.1 The vorticity equation; 7.2.2 Vortex filaments; 7.2.4 Vortex tubes; the page also discusses Kelvin's circulation theorem); 7.3 Viscous effects (7.3.1 The vorticity equation for a viscous fluid; 7.3.2 The Burgers vortex)
  - https://eng.libretexts.org/Bookshelves/Civil_Engineering/All_Things_Flow_-_Fluid_Mechanics_for_the_Natural_Sciences_(Smyth)/07%3A_Vortices/7.02%3A_Vortex_dynamics_in_a_homogeneous%2C_inviscid_fluid (loaded); .../7.03%3A_Viscous_effects (loaded)
- MIT 16.01 fluids F11 "Vorticity and Circulation"; 16.03–16.04 F6 "3D Vortex Filaments, Lifting Line Theory"
- Extract: ω = ∇×u, Γ = ∮u·dl = ∬ω·n dA; Kelvin's theorem (DΓ/Dt = 0 for inviscid, barotropic, conservative forces) as the reason starting vortices exist; Biot–Savart for a straight segment (the vortex-filament calculator, reused in P10.11 VLM); Lamb–Oseen core growth r_c ∝ √(νt).
- Partial GAP: no verified free source derives the Lamb–Oseen vortex as a named result with its core-radius definition (2.25 poses it as an exercise; Smyth treats the Burgers vortex instead); the derivation and the wake-vortex-decay example are original.

## 08-10 Laminar Boundary Layers

Sources: `@leishman-intro-aerospace-flight-vehicles`, `@mit-ocw-16-100-aerodynamics`, `@libretexts-liburdy-intermediate-fluid-mechanics`, `@mit-ocw-2-25-advanced-fluid-mechanics`, `@nasa-glenn-beginners-guide-aeronautics`

- Leishman, "Boundary Layer Flows": Laminar & Turbulent Boundary Layers; Shear Stresses & Skin Friction; Development of a Boundary Layer; Other Properties of Boundary Layers; Pressure Gradients & Flow Separation; Boundary Layer Equations; Laminar Boundary Layer – Blasius Solution
  - .../chapter/boundary-layers/ (loaded)
- MIT 16.100: Lecture 29 "Laminar Boundary Layer Order of Magnitude Analysis"; Lecture 30 "Solutions of the Laminar Boundary Layer Equations" (Flat Plate (Laminar): Blasius Solution; Falkner-Skan Solutions); Lecture 25 "Falkner-Skan Flows"
  - https://ocw.mit.edu/courses/16-100-aerodynamics-fall-2005/resources/16100lectre30_cj/ (loaded; PDF loaded)
- Liburdy Ch. 9 "Boundary Layer Flows": Physical Attributes of a Boundary Layer; Boundary Layer Equations; Laminar Boundary Layer Solution; Pressure Gradient Effects
  - https://eng.libretexts.org/Bookshelves/Civil_Engineering/Intermediate_Fluid_Mechanics_(Liburdy)/09%3A_Boundary_Layer_Flows (loaded)
- McKinley, "The Generalized Boundary Layer Equations" (MIT 2.25, November 2004)
  - https://ocw.mit.edu/courses/2-25-advanced-fluid-mechanics-fall-2013/resources/mit2_25f13_generalizedeqn/ (loaded; PDF loaded)
- NASA Glenn "Boundary Layer" — https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/boundary-layer/ (loaded)
- Extract: Prandtl's scaling δ/x ~ Re_x^(−1/2); similarity variable η = y√(U/νx); f‴ + ½ff″ = 0 solved by shooting (P02 BVP methods); f″(0) = 0.332, δ99 ≈ 5.0x/√Re_x, δ* = 1.72x/√Re_x, θ = 0.664x/√Re_x as the solver's checks.

## 08-11 Integral Boundary Layer Methods

Sources: `@mit-ocw-16-100-aerodynamics`, `@libretexts-liburdy-intermediate-fluid-mechanics`, `@leishman-intro-aerospace-flight-vehicles`

- MIT 16.100: Lecture 22 "Integral Boundary Layer Equations" (displacement thickness, momentum thickness, the integral momentum equation); Lecture 23 "Correlation Methods for Integral Boundary Layers" (Thwaites 1949: λ = θ²(due/dx)/ν, the correlation assumption and the resulting quadrature for θ²); Lecture 24 "Method of Assumed Profiles"
  - https://ocw.mit.edu/courses/16-100-aerodynamics-fall-2005/resources/16100lectre22_cg/ (loaded; PDF loaded); .../16100lectre23/ (loaded; PDF loaded); .../16100lectre24_cg/ (loaded; PDF loaded)
- Liburdy Ch. 5 "Integral Boundary Layer Relationships": Scaling Parameters; Boundary layer thickness; General Integral Boundary Layer Equations; Solution Method
  - https://eng.libretexts.org/Bookshelves/Civil_Engineering/Intermediate_Fluid_Mechanics_(Liburdy)/05%3A_Integral_Boundary_Layer_Relationships (loaded)
- Leishman, "Boundary Layer Flows": von Kármán Momentum Integral Equation
- Extract: dθ/dx + (2 + H)(θ/ue)due/dx = Cf/2; Thwaites θ² = (0.45ν/ue⁶)∫ue⁵dx; separation at λ ≈ −0.09; flat-plate check against 08-10 Blasius; feeds the P10.06 viscous–inviscid coupling.
- Partial GAP: Lecture 23 gives Thwaites's integral but not the l(λ) and H(λ) correlation fits needed for Cf and δ*, and no verified free source tabulates them; those fits (cited to Thwaites 1949 / later texts as published data) and the marching implementation are original.

## 08-12 Transition

Sources: `@ntrs-malik-1987-stability-theory-lfc` (supplementary), `@leishman-intro-aerospace-flight-vehicles`, `@libretexts-liburdy-intermediate-fluid-mechanics`, `@mit-ocw-2-27-turbulent-flow-transport` (supplementary)

- Malik, "Stability Theory Applications to Laminar-Flow Control" (NASA Langley NLF/LFC conference, 1987): NEED FOR TRANSITION PREDICTION METHODS; AVAILABLE METHODS; CORRELATIONS; PREDICTION METHODS BASED UPON BOUNDARY-LAYER STABILITY THEORY (Orr-Sommerfeld equation; e^N method; N-factor envelopes)
  - https://ntrs.nasa.gov/citations/19900003197 (NTRS record and text loaded)
- Leishman, "Boundary Layer Flows": Laminar & Turbulent Boundary Layers; Surface Roughness; Laminar Separation Bubbles; Boundary Layer Flow Control
- Liburdy Ch. 10 "Introduction to Turbulence Effects": Initial Laminar Starting Length (the chapter discusses the critical Reynolds number and transition)
  - https://eng.libretexts.org/Bookshelves/Civil_Engineering/Intermediate_Fluid_Mechanics_(Liburdy)/10%3A_Introduction_to_Turbulence_Effects (loaded)
- MIT 2.27 outline "Concepts in Turbulence": 3.1 Comments on laminar flow, its stability, and the transition to turbulent flow
  - https://ocw.mit.edu/courses/2-27-turbulent-flow-and-transport-spring-2002/resources/concepts/ (loaded; PDF loaded)
- Extract: normal-mode ansatz → Orr–Sommerfeld eigenproblem (P01.16 eigenvalues; solved numerically with Chebyshev or finite differences); neutral curve and Re_crit; amplification N = ∫(−αi)dx and the e^N criterion with N ≈ 9 for low-disturbance flight; free-stream turbulence and roughness as bypass routes.
- Partial GAP: no verified free source gives an introductory derivation of the Orr–Sommerfeld equation or a worked e^N calculation for a Blasius layer; Malik shows the method in use but is not a tutorial. The derivation, the eigenvalue solver and the e^N marching on the 08-10 profile are original.

## 08-13 Turbulence Fundamentals

Sources: `@leishman-intro-aerospace-flight-vehicles`, `@mit-ocw-16-100-aerodynamics`, `@mit-ocw-2-29-numerical-fluid-mechanics`, `@libretexts-liburdy-intermediate-fluid-mechanics`, `@mit-ocw-2-27-turbulent-flow-transport` (supplementary), `@libretexts-smyth-all-things-flow` (supplementary), `@jhtdb-johns-hopkins-turbulence-database` (supplementary)

- Leishman, "Turbulent Flows": Visualizing Turbulent Flows; Reynolds Decomposition; Navier-Stokes Equations; Including Turbulent Flow Effects; Turbulence Scales & Turbulent Kinetic Energy (Kolmogorov scales)
  - .../chapter/turbulent-flows/ (loaded)
- MIT 16.100 Lecture 26 "Effect of Turbulent Fluctuations on Mean Flow: Reynolds-Averaging"
  - https://ocw.mit.edu/courses/16-100-aerodynamics-fall-2005/resources/16100lectre26_cg/ (loaded; PDF loaded)
- MIT 2.29 Lecture 24 ("Turbulent Flows and their Numerical Modeling — Properties of Turbulent Flows; Stirring and Mixing; Energy Cascade and Scales") and Lecture 25 (Kolmogorov microscale; Direct Numerical Simulations (DNS); Reynolds-averaged Navier-Stokes (RANS); Large-Eddy Simulations (LES))
  - https://ocw.mit.edu/courses/2-29-numerical-fluid-mechanics-spring-2015/resources/mit2_29s15_lecture24/ (loaded; PDF loaded); .../mit2_29s15_lecture25/ (loaded; PDF loaded)
- Liburdy Ch. 10: Reynolds Averaging
- MIT 2.27 outline "Concepts in Turbulence": 3.2 Features of turbulent flows; 3.3 The energy cascade and the smallest scales of turbulent motion (the Kolmogorov microscale); 3.5 Statistical averages of random quantities in turbulent flow; 3.6 The mean-flow equations … The closure problem
- Smyth 5.4 Pancakes and noodles- the geometry of turbulence
- Data (supplementary, ODC-By): JHTDB forced isotropic turbulence via the cutout service or Python client, for a measured E(k) to compare with the −5/3 law — https://turbulence.idies.jhu.edu/home (loaded); https://turbulence.idies.jhu.edu/citing (loaded)
- Extract: u = ū + u′; mean, variance, autocorrelation and integral time scale from a signal; Taylor's frozen-turbulence hypothesis to get E(k) from a time series by FFT (P03 spectral tools); η = (ν³/ε)^(1/4); the Reynolds-stress closure problem that motivates 08-15.
- The bundled signals are synthetic (generated with a prescribed spectrum); JHTDB data are optional for learners to fetch themselves.

## 08-14 Turbulent Boundary Layers

Sources: `@leishman-intro-aerospace-flight-vehicles`, `@libretexts-liburdy-intermediate-fluid-mechanics`, `@mit-ocw-2-27-turbulent-flow-transport` (supplementary), `@lee-moser-2015-channel-dns-data` (supplementary), `@kth-schlatter-orlu-tbl-dns-data` (supplementary), `@jhtdb-johns-hopkins-turbulence-database` (supplementary)

- Leishman, "Turbulent Flows": Wall-Bounded Turbulent Flows (law of the wall, log law). "Boundary Layer Flows": Turbulent Boundary Layers; Skin Friction Drag on an Airfoil (the chapter uses the 1/7-power turbulent profile)
  - .../chapter/turbulent-flows/ (loaded); .../chapter/boundary-layers/ (loaded)
- Liburdy Ch. 10 "Introduction to Turbulence Effects": Turbulent Boundary Layer; Scaling for Turbulent Boundary Layer Flows; Initial Laminar Starting Length; Universal Velocity Profile
- MIT 2.27 outlines: "Bounded Flows I: General Near-Wall Scaling Laws; Flows in Pipes, Channels, etc." (5.1 inner/outer layers, Prandtl's universal law of the wall, viscous and logarithmic sublayers, velocity-defect law); "Introduction to Turbulent Boundary Layers" (6.4 Coles's mean-velocity profile with wake function; 6.6 smooth and rough flat plates)
  - https://ocw.mit.edu/courses/2-27-turbulent-flow-and-transport-spring-2002/resources/simple_turb_wall_flows/ (loaded; PDF loaded); .../resources/boundary_layers/ (loaded; PDF loaded)
- Data (supplementary; citation terms, no license — learners download, nothing redistributed): Lee & Moser channel DNS mean profiles (LM_Channel_<Re>_mean_prof.dat, Re_tau 180–5200) — https://turbulence.oden.utexas.edu/channel2015/content/Data_2015_5200.html (loaded); Schlatter & Örlü zero-pressure-gradient TBL profiles vel_<Re_theta>.prof, Re_theta 670–4060 — https://www.mech.kth.se/~pschlatt/DATA/ (loaded)
- Extract: u⁺ = y⁺ in the viscous sublayer, u⁺ = (1/κ)ln y⁺ + B in the log layer (κ ≈ 0.41, B ≈ 5.0); fitting κ and B to a DNS profile; Cf correlations for a flat plate (e.g. Cf = 0.0592 Re_x^(−1/5) and Prandtl–Schlichting forms) compared with laminar Blasius; feeds the P10.13 drag build-up.
- Lesson-bundled profiles are synthetic (built from the composite law of the wall); comparison with Lee & Moser or KTH data is a learner-download exercise.

## 08-15 Turbulence Modeling

Sources: `@leishman-intro-aerospace-flight-vehicles`, `@mit-ocw-2-29-numerical-fluid-mechanics`, `@mit-ocw-2-27-turbulent-flow-transport` (supplementary)

- Leishman, "Turbulent Flows": Eddy Viscosity Models; Turbulence Models (the chapter discusses mixing length, k-epsilon and Spalart–Allmaras); Numerical Solution of the RANS Equations
- MIT 2.29 Lecture 25: Reynolds-averaged Navier-Stokes (RANS); Turbulence closures: Eddy viscosity and diffusivity, Mixing-length Models, k-ε Models
  - https://ocw.mit.edu/courses/2-29-numerical-fluid-mechanics-spring-2015/resources/mit2_29s15_lecture25/ (loaded; PDF loaded)
- MIT 2.27 outlines: "Concepts in Turbulence" 3.8 Prandtl's (momentum-transfer) mixing length hypothesis; "Introduction to Turbulence Models": 8.1 Approaches to closure; 8.2 Equations for … the turbulence intensity k; 8.3 Modeling of the terms in the k-equation; 8.4 A calibration of the empirical constants … Wall boundary condition on k; 8.5 An equation for ε; 8.6 Application of the k-epsilon model in some simple cases
  - https://ocw.mit.edu/courses/2-27-turbulent-flow-and-transport-spring-2002/resources/8_turbulence_models/ (loaded; PDF loaded)
- Extract: νt = ℓm²|dū/dy| with ℓm = κy (and van Driest damping) reproduces the 08-14 log law in a 1D channel/Couette RANS solve; k-ε transport equations and standard constants at overview level; why wall functions exist.
- Partial GAP: no verified free source works a mixing-length channel-flow solution step by step or lists the standard k-ε constants with their calibration; the 1D RANS solver and the constants table (cited to Launder & Spalding 1974 as published values) are original.

## 08-16 Flow Separation & Drag

Sources: `@leishman-intro-aerospace-flight-vehicles`, `@nasa-glenn-beginners-guide-aeronautics`, `@mit-ocw-1-060-engineering-mechanics-ii`, `@mit-ocw-2-25-advanced-fluid-mechanics`

- Leishman, "Bluff Body Flows": Drag Coefficients; Drag Coefficients of Simple Shapes; Flow Patterns Around a Circular Cylinder; Drag of Circular Cylinders & Spheres; Roughness Effects on Cylinders & Spheres; Why Do Golf Balls Have Dimples?; Vortex-Induced Vibration; Streamlining. "Time-Dependent Flows": Strouhal Number. "Boundary Layer Flows": Pressure Gradients & Flow Separation; Flow Separation on an Airfoil
  - .../chapter/bluff-body-flows/ (loaded); .../chapter/unsteady-flows/ (listed on the loaded contents page)
- NASA Glenn: "Drag of a Sphere" (Cd vs. Re figure with numbered cases of flow past a cylinder and sphere, from ideal flow to a turbulent wake); "Shape Effects on Drag"
  - https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/drag-of-a-sphere/ (loaded); .../shape-effects-on-drag/ (loaded)
- MIT 1.060 L20 "Drag Forces on Immersed Bodies"; L21 "Lift Forces on Immersed Bodies"
- MIT 2.25 "Boundary Layers, Separation, and Drag" (lectures 21–23; readings in Panton Ch. 20 and Kundu & Cohen Ch. 9 are paid)
  - https://ocw.mit.edu/courses/2-25-advanced-fluid-mechanics-fall-2013/pages/boundary-layers-separation-and-drag/ (loaded)
- Extract: adverse pressure gradient → separation (link to Thwaites λ in 08-11); pressure vs. friction drag; Cd(Re) regimes for cylinder and sphere and the drag crisis; St = fD/U ≈ 0.2 for vortex shedding; d'Alembert's paradox as the potential-flow (08-08) contrast.
- Partial GAP: no verified free source provides tabulated cylinder/sphere Cd(Re) data with a license; the explorer uses published correlations (e.g. Stokes, Oseen and empirical fits) written up originally, and Leishman's/NASA's figures are linked, not reproduced.

## 08-17 Intro to Compressible Flow

Sources: `@nasa-glenn-beginners-guide-aeronautics`, `@mit-ocw-16-100-aerodynamics`, `@mit-ocw-16-01-unified-engineering-fluids`, `@libretexts-bar-meir-fluid-mechanics`, `@ntrs-naca-tr-1135-compressible-flow`, `@leishman-intro-aerospace-flight-vehicles`, `@ntrs-us-standard-atmosphere-1976`

- NASA Glenn: "Speed of Sound"; "Mach Number"; "Isentropic Flow Equations"
  - https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/speed-of-sound-interactive/ (loaded); https://www.grc.nasa.gov/www/k-12/airplane/mach.html (loaded); https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html (loaded)
- MIT 16.100: Lecture 39 "Derivation of Sound Wave Properties"; Lecture 38 "Behavior of Isentropic Flow in Quasi-1D"; Lecture 37 "Critical Mach Number"
  - https://ocw.mit.edu/courses/16-100-aerodynamics-fall-2005/resources/16100lectre39_cj/ (loaded; PDF loaded); .../16100lectre38_cj/ (loaded; PDF loaded)
- MIT 16.03–16.04 fluids: F11 "Introduction to Compressible Flows, Thermodynamic Concepts"; F12 "Energy Conservation"; F13 "Stagnation Quantities, Introduction to Shock Waves"; F14 "Normal Shock Waves, Speed of Sound"; F15 "Mach Number Relations, Normal Shock Properties"; F16 "Shock Losses, Compressible Flow Pitot-Tube"
- Bar-Meir Ch. 11 "Compressible Flow One Dimensional": 11.1 What is Compressible Flow?; 11.3 Speed of Sound (11.3.2 Speed of Sound in Ideal and Perfect Gases); 11.4 Isentropic Flow (11.4.1 Stagnation State for Ideal Gas Model; 11.4.2 Isentropic Converging–Diverging Flow in Cross Section; 11.4.3.1 The pressure Mach number relationship; 11.4.3.2 Relationship Between the Mach Number and Cross Section Area)
  - .../Fluid_Mechanics_(Bar-Meir)/11%3A_Compressible_Flow_One_Dimensional/11.4_Isentropic_Flow (loaded)
- NACA TR 1135 (1953): equations for continuous one-dimensional flow and the tables of dimensionless ratios vs. Mach number for air as a perfect gas (γ = 7/5), used as the check table
  - https://ntrs.nasa.gov/citations/19930091059 (NTRS record and text loaded)
- Leishman, "Mach Number & Reynolds Number": Definition of Mach Number; Significance of the Speed of Sound; Mach Number Regimes. "Airspeed Definitions & Measurement": Airspeed Measurement at Higher Mach Numbers
- Speed of sound vs. altitude: USSA-1976 1.3.10 Speed of sound and Table III
- Extract: a = √(γRT); T0/T, p0/p, ρ0/ρ as functions of M; ρ changes < 5 % below M ≈ 0.3 (when incompressible Bernoulli is acceptable); compressible pitot relation; the calculator is checked against the NACA 1135 tables and becomes P11.01's starting point.

## 08-18 Incompressible CFD from Scratch

Sources: `@barba-cfd-python-12-steps`, `@mit-ocw-2-29-numerical-fluid-mechanics`, `@nasa-nparc-cfd-vv-tutorial`, `@ghia-ghia-shin-1982-cavity` (optional-paid, data citation only)

- Barba & Forsyth, *CFD Python*: 13_Step_10.ipynb "Step 10: 2D Poisson Equation"; 14_Step_11.ipynb "Step 11: Cavity Flow with Navier–Stokes" (Discretized equations); 15_Step_12.ipynb "Step 12: Channel Flow with Navier–Stokes" (Steps 9 Laplace and 7 Diffusion as needed)
  - https://github.com/barbagroup/CFDPython (loaded); https://raw.githubusercontent.com/barbagroup/CFDPython/master/lessons/14_Step_11.ipynb (loaded)
- MIT 2.29: Lecture 17 "Numerical Methods for the Navier-Stokes Equations" (Choice of Variable Arrangement on the Grid; Calculation of the Pressure; Pressure Correction Methods; A Simple Explicit Scheme); Lecture 19 (Projection Methods – Non-Incremental, Incremental and Rotational-incremental Schemes; Fractional Step Methods; Streamfunction-Vorticity Methods: scheme and boundary conditions; Artificial Compressibility Methods; Boundary Conditions)
  - https://ocw.mit.edu/courses/2-29-numerical-fluid-mechanics-spring-2015/resources/mit2_29s15_lecture17/ (loaded; PDF loaded); .../mit2_29s15_lecture19/ (loaded; PDF loaded)
- NPARC tutorial: "Examining Spatial (Grid) Convergence" (observed order, Richardson extrapolation, GCI)
  - https://www.grc.nasa.gov/www/wind/valid/tutorial/tutorial.html (loaded); .../spatconv.html (pattern)
- Benchmark: Ghia, Ghia & Shin (1982), J. Comput. Phys. 48, 387–411 — paywalled; the lesson embeds the published Re = 100 centreline values as "published data from Ghia et al. 1982, Table I (u on the vertical centreline) and Table II (v on the horizontal centreline)" — https://doi.org/10.1016/0021-9991(82)90058-4 (DOI resolved; article page 403)
- Extract: fractional-step/projection on a staggered or collocated grid with a pressure Poisson solve (P03 linear solvers); CFL and diffusion time-step limits; lid boundary conditions and corner singularities; centreline profiles vs. Ghia at Re = 100 on three grids with observed order and GCI.
- Partial GAP: Step 11 does not compare with Ghia, uses a collocated grid with no discussion of checkerboarding, and stops short of grid-convergence study; the validation workflow is original.

## 08-19 Home Lab: Flow & Pressure Measurement

Sources: `@libretexts-ahmari-kabir-fluid-mechanics-lab-manual`, `@openstax-university-physics-vol1`, `@leishman-intro-aerospace-flight-vehicles`, `@nasa-glenn-beginners-guide-aeronautics`, `@mit-ocw-1-060-engineering-mechanics-ii`

- Ahmari & Kabir: 1.6 "Experiment #6: Orifice and Free Jet Flow" (7. Theory — 7.1 Determination of the Coefficient of Velocity; 7.2 Determination of the Coefficient of Discharge; 8. Experimental Procedure; 9. Results and Calculations); 1.7 "Experiment #7: Osborne Reynolds' Demonstration"; 1.1 "Experiment #1: Hydrostatic Pressure"
  - https://eng.libretexts.org/Bookshelves/Mechanical_Engineering/Applied_Fluid_Mechanics_Lab_Manual_(Ahmari_and_Kabir)/01%3A_Lab_Manual/01.6%3A_Experiment_%236%3A_Orifice_and_Free_Jet_Flow (loaded); .../01.7%3A_Experiment_%237%3A_Osborne_Reynolds'_Demonstration (loaded during contents crawl)
- OpenStax Vol. 1 14.2 Measuring Pressure (Gauge Pressure vs. Absolute Pressure; Manometers; Example "Fluid Heights in an Open U-Tube"); 14.6 Bernoulli's Equation (Analyzing Bernoulli's Equation — the basis of the Torricelli efflux speed); 14.7 Viscosity and Turbulence (Measuring Turbulence — Reynolds number)
- Leishman, "Fluid Statics & the Hydrostatic Equation": Piezometer Tube; U-Tube Manometer. "Applications of the Conservation Laws": Pitot Tube & Pitot-Static Tube. "Mach Number & Reynolds Number": Definition of Reynolds Number
- NASA Glenn "Pitot – Static Tube – Speedometer" (loaded, 08-04); legacy "Reynolds Number" — https://www.grc.nasa.gov/www/k-12/airplane/reynolds.html (loaded)
- MIT 1.060 L4 "Manometry, Stability of Dams and Gates"; L31 "Discharge From Lakes and Under Gates"
- Extract: Torricelli v = Cv√(2gh) and Q = Cd·A√(2gh); draining-time ODE dh/dt = −(Cd·a/A)√(2gh) fitted to timed height readings to estimate Cd (≈ 0.6 for a sharp-edged hole); U-tube Δh → Δp = ρw g Δh → V = √(2Δp/ρair) for a fan with a home-made pitot probe; Re = VD/ν for both experiments; uncertainty propagation from P00.14.
- Partial GAP: all verified procedures use lab benches; the home procedure (bottle or bucket with a drilled orifice, phone-video timing, clear-tubing water manometer, drinking-straw pitot probe on a household fan, safety notes on water and electricity) is original.


## Optional textbooks for deeper study

Paid; never required. Only `@ghia-ghia-shin-1982-cavity` is recorded in `sources.yaml` (as optional-paid), because P08.18 cites its tabulated data; the books below are not recorded.

- Y. A. Çengel and J. M. Cimbala, *Fluid Mechanics: Fundamentals and Applications* (McGraw-Hill) — the most complete worked-example treatment of 08-01 … 08-07 and 08-16, including Moody/Colebrook and minor-loss tables.
- F. M. White, *Fluid Mechanics* and *Viscous Fluid Flow* (McGraw-Hill) — the readings MIT 2.06 and 2.29 assign; *Viscous Fluid Flow* covers exact N-S solutions, Thwaites, Orr–Sommerfeld/e^N and turbulent-BL correlations (08-06, 08-10 … 08-14).
- J. D. Anderson, *Fundamentals of Aerodynamics* (McGraw-Hill) — the text keyed to MIT Unified Engineering's fluids notes; potential flow, vorticity and compressible flow from the aerospace side (08-02, 08-08, 08-09, 08-17).
- S. B. Pope, *Turbulent Flows* (Cambridge University Press) — statistics, spectra, Kolmogorov theory, wall flows and RANS modeling (08-13 … 08-15).
- H. Schlichting and K. Gersten, *Boundary-Layer Theory* (Springer) — the reference for laminar and turbulent boundary layers, stability theory and transition (08-10 … 08-14).
