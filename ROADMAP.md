# Aerospace Engineering from Scratch: Roadmap

Draft 2 (2026-10-11). Every phase and lesson in the curriculum, in teaching order. Draft 2 applies the review in [docs/RECOMMENDATIONS-2026-10-11.md](docs/RECOMMENDATIONS-2026-10-11.md): bridge modules, home experiments from Phase 00 onward, a six-part lesson model, mastery-based assessment, repaired prerequisites, and honest claims about time and equivalence.

**Goal.** Someone who finishes a learning path can do the work expected of a graduate from a top aerospace program: establish requirements, identify the governing physics, make justified assumptions, solve and simulate, design a practical test, evaluate risk and uncertainty, improve the solution, and defend it under technical scrutiny.

> **Learn aerospace engineering from first principles. Design it on your computer. Build it at home. Test it against reality.**

Every major concept connects to something the learner can calculate, simulate, measure, build, or observe. A motivated self-learner with a computer, basic tools, a modest budget, and no university lab must be able to complete every *required* step. Lower the equipment and institutional barriers, never the standard of reasoning, documentation, and evidence.

**How a lesson works.** Each lesson follows a six-part arc: *1. The engineering problem → 2. Understand and derive → 3. Build the computational model → 4. Experiment at home or investigate real data → 5. Validate and improve → 6. Defend your work.* "Build" means writing the method yourself in Python, testing it, characterizing its error, and only then comparing it with a library or professional tool. You write the panel method before you open XFOIL, and the orbit propagator before you touch GMAT. "From scratch" means understanding the assumptions, governing equations, numerical method, and validation; it does not mean rebuilding every production tool. The **Build it** column below names what you end up with. Not every lesson needs hardware, but every technical phase needs applied evidence, and every physical experiment needs a theoretical basis.

**Practical modes.** Each practical activity is labelled as one of: **physical experiment** (you measure a real system), **computational experiment** (you study a model or numerical method, including sensitivity and validation), **data-based investigation** (you analyze authentic measurements collected by others), or **design/build activity** (you fabricate or configure an artifact). They complement each other and are not interchangeable: a simulation does not teach calibration, tolerance, noise, or wiring.

**Competence, not completion.** Each major lesson assesses four things separately: mathematical and physical reasoning, analytical problem solving, computational competence, and engineering interpretation. A passing code test does not pass the lesson. See [Assessment](#assessment).

**How to read the tables.** Lessons inside a phase are sequential. The **Needs** column lists only prerequisites from *other* phases, written `P<phase>.<lesson>` (for example `P02.04` is Phase 2, lesson 4). A blank means the earlier lessons in the same phase are enough.

**Legend:** ✅ Complete · 🚧 In progress · ⬚ Planned · **★** graduate-level depth, optional on every path · 🔬 lesson whose practical activity is a physical home experiment.

---

## Benchmark: what "top-program graduate" means

The scope is the union of three references, so a finished path covers everything each one requires. Three claims are kept separate throughout:

1. **Curriculum scope comparable:** the subject matter maps to the benchmarks below. This is what the tables in this file establish.
2. **Demonstrated competence comparable:** a learner passes independent, demanding assessments and completes validated projects comparable to undergraduate engineering work. This is a goal each learner demonstrates, not a status earned by finishing pages.
3. **Credential equivalence:** none. A self-directed roadmap grants no accredited degree, professional authorization, or engineering license, and never will claim to.

| Reference | What it requires | Where it lands here |
|---|---|---|
| **ABET program criteria, aerospace** | Aeronautical: aerodynamics, aerospace materials, structures, propulsion, flight mechanics, stability and control. Astronautical: orbital mechanics, space environment, attitude determination and control, telecommunications, space structures, rocket propulsion. Design that integrates these topics. | Phases 10–23, design in 24–25 |
| **ABET student outcomes 1–7** | Solve complex problems; design to meet needs; communicate; ethics and professional responsibility; teamwork; experimentation and data analysis; learn new things independently. | Mapped in [Outcome coverage](#outcome-coverage) |
| **MIT Course 16 (benchmark degree)** | Core: programming, Unified Engineering (materials and structures, signals and systems, fluids, thermo and propulsion), dynamics, automatic control, probability and statistics, differential equations. Four professional areas from aerodynamics, structures, propulsion, computation, feedback control, real-time software, communications, human systems, autonomy. A design capstone (aircraft, spacecraft, or autonomous vehicle) and a build-and-test lab. | Core = Phases 0–9; professional areas = 10–23; capstone and lab = 23, 26 |

To be cross-checked against the current ABET criteria and the Stanford, Caltech GALCIT, Georgia Tech, and Purdue course listings during step 2, when each phase gets its reference list; the mapping is versioned and re-verified rather than treated as permanent. Graduate-level depth (first-year master's) is included where top undergraduates routinely reach it, marked **★**. No mandatory lesson depends on a ★ lesson.

---

## Phase overview

| Phase | Title | Lessons | Est. hours |
|---|---|---|---|
| B1 | Bridge: Precalculus & Trigonometry (optional, placement-tested) | 8 | 20 |
| B2 | Bridge: General Chemistry (optional, placement-tested) | 6 | 15 |
| 00 | Setup, Engineering Computing & First Measurements ✅ | 15 | 24 |
| 01 | Mathematics I: Calculus & Linear Algebra ✅ | 18 | 45 |
| 02 | Mathematics II: Differential Equations, Transforms & Probability ✅ | 20 | 55 |
| 03 | Numerical Methods & Scientific Computing ✅ | 18 | 45 |
| 04 | Physics Foundations ✅ | 15 | 38 |
| 05 | Statics, Mechanics of Materials & CAD ✅ | 18 | 46 |
| 06 | Dynamics & Vibrations ✅ | 19 | 48 |
| 07 | Thermodynamics ✅ | 15 | 38 |
| 08 | Fluid Mechanics ✅ | 19 | 48 |
| 09 | Signals, Systems & Control ✅ | 21 | 53 |
| 10 | Aerodynamics | 23 | 64 |
| 11 | Compressible Flow, Hypersonics & CFD | 18 | 50 |
| 12 | Aerospace Materials, Manufacturing & Drawings | 16 | 35 |
| 13 | Aerospace Structures | 21 | 59 |
| 14 | Air-Breathing Propulsion | 16 | 40 |
| 15 | Rocket & Space Propulsion | 16 | 40 |
| 16 | Aircraft Performance | 14 | 30 |
| 17 | Flight Dynamics, Stability & Control | 19 | 50 |
| 18 | Estimation, Navigation, Guidance & Autonomy | 20 | 55 |
| 19 | Orbital Mechanics & Astrodynamics | 20 | 50 |
| 20 | Spacecraft Systems & the Space Environment | 21 | 54 |
| 21 | Avionics, Embedded & Flight Software | 17 | 44 |
| 22 | Experimental Methods & Testing | 14 | 35 |
| 23 | Systems Engineering & Professional Practice | 16 | 30 |
| 24 | Aircraft Design | 14 | 45 |
| 25 | Spacecraft & Mission Design | 14 | 45 |
| 26 | Capstone Projects | 10 | 200+ |
| | **Total (excluding bridge modules)** | **467** | **~1,170 lesson hours + capstones** |

**About the hours.** The estimates above cover only reading and working through the lesson itself, and none of them has been validated with real learners. A full path also includes derivations and problem sets, programming and debugging, equipment preparation and experiments, data analysis and reporting, exams and remediation, and the capstones. Each of those is tracked as its own estimate once pilot learners have measured it; until then this file deliberately gives no total, and the standard for finishing is demonstrated competence, not hours spent. A novice will take far longer than an experienced programmer on the same assignment, so published estimates will be ranges by starting level.

Every phase ends with a **phase exam** (an integrated test with a suggested timed, closed-reference self-assessment mode) and a **mini-project** (a rubric-graded engineering investigation that integrates several lessons). From the first mini-project onward, every project states requirements, margins, risks, and a verification plan; Phase 23 formalizes what the learner has been practising since Phase 00.

---

## Bridge modules (optional, before Phase 01) ⬚

A short diagnostic test at the start of Phase 00 decides whether a learner needs these. They are never assumed: a true beginner is not silently expected to know algebra, trigonometry, or chemistry.

### B1: Precalculus & Trigonometry ⬚ (~20 h)

| # | Lesson | Build it | Needs |
|---|---|---|---|
| 01 | Algebra & Units Refresher | Unit-consistent rearrangement drills | |
| 02 | Functions & Graphs | Function plotter | |
| 03 | Exponentials & Logarithms | Atmospheric pressure decay | |
| 04 | Trigonometry & the Unit Circle | Angle-of-attack geometry | |
| 05 | Vectors in the Plane | Force resolution | |
| 06 | Geometry & Coordinate Systems | Wing planform area and MAC | |
| 07 | Complex Numbers Intro | Rotation by multiplication | |
| 08 | Bridge Exam & Placement | Pass to enter Phase 01 | |

### B2: General Chemistry ⬚ (~15 h)

| # | Lesson | Build it | Needs |
|---|---|---|---|
| 01 | Atoms, Moles & Stoichiometry | Propellant mixture ratio calculator | |
| 02 | Bonding & Molecular Structure | Bond-energy estimator | |
| 03 | Gases & Solutions | Ideal gas law lab with a syringe 🔬 | |
| 04 | Thermochemistry | Heat of combustion from bond energies | |
| 05 | Reaction Rates & Equilibrium Intro | Simple rate model | |
| 06 | Electrochemistry & Corrosion Basics | Battery cell voltage, galvanic series | |

## Phase 00: Setup, Engineering Computing & First Measurements ✅ (~24 h)

| # | Lesson | Build it | Needs |
|---|---|---|---|
| 01 | Dev Environment for Engineers | Python, uv, VS Code, Jupyter set up and verified | |
| 02 | Git & Collaboration | Fork, branch, PR on this repo | |
| 03 | NumPy as an Engineering Calculator | Vectorized ISA atmosphere table | |
| 04 | Plotting That Engineers Trust | Publication-quality plot style module | |
| 05 | Units & Dimensional Consistency | Unit-checked calculations with `pint`; Mars Climate Orbiter post-mortem | |
| 06 | SciPy Tour | Integrate, optimize, interpolate on a real dataset | |
| 07 | Testing Engineering Code | pytest suite with tolerance-based assertions | |
| 08 | Reading Data: CSV, HDF5, Telemetry Logs | Parser for a flight log | |
| 09 | Symbolic Math with SymPy | Derive and lambdify an equation of motion | |
| 10 | Writing Technical Reports | Report template (Markdown → PDF via LaTeX) | |
| 11 | Reading a Technical Paper | Annotated summary of a NASA technical report | |
| 12 | Engineering Judgment & Sanity Checks | Order-of-magnitude checker for an airliner | |
| 13 | Safety & Standards for Home Experiments | Hazard checklist and stop-criteria template used by every 🔬 lesson | |
| 14 | Measurement, Calibration & Uncertainty | Calibrate a kitchen scale and a phone sensor; report uncertainty 🔬 | |
| 15 | First Experiment: Predict, Measure, Compare | Time a falling object with phone video; compare with prediction 🔬 | |

## Phase 01: Mathematics I: Calculus & Linear Algebra ✅ (~45 h)

| # | Lesson | Build it | Needs |
|---|---|---|---|
| 01 | Limits, Continuity & Derivatives | Numerical vs symbolic derivative comparison | P00.03 |
| 02 | Applications of Derivatives | Max range glide from a drag polar | |
| 03 | Integrals & the Fundamental Theorem | Area and centroid of an airfoil section | |
| 04 | Techniques of Integration | Moment of inertia of a tapered beam | |
| 05 | Sequences, Series & Taylor Expansions | Small-angle approximation error explorer | |
| 06 | Parametric Curves & Polar Coordinates | Plot a conic orbit from r(θ) | |
| 07 | Vectors, Dot & Cross Products | Vector library for 3D mechanics | |
| 08 | Multivariable Functions & Partial Derivatives | Sensitivity of lift to α and V | |
| 09 | Gradients, Directional Derivatives & Optimization | Minimum-drag speed by gradient search | |
| 10 | Multiple Integrals | Mass properties of a wing box | |
| 11 | Vector Fields, Line & Surface Integrals | Circulation around a vortex | |
| 12 | Divergence, Stokes & Green's Theorems | Verify divergence theorem numerically | |
| 13 | Systems of Linear Equations | Gaussian elimination from scratch | |
| 14 | Vector Spaces, Basis & Rank | Detect a mechanism in a truss | |
| 15 | Matrix Transformations & Rotation Matrices | Direction cosine matrix tools | |
| 16 | Eigenvalues & Eigenvectors | Principal axes of an inertia tensor | |
| 17 | Orthogonality, Least Squares & QR | Fit a drag polar to wind tunnel data | |
| 18 | SVD & Matrix Norms | Conditioning of a structural stiffness matrix | |

## Phase 02: Mathematics II: Differential Equations, Transforms & Probability ✅ (~55 h)

| # | Lesson | Build it | Needs |
|---|---|---|---|
| 01 | First-Order ODEs | Rocket mass depletion model | P01.04 |
| 02 | Second-Order Linear ODEs | Damped spring-mass solver | |
| 03 | Forced Response & Resonance | Resonance explorer | |
| 04 | Systems of ODEs & Phase Portraits | Phase portrait plotter | P01.16 |
| 05 | Matrix Exponential & Linear Systems | State transition matrix | |
| 06 | Nonlinear ODEs, Equilibria & Linearization | Pendulum linearization | |
| 07 | Complex Numbers & Phasors | Phasor calculator | |
| 08 | Laplace Transform | Transfer-function solver for ODEs | |
| 09 | Fourier Series | Square-wave synthesis | |
| 10 | Fourier Transform & the FFT | FFT of vibration data | |
| 11 | PDEs: Heat, Wave & Laplace Equations | Separation-of-variables solutions | P01.12 |
| 12 | Complex Variables & Conformal Maps | Joukowski airfoil from a circle | |
| 13 | Calculus of Variations | Brachistochrone and minimum-time paths | |
| 14 | Tensor Notation & Index Gymnastics | Stress tensor transformations | P01.15 |
| 15 | Probability Fundamentals | Monte Carlo dice-to-dispersion | |
| 16 | Random Variables & Distributions | Landing dispersion ellipse | |
| 17 | Statistics & Hypothesis Testing | Is this material batch in spec? | |
| 18 | Regression & Curve Fitting | Thrust curve fit with confidence bands | P01.17 |
| 19 | Random Processes & Power Spectra | Gust spectrum (Dryden) generator | |
| 20 | Bayesian Inference | Update a failure-rate estimate | |

## Phase 03: Numerical Methods & Scientific Computing ✅ (~45 h)

| # | Lesson | Build it | Needs |
|---|---|---|---|
| 01 | Floating Point & Numerical Error | Catastrophic cancellation demo | P00.03 |
| 02 | Root Finding | Newton, secant, bisection; solve Kepler's equation | P01.05 |
| 03 | Interpolation & Splines | Atmospheric table interpolator | |
| 04 | Numerical Differentiation & Integration | Simpson, Gauss quadrature | |
| 05 | ODE Integrators I: Euler to RK4 | RK4 from scratch | P02.04 |
| 06 | ODE Integrators II: Adaptive & Stiff | RK45 with step control, implicit Euler | |
| 07 | Symplectic Integrators | Energy-conserving orbit propagation | |
| 08 | Direct Linear Solvers | LU, Cholesky from scratch | P01.13 |
| 09 | Iterative & Sparse Solvers | Conjugate gradient on a 2D Poisson problem | |
| 10 | Nonlinear Systems | Newton-Raphson for trim | |
| 11 | Finite Differences for PDEs | 1D heat equation, stability analysis | P02.11 |
| 12 | Finite Volume Methods | 1D advection with upwinding | |
| 13 | Finite Element Basics | 1D bar FEM solver | |
| 14 | Unconstrained Optimization | BFGS from scratch | P01.09 |
| 15 | Constrained Optimization | SQP for a sizing problem | |
| 16 | Multidisciplinary Design Optimization | Wing aero-structural MDO with OpenMDAO | |
| 17 | Uncertainty Quantification | Monte Carlo and polynomial chaos | P02.16 |
| 18 | Verification & Validation | Method of manufactured solutions | |

## Phase 04: Physics Foundations ✅ (~35 h)

| # | Lesson | Build it | Needs |
|---|---|---|---|
| 01 | Kinematics & Newton's Laws | Projectile with drag | P01.07 |
| 02 | Work, Energy & Power | Energy-height of an aircraft | |
| 03 | Momentum & Collisions | Variable-mass system intuition | |
| 04 | Rotation & Angular Momentum | Spinning top simulation | |
| 05 | Gravitation & Kepler's Laws | Two-body simulation | P03.05 |
| 06 | Oscillations & Waves | Wave propagation visualizer | P02.02 |
| 07 | Kinetic Theory of Gases | Molecular speed distributions | |
| 08 | Electrostatics & Circuits | DC circuit solver | P01.13 |
| 09 | Magnetism & Induction | Magnetorquer field model | |
| 10 | AC Circuits & Electronics Basics | RC filter design | P02.07 |
| 11 | Electromagnetic Waves & Antennas | Dipole radiation pattern | |
| 12 | Optics & Sensors | Star tracker optics budget | |
| 13 | Special Relativity for GPS | Clock corrections in GNSS | |
| 14 | Atomic & Radiation Physics | Radiation dose basics | |
| 15 | Electronics Bench Lab | Breadboard, multimeter, microcontroller: measure voltage, current, and a sensor 🔬 | P00.14 |

## Phase 05: Statics & Mechanics of Materials ✅ (~40 h)

| # | Lesson | Build it | Needs |
|---|---|---|---|
| 01 | Forces, Moments & Free-Body Diagrams | FBD solver for a landing gear | P04.01 |
| 02 | Equilibrium of Rigid Bodies | Reaction solver | P01.13 |
| 03 | Trusses: Joints & Sections | Truss solver from scratch | |
| 04 | Frames, Machines & Friction | Control linkage analysis | |
| 05 | Centroids & Moments of Inertia | Section property calculator | P01.10 |
| 06 | Stress & Strain | Stress-strain curve analyzer | |
| 07 | Axial Loading & Thermal Stress | Bolted joint preload | |
| 08 | Torsion of Shafts | Turbine shaft sizing | |
| 09 | Shear & Bending Moment Diagrams | V–M diagram generator | |
| 10 | Bending Stress in Beams | Wing spar stress | |
| 11 | Shear Stress in Beams | Shear flow in an I-beam | |
| 12 | Beam Deflection | Deflection by integration and superposition | P02.02 |
| 13 | Stress Transformation & Mohr's Circle | Mohr's circle tool | P02.14 |
| 14 | Failure Theories | von Mises, Tresca margin calculator | |
| 15 | Column Buckling | Euler and Johnson buckling | |
| 16 | Energy Methods | Castigliano deflection solver | |
| 17 | CAD & Technical Drawings I | Model and dimension a bracket; produce a drawing | |
| 18 | Home Lab: Cantilever Beam Deflection | Load a clamped beam, measure deflection, compare with theory 🔬 | P00.15 |

## Phase 06: Dynamics & Vibrations ✅ (~45 h)

| # | Lesson | Build it | Needs |
|---|---|---|---|
| 01 | Particle Kinematics in Many Coordinates | Polar, cylindrical, path coordinates | P04.01 |
| 02 | Particle Kinetics | Pull-up maneuver loads | |
| 03 | Work-Energy & Impulse-Momentum | Arrestor-hook landing | |
| 04 | Relative Motion & Rotating Frames | Coriolis on a rotating Earth | P01.15 |
| 05 | Rigid Body Kinematics in 2D | Four-bar linkage | |
| 06 | Rigid Body Kinetics in 2D | Rolling, sliding wheel | |
| 07 | 3D Rotations: Euler Angles & DCMs | Attitude kinematics library | P01.15 |
| 08 | Quaternions | Quaternion library from scratch | |
| 09 | 3D Rigid Body Dynamics: Euler's Equations | Tumbling body simulator (Dzhanibekov effect) | P03.05 |
| 10 | Gyroscopes & Spin Stability | Spin-stabilized satellite | |
| 11 | Lagrangian Mechanics | Double pendulum from Lagrange | P02.13 |
| 12 | Hamiltonian Mechanics ★ | Phase-space orbits | |
| 13 | Single-DOF Vibration | Free and forced response | P02.03 |
| 14 | Multi-DOF Vibration & Modes | Modal analysis of a 3-mass system | P01.16 |
| 15 | Continuous Systems | Beam bending modes | P02.11 |
| 16 | Random Vibration | Launch vibration response (PSD to Grms) | P02.19 |
| 17 | Vibration Isolation | Isolator design for an IMU | |
| 18 | Multibody Simulation | Landing gear multibody sim | |
| 19 | Home Lab: Pendulum System Identification | Log a pendulum with phone or IMU; fit damping and frequency 🔬 | P02.18 |

## Phase 07: Thermodynamics ✅ (~35 h)

| # | Lesson | Build it | Needs |
|---|---|---|---|
| 01 | Systems, Properties & State | Ideal-gas property calculator | P04.07 |
| 02 | The First Law | Tank blowdown model | P02.01 |
| 03 | Properties of Pure Substances | Steam/propellant table lookup | P03.03 |
| 04 | The Second Law & Entropy | Entropy generation tracker | |
| 05 | Isentropic Processes & Efficiencies | Compressor efficiency calculator | |
| 06 | Power Cycles: Brayton | Ideal Brayton cycle analyzer | |
| 07 | Real Brayton Cycles | Intercooled, recuperated cycle | |
| 08 | Exergy | Exergy loss audit of an engine | |
| 09 | Gas Mixtures & Psychrometrics | Cabin air conditioning | |
| 10 | Chemical Reactions & Combustion | Adiabatic flame temperature | |
| 11 | Chemical Equilibrium | Equilibrium solver (Gibbs minimization) | P03.15 |
| 12 | Heat Transfer: Conduction | 2D conduction solver | P03.11 |
| 13 | Heat Transfer: Convection | Convective cooling correlations | |
| 14 | Heat Transfer: Radiation | View factor and radiative balance | |
| 15 | Home Lab: Cooling Curve | Measure cooling of a small object; fit lumped and distributed models 🔬 | P00.14 |

## Phase 08: Fluid Mechanics ✅ (~45 h)

| # | Lesson | Build it | Needs |
|---|---|---|---|
| 01 | Fluid Properties & Hydrostatics | Standard atmosphere from scratch | P07.01 |
| 02 | Fluid Kinematics | Streamline, pathline, streakline plotter | P01.11 |
| 03 | Control Volume: Mass & Momentum | Jet thrust from momentum balance | |
| 04 | Bernoulli & Energy Equation | Pitot-static airspeed | |
| 05 | Dimensional Analysis & Similarity | Buckingham Pi tool; scaling a wind tunnel model | |
| 06 | Navier-Stokes Equations | Derivation and exact solutions | P02.11 |
| 07 | Pipe Flow & Losses | Propellant feed line pressure drop | |
| 08 | Potential Flow | Superposition of elementary flows | P02.12 |
| 09 | Vorticity & Circulation | Kelvin's theorem demo | |
| 10 | Boundary Layers: Laminar | Blasius solution | P03.05 |
| 11 | Integral Boundary Layer Methods | Thwaites method | |
| 12 | Transition | e^N method intro | |
| 13 | Turbulence Fundamentals | Turbulence statistics from data | P02.19 |
| 14 | Turbulent Boundary Layers | Log-law and skin friction | |
| 15 | Turbulence Modeling | Mixing-length and k-ε overview | |
| 16 | Flow Separation & Drag | Drag of bluff bodies | |
| 17 | Intro to Compressible Flow | Speed of sound and Mach number | P07.05 |
| 18 | Incompressible CFD from Scratch | 2D lid-driven cavity solver | P03.12 |
| 19 | Home Lab: Flow & Pressure Measurement | Low-energy flow rig; Reynolds reasoning and calibrated pressure 🔬 | P00.14 |

## Phase 09: Signals, Systems & Control ✅ (~50 h)

| # | Lesson | Build it | Needs |
|---|---|---|---|
| 01 | Signals & Systems | Signal toolkit | P02.08 |
| 02 | Convolution & Impulse Response | Discrete convolution | |
| 03 | Transfer Functions & Block Diagrams | Block diagram reducer | |
| 04 | Time Response & Poles | Step response analyzer | |
| 05 | Frequency Response & Bode Plots | Bode plotter from scratch | P02.10 |
| 06 | Sampling & Aliasing | Anti-aliasing filter design | |
| 07 | Digital Filters | FIR and IIR from scratch | |
| 08 | Feedback & Steady-State Error | Cruise controller | |
| 09 | Root Locus | Root locus plotter | |
| 10 | Nyquist & Stability Margins | Nyquist plotter, gain and phase margins | |
| 11 | PID Design & Tuning | Pitch attitude PID with anti-windup | |
| 12 | Lead-Lag Compensation | Loop shaping a servo | |
| 13 | State-Space Models | State-space from ODEs | P02.05 |
| 14 | Controllability & Observability | Rank tests | P01.14 |
| 15 | Pole Placement & Observers | Full-state feedback + Luenberger observer | |
| 16 | LQR | Riccati solver, LQR design | |
| 17 | Digital Control | Discretization, z-transform | |
| 18 | Robustness & Uncertainty ★ | Sensitivity functions, H∞ intro | |
| 19 | Nonlinear Control Basics | Lyapunov, feedback linearization | P02.06 |
| 20 | Model Predictive Control ★ | Linear MPC from scratch | P03.15 |
| 21 | Home Lab: Servo Feedback Control | Close the loop on a servo or small motor; measure transient response 🔬 | P04.15 |

## Phase 10: Aerodynamics ⬚ (~60 h)

| # | Lesson | Build it | Needs |
|---|---|---|---|
| 01 | Forces, Moments & Coefficients | Coefficient calculator from balance data | P08.05 |
| 02 | Airfoil Geometry | NACA 4/5-digit generator | |
| 03 | The Kutta Condition & Circulation | Kutta-Joukowski lift | P08.08 |
| 04 | Thin Airfoil Theory | Thin airfoil solver (camber line) | P02.09 |
| 05 | Vortex Panel Method | 2D panel method from scratch | P01.13 |
| 06 | Viscous-Inviscid Coupling | Panel method + Thwaites/Head boundary layer | P08.11 |
| 07 | Airfoil Design & Selection | Inverse design exercise | |
| 08 | High-Lift Devices | Flap and slat effects | |
| 09 | Finite Wings & Induced Drag | Downwash and induced drag | |
| 10 | Prandtl Lifting-Line Theory | Lifting-line solver | |
| 11 | Vortex Lattice Method | VLM from scratch | P03.08 |
| 12 | Wing Planform Effects | Taper, sweep, twist optimizer | P03.14 |
| 13 | Drag Build-Up | Component drag estimation | P08.14 |
| 14 | Low Reynolds Number Aerodynamics | Drone and glider airfoils | P08.12 |
| 15 | Stall & Post-Stall | Stall progression model | |
| 16 | Fuselage & Body Aerodynamics | Slender body theory | |
| 17 | Propeller & Rotor Theory | Blade element momentum code | |
| 18 | Helicopter Aerodynamics | Hover and forward flight performance | |
| 19 | Unsteady Aerodynamics | Theodorsen, Wagner functions | P02.08 |
| 20 | Ground Effect & Interference | Wing-in-ground-effect estimate | |
| 21 | Using XFOIL, AVL & OpenVSP | Validate your codes against the tools | |
| 22 | Aerodynamic Validation Project | Compare your VLM to wind tunnel data | |
| 23 | Home Lab: Airfoil Lift & Drag on a Bench Rig | Guarded tabletop airflow rig and simple balance; theory vs computation vs measurement 🔬 | P08.19 |

## Phase 11: Compressible Flow, Hypersonics & CFD ⬚ (~50 h)

| # | Lesson | Build it | Needs |
|---|---|---|---|
| 01 | Isentropic Flow Relations | Compressible flow calculator | P08.17 |
| 02 | Normal Shocks | Normal shock solver | |
| 03 | Oblique Shocks | θ-β-M solver | P03.02 |
| 04 | Prandtl-Meyer Expansion | Expansion fan calculator | |
| 05 | Shock-Expansion Theory | Supersonic diamond airfoil | |
| 06 | Quasi-1D Nozzle Flow | Converging-diverging nozzle solver | |
| 07 | Method of Characteristics | Minimum-length nozzle design | |
| 08 | Linearized Supersonic Flow | Ackeret theory | |
| 09 | Transonic Flow & Supercritical Airfoils | Critical Mach estimator | |
| 10 | Swept Wings & Area Rule | Wave drag via area rule | |
| 11 | Fanno & Rayleigh Flow | Duct flow with friction and heating | |
| 12 | Hypersonic Flow Fundamentals | Newtonian impact theory | |
| 13 | Aerothermodynamics & Heating | Stagnation heating (Fay-Riddell, Sutton-Graves) | P07.13 |
| 14 | Real-Gas & High-Temperature Effects ★ | Equilibrium air model | P07.11 |
| 15 | Euler Equations & Riemann Problems | Godunov solver for shock tube | P03.12 |
| 16 | 2D Compressible CFD | 2D Euler solver (supersonic wedge) | |
| 17 | RANS CFD in Practice | Run SU2 on an airfoil, grid convergence | P03.18 |
| 18 | Mesh Generation | Structured mesh generator | |

## Phase 12: Aerospace Materials & Manufacturing ⬚ (~30 h)

| # | Lesson | Build it | Needs |
|---|---|---|---|
| 01 | Atomic Structure & Bonding | Crystal structure visualizer | P04.14 |
| 02 | Mechanical Properties & Testing | Tensile test data reducer | P05.06 |
| 03 | Aluminum & Titanium Alloys | Material selection database | |
| 04 | Steels & Superalloys | Creep life estimate (Larson-Miller) | |
| 05 | Material Selection (Ashby Methods) | Ashby chart from scratch | P03.14 |
| 06 | Polymers & Composites Overview | Rule of mixtures | |
| 07 | Ceramics & Thermal Protection Materials | TPS material trade | |
| 08 | Fatigue | S-N curve and Miner's rule | P02.16 |
| 09 | Fracture Mechanics | Crack growth (Paris law) | |
| 10 | Corrosion & Environmental Degradation | Corrosion allowance planner | |
| 11 | Conventional Manufacturing | Machining cost and tolerance stack-up | |
| 12 | Composite Manufacturing | Layup, cure cycle | |
| 13 | Additive Manufacturing | AM design rules checker | |
| 14 | Joining: Fasteners, Bonding, Welding | Fastener joint analysis | |
| 15 | GD&T & Tolerance Analysis | Tolerance stack-up for a fitted assembly | P05.17 |
| 16 | CAD II: Assemblies & Drawings for Manufacture | Assembly with a bill of materials and inspection drawing | |

## Phase 13: Aerospace Structures ⬚ (~55 h)

| # | Lesson | Build it | Needs |
|---|---|---|---|
| 01 | Loads on Aircraft & Spacecraft | V-n diagram, launch load factors | P06.02 |
| 02 | Elasticity Fundamentals | 2D elasticity solutions | P02.14 |
| 03 | Thin-Walled Beams: Bending | Unsymmetric bending of a stringer section | P05.10 |
| 04 | Shear Flow in Open Sections | Shear center finder | |
| 05 | Shear Flow in Closed & Multi-Cell Sections | Wing box torsion solver | |
| 06 | Idealized Structures | Boom idealization of a wing box | |
| 07 | Wing & Fuselage Structural Analysis | Wing box sizing tool | |
| 08 | Plates & Shells | Plate bending solver | P02.11 |
| 09 | Buckling of Plates & Stiffened Panels | Panel buckling margins | P05.15 |
| 10 | Matrix Structural Analysis | Direct stiffness method (frames) | P03.08 |
| 11 | Finite Element Method: 2D/3D | 2D plane-stress FEM from scratch | P03.13 |
| 12 | FEM in Practice | CalculiX validation; mesh convergence | |
| 13 | Composite Laminate Theory | Classical lamination theory code | P12.06 |
| 14 | Composite Failure Criteria | Tsai-Wu, Hashin failure envelopes | |
| 15 | Damage Tolerance & Fail-Safe Design | Inspection interval planning | P12.09 |
| 16 | Structural Dynamics | Modal FEM of a wing | P06.14 |
| 17 | Static Aeroelasticity | Divergence and control reversal | P10.10 |
| 18 | Dynamic Aeroelasticity: Flutter | p-k flutter solver | P10.19 |
| 19 | Thermal Structures | Thermal stress in a hot structure | P07.12 |
| 20 | Margins of Safety & Stress Reports | Formal stress report | |
| 21 | Home Lab: Instrumented Structure | Strain-gauged built-up beam or coupon; stiffness, failure mode, margins 🔬 | P05.18 |

## Phase 14: Air-Breathing Propulsion ⬚ (~40 h)

| # | Lesson | Build it | Needs |
|---|---|---|---|
| 01 | Thrust, Efficiency & Specific Fuel Consumption | Thrust equation explorer | P08.03 |
| 02 | Ideal Ramjet | Ramjet cycle | P11.01 |
| 03 | Ideal Turbojet | Turbojet cycle | P07.06 |
| 04 | Turbofans | Bypass ratio trade | |
| 05 | Turboprops & Turboshafts | Turboprop cycle | |
| 06 | Real Cycle Analysis | Component-efficiency cycle deck | |
| 07 | Off-Design Performance | Component matching | P03.10 |
| 08 | Inlets & Diffusers | Supersonic inlet design | P11.03 |
| 09 | Axial Compressors | Velocity triangles, stage loading | |
| 10 | Centrifugal Compressors | Impeller sizing | |
| 11 | Combustors | Combustor sizing, emissions | P07.10 |
| 12 | Turbines & Cooling | Turbine stage and blade cooling | |
| 13 | Nozzles | Nozzle area scheduling | P11.06 |
| 14 | Propellers & Electric Propulsion for Aircraft | Motor-prop matching | P10.17 |
| 15 | Scramjets ★ | Scramjet cycle | |
| 16 | Engine Integration & Noise | Installation losses | |

## Phase 15: Rocket & Space Propulsion ⬚ (~40 h)

| # | Lesson | Build it | Needs |
|---|---|---|---|
| 01 | The Rocket Equation & Staging | Staging optimizer | P04.03 |
| 02 | Thrust, Isp & Nozzle Performance | Ideal rocket performance | P11.06 |
| 03 | Thermochemistry of Propellants | Mini-CEA equilibrium code | P07.11 |
| 04 | Nozzle Design | Bell nozzle contour (Rao) | P11.07 |
| 05 | Liquid Rocket Engines | Engine cycle comparison | |
| 06 | Feed Systems & Turbopumps | Pressure-fed vs pump-fed sizing | P08.07 |
| 07 | Injectors & Combustion Instability | Injector design | |
| 08 | Thrust Chamber Cooling | Regenerative cooling model | P07.13 |
| 09 | Solid Rocket Motors | Grain regression simulator | |
| 10 | Hybrid Rockets | Hybrid motor model | |
| 11 | Electric Propulsion I: Electrothermal & Electrostatic | Ion thruster performance | P04.08 |
| 12 | Electric Propulsion II: Hall & Electromagnetic | Hall thruster scaling | |
| 13 | Cold Gas & Monopropellant Systems | RCS thruster sizing | |
| 14 | Propulsion System Design | Full propulsion system budget | |
| 15 | Test Stands & Static Fires | Test data reduction | |
| 16 | Advanced Concepts ★ | Nuclear thermal, solar sail estimates | |

## Phase 16: Aircraft Performance ⬚ (~30 h)

| # | Lesson | Build it | Needs |
|---|---|---|---|
| 01 | The Atmosphere & Airspeeds | IAS/CAS/EAS/TAS converter | P08.04 |
| 02 | Drag Polars & Lift-to-Drag | Polar builder | P10.13 |
| 03 | Steady Level Flight | Thrust/power required curves | |
| 04 | Climb & Descent | Rate of climb, ceiling | |
| 05 | Gliding Flight | Glide polar, best glide | |
| 06 | Range & Endurance: Jets | Breguet range | |
| 07 | Range & Endurance: Props & Electric | Battery-electric range | P14.14 |
| 08 | Takeoff & Landing | Field length calculator | |
| 09 | Turning Flight & Maneuvers | Turn performance, V-n | |
| 10 | Energy Methods | Specific excess power, min time to climb | |
| 11 | Mission Analysis | Mission segment simulator | |
| 12 | Payload-Range Diagrams | Airliner payload-range | |
| 13 | Rotorcraft Performance | Helicopter power curve | P10.18 |
| 14 | Performance Flight Test Reduction | Reduce real flight test data | |

## Phase 17: Flight Dynamics, Stability & Control ⬚ (~45 h)

| # | Lesson | Build it | Needs |
|---|---|---|---|
| 01 | Reference Frames for Flight | Frame transform library | P06.07 |
| 02 | 6-DOF Equations of Motion | 6-DOF aircraft simulator | P06.09 |
| 03 | Static Longitudinal Stability | Neutral point, static margin | P10.12 |
| 04 | Trim & Control Power | Trim solver | P03.10 |
| 05 | Static Lateral-Directional Stability | Dihedral and weathercock | |
| 06 | Stability Derivatives | Estimate derivatives with your VLM | P10.11 |
| 07 | Linearization & Small Perturbations | Linear models from the 6-DOF sim | P02.06 |
| 08 | Longitudinal Modes | Short period and phugoid | P09.13 |
| 09 | Lateral-Directional Modes | Roll, spiral, Dutch roll | |
| 10 | Handling Qualities | Cooper-Harper and MIL-F-8785 checks | |
| 11 | Stability Augmentation | Yaw damper design | P09.10 |
| 12 | Autopilots | Altitude and heading hold | P09.11 |
| 13 | Modern Flight Control | LQR-based inner loop | P09.16 |
| 14 | Fly-by-Wire & Nonlinear Control | Dynamic inversion | P09.19 |
| 15 | Atmospheric Disturbances | Dryden/von Kármán turbulence in the sim | P02.19 |
| 16 | Multirotor Dynamics & Control | Quadrotor sim and controller | |
| 17 | Spin, Departure & High-α ★ | Departure prediction | |
| 18 | System Identification from Flight Data | Output-error parameter estimation | P02.18 |
| 19 | Home Lab: Instrumented Glider Flight Test | Build, trim, and fly a foam glider with a logger; glide and stability data 🔬 | P10.23, P16.05 |

## Phase 18: Estimation, Navigation, Guidance & Autonomy ⬚ (~55 h)

| # | Lesson | Build it | Needs |
|---|---|---|---|
| 01 | Estimation Fundamentals | Least-squares estimator | P01.17 |
| 02 | Recursive Estimation | Recursive least squares | P02.16 |
| 03 | The Kalman Filter | Kalman filter from scratch | P09.13 |
| 04 | Extended Kalman Filter | EKF for a falling body | |
| 05 | Unscented & Particle Filters | UKF and particle filter | |
| 06 | Inertial Sensors | IMU error model and simulator | P06.07 |
| 07 | Strapdown Inertial Navigation | INS mechanization | |
| 08 | GNSS Positioning | Pseudorange solver | P04.13 |
| 09 | GNSS/INS Integration | Loosely coupled GNSS/INS EKF | |
| 10 | Attitude Estimation | Complementary filter, MEKF | P06.08 |
| 11 | Air Data & Other Sensors | Air data computer | |
| 12 | Vision-Based Navigation | Visual odometry basics | |
| 13 | Guidance Laws | Proportional navigation | |
| 14 | Trajectory Optimization | Direct collocation solver | P03.15 |
| 15 | Optimal Control Theory ★ | Pontryagin, shooting methods | P02.13 |
| 16 | Powered Descent Guidance ★ | Convex-optimization landing guidance | |
| 17 | Path Planning | A*, RRT for a UAV | |
| 18 | Decision Making Under Uncertainty | MDPs for autonomy | |
| 19 | Sensor Fusion Architectures | Fault detection and isolation | |
| 20 | GNC Integration Project | Closed-loop GNC on a simulated UAV or spacecraft | P17.16 or P20.08 |

## Phase 19: Orbital Mechanics & Astrodynamics ⬚ (~50 h)

| # | Lesson | Build it | Needs |
|---|---|---|---|
| 01 | The Two-Body Problem | Two-body propagator | P04.05 |
| 02 | Conic Sections & Orbit Geometry | Orbit plotter | P01.06 |
| 03 | Orbital Elements | State vector ↔ elements converter | |
| 04 | Kepler's Equation & Time of Flight | Kepler solver | P03.02 |
| 05 | Coordinate Systems & Time | ECI/ECEF/topocentric, time scales | P06.07 |
| 06 | Ground Tracks & Visibility | Ground track and pass predictor | |
| 07 | Orbital Maneuvers | Hohmann, bi-elliptic, plane changes | |
| 08 | Lambert's Problem | Lambert solver | |
| 09 | Orbit Perturbations | J2, drag, SRP, third-body | P03.06 |
| 10 | Sun-Synchronous & Special Orbits | Orbit design tool | |
| 11 | Relative Motion & Rendezvous | Clohessy-Wiltshire, rendezvous planner | |
| 12 | Interplanetary Transfers | Patched conics, porkchop plots | |
| 13 | Gravity Assists | Flyby calculator | |
| 14 | The Three-Body Problem | CR3BP, Lagrange points | |
| 15 | Halo Orbits & Low-Energy Transfers ★ | Differential correction for halo orbits | |
| 16 | Launch Trajectories | Ascent trajectory simulator | P15.01 |
| 17 | Atmospheric Entry | Entry trajectory and heating | P11.13 |
| 18 | Low-Thrust Trajectories | Low-thrust spiral | P15.11 |
| 19 | Orbit Determination | Batch least squares OD | P18.01 |
| 20 | Professional Tools | Validate against GMAT and TLE/SGP4 | |

## Phase 20: Spacecraft Systems & the Space Environment ⬚ (~50 h)

| # | Lesson | Build it | Needs |
|---|---|---|---|
| 01 | Spacecraft Architecture | Subsystem block diagram | |
| 02 | The Space Environment | Environment model: vacuum, AO, plasma | P19.09 |
| 03 | Radiation Effects | Dose-depth and SEE estimates | P04.14 |
| 04 | Attitude Representation & Kinematics | Attitude kinematics library | P06.08 |
| 05 | Attitude Dynamics & Disturbance Torques | Gravity gradient, SRP, magnetic torques | P06.09 |
| 06 | Attitude Determination | TRIAD, QUEST, q-method | P18.10 |
| 07 | Attitude Actuators | Reaction wheels, magnetorquers, CMGs | |
| 08 | Attitude Control | Detumble and pointing controller | P09.15 |
| 09 | Momentum Management | Wheel desaturation | |
| 10 | Electrical Power: Solar & Batteries | Power budget simulator | P04.08 |
| 11 | Power Distribution & Regulation | EPS architecture trade | |
| 12 | Spacecraft Thermal Control | Lumped-node thermal model | P07.14 |
| 13 | Communications & Link Budgets | Link budget calculator | P04.11 |
| 14 | Modulation, Coding & Ground Segments | BER curves, coding gain | P09.06 |
| 15 | Command & Data Handling | Telemetry frame encoder | |
| 16 | Spacecraft Structures & Mechanisms | Launch-load structural check | P13.01 |
| 17 | Spacecraft Propulsion Integration | Δv and propellant budget | P15.14 |
| 18 | Payloads & Remote Sensing | Imaging resolution and swath | P04.12 |
| 19 | Reliability & Redundancy in Space | Fault tree for an EPS | |
| 20 | CubeSat Systems | CubeSat subsystem sizing | |
| 21 | Home Lab: Desktop Attitude Demonstrator | Guarded low-energy reaction-wheel platform with IMU; estimate and control attitude 🔬 | P09.21, P18.10 |

## Phase 21: Avionics, Embedded & Flight Software ⬚ (~40 h)

| # | Lesson | Build it | Needs |
|---|---|---|---|
| 01 | Digital Logic & Computer Architecture | Logic simulator | P04.08 |
| 02 | Microcontrollers & C for Embedded | Blink to sensor read on a microcontroller | |
| 03 | Communication Buses | I2C, SPI, UART, CAN drivers | |
| 04 | Avionics Buses | MIL-STD-1553 and ARINC 429 decoders | |
| 05 | Real-Time Operating Systems | Rate-monotonic scheduler | |
| 06 | Flight Software Architecture | Component-based FSW (F´ style) | |
| 07 | Fault Management & Watchdogs | Safe-mode state machine | |
| 08 | Telemetry & Command | Telemetry frame format, ground station decoder | |
| 09 | Software-in-the-Loop Simulation | SIL with your 6-DOF aircraft or spacecraft attitude sim | P17.02 or P20.08 |
| 10 | Hardware-in-the-Loop Testing | HIL rig design | |
| 11 | Safety-Critical Software (DO-178C) | Requirements-to-test trace | |
| 12 | Formal Methods & Static Analysis ★ | Model checking a mode logic | |
| 13 | Autopilot Firmware: PX4/ArduPilot | Modify and fly in simulation | P17.12 |
| 14 | Radiation-Tolerant Computing | TMR and EDAC | P20.03 |
| 15 | Avionics Integration | Integrated avionics architecture | |
| 16 | Cybersecurity for Aerospace Systems | Threat model for a UAV link | |
| 17 | Home Lab: CubeSat Bench Prototype | Non-flight power, sensor, and telemetry model with fault injection 🔬 | P20.10, P20.15 |

## Phase 22: Experimental Methods & Testing ⬚ (~35 h)

| # | Lesson | Build it | Needs |
|---|---|---|---|
| 01 | Measurement Systems & Uncertainty | Uncertainty propagation tool | P02.17 |
| 02 | Sensors & Transducers | Strain gauge bridge | P04.10 |
| 03 | Data Acquisition | DAQ pipeline (Arduino or similar) | P09.06 |
| 04 | Design of Experiments | Factorial test plan | |
| 05 | Wind Tunnel Testing | Wall and blockage corrections | P10.01 |
| 06 | Flow Visualization & PIV | PIV cross-correlation from scratch | |
| 07 | Pressure & Force Measurement | Balance calibration | |
| 08 | Structural Testing | Strain survey on a beam | P05.10 |
| 09 | Vibration & Modal Testing | Modal parameters from tap test | P06.14 |
| 10 | Propulsion Testing | Thrust stand design; instrumented electric-fan thrust test 🔬 | P14.01 or P15.02 |
| 11 | Environmental Testing | Thermal vac and vibe test plans | |
| 12 | Flight Testing | Flight test card and data reduction | P16.14 |
| 13 | Model Validation | Sim vs flight data correlation | P03.18 |
| 14 | Lab Report | Full lab report on a home experiment | P00.10 |

## Phase 23: Systems Engineering & Professional Practice ⬚ (~30 h)

| # | Lesson | Build it | Needs |
|---|---|---|---|
| 01 | The Systems Engineering Life Cycle | Life-cycle plan (NASA SE Handbook) | |
| 02 | Requirements Engineering | Requirements set with verification methods | |
| 03 | Functional Analysis & Architecture | Functional flow block diagram | |
| 04 | Trade Studies | Weighted trade with sensitivity | P02.16 |
| 05 | Interfaces & Configuration Management | ICD template | |
| 06 | Model-Based Systems Engineering | SysML v2 model of a CubeSat | |
| 07 | Reliability Engineering | Reliability block diagrams | P02.15 |
| 08 | Safety Analysis | FMEA, FTA, hazard analysis | |
| 09 | Certification & Regulation | FAA Part 23/25 and Part 107 walkthrough | |
| 10 | Space Regulation & Debris Mitigation | Debris compliance check | |
| 11 | Human Factors | Cockpit workload assessment | |
| 12 | Cost Estimation | Parametric cost model | |
| 13 | Project Management | Schedule and risk register | |
| 14 | Ethics in Aerospace Engineering | Case studies: Challenger, Columbia, 737 MAX | |
| 15 | Teamwork & Design Reviews | Run a mock PDR | |
| 16 | Technical Communication | Briefing deck and design review package | |

## Phase 24: Aircraft Design ⬚ (~45 h)

| # | Lesson | Build it | Needs |
|---|---|---|---|
| 01 | The Design Process & Requirements | Requirements for a design brief | P23.02 |
| 02 | Initial Sizing | Weight fraction sizing | P16.06 |
| 03 | Constraint Analysis | T/W–W/S constraint diagram | P16.10 |
| 04 | Configuration Layout | OpenVSP layout | |
| 05 | Wing Design | Wing sizing and airfoil selection | P10.12 |
| 06 | Tail Sizing & Stability | Tail volume and static margin | P17.03 |
| 07 | Propulsion Selection & Integration | Engine selection and scaling | P14.06 |
| 08 | Landing Gear & Subsystems | Gear placement and tip-over check | |
| 09 | Weight & Balance | Component weight build-up and CG envelope | |
| 10 | Structural Concept | Wing box concept and sizing | P13.07 |
| 11 | Performance & Mission Verification | Mission sim check | P16.11 |
| 12 | Cost & Operations | Direct operating cost | P23.12 |
| 13 | Design Optimization | Aircraft-level MDO | P03.16 |
| 14 | Design Report | Conceptual design report | |

## Phase 25: Spacecraft & Mission Design ⬚ (~45 h)

| # | Lesson | Build it | Needs |
|---|---|---|---|
| 01 | Mission Objectives & Concept of Operations | ConOps document | P23.02 |
| 02 | Mission Architecture Trades | Architecture trade space | P23.04 |
| 03 | Orbit & Constellation Design | Coverage analysis tool | P19.10 |
| 04 | Δv Budget & Propulsion Sizing | Mission Δv budget | P15.14 |
| 05 | Mass & Power Budgets | Integrated budget spreadsheet-as-code | P20.10 |
| 06 | Payload Sizing | Instrument sizing | P20.18 |
| 07 | ADCS Sizing | Actuator and sensor selection | P20.08 |
| 08 | Thermal Design | Hot and cold case analysis | P20.12 |
| 09 | Communications Architecture | End-to-end link and ground network | P20.13 |
| 10 | Launch Vehicle Selection & Integration | Launch vehicle trade | |
| 11 | Operations & Ground Systems | Ops plan | |
| 12 | Reliability, Risk & Cost | Mission risk and cost estimate | P23.07 |
| 13 | Concurrent Engineering Sessions | Run a concurrent design study | |
| 14 | Mission Design Report | Mission concept report (PDR level) | |

## Phase 26: Capstone Projects ⬚ (~200+ h)

Each capstone mirrors a university senior design or lab course: written requirements, a verification plan, a design review, analysis, a build or high-fidelity simulation, test, uncertainty and failure analysis, and a final report with an oral or video defense. Every path requires **one design capstone (01–04, 09) and one build-and-test capstone (05–08)**.

Every build-and-test capstone ships with three equipment tiers (minimum-cost, recommended, enhanced), an itemized BOM with substitutions and a reuse map, full build and safety documentation, and a computational or open-data fallback. The fallback earns modeling and simulation competencies, not physical experimentation competencies; to claim a full design-plus-build completion the learner performs at least one safe, accessible physical build-and-measure project that meets the common experiment rubric. No capstone requires at-home combustion, propellant fabrication, pressurized firing equipment, or autonomous outdoor flight; those happen only at qualified supervised clubs or facilities, as optional extensions.

| # | Capstone | Type | Needs |
|---|---|---|---|
| 01 | Design a Regional Electric Aircraft | Design | P24 |
| 02 | Design a CubeSat Earth-Observation Mission | Design | P25 |
| 03 | Design a Small Launch Vehicle | Design | P15, P19.16 |
| 04 | Design a Lunar Lander and Its Descent Guidance | Design | P18.16, P19.11, P20.08; P25 for the full-mission version |
| 05 | Build and Fly a Fixed-Wing UAV with Your Own Autopilot | Build & test | P17, P21 (01–16) |
| 06 | Propulsion Test Campaign: Instrumented Electric Thrust Stand at Home; Motor Firing Only at a Supervised Club | Build & test | P15, P22.10 |
| 07 | Build a Low-Speed Wind Tunnel and Validate Your Aero Codes | Build & test | P10 (01–13, 23), P22.05 |
| 08 | Build a Desktop Attitude-Control Demonstrator and Validate Your ADCS Models | Build & test | P20 (04–10), P21.17 |
| 09 | Design a Composite Wing Box to Certification Loads | Design | P12, P13, P23 |
| 10 | Open Research Capstone (reproduce a published result, then extend it) | Research | Any |

---

## Learning paths

Every path shares the **Core** (Phases 00–09, plus the bridge modules when the placement test says so). Paths then add professional-area phases and end in capstones. A phase listed without a range means all of its lessons; `20 (03)` means only lesson 03 of Phase 20, included as a bridge. Each path contains the transitive closure of every mandatory prerequisite of every lesson in it; `scripts/check_paths.py` verifies this from the tables in this file and must pass before any change to the roadmap is merged.

| Path | Phases after Core | Capstones | Scope |
|---|---|---|---|
| **Aerospace Engineer (full)** | 10–25 | One of 01–04 or 09, plus one of 05–08 | Full ABET aeronautical and astronautical scope |
| **Aeronautical Engineer** | 10, 11, 12, 13, 14, 16, 17, 18 (01–13), 20 (03), 21 (01–16), 22, 23, 24 | 01, plus 05 or 07 | ABET aeronautical scope |
| **Astronautical Engineer** | 10 (01), 11 (01–07, 12–13), 12, 13 (01–16, 19–21), 15, 18, 19, 20, 21 (01–12, 14–17), 22 (01–11, 13–14), 23, 25 | 02, 03, or 04, plus 06 or 08 | ABET astronautical scope |
| **Propulsion Specialist** | 10 (01–13, 17), 11, 12, 14, 15, 19 (01–04, 16), 22 (01–11, 13–14), 23 | 03, plus 06 | Specialization, not a degree scope |
| **GNC & Autonomy Specialist** | 10 (01–05, 09–13, 23), 16 (01–05), 17, 18, 19 (01–11), 20 (03–10, 15), 21 | 04, plus 05 or 08 | Specialization, not a degree scope |
| **Structures & Materials Specialist** | 10 (01–13, 19, 23), 12, 13, 22 (01–09, 11, 13–14), 23 | 09, plus 07 | Specialization, not a degree scope |

"Scope" means the subject matter matches the benchmark. It is claim 1 of the three in [Benchmark](#benchmark-what-top-program-graduate-means), never claim 2 or 3. The specialist paths are deliberately narrower than a degree, the way a master's student who already has the core would study.

---

## Home laboratory

Equipment is bought at the point of use, never as an upfront kit, and every 🔬 lesson offers a minimum-cost, a recommended, and an enhanced method, with household and phone-sensor substitutions where they are scientifically sound. No lesson requires a 3D printer, a machine shop, proprietary software, or a university lab. Every lab guide opens with a safety box (hazards, controls, supervision, operating envelope, stop criteria), and every instrument gets a calibration plan: a cheap instrument with known uncertainty teaches more than an uncalibrated expensive one.

| Level | Typical equipment | Introduced | Used for |
|---|---|---|---|
| **1. Basic engineering lab** | Ruler, tape, digital scale, weights, hand tools, safety glasses, clamps, breadboard, multimeter, microcontroller and a few sensors, phone camera | Phases 00–04 | Mechanics, beam deflection, pendulum dynamics, motion measurement, circuits, calibration, logging, first feedback control |
| **2. Prototyping lab** (optional) | Soldering tools with ventilation, small servos and motors, IMU, load cell, simple fixtures; 3D printer *or* cardboard, foam, pre-cut parts | Phases 05–13 | Materials testing, instrumented structures, actuators, embedded systems, controls, bench aerodynamics |
| **3. Aerospace experimentation** (optional) | Guarded low-speed bench airflow rig, instrumented foam glider, flight-data logger, desktop attitude platform, CubeSat bench prototype | Phases 10–21 | Lift and drag, structural validation, flight performance, sensor fusion, closed-loop control, subsystem tests |

### Integrated home projects

Each project is a package: engineering question, theory, unit-tested model, tiered BOM, safety box, procedure, measurement protocol, analysis, deliverables, and rubric. Matching one data point does not pass; the learner must explain the discrepancies.

| Project | Lesson(s) | Build or investigate | Evidence required |
|---|---|---|---|
| 01 Cantilever beam | P05.18 | Load a small beam, measure deflection | Predicted vs measured, load curve, uncertainty, error explanation |
| 02 Dynamic system identification | P06.19 | Pendulum or spring-mass, logged motion | Equations of motion, fitted damping, sensitivity |
| 03 Heat-transfer investigation | P07.15 | Cooling of a safe small object | Lumped vs distributed assumptions, fitted parameters |
| 04 Simple fluid experiment | P08.19 | Low-energy flow and pressure rig | Similarity reasoning, calibrated measurement, model comparison |
| 05 Low-speed bench airflow rig | P10.23, P22.05 | Guarded tabletop tunnel or equivalent | Flow characterization, blockage, repeatability |
| 06 Airfoil geometry and test | P10.23 | Hand-built sections tested on the rig | Theory vs computation vs experiment; stall and viscous limits |
| 07 Instrumented glider | P17.19 | Design, build, fly a foam glider with a logger | CG, stability, glide performance, structures, flight data |
| 08 Materials and structure comparison | P13.21 | Coupons, joints, built-up beams at low load | Stiffness, failure mode, variability, margins |
| 09 Feedback-control demonstrator | P09.21 | Servo or small motor under closed-loop control | Controller model, measured transient, robustness |
| 10 Satellite attitude demonstrator | P20.21 | Desktop platform with IMU and guarded reaction wheel | Attitude estimate, torque budget, tracking error |
| 11 CubeSat bench prototype | P21.17 | Non-flight power, telemetry, sensors, software | Power budget, fault cases, interface verification |
| 12 Autonomous flight in simulation | P18.20, P21.09 | Closed-loop UAV simulation; optional safe hardware extension | Stability, navigation error, mode transitions, failure injection |
| 13 Propulsion modeling and safe analogs | P22.10, P15.15 | Cycle and nozzle models, electric-fan thrust stand, published motor data | Thrust and efficiency predictions, validation |
| 14 Mission systems trade study | P25.02 | Constrained spacecraft mission with existing tools and data | Requirements traceability, budgets, risks |

### Flagship tracks

Two integrated systems that the learner revisits as each new discipline arrives, so the course is a spiral rather than a line. After the Core, a learner picks one, does both, or substitutes an equivalent supervised or club project; nobody is made to buy two hardware kits.

- **Track A, progressive aircraft:** a foam glider that grows through geometry and loads (P05), lift and drag estimates checked on the bench rig (P10), spar sizing checked against measured deflection (P13), glide, trim, and CG sensitivity (P16–17), a lightweight logger and a simulation-tested controller (P18, P21), controlled low-risk flight tests, and a final defense. Sophistication comes from the engineering process, not the cost of the airframe.
- **Track B, progressive spacecraft:** a notional mission whose orbit and environment are computed early (P19–20), subsystem budgets (P20), a desktop sensor, power, and telemetry prototype with attitude-control demonstration where practical (P20.21, P21.17), integration and fault handling (P21), subsystem tests, and a mission-design package plus prototype verification report (P25).

---

## Assessment

Four competencies are assessed separately for every major technical lesson: **mathematical and physical reasoning**, **analytical problem solving**, **computational competence**, and **engineering interpretation** (what the result means, when the model is valid, what could fail, how to validate it). A learner who produces plausible numbers without understanding the physics does not pass.

| Layer | What it is |
|---|---|
| Quick checks | Auto-graded questions on terminology, units, limiting cases, and intuition, each with feedback |
| Problem sets | Original derivations and unfamiliar engineering scenarios, with hint → method → full worked solution scaffolding |
| Implementation tests | Public unit tests plus unseen validation cases; residuals, conservation laws, and convergence rather than exact output matching |
| Experiment evaluations | Calibration, raw data, uncertainty, reproducibility, quality of conclusions |
| Phase exam | Integrated analytic, conceptual, computational, and interpretation questions, with a timed self-assessment mode |
| Phase mini-project | Rubric-graded investigation integrating several lessons; the phase is passed on exam *and* project |
| Design reviews | Written package and oral or video defense against a PDR/CDR-style checklist (capstones) |
| External challenge | Optional peers, mentors, clubs, or reviewers assess the reasoning |

Project rubric (proposed weights, to be calibrated in pilots): correct physics and assumptions 25%, calculation and numerical correctness 20%, experiment design and data quality 20%, verification and validation 20%, communication, reproducibility, and safety documentation 15%. Minimum thresholds apply in physics and safety; presentation quality cannot compensate for either. Every major lesson ends with adversarial self-checks: what if the load doubles, which assumption fails first, what measurement would falsify this, does the result survive a change of mesh or time step.

Achievement is reported in three distinct levels: **practice completed**, **self-assessed mastery**, and **externally reviewed mastery**. None of them is a credential.

---

## Outcome coverage

How the curriculum hits each ABET student outcome. Self-learners can't get a team or a grader by default, so the outcomes that need other people get explicit mechanisms.

| ABET outcome | Where it's built | Mechanism |
|---|---|---|
| 1. Solve complex problems | Every phase | Problem sets and phase exams with worked solutions |
| 2. Design to meet needs | 24, 25, 26; every mini-project from Phase 00 | Design briefs with requirements, margins, risks, and a verification plan |
| 3. Communicate | 00.10, 22.14, 23.16, every capstone | Report and briefing rubrics; public write-ups |
| 4. Ethics & professional responsibility | 23.09, 23.10, 23.14 | Case studies; ethics section required in each design report |
| 5. Teamwork | 23.15, capstones | Optional cohort and peer-review program via GitHub Discussions |
| 6. Experimentation & data analysis | 🔬 lessons from 00.14 onward, 22 as the formal methods phase | Home experiments with calibration and uncertainty, public datasets (NASA, flight logs) |
| 7. Learn independently | 00.11, open research capstone | Paper-reading lessons; ★ lessons point to the literature |

---

## Decisions for step 2

Settled on 2026-10-11 (full text in [docs/RECOMMENDATIONS-2026-10-11.md](docs/RECOMMENDATIONS-2026-10-11.md), section 20).

- **Hardware.** Every hands-on lesson and capstone publishes three budget-tier BOMs with substitutions and a reuse map, full build and safety documentation, and a simulation or open-data fallback, and states explicitly which outcomes were simulated and which were physically measured. Prices are verified at authoring time and dated, never treated as permanent. Motor firing happens only in qualified supervised settings; accessible physical projects (bench rig, instrumented structure, attitude demonstrator) keep the build-and-test route open to everyone.
- **Exams.** Auto-graded checks, original rigorous problem sets, phase exams, and mini-projects, with hint → method → full-solution scaffolding. Full solutions live in the public repo in a separate `solutions/` tree, hidden by default in the website but freely accessible, with no payment, registration, or time locks. No private instructor folder: in a public repo that is not access control. Public unsupervised exams never claim credential equivalence.
- **Textbook anchors.** Each phase gets a verified map from lessons to exact sections of free-to-access references (seeded with MIT OpenCourseWare, OpenStax, Leishman's *Introduction to Aerospace Flight Vehicles*, NASA NTRS, and the NASA Systems Engineering Handbook), with paid textbooks listed only as clearly optional enrichment and never as prerequisites. A machine-readable manifest records title, creator, version, URL, license, access tier, mapped lessons, and verification date, with link checking in CI and periodic human review. Coverage gaps are filled with original instruction, never with a paid book.

**Repository layout (decided):** `curriculum/phase-XX/lesson-YY/{lesson.md, problems.md, starter/, tests/}`, `solutions/phase-XX/lesson-YY/{hints.md, worked-solutions.md}`, `assessments/phase-XX/{exam.md, grading-rubric.md, parameterized-cases/}`, `projects/<id>/{bom.yaml, safety.md, procedure.md, data/, report-template.md}`, `references/{sources.yaml, phase-XX.md}`, `scripts/check_paths.py`.

**Definition of done for step 2.** Each phase has a verified reference map and an original assessment plan; every physical project has a safety-reviewed, cost-tiered build package and an honest accessible alternative; the lesson template encodes the six-part arc and the metadata block (outcomes, prerequisites, engineering question, assumptions, analytical work, computational deliverable, practical activity, equipment, safety, evidence, assessment, references, extensions); and another author can write a lesson that meets the same standard.

**Still open, deferred to pilots.** Time estimates per activity type; the diagnostic entry test; which three lessons (beam, pendulum, first data collection) to pilot with real beginners before anything scales.
