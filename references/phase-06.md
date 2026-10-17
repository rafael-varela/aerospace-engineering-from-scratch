# Reference map: Phase 06: Dynamics & Vibrations

**Note 2026-10-17:** first mapping of Phase 06. **No live page could be loaded on 2026-10-17.** WebFetch failed with "getaddrinfo ENOTFOUND" for every host tried: ocw.mit.edu, openstax.org, ntrs.nasa.gov, phyphox.org and google.com. A direct HTTPS request through the session's egress proxy was refused ("CONNECT tunnel failed, response 403", connect_rejected) for ocw.mit.edu. Status 000 came back for openstax.org, ntrs.nasa.gov, standards.nasa.gov, phys.libretexts.org, phyphox.org, physlets.org, web.archive.org and dspace.mit.edu. Only github.com answered, and repository access there was not enabled for this session. So this map adds **no new source ids**. Every section title below was copied from this repository's earlier reference maps (phase-01 to phase-05), which recorded those titles from pages loaded between 2026-10-11 and 2026-10-17. The existing manifest entries were extended to name the Phase 06 lessons; their `verified` dates were not changed. **Rejected, for now, because nothing could be loaded:**
- MIT OCW 16.07 *Dynamics* (Fall 2009); 2.003SC *Engineering Dynamics* (Fall 2011); 2.003J/1.053J *Dynamics and Vibration* (Spring 2007); 16.61 *Aerospace Dynamics*.
- OpenStax *University Physics* Vol. 1 sections not recorded in an earlier map (for example 11.1 and the Coriolis subsection of 6.3).
- LibreTexts classical-mechanics texts (Cline; Fowler; Tatum).
- NASA-HDBK-7005, NASA-STD-7001B and NASA-HDBK-7004.
- NTRS NASA-TM-74839 and NACA landing-gear reports.
- phyphox and Tracker.

A web *search* (not a page load) was used only to get the names and URLs of these items right. Each one appears below as a "candidate (unverified)" on a GAP line, with no `@` id. A candidate becomes citable only after a page load verifies its contents and licence and it is added to `sources.yaml`. Lesson authors must not cite candidates.

Section numbers and titles are copied exactly as the earlier maps recorded them from the publishers' own tables of contents, section pages or PDF outlines. Every source here is **required-free** unless marked otherwise; sources are identified by an `@`-prefixed identifier and recorded in [sources.yaml](sources.yaml) with license, version, and caveats. **GAP** marks a lesson that no listed free source covers adequately: the lesson must carry original instruction there, and any "candidate" named is unverified until it is added to the manifest.

Conventions: OpenStax section URLs follow `https://openstax.org/books/<book-slug>/pages/<n>-<m>-<title-slug>`; a URL written as `.../pages/...` continues the last full URL above it. MIT OCW lecture URLs follow `https://ocw.mit.edu/courses/<course-slug>/resources/lecture-<n>-<slug>/`. "(loaded YYYY-MM-DD)" means an earlier phase's map records that exact URL as loaded on that date. "(pattern)" means the URL was built from a verified slug scheme and has not been loaded. Nothing in this map was loaded on 2026-10-17. Licences:
- **CC BY-NC-SA 4.0:** OpenStax, MIT OCW and LibreTexts.
- **CC BY 4.0:** Langtangen & Linge.
- **BSD 3-Clause:** SciPy and SymPy.
- Free to read, not to redistribute: JCGM 100:2008.
- US Government work: NASA *Basics of Space Flight*.

The curriculum links to these sources and writes its own explanations; nothing is copied or adapted.


## 06-01 Particle Kinematics in Many Coordinates

Sources: `@openstax-university-physics-vol1`, `@openstax-calculus-volume-3`, `@openstax-calculus-volume-2`

- OpenStax University Physics Vol. 1 Ch. 4 "Motion in Two and Three Dimensions": 4.3 Projectile Motion; 4.4 Uniform and Nonuniform Circular Motion
  - https://openstax.org/books/university-physics-volume-1/pages/4-4-uniform-and-nonuniform-circular-motion (pattern)
- OpenStax Calculus Vol. 3: 3.1 Vector-Valued Functions and Space Curves; 3.2 Calculus of Vector-Valued Functions; 3.4 Motion in Space
  - https://openstax.org/books/calculus-volume-3/pages/3-4-motion-in-space (pattern)
- OpenStax Calculus Vol. 2: 7.1 Parametric Equations; 7.3 Polar Coordinates
  - https://openstax.org/books/calculus-volume-2/pages/7-3-polar-coordinates (pattern)
- Extract: velocity and acceleration as derivatives of a vector-valued position; tangential and normal components of acceleration (3.4) as the path-coordinate description; centripetal and tangential acceleration in circular motion (4.4); polar coordinates as the geometric starting point.
- Partial GAP: none of the recorded sections derive the polar or cylindrical unit-vector derivatives (ė_r = θ̇ e_θ, ė_θ = −θ̇ e_r) or the radial and transverse acceleration components (r̈ − rθ̇², rθ̈ + 2ṙθ̇). That derivation, and the conversions between the three descriptions, are original. Candidate (unverified): MIT OCW 16.07 *Dynamics* (Fall 2009) lecture notes, https://ocw.mit.edu/courses/16-07-dynamics-fall-2009/pages/lecture-notes/

## 06-02 Particle Kinetics

Sources: `@openstax-university-physics-vol1`, `@scipy-userguide`

- OpenStax Vol. 1 Ch. 5 "Newton's Laws of Motion": 5.3 Newton's Second Law; 5.7 Drawing Free-Body Diagrams
- Ch. 6 "Applications of Newton's Laws": 6.1 Solving Problems with Newton's Laws; 6.3 Centripetal Force; 6.4 Drag Force and Terminal Speed
  - https://openstax.org/books/university-physics-volume-1/pages/6-3-centripetal-force (pattern); .../pages/6-4-drag-force-and-terminal-speed (pattern)
- SciPy User Guide: "Integration (scipy.integrate)" > "Ordinary differential equations (solve_ivp)"
  - https://docs.scipy.org/doc/scipy/tutorial/index.html (pattern)
- Extract: ΣF = m a resolved in path (tangential–normal) coordinates; centripetal force for a curved path; integrating a point-mass flight path with `solve_ivp`.
- Partial GAP: no recorded free section treats load factor (n = L/W, n = 1 + V²/(gR) at the bottom of a pull-up), seat and wing loads, or g-onset. The aircraft application is original. Candidate (unverified): MIT OCW 16.07 lecture notes (particle kinetics in normal–tangential coordinates).

## 06-03 Work-Energy & Impulse-Momentum

Sources: `@openstax-university-physics-vol1`, `@scipy-userguide`

- OpenStax Vol. 1 Ch. 7 "Work and Kinetic Energy": 7.1 Work; 7.2 Kinetic Energy; 7.3 Work-Energy Theorem; 7.4 Power
- Ch. 8 "Potential Energy and Conservation of Energy": 8.1 Potential Energy of a System; 8.2 Conservative and Non-Conservative Forces; 8.3 Conservation of Energy
  - https://openstax.org/books/university-physics-volume-1/pages/8-3-conservation-of-energy (pattern)
- Ch. 9 "Linear Momentum and Collisions": 9.1 Linear Momentum; 9.2 Impulse and Collisions; 9.3 Conservation of Linear Momentum
  - .../pages/9-2-impulse-and-collisions (pattern)
- SciPy "Ordinary differential equations (solve_ivp)" (the independent time-integration check)
- Extract: W = ΔT along the stopping stroke; the impulse ∫F dt = Δp; why a force law that depends on velocity (a damper) needs both views, plus time integration, to give the peak deceleration.
- Partial GAP: no verified free source describes carrier arresting gear (the cable and its energy absorber) or gives realistic arrestment numbers. The force law and its numbers are original and labelled *illustrative*.

## 06-04 Relative Motion & Rotating Frames

Sources: `@openstax-university-physics-vol1`, `@mit-ocw-18-06-linear-algebra`

- OpenStax Vol. 1 Ch. 4: 4.5 Relative Motion in One and Two Dimensions
  - https://openstax.org/books/university-physics-volume-1/pages/4-5-relative-motion-in-one-and-two-dimensions (pattern)
- OpenStax Vol. 1 Ch. 6: 6.3 Centripetal Force (the earlier map recorded only the section title; whether the section discusses non-inertial frames was not checked)
- MIT 18.06 Lecture 31: Change of basis; image compression (coordinate changes between frames)
  - https://ocw.mit.edu/courses/18-06-linear-algebra-spring-2010/resources/lecture-31-change-of-basis-image-compression/ (pattern)
- Extract: relative velocity in translating frames (4.5); a change of basis as the algebra behind expressing one vector in two frames.
- GAP: no verified free source derives the rotating-frame transport theorem (d/dt)_I = (d/dt)_B + ω×, the five-term acceleration (including Coriolis 2ω×v_rel and centripetal ω×(ω×r)), Coriolis deflection on a rotating Earth, or the Foucault pendulum. All of this is original. Candidates (unverified): MIT OCW 16.07 lecture notes on relative motion and rotating frames; MIT OCW 2.003SC lecture "Newton's Laws, Vectors and Reference Frames" (resource MIT2_003SCF11Kinematic.pdf seen in search results).

## 06-05 Rigid Body Kinematics in 2D

Sources: `@openstax-university-physics-vol1`, `@scipy-userguide`, `@ocw-16-90`

- OpenStax Vol. 1 Ch. 10 "Fixed-Axis Rotation": 10.1 Rotational Variables; 10.3 Relating Angular and Translational Quantities
  - https://openstax.org/books/university-physics-volume-1/pages/10-3-relating-angular-and-translational-quantities (pattern)
- SciPy "Optimization (scipy.optimize)" > "Root finding" > "Sets of equations" (loop-closure equations solved numerically)
- MIT 16.90 Unit 1: "1.7.4 Newton-Raphson Implement Implicit Methods on Nonlinear Problems" (titles as published; Newton's method for a nonlinear system)
  - https://ocw.mit.edu/courses/16-90-computational-methods-in-aerospace-engineering-spring-2014/ (pattern)
- Extract: v = ω × r for a point on a rotating body; Newton's method on the two loop-closure equations of a four-bar linkage.
- GAP: no verified free source covers planar rigid-body kinematics for linkages: the relative-velocity and relative-acceleration equations v_B = v_A + ω × r_B/A, instant centres, Grashof's criterion, and four-bar position, velocity and acceleration analysis. All of this is original. Candidates (unverified): MIT OCW 16.07 lecture notes (2-D rigid-body kinematics); MIT OCW 2.003SC lectures.

## 06-06 Rigid Body Kinetics in 2D

Sources: `@openstax-university-physics-vol1`, `@libretexts-baker-haynes-engineering-statics`

- OpenStax Vol. 1 Ch. 6: 6.2 Friction
- Ch. 10: 10.4 Moment of Inertia and Rotational Kinetic Energy; 10.6 Torque; 10.7 Newton's Second Law for Rotation
  - https://openstax.org/books/university-physics-volume-1/pages/10-7-newtons-second-law-for-rotation (pattern)
- Baker & Haynes Ch. 9 "Friction" (engineeringstatics.org outline): 9.1 Dry Friction (9.1.1 Coulomb Friction)
  - https://engineeringstatics.org/Chapter_09-dry-friction.html (as linked from the outline loaded 2026-10-17)
- Baker & Haynes Ch. 10 "Moments of Inertia": 10.8 Mass Moment of Inertia
  - https://eng.libretexts.org/Bookshelves/Mechanical_Engineering/Engineering_Statics%3A_Open_and_Interactive_(Baker_and_Haynes)/10%3A_Moments_of_Inertia (loaded 2026-10-17)
- Extract: ΣF = m a_G together with ΣM_G = I_G α; kinetic Coulomb friction μ_k N opposing the slip velocity; the slip condition v_G ≠ ωr, which ends when rolling begins.
- Partial GAP: rolling motion (OpenStax 11.1) was not recorded in an earlier map, so it is not cited. No verified free source treats wheel spin-up at touchdown or the friction energy lost while the tyre slides, so the sliding-to-rolling transition and its two-way energy check are original.

## 06-07 3D Rotations: Euler Angles & DCMs

Sources: `@mit-ocw-18-06-linear-algebra`, `@nasa-basics-of-space-flight`

- MIT 18.06 Lecture 30: Linear transformations and their matrices; Lecture 31: Change of basis; image compression; Lecture 17: Orthogonal matrices and Gram-Schmidt; Lecture 18: Properties of determinants
  - https://ocw.mit.edu/courses/18-06-linear-algebra-spring-2010/resources/lecture-17-orthogonal-matrices-and-gram-schmidt/ (pattern); .../resources/lecture-30-linear-transformations-and-their-matrices/ (pattern)
- Context only: NASA *Basics of Space Flight* Chapter 11 "Onboard Systems" page 2, "Attitude and Articulation Control Subsystems (AACS)", qualitative
  - https://science.nasa.gov/learn/basics-of-space-flight/chapter11-2/ (loaded 2026-10-16)
- Extract: QᵀQ = I and det Q = +1 define a rotation; Gram–Schmidt as the re-orthonormalisation step for a DCM that has drifted during integration.
- GAP (as already noted in the 01-15 map): Euler-angle sequences (3-2-1 and others), DCM construction from elementary rotations, Euler-angle rates from body rates, gimbal lock at θ = ±90°, and the DCM kinematic equation Ċ = −[ω×]C are all original. Candidate (unverified): NTRS 19770024290, NASA-TM-74839, D. M. Henderson (McDonnell Douglas Technical Services), 1977, on the relationships between Euler angles, quaternions and transformation matrices, https://ntrs.nasa.gov/citations/19770024290. A companion contractor report, NTRS 19770019231, gives the matching utility subroutines.

## 06-08 Quaternions

Sources: `@mit-ocw-18-06-linear-algebra`, `@nasa-basics-of-space-flight`

- MIT 18.06 Lecture 17: Orthogonal matrices and Gram-Schmidt (the DCM side of every quaternion ↔ DCM check)
- Context only: Basics of Space Flight Chapter 11 page 2, AACS (run-in labels "Celestial Reference" and "Inertial Reference": star trackers and inertial references, qualitative)
  - https://science.nasa.gov/learn/basics-of-space-flight/chapter11-2/ (loaded 2026-10-16)
- GAP: no verified free source covers quaternion algebra (product, conjugate, rotating a vector by q v q*), quaternion ↔ DCM ↔ Euler conversions, the quaternion kinematic equation q̇ = ½ q ⊗ (0, ω), norm drift and renormalisation, or SLERP. The whole lesson is original. Candidate (unverified): NASA-TM-74839 (see 06-07).

## 06-09 3D Rigid Body Dynamics: Euler's Equations

Sources: `@openstax-university-physics-vol1`, `@libretexts-baker-haynes-engineering-statics`, `@mit-ocw-18-06-linear-algebra`, `@ocw-16-90`, `@scipy-userguide`

- OpenStax Vol. 1 Ch. 11 "Angular Momentum": 11.2 Angular Momentum; 11.3 Conservation of Angular Momentum
  - https://openstax.org/books/university-physics-volume-1/pages/11-3-conservation-of-angular-momentum (pattern)
- Baker & Haynes Ch. 10: 10.7 Products of Inertia; 10.8 Mass Moment of Inertia
- MIT 18.06 Lecture 25: Symmetric matrices and positive definiteness (the inertia tensor's principal axes are the eigenvectors of a symmetric matrix)
  - https://ocw.mit.edu/courses/18-06-linear-algebra-spring-2010/resources/lecture-25-symmetric-matrices-and-positive-definiteness/ (pattern)
- MIT 16.90 Unit 1: "1.9 Runge-Kutta Methods"; SciPy "Ordinary differential equations (solve_ivp)"
- Extract: L = I ω with a full inertia tensor; principal axes from the eigenproblem; integrating a nonlinear ODE system and monitoring conserved quantities.
- Partial GAP: no verified free source derives Euler's equations I₁ω̇₁ = (I₂ − I₃)ω₂ω₃ (and the cyclic pair), the instability of rotation about the intermediate axis (the Dzhanibekov effect), or the energy and momentum ellipsoids. That theory is original. Candidate (unverified): MIT OCW 16.07 Lecture 28, "3D Rigid Body Dynamics: Equations of Motion; Euler's Equations" (resource mit16_07f09_lec28).

## 06-10 Gyroscopes & Spin Stability

Sources: `@openstax-university-physics-vol1`, `@nasa-basics-of-space-flight`

- OpenStax Vol. 1 Ch. 11: 11.4 Precession of a Gyroscope
  - https://openstax.org/books/university-physics-volume-1/pages/11-4-precession-of-a-gyroscope (pattern)
- Basics of Space Flight Chapter 11 page 2, AACS: run-in labels "Spin" and "3-Axis" (spin-stabilised compared with three-axis-stabilised spacecraft, qualitative)
  - https://science.nasa.gov/learn/basics-of-space-flight/chapter11-2/ (loaded 2026-10-16)
- Extract: gyroscopic precession rate Ω = τ/L; why a spinning spacecraft holds its attitude.
- Partial GAP: no verified free source covers the major/minor-axis rule, the nutation frequency of an axisymmetric spinner, or the energy-sink argument that turns a minor-axis spin into a flat spin. The Explorer 1 history (1958, whip antennas dissipating energy) also needs a verified primary source. All of this is original, and the history is stated without a citation until one is verified. Candidates (unverified): MIT OCW 16.07 Lecture 29 (Euler angles and the free motion of a rotating body) and Lecture 30 (tops and gyroscopes) (resources mit16_07f09_lec29, mit16_07f09_lec30); MIT OCW 16.61 *Aerospace Dynamics*.

## 06-11 Lagrangian Mechanics

Sources: `@mit-ocw-8-033-relativity` (supplementary-free), `@mit-ocw-16-323-principles-of-optimal-control`, `@sympy-tutorial`, `@mit-ocw-18-03-differential-equations`, `@scipy-userguide`

- MIT 8.033 Session 8 "Variational Calculus" slides (Euler–Lagrange equations and constants of motion; brachistochrone)
  - https://ocw.mit.edu/courses/8-033-relativity-fall-2006/resources/lecture8_euler/ (pattern)
- MIT 16.323 Lecture 5: Calculus of variations (functionals, first variation, Euler–Lagrange equation); Lecture 2: Nonlinear optimization: constrained nonlinear optimization, Lagrange multipliers
  - https://ocw.mit.edu/courses/16-323-principles-of-optimal-control-spring-2008/resources/lec5/ (pattern)
- SymPy 1.14.0 "Introductory Tutorial": "SymPy Features" > "Calculus", "Solvers", "Matrices" (deriving and solving the equations of motion symbolically, then lambdify)
  - https://docs.sympy.org/latest/tutorials/intro-tutorial/index.html (pattern)
- MIT 18.03 notes 38 "Linearization near equilibria; the nonlinear pendulum" (small-angle linearisation before the normal-mode step)
  - https://ocw.mit.edu/courses/18-03-differential-equations-spring-2010/pages/lecture-notes/ (pattern)
- Extract: the Euler–Lagrange equation from P02.13, applied to L = T − V in generalised coordinates; symbolic derivation; linearising about equilibrium.
- Partial GAP (as noted in the 02-13 map): neither source applies Hamilton's principle to mechanical systems. Generalised coordinates, L = T − V for a double pendulum, generalised forces, and the check against Newton–Euler are original. Candidates (unverified): MIT OCW 2.003SC Lecture 15 "Introduction to Lagrange With Examples", Lecture 16 "Kinematic Approach to Finding Generalized Forces" and Lecture 17 "Practice Finding EOM Using Lagrange Equations" (titles seen in search results for the course's resource index, https://ocw.mit.edu/courses/2-003sc-engineering-dynamics-fall-2011/pages/resource-index/); MIT OCW 16.07 "Energy methods: Lagrange's equations"; D. Cline, *Variational Principles in Classical Mechanics* (LibreTexts; the licence was reported as CC BY-NC-SA 4.0 but not checked on the page).

## 06-12 Hamiltonian Mechanics ★

Sources: `@langtangen-linge-fd-computing-pdes`, `@mit-ocw-18-03-differential-equations`, `@openstax-university-physics-vol1`, `@ocw-18-085`

- Langtangen & Linge, *Finite Difference Computing with PDEs*, Ch. 1 "Vibration ODEs": 1.6 Energy Considerations (1.6.1 Derivation of the Energy Expression; 1.6.2 An Error Measure Based on Energy); 1.7 The Euler-Cromer Method (1.7.4 The Störmer-Verlet Algorithm); 1.5 Alternative Schemes Based on 1st-Order Equations (1.5.1 The Forward Euler Scheme)
  - https://link.springer.com/chapter/10.1007/978-3-319-55456-3_1 (pattern)
- MIT 18.03 notes 35 "Qualitative behavior of linear systems; phase plane"; notes 38 "Linearization near equilibria; the nonlinear pendulum"
- OpenStax Vol. 1 Ch. 13 "Gravitation": 13.3 Gravitational Potential Energy and Total Energy; 13.5 Kepler's Laws of Planetary Motion
  - https://openstax.org/books/university-physics-volume-1/pages/13-5-keplers-laws-of-planetary-motion (pattern)
- MIT 18.085 "Lecture 10: Finite Differences in Time"
  - https://ocw.mit.edu/courses/18-085-computational-science-and-engineering-i-fall-2008/resources/lecture-10-finite-differences-in-time/ (pattern)
- Extract: pendulum phase portraits (separatrix, libration and rotation); energy as the conserved quantity that labels each orbit; why Störmer–Verlet keeps the energy error bounded while forward Euler's grows.
- Partial GAP (the same gap as 03-07): no verified free source covers the Legendre transform, Hamilton's equations, canonical momenta, Liouville's theorem (phase-space area preservation) or symplectic maps. These are original.

## 06-13 Single-DOF Vibration

Sources: `@mit-ocw-18-03-differential-equations`, `@openstax-university-physics-vol1`, `@langtangen-linge-fd-computing-pdes`, `@scipy-userguide`

- MIT 18.03 L9 Solving Second-order Linear ODE's with Constant Coefficients; L10 Continuation: Complex Characteristic Roots; L13 Finding Particular Solutions to Inhomogeneous ODEs; L14 Interpretation of the Exceptional Case: Resonance
  - https://ocw.mit.edu/courses/18-03-differential-equations-spring-2010/resources/lecture-14-interpretation-of-the-exceptional-case-resonance/ (pattern)
- MIT 18.03 notes 11 "Modes and the characteristic polynomial"; 12 "Good vibrations, damping conditions"; 13 "Exponential response formula, spring drive"; 14 "Complex gain, dashpot drive"; 16 "Frequency response"
  - https://ocw.mit.edu/courses/18-03-differential-equations-spring-2010/pages/lecture-notes/ (pattern)
- OpenStax Vol. 1 Ch. 15 "Oscillations": 15.1 Simple Harmonic Motion; 15.2 Energy in Simple Harmonic Motion; 15.5 Damped Oscillations; 15.6 Forced Oscillations
  - https://openstax.org/books/university-physics-volume-1/pages/15-6-forced-oscillations (pattern)
- Langtangen & Linge Ch. 1: 1.4 Analysis of the Numerical Scheme (1.4.2 The Error in the Numerical Frequency)
- Extract: ω_n, ζ and the under-, critically and over-damped cases; the amplitude and phase of the response to harmonic forcing (complex gain); resonance; checking the closed form against `solve_ivp`.
- Partial GAP: the 18.03 notes 14 title shows a dashpot-driven input, but no verified section treats base excitation and the transmissibility ratio |X/Y| explicitly, or the logarithmic decrement δ = ln(x_k/x_{k+1}). These are original. Candidates (unverified): MIT OCW 2.003SC Lecture 19 "Introduction to Mechanical Vibration" and Lecture 20 "Linear System Modeling a Single Degree of Freedom Oscillator"; MIT OCW 2.003J (Spring 2007).

## 06-14 Multi-DOF Vibration & Modes

Sources: `@mit-ocw-18-03-differential-equations`, `@mit-ocw-18-06-linear-algebra`, `@scipy-userguide`

- MIT 18.03 L29 Matrix Exponentials; L30 Decoupling Linear Systems with Constant Coefficients; notes 36 "Normal modes and the matrix exponential"
  - https://ocw.mit.edu/courses/18-03-differential-equations-spring-2010/resources/lecture-29-matrix-exponentials/ (pattern)
- MIT 18.06 Lecture 21: Eigenvalues and eigenvectors; Lecture 22: Diagonalization and powers of A; Lecture 25: Symmetric matrices and positive definiteness
  - https://ocw.mit.edu/courses/18-06-linear-algebra-spring-2010/resources/lecture-21-eigenvalues-and-eigenvectors/ (pattern)
- SciPy "Linear Algebra (scipy.linalg)" (the generalised symmetric eigenproblem solver)
- Extract: normal modes as eigenvectors; decoupling by a change of coordinates (modal coordinates); real eigenvalues and orthogonal eigenvectors for symmetric matrices.
- Partial GAP: the recorded sources treat the standard eigenproblem A x = λ x. None covers the generalised problem K φ = ω² M φ, mass normalisation (φᵀ M φ = 1), M- and K-orthogonality, modal superposition of a forced response, or proportional (Rayleigh) damping. These are original. Candidates (unverified): MIT OCW 2.003SC Lecture 22 "Finding Natural Frequencies and Mode Shapes of a 2 DOF System", Lecture 23 "Vibration by Mode Superposition" and Lecture 24 "Modal Analysis: Orthogonality, Mass Stiffness, Damping Matrix".

## 06-15 Continuous Systems

Sources: `@ocw-18-085`, `@ocw-16-90`, `@mit-ocw-3-11-mechanics-of-materials`, `@mit-ocw-18-03-differential-equations`, `@scipy-userguide`

- MIT 18.085 "Lecture 17: Finite Elements in 1D (part 1)", "Lecture 18: Finite Elements in 1D (part 2)", "Lecture 19: Quadratic/Cubic Elements", "Lecture 20: Element Matrices; 4th Order Bending Equations"
  - https://ocw.mit.edu/courses/18-085-computational-science-and-engineering-i-fall-2008/ (pattern)
- MIT 16.90 Unit 2: "2.9 Introduction to Finite Elements"; "2.10 More on Finite Element Methods"
- Roylance, MIT 3.11 module "Beam Displacements" (PDF dated November 30, 2000): Multiple integration (the static EI v'''' = q equation behind the dynamic one)
  - https://ocw.mit.edu/courses/3-11-mechanics-of-materials-fall-1999/resources/mit3_11f99_bdisp/ (loaded 2026-10-17)
- MIT 18.03 notes 20 "Fourier series" (separation of variables and eigenfunction expansions; see also the 02-11 map)
- SciPy "Linear Algebra (scipy.linalg)"; "Sparse eigenvalue problems with ARPACK"
- Extract: the fourth-order bending operator and its element matrices; assembling stiffness (and mass) matrices; convergence of a discrete eigenproblem.
- Partial GAP: no verified free source derives the Euler–Bernoulli beam vibration equation EI w'''' + ρA ẅ = 0, the cantilever frequency equation cos βL cosh βL = −1, the pinned–pinned frequencies, the consistent mass matrix of a beam element, or the first bending modes of a wing. These are original.

## 06-16 Random Vibration

Sources: `@mit-ocw-6-011-signals-systems-and-inference`, `@mit-ocw-6-003-signals-and-systems`, `@mit-ocw-18-03-differential-equations`, `@scipy-userguide`

- MIT 6.011 Lecture 16: Wide-sense Stationary Processes, LTI Filtering of WSS Processes; Lecture 17: LTI Filtering of WSS Processes; Lecture 18: Power Spectral Density; Lecture 19: Einstein-Wiener-Khinchin Theorem, PSD Applications, Modeling Filters
  - https://ocw.mit.edu/courses/6-011-signals-systems-and-inference-spring-2018/pages/lecture-notes/ (pattern)
- MIT 6.003 Lecture 21: Sampling (sampling a synthesised random signal; aliasing)
  - https://ocw.mit.edu/courses/6-003-signals-and-systems-fall-2011/resources/lecture-21-sampling/ (pattern)
- MIT 18.03 notes 16 "Frequency response" (the |H(ω)| that shapes the output PSD)
- SciPy "Integration (scipy.integrate)" (integrating a PSD for the mean square)
- Extract: S_yy(ω) = |H(ω)|² S_xx(ω) for an LTI system driven by a WSS input; mean square equals the area under the PSD; Wiener–Khinchin.
- Partial GAP: no verified free source covers acceleration PSDs in g²/Hz, Grms, the piecewise log–log test spectra used in launch testing, Miles' equation Grms = √(π/2 · f_n · Q · W(f_n)), or the 3σ peak convention. These are original, and any spectrum is labelled a *NASA-STD-7001-style illustrative* input, not the standard's values. Candidates (unverified): NASA-STD-7001B *Payload Vibroacoustic Test Criteria*, https://standards.nasa.gov/standard/NASA/NASA-STD-7001; NASA-HDBK-7005 *Dynamic Environmental Criteria* (search results pointed only to third-party copies, so no official URL was confirmed).

## 06-17 Vibration Isolation

Sources: `@mit-ocw-18-03-differential-equations`, `@mit-ocw-6-011-signals-systems-and-inference`

- MIT 18.03 notes 14 "Complex gain, dashpot drive"; 16 "Frequency response"; 30 "The transfer function and frequency response"
  - https://ocw.mit.edu/courses/18-03-differential-equations-spring-2010/pages/lecture-notes/ (pattern)
- MIT 6.011 Lecture 17: LTI Filtering of WSS Processes (the random-input response that 06-16 builds on)
- Extract: an isolator as a low-pass filter from base motion to payload motion; gain above and below the corner frequency; how damping trades resonant peak against high-frequency attenuation.
- Partial GAP: no verified free source treats isolator design: transmissibility curves and the √2 crossover, choosing stiffness for a target corner frequency, sway space (rattle space), or IMU isolation. These are original. Candidate (unverified): NASA-HDBK-7005 (see 06-16).

## 06-18 Multibody Simulation

Sources: `@ocw-16-90`, `@scipy-userguide`, `@openstax-university-physics-vol1`

- MIT 16.90 Unit 1: "1.6 Systems of ODE's and Eigenvalue Stability"; "1.7 Stiffness and Implicit Methods" (1.7.1 Stiffness)
  - https://ocw.mit.edu/courses/16-90-computational-methods-in-aerospace-engineering-spring-2014/ (pattern)
- SciPy "Ordinary differential equations (solve_ivp)" (RK45 compared with Radau and BDF for a stiff strut model; events for touchdown and liftoff)
- OpenStax Vol. 1 Ch. 8: 8.2 Conservative and Non-Conservative Forces (energy accounting with a dissipative damper)
- Extract: why a stiff tyre spring combined with a soft gas spring makes the system stiff, and choosing an implicit integrator; an energy budget (kinetic + gravitational + gas-spring = initial − damper − tyre losses).
- GAP: no verified free source covers constrained multibody formulations (constraint forces, Lagrange multipliers, index reduction) or oleo-pneumatic strut modelling (polytropic gas spring p V^n = const, orifice damping proportional to ẋ|ẋ|, drop-test load factor). The model and its numbers are original and labelled *illustrative*. Candidates (unverified): NACA Report 1154, Milwitzky and Cook, "Analysis of Landing-Gear Behavior" (1953); NACA-TN-2743, Flügge, "Landing-gear impact" (1952). Both appeared in NTRS search results; neither page was loaded.

## 06-19 Home Lab: Pendulum System Identification

Sources: `@openstax-university-physics-vol1`, `@mit-ocw-18-03-differential-equations`, `@jcgm-100-2008`, `@openstax-introductory-statistics-2e`, `@scipy-userguide`, `@mit-ocw-6-003-signals-and-systems`

- OpenStax Vol. 1 Ch. 15: 15.4 Pendulums; 15.5 Damped Oscillations
  - https://openstax.org/books/university-physics-volume-1/pages/15-5-damped-oscillations (pattern)
- MIT 18.03 notes 12 "Good vibrations, damping conditions"; notes 38 "Linearization near equilibria; the nonlinear pendulum"
- JCGM 100:2008 (GUM): Type A and Type B evaluation; combined standard uncertainty; coverage factor
  - https://www.bipm.org/documents/20126/2071204/JCGM_100_2008_E.pdf (pattern)
- OpenStax Introductory Statistics 2e: 12.3 The Regression Equation (residuals); 12.6 Outliers
- SciPy "Optimization (scipy.optimize)" (nonlinear least squares). The `curve_fit` reference page was not recorded in an earlier map, so cite the chapter only.
- MIT 6.003 Lecture 21: Sampling; Lecture 22: Sampling and Quantization (phone frame rates and IMU sample rates; aliasing of the swing)
- Extract: θ(t) = A e^(−ζω_n t) cos(ω_d t + φ) as the model; logarithmic decrement compared with a least-squares fit; parameter uncertainty from the covariance compared with the scatter of repeat runs; GUM-style reporting.
- Partial GAP: no verified free source distinguishes viscous, Coulomb (linear in time) and quadratic (amplitude-dependent) damping from a decay envelope, or covers bootstrap uncertainty for curve-fit parameters. These are original. The procedure, bill of materials and safety box live in `projects/p06-pendulum-sysid/` and are original. Candidates (unverified; neither tool's page or licence could be loaded): phyphox (RWTH Aachen University; a search result reported the code as GPL-3.0), https://phyphox.org/; Tracker video analysis (Open Source Physics; a search result reported GPL-3.0), https://physlets.org/tracker/. Until these are verified, name any free sensor-logging or video-tracking app neutrally, as the author brief requires.


## Optional textbooks for deeper study

Paid; never required. Not recorded in `sources.yaml`.

- J. L. Meriam, L. G. Kraige and J. N. Bolton, *Engineering Mechanics: Dynamics* (Wiley): the standard problem sets for particle and planar rigid-body dynamics (06-01 to 06-06), and the main reference for MIT 16.07.
- R. C. Hibbeler, *Engineering Mechanics: Dynamics* (Pearson): an alternative to Meriam for 06-01 to 06-06, with many worked linkage and rolling examples.
- D. T. Greenwood, *Principles of Dynamics* and *Advanced Dynamics* (Prentice Hall / Cambridge University Press): rotating frames, rigid-body 3-D dynamics, Lagrange and Hamilton (06-04, 06-09 to 06-12).
- D. J. Inman, *Engineering Vibration* (Pearson): SDOF, MDOF, continuous systems and isolation (06-13 to 06-17).
- W. T. Thomson and M. D. Dahleh, *Theory of Vibration with Applications* (Prentice Hall): a classic alternative to Inman, with random vibration (06-13 to 06-17).
- J. R. Wertz (ed.), *Spacecraft Attitude Determination and Control* (Kluwer/Springer): attitude representations and spin dynamics (06-07 to 06-10).
- F. L. Markley and J. L. Crassidis, *Fundamentals of Spacecraft Attitude Determination and Control* (Springer): the modern reference for quaternions and attitude kinematics (06-07, 06-08).
- P. H. Wirsching, T. L. Paez and K. Ortiz, *Random Vibrations: Theory and Practice* (Wiley / Dover): PSDs, response of SDOF and MDOF systems to random input, and Miles' equation (06-16, 06-17).
