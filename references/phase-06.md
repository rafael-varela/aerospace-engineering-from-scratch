# Reference map: Phase 06: Dynamics & Vibrations

**Note 2026-10-17 (update):** the first pass of this map, earlier on 2026-10-17, could load no live page. That pass reused section titles that the phase-01 to phase-05 maps had recorded, and it extended those ids. Once the network opened, the candidates were loaded again with `curl -sSL` through the session's HTTPS proxy. WebFetch still failed with "getaddrinfo ENOTFOUND". HTML pages were reduced to text with a small `html.parser` script; PDFs were read with `pdftotext` or as page images.

**Loaded and added** (nine new ids):
- MIT OCW 16.07 *Dynamics* (Fall 2009): the lecture-notes page, plus the PDFs of L5, L6, L8, L19–L21 and L24–L30.
- MIT OCW 2.003SC *Engineering Dynamics* (Fall 2011): the resource index, plus the lecture pages for Lectures 1, 3, 4, 13, 15, 19, 20, 21, 24 and 27.
- MIT OCW 16.61 *Aerospace Dynamics* (Spring 2003): the lecture-notes page, plus the PDFs of lectures 11, 14 and 15. These are handwritten scans, so they are cited by lecture title only.
- NASA-STD-7001C *Payload Vibroacoustic Test Criteria*: the standards.nasa.gov page and the full PDF. Revision C, approved 2026-07-14, supersedes 7001B.
- NTRS 19770024290 (NASA-TM-74839) and NTRS 19930092181 (NACA-TR-1154): the NTRS citation API records and both PDFs.
- LibreTexts, Cline, *Variational Principles in Classical Mechanics*: the book page, the chapter pages for Chapters 3, 4, 6, 8 and 12–16, and four section pages.
- The phyphox home, "Source code" and FAQ pages.
- Tracker: physlets.org/tracker, which redirects to its GitHub Pages site, and that site.

OpenStax *University Physics* Vol. 1 pages 4.4, 6.3, 11.1, 11.4 and 15.4–15.6 were loaded, and `@openstax-university-physics-vol1` now cites 11.1 and the 6.3 Coriolis subsection.

**Rejected:**
- NASA-HDBK-7005 and NASA-HDBK-7004: their standards.nasa.gov pages returned "Access denied" (HTTP 403). The NTRS record 20190026820 is only a 2017 JPL conference-paper record about revising HDBK-7005, with no download. Third-party copies of the handbook are not cited.
- NASA-STD-7001B: superseded by 7001C. The Rev B PDF is listed only under the standard's document history.
- MIT OCW 2.003J (Spring 2007): the lecture-notes URL tried returned 404, and 2.003SC covers the same syllabus.
- NACA-TN-2743 (Flügge, "Landing-gear impact", 1952): its NTRS record was verified as public use permitted, but it is not added, because NACA-TR-1154 covers the lesson and has a verified table of contents.
- LibreTexts Fowler and Tatum: not loaded, because Cline covers every topic they would have covered.

The section titles carried over from the first pass (OpenStax Calculus, 18.03, 18.06, 6.011, 6.003, 16.90, 18.085, Langtangen & Linge, SciPy, SymPy, JCGM) are still copied from the earlier maps and were not reloaded today.

Verified 2026-10-17 by loading the live pages, except where a line says its titles come from an earlier map. Section numbers and titles are copied from the publishers' own tables of contents, section pages or PDF outlines. Every source here is **required-free** unless marked otherwise; sources are identified by an `@`-prefixed identifier and recorded in [sources.yaml](sources.yaml) with license, version, and caveats. **GAP** marks a lesson that no listed free source covers adequately: the lesson must carry original instruction there, and any "candidate" named is unverified until it is added to the manifest.

Conventions: OpenStax section URLs follow `https://openstax.org/books/<book-slug>/pages/<n>-<m>-<title-slug>`; a URL written as `.../pages/...` continues the last full URL above it. Other URL schemes:
- MIT OCW 16.07 lecture PDFs: `https://ocw.mit.edu/courses/16-07-dynamics-fall-2009/resources/mit16_07f09_lec<NN>/`.
- MIT OCW 2.003SC lecture pages: `https://ocw.mit.edu/courses/2-003sc-engineering-dynamics-fall-2011/resources/lecture-<n>-<slug>-1/` (slugs as listed in the resource index).
- MIT OCW 16.61 lectures: `https://ocw.mit.edu/courses/16-61-aerospace-dynamics-spring-2003/resources/lecture<n>/`. Lecture 1 is `lecture1a` and lecture 15 is `lecture_15`.
- Cline chapters: `https://phys.libretexts.org/Bookshelves/Classical_Mechanics/Variational_Principles_in_Classical_Mechanics_(Cline)/<NN>%3A_<Chapter_Title>`.

"(loaded)" means that URL was loaded on 2026-10-17. "(loaded YYYY-MM-DD)" means an earlier map records it as loaded on that date. "(pattern)" means the URL was built from a verified scheme. Licences:
- **CC BY-NC-SA 4.0:** MIT OCW, OpenStax and LibreTexts (Cline).
- **CC BY 4.0:** Langtangen & Linge.
- **US Government works:** NASA-STD-7001C (marked "Approved for public release—distribution is unlimited"), and the NTRS reports (determination "GOV_PUBLIC_USE_PERMITTED"). NASA *Basics of Space Flight* is also a US Government work.
- **GNU GPL:** phyphox source code (version not stated on the page; the name is a registered trademark) and Tracker (GPL version 3).
- **BSD 3-Clause:** SciPy and SymPy.
- Free to read, not to redistribute: JCGM 100:2008.

The curriculum links to these sources and writes its own explanations; nothing is copied or adapted.


## 06-01 Particle Kinematics in Many Coordinates

Sources: `@mit-ocw-16-07-dynamics`, `@mit-ocw-2-003sc-engineering-dynamics`, `@openstax-university-physics-vol1`, `@openstax-calculus-volume-3`, `@openstax-calculus-volume-2`

- MIT 16.07 L4 "Curvilinear motion; Cartesian coordinates"; L5 "Other coordinate systems" (PDF headings: Polar Coordinates (r − θ); Velocity vector; Acceleration vector; Change of basis; Equations of Motion; Cylindrical Coordinates (r − θ − z); Spherical Coordinates (r − θ − φ); Application Examples); L6 "Intrinsic coordinates" (PDF headings: Intrinsic coordinates: Tangential, Normal and Binormal components; Equations of Motion in Intrinsic Coordinates)
  - https://ocw.mit.edu/courses/16-07-dynamics-fall-2009/pages/lecture-notes/ (loaded); .../resources/mit16_07f09_lec05/ (loaded; PDF loaded); .../resources/mit16_07f09_lec06/ (loaded; PDF loaded)
- MIT 2.003SC Lecture 4: Movement of a Particle in Circular Motion w/ Polar Coordinates (video; the page says it covers velocity and acceleration in polar and cylindrical coordinates and the definition of normal and tangential coordinates)
  - https://ocw.mit.edu/courses/2-003sc-engineering-dynamics-fall-2011/resources/lecture-4-movement-of-a-particle-in-circular-motion-w-polar-coordinates-1/ (loaded)
- OpenStax Vol. 1 4.4 Uniform and Nonuniform Circular Motion (subsections: Centripetal Acceleration; Nonuniform Circular Motion; Total Acceleration during Circular Motion)
  - https://openstax.org/books/university-physics-volume-1/pages/4-4-uniform-and-nonuniform-circular-motion (loaded)
- OpenStax Calculus Vol. 3: 3.1 Vector-Valued Functions and Space Curves; 3.4 Motion in Space. Calculus Vol. 2: 7.3 Polar Coordinates (titles from the phase-01 map)
  - https://openstax.org/books/calculus-volume-3/pages/3-4-motion-in-space (pattern)
- Extract: the unit-vector derivatives and the r̈ − rθ̇², rθ̈ + 2ṙθ̇ components (16.07 L5); a_t = v̇ and a_n = v²/ρ (L6); converting between the three descriptions with a change of basis (L5 "Change of basis"). 16.07 L5's helix-descent aircraft example is the aerospace hook.

## 06-02 Particle Kinetics

Sources: `@mit-ocw-16-07-dynamics`, `@openstax-university-physics-vol1`, `@scipy-userguide`

- MIT 16.07 L6 "Intrinsic coordinates": Equations of Motion in Intrinsic Coordinates; L2 "Degrees of freedom and constraints, rectilinear motion"
  - https://ocw.mit.edu/courses/16-07-dynamics-fall-2009/resources/mit16_07f09_lec06/ (loaded; PDF loaded)
- OpenStax Vol. 1 Ch. 6: 6.1 Solving Problems with Newton's Laws; 6.3 Centripetal Force (subsections: Banked Curves; Inertial Forces and Noninertial (Accelerated) Frames: The Coriolis Force); 6.4 Drag Force and Terminal Speed
  - https://openstax.org/books/university-physics-volume-1/pages/6-3-centripetal-force (loaded)
- SciPy "Integration (scipy.integrate)" > "Ordinary differential equations (solve_ivp)" (title from the phase-03 map)
- Extract: ΣF_t = m v̇ and ΣF_n = m v²/ρ; a banked curve as the model for a coordinated turn; integrating a point-mass flight path.
- Partial GAP: no verified source defines aircraft load factor n = L/W, the pull-up relation n = 1 + V²/(gR), or g-onset. The aircraft application is original.

## 06-03 Work-Energy & Impulse-Momentum

Sources: `@mit-ocw-16-07-dynamics`, `@openstax-university-physics-vol1`, `@scipy-userguide`

- MIT 16.07 L9 "Linear impulse and momentum; collisions"; L12 "Work and energy"; L13 "Conservative internal forces and potential energy"
  - https://ocw.mit.edu/courses/16-07-dynamics-fall-2009/resources/mit16_07f09_lec12/ (pattern)
- OpenStax Vol. 1 Ch. 7: 7.1 Work; 7.3 Work-Energy Theorem; Ch. 8: 8.2 Conservative and Non-Conservative Forces; 8.3 Conservation of Energy; Ch. 9: 9.2 Impulse and Collisions (titles from the phase-04 map)
  - https://openstax.org/books/university-physics-volume-1/pages/9-2-impulse-and-collisions (pattern)
- SciPy "Ordinary differential equations (solve_ivp)"
- Extract: W = ΔT along the stopping stroke; ∫F dt = Δp; why a damper whose force depends on velocity needs time integration to give the peak deceleration.
- Partial GAP: no verified source describes carrier arresting gear (the cable and its energy absorber). The force law and its numbers are original and labelled *illustrative*.

## 06-04 Relative Motion & Rotating Frames

Sources: `@mit-ocw-16-07-dynamics`, `@libretexts-cline-variational-principles`, `@mit-ocw-2-003sc-engineering-dynamics`, `@mit-ocw-16-61-aerospace-dynamics`, `@openstax-university-physics-vol1`, `@mit-ocw-18-06-linear-algebra`

- MIT 16.07 L7 "Relative motion using translating axes"; L8 "Relative motion using rotating axes" (PDF headings: Time derivative of a vector in a rotating frame: Coriolis' theorem; Position vector; Velocity vector; Acceleration vector; Acceleration in Rotating Cartesian Coordinates)
  - https://ocw.mit.edu/courses/16-07-dynamics-fall-2009/resources/mit16_07f09_lec08/ (loaded; PDF loaded)
- Cline Ch. 12 "Non-inertial Reference Frames": 12.3 Rotating Reference Frame; 12.4 Reference Frame Undergoing Rotation Plus Translation; 12.5 Newton's Law of Motion in a Non-Inertial Frame; 12.7 Centrifugal Force; 12.8 Coriolis Force; 12.10 Effective gravitational force near the surface of the Earth; 12.11 Free Motion on the Earth (it treats projectiles); 12.13 Foucault pendulum
  - https://phys.libretexts.org/Bookshelves/Classical_Mechanics/Variational_Principles_in_Classical_Mechanics_(Cline)/12%3A_Non-inertial_Reference_Frames (loaded); .../12.11%3A_Free_Motion_on_the_Earth (loaded)
- MIT 2.003SC Lecture 1: History of Dynamics; Motion in Moving Reference Frames; Lecture 3: Motion of Center of Mass; Acceleration in Rotating Ref. Frames; Recitation 3 Notes: V and A of a Point in a Moving Frame (PDF); Lecture 8: Fictitious Forces and Rotating Mass
  - https://ocw.mit.edu/courses/2-003sc-engineering-dynamics-fall-2011/resources/lecture-3-motion-of-center-of-mass-acceleration-in-rotating-ref-1/ (loaded); .../resources/mit2_003scf11_rec3notes/ (pattern)
- MIT 16.61 lecture 2 "Coriolis "Demystified"" (handwritten PDF); lecture 4 "Introduction to Multiple Intermediate Frames"
  - https://ocw.mit.edu/courses/16-61-aerospace-dynamics-spring-2003/resources/lecture2/ (pattern)
- OpenStax Vol. 1 6.3 Centripetal Force, subsection "Inertial Forces and Noninertial (Accelerated) Frames: The Coriolis Force" (qualitative)
  - https://openstax.org/books/university-physics-volume-1/pages/6-3-centripetal-force (loaded)
- Extract: the transport theorem (d/dt)_I = (d/dt)_B + Ω×; the five-term acceleration (16.07 L8); Earth's effective g, deflection of a projectile and the Foucault precession rate Ω sin λ (Cline 12.10, 12.11, 12.13). Check rotating-frame results against inertial-frame integration.

## 06-05 Rigid Body Kinematics in 2D

Sources: `@mit-ocw-16-07-dynamics`, `@scipy-userguide`, `@ocw-16-90`, `@openstax-university-physics-vol1`

- MIT 16.07 L21 "2D rigid body dynamics" (PDF headings: Introduction; Kinematics of Two-Dimensional Rigid Body Motion; Parallel Axes; Body Axes; Instantaneous Center of Rotation; the "Rolling Cylinder" example)
  - https://ocw.mit.edu/courses/16-07-dynamics-fall-2009/resources/mit16_07f09_lec21/ (loaded; PDF loaded)
- SciPy "Optimization (scipy.optimize)" > "Root finding" > "Sets of equations"; MIT 16.90 Unit 1 "1.7.4 Newton-Raphson Implement Implicit Methods on Nonlinear Problems" (titles from the phase-03 map)
- OpenStax Vol. 1 10.3 Relating Angular and Translational Quantities (title from the phase-04 map)
- Extract: v_P = v_O + ω × r_P/O and the matching acceleration; instantaneous centre, space centrode and body centrode (16.07 L21); Newton's method on the loop-closure equations.
- Partial GAP: no verified source treats four-bar linkages (loop-closure position analysis, transmission angle, Grashof classification). The linkage is original.

## 06-06 Rigid Body Kinetics in 2D

Sources: `@mit-ocw-16-07-dynamics`, `@openstax-university-physics-vol1`, `@libretexts-cline-variational-principles`, `@mit-ocw-2-003sc-engineering-dynamics`, `@libretexts-baker-haynes-engineering-statics`

- MIT 16.07 L21 "2D rigid body dynamics" (PDF headings: Conservation of Angular Momentum about the Center of Mass; Radius of Gyration; Parallel Axis Theorem; Summary: Governing Equations; Work; Work done by couples); L22 "2D rigid body dynamics: work and energy"; L23 "2D rigid body dynamics: impulse and momentum"
  - https://ocw.mit.edu/courses/16-07-dynamics-fall-2009/resources/mit16_07f09_lec21/ (loaded; PDF loaded); .../resources/mit16_07f09_lec22/ (pattern)
- OpenStax Vol. 1 11.1 Rolling Motion (subsections: Rolling Motion without Slipping; Rolling Motion with Slipping; Rolling Down an Inclined Plane with Slipping)
  - https://openstax.org/books/university-physics-volume-1/pages/11-1-rolling-motion (loaded)
- Cline 13.24 The Rolling Wheel
- MIT 2.003SC Recitation 4 Notes: Torque and Angular Momentum, Pendulum with Torsional Spring, Rolling Pipe on Moving Truck (PDF)
  - https://ocw.mit.edu/courses/2-003sc-engineering-dynamics-fall-2011/resources/mit2_003scf11_rec4notes1/ (pattern)
- Baker & Haynes 9.1 Dry Friction (9.1.1 Coulomb Friction); 10.8 Mass Moment of Inertia (titles from the phase-05 map)
- Extract: ΣF = m a_G together with ΣM_G = I_G α; kinetic friction while slipping; the end of slip when v_G = ωR (OpenStax 11.1 "Rolling Motion with Slipping"); friction work, checked by 16.07 L22 and by L23.
- Partial GAP: wheel spin-up at landing-gear touchdown (tyre rolling radius, the realistic μ, the friction energy dissipated in the tyre contact patch) is original and labelled *illustrative*.

## 06-07 3D Rotations: Euler Angles & DCMs

Sources: `@ntrs-nasa-tm-74839-euler-angles-quaternions`, `@mit-ocw-16-07-dynamics`, `@libretexts-cline-variational-principles`, `@mit-ocw-16-61-aerospace-dynamics`, `@mit-ocw-18-06-linear-algebra`, `@nasa-basics-of-space-flight`

- NASA-TM-74839 (JSC-12960, 77-FM-37), Henderson, "Euler Angles, Quaternions, and Transformation Matrices: Working Relationships", July 1977: 2.1 Euler Angle Transformation Matrices; Appendix A Relationships for the Three-Axis Euler Angle Rotation Sequences (all twelve sequences, and the Euler angles in terms of the matrix elements); Appendix B Computer Subroutines for the Relationships
  - https://ntrs.nasa.gov/citations/19770024290 (record loaded through the NTRS API; PDF loaded)
- MIT 16.07 L3 "Vectors, matrices and coordinate transformations"; L25 "3D rigid body kinematics" (PDF headings: Rotation About a Fixed Point; Angular Velocity About a Fixed Point; Angular Acceleration; Instantaneous Axis of Rotation; General Motion); L29 "3D rigid body dynamics" (Euler angles as spin, nutation and precession)
  - https://ocw.mit.edu/courses/16-07-dynamics-fall-2009/resources/mit16_07f09_lec25/ (loaded; PDF loaded); .../resources/mit16_07f09_lec29/ (loaded; PDF loaded)
- Cline 13.13 Euler Angles; 13.14 Angular Velocity
- MIT 16.61 lecture 11 "Kinematics of Rigid Bodies" (handwritten PDF)
- MIT 18.06 Lecture 17: Orthogonal matrices and Gram-Schmidt (title from the phase-01 map; re-orthonormalisation)
- Context: Basics of Space Flight Chapter 11 page 2, AACS (loaded 2026-10-16)
- Extract: building a DCM from three elementary rotations, using all twelve sequences as a check table (TM-74839 App. A); the Euler-angle rate equations from ω (Cline 13.14); gimbal lock where those rates become singular.
- Caveat: TM-74839's matrix M maps vectors from the rotated frame into the original frame (the transpose of the usual "inertial-to-body" DCM), and the scan is a typescript with OCR noise. Check every convention before using it. Partial GAP: integrating the DCM kinematic equation Ċ = −[ω×]C with re-orthonormalisation, and the 3-2-1 aircraft yaw-pitch-roll convention, are original.

## 06-08 Quaternions

Sources: `@ntrs-nasa-tm-74839-euler-angles-quaternions`, `@mit-ocw-18-06-linear-algebra`, `@nasa-basics-of-space-flight`

- NASA-TM-74839: 2.2 Transformation Matrices Using the Hamilton Quaternion; 2.3 Euler Angle and Quaternion Relationships; Appendix A (the quaternion as a function of the Euler angles, for each of the twelve sequences)
  - https://ntrs.nasa.gov/citations/19770024290 (record loaded; PDF loaded). Companion: NASA-CR-151435 (DN-1.4-8-020), "Euler angles, quaternions, and transformation matrices for space shuttle analysis", June 1977, https://ntrs.nasa.gov/citations/19770019231 (record loaded; not added as a separate id)
- MIT 18.06 Lecture 17: Orthogonal matrices and Gram-Schmidt (the DCM side of the quaternion ↔ DCM check)
- Context: Basics of Space Flight Chapter 11 page 2, AACS ("Celestial Reference", "Inertial Reference")
- Extract: the quaternion ↔ DCM ↔ Euler conversions, cross-checked against the lesson 07 library (TM-74839 §2.2–2.3, App. A).
- Partial GAP: quaternion algebra from first principles (the product, the conjugate, rotating a vector by q v q*), the quaternion kinematic equation q̇ = ½ q ⊗ (0, ω) with renormalisation, sign ambiguity (q and −q), and SLERP are not titled sections of any verified source. These are original.

## 06-09 3D Rigid Body Dynamics: Euler's Equations

Sources: `@mit-ocw-16-07-dynamics`, `@libretexts-cline-variational-principles`, `@mit-ocw-2-003sc-engineering-dynamics`, `@mit-ocw-16-61-aerospace-dynamics`, `@openstax-university-physics-vol1`, `@libretexts-baker-haynes-engineering-statics`, `@mit-ocw-18-06-linear-algebra`, `@ocw-16-90`, `@scipy-userguide`

- MIT 16.07 L26 "3D rigid body dynamics: the inertia tensor" (PDF headings: Angular Momentum; Principal Axes of Inertia; Parallel Axis Theorem; Rotation of Axes; Inertia Tensor of a Cube; Products of Inertia in a Plane); L27 "3D rigid body dynamics: kinetic energy, instability, equations of motion" (PDF headings: Kinetic Energy for Systems of Particles; Governing Equations for Rotational Motion of a Three-Dimensional Body); L28 "3D rigid body dynamics: equations of motion; Euler's equations" (PDF headings: Body-Fixed Axis; Stability of Free Motion about a Principal Axis; Stability of a Gyrostat)
  - https://ocw.mit.edu/courses/16-07-dynamics-fall-2009/resources/mit16_07f09_lec28/ (loaded; PDF loaded)
- Cline Ch. 13 "Rigid-body Rotation": 13.4 Inertia Tensor; 13.6 Principal Axis System; 13.11 Angular Momentum and Angular Velocity Vectors; 13.12 Kinetic Energy of Rotating Rigid Body; 13.17 Euler's equations of motion for rigid-body rotation; 13.21 Torque-free rotation of an asymmetric rigid rotor; 13.22 Stability of torque-free rotation of an asymmetric body
  - https://phys.libretexts.org/Bookshelves/Classical_Mechanics/Variational_Principles_in_Classical_Mechanics_(Cline)/13%3A_Rigid-body_Rotation (loaded); .../13.22%3A_Stability_of_torque-free_rotation_of_an_asymmetric_body (loaded)
- MIT 2.003SC Lecture 10: Equations of Motion, Torque, Angular Momentum of Rigid Bodies; Lecture 11: Mass Moment of Inertia of Rigid Bodies; Lecture 13: Four Classes of Problems With Rotational Motion (the page says it covers "an alternative form of Euler's equation")
  - https://ocw.mit.edu/courses/2-003sc-engineering-dynamics-fall-2011/resources/lecture-13-four-classes-of-problems-with-rotational-motion-1/ (loaded)
- MIT 16.61 lecture 12 "Rigid Body Dynamics" (handwritten PDF)
- Earlier-map titles: OpenStax 11.2 Angular Momentum, 11.3 Conservation of Angular Momentum; Baker & Haynes 10.7 Products of Inertia; 18.06 Lecture 25: Symmetric matrices and positive definiteness; 16.90 "1.9 Runge-Kutta Methods"; SciPy solve_ivp
- Extract: Euler's equations in principal axes; the intermediate-axis instability by linearising (16.07 L28 "Stability of Free Motion about a Principal Axis"; Cline 13.22); conservation of T and |H| as the integration check.

## 06-10 Gyroscopes & Spin Stability

Sources: `@mit-ocw-16-07-dynamics`, `@libretexts-cline-variational-principles`, `@mit-ocw-16-61-aerospace-dynamics`, `@openstax-university-physics-vol1`, `@nasa-basics-of-space-flight`

- MIT 16.07 L27 (an Explorer 1 passage: launched 31 January 1958, spin-stabilised about its minimum-inertia axis, flexible antennas dissipated energy, and it tipped over into a spin about the maximum-inertia axis); L28 "Stability of a Gyrostat"; L29 "Description of Free Motions of a Rotating Body Using Euler Angles"; "Extreme Aircraft Dynamics"; L30 "3D rigid body dynamics: tops and gyroscopes" (PDF headings: Steady Precession: Gyroscopic Motion; Unsteady Precession of a Top: Integrals of Motion; Spinning Top by Lagrange's Equation; Extreme Tops: Gyroscopes); L31 "Inertial instruments and inertial navigation"
  - https://ocw.mit.edu/courses/16-07-dynamics-fall-2009/resources/mit16_07f09_lec27/ (loaded; PDF loaded); .../resources/mit16_07f09_lec30/ (loaded; PDF loaded)
- Cline 13.20 Torque-free rotation of an inertially-symmetric rigid rotor; 13.22 Stability of torque-free rotation of an asymmetric body (it mentions Explorer 1 and energy dissipation); 13.23 Symmetric rigid rotor subject to torque about a fixed point
- MIT 16.61 lecture 13 "Axisymmetric Rotations"; lecture 14 "Gyroscopes"; lecture 15 "Spacecraft Attitude Dynamics" (handwritten PDFs; lecture 15 loaded)
  - https://ocw.mit.edu/courses/16-61-aerospace-dynamics-spring-2003/resources/lecture_15/ (loaded; PDF loaded)
- OpenStax Vol. 1 11.4 Precession of a Gyroscope
  - https://openstax.org/books/university-physics-volume-1/pages/11-4-precession-of-a-gyroscope (loaded)
- Basics of Space Flight Chapter 11 page 2, AACS ("Spin", "3-Axis"; loaded 2026-10-16)
- Extract: the body-frame nutation rate of an axisymmetric spinner (Cline 13.20); the major/minor-axis rule with energy dissipation (16.07 L27–L28; Cline 13.22); steady precession Ω = τ/(I ω_s) (16.07 L30; OpenStax 11.4). The Explorer 1 history can now be cited to 16.07 L27.

## 06-11 Lagrangian Mechanics

Sources: `@mit-ocw-16-07-dynamics`, `@mit-ocw-2-003sc-engineering-dynamics`, `@libretexts-cline-variational-principles`, `@mit-ocw-16-61-aerospace-dynamics`, `@mit-ocw-8-033-relativity` (supplementary-free), `@mit-ocw-16-323-principles-of-optimal-control`, `@sympy-tutorial`, `@mit-ocw-18-03-differential-equations`, `@scipy-userguide`

- MIT 16.07 L20 "Energy methods: Lagrange's equations" (PDF headings: Derivation of Lagrange's Equations in Cartesian Coordinates; Mass-spring System; Extension to General Coordinate Systems; Simple Pendulum by Lagrange's Equations; Double Pendulum by Lagrange's Equations; Derivation of Lagrange's Equation for General Coordinate Systems)
  - https://ocw.mit.edu/courses/16-07-dynamics-fall-2009/resources/mit16_07f09_lec20/ (loaded; PDF loaded)
- MIT 2.003SC Lecture 15: Introduction to Lagrange With Examples (generalised coordinates; complete, independent and holonomic); Lecture 16: Kinematic Approach to Finding Generalized Forces; Lecture 17: Practice Finding EOM Using Lagrange Equations; Recitation 8 Notes: Cart and Pendulum (Lagrange) (PDF); Recitation 9 Notes: Generalized Forces with Double Pendulum Example (PDF)
  - https://ocw.mit.edu/courses/2-003sc-engineering-dynamics-fall-2011/resources/lecture-15-introduction-to-lagrange-with-examples-1/ (loaded); .../resources/mit2_003scf11_rec9notes1/ (pattern)
- Cline Ch. 6 "Lagrangian Dynamics": 6.3 Lagrange Equations from d'Alembert's Principle; 6.4 Lagrange equations from Hamilton's Principle; 6.5 Constrained Systems; 6.8 Applications to systems involving holonomic constraints; 6.13 The Lagrangian versus the Newtonian approach to classical mechanics. Ch. 4 "Nonlinear Systems and Chaos": 4.6 Differentiation Between Ordered and Chaotic Motion
  - https://phys.libretexts.org/Bookshelves/Classical_Mechanics/Variational_Principles_in_Classical_Mechanics_(Cline)/06%3A_Lagrangian_Dynamics (loaded)
- MIT 16.61 lecture 7 "Lagrange's Equations"; lecture 8 "Examples Using Lagrange's Equations"; lecture 9 "Virtual Work and the Derivation of Lagrange's Equations"; lecture 10 "Friction in Lagrange's Equations"
- Earlier-map titles: 8.033 Session 8 "Variational Calculus"; 16.323 Lecture 5 "Calculus of variations"; SymPy tutorial "Calculus", "Solvers", "Matrices"; 18.03 notes 38 "Linearization near equilibria; the nonlinear pendulum"
- Extract: generalised coordinates and forces; the double pendulum derived exactly as in 16.07 L20 (the independent check against Newton–Euler is the lesson's own); small-angle normal modes; chaos at large amplitude (Cline 4.6).

## 06-12 Hamiltonian Mechanics ★

Sources: `@libretexts-cline-variational-principles`, `@langtangen-linge-fd-computing-pdes`, `@mit-ocw-18-03-differential-equations`, `@openstax-university-physics-vol1`, `@ocw-18-085`

- Cline Ch. 8 "Hamiltonian Mechanics": 8.2 Legendre Transformation between Lagrangian and Hamiltonian mechanics; 8.3 Hamilton's Equations of Motion; 8.4 Hamiltonian in Different Coordinate Systems; 8.5 Applications of Hamiltonian Dynamics. Ch. 3: 3.4 Geometrical Representations of Dynamical Motion (configuration, state and phase space; plane-pendulum phase portrait). Ch. 15 "Advanced Hamiltonian Mechanics": 15.3 Canonical Transformations in Hamiltonian Mechanics; 15.7 Symplectic Representation. Ch. 11 "Conservative two-body Central Forces"
  - https://phys.libretexts.org/Bookshelves/Classical_Mechanics/Variational_Principles_in_Classical_Mechanics_(Cline)/08%3A_Hamiltonian_Mechanics (loaded); .../03%3A_Linear_Oscillators/3.04%3A_Geometrical_Representations_of_Dynamical_Motion (loaded); .../15.07%3A_Symplectic_Representation (loaded)
- Langtangen & Linge Ch. 1: 1.6 Energy Considerations; 1.7 The Euler-Cromer Method (1.7.4 The Störmer-Verlet Algorithm) (titles from the phase-03 map)
- 18.03 notes 35, 38; OpenStax 13.3, 13.5; 18.085 Lecture 10 (titles from the earlier maps)
- Extract: H = T + V through the Legendre transform; Hamilton's equations as a first-order system; pendulum and Kepler phase portraits; bounded energy error under Störmer–Verlet.
- Partial GAP: no verified section is titled for Liouville's theorem (none of the section pages loaded contained the word), so phase-space area preservation, and its numerical demonstration with symplectic and non-symplectic integrators, are original.

## 06-13 Single-DOF Vibration

Sources: `@mit-ocw-2-003sc-engineering-dynamics`, `@mit-ocw-16-07-dynamics`, `@libretexts-cline-variational-principles`, `@mit-ocw-18-03-differential-equations`, `@openstax-university-physics-vol1`, `@langtangen-linge-fd-computing-pdes`, `@scipy-userguide`

- MIT 2.003SC Lecture 19: Introduction to Mechanical Vibration (the page lists SDOF free decay, natural frequencies and damping ratios); Lecture 20: Linear System Modeling a Single Degree of Freedom Oscillator (the page lists experimental determination of damping ratio and the steady-state response to harmonic input); Recitation 10 Notes: Natural Frequency From Deflection & Frequency Response (PDF)
  - https://ocw.mit.edu/courses/2-003sc-engineering-dynamics-fall-2011/resources/lecture-19-introduction-to-mechanical-vibration-1/ (loaded); .../resources/lecture-20-linear-system-modeling-a-single-degree-of-freedom-oscillator-1/ (loaded)
- MIT 16.07 L19 "Vibration, normal modes, natural frequencies, instability" (PDF heading: Vibration, Instability)
- Cline 3.5 Linearly-damped Free Linear Oscillator; 3.6 Sinusoidally-driven, linearly-damped, linear oscillator
  - https://phys.libretexts.org/Bookshelves/Classical_Mechanics/Variational_Principles_in_Classical_Mechanics_(Cline)/03%3A_Linear_Oscillators (loaded)
- OpenStax Vol. 1 15.5 Damped Oscillations; 15.6 Forced Oscillations
  - https://openstax.org/books/university-physics-volume-1/pages/15-6-forced-oscillations (loaded)
- 18.03 notes 11–16 and L9, L10, L13, L14; Langtangen & Linge 1.4.2 (titles from the earlier maps)
- Extract: ω_n, ζ and the logarithmic decrement (2.003SC L20, "experimental determination of damping ratio"); complex frequency response (2.003SC L20, Recitation 10).
- Partial GAP: base excitation and the transmissibility ratio are not named on the loaded 2.003SC pages for Lectures 19 and 20. Derive them here; 06-17 cites 2.003SC Lecture 21 for isolation.

## 06-14 Multi-DOF Vibration & Modes

Sources: `@mit-ocw-2-003sc-engineering-dynamics`, `@mit-ocw-16-07-dynamics`, `@libretexts-cline-variational-principles`, `@mit-ocw-18-03-differential-equations`, `@mit-ocw-18-06-linear-algebra`, `@scipy-userguide`

- MIT 2.003SC Lecture 22: Finding Natural Frequencies and Mode Shapes of a 2 DOF System; Lecture 23: Vibration by Mode Superposition; Lecture 24: Modal Analysis: Orthogonality, Mass Stiffness, Damping Matrix (the page lists the modal expansion theorem, uncoupled equations, modal initial conditions, Rayleigh damping); Lecture 25: Modal Analysis: Response to IC's and to Harmonic Forces; Lecture 26: Response of 2-DOF Systems by the Use of Transfer Functions; Recitation 12 Notes: Modal Analysis (PDF)
  - https://ocw.mit.edu/courses/2-003sc-engineering-dynamics-fall-2011/resources/lecture-24-modal-analysis-orthogonality-mass-stiffness-damping-matrix-1/ (loaded)
- MIT 16.07 L19 (PDF headings: Vibration of Multi-Degree of Freedom Systems; Normal Modes of Multi-Degree of Freedom Systems; Characteristic Value Problem; Solution for Unequal Masses)
  - https://ocw.mit.edu/courses/16-07-dynamics-fall-2009/resources/mit16_07f09_lec19/ (loaded; PDF loaded)
- Cline Ch. 14 "Coupled Linear Oscillators": 14.3 Normal Modes; 14.6 General Analytic Theory for Coupled Linear Oscillators; 14.8 Three-body coupled linear oscillator systems; 14.11 Damped Coupled Linear Oscillators
  - https://phys.libretexts.org/Bookshelves/Classical_Mechanics/Variational_Principles_in_Classical_Mechanics_(Cline)/14%3A_Coupled_Linear_Oscillators (loaded)
- 18.03 L30 and notes 36; 18.06 Lectures 21, 22, 25; SciPy "Linear Algebra (scipy.linalg)" (titles from the earlier maps)
- Extract: K φ = ω² M φ; mass normalisation and orthogonality; modal superposition of forced response; Rayleigh damping (2.003SC L24). Cline 14.8 is the three-mass case.

## 06-15 Continuous Systems

Sources: `@mit-ocw-2-003sc-engineering-dynamics`, `@libretexts-cline-variational-principles`, `@ocw-18-085`, `@ocw-16-90`, `@mit-ocw-3-11-mechanics-of-materials`, `@mit-ocw-18-03-differential-equations`, `@scipy-userguide`

- MIT 2.003SC Lecture 27: Vibration of Continuous Structures: Strings, Beams, Rods, etc. (the page lists the wave equation for strings and rods, and the vibration of beams with dispersion in wave propagation)
  - https://ocw.mit.edu/courses/2-003sc-engineering-dynamics-fall-2011/resources/lecture-27-vibration-of-continuous-structures-strings-beams-rods-etc-1/ (loaded)
- Cline 14.10 Discrete Lattice Chain; 16.2 The Continuous Uniform Linear Chain (the limit from a discrete model to a continuous one)
  - https://phys.libretexts.org/Bookshelves/Classical_Mechanics/Variational_Principles_in_Classical_Mechanics_(Cline)/16%3A_Analytical_Formulations_for_Continuous_Systems (loaded)
- Earlier-map titles: 18.085 Lectures 17–20 ("Lecture 20: Element Matrices; 4th Order Bending Equations"); 16.90 "2.9 Introduction to Finite Elements"; Roylance "Beam Displacements" (loaded 2026-10-17); 18.03 notes 20; SciPy "Sparse eigenvalue problems with ARPACK"
- Extract: the bending wave equation and its dispersion (2.003SC L27); refining a discrete chain toward the continuum (Cline 16.2); fourth-order element matrices (18.085 L20).
- Partial GAP: no verified titled section gives the cantilever frequency equation cos βL cosh βL = −1, the pinned–pinned frequencies (nπ/L)²√(EI/ρA), or the consistent mass matrix of a beam element. These are original.

## 06-16 Random Vibration

Sources: `@nasa-std-7001c-payload-vibroacoustic-test-criteria`, `@mit-ocw-6-011-signals-systems-and-inference`, `@mit-ocw-6-003-signals-and-systems`, `@mit-ocw-18-03-differential-equations`, `@scipy-userguide`

- NASA-STD-7001C, *Payload Vibroacoustic Test Criteria* (approved 2026-07-14; supersedes 7001B): 4.1 Summary of Verification Test Requirements (Table 1—Summary of Verification Test Requirement Levels and Durations); 4.3 Test Levels (4.3.1 Qualification and Protoflight Tests; 4.3.2 Acceptance Testing; 4.3.3 Workmanship); 4.4.2 Random Vibration Tests; 4.4.3 Test Duration; 4.4.4 Test Control Tolerances; Appendix C.1.3 Mass Attenuation of the Minimum Workmanship Vibration Level
  - https://standards.nasa.gov/standard/NASA/NASA-STD-7001 (loaded); PDF "NASA-STD-7001C - Final 7-14-2026.pdf" linked from that page (loaded)
  - Table 2—Component Minimum Workmanship Random Vibration Test Levels (recorded from the PDF): 20 Hz @ 0.01 g²/Hz; 20 to 80 Hz @ +3 dB/oct; 80 to 500 Hz @ 0.04 g²/Hz; 500 to 2000 Hz @ −3 dB/oct; 2000 Hz @ 0.01 g²/Hz; overall 6.8 grms (components under 50 kg). Table 1: qualification = flight acceptance + 3 dB; minimum durations 2 min (qualification, single mission), 1 min (protoflight, flight acceptance). Lesson authors may use these levels with the citation.
- 6.011 Lectures 16–19; 6.003 Lecture 21: Sampling; 18.03 notes 16 (titles from the phase-02 map)
- Extract: integrating a log–log PSD to Grms (the Table 2 spectrum should integrate to about 6.8 grms, a built-in check); S_yy = |H|² S_xx; test levels relative to the maximum expected flight level.
- Partial GAP: NASA-STD-7001C does not state Miles' equation or the 3σ peak convention (no occurrence of "Miles" in the text). Both are original.

## 06-17 Vibration Isolation

Sources: `@mit-ocw-2-003sc-engineering-dynamics`, `@nasa-std-7001c-payload-vibroacoustic-test-criteria`, `@mit-ocw-18-03-differential-equations`, `@mit-ocw-6-011-signals-systems-and-inference`

- MIT 2.003SC Lecture 21: Vibration Isolation (the page lists the Hx/F transfer function with complex numbers, vibration isolation system design, predicting natural frequency by √(g/δ), and isolation of the source); Recitation 10 Notes: Natural Frequency From Deflection & Frequency Response (PDF)
  - https://ocw.mit.edu/courses/2-003sc-engineering-dynamics-fall-2011/resources/lecture-21-vibration-isolation-1/ (loaded); .../resources/mit2_003scf11_rec10note1/ (pattern)
- NASA-STD-7001C 4.3.3 Workmanship (Table 2 as the input spectrum for an IMU-sized component); 4.4.6 Test Tailoring Methods
- 18.03 notes 14, 16, 30; 6.011 Lecture 17 (titles from the earlier maps)
- Extract: transmissibility and the √2 frequency ratio; static deflection fixing f_n; how damping trades the resonant peak against attenuation at high frequency; the random response to the Table 2 input.
- Partial GAP: sway space and IMU-specific isolator sizing are original and labelled *illustrative*.

## 06-18 Multibody Simulation

Sources: `@ntrs-naca-tr-1154-landing-gear-behavior`, `@mit-ocw-16-61-aerospace-dynamics`, `@ocw-16-90`, `@scipy-userguide`, `@openstax-university-physics-vol1`

- NACA-TR-1154, Milwitzky and Cook, "Analysis of Landing-Gear Behavior" (1953), contents:
  - Mechanics of Landing Gear (Dynamics of System; Forces in Shock Strut: Hydraulic force, Pneumatic force, Internal friction force; Forces on Tire).
  - Equations of Motion (Motion Prior to Shock-Strut Deflection; Motion Subsequent to Beginning of Shock-Strut Deflection).
  - Solution of Equations of Motion (Numerical Integration Procedures; Use of Tire Force-Deflection Characteristics; Effect of Drag Loads).
  - Evaluation of Analysis by Comparison of Calculated Results with Experimental Data (Normal Impact; Impact With Tire Bottoming).
  - Parameter Studies (Representation of Tire Force-Deflection Characteristics; Effect of Orifice Discharge Coefficient; Effect of Air-Compression Process).
  - Simplification of Equations of Motion.
  - Appendix A—Numerical Integration Procedures (Runge-Kutta Procedure).
  - Appendix B—Source of Experimental Data.
  - URLs: https://ntrs.nasa.gov/citations/19930092181 (record loaded through the NTRS API; PDF loaded, contents page read)
- MIT 16.61 lecture 6 "Numerical Solution of Nonlinear Differential Equations"; lecture 10 "Friction in Lagrange's Equations"
- Earlier-map titles: 16.90 "1.7 Stiffness and Implicit Methods"; SciPy solve_ivp; OpenStax 8.2
- Extract: the two-mass landing-gear model (upper mass and lower mass), with a pneumatic force from a polytropic air spring and a hydraulic force ∝ ẋ² through an orifice; tyre force–deflection; checking a drop test against its data (TR-1154 "Evaluation of Analysis"); an energy budget; a stiff integrator.
- Partial GAP: constraint-force formulations for general multibody systems (Lagrange multipliers, index reduction) are not a titled section of any verified source. They are original.

## 06-19 Home Lab: Pendulum System Identification

Sources: `@openstax-university-physics-vol1`, `@mit-ocw-16-07-dynamics`, `@mit-ocw-2-003sc-engineering-dynamics`, `@mit-ocw-18-03-differential-equations`, `@jcgm-100-2008`, `@openstax-introductory-statistics-2e`, `@scipy-userguide`, `@mit-ocw-6-003-signals-and-systems`, `@phyphox-app` (supplementary-free), `@tracker-video-analysis` (supplementary-free)

- OpenStax Vol. 1 15.4 Pendulums (subsections: The Simple Pendulum; Physical Pendulum; Torsional Pendulum); 15.5 Damped Oscillations
  - https://openstax.org/books/university-physics-volume-1/pages/15-4-pendulums (loaded); .../pages/15-5-damped-oscillations (loaded)
- MIT 16.07 L24 "Pendulums" (PDF headings: Simple Pendulum; Small amplitude approximation; Compound Pendulum; Center of Percussion; Reversibility)
  - https://ocw.mit.edu/courses/16-07-dynamics-fall-2009/resources/mit16_07f09_lec24/ (loaded; PDF loaded)
- MIT 2.003SC Lecture 20 (experimental determination of damping ratio)
- Earlier-map titles: 18.03 notes 12 and 38; JCGM 100:2008; Introductory Statistics 12.3 and 12.6; SciPy "Optimization (scipy.optimize)"; 6.003 Lectures 21–22
- Logging tools (supplementary-free): phyphox (RWTH Aachen University), "Download for free" on the home page, with a FAQ section "Export formats"; Tracker video analysis and modeling tool (Open Source Physics), GPL v3
  - https://phyphox.org/ (loaded); https://phyphox.org/source/ (loaded; licence statement); https://phyphox.org/faq/ (loaded); https://opensourcephysics.github.io/tracker-website/ (loaded; physlets.org/tracker redirects here)
- Extract: a physical pendulum model (16.07 L24, OpenStax 15.4); the decay envelope and the logarithmic decrement; a curve fit with covariance-based uncertainty and GUM reporting; logging the phone's gyroscope (phyphox) or tracking a video (Tracker).
- Partial GAP: no verified source tells viscous, Coulomb and quadratic damping apart by the shape of the decay envelope, or covers bootstrap uncertainty for fit parameters. These are original. The procedure, bill of materials and safety box live in `projects/p06-pendulum-sysid/` and are original. The lesson may name phyphox or Tracker as one option, but must not require either one.


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
