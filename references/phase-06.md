# Reference map: Phase 06: Dynamics & Vibrations

**Note 2026-10-17 (third pass):** Pages were loaded with curl through the session proxy because WebFetch could not resolve hosts in this environment. Every source cited below was downloaded on 2026-10-17 with `curl -sSL --max-time 60`, into a scratch folder used only for downloads. HTML was reduced to text and headings with an `html.parser` script, and PDFs with `pdftotext`. Scanned PDFs were read as page images. Every section number, title and licence line below was read from those downloaded files; no title is carried over from an earlier map.

**Reloaded in this pass** (titles from the earlier phase maps now replaced by what the pages show):
- OpenStax *University Physics* Vol. 1: chapter-introduction pages for Ch. 4, 6–11, 13 and 15, plus 25 section pages.
- OpenStax *Calculus* Vol. 2 7.3 and Vol. 3 3.1 and 3.4.
- OpenStax *Introductory Statistics 2e* 12.3 and 12.6.
- MIT OCW: 18.03 (Lecture Notes page and lecture pages), 18.06, 6.003 and 18.085 (video lists), 6.011 (Lecture Notes page), 3.11 (Roylance "Beam Displacements" page and PDF).
- Langtangen & Linge, Ch. 1 chapter page and the book page.
- SciPy User Guide: index, integrate, optimize and linalg pages.
- SymPy Introductory Tutorial.
- The JCGM 100:2008 PDF.
- NASA *Basics of Space Flight* Chapter 11 page 2.
- Baker & Haynes Ch. 10 (LibreTexts) and engineeringstatics.org §9.1.

**Loaded in the second pass the same day:**
- MIT OCW 16.07 (lecture list and 13 lecture PDFs).
- MIT OCW 2.003SC (resource index and lecture pages).
- MIT OCW 16.61 (lecture list and three PDFs).
- Cline on LibreTexts.
- NASA-STD-7001C (page and PDF).
- NTRS 19770024290 and 19930092181 (API records and PDFs).
- The phyphox and Tracker sites.

**Not verified, and not cited:**
- MIT OCW 16.90: the course URL in the manifest returned 404. The DSpace record page loaded, but it shows no section titles, and its archive download returned HTTP 429 twice.
- NASA-HDBK-7005 and NASA-HDBK-7004: standards.nasa.gov returned "Access denied" (HTTP 403 from the site, not from the proxy).

**Rejected or removed:**
- NASA-STD-7001B: superseded by 7001C.
- MIT 2.003J: lecture-notes URL 404.
- NACA-TN-2743: record verified but not needed.
- MIT 8.033 and 16.323: removed from 06-11. The 8.033 resource page shows only a file name; the 16.323 "Calculus of Variations" page loaded, but both sources are redundant with 16.07 L20 and Cline Ch. 6.

Verified 2026-10-17 by loading the live pages; section numbers and titles are copied from the publishers' own tables of contents, section pages or PDF outlines. Every source here is **required-free** unless marked otherwise; sources are identified by an `@`-prefixed identifier and recorded in [sources.yaml](sources.yaml) with license, version, and caveats. **GAP** marks a lesson that no listed free source covers adequately: the lesson must carry original instruction there, and any "candidate" named is unverified until it is added to the manifest.

Conventions: OpenStax section URLs follow `https://openstax.org/books/<book-slug>/pages/<n>-<m>-<title-slug>`; a URL written as `.../pages/...` continues the last full URL above it. OpenStax chapter names come from each book's table of contents, which is embedded in the loaded pages. Other URL schemes:
- MIT OCW 16.07 lecture PDFs: `https://ocw.mit.edu/courses/16-07-dynamics-fall-2009/resources/mit16_07f09_lec<NN>/`.
- MIT OCW 2.003SC lecture pages: `https://ocw.mit.edu/courses/2-003sc-engineering-dynamics-fall-2011/resources/lecture-<n>-<slug>-1/`.
- MIT OCW 16.61 lectures: `https://ocw.mit.edu/courses/16-61-aerospace-dynamics-spring-2003/resources/lecture<n>/`. Lecture 1 is `lecture1a` and lecture 15 is `lecture_15`.
- Cline chapters: `https://phys.libretexts.org/Bookshelves/Classical_Mechanics/Variational_Principles_in_Classical_Mechanics_(Cline)/<NN>%3A_<Chapter_Title>`.

"(loaded)" means that URL was downloaded on 2026-10-17. "(pattern)" means it was built from a verified scheme and not downloaded. "(listed on loaded page)" means the title and link were read on a loaded index page. Licences, as stated on the loaded pages:
- **CC BY-NC-SA 4.0:** MIT OCW (footer link to by-nc-sa/4.0), OpenStax ("licensed under a Creative Commons Attribution-NonCommercial-ShareAlike License") and LibreTexts.
- **CC BY 4.0:** Langtangen & Linge ("This book is open access under a CC BY 4.0 license").
- **US Government works:** NASA-STD-7001C ("APPROVED FOR PUBLIC RELEASE—DISTRIBUTION IS UNLIMITED") and the NTRS reports ("GOV_PUBLIC_USE_PERMITTED").
- **GNU GPL:** phyphox source code and Tracker (GPL version 3).
- **Copyright JCGM member organizations, free to read:** JCGM 100:2008 ("Copyright of this document is shared jointly by the JCGM member organizations").
- **BSD 3-Clause:** SciPy and SymPy (licence recorded in the manifest; not stated on the docs pages).

The curriculum links to these sources and writes its own explanations; nothing is copied or adapted.


## 06-01 Particle Kinematics in Many Coordinates

Sources: `@mit-ocw-16-07-dynamics`, `@mit-ocw-2-003sc-engineering-dynamics`, `@openstax-university-physics-vol1`, `@openstax-calculus-volume-3`, `@openstax-calculus-volume-2`

- MIT 16.07 L4 "Curvilinear motion; Cartesian coordinates"; L5 "Other coordinate systems" (PDF headings: Polar Coordinates (r − θ); Velocity vector; Acceleration vector; Change of basis; Equations of Motion; Cylindrical Coordinates (r − θ − z); Spherical Coordinates (r − θ − φ); Application Examples); L6 "Intrinsic coordinates" (PDF headings: Intrinsic coordinates: Tangential, Normal and Binormal components; Equations of Motion in Intrinsic Coordinates)
  - https://ocw.mit.edu/courses/16-07-dynamics-fall-2009/pages/lecture-notes/ (loaded); .../resources/mit16_07f09_lec05/ (loaded; PDF loaded); .../resources/mit16_07f09_lec06/ (loaded; PDF loaded)
- MIT 2.003SC Lecture 4: Movement of a Particle in Circular Motion w/ Polar Coordinates (video; the page says it covers velocity and acceleration in polar and cylindrical coordinates and the definition of normal and tangential coordinates)
  - https://ocw.mit.edu/courses/2-003sc-engineering-dynamics-fall-2011/resources/lecture-4-movement-of-a-particle-in-circular-motion-w-polar-coordinates-1/ (loaded)
- OpenStax University Physics Vol. 1 Ch. 4 "Motion in Two and Three Dimensions": 4.4 Uniform and Nonuniform Circular Motion (subsections: Centripetal Acceleration; Equations of Motion for Uniform Circular Motion; Nonuniform Circular Motion; worked example "Total Acceleration during Circular Motion")
  - https://openstax.org/books/university-physics-volume-1/pages/4-4-uniform-and-nonuniform-circular-motion (loaded)
- OpenStax Calculus Vol. 3 Ch. 3 "Vector-Valued Functions": 3.1 Vector-Valued Functions and Space Curves; 3.4 Motion in Space (subsections: Motion Vectors in the Plane and in Space; Components of the Acceleration Vector; Projectile Motion; Kepler's Laws; student project "Navigating a Banked Turn")
  - https://openstax.org/books/calculus-volume-3/pages/3-1-vector-valued-functions-and-space-curves (loaded); .../pages/3-4-motion-in-space (loaded)
- OpenStax Calculus Vol. 2 Ch. 7 "Parametric Equations and Polar Coordinates": 7.3 Polar Coordinates (subsections: Defining Polar Coordinates; Polar Curves; Symmetry in Polar Coordinates)
  - https://openstax.org/books/calculus-volume-2/pages/7-3-polar-coordinates (loaded)
- Extract: the unit-vector derivatives and the r̈ − rθ̇², rθ̈ + 2ṙθ̇ components (16.07 L5); a_t = v̇ and a_n = v²/ρ (16.07 L6; Calculus Vol. 3 3.4 "Components of the Acceleration Vector"); converting between the three descriptions with a change of basis (16.07 L5).

## 06-02 Particle Kinetics

Sources: `@mit-ocw-16-07-dynamics`, `@openstax-university-physics-vol1`, `@scipy-userguide`

- MIT 16.07 L6 "Intrinsic coordinates": Equations of Motion in Intrinsic Coordinates; L2 "Degrees of freedom and constraints, rectilinear motion"
  - https://ocw.mit.edu/courses/16-07-dynamics-fall-2009/resources/mit16_07f09_lec06/ (loaded; PDF loaded)
- OpenStax Vol. 1 Ch. 6 "Applications of Newton's Laws": 6.1 Solving Problems with Newton's Laws (subsections: Problem-Solving Strategies; Particle Equilibrium; Particle Acceleration; Newton's Laws of Motion and Kinematics); 6.3 Centripetal Force (subsections: Banked Curves; Inertial Forces and Noninertial (Accelerated) Frames: The Coriolis Force); 6.4 Drag Force and Terminal Speed
  - https://openstax.org/books/university-physics-volume-1/pages/6-introduction (loaded); .../pages/6-1-solving-problems-with-newtons-laws (loaded); .../pages/6-3-centripetal-force (loaded); .../pages/6-4-drag-force-and-terminal-speed (loaded)
- SciPy User Guide: "Integration (scipy.integrate)" > "Ordinary differential equations (solve_ivp)"
  - https://docs.scipy.org/doc/scipy/tutorial/integrate.html (loaded; SciPy v1.18.0 Manual)
- Extract: ΣF_t = m v̇ and ΣF_n = m v²/ρ; a banked curve as the model for a coordinated turn; integrating a point-mass flight path.
- Partial GAP: no verified source defines aircraft load factor n = L/W, the pull-up relation n = 1 + V²/(gR), or g-onset. The aircraft application is original.

## 06-03 Work-Energy & Impulse-Momentum

Sources: `@mit-ocw-16-07-dynamics`, `@openstax-university-physics-vol1`, `@scipy-userguide`

- MIT 16.07 L9 "Linear impulse and momentum; collisions"; L12 "Work and energy"; L13 "Conservative internal forces and potential energy" (titles listed on loaded page)
  - https://ocw.mit.edu/courses/16-07-dynamics-fall-2009/pages/lecture-notes/ (loaded)
- OpenStax Vol. 1 Ch. 7 "Work and Kinetic Energy": 7.1 Work (subsections: Work Done by a Force; Work Done by Forces that Vary); 7.3 Work-Energy Theorem (example "Determining a Stopping Force"). Ch. 8 "Potential Energy and Conservation of Energy": 8.2 Conservative and Non-Conservative Forces; 8.3 Conservation of Energy. Ch. 9 "Linear Momentum and Collisions": 9.2 Impulse and Collisions (subsections: Impulse; Impulse-Momentum Theorem; Momentum and Force)
  - https://openstax.org/books/university-physics-volume-1/pages/7-3-work-energy-theorem (loaded); .../pages/8-3-conservation-of-energy (loaded); .../pages/9-2-impulse-and-collisions (loaded)
- SciPy "Ordinary differential equations (solve_ivp)"
- Extract: W = ΔT along the stopping stroke; ∫F dt = Δp; why a damper whose force depends on velocity needs time integration to give the peak deceleration.
- Partial GAP: no verified source describes carrier arresting gear (the cable and its energy absorber). The force law and its numbers are original and labelled *illustrative*.

## 06-04 Relative Motion & Rotating Frames

Sources: `@mit-ocw-16-07-dynamics`, `@libretexts-cline-variational-principles`, `@mit-ocw-2-003sc-engineering-dynamics`, `@mit-ocw-16-61-aerospace-dynamics`, `@openstax-university-physics-vol1`, `@mit-ocw-18-06-linear-algebra`

- MIT 16.07 L7 "Relative motion using translating axes"; L8 "Relative motion using rotating axes" (PDF headings: Time derivative of a vector in a rotating frame: Coriolis' theorem; Position vector; Velocity vector; Acceleration vector; Acceleration in Rotating Cartesian Coordinates)
  - https://ocw.mit.edu/courses/16-07-dynamics-fall-2009/resources/mit16_07f09_lec08/ (loaded; PDF loaded)
- Cline Ch. 12 "Non-inertial Reference Frames": 12.3 Rotating Reference Frame; 12.4 Reference Frame Undergoing Rotation Plus Translation; 12.5 Newton's Law of Motion in a Non-Inertial Frame; 12.7 Centrifugal Force; 12.8 Coriolis Force; 12.10 Effective gravitational force near the surface of the Earth; 12.11 Free Motion on the Earth (it treats projectiles); 12.13 Foucault pendulum
  - https://phys.libretexts.org/Bookshelves/Classical_Mechanics/Variational_Principles_in_Classical_Mechanics_(Cline)/12%3A_Non-inertial_Reference_Frames (loaded); .../12.11%3A_Free_Motion_on_the_Earth (loaded)
- MIT 2.003SC Lecture 1: History of Dynamics; Motion in Moving Reference Frames; Lecture 3: Motion of Center of Mass; Acceleration in Rotating Ref. Frames; Recitation 3 Notes: V and A of a Point in a Moving Frame; Lecture 8: Fictitious Forces and Rotating Mass
  - https://ocw.mit.edu/courses/2-003sc-engineering-dynamics-fall-2011/resources/lecture-3-motion-of-center-of-mass-acceleration-in-rotating-ref-1/ (loaded); .../resources/mit2_003scf11_rec3notes/ (listed on loaded page)
- MIT 16.61 lecture 2 "Coriolis "Demystified""; lecture 4 "Introduction to Multiple Intermediate Frames" (titles listed on loaded page; handwritten PDFs)
  - https://ocw.mit.edu/courses/16-61-aerospace-dynamics-spring-2003/pages/lecture-notes/ (loaded)
- OpenStax Vol. 1 6.3 Centripetal Force, subsection "Inertial Forces and Noninertial (Accelerated) Frames: The Coriolis Force" (qualitative)
  - https://openstax.org/books/university-physics-volume-1/pages/6-3-centripetal-force (loaded)
- MIT 18.06 Lecture 31: Change of basis; image compression (coordinate changes between frames)
  - https://ocw.mit.edu/courses/18-06-linear-algebra-spring-2010/video_galleries/video-lectures/ (loaded; lecture listed there)
- Extract: the transport theorem (d/dt)_I = (d/dt)_B + Ω×; the five-term acceleration (16.07 L8); Earth's effective g, deflection of a projectile and the Foucault precession rate (Cline 12.10, 12.11, 12.13). Check rotating-frame results against inertial-frame integration.

## 06-05 Rigid Body Kinematics in 2D

Sources: `@mit-ocw-16-07-dynamics`, `@scipy-userguide`, `@openstax-university-physics-vol1`

- MIT 16.07 L21 "2D rigid body dynamics" (PDF headings: Introduction; Kinematics of Two-Dimensional Rigid Body Motion; Parallel Axes; Body Axes; Instantaneous Center of Rotation; the "Rolling Cylinder" example)
  - https://ocw.mit.edu/courses/16-07-dynamics-fall-2009/resources/mit16_07f09_lec21/ (loaded; PDF loaded)
- SciPy User Guide: "Optimization (scipy.optimize)" > "Root finding" > "Sets of equations"
  - https://docs.scipy.org/doc/scipy/tutorial/optimize.html (loaded)
- OpenStax Vol. 1 Ch. 10 "Fixed-Axis Rotation": 10.3 Relating Angular and Translational Quantities (subsections: Angular vs. Linear Variables; Relationships between Rotational and Translational Motion)
  - https://openstax.org/books/university-physics-volume-1/pages/10-3-relating-angular-and-translational-quantities (loaded)
- Extract: v_P = v_O + ω × r_P/O and the matching acceleration; instantaneous centre, space centrode and body centrode (16.07 L21); Newton's method on the loop-closure equations.
- Partial GAP: no verified source treats four-bar linkages (loop-closure position analysis, transmission angle, Grashof classification). The linkage is original.

## 06-06 Rigid Body Kinetics in 2D

Sources: `@mit-ocw-16-07-dynamics`, `@openstax-university-physics-vol1`, `@libretexts-cline-variational-principles`, `@mit-ocw-2-003sc-engineering-dynamics`, `@libretexts-baker-haynes-engineering-statics`

- MIT 16.07 L21 "2D rigid body dynamics" (PDF headings: Conservation of Angular Momentum about the Center of Mass; Radius of Gyration; Parallel Axis Theorem; Summary: Governing Equations; Work; Work done by couples); L22 "2D rigid body dynamics: work and energy"; L23 "2D rigid body dynamics: impulse and momentum" (L22, L23 titles listed on loaded page)
  - https://ocw.mit.edu/courses/16-07-dynamics-fall-2009/resources/mit16_07f09_lec21/ (loaded; PDF loaded)
- OpenStax Vol. 1 Ch. 11 "Angular Momentum": 11.1 Rolling Motion (subsections: Rolling Motion without Slipping; Rolling Motion with Slipping; Conservation of Mechanical Energy in Rolling Motion; example "Rolling Down an Inclined Plane with Slipping")
  - https://openstax.org/books/university-physics-volume-1/pages/11-1-rolling-motion (loaded)
- Cline 13.24 The Rolling Wheel (listed on the loaded Ch. 13 page)
- MIT 2.003SC Recitation 4 Notes: Torque and Angular Momentum, Pendulum with Torsional Spring, Rolling Pipe on Moving Truck (listed on loaded resource index)
- Baker & Haynes, engineeringstatics.org Section 9.1 Dry Friction (Subsection 9.1.1 Coulomb Friction; 9.1.4 Coulomb Friction Examples, with "Kinetic friction"); LibreTexts Ch. 10 "Moments of Inertia": 10.8 Mass Moment of Inertia
  - https://engineeringstatics.org/Chapter_09-dry-friction.html (loaded); https://eng.libretexts.org/Bookshelves/Mechanical_Engineering/Engineering_Statics%3A_Open_and_Interactive_(Baker_and_Haynes)/10%3A_Moments_of_Inertia (loaded)
- Extract: ΣF = m a_G together with ΣM_G = I_G α; kinetic friction while slipping; the end of slip when v_G = ωR (OpenStax 11.1 "Rolling Motion with Slipping"); friction work, checked by 16.07 L22 and by L23.
- Partial GAP: wheel spin-up at landing-gear touchdown (tyre rolling radius, the realistic μ, the friction energy dissipated in the tyre contact patch) is original and labelled *illustrative*.

## 06-07 3D Rotations: Euler Angles & DCMs

Sources: `@ntrs-nasa-tm-74839-euler-angles-quaternions`, `@mit-ocw-16-07-dynamics`, `@libretexts-cline-variational-principles`, `@mit-ocw-16-61-aerospace-dynamics`, `@mit-ocw-18-06-linear-algebra`, `@nasa-basics-of-space-flight`

- NASA-TM-74839 (JSC-12960, 77-FM-37), Henderson, "Euler Angles, Quaternions, and Transformation Matrices: Working Relationships", July 1977: 2.1 Euler Angle Transformation Matrices; Appendix A Relationships for the Three-Axis Euler Angle Rotation Sequences; Appendix B Computer Subroutines for the Relationships
  - https://ntrs.nasa.gov/citations/19770024290 (record loaded through the NTRS API; PDF loaded)
- MIT 16.07 L3 "Vectors, matrices and coordinate transformations" (title listed on loaded page); L25 "3D rigid body kinematics" (PDF headings: Rotation About a Fixed Point; Angular Velocity About a Fixed Point; Angular Acceleration; Instantaneous Axis of Rotation; General Motion); L29 "3D rigid body dynamics" (Euler angles as spin, nutation and precession)
  - https://ocw.mit.edu/courses/16-07-dynamics-fall-2009/resources/mit16_07f09_lec25/ (loaded; PDF loaded); .../resources/mit16_07f09_lec29/ (loaded; PDF loaded)
- Cline 13.13 Euler Angles; 13.14 Angular Velocity (listed on the loaded Ch. 13 page)
- MIT 16.61 lecture 11 "Kinematics of Rigid Bodies" (handwritten PDF loaded)
- MIT 18.06 Lecture 17: Orthogonal matrices and Gram-Schmidt; Lecture 18: Properties of determinants; Lecture 30: Linear transformations and their matrices
  - https://ocw.mit.edu/courses/18-06-linear-algebra-spring-2010/video_galleries/video-lectures/ (loaded); .../resources/lecture-17-orthogonal-matrices-and-gram-schmidt/ (listed on loaded page)
- Context: NASA *Basics of Space Flight*, "Chapter 11: Onboard Systems" (page 2), section "Attitude and Articulation Control Subsystems (AACS)", run-in paragraphs "Spin", "3-Axis", "Articulation", "Celestial Reference", "Inertial Reference" (qualitative)
  - https://science.nasa.gov/learn/basics-of-space-flight/chapter11-2/ (loaded)
- Extract: building a DCM from three elementary rotations, using all twelve sequences as a check table (TM-74839 App. A); the Euler-angle rate equations from ω (Cline 13.14); gimbal lock where those rates become singular; QᵀQ = I and det = +1 (18.06 L17, L18).
- Caveat: TM-74839's matrix M maps vectors from the rotated frame into the original frame (the transpose of the usual "inertial-to-body" DCM), and the scan is a typescript with OCR noise. Check every convention before using it. Partial GAP: integrating the DCM kinematic equation Ċ = −[ω×]C with re-orthonormalisation, and the 3-2-1 aircraft yaw-pitch-roll convention, are original.

## 06-08 Quaternions

Sources: `@ntrs-nasa-tm-74839-euler-angles-quaternions`, `@mit-ocw-18-06-linear-algebra`, `@nasa-basics-of-space-flight`

- NASA-TM-74839: 2.2 Transformation Matrices Using the Hamilton Quaternion; 2.3 Euler Angle and Quaternion Relationships; Appendix A Relationships for the Three-Axis Euler Angle Rotation Sequences
  - https://ntrs.nasa.gov/citations/19770024290 (record loaded; PDF loaded). Companion: NASA-CR-151435 (DN-1.4-8-020), "Euler angles, quaternions, and transformation matrices for space shuttle analysis", June 1977, https://ntrs.nasa.gov/citations/19770019231 (record loaded; not added as a separate id)
- MIT 18.06 Lecture 17: Orthogonal matrices and Gram-Schmidt (the DCM side of the quaternion ↔ DCM check; listed on the loaded video list)
- Context: *Basics of Space Flight* "Chapter 11: Onboard Systems", "Attitude and Articulation Control Subsystems (AACS)": "Celestial Reference" and "Inertial Reference" paragraphs
  - https://science.nasa.gov/learn/basics-of-space-flight/chapter11-2/ (loaded)
- Extract: the quaternion ↔ DCM ↔ Euler conversions, cross-checked against the lesson 07 library (TM-74839 §2.2–2.3, App. A).
- Partial GAP: quaternion algebra from first principles (the product, the conjugate, rotating a vector by q v q*), the quaternion kinematic equation q̇ = ½ q ⊗ (0, ω) with renormalisation, sign ambiguity (q and −q), and SLERP are not titled sections of any verified source. These are original.

## 06-09 3D Rigid Body Dynamics: Euler's Equations

Sources: `@mit-ocw-16-07-dynamics`, `@libretexts-cline-variational-principles`, `@mit-ocw-2-003sc-engineering-dynamics`, `@mit-ocw-16-61-aerospace-dynamics`, `@openstax-university-physics-vol1`, `@libretexts-baker-haynes-engineering-statics`, `@mit-ocw-18-06-linear-algebra`, `@scipy-userguide`

- MIT 16.07 L26 "3D rigid body dynamics: the inertia tensor" (PDF headings: Angular Momentum; Principal Axes of Inertia; Parallel Axis Theorem; Rotation of Axes; Inertia Tensor of a Cube; Products of Inertia in a Plane); L27 "3D rigid body dynamics: kinetic energy, instability, equations of motion" (PDF headings: Kinetic Energy for Systems of Particles; Governing Equations for Rotational Motion of a Three-Dimensional Body); L28 "3D rigid body dynamics: equations of motion; Euler's equations" (PDF headings: Body-Fixed Axis; Stability of Free Motion about a Principal Axis; Stability of a Gyrostat)
  - https://ocw.mit.edu/courses/16-07-dynamics-fall-2009/resources/mit16_07f09_lec28/ (loaded; PDF loaded)
- Cline Ch. 13 "Rigid-body Rotation": 13.4 Inertia Tensor; 13.6 Principal Axis System; 13.11 Angular Momentum and Angular Velocity Vectors; 13.12 Kinetic Energy of Rotating Rigid Body; 13.17 Euler's equations of motion for rigid-body rotation; 13.21 Torque-free rotation of an asymmetric rigid rotor; 13.22 Stability of torque-free rotation of an asymmetric body
  - https://phys.libretexts.org/Bookshelves/Classical_Mechanics/Variational_Principles_in_Classical_Mechanics_(Cline)/13%3A_Rigid-body_Rotation (loaded); .../13.22%3A_Stability_of_torque-free_rotation_of_an_asymmetric_body (loaded)
- MIT 2.003SC Lecture 10: Equations of Motion, Torque, Angular Momentum of Rigid Bodies; Lecture 11: Mass Moment of Inertia of Rigid Bodies; Lecture 13: Four Classes of Problems With Rotational Motion (the page says it covers "an alternative form of Euler's equation")
  - https://ocw.mit.edu/courses/2-003sc-engineering-dynamics-fall-2011/resources/lecture-13-four-classes-of-problems-with-rotational-motion-1/ (loaded)
- MIT 16.61 lecture 12 "Rigid Body Dynamics" (title listed on loaded page; handwritten PDF)
- OpenStax Vol. 1 Ch. 11 "Angular Momentum": 11.2 Angular Momentum (subsections include Angular Momentum of a Rigid Body); 11.3 Conservation of Angular Momentum
  - https://openstax.org/books/university-physics-volume-1/pages/11-2-angular-momentum (loaded); .../pages/11-3-conservation-of-angular-momentum (loaded)
- Baker & Haynes Ch. 10 "Moments of Inertia": 10.7 Products of Inertia; 10.8 Mass Moment of Inertia (listed on loaded chapter page)
- MIT 18.06 Lecture 25: Symmetric matrices and positive definiteness (listed on loaded video list)
- SciPy "Integration (scipy.integrate)" > "Ordinary differential equations (solve_ivp)"
- Extract: Euler's equations in principal axes; the intermediate-axis instability by linearising (16.07 L28 "Stability of Free Motion about a Principal Axis"; Cline 13.22); conservation of T and |H| as the integration check.

## 06-10 Gyroscopes & Spin Stability

Sources: `@mit-ocw-16-07-dynamics`, `@libretexts-cline-variational-principles`, `@mit-ocw-16-61-aerospace-dynamics`, `@openstax-university-physics-vol1`, `@nasa-basics-of-space-flight`

- MIT 16.07 L27 (an Explorer 1 passage: launched 31 January 1958, spin-stabilised about its minimum-inertia axis, flexible antennas dissipated energy, and it tipped over into a spin about the maximum-inertia axis); L28 "Stability of a Gyrostat"; L29 "Description of Free Motions of a Rotating Body Using Euler Angles"; "Extreme Aircraft Dynamics"; L30 "3D rigid body dynamics: tops and gyroscopes" (PDF headings: Steady Precession: Gyroscopic Motion; Unsteady Precession of a Top: Integrals of Motion; Spinning Top by Lagrange's Equation; Extreme Tops: Gyroscopes); L31 "Inertial instruments and inertial navigation"
  - https://ocw.mit.edu/courses/16-07-dynamics-fall-2009/resources/mit16_07f09_lec27/ (loaded; PDF loaded); .../resources/mit16_07f09_lec30/ (loaded; PDF loaded)
- Cline 13.20 Torque-free rotation of an inertially-symmetric rigid rotor; 13.22 Stability of torque-free rotation of an asymmetric body (it mentions Explorer 1 and energy dissipation); 13.23 Symmetric rigid rotor subject to torque about a fixed point
- MIT 16.61 lecture 13 "Axisymmetric Rotations"; lecture 14 "Gyroscopes"; lecture 15 "Spacecraft Attitude Dynamics" (handwritten PDFs; 14 and 15 loaded)
  - https://ocw.mit.edu/courses/16-61-aerospace-dynamics-spring-2003/resources/lecture_15/ (loaded; PDF loaded)
- OpenStax Vol. 1 11.4 Precession of a Gyroscope (example "Period of Precession")
  - https://openstax.org/books/university-physics-volume-1/pages/11-4-precession-of-a-gyroscope (loaded)
- Context: *Basics of Space Flight* "Chapter 11: Onboard Systems", "Attitude and Articulation Control Subsystems (AACS)": "Spin" and "3-Axis" paragraphs
  - https://science.nasa.gov/learn/basics-of-space-flight/chapter11-2/ (loaded)
- Extract: the body-frame nutation rate of an axisymmetric spinner (Cline 13.20); the major/minor-axis rule with energy dissipation (16.07 L27–L28; Cline 13.22); steady precession Ω = τ/(I ω_s) (16.07 L30; OpenStax 11.4). The Explorer 1 history is cited to 16.07 L27.

## 06-11 Lagrangian Mechanics

Sources: `@mit-ocw-16-07-dynamics`, `@mit-ocw-2-003sc-engineering-dynamics`, `@libretexts-cline-variational-principles`, `@mit-ocw-16-61-aerospace-dynamics`, `@sympy-tutorial`, `@mit-ocw-18-03-differential-equations`, `@scipy-userguide`

- MIT 16.07 L20 "Energy methods: Lagrange's equations" (PDF headings: Derivation of Lagrange's Equations in Cartesian Coordinates; Mass-spring System; Extension to General Coordinate Systems; Simple Pendulum by Lagrange's Equations; Double Pendulum by Lagrange's Equations; Derivation of Lagrange's Equation for General Coordinate Systems)
  - https://ocw.mit.edu/courses/16-07-dynamics-fall-2009/resources/mit16_07f09_lec20/ (loaded; PDF loaded)
- MIT 2.003SC Lecture 15: Introduction to Lagrange With Examples (generalised coordinates; complete, independent and holonomic); Lecture 16: Kinematic Approach to Finding Generalized Forces; Lecture 17: Practice Finding EOM Using Lagrange Equations; Recitation 8 Notes: Cart and Pendulum (Lagrange); Recitation 9 Notes: Generalized Forces with Double Pendulum Example
  - https://ocw.mit.edu/courses/2-003sc-engineering-dynamics-fall-2011/resources/lecture-15-introduction-to-lagrange-with-examples-1/ (loaded); .../pages/resource-index/ (loaded)
- Cline Ch. 6 "Lagrangian Dynamics": 6.3 Lagrange Equations from d'Alembert's Principle; 6.4 Lagrange equations from Hamilton's Principle; 6.5 Constrained Systems; 6.8 Applications to systems involving holonomic constraints; 6.13 The Lagrangian versus the Newtonian approach to classical mechanics. Ch. 4 "Nonlinear Systems and Chaos": 4.6 Differentiation Between Ordered and Chaotic Motion
  - https://phys.libretexts.org/Bookshelves/Classical_Mechanics/Variational_Principles_in_Classical_Mechanics_(Cline)/06%3A_Lagrangian_Dynamics (loaded)
- MIT 16.61 lecture 7 "Lagrange's Equations"; lecture 8 "Examples Using Lagrange's Equations"; lecture 9 "Virtual Work and the Derivation of Lagrange's Equations"; lecture 10 "Friction in Lagrange's Equations" (titles listed on loaded page)
- SymPy 1.14.0 "Introductory Tutorial": Preliminaries; Introduction; Gotchas; "SymPy Features" (Basic Operations; Printing; Simplification; Calculus; Solvers; Matrices; Advanced Expression Manipulation); What's Next
  - https://docs.sympy.org/latest/tutorials/intro-tutorial/index.html (loaded)
- MIT 18.03 lecture notes 37 "Nonlinear systems"; 38 "Linearization near equilibria; the nonlinear pendulum"
  - https://ocw.mit.edu/courses/18-03-differential-equations-spring-2010/pages/lecture-notes/ (loaded)
- SciPy "Integration (scipy.integrate)" > "Ordinary differential equations (solve_ivp)"
- Extract: generalised coordinates and forces; the double pendulum derived exactly as in 16.07 L20 (the independent check against Newton–Euler is the lesson's own); small-angle normal modes (18.03 note 38); chaos at large amplitude (Cline 4.6).

## 06-12 Hamiltonian Mechanics ★

Sources: `@libretexts-cline-variational-principles`, `@langtangen-linge-fd-computing-pdes`, `@mit-ocw-18-03-differential-equations`, `@openstax-university-physics-vol1`, `@ocw-18-085`

- Cline Ch. 8 "Hamiltonian Mechanics": 8.2 Legendre Transformation between Lagrangian and Hamiltonian mechanics; 8.3 Hamilton's Equations of Motion; 8.4 Hamiltonian in Different Coordinate Systems; 8.5 Applications of Hamiltonian Dynamics. Ch. 3 "Linear Oscillators": 3.4 Geometrical Representations of Dynamical Motion (headings: Configuration space; State space; Phase space; Plane pendulum). Ch. 15 "Advanced Hamiltonian Mechanics": 15.3 Canonical Transformations in Hamiltonian Mechanics; 15.7 Symplectic Representation. Ch. 11 "Conservative two-body Central Forces"
  - https://phys.libretexts.org/Bookshelves/Classical_Mechanics/Variational_Principles_in_Classical_Mechanics_(Cline)/08%3A_Hamiltonian_Mechanics (loaded); .../03%3A_Linear_Oscillators/3.04%3A_Geometrical_Representations_of_Dynamical_Motion (loaded); .../15.07%3A_Symplectic_Representation (loaded)
- Langtangen & Linge, *Finite Difference Computing with PDEs*, Ch. 1 "Vibration ODEs": 1.5 Alternative Schemes Based on 1st-Order Equations (1.5.1 The Forward Euler Scheme; 1.5.5 Runge-Kutta Methods); 1.6 Energy Considerations (1.6.1 Derivation of the Energy Expression; 1.6.2 An Error Measure Based on Energy); 1.7 The Euler-Cromer Method (1.7.1 Forward-Backward Discretization; 1.7.4 The Störmer-Verlet Algorithm)
  - https://link.springer.com/chapter/10.1007/978-3-319-55456-3_1 (loaded); https://link.springer.com/book/10.1007/978-3-319-55456-3 (loaded; licence)
- MIT 18.03 lecture notes 35 "Qualitative behavior of linear systems; phase plane"; 38 "Linearization near equilibria; the nonlinear pendulum"; 39 "Limitations of the linear: limit cycles and chaos"
  - https://ocw.mit.edu/courses/18-03-differential-equations-spring-2010/pages/lecture-notes/ (loaded)
- OpenStax Vol. 1 Ch. 13 "Gravitation": 13.3 Gravitational Potential Energy and Total Energy (subsections: Gravitational Potential Energy beyond Earth; Conservation of Energy); 13.5 Kepler's Laws of Planetary Motion (subsections: Kepler's First Law; Orbital Transfers; Kepler's Second Law; Kepler's Third Law)
  - https://openstax.org/books/university-physics-volume-1/pages/13-3-gravitational-potential-energy-and-total-energy (loaded); .../pages/13-5-keplers-laws-of-planetary-motion (loaded)
- MIT 18.085 Lecture 10: Finite Differences in Time
  - https://ocw.mit.edu/courses/18-085-computational-science-and-engineering-i-fall-2008/video_galleries/video-lectures/ (loaded; lecture listed there)
- Extract: H = T + V through the Legendre transform; Hamilton's equations as a first-order system; pendulum and Kepler phase portraits; bounded energy error under Störmer–Verlet.
- Partial GAP: no verified section is titled for Liouville's theorem (none of the section pages loaded contained the word), so phase-space area preservation, and its numerical demonstration with symplectic and non-symplectic integrators, are original.

## 06-13 Single-DOF Vibration

Sources: `@mit-ocw-2-003sc-engineering-dynamics`, `@mit-ocw-16-07-dynamics`, `@libretexts-cline-variational-principles`, `@mit-ocw-18-03-differential-equations`, `@openstax-university-physics-vol1`, `@langtangen-linge-fd-computing-pdes`, `@scipy-userguide`

- MIT 2.003SC Lecture 19: Introduction to Mechanical Vibration (the page lists SDOF free decay, natural frequencies and damping ratios); Lecture 20: Linear System Modeling a Single Degree of Freedom Oscillator (the page lists experimental determination of damping ratio and the steady-state response to harmonic input); Recitation 10 Notes: Natural Frequency From Deflection & Frequency Response
  - https://ocw.mit.edu/courses/2-003sc-engineering-dynamics-fall-2011/resources/lecture-19-introduction-to-mechanical-vibration-1/ (loaded); .../resources/lecture-20-linear-system-modeling-a-single-degree-of-freedom-oscillator-1/ (loaded)
- MIT 16.07 L19 "Vibration, normal modes, natural frequencies, instability" (PDF heading: Vibration, Instability)
- Cline Ch. 3 "Linear Oscillators": 3.5 Linearly-damped Free Linear Oscillator; 3.6 Sinusoidally-driven, linearly-damped, linear oscillator
  - https://phys.libretexts.org/Bookshelves/Classical_Mechanics/Variational_Principles_in_Classical_Mechanics_(Cline)/03%3A_Linear_Oscillators (loaded)
- MIT 18.03 lecture notes 11 "Modes and the characteristic polynomial"; 12 "Good vibrations, damping conditions"; 13 "Exponential response formula, spring drive"; 14 "Complex gain, dashpot drive"; 15 "Operators, undetermined coefficients, resonance"; 16 "Frequency response"; 18.03 Lecture 9: Solving Second-order Linear ODE's with Constant Coefficients; Lecture 10: Continuation: Complex Characteristic Roots; Lecture 13: Finding Particular Solutions to Inhomogeneous ODEs; Lecture 14: Interpretation of the Exceptional Case: Resonance
  - https://ocw.mit.edu/courses/18-03-differential-equations-spring-2010/pages/lecture-notes/ (loaded); .../resources/lecture-14-interpretation-of-the-exceptional-case-resonance/ (loaded); .../video_galleries/video-lectures/ (loaded; Lecture 13 listed there)
- OpenStax Vol. 1 Ch. 15 "Oscillations": 15.5 Damped Oscillations; 15.6 Forced Oscillations
  - https://openstax.org/books/university-physics-volume-1/pages/15-5-damped-oscillations (loaded); .../pages/15-6-forced-oscillations (loaded)
- Langtangen & Linge Ch. 1: 1.4 Analysis of the Numerical Scheme (1.4.2 The Error in the Numerical Frequency)
- SciPy "Integration (scipy.integrate)" > "Ordinary differential equations (solve_ivp)"
- Extract: ω_n, ζ and the logarithmic decrement (2.003SC L20, "experimental determination of damping ratio"); complex frequency response (18.03 notes 14 and 16; 2.003SC L20, Recitation 10).
- Partial GAP: base excitation and the transmissibility ratio are not named on the loaded 2.003SC pages for Lectures 19 and 20. Derive them here; 06-17 cites 2.003SC Lecture 21 for isolation.

## 06-14 Multi-DOF Vibration & Modes

Sources: `@mit-ocw-2-003sc-engineering-dynamics`, `@mit-ocw-16-07-dynamics`, `@libretexts-cline-variational-principles`, `@mit-ocw-18-03-differential-equations`, `@mit-ocw-18-06-linear-algebra`, `@scipy-userguide`

- MIT 2.003SC Lecture 22: Finding Natural Frequencies and Mode Shapes of a 2 DOF System; Lecture 23: Vibration by Mode Superposition; Lecture 24: Modal Analysis: Orthogonality, Mass Stiffness, Damping Matrix (the page lists the modal expansion theorem, uncoupled equations, modal initial conditions, Rayleigh damping); Lecture 25: Modal Analysis: Response to IC's and to Harmonic Forces; Lecture 26: Response of 2-DOF Systems by the Use of Transfer Functions; Recitation 12 Notes: Modal Analysis
  - https://ocw.mit.edu/courses/2-003sc-engineering-dynamics-fall-2011/resources/lecture-24-modal-analysis-orthogonality-mass-stiffness-damping-matrix-1/ (loaded)
- MIT 16.07 L19 (PDF headings: Vibration of Multi-Degree of Freedom Systems; Normal Modes of Multi-Degree of Freedom Systems; Characteristic Value Problem; Solution for Unequal Masses)
  - https://ocw.mit.edu/courses/16-07-dynamics-fall-2009/resources/mit16_07f09_lec19/ (loaded; PDF loaded)
- Cline Ch. 14 "Coupled Linear Oscillators": 14.3 Normal Modes; 14.6 General Analytic Theory for Coupled Linear Oscillators; 14.8 Three-body coupled linear oscillator systems; 14.11 Damped Coupled Linear Oscillators
  - https://phys.libretexts.org/Bookshelves/Classical_Mechanics/Variational_Principles_in_Classical_Mechanics_(Cline)/14%3A_Coupled_Linear_Oscillators (loaded)
- MIT 18.03 Lecture 29: Matrix Exponentials; Lecture 30: Decoupling Linear Systems with Constant Coefficients; lecture notes 36 "Normal modes and the matrix exponential"
  - https://ocw.mit.edu/courses/18-03-differential-equations-spring-2010/resources/lecture-30-decoupling-linear-systems-with-constant-coefficients/ (loaded)
- MIT 18.06 Lecture 21: Eigenvalues and eigenvectors; Lecture 22: Diagonalization and powers of A; Lecture 25: Symmetric matrices and positive definiteness (listed on loaded video list)
- SciPy "Linear Algebra (scipy.linalg)" > "Decompositions" > "Eigenvalues and eigenvectors"
  - https://docs.scipy.org/doc/scipy/tutorial/linalg.html (loaded)
- Extract: K φ = ω² M φ; mass normalisation and orthogonality; modal superposition of forced response; Rayleigh damping (2.003SC L24). Cline 14.8 is the three-mass case.

## 06-15 Continuous Systems

Sources: `@mit-ocw-2-003sc-engineering-dynamics`, `@libretexts-cline-variational-principles`, `@ocw-18-085`, `@mit-ocw-3-11-mechanics-of-materials`, `@scipy-userguide`

- MIT 2.003SC Lecture 27: Vibration of Continuous Structures: Strings, Beams, Rods, etc. (the page lists the wave equation for strings and rods, and the vibration of beams with dispersion in wave propagation)
  - https://ocw.mit.edu/courses/2-003sc-engineering-dynamics-fall-2011/resources/lecture-27-vibration-of-continuous-structures-strings-beams-rods-etc-1/ (loaded)
- Cline 14.10 Discrete Lattice Chain; 16.2 The Continuous Uniform Linear Chain (the limit from a discrete model to a continuous one)
  - https://phys.libretexts.org/Bookshelves/Classical_Mechanics/Variational_Principles_in_Classical_Mechanics_(Cline)/16%3A_Analytical_Formulations_for_Continuous_Systems (loaded)
- MIT 18.085 Lecture 17: Finite Elements in 1D (part 1); Lecture 18: Finite Elements in 1D (part 2); Lecture 19: Quadratic/Cubic Elements; Lecture 20: Element Matrices; 4th Order Bending Equations
  - https://ocw.mit.edu/courses/18-085-computational-science-and-engineering-i-fall-2008/video_galleries/video-lectures/ (loaded; lectures listed there); .../resources/lecture-20-element-matrices-4th-order-bending-equations/ (listed on loaded page)
- Roylance, MIT 3.11 module "Beam Displacements" (PDF dated November 30, 2000): Introduction; Multiple integration; Energy method
  - https://ocw.mit.edu/courses/3-11-mechanics-of-materials-fall-1999/resources/mit3_11f99_bdisp/ (loaded; PDF loaded)
- SciPy "Linear Algebra (scipy.linalg)" > "Decompositions" > "Eigenvalues and eigenvectors"
- Extract: the bending wave equation and its dispersion (2.003SC L27); refining a discrete chain toward the continuum (Cline 16.2); fourth-order element matrices (18.085 L20); EI w'''' from Roylance's "Multiple integration".
- Partial GAP: no verified titled section gives the cantilever frequency equation cos βL cosh βL = −1, the pinned–pinned frequencies (nπ/L)²√(EI/ρA), or the consistent mass matrix of a beam element. These are original.

## 06-16 Random Vibration

Sources: `@nasa-std-7001c-payload-vibroacoustic-test-criteria`, `@mit-ocw-6-011-signals-systems-and-inference`, `@mit-ocw-6-003-signals-and-systems`, `@mit-ocw-18-03-differential-equations`, `@scipy-userguide`

- NASA-STD-7001C, *Payload Vibroacoustic Test Criteria* (approved 2026-07-14; supersedes 7001B): 4.1 Summary of Verification Test Requirements (Table 1—Summary of Verification Test Requirement Levels and Durations); 4.3 Test Levels (4.3.1 Qualification and Protoflight Tests; 4.3.2 Acceptance Testing; 4.3.3 Workmanship); 4.4.2 Random Vibration Tests; 4.4.3 Test Duration; 4.4.4 Test Control Tolerances; Appendix C.1.3 Mass Attenuation of the Minimum Workmanship Vibration Level
  - https://standards.nasa.gov/standard/NASA/NASA-STD-7001 (loaded); PDF "NASA-STD-7001C - Final 7-14-2026.pdf" linked from that page (loaded)
  - Table 2—Component Minimum Workmanship Random Vibration Test Levels (recorded from the PDF): 20 Hz @ 0.01 g²/Hz; 20 to 80 Hz @ +3 dB/oct; 80 to 500 Hz @ 0.04 g²/Hz; 500 to 2000 Hz @ −3 dB/oct; 2000 Hz @ 0.01 g²/Hz; overall 6.8 grms (components under 50 kg). Table 1: qualification = flight acceptance + 3 dB; minimum durations 2 min (qualification, single mission), 1 min (protoflight, flight acceptance). Lesson authors may use these levels with the citation.
- MIT 6.011 lecture notes (resource titles as shown on the page): "6.011 Signals, Systems and Inference, Lecture 16 Wide-sense Stationary Processes, LTI Filtering of WSS Processes"; "… Lecture 17 LTI Filtering of WSS Processes"; "… Lecture 18 Power Spectral Density"; "… Lecture 19 Einstein-Wiener-Khinchin Theorem, PSD Applications, Modeling Filters"
  - https://ocw.mit.edu/courses/6-011-signals-systems-and-inference-spring-2018/pages/lecture-notes/ (loaded)
- MIT 6.003 Lecture 21: Sampling (sampling a synthesised random signal; aliasing)
  - https://ocw.mit.edu/courses/6-003-signals-and-systems-fall-2011/video_galleries/lecture-videos/ (loaded; lecture listed there)
- MIT 18.03 lecture notes 16 "Frequency response" (the |H(ω)| that shapes the output PSD)
- SciPy "Integration (scipy.integrate)" (integrating a PSD for the mean square)
  - https://docs.scipy.org/doc/scipy/tutorial/integrate.html (loaded)
- Extract: integrating a log–log PSD to Grms (the Table 2 spectrum integrates to 6.78 grms, a built-in check); S_yy = |H|² S_xx; test levels relative to the maximum expected flight level.
- Partial GAP: NASA-STD-7001C does not state Miles' equation or the 3σ peak convention (no occurrence of "Miles" in the text). Both are original.

## 06-17 Vibration Isolation

Sources: `@mit-ocw-2-003sc-engineering-dynamics`, `@nasa-std-7001c-payload-vibroacoustic-test-criteria`, `@mit-ocw-18-03-differential-equations`, `@mit-ocw-6-011-signals-systems-and-inference`

- MIT 2.003SC Lecture 21: Vibration Isolation (the page lists the Hx/F transfer function with complex numbers, vibration isolation system design, predicting natural frequency by √(g/δ), and isolation of the source); Recitation 10 Notes: Natural Frequency From Deflection & Frequency Response
  - https://ocw.mit.edu/courses/2-003sc-engineering-dynamics-fall-2011/resources/lecture-21-vibration-isolation-1/ (loaded); .../resources/mit2_003scf11_rec10note1/ (listed on loaded page)
- NASA-STD-7001C 4.3.3 Workmanship (Table 2 as the input spectrum for an IMU-sized component); 4.4.6 Test Tailoring Methods
- MIT 18.03 lecture notes 14 "Complex gain, dashpot drive"; 16 "Frequency response"; 30 "The transfer function and frequency response"
  - https://ocw.mit.edu/courses/18-03-differential-equations-spring-2010/pages/lecture-notes/ (loaded)
- MIT 6.011 "6.011 Signals, Systems and Inference, Lecture 17 LTI Filtering of WSS Processes" (the random-input response that 06-16 builds on)
- Extract: transmissibility and the √2 frequency ratio; static deflection fixing f_n; how damping trades the resonant peak against attenuation at high frequency; the random response to the Table 2 input.
- Partial GAP: sway space and IMU-specific isolator sizing are original and labelled *illustrative*.

## 06-18 Multibody Simulation

Sources: `@ntrs-naca-tr-1154-landing-gear-behavior`, `@mit-ocw-16-07-dynamics`, `@mit-ocw-16-61-aerospace-dynamics`, `@scipy-userguide`, `@openstax-university-physics-vol1`

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
- MIT 16.07 L20 "Energy methods: Lagrange's equations" (generalised coordinates and forces for the two-mass model)
  - https://ocw.mit.edu/courses/16-07-dynamics-fall-2009/resources/mit16_07f09_lec20/ (loaded; PDF loaded)
- MIT 16.61 lecture 6 "Numerical Solution of Nonlinear Differential Equations"; lecture 10 "Friction in Lagrange's Equations" (titles listed on loaded page)
- SciPy "Integration (scipy.integrate)" > "Ordinary differential equations (solve_ivp)", including "Solving a system with a banded Jacobian matrix" (implicit methods for stiff systems)
  - https://docs.scipy.org/doc/scipy/tutorial/integrate.html (loaded)
- OpenStax Vol. 1 Ch. 8: 8.2 Conservative and Non-Conservative Forces (energy accounting with a dissipative damper)
  - https://openstax.org/books/university-physics-volume-1/pages/8-2-conservative-and-non-conservative-forces (loaded)
- Extract: the two-mass landing-gear model (upper mass and lower mass), with a pneumatic force from a polytropic air spring and a hydraulic force ∝ ẋ² through an orifice; tyre force–deflection; checking a drop test against its data (TR-1154 "Evaluation of Analysis"); an energy budget; a stiff integrator.
- Partial GAP: constraint-force formulations for general multibody systems (Lagrange multipliers, index reduction) are not a titled section of any verified source, and no verified source explains stiffness as such (16.90 could not be verified). Both are original.

## 06-19 Home Lab: Pendulum System Identification

Sources: `@openstax-university-physics-vol1`, `@mit-ocw-16-07-dynamics`, `@mit-ocw-2-003sc-engineering-dynamics`, `@mit-ocw-18-03-differential-equations`, `@jcgm-100-2008`, `@openstax-introductory-statistics-2e`, `@scipy-userguide`, `@mit-ocw-6-003-signals-and-systems`, `@phyphox-app` (supplementary-free), `@tracker-video-analysis` (supplementary-free)

- OpenStax Vol. 1 Ch. 15 "Oscillations": 15.4 Pendulums (subsections: The Simple Pendulum; Physical Pendulum; Torsional Pendulum); 15.5 Damped Oscillations
  - https://openstax.org/books/university-physics-volume-1/pages/15-4-pendulums (loaded); .../pages/15-5-damped-oscillations (loaded)
- MIT 16.07 L24 "Pendulums" (PDF headings: Simple Pendulum; Small amplitude approximation; Compound Pendulum; Center of Percussion; Reversibility)
  - https://ocw.mit.edu/courses/16-07-dynamics-fall-2009/resources/mit16_07f09_lec24/ (loaded; PDF loaded)
- MIT 2.003SC Lecture 20: Linear System Modeling a Single Degree of Freedom Oscillator (experimental determination of damping ratio)
- MIT 18.03 lecture notes 12 "Good vibrations, damping conditions"; 38 "Linearization near equilibria; the nonlinear pendulum"
- JCGM 100:2008 (GUM): 4 Evaluating standard uncertainty (4.2 Type A evaluation of standard uncertainty; 4.3 Type B evaluation of standard uncertainty); 5 Determining combined standard uncertainty (5.1 Uncorrelated input quantities); 6 Determining expanded uncertainty (6.3 Choosing a coverage factor)
  - https://www.bipm.org/documents/20126/2071204/JCGM_100_2008_E.pdf (loaded)
- OpenStax Introductory Statistics 2e Ch. 12 "Linear Regression and Correlation": 12.3 The Regression Equation (subsections: Residuals Plots; Least Squares Criteria for Best Fit; The Correlation Coefficient r); 12.6 Outliers
  - https://openstax.org/books/introductory-statistics-2e/pages/12-3-the-regression-equation (loaded); .../pages/12-6-outliers (loaded)
- SciPy User Guide: "Optimization (scipy.optimize)" > "Least-squares minimization (least_squares)"
  - https://docs.scipy.org/doc/scipy/tutorial/optimize.html (loaded). The `curve_fit` reference page was not loaded.
- MIT 6.003 Lecture 21: Sampling; Lecture 22: Sampling and Quantization (phone frame rates and IMU sample rates; aliasing of the swing)
- Logging tools (supplementary-free): phyphox (RWTH Aachen University), "Download for free" on the home page, with a FAQ section "Export formats"; Tracker video analysis and modeling tool (Open Source Physics), GPL v3
  - https://phyphox.org/ (loaded); https://phyphox.org/source/ (loaded; licence statement); https://phyphox.org/faq/ (loaded); https://opensourcephysics.github.io/tracker-website/ (loaded; physlets.org/tracker redirects here)
- Extract: a physical pendulum model (16.07 L24, OpenStax 15.4); the decay envelope and the logarithmic decrement; a nonlinear least-squares fit with covariance-based uncertainty and GUM reporting; logging the phone's gyroscope (phyphox) or tracking a video (Tracker).
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
