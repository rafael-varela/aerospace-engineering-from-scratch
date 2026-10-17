# Reference map: Phase 05: Statics & Mechanics of Materials

**Note 2026-10-17:** first mapping of Phase 05. Verified 2026-10-17 by loading the live pages of LibreTexts (Baker & Haynes *Engineering Statics: Open and Interactive*; Wierzbicki *Structural Mechanics*), engineeringstatics.org, TU Delft OPEN Books and the DelftX front-matter PDF of Alderliesten, MIT OCW 3.11 (module resource pages and the module PDFs cited below), MIT OCW 16.001 (lecture-notes page, buckling notes, Lab 01 and Lab 02 PDFs), OpenStax University Physics Vol. 1 Ch. 12 and Vol. 2 §1.3, NTRS (Astronautic Structures Manual Vols. I–III, NASA RP-1228), the FreeCAD manual repository on GitHub, and Ford's *Engineering Graphics and Design* on Pressbooks. Eight ids were added to the manifest; `@mit-ocw-3-11-mechanics-of-materials`, `@openstax-university-physics-vol1` and `@openstax-university-physics-vol2` were extended. Not cited: the FreeCAD wiki (wiki.freecad.org returned an Anubis bot-protection "Access Denied" page for TechDraw Workbench, PartDesign Workbench and Manual:Generating 2D drawings, so neither content nor license could be verified); the LibreTexts book shelved as "Introduction to Aerospace Structures and Materials (Alderliesten)" (its units — live loads, arches and cables, influence lines, moment distribution — do not match the TU Delft table of contents and appear to be a different civil-engineering text under the wrong name); Wikipedia "Johnson's parabolic formula" (the fetched text lost its equations and showed no license footer); MatWeb (not openly licensed).

Verified 2026-10-17 by loading the live pages; section numbers and titles are copied from the publishers' own tables of contents, section pages or PDF outlines. Every source here is **required-free** unless marked otherwise; sources are identified by an `@`-prefixed identifier and recorded in [sources.yaml](sources.yaml) with license, version, and caveats. **GAP** marks a lesson that no listed free source covers adequately: the lesson must carry original instruction there, and any "candidate" named is unverified until it is added to the manifest.

Conventions: Baker & Haynes section URLs follow `https://eng.libretexts.org/Bookshelves/Mechanical_Engineering/Engineering_Statics%3A_Open_and_Interactive_(Baker_and_Haynes)/<NN>%3A_<Chapter_Title>/<N>.<MM>%3A_<Section_Title>`; a URL written as `.../<NN>%3A_...` continues that prefix. Section titles are given as LibreTexts shows them; the PreTeXt original at engineeringstatics.org words a few differently (e.g. 4.1 "Moment of Force" vs. LibreTexts "Direction of a Moment"; Ch. 8 "Internal Forces" vs. "Internal Loadings"). MIT OCW 3.11 module pages follow `https://ocw.mit.edu/courses/3-11-mechanics-of-materials-fall-1999/resources/mit3_11f99_<code>/`, each with one Roylance PDF. MIT OCW 16.001 resources follow `https://ocw.mit.edu/courses/16-001-unified-engineering-materials-and-structures-fall-2021/resources/mit16_001_f21_<code>_pdf/`. Licenses: LibreTexts, MIT OCW, OpenStax and Ford are **CC BY-NC-SA 4.0**; Alderliesten is CC BY-NC-SA 4.0 on the TU Delft catalogue page but CC BY-NC-ND 4.0 in the book's own front matter, so treat it as link-only; NASA reports are US Government works ("Public Use Permitted"); the FreeCAD manual is CC BY 4.0. The curriculum links to these sources and writes its own explanations; nothing is copied or adapted.


## 05-01 Forces, Moments & Free-Body Diagrams

Sources: `@libretexts-baker-haynes-engineering-statics`, `@openstax-university-physics-vol1`, `@tudelft-alderliesten-aerospace-structures-materials`

- Baker & Haynes Ch. 1 "Introduction to Statics" (engineeringstatics.org outline): 1.1 Newton's Laws of Motion; 1.2 Units; 1.3 Forces; 1.4 Problem Solving
  - https://engineeringstatics.org/Chapter_01.html (loaded)
- Ch. 2 "Forces and Other Vectors": 2.1 Vectors; 2.3 Two Dimensional Coordinate Systems; 2.4 Three Dimensional Coordinate Systems; 2.5 Unit Vectors; 2.6 Vector Addition; 2.7 Dot Products; 2.8 Cross Products
  - https://eng.libretexts.org/Bookshelves/Mechanical_Engineering/Engineering_Statics%3A_Open_and_Interactive_(Baker_and_Haynes)/02%3A_Forces_and_Other_Vectors (loaded)
- Ch. 4 "Moments and Static Equivalence": 4.1 Direction of a Moment; 4.2 Magnitude of a Moment; 4.3 Varignon's Theorem; 4.4 Moments in Three Dimensions; 4.5 Couples; 4.6 Equivalent Transformations; 4.7 Statically Equivalent Systems
  - .../04%3A_Moments_and_Static_Equivalence (loaded)
- Ch. 5 "Rigid Body Equilibrium": 5.2 Free Body Diagrams
  - .../05%3A_Rigid_Body_Equilibrium/5.02%3A_Free_Body_Diagrams (pattern)
- OpenStax Vol. 1 Ch. 12: 12.1 Conditions for Static Equilibrium (review of the Phase 04 FBD habit, now with torque)
  - https://openstax.org/books/university-physics-volume-1/pages/12-1-conditions-for-static-equilibrium (loaded)
- Aerospace context (Alderliesten): Ch. 6 "Aircraft & spacecraft loads": 6.2 Externally loaded airframe; 6.3 Load path; 6.4 Loads and load paths in an airframe
- Extract: force and moment as vectors (r × F), couples, reduction of a load system to a force-couple at a point; drawing an FBD with correct support symbols.
- Partial GAP: no verified free source treats landing-gear ground loads (strut, drag brace, side brace, tyre reactions); the FBD-solver build and its load cases are original.

## 05-02 Equilibrium of Rigid Bodies

Sources: `@libretexts-baker-haynes-engineering-statics`, `@openstax-university-physics-vol1`, `@mit-ocw-16-001-unified-materials-structures`

- Baker & Haynes Ch. 3 "Equilibrium of Particles" (engineeringstatics.org outline): 3.1 Equilibrium; 3.2 Particles; 3.3 1D Particle Equilibrium; 3.4 2D Particle Equilibrium; 3.5 3D Particle Equilibrium
- Ch. 5 "Rigid Body Equilibrium": 5.1 Degree of Freedom; 5.2 Free Body Diagrams; 5.3 Equations of Equilibrium; 5.4 2D Rigid Body Equilibrium; 5.5 3D Rigid Body Equilibrium; 5.6 Stability and Determinacy; 5.7 Equilibrium Examples
  - .../05%3A_Rigid_Body_Equilibrium (loaded); .../5.06%3A_Stability_and_Determinacy (pattern)
- OpenStax Vol. 1 Ch. 12: 12.1 Conditions for Static Equilibrium; 12.2 Examples of Static Equilibrium
  - https://openstax.org/books/university-physics-volume-1/pages/12-2-examples-of-static-equilibrium (loaded)
- MIT 16.001 lecture notes: "Structural Supports and Reactions; Requirements for Structural Equilibrium" (resource mit16_001_f21_lec05_pdf)
  - https://ocw.mit.edu/courses/16-001-unified-engineering-materials-and-structures-fall-2021/pages/lecture-notes/ (loaded)
- Extract: support reactions by type; writing ΣF = 0, ΣM = 0 as a linear system A·x = b (ties to P01.13); counting unknowns vs. equations to classify determinate, indeterminate and improperly constrained bodies.

## 05-03 Trusses: Joints & Sections

Sources: `@libretexts-baker-haynes-engineering-statics`, `@mit-ocw-3-11-mechanics-of-materials`, `@mit-ocw-16-001-unified-materials-structures`

- Baker & Haynes Ch. 6 "Equilibrium of Structures": 6.1 Structures; 6.2 Interactions between members; 6.3 Trusses; 6.4 Method of Joints; 6.5 Method of Sections
  - .../06%3A_Equilibrium_of_Structures (loaded); .../6.04%3A_Method_of_Joints (pattern)
- Roylance, MIT 3.11 module "Trusses" (PDF dated June 8, 2000): Introduction; Statics analysis of forces; Deflections; Geometrical approach; Energy approach; Statically indeterminate trusses; Matrix analysis of trusses; Stiffness matrix for a single truss element; Assembly of multiple element contributions
  - https://ocw.mit.edu/courses/3-11-mechanics-of-materials-fall-1999/resources/mit3_11f99_truss/ (loaded; PDF loaded)
- MIT 16.001 lecture notes: "Truss Analysis" (mit16_001_f21_lec06lec07_pdf); "Statically Indeterminate Systems" (mit16_001_f21_lec08_pdf); Lab 01 §3 "Truss Lab" (Howe truss, measured vs. predicted member loads)
- Extract: two-force members; joint equilibrium assembled into one matrix (the from-scratch solver); method of sections for a single member; Roylance's matrix section previews the direct stiffness method.

## 05-04 Frames, Machines & Friction

Sources: `@libretexts-baker-haynes-engineering-statics`

- Baker & Haynes Ch. 6 "Equilibrium of Structures": 6.2 Interactions between members; 6.6 Frames and Machines
  - .../06%3A_Equilibrium_of_Structures/6.06%3A_Frames_and_Machines (pattern)
- Ch. 9 "Friction" (engineeringstatics.org outline): 9.1 Dry Friction (9.1.1 Coulomb Friction; 9.1.2 Friction Angle and Friction Resultant; 9.1.3 Normal Forces; 9.1.4 Coulomb Friction Examples); 9.2 Slipping vs. Tipping; 9.3 Wedges; 9.4 Screw Threads; 9.5 Flexible Belts; 9.6 Journal Bearings; 9.7 Rotating Discs
  - https://engineeringstatics.org/Chapter_09-dry-friction.html (as linked from the loaded outline)
  - LibreTexts packs the whole chapter into one page, "9.1: Friction", with subsections 9.1.1 Dry Friction … 9.1.7 Disc Friction, plus a stray "9.2: Untitled Chapter 5": .../09%3A_Friction/9.01%3A_Friction (loaded)
- Extract: dismembering a frame into multi-force members with equal-and-opposite pin forces; mechanical advantage of a bellcrank/pushrod linkage; Coulomb friction, impending motion, slip vs. tip; screw-thread friction feeds 05-07 (bolt torque).
- Partial GAP: no verified free source analyses an aircraft control linkage (stick, bellcranks, pushrods, cable quadrants); the linkage geometry and loads are original.

## 05-05 Centroids & Moments of Inertia

Sources: `@libretexts-baker-haynes-engineering-statics`, `@mit-ocw-16-001-unified-materials-structures`

- Baker & Haynes Ch. 7 "Centroids and Centers of Gravity": 7.1 Weighted Averages; 7.2 Center of Gravity; 7.3 Center of Mass; 7.4 Centroids; 7.5 Centroids using Composite Parts; 7.7 Centroids using Integration; 7.8 Distributed Loads
  - .../07%3A_Centroids_and_Centers_of_Gravity (loaded)
- Ch. 10 "Moments of Inertia": 10.1 Integral Properties of Shapes; 10.2 Moments of Inertia of Common Shapes; 10.3 Parallel Axis Theorem; 10.4 Moment of Inertia of Composite Shapes; 10.5 Polar Moment of Inertia; 10.6 Radius of Gyration; 10.7 Products of Inertia; 10.8 Mass Moment of Inertia
  - .../10%3A_Moments_of_Inertia (loaded); .../10.03%3A_Parallel_Axis_Theorem (pattern)
- MIT 16.001 lecture notes: "Moment of Intertia [sic] of Beam Cross Section, part 01" and "part 02" (mit16_001_f21_lec29part01_pdf, ..._lec29part02_pdf)
- Extract: composite-area tables (A, ȳ, Aȳ, I_c, Ad²); parallel-axis theorem; products of inertia and principal axes (link to eigenvalues, P01.10); 7.8 turns a distributed load into its resultant for 05-09.

## 05-06 Stress & Strain

Sources: `@openstax-university-physics-vol1`, `@mit-ocw-3-11-mechanics-of-materials`, `@tudelft-alderliesten-aerospace-structures-materials`, `@mit-ocw-16-001-unified-materials-structures`

- OpenStax Vol. 1 Ch. 12: 12.3 Stress, Strain, and Elastic Modulus (Tensile or Compressive Stress, Strain, and Young's Modulus; Bulk Stress, Strain, and Modulus; Shear Stress, Strain, and Modulus); 12.4 Elasticity and Plasticity
  - https://openstax.org/books/university-physics-volume-1/pages/12-3-stress-strain-and-elastic-modulus (loaded); .../12-4-elasticity-and-plasticity (pattern)
- Roylance, MIT 3.11 modules: "Introduction to Elastic Response" (PDF title "Introduction to Elasticity", January 21, 2000: Tensile strength and tensile stress; Stiffness); "Stress-Strain Curves" (mit3_11f99_ss); "Constitutive Relations" (PDF title "Constitutive Equations", October 4, 2000: Isotropic elastic materials; Hydrostatic and distortional components; Anisotropic materials)
  - https://ocw.mit.edu/courses/3-11-mechanics-of-materials-fall-1999/resources/mit3_11f99_elas_1/ (loaded; PDF loaded); .../mit3_11f99_const/ (loaded; PDF loaded); .../mit3_11f99_ss/ (pattern)
- Alderliesten Ch. 1 "Material physics & properties": 1.2 Stress-strain; 1.3 Loading modes; 1.4 Engineering terminology; 1.5 Normal stress; 1.6 Shear stress; 1.8 Stiffness and apparent stiffness; 1.10 Toughness
- MIT 16.001: "Concept and Uses of Stress at a Point"; "Strain; Transformation of Strain Components"; "Constitutive Material Response"; Lab 01 §2 "Material Tensile Testing Lab" with Appendix B "Analyzing Material Tensile Strength Testing Data"
- Extract: engineering vs. true stress/strain; E, ν, G and G = E/2(1+ν); yield (0.2 % offset), ultimate, elongation read off a curve — the analyzer's tasks; Lab 01 Appendix B shows how real test data are reduced.

## 05-07 Axial Loading & Thermal Stress

Sources: `@openstax-university-physics-vol2`, `@nasa-rp-1228-fastener-design-manual`, `@nasa-tm-x-73305-astronautic-structures-manual`, `@mit-ocw-3-11-mechanics-of-materials`, `@tudelft-alderliesten-aerospace-structures-materials`

- OpenStax Vol. 2 1.3 Thermal Expansion: subsections "Linear Thermal Expansion" and "Thermal Stress" (Example 1.4 Calculating Thermal Stress)
  - https://openstax.org/books/university-physics-volume-2/pages/1-3-thermal-expansion (loaded)
- Roylance "Trusses": "Deflections" and "Statically indeterminate trusses" (axial δ = PL/AE and compatibility for redundant bars; the module has no thermal content)
- NASA RP-1228 *Fastener Design Manual* (Barrett, 1990): "General Design Information" → "Fastener Torque" (torque coefficient relation T = KFd and torque-measurement methods); "Fatigue-Resistant Bolts" → "Fatigue (Cyclic) Loading of Bolts" (bolt and clamped-member spring constants, stiffness ratio, load-elongation diagram showing preload, external load and separation); "Design Criteria"
  - https://ntrs.nasa.gov/citations/19900009424 (loaded; PDF loaded)
- NASA TM X-73307 *Astronautic Structures Manual* Vol. III, Section D "THERMAL STRESSES": 3.0 STRENGTH OF MATERIALS SOLUTIONS (3.0.1 Unrestrained Beam -- Thermal Loads Only; 3.0.2 Restrained Beam -- Thermal Loads Only; 3.0.6 Trusses)
  - https://ntrs.nasa.gov/citations/19760071127 (loaded; PDF partly loaded)
- Alderliesten: 2.2 The effect of ambient temperature; 11.2 Mechanically fastened joints
- Extract: δ = PL/AE and compatibility for statically indeterminate bars; σ = −EαΔT for a fully restrained bar; a bolted joint as two springs in series/parallel — preload, load share, separation — with T = KFd linking wrench torque to preload.

## 05-08 Torsion of Shafts

Sources: `@mit-ocw-3-11-mechanics-of-materials`, `@mit-ocw-16-001-unified-materials-structures`, `@tudelft-alderliesten-aerospace-structures-materials`, `@libretexts-baker-haynes-engineering-statics`

- Roylance, MIT 3.11 module "Shear and Torsion" (PDF dated June 23, 2000): Shearing stresses and strains; Statics - Twisting Moments; Torsional Stresses and Displacements; Energy method for rotational displacement; Noncircular sections: the Prandtl membrane analogy
  - https://ocw.mit.edu/courses/3-11-mechanics-of-materials-fall-1999/resources/mit3_11f99_torsion/ (loaded; PDF loaded)
- MIT 16.001 lecture notes: "Torsion" (mit16_001_f21_lec30lec31lec32_pdf)
- Baker & Haynes 10.5 Polar Moment of Inertia (J for solid and hollow sections)
- Alderliesten Ch. 7 "Translating loads to stresses": 7.3 Torsional loading of a fuselage structure (closed thin-walled section context)
- Extract: τ = Tr/J, θ = TL/GJ; power P = Tω for shaft sizing (Roylance Example 4 sizes a shaft from 100 hp at 1800 rpm); solid vs. hollow shaft mass efficiency.
- Partial GAP: Roylance does not treat thin-walled closed tubes (Bredt–Batho); whether Alderliesten 7.3 derives it was not verified, so the thin-wall formula is written originally if the lesson uses it.

## 05-09 Shear & Bending Moment Diagrams

Sources: `@libretexts-baker-haynes-engineering-statics`, `@mit-ocw-3-11-mechanics-of-materials`

- Baker & Haynes Ch. 8 "Internal Loadings": 8.1 Internal Loads; 8.2 Sign Conventions; 8.3 Internal Loads at a Point; 8.4 Shear and Bending Moment Diagrams; 8.5 Section Cut Method; 8.6 Relation Between Loading, Shear and Moment; 8.7 Graphical Method; 8.8 Integration Method; 8.9 Geogebra Interactives
  - .../08%3A_Internal_Loadings (loaded); .../8.06%3A_Relation_Between_Loading_Shear_and_Moment (pattern)
- Baker & Haynes 7.8 Distributed Loads (resultant and its line of action)
- Roylance, MIT 3.11 module "Shear and Bending Moment Diagrams" (resource mit3_11f99_statics)
  - https://ocw.mit.edu/courses/3-11-mechanics-of-materials-fall-1999/resources/mit3_11f99_statics/ (pattern; listed on the loaded Modules page)
- Extract: dV/dx = −w, dM/dx = V; singularity (Macaulay) functions or numerical integration for the generator; check that V and M close to zero at a free end.

## 05-10 Bending Stress in Beams

Sources: `@mit-ocw-3-11-mechanics-of-materials`, `@mit-ocw-16-001-unified-materials-structures`, `@tudelft-alderliesten-aerospace-structures-materials`

- Roylance, MIT 3.11 module "Stresses in Beams" (PDF dated November 21, 2000): Introduction; Normal Stresses; Example 1
  - https://ocw.mit.edu/courses/3-11-mechanics-of-materials-fall-1999/resources/mit3_11f99_bstress/ (loaded; PDF loaded)
- MIT 16.001 lecture notes: "Transverse Loading of Slender Structural Elements; Beam Theory" (mit16_001_f21_lec22lec23lec24_pdf)
- Alderliesten Ch. 7 "Translating loads to stresses": 7.4 Bending of a wing structure; 7.5 Case study: bending of wing spar; Ch. 5 "Aircraft & spacecraft structures": 5.5 Wing structures
- Extract: σ = −My/I from plane sections; neutral axis through the centroid (05-05); section modulus; Alderliesten 7.5 is the direct aerospace model for the wing-spar build.

## 05-11 Shear Stress in Beams

Sources: `@mit-ocw-3-11-mechanics-of-materials`, `@tudelft-alderliesten-aerospace-structures-materials`

- Roylance "Stresses in Beams": Shear stresses; Example 2 (T-beam; notes the web–flange junction of I and T beams); Example 3
- Alderliesten 5.6 Torsion box; 7.3 Torsional loading of a fuselage structure (shear in thin-walled skins)
- Extract: τ = VQ/(Ib); parabolic distribution in a rectangle with τ_max = 3V/2A; why the web carries almost all the shear in an I-section.
- Partial GAP: Roylance derives VQ/Ib but never uses the term "shear flow" and does not work an I-beam distribution; no verified free source treats q = VQ/I in thin-walled open sections. Shear flow in the I-beam (flanges and web, sum to V) is original.

## 05-12 Beam Deflection

Sources: `@mit-ocw-3-11-mechanics-of-materials`, `@libretexts-wierzbicki-structural-mechanics`, `@mit-ocw-16-001-unified-materials-structures`

- Roylance, MIT 3.11 module "Beam Displacements" (PDF dated November 30, 2000): Introduction; Multiple integration; Energy method; Superposition
  - https://ocw.mit.edu/courses/3-11-mechanics-of-materials-fall-1999/resources/mit3_11f99_bdisp/ (loaded; PDF loaded)
- Wierzbicki, *Structural Mechanics*: Ch. 4 "Solution Method for Beam Deflections"
  - https://eng.libretexts.org/Bookshelves/Mechanical_Engineering/Structural_Mechanics_(Wierzbicki)/04%3A_Solution_Method_for_Beam_Deflections (listed on the loaded book page)
- MIT 16.001 "Transverse Loading of Slender Structural Elements; Beam Theory"
- Extract: EI v'' = M integrated twice with boundary conditions; superposition from a table of standard cases (Roylance's Fig. 6 table of cantilever and simply supported cases); P02.02 numerical ODE solution as a cross-check.

## 05-13 Stress Transformation & Mohr's Circle

Sources: `@mit-ocw-3-11-mechanics-of-materials`, `@mit-ocw-16-001-unified-materials-structures`, `@tudelft-alderliesten-aerospace-structures-materials`, `@nasa-tm-x-73305-astronautic-structures-manual`

- Roylance, MIT 3.11 module "Tensor Transformations" (PDF title "Transformation of Stresses and Strains", May 14, 2001): Introduction; Direct approach; Mohr's circle; General approach; Principal stresses and planes in three dimensions
  - https://ocw.mit.edu/courses/3-11-mechanics-of-materials-fall-1999/resources/mit3_11f99_trans/ (loaded; PDF loaded)
- MIT 16.001 lecture notes: "Stress Equilibrium"; "Transformation of Stress Components" (mit16_001_f21_lec11lec12_pdf)
- Alderliesten 1.7 Bi-axial loading
- Astronautic Structures Manual Vol. I, Section A3 "COMBINED STRESSES": 3.1.0 Combined Stresses
  - https://ntrs.nasa.gov/citations/19760071125 (loaded; PDF partly loaded)
- Extract: transformation equations as a rotation of a symmetric tensor (P02.14); principal values as eigenvalues; Mohr's circle construction and reading τ_max. The 3.11 "Mohr's circle (Java applet)" and "strs3d" tools listed on the Modules page are Java applets that no longer run in browsers — the lesson's own tool replaces them.

## 05-14 Failure Theories

Sources: `@mit-ocw-3-11-mechanics-of-materials`, `@nasa-tm-x-73305-astronautic-structures-manual`, `@tudelft-alderliesten-aerospace-structures-materials`

- Roylance, MIT 3.11 module "Yield and Plastic Flow" (PDF dated October 15, 2001): Introduction; Multiaxial stress states (Tresca maximum-shear and von Mises distortion-energy criteria); Effect of hydrostatic pressure
  - https://ocw.mit.edu/courses/3-11-mechanics-of-materials-fall-1999/resources/mit3_11f99_yield/ (loaded; PDF loaded)
- Astronautic Structures Manual Vol. I, Section A1 "STRESS AND STRAIN": 1.4.0 Theories of Failure (1.4.1 Elastic Failure; 1.4.2 Interaction Curves); Section A3: 3.2.0 Stress Ratios, Interaction Curves, and Factor of Safety
- Alderliesten Ch. 8 "Considering strength & stiffness": 8.2 Structural performance; 8.3 Selecting the appropriate criterion; Ch. 9 "Design & certification": 9.2 Safety, regulations and specifications
- Extract: σ_vM and τ_max from principal stresses; the two yield envelopes in σ1–σ2 space; margin of safety MS = allowable/(FS·applied) − 1 as aerospace practice reports it (ASM A3.2.0).

## 05-15 Column Buckling

Sources: `@mit-ocw-16-001-unified-materials-structures`, `@libretexts-wierzbicki-structural-mechanics`, `@mit-ocw-3-11-mechanics-of-materials`, `@nasa-tm-x-73305-astronautic-structures-manual`

- MIT 16.001 "Buckling of Beams" (lecture notes 33–35): Introduction to bifurcation of equilibrium and structural instability; A single degree-of-freedom example; Columns; Column buckling (Euler buckling); Example: Simply-supported column; Other boundary conditions; Effects of initial imperfections; Failure of columns
  - https://ocw.mit.edu/courses/16-001-unified-engineering-materials-and-structures-fall-2021/resources/mit16_001_f21_lec33lec34lec35_pdf/ (loaded; PDF loaded)
- Wierzbicki, *Structural Mechanics* Ch. 8 "Stability of Elastic Structures": 8.3 Stability of Elastic Column Using the Energy Method; 8.4 Effect of Structural Imperfections; 8.6 Plastic Buckling of Columns (tangent modulus E_t and critical slenderness β_c = π√(E/σ_y)); Ch. 9 "Advanced Topic in Column Buckling": 9.1 The Tallest Column
  - https://eng.libretexts.org/Bookshelves/Mechanical_Engineering/Structural_Mechanics_(Wierzbicki)/08%3A_Stability_of_Elastic_Structures (loaded); .../8.06%3A_Plastic_Buckling_of_Columns (loaded)
- Roylance "Stresses in Beams": section "Buckling"
- Astronautic Structures Manual: Section C1 "COLUMNS" (listed in the Vol. I subject index; its body was not reached through the verification proxy, so subsection titles are not recorded); A3.3.0 Interaction for Beam-Columns
- Extract: P_cr = π²EI/(KL)² with effective-length factors; slenderness ratio and the Euler/yield crossover; imperfection sensitivity.
- Partial GAP: no verified free source states the Johnson parabolic formula σ_cr = σ_y − (σ_y/2π)²(KL/r)²/E or its tangent point at KL/r = π√(2E/σ_y); that part is original (Wierzbicki 8.6's tangent-modulus discussion is the closest verified material). Candidate: ASM Section C1, once its pages are verified.

## 05-16 Energy Methods

Sources: `@mit-ocw-3-11-mechanics-of-materials`, `@libretexts-wierzbicki-structural-mechanics`

- Roylance "Trusses": Energy approach (strain energy and Castigliano's theorem, Example 3)
- Roylance "Beam Displacements": Energy method (Examples 3 and 4: three-point bend and tapered cantilever)
- Roylance "Shear and Torsion": Energy method for rotational displacement
- Wierzbicki, *Structural Mechanics* Ch. 7 "Energy Methods in Elasticity": 7.1 The Concept of Potential Energy; 7.2 Equivalence of the Minimum Potential Energy and Principle of Virtual Work; 7.6 Castigliano Theorem
  - https://eng.libretexts.org/Bookshelves/Mechanical_Engineering/Structural_Mechanics_(Wierzbicki)/07%3A_Energy_Methods_in_Elasticity (loaded)
- Extract: U for axial, bending and torsion members; δ = ∂U/∂P with a dummy load at the point of interest; checking the 05-12 deflections and a 05-03 truss joint deflection by Castigliano.

## 05-17 CAD & Technical Drawings I

Sources: `@freecad-manual-van-havre`, `@ford-engineering-graphics-and-design`

- FreeCAD manual (van Havre), "Discovering FreeCAD": The FreeCAD interface; Navigating in the 3D view; Parametric objects. "Working with FreeCAD": Modeling for product design (Part Design workbench, Sketcher constraints, Pad, Pocket, Linear Pattern); Generating 2D drawings (TechDraw page from a template, inserting views, dimensions and title block, export to SVG/DXF/PDF); Using spreadsheets
  - https://github.com/yorikvanhavre/FreeCAD-manual/blob/master/SUMMARY.md (loaded); .../working_with_freecad/generating_2d_drawings.md (loaded); .../working_with_freecad/modeling_for_product_design.md (loaded)
- Ford, *Engineering Graphics and Design*: 1. Orthographic Projection; 2. Multiview Layout and Orientation; 4. Section Views; 7. Dimensioning Standards; 8. Conventional Tolerancing; 9. Standard Fits (10. Intro to GD&T and Datum Frames is optional)
  - https://uw.pressbooks.pub/enggraphics/ (loaded)
- Extract: a fully constrained sketch, design intent and parameters (spreadsheet-driven dimensions); third-angle multiview layout; dimensioning rules (no redundant dimensions, dimension features from datums); a title block with material, units and tolerances.
- Partial GAP: the FreeCAD manual's first edition targets FreeCAD 0.16 and its 2D-drawing chapter adds dimensions with the third-party Drawing Dimensioning add-on rather than TechDraw's own dimension tools; the FreeCAD wiki pages that document current (1.x) TechDraw and PartDesign could not be loaded (bot protection). Step-by-step instructions for the current FreeCAD release are original; the official wiki is the candidate to add once verifiable.

## 05-18 Home Lab: Cantilever Beam Deflection

Sources: `@mit-ocw-3-11-mechanics-of-materials`, `@mit-ocw-16-001-unified-materials-structures`, `@libretexts-baker-haynes-engineering-statics`

- Theory: Roylance "Beam Displacements" — Multiple integration (cantilever tip deflection PL³/3EI); Superposition; Roylance "Stresses in Beams" — Normal Stresses (check that the beam stays elastic)
- Section properties: Baker & Haynes 10.2 Moments of Inertia of Common Shapes (I = bh³/12 for a rectangular strip)
- Model lab: MIT 16.001 Lab 02 "Unified Materials and Structures Lab 2" §1 "Deflections of a cantilevered beam subject to concentrated loads" (tip deflection under hanging weights compared with beam theory) and §2 "Measurement of buckling loads in a slender beam subject to uniaxial compression"
  - https://ocw.mit.edu/courses/16-001-unified-engineering-materials-and-structures-fall-2021/resources/mit16_001s_f21_lab02_pdf/ (loaded; PDF loaded)
- E values: Roylance appendix "Material Properties" (mit3_11f99_props, pattern) for aluminium, steel and wood handbook values to compare against
- Extract: plotting δ vs. P and fitting the slope to back out EI (then E); sources of error — clamp compliance, measured vs. nominal thickness (h³ sensitivity), large-deflection limits; uncertainty propagation from P00.15.
- The home procedure (household clamp, ruler or 3D-printed beam, coins or water bottles as weights, phone-camera deflection reading, safety notes) is original; the MIT lab uses lab equipment and is cited as the model, not as the procedure.


## Optional textbooks for deeper study

Paid; never required. Not recorded in `sources.yaml`.

- R. C. Hibbeler, *Engineering Mechanics: Statics* and *Mechanics of Materials* (Pearson) — the standard US problem sets for 05-01 … 05-16.
- J. M. Gere and B. J. Goodno, *Mechanics of Materials* (Cengage) — thorough treatment of thermal stress, shear flow and inelastic columns (05-07, 05-11, 05-15).
- T. H. G. Megson, *Aircraft Structures for Engineering Students* (Butterworth-Heinemann) — thin-walled torsion, open-section shear flow and aircraft-specific buckling; the natural next text after this phase.
- E. F. Bruhn, *Analysis and Design of Flight Vehicle Structures* (Jacobs Publishing) — industry handbook for columns (including Johnson), joints and margins of safety.
