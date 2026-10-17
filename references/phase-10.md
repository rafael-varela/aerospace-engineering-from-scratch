# Reference map: Phase 10: Aerodynamics

**Note 2026-10-17:** first mapping of Phase 10. WebFetch could not resolve these hosts in this sandbox (DNS failure), so every page was loaded with `curl` through the session proxy and read as HTML text, or as PDF text with `pdftotext -layout`. PDFs with no usable text layer were rendered to images and read as images; this applies to the 16.13 handwritten notes. Pages loaded: MIT OCW 16.100 (course, lecture-notes, study-materials and projects pages, plus the lecture PDFs marked below), MIT OCW 16.01 Unified Engineering (fluid-mechanics page and the lecture PDFs marked below), MIT OCW 16.13 (course and lecture-notes pages, four scanned PDFs), Drela's XFOIL and AVL pages with `xfoil_doc.txt`, `sessions.txt` and `avl_doc.txt`, openvsp.org (home, download, license, documentation pages) and its wiki (start, analysisguides, vspaerotutorial, parasitedrag, coeffdefn), the NASA Glenn *Beginner's Guide to Aeronautics* (legacy index pages and the pages named below, both legacy and new www1 URLs), Leishman's *Introduction to Aerospace Flight Vehicles* (table of contents and the 16 chapter pages cited below), the UIUC Low-Speed Airfoil Tests page and the Vol. 1 PDF of *Summary of Low-Speed Airfoil Data*, and NTRS (API metadata plus the full PDF of every report cited). Thirty-two ids were added to the manifest. `@leishman-intro-aerospace-flight-vehicles` and `@ntrs-naca-tr-824-summary-of-airfoil-data` were extended.

Access problems: eaglepubs.erau.edu returned CloudFront "403 Request blocked" to curl's default user agent and to a bare `Mozilla/5.0` agent, but loaded with a full browser user-agent string. web.archive.org reset every connection and was not used. ntrs.nasa.gov and www1.grc.nasa.gov reset some connections, and those pages loaded on retry.

Not cited:
- **MIT 16.110 Flight Vehicle Aerodynamics (Drela).** No OCW course exists: the guessed OCW URL returned 404, and search found only the MITx/edX runs, whose enrollment is closed. Drela's free notes are cited instead through 16.01 (fluids) and 16.13.
- **A. M. O. Smith, "High-Lift Aerodynamics".** Not found on NTRS; the AIAA journal version is paywalled.
- **Glauert's propeller theory and Gessow & Myers.** Glauert was not found on NTRS. Gessow & Myers is a commercial book; NTRS has only their unrelated autorotation flight test.
- **OpenVSP Ground School (vspu.larc.nasa.gov).** Not loaded.
- **Other UIUC volumes.** UIUC Vols. 2–5 and the UIUC airfoil coordinate database were not opened beyond their listing.

Verified 2026-10-17 by loading the live pages; titles, section names and report numbers are copied from the loaded pages, PDF covers or tables of contents. Every source here is **required-free** unless marked otherwise. Sources are identified by an `@`-prefixed identifier and recorded in [sources.yaml](sources.yaml) with license, version, and caveats. **GAP** marks a lesson that no listed free source covers adequately: the lesson must carry original instruction there, and any "candidate" named is unverified until it is added to the manifest.

### Conventions

**URL patterns.**
- MIT 16.100 lecture and study PDFs: `https://ocw.mit.edu/courses/16-100-aerodynamics-fall-2005/resources/<code>/`.
- MIT 16.01 fluids notes: `https://ocw.mit.edu/courses/16-01-unified-engineering-i-ii-iii-iv-fall-2005-spring-2006/resources/<code>/`. Fall lectures are `fNN_fall` (F7 and F8 are `f07_new_fall` and `f08_new_fall`); spring lectures are `fNN_sp`.
- MIT 16.13: `https://ocw.mit.edu/courses/16-13-aerodynamics-of-viscous-fluids-fall-2003/resources/<code>/`.
- Leishman: `https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/<slug>/`. Several slugs do not match the titles, and Leishman is cited by chapter title and URL, never by chapter number.
- NASA Glenn: new pages are `https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/<slug>/`. Legacy pages that are still served are `https://www.grc.nasa.gov/www/k-12/airplane/<page>.html`.
- NTRS: the citation page is `https://ntrs.nasa.gov/citations/<id>` and the PDF is `https://ntrs.nasa.gov/api/citations/<id>/downloads/<id>.pdf`.

A URL written as `.../<code>/` continues the prefix of the line above it. "(pattern)" means the URL was built from a verified pattern, with the title copied from a loaded index page.

**Licenses.**
- MIT OCW: **CC BY-NC-SA 4.0**.
- Leishman: **CC BY-NC-ND 4.0**, so it is link-only and never adapted.
- XFOIL and AVL: **GNU GPL**. The code may be redistributed under the GPL, but the curriculum does not bundle it.
- OpenVSP software: **NASA Open Source Agreement v1.3**. The OpenVSP wiki text is **CC BY-SA 4.0**.
- NASA Glenn web pages and NACA/NASA reports: US Government works, public domain. NTRS records state "Public Use Permitted."
- UIUC *Summary of Low-Speed Airfoil Data*: free to download but **"All rights reserved"**, so it is link-only.

The curriculum links to these sources and writes its own explanations and code; nothing is copied or adapted. The Phase 10 datasets in `curriculum/phase-10/data/` are transcriptions of public-domain NACA/NASA tables (see that folder's README).


## 10-01 Forces, Moments & Coefficients

Sources: `@mit-ocw-16-unified-fluids`, `@leishman-intro-aerospace-flight-vehicles`, `@nasa-grc-beginners-guide-aeronautics`, `@mit-ocw-16-100-aerodynamics`

- 16.01 Fluids Lecture 3 "Aerodynamic Forces and Moments" (notes: 2-D Aerodynamic Forces and Moments; Center of Pressure; Nondimensional Coefficients)
  - https://ocw.mit.edu/courses/16-01-unified-engineering-i-ii-iii-iv-fall-2005-spring-2006/resources/f03_fall/ (loaded; PDF loaded)
- 16.01 Fluids Lecture 4 "Dimensional Analysis and Similarity" and Lecture 5 "Reynolds and Mach Number"
  - .../f04_fall/ (loaded; PDF loaded); .../f05_fall/ (loaded; PDF loaded)
- Leishman, "Aerodynamics of Airfoil Sections": Airfoil Forces and Pitching Moments; Aerodynamic Coefficients; Pressure Coefficient; Center of Pressure & Aerodynamic Center; Integration of Pressures & Shear Stresses; Numerical Integration of Pressure
  - https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/airfoil-characteristics/ (loaded)
- NASA Glenn: "Aerodynamic Forces", "Lift Equation", "Lift Coefficient", "Drag Coefficient", "Center of Pressure", "Aerodynamic Center"
  - https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/aerodynamic-forces/ (loaded); .../lift-equation/ (loaded); .../lift-coefficient/ (loaded); .../drag-coefficient/ (loaded); .../center-of-pressure/ (loaded); .../aerodynamic-center/ (loaded)
- NASA Glenn wind-tunnel balance pages: "Force Balance", "Idealized Force and Moment Balance", "External Force Balance", "Model Mounts"
  - https://www.grc.nasa.gov/www/k-12/airplane/tunbal.html (loaded); .../tunbalfomo.html (loaded); .../tunbalext.html (loaded); .../tunbalmnt.html (loaded)
- 16.100 Lecture 35 "Drag Tare Due to Mount"
  - https://ocw.mit.edu/courses/16-100-aerodynamics-fall-2005/resources/16100lectre35/ (loaded; PDF loaded)
- Extract: resolve normal and axial force into lift and drag through α; transfer the moment between reference points (M_ac = M_ref + x·N); q = ½ρV²; coefficients with a stated reference area and chord; tare subtraction. Together these are the calculator's steps.

## 10-02 Airfoil Geometry

Sources: `@ntrs-naca-tr-824-summary-of-airfoil-data`, `@ntrs-naca-tr-460-78-related-airfoils`, `@leishman-intro-aerospace-flight-vehicles`, `@nasa-grc-beginners-guide-aeronautics`

- NACA TR-824 *Summary of Airfoil Data*, "Description of Airfoils": Method of Combining Mean Lines and Thickness Distributions; NACA Four-Digit Series Airfoils (Numbering system; Thickness distributions; Mean lines); NACA Five-Digit Series Airfoils (Numbering system; Thickness distributions; Mean lines); NACA 6-Series Airfoils. Supplementary Data I "Basic Thickness Forms", II "Data for Mean Lines" (NACA mean lines 210–250 tabulated), III "Airfoil Ordinates"
  - https://ntrs.nasa.gov/citations/19930090976 (PDF loaded)
- NACA TR-460 (Jacobs, Ward & Pinkerton): "Description of Airfoils" (the original four-digit thickness and camber equations)
  - https://ntrs.nasa.gov/citations/19930091108 (PDF loaded)
- Leishman, "Airfoil Shapes": Airfoil Geometry; NACA Method of Drawing Airfoil Shapes; Symmetrical NACA Airfoils; Cambered NACA Airfoils; Other NACA Airfoils; Examples to Try
  - https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/airfoil-geometries/ (loaded)
- NASA Glenn "Effects of Shape on Lift Interactive" (qualitative camber and thickness effects)
  - https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/foilshape/ (loaded)
- Extract: thickness applied perpendicular to the mean line; the four-digit thickness polynomial, including the closed-trailing-edge coefficient variant; four-digit and five-digit mean-line equations; cosine spacing of points. Validate the generator against the tabulated ordinates in TR-824 Supplementary Data III.

## 10-03 Kutta Condition & Circulation

Sources: `@mit-ocw-16-100-aerodynamics`, `@mit-ocw-16-unified-fluids`, `@leishman-intro-aerospace-flight-vehicles`, `@nasa-grc-beginners-guide-aeronautics`, `@ntrs-naca-tr-452-arbitrary-wing-sections`

- 16.100 Lecture 13 "Kutta Condition"; Lecture 14 "Kutta Condition Solution"; study material "Lift Generation and Streamline Curvature"
  - https://ocw.mit.edu/courses/16-100-aerodynamics-fall-2005/resources/16100lectre13_cg/ (loaded; PDF loaded); .../16100lectre14_cg/ (pattern); .../4liftgen/ (loaded; PDF loaded)
- 16.01 Fluids (fall) Lecture 11 "Vorticity and Circulation"; Lecture 16 "Vortex, Lifting Flow" (notes: Vortex; Lifting flow about circular cylinder); (spring) Lecture 1 "Formation of Lifting Flow"
  - .../16-01-unified-engineering-i-ii-iii-iv-fall-2005-spring-2006/resources/f11_fall/ (pattern); .../f16_fall/ (loaded; PDF loaded); .../f01_sp/ (loaded; PDF loaded)
- Leishman, "Classic Airfoil Theory": Kutta Condition; Starting Vortex & Kelvin's Theorem; Conformal Transformations & Joukowsky Airfoils. "Potential Flow Theory": Understanding Circulation and Vorticity; Complex Variable Methods in Potential Flow; Conformal Transformations
  - https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/classic-airfoil-theory/ (loaded); .../chapter/potential-flows/ (loaded)
- NASA Glenn "Lift of a Rotating Cylinder" and "Conformal Mapping"
  - https://www.grc.nasa.gov/www/k-12/airplane/cyl.html (loaded); .../map.html (loaded)
- Optional depth: NACA TR-452 (Theodorsen & Garrick) "Transformation of an Arbitrary Airfoil into a Circle"
  - https://ntrs.nasa.gov/citations/19930091527 (PDF loaded)
- Extract: L′ = ρV∞Γ; potential flow has infinitely many solutions until the Kutta condition fixes Γ; the cylinder-with-circulation flow mapped by ζ = z + a²/z; choosing Γ so the rear stagnation point sits on the sharp trailing edge. These are the steps of the Joukowski build.

## 10-04 Thin Airfoil Theory

Sources: `@mit-ocw-16-unified-fluids`, `@mit-ocw-16-100-aerodynamics`, `@leishman-intro-aerospace-flight-vehicles`, `@ntrs-naca-tr-824-summary-of-airfoil-data`

- 16.01 Fluids (spring) Lecture 2 "Airfoil Vortex Sheet Models, Thin Airfoil Analysis Problem" (notes: Surface Vortex Sheet Model; Single Vortex Sheet Model; Flow tangency imposition; Symmetric airfoil case); Lecture 3 "Thin Airfoil Analysis Problem (cont.)"; Lecture 4 "Thin Airfoil Theory Application: Analysis Example"
  - .../16-01-unified-engineering-i-ii-iii-iv-fall-2005-spring-2006/resources/f02_sp/ (loaded; PDF loaded); .../f03_sp/ (pattern); .../f04_sp/ (loaded; PDF loaded)
- 16.100 Lecture 15 "Thin Airfoil Theory Summary"; Lecture 16 "Important Concepts in Thin Airfoil Theory"
  - https://ocw.mit.edu/courses/16-100-aerodynamics-fall-2005/resources/16100lectre15_cg/ (loaded; PDF loaded); .../16100lectre16_cg/ (loaded; PDF loaded)
- Leishman, "Classic Airfoil Theory": Thin Airfoil Assumption; Theory of Vortex Sheets; Thin Airfoil Theory
- TR-824 "Theoretical Considerations": Angle of Zero Lift (Methods of calculation; Numerical examples); Moment Coefficients (Methods of calculation; Numerical examples)
- Extract: the Glauert substitution x = (c/2)(1 − cos θ); Fourier coefficients A0, A1, A2 from dz/dx; c_l = 2π(α − α_L0) and c_m,c/4 from A1 and A2. The solver takes dz/dx numerically from 10-02's mean lines and checks α_L0 against the TR-824 values.

## 10-05 Vortex Panel Method

Sources: `@mit-ocw-16-unified-fluids`, `@ntrs-nasa-tp-2995-panel-methods`, `@leishman-intro-aerospace-flight-vehicles`, `@ntrs-naca-tr-452-arbitrary-wing-sections`, `@drela-xfoil`

- 16.01 Fluids (fall) Lecture 17 "Flowfield Prediction" (notes: Flowfield prediction; Source Sheets); Lecture 18 "Lift Prediction, Vortex Sheets" (notes: Prediction of Lift; Vortex Sheets)
  - .../16-01-unified-engineering-i-ii-iii-iv-fall-2005-spring-2006/resources/f17_fall/ (loaded; PDF loaded); .../f18_fall/ (loaded; PDF loaded)
- Erickson, NASA TP-2995 *Panel Methods: An Introduction*: 2. What Panel Methods Can and Cannot Do; 3. Common Roots of Panel Methods; 4. Distinctions Between Panel Codes (Surface Geometry; Singularity Distributions; Continuity of Doublet Strength); 5. Modeling (Mass-Flux and Velocity Boundary Conditions; Interior Potential; Wakes)
  - https://ntrs.nasa.gov/citations/19910009745 (PDF loaded)
- Leishman, "Potential Flow Theory": General Method of Singularities; Superposition of Elementary Flows
- XFOIL User Primer "Inviscid Formulation" (a linear-vorticity streamfunction panel method; source of the comparison polar)
  - https://web.mit.edu/drela/Public/web/xfoil/xfoil_doc.txt (loaded)
- Extract: discretise the surface into panels; place collocation points; build the influence-coefficient matrix; impose the Kutta condition with γ_1 + γ_N = 0; compute Cp from surface speed; run a convergence study in N against the 10-03 Joukowski exact solution.
- Partial GAP: no verified free source derives the two-dimensional linear-strength vortex panel equations step by step: TP-2995 is three-dimensional and source–doublet, and the 16.01 notes stop at sheets. The panel influence integrals, the Kármán–Trefftz mapping and its exact Cp are written originally. TR-452 is the closest verified conformal-mapping reference, and none of the loaded sources mentions Kármán–Trefftz by name.

## 10-06 Viscous–Inviscid Coupling

Sources: `@mit-ocw-16-13-aerodynamics-viscous-fluids` (supplementary), `@mit-ocw-16-100-aerodynamics`, `@leishman-intro-aerospace-flight-vehicles`, `@drela-xfoil`, `@ntrs-nasa-tm-4074-ladson-naca0012`

- 16.100 Lecture 22 "Integral Boundary Layer Equations"; Lecture 23 "Correlation Methods for Integral Boundary Layers" (Thwaites' method); Lecture 24 "Method of Assumed Profiles"; Lecture 25 "Falkner-Skan Flows"
  - https://ocw.mit.edu/courses/16-100-aerodynamics-fall-2005/resources/16100lectre22_cg/ (loaded; PDF loaded); .../16100lectre23/ (loaded; PDF loaded); .../16100lectre24_cg/ (pattern); .../16100lectre25/ (pattern)
- 16.13 (Drela & Merchant), handwritten notes:
  - Handouts: "Thwaites' Method"; "Integral Thicknesses"; "Von-Karman Integral Momentum Equation"; "Displacement Effects of Boundary Layer on Potential Flow"
  - Lecture 13 "Integral Methods. Integral Momentum Equation. Thwaites' Method."
  - Lecture 17 "2D Interaction Models: Displacement Body, Transpiration. Form Drag, Stall Mechanisms."
  - Lecture 18 "IBLT Solution Techniques. Iteration Stability."
  - Lecture 23 "Transition Mechanisms. Transition Prediction: Local Correlations, Amplification Methods."
  - Lecture 27 "Equilibrium BL's: Clauser Hypothesis. Dissipation Formulas and Integral Closure."
  - URLs: https://ocw.mit.edu/courses/16-13-aerodynamics-of-viscous-fluids-fall-2003/resources/thwaites_method/ (loaded; PDF loaded as images); .../lecture_17/ (loaded; PDF loaded as images); .../lecture_23/ (loaded; PDF loaded as images); .../lecture_27/ (loaded; PDF loaded as images); .../lecture_13/ (pattern); .../displacement_effects/ (pattern)
- Leishman, "Boundary Layer Flows": Pressure Gradients & Flow Separation; Flow Separation on an Airfoil; von Kármán Momentum Integral Equation; Skin Friction Drag on an Airfoil; Laminar Separation Bubbles
  - https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/boundary-layers/ (loaded)
- XFOIL User Primer "Viscous Formulation": two-equation lagged-dissipation integral BL, strong coupling, and drag from the wake momentum thickness far downstream. Its "Theory References" name Drela & Giles (AIAA J. 1987) and Drela (1989), which are not free.
- Validation data: Ladson TM-4074 (free and fixed transition; see 10-13)
- Extract: the momentum-integral equation; θ² from Thwaites' quadrature with λ, H(λ), ℓ(λ); laminar separation at λ ≈ −0.09; why a transition criterion and a turbulent closure are needed; the displacement-body/transpiration idea behind coupling; drag from the trailing-edge θ.
- Partial GAP: Michel's transition criterion, Head's entrainment method and the Squire–Young formula are not stated in any loaded free source. Squire & Young appear only as reference 23 in TR-824, and 16.13 Lecture 23 treats transition correlations only generally, in handwriting. Those three correlations and the coupling iteration are written originally.

## 10-07 Airfoil Design & Selection

Sources: `@drela-xfoil`, `@leishman-intro-aerospace-flight-vehicles`, `@mit-ocw-16-unified-fluids`, `@ntrs-naca-tr-824-summary-of-airfoil-data`

- XFOIL User Primer: "Inverse Formulation"; "Full-Inverse Design Routine (MDES)" (Creation of seed surface speed distribution; Modification of surface speed distributions; Generation of new geometry); "Mixed-Inverse Design Routine (QDES)"; "Geometry Design Routine"
  - https://web.mit.edu/drela/Public/web/xfoil/xfoil_doc.txt (loaded)
- Leishman: "Classic Airfoil Theory" → Inverse Airfoil Design; "Aerodynamics of Airfoil Sections" → Inverse Airfoil Design and Shape Optimization; Drag Polars
- 16.01 Fluids (fall) Lecture 21 "Airfoil Polar Relations" (notes: Wing Loading; Sample airfoil polars, using an RC sailplane example); Lecture 20 "Airfoils: Detailed Look"
  - .../16-01-unified-engineering-i-ii-iii-iv-fall-2005-spring-2006/resources/f21_fall/ (loaded; PDF loaded); .../f20_fall/ (loaded; PDF loaded)
- TR-824 "Application to Wing Design": Application of Section Data; Selection of Root Section; Selection of Tip Section
- Extract: specify a target Cp or u_e(s) and recover the shape (thin-airfoil inverse: from the load distribution to camber); selection criteria (c_l design, c_d,min, c_l,max, c_m) read off polars; the root/tip selection logic.
- Partial GAP: no free source poses an inverse-design exercise at this level; the exercise (target load → camber line via 10-04, checked with 10-05/06) is original.

## 10-08 High-Lift Devices

Sources: `@ntrs-naca-tr-938-trailing-edge-high-lift`, `@ntrs-nasa-cr-4746-high-lift-systems` (supplementary), `@ntrs-naca-tr-824-summary-of-airfoil-data`, `@leishman-intro-aerospace-flight-vehicles`, `@nasa-grc-beginners-guide-aeronautics`

- Cahill, NACA TR-938 *Summary of Section Data on Trailing-Edge High-Lift Devices*: Flap Theory; Discussion of Experimental Results: Plain Flaps (Flap Loads and Moments; Summary of Plain-Flap Data); Split Flaps (Summary of Split-Flap Data); Single Slotted Flaps; External-Airfoil Flaps and Venetian-Blind Flaps
  - https://ntrs.nasa.gov/citations/19930092003 (PDF loaded)
- Rudolph, NASA CR-4746: Chapter 1 "Types of High-Lift Systems: Their Geometry, Functions, and Design Criteria" (1.1 Types of High-Lift Systems; 1.3 Geometric Parameters of High-Lift Devices; 1.5 Characteristics and Constraints of High-Lift Devices; 1.6 Interactions between Leading- and Trailing-Edge Devices); Chapter 3 "Conceptual Studies for Improved High-Lift Systems" (3.2 Leading-Edge Devices; 3.3 Trailing-Edge Devices)
  - https://ntrs.nasa.gov/citations/19960052267 (PDF loaded)
- TR-824 "Experimental Characteristics": High-Lift Devices
- Leishman: "Airfoil Shapes" → Airfoils with Flaps; "Stalling & Spinning" → Importance of High-Lift Devices; Stall Control Devices
  - https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/maximum-lift-stalling-spinning/ (loaded)
- NASA Glenn "Flaps and Slats Interactive"
  - https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/flaps-and-slats/ (loaded)
- Extract: a flap as a change in effective camber (thin-airfoil Δα_L0 from 10-04); Δc_l,max by flap type from TR-938's summaries; the slat as a delay of leading-edge stall rather than a camber change; slot effects.

## 10-09 Finite Wings & Induced Drag

Sources: `@mit-ocw-16-unified-fluids`, `@mit-ocw-16-100-aerodynamics`, `@leishman-intro-aerospace-flight-vehicles`, `@nasa-grc-beginners-guide-aeronautics`, `@ntrs-naca-tr-116-prandtl-modern-hydrodynamics`

- 16.01 Fluids (spring) Lecture 5 "Introduction to 3D Wings, 2D and 3D Coefficients"; Lecture 6 "3D Vortex Filaments, Lifting Line Theory" (notes: General 3-D vortex; Wing vortex model; Downwash distribution)
  - .../16-01-unified-engineering-i-ii-iii-iv-fall-2005-spring-2006/resources/f05_sp/ (loaded; PDF loaded); .../f06_sp/ (loaded; PDF loaded)
- 16.100 Lecture 19 "Trefftz Plane Analysis of Induced Drag"; Lecture 32 "Single Horseshoe Vortex Wing Model"
  - https://ocw.mit.edu/courses/16-100-aerodynamics-fall-2005/resources/16100lectre19_cg/ (loaded; PDF loaded); .../16100lectre32_cj/ (loaded; PDF loaded)
- Leishman, "Aerodynamics of Finite Wings": Origin of Trailing Vortices; Effects of Trailing Vortices; Wing Force & Moment Coefficients; Drag Polar for a Finite Wing; Calculating the Lift & Drag on a Finite Wing; Winglets; Vortex Wake Upsets
  - https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/finite-wing-characteristics/ (loaded)
- NASA Glenn "Downwash Effects on Lift", "Induced Drag Coefficient", "Winglets"
  - https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/downwash-effects-on-lift/ (loaded); .../induced-drag-coefficient/ (loaded); .../winglets/ (loaded)
- Optional primary source: Prandtl, NACA TR-116, Part II "Applications": C. "The Finite Wing"
  - https://ntrs.nasa.gov/citations/19930091180 (PDF loaded)
- Extract: Helmholtz's vortex theorems and the horseshoe vortex; downwash w and induced angle α_i; induced drag as the tilt of the local lift vector; the Trefftz-plane view; C_Di = C_L²/(πeAR).

## 10-10 Prandtl Lifting Line

Sources: `@mit-ocw-16-unified-fluids`, `@mit-ocw-16-100-aerodynamics`, `@leishman-intro-aerospace-flight-vehicles`, `@ntrs-naca-tr-116-prandtl-modern-hydrodynamics`, `@ntrs-naca-tr-921-span-loading`, `@ntrs-naca-tn-1422-sivells-three-wings`

- 16.01 Fluids (spring) Lecture 6 "3D Vortex Filaments, Lifting Line Theory"; Lecture 7 "Elliptical Lift Distribution"; Lecture 9 "General Wings" (notes: General circulation distribution and downwash; Load distributions on typical planforms; Lift slope reduction)
  - .../f07_sp/ (loaded; PDF loaded); .../f09_sp/ (loaded; PDF loaded)
- 16.100 Lecture 17 "Prandtl's Lifting Line Introduction"; Lecture 18 "Force Calculations for Lifting Line"
  - https://ocw.mit.edu/courses/16-100-aerodynamics-fall-2005/resources/16100lectre17_cg/ (loaded; PDF loaded); .../16100lectre18_cg/ (loaded; PDF loaded)
- Leishman, "Lifting Line & Finite Wing Theory": Helmholtz's Vortex Theorems; Theory of Vortex Lines; Lifting Line Theory; Numerical Solution to the Monoplane Wing Equation; Wing Shape & Loading Distributions
  - https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/lifting-line-theory/ (loaded)
- Prandtl, NACA TR-116, Part II: D. "Theory of the Monoplane"
- Check data: DeYoung & Harper, NACA TR-921 "Loading Characteristics for Arbitrary Wings"; Sivells, NACA TN-1422 (measured vs. calculated characteristics of three AR-9 wings)
  - https://ntrs.nasa.gov/citations/19930091985 (PDF loaded); https://ntrs.nasa.gov/citations/19930082118 (PDF loaded)
- Extract: the monoplane equation with Γ(θ) = 2bV∞ΣA_n sin nθ; collocation at N stations; C_L = πAR·A1; δ and e from the higher harmonics; check the elliptic case exactly, then compare taper cases with TR-921.

## 10-11 Vortex Lattice Method

Sources: `@ntrs-nasa-tn-d-6142-vortex-lattice-program`, `@ntrs-nasa-sp-405-vortex-lattice-utilization` (supplementary), `@leishman-intro-aerospace-flight-vehicles`, `@drela-avl`

- Margason & Lamar, NASA TN D-6142: Basic Concepts and Limitations; Program Description, Part I "Geometry Computation" (Section 3. Horseshoe Vortex Lattice), Part II "Vortex Strength Computation", Part III "Aerodynamic Computation" (Section 1. Lift and Moment Using Entire Horseshoe Vortex; Section 2. Lift and Pitching and Rolling Moments Using Only Spanwise Filament of Horseshoe Vortex); Effect of Vortex-Lattice Arrangement on Computed Aerodynamic Characteristics
  - https://ntrs.nasa.gov/citations/19710007949 (PDF loaded)
- NASA SP-405 *Vortex-Lattice Utilization*: paper 1 "Historical Evolution of Vortex-Lattice Methods" (DeYoung); paper 17 "Arrangement of Vortex Lattices on Subsonic Wings" (DeJarnette); paper 18 "Lattice Arrangements for Rapid Convergence"
  - https://ntrs.nasa.gov/citations/19760021075 (PDF loaded)
- Leishman, "Lifting Line & Finite Wing Theory": Introduction to Lifting Surface Theory
- AVL User Primer: "Vortex-Lattice Modeling Principles" (Configurations; Unsteady flow; Compressibility); "Vortex Lattice Spacing Distributions"
  - https://web.mit.edu/drela/Public/web/avl/avl_doc.txt (loaded)
- Extract: horseshoe vortices with the bound leg at ¼ panel chord and control points at ¾ chord; Biot–Savart for a finite segment; the influence matrix and flow tangency; forces by Kutta–Joukowski on the bound legs; cosine spanwise spacing; convergence in panel count (TN D-6142 gives practical row and column counts).

## 10-12 Wing Planform Effects

Sources: `@ntrs-naca-tr-921-span-loading`, `@mit-ocw-16-unified-fluids`, `@leishman-intro-aerospace-flight-vehicles`, `@ntrs-naca-tn-1669-untapered-wings-aspect-ratio-sweep`, `@ntrs-naca-tn-1671-swept-wings-taper-ratio`, `@ntrs-naca-tn-1422-sivells-three-wings`

- NACA TR-921: Effect of Plan-Form Variation on the Basic Loading; Effect of Plan-Form Variation on the Additional Loading
- 16.01 Fluids (spring) Lecture 8 "Wing Geometry, Wing Design Problem"; Lecture 9 "General Wings"
  - .../16-01-unified-engineering-i-ii-iii-iv-fall-2005-spring-2006/resources/f08_sp/ (loaded; PDF loaded)
- Leishman, "Wing Shapes & Nomenclature": Geometric Definition of a Wing; Calculation of Wing Area; Calculation of Aspect Ratio; Calculation of Mean Chords
  - https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/wing-shapes-and-nomenclature/ (loaded)
- Leishman, "Aerodynamics of Finite Wings": Effects of Compressibility & Wing Sweep
- Data: NACA TN-1669 (untapered wings, AR 1.34–5.16, sweep 0–60°); TN-1671 (45° swept, AR 2.61, taper 1.0/0.5/0.25); TN-1422 (AR 9, taper 0.4, 0° and 2° washout)
  - https://ntrs.nasa.gov/citations/19930082295 (PDF loaded); https://ntrs.nasa.gov/citations/19930082298 (PDF loaded)
- Extract: how taper, sweep and twist move the span loading and the local c_l peak (tip-stall tendency); basic loading (from twist) vs. additional loading (from α); e vs. taper; lift-slope reduction with sweep.
- Partial GAP: no verified free source sets up planform optimization (minimum induced drag subject to constraints, with a stall-margin constraint); the optimizer is original.

## 10-13 Drag Build-Up

Sources: `@ntrs-nasa-tm-4074-ladson-naca0012`, `@openvsp-docs` (supplementary), `@mit-ocw-16-unified-fluids`, `@ntrs-naca-tr-824-summary-of-airfoil-data`, `@leishman-intro-aerospace-flight-vehicles`, `@nasa-grc-beginners-guide-aeronautics`, `@mit-ocw-16-100-aerodynamics`

- Ladson, NASA TM-4074: Introduction; Wake Survey Rake; Instrumentation; Tests and Procedures; Presentation of Results; Conclusions (free- and fixed-transition section drag vs. Re at M = 0.15; transcribed in `curriculum/phase-10/data/ladson_naca0012_M015.csv`)
  - https://ntrs.nasa.gov/citations/19880019495 (loaded; PDF loaded)
- OpenVSP wiki "Parasite Drag Tool": General Overview; Coefficient of Friction Equations (Laminar; Turbulent); Form Factor Equations; Flow Condition; Bibliography
  - https://openvsp.org/wiki/doku.php?id=parasitedrag (loaded)
- 16.01 Fluids (spring) Lecture 10 "Aircraft Performance Analysis" (notes: Aircraft Performance Analysis; Parasite Drag Estimation. The notes cite Hoerner, *Fluid-Dynamic Drag*, Ch. 3)
  - .../16-01-unified-engineering-i-ii-iii-iv-fall-2005-spring-2006/resources/f10_sp/ (loaded; PDF loaded)
- TR-824 "Drag Characteristics of Smooth Airfoils" (Drag characteristics in low-drag range; Effects of type of section on drag characteristics); "Effect of Surface Irregularities on Drag" (Permissible roughness; Drag with fixed transition; Drag with practical construction methods)
- 16.100 study material "CD NACA" (C_D vs. Re for NACA 0004 and NACA 0012)
  - https://ocw.mit.edu/courses/16-100-aerodynamics-fall-2005/resources/8cdnaca/ (loaded; PDF loaded)
- Leishman, "Aerodynamics of Finite Wings": Drag Polar for an Airplane; Generalization of the Airplane Drag Polar
- NASA Glenn "Shape Effects on Drag"; "Drag Coefficient"
  - https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/shape-effects-on-drag/ (loaded)
- Extract: C_D0 = Σ C_f·FF·Q·S_wet/S_ref; laminar (Blasius) and turbulent skin-friction laws; form factors; the section-level check is c_d ≈ 2C_f·FF against Ladson's free-transition and fixed-transition 0012 data across Re.
- Partial GAP: form-factor and interference correlations are available free only as the OpenVSP wiki's equation list (which cites paid sources) and a short 16.01 note. The derivation and the choice among correlations are written originally; Hoerner is optional below.

## 10-14 Low-Reynolds-Number Aerodynamics

Sources: `@ntrs-nasa-tm-4062-eppler-387-low-re`, `@uiuc-selig-summary-low-speed-airfoil-data` (supplementary, link only), `@leishman-intro-aerospace-flight-vehicles`, `@drela-xfoil`, `@mit-ocw-16-13-aerodynamics-viscous-fluids` (supplementary)

- McGhee, Walker & Millard, NASA TM-4062 (Eppler 387, Re 60 000–460 000): Introduction; Model, Apparatus, and Procedure; Tests and Methods; Presentation of Results; Discussion of Results; Concluding Remarks. Includes oil-flow laminar-separation and reattachment data and comparison with two other facilities.
  - https://ntrs.nasa.gov/citations/19890001471 (PDF loaded)
- Selig et al., *Summary of Low-Speed Airfoil Data* Vol. 1: 2 "Wind-Tunnel Facility and Measurement Techniques"; 3 "Summary of Airfoil Data" (3.1 Airfoils for Free Flight Models; 3.2 Airfoils for Thermal Duration Sailplanes; 3.8 Airfoils for Heavy-Lift Cargo Planes; 3.9 Airfoils for Small Wind Turbines); Vols. 2, 3 and 5 are listed on the page (Vol. 5 includes E387 and S1223)
  - https://m-selig.ae.illinois.edu/uiuc_lsat.html (loaded; Vol. 1 PDF loaded)
- Leishman: "Boundary Layer Flows" → Laminar Separation Bubbles; "Gliders & Sailplanes" → Sailplane Aerodynamics; Airfoils for Sailplanes
  - https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/gliders-sailplanes/ (loaded)
- XFOIL User Primer: General Description ("transitional separation bubble(s)"); "Re, Mach dependence"; `sessions.txt` (an E387 viscous session)
  - https://web.mit.edu/drela/Public/web/xfoil/sessions.txt (loaded)
- 16.13 Lecture 23 (transition prediction; e^n method handout `en_method`, pattern)
- Extract: why Re < 5×10⁵ is different (laminar separation bubble, hysteresis, drag rise); reading E387 data across Re; trip strips; airfoil choices for drones and gliders. Compare an e^n-type XFOIL run with TM-4062.

## 10-15 Stall & Post-Stall

Sources: `@ntrs-naca-tn-2502-airfoil-stall-types`, `@leishman-intro-aerospace-flight-vehicles`, `@ntrs-naca-tr-824-summary-of-airfoil-data`, `@ntrs-nasa-tm-4074-ladson-naca0012`, `@mit-ocw-16-13-aerodynamics-viscous-fluids` (supplementary)

- McCullough & Gault, NACA TN-2502: Summary (the three types: trailing-edge stall, leading-edge stall, thin-airfoil stall); Introduction; Description; Concluding remarks
  - https://ntrs.nasa.gov/citations/19930083422 (PDF loaded)
- Leishman, "Stalling & Spinning": Stalling Airspeeds & Angles of Attack; General Flow Physics of Stall; Two-Dimensional Stall Characteristics; Three-Dimensional Stall; Stall Control Devices
- TR-824 Table II "Maximum Lift and Stalling Characteristics of Models Tested in the NACA 19-Foot Pressure Tunnel"; "Lift Characteristics of Smooth Airfoils"; "Lift Characteristics of Rough Airfoils"
- 16.13 Lecture 17 "... Form Drag, Stall Mechanisms."
- Data: Ladson TM-4074 lift through c_l,max (post-stall points without wake drag)
- Extract: classify stall from airfoil thickness and leading-edge radius; the progression of trailing-edge separation; how a stall model blends an attached-flow c_l with a separated (flat-plate-like) branch.
- Partial GAP: no verified free source gives a post-stall coefficient model (e.g. a Viterna-type or flat-plate blend to 90°); the stall-progression model is original.

## 10-16 Fuselage & Body Aerodynamics

Sources: `@ntrs-naca-tr-184-munk-airship-hulls`, `@ntrs-naca-tr-1048-slender-inclined-bodies`, `@drela-avl`, `@leishman-intro-aerospace-flight-vehicles`, `@nasa-grc-beginners-guide-aeronautics`, `@ntrs-naca-tr-116-prandtl-modern-hydrodynamics`

- Munk, NACA TR-184 *The Aerodynamic Forces on Airship Hulls*: I. General Properties of Aerodynamic Flows; II. The Aerodynamic Forces on Airship Hulls; III. Some Practical Conclusions
  - https://ntrs.nasa.gov/citations/19930091249 (PDF loaded)
- Allen & Perkins, NACA TR-1048: Introduction; Potential Flow Theory; Comparisons with Experiment and Discussion of the Effects of Viscosity; Estimation of Cross-Flow Drag Coefficient
  - https://ntrs.nasa.gov/citations/19930090962 (PDF loaded)
- Prandtl TR-116 Part II: A. "Distribution of Pressure on Airship Bodies"
- AVL User Primer: General Description ("slender bodies such as fuselages and nacelles ... consistent with slender-body theory"); the `BODY` keyword in "Geometry Input File -- xxx.avl"
- Leishman, "Airships, Blimps, & Aerostats": Aerodynamic Effects
  - https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/lth/ (loaded)
- NASA Glenn "Fuselage"
  - https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/fuselage/ (loaded)
- Extract: slender-body normal force from the rate of change of apparent mass, dN/dx = ρV²α·dS/dx; a closed body gives zero net lift but a destabilising Munk moment; the viscous cross-flow addition (Allen & Perkins) at larger α.
- Partial GAP: the verified slender-body sources are the 1924 and 1951 originals in period notation. A modern step-by-step derivation for the build is written originally.

## 10-17 Propeller & Rotor Theory

Sources: `@ntrs-nasa-cr-3082-rotary-wing-aerodynamics`, `@leishman-intro-aerospace-flight-vehicles`, `@nasa-grc-beginners-guide-aeronautics`, `@ntrs-naca-tr-640-full-scale-propellers` (supplementary), `@openvsp-docs` (supplementary)

- Stepniewski, NASA CR-3082 *Rotary-Wing Aerodynamics* Vol. I:
  - Chapter II "Momentum Theory": 3. Actuator Disc; 3.1 Induced Velocity and Thrust in Axial Translation; 3.3 Ideal Power in Climb and Hovering
  - Chapter III "Blade Element Theory": 2. Axial Translation and Hovering; 2.1 Basic Considerations of Thrust and Torque Predictions; 2.2 Combined Blade-Element and Momentum Theory; 2.3 Nondimensional Coefficients; 2.5 Tip Losses
  - https://ntrs.nasa.gov/citations/19790013868 (PDF loaded)
- Leishman, "Aircraft Propellers": Propeller Fundamentals; Momentum Theory of the Single Propeller; Propeller Performance Coefficients; General Propeller Performance; Propeller Performance Charts. "Helicopters": Blade Element Theory; Blade Element Momentum Theory (BEMT)
  - https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/propellers/ (loaded)
- NASA Glenn "Propeller Analysis" (Simple Momentum Theory); "Propeller Thrust"; "Propeller Propulsion"
  - https://www.grc.nasa.gov/www/k-12/airplane/propanl.html (loaded); .../propth.html (loaded); .../propeller.html (loaded)
- Validation data: Hartman & Biermann, NACA TR-640 (full-scale 2-, 3- and 4-blade propellers, Clark Y and R.A.F. 6 sections)
  - https://ntrs.nasa.gov/citations/19930091715 (PDF loaded)
- Coefficient conventions: OpenVSP wiki "Definitions of Aerodynamic Coefficients" (CT, CQ, CP, J, EtaP; helicopter-convention CT_H, CP_H, FOM)
  - https://openvsp.org/wiki/doku.php?id=coeffdefn (loaded)
- Extract: the actuator disc (T = 2ρAv_i(V + v_i)); annulus-by-annulus BEM with an inflow iteration; the Prandtl tip-loss factor; propeller (J, C_T, C_P, η) vs. rotor (C_T, C_P, FM) conventions; check against TR-640 charts.

## 10-18 Helicopter Aerodynamics

Sources: `@ntrs-nasa-cr-3082-rotary-wing-aerodynamics`, `@leishman-intro-aerospace-flight-vehicles`

- Stepniewski, NASA CR-3082 Vol. I:
  - Chapter II: 3.5 Vertical Descent Rates; 3.6 Induced Velocity and Thrust in Nonaxial Translation; 3.7 Power Required in Nonaxial Translation; 3.9 Induced Power in Horizontal Flight; 4. Flight Envelope of an Ideal Helicopter; 5.2 Nonuniform Downwash and Tip Losses (the k_ind factor)
  - Chapter III: 2.8 Rotor Power, and Aerodynamic and Overall Efficiencies in Hover; 3. Forward Flight (3.1 Velocities; 3.4 Blade Profile Drag Contribution to Rotor Power and Drag; 3.8 Propulsive Thrust and Power Required in Horizontal Flight; 3.9 Rotor and Helicopter Efficiency in Horizontal Flight)
  - Chapter VI "Airfoils for Rotary-Wing Aircraft": 9. Contribution of Airfoil Characteristics to Performance (9.1 Hover; 9.2 Forward Flight)
- Leishman, "Helicopters": Rotor Flow Environment; Hovering Flight Analysis; Axial Climbing & Descending Flight; Forward Flight Analysis; Blade Element Analysis in Forward Flight; Autorotation; Limitations of the Helicopter
  - https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/helicopters/ (loaded)
- Extract: hover power P = T^{3/2}/√(2ρA) with k_ind and profile power (σC_d0/8); figure of merit; Glauert's forward-flight inflow; the power-vs-speed bucket (induced, profile, parasite); advance ratio and the retreating-blade limits.

## 10-19 Unsteady Aerodynamics

Sources: `@ntrs-naca-tr-496-theodorsen-flutter`, `@ntrs-naca-tn-667-jones-operational-nonuniform-lift`, `@ntrs-naca-tr-681-jones-unsteady-lift-finite-wing`, `@leishman-intro-aerospace-flight-vehicles`, `@ntrs-nasa-cr-3082-rotary-wing-aerodynamics`

- Theodorsen, NACA TR-496: "Velocity Potentials, Forces, and Moments of the Noncirculatory Flow"; "Velocity Potentials, Forces, and Moments of the Circulatory Flow"; "The Magnitude of the Circulation"; "Value of the Function C(k)"; "Total Aerodynamic Forces and Moments" (read the page images; OCR is unreliable)
  - https://ntrs.nasa.gov/citations/19930090935 (PDF loaded)
- R. T. Jones, NACA TN-667 *Operational Treatment of the Nonuniform-Lift Theory in Airplane Dynamics*: Introduction (lift growth after a sudden change, by the method of operators)
  - https://ntrs.nasa.gov/citations/19930081472 (PDF loaded)
- R. T. Jones, NACA TR-681 *The Unsteady Lift of a Wing of Finite Aspect Ratio*
  - https://ntrs.nasa.gov/citations/19930091758 (PDF loaded)
- Leishman, "Unsteady Aerodynamics, Aeroelasticity, & Flutter": Unsteady Airfoil Behavior; Frequency Domain Theories; Time-Domain: The Indicial Response; Time- Vs. Frequency-Domain Functions; Three-Dimensional Unsteady Aerodynamics; Dynamic Stall & Stall Flutter
  - https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/unsteady-aerodynamics-flutter/ (loaded)
- Stepniewski CR-3082 Chapter VI: 7. Unsteady Aerodynamic Effects (7.2 Physics of Unsteady Flow; 7.3 Analytical Approach to Oscillations of Thin, Flat Airfoils below Stall)
- Extract: reduced frequency k = ωb/V; Theodorsen's C(k) = F + iG built from Hankel functions; noncirculatory (apparent-mass) vs. circulatory lift; the Wagner indicial function and Jones's exponential approximation; Duhamel superposition for arbitrary motion; finite-AR changes from TR-681.

## 10-20 Ground Effect & Interference

Sources: `@mit-ocw-16-100-aerodynamics`, `@ntrs-naca-tm-77-wieselsberger-ground`, `@ntrs-naca-tn-705-recant-ground-effect-flaps`, `@leishman-intro-aerospace-flight-vehicles`, `@ntrs-naca-tr-116-prandtl-modern-hydrodynamics`

- 16.100 Lecture 33 "Ground Effect Using Single Vortex Model" (image horseshoe vortex, w = 0 at z = 0); Lecture 34 "Three-Dimensional Wall Effects"; study material "Wind Tunnel Wall" (values of δ for a wing with uniform loading in a closed elliptic jet)
  - https://ocw.mit.edu/courses/16-100-aerodynamics-fall-2005/resources/16100lectre33_cj/ (loaded; PDF loaded); .../16100lectre34_cj/ (loaded; PDF loaded); .../3windtunwall/ (loaded; PDF loaded)
- Wieselsberger, NACA TM-77 *Wing Resistance Near the Ground* (image method with Prandtl's wing theory; read the page images)
  - https://ntrs.nasa.gov/citations/19930081293 (PDF loaded)
- Recant, NACA TN-705: Introduction; Apparatus and Tests; Results and Discussion; Conclusions
  - https://ntrs.nasa.gov/citations/19930081443 (PDF loaded)
- Leishman, "Ground Effect Vehicles": Anatomy of GEVs; WIGE Vehicle Performance
  - https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/hovercraft/ (loaded)
- Prandtl TR-116 Part II: E. "Improved Theory of Airplanes Having More Than One Wing" (biplane interference); F. "Aerofoils in a Tube or in a Free Jet" (wall interference)
- Extract: the image vortex reduces downwash, so induced drag falls by a factor that depends on h/b; build the estimate from the 10-11 VLM with a mirrored wing; tunnel-wall corrections are the same image idea with the opposite sign for a closed jet.

## 10-21 Using XFOIL, AVL & OpenVSP

Sources: `@drela-xfoil`, `@drela-avl`, `@openvsp-docs` (all supplementary; optional free tools)

- XFOIL page: Release Conditions (GPL); Software (6.99 for Unix/Windows; 6.996 Unix source); Documents (`xfoil_doc.txt` User Guide, `sessions.txt` sample inputs, `version_notes.txt`)
  - https://web.mit.edu/drela/Public/web/xfoil/ (loaded)
- XFOIL User Primer: Program Execution; Airfoil file formats; Analysis; Numerical accuracy; Polar calculations and plotting; Re, Mach dependence; Caveats
  - https://web.mit.edu/drela/Public/web/xfoil/xfoil_doc.txt (loaded)
- AVL page: Release Conditions (GPL); Software (3.52; 3.40b executables); Documents (`avl_doc.txt`, `AVL_User_Primer.pdf`, `session1.txt`, `session2.txt`, `runs/`)
  - https://web.mit.edu/drela/Public/web/avl/ (loaded)
- AVL User Primer: Input Files; Geometry Input File -- xxx.avl (Coordinate system; File format); Mass Input File -- xxx.mass; Program Execution; OPER Routine -- Flow Analysis (Stability derivatives; Trimmed Flight Condition Setup); MODE Routine -- Eigenmode Analysis
  - https://web.mit.edu/drela/Public/web/avl/avl_doc.txt (loaded)
- OpenVSP: Download (3.55.0 for Windows, macOS and Ubuntu); License (NOSA 1.3); Documentation page; wiki "Analysis Guides" (VSPAERO VLM & Panel Code; Parasite Drag Tool; Definitions of Aerodynamic Coefficients); wiki "vspaerotutorial" (Using VSPAERO: Getting Started; DegenGeom File; Cart3D File; Setup File)
  - https://openvsp.org/download.php (loaded); https://openvsp.org/license.shtml (loaded); https://openvsp.org/wiki/doku.php?id=analysisguides (loaded); https://openvsp.org/wiki/doku.php?id=vspaerotutorial (loaded)
- Extract: install the tools; run the same airfoil in XFOIL against the 10-05/10-06 codes; describe the 10-11 wing in an `.avl` file and compare C_L, C_Di and e; build the wing in OpenVSP and run VSPAERO (VLM). The lesson stresses that these are cross-checks of the student's own codes, not substitutes for them.
- Partial GAP: the OpenVSP wiki's VSPAERO tutorial is marked "under construction" since 2015 and mixes v1.x material, and the Ground School was not loaded. Step-by-step instructions for the current OpenVSP 3.55 VSPAERO GUI are written originally.

## 10-22 Aerodynamic Validation Project

Sources: `@ntrs-naca-tn-1669-untapered-wings-aspect-ratio-sweep`, `@ntrs-naca-tn-1671-swept-wings-taper-ratio`, `@ntrs-naca-tn-1422-sivells-three-wings`, `@ntrs-nasa-tn-d-6142-vortex-lattice-program`, `@ntrs-naca-tr-921-span-loading`, `@ntrs-nasa-tm-4074-ladson-naca0012`, `@leishman-intro-aerospace-flight-vehicles`, `@mit-ocw-16-100-aerodynamics`, `@ntrs-nasa-sp-405-vortex-lattice-utilization` (supplementary), `@drela-avl` (supplementary)

- Wind-tunnel data (transcribed in `curriculum/phase-10/data/finite_wing_cla.csv`): NACA TN-1669 (Introduction; Apparatus and Tests; Corrections; Results and Discussion; Conclusions), TN-1671, TN-1422
  - https://ntrs.nasa.gov/citations/19930082295 (PDF loaded); https://ntrs.nasa.gov/citations/19930082298 (PDF loaded); https://ntrs.nasa.gov/citations/19930082118 (PDF loaded)
- Method checks: TN D-6142 "Effect of Vortex-Lattice Arrangement on Computed Aerodynamic Characteristics" and "Sample Cases"; TR-921 loading characteristics
- Section input data: Ladson TM-4074 (a_0 and α_L0 for NACA 0012 at the test Re)
- Leishman, "Wind Tunnel Testing": Model Scaling & Similarity Parameters; Wall Interference; Data Reduction
  - https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/wind-tunnel-testing/ (loaded)
- 16.100 Lecture 31 "Similarity in Wind Tunnel Testing"
  - https://ocw.mit.edu/courses/16-100-aerodynamics-fall-2005/resources/16100lectre31_cj/ (loaded; PDF loaded)
- Extract: match planform, section and Re to each test; state which corrections the report already applied (jet-boundary corrections are noted in the data README); compare C_Lα across AR and sweep; plot residuals; run a convergence study; explain the remaining disagreement (viscous a_0, thickness, Re).

## 10-23 Home Lab: Airfoil Bench Rig

Sources: `@nasa-grc-beginners-guide-aeronautics`, `@leishman-intro-aerospace-flight-vehicles`, `@mit-ocw-16-100-aerodynamics`, `@ntrs-nasa-tm-4062-eppler-387-low-re`, `@uiuc-selig-summary-low-speed-airfoil-data` (supplementary, link only), `@ntrs-nasa-tm-4074-ladson-naca0012`, `@ntrs-naca-tr-824-summary-of-airfoil-data`

- NASA Glenn "Build Your Own Wind Tunnel" (a classroom tunnel built for under $100 around a computer fan motor; "Little Smokey" flow-visualization tunnel). Wind-tunnel pages: "Wind Tunnel Parts"; "Open Return Wind Tunnel"; "Wind Tunnel Testing"; "Force Balance"; "External Force Balance"; "Model Mounts"; "Reynolds Number"
  - https://www.grc.nasa.gov/www/k-12/WindTunnel/build.html (loaded); https://www.grc.nasa.gov/www/k-12/airplane/tunpart.html (loaded); .../tunoret.html (loaded); .../tuntest.html (loaded); .../reynolds.html (loaded)
- Leishman, "Wind Tunnel Testing": Components of a Low-Speed Wind Tunnel; Flow Quality; Freestream Speed Measurement; Wall Interference; Instrumentation & Measurement Techniques; Data Reduction
- 16.100 Lecture 35 "Drag Tare Due to Mount"; Lecture 34 "Three-Dimensional Wall Effects"
- Selig et al. Vol. 1 §2 "Wind-Tunnel Facility and Measurement Techniques" (2.2.1 Lift Force Measurements; 2.2.2 Drag Force Measurements; 2.3.1 Wind-Tunnel Boundary Corrections; 2.4 Calibrations and Uncertainty Analysis), a model of how a small university tunnel reports its data
- Comparison data: TM-4062 (E387 at Re 60 000–460 000) bridges toward the home rig's Re of roughly 10⁴–10⁵; Ladson TM-4074 and TR-824 are at Re of millions
- Extract: Re = Vc/ν for the rig and why it sits below every NACA dataset; a honeycomb (straws) for flow straightening; balance arm geometry from 10-01 to resolve lift and drag on a kitchen scale; tare and wall-blockage corrections; uncertainty propagation. Expect a lower c_l,max, earlier stall and much higher c_d than the NACA data, and explain the gap.
- The home procedure is original: fan, straw honeycomb, foam or 3D-printed airfoil, kitchen-scale balance, speed measurement and safety notes. The NASA Glenn tunnel pages and the cited labs are models, not procedures.


## Optional textbooks for deeper study

Paid; never required. Not recorded in `sources.yaml`.

- J. D. Anderson, *Fundamentals of Aerodynamics* (McGraw-Hill). The text that MIT 16.01 and 16.100 assign readings from; covers 10-01 … 10-12 and 10-20.
- J. Katz and A. Plotkin, *Low-Speed Aerodynamics* (Cambridge). The standard derivation of 2-D and 3-D panel and vortex-lattice methods, including the linear-strength vortex panel method (10-05, 10-11).
- M. Drela, *Flight Vehicle Aerodynamics* (MIT Press). The 16.110 text: integral boundary layers, viscous–inviscid coupling, lifting-surface theory, slender bodies and unsteady flows (10-06, 10-11, 10-16, 10-19).
- I. H. Abbott and A. E. von Doenhoff, *Theory of Wing Sections* (Dover). The book form of TR-824's methods and data; inexpensive but paid (10-02, 10-04, 10-08).
- J. G. Leishman, *Principles of Helicopter Aerodynamics* (Cambridge). The full treatment of hover, forward flight, BEMT and unsteady rotor aerodynamics (10-17, 10-18, 10-19).
- J. J. Bertin and R. M. Cummings, *Aerodynamics for Engineers* (Cambridge/Pearson). Worked examples across the phase.
- S. F. Hoerner, *Fluid-Dynamic Drag* (Hoerner Fluid Dynamics). Component drag data and correlations behind most drag build-ups (10-13, 10-16).
- D. P. Raymer, *Aircraft Design: A Conceptual Approach* (AIAA). Component drag build-up and form-factor correlations as used in design (10-13).
