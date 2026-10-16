# Reference map: Phase 03: Numerical Methods & Scientific Computing

Verified 2026-10-11 (03-06, 03-07, 03-17, 03-18 fills added and verified 2026-10-16) by loading the live pages; section numbers and titles are copied from the publishers' own chapter outlines or section pages. Every source here is **required-free** unless marked otherwise; sources are identified by an `@`-prefixed identifier and recorded in [sources.yaml](sources.yaml) with license, version, and caveats. **GAP** marks a lesson that no listed free source covers adequately: the lesson must carry original instruction there, and any "candidate" named is unverified until it is added to the manifest.

Conventions: OpenStax section URLs follow `https://openstax.org/books/<book-slug>/pages/<n>-<m>-<title-slug>`; a URL written as `.../pages/...` continues the last full URL above it. MIT OCW video URLs follow `https://ocw.mit.edu/courses/<course-slug>/resources/lecture-<n>-<slug>/`. OpenStax and MIT OCW are **CC BY-NC-SA 4.0**, Leishman is **CC BY-NC-ND 4.0**, Langtangen & Linge *Finite Difference Computing with PDEs* is **CC BY 4.0** and Langtangen *Finite Difference Computing with Exponential Decay Models* is **CC BY-NC 4.0** (Springer Open; cite by chapter DOI): the curriculum links to them and writes its own explanations; nothing is copied or adapted.


## 03-01 Floating Point & Numerical Error

Sources: `@ocw-18-330`, `@ocw-18-335j`, `@goldberg-1991`

- 18.335J Week 1 "Lecture 2: Floating-Point Arithmetic" (handout "Floating-Point Arithmetic"; "Some Myths
  about Floating-Point Arithmetic (PDF)"); Week 2 "Lecture 3: Floating-Point Summation and Backwards
  Stability", "Lecture 5: Condition Numbers". https://ocw.mit.edu/courses/18-335j-introduction-to-numerical-methods-spring-2019/
  Extract: IEEE 754 formats, machine epsilon, backward error, conditioning.
- D. Goldberg, "What Every Computer Scientist Should Know About Floating-Point Arithmetic" (Computing
  Surveys, Mar 1991; reprinted as Appendix D of Oracle Numerical Computation Guide).
  https://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html
  Extract: "Rounding Error" > "Relative Error and Ulps", "Cancellation"; "The IEEE Standard" > "Special Quantities".
- 18.330 "Lecture notes 1–2: Series and Sequences" (1.2 The big-O notation) for asymptotic notation only
  (it does not cover floating point).

## 03-02 Root Finding

Sources: `@ocw-18-330`, `@ocw-18-335j`, `@scipy-userguide`

- 18.330 "Lecture notes 9–10: Nonlinear Equations". https://ocw.mit.edu/courses/18-330-introduction-to-numerical-analysis-spring-2012/pages/lecture-notes/
  Extract: bisection, Newton, secant; convergence orders.
- 18.335J Week 1 "Lecture 1: Course Overview, Newton's Method for Root-Finding" (handout "Square Roots via Newton's Method").
- SciPy User Guide "Optimization (scipy.optimize)" > "Root finding" > "Scalar functions". https://docs.scipy.org/doc/scipy/tutorial/optimize.html

## 03-03 Interpolation & Splines

Sources: `@ocw-18-330`, `@ocw-18-335j`, `@scipy-userguide`

- 18.330 "Lecture notes 5–8: Interpolation". Extract: Lagrange/Newton forms, Runge phenomenon, splines.
- 18.330 "Lecture notes 21–25: Spectral Interpolation, Differentiation, Quadrature" (Chebyshev).
- 18.335J Week 13 "Lecture 35: Chebyshev Approximation".
- SciPy User Guide "Interpolation (scipy.interpolate)". https://docs.scipy.org/doc/scipy/tutorial/interpolate.html

## 03-04 Numerical Differentiation & Integration

Sources: `@ocw-18-330`, `@ocw-18-335j`, `@scipy-userguide`

- 18.330 "Lecture notes 3–4: Integrals as Sums and Derivatives as Differences". Extract: FD stencils,
  truncation vs roundoff trade-off, Newton–Cotes.
- 18.335J Week 12 "Lecture 33: Numerical Integration and the Convergence of the Trapezoidal Rule",
  "Lecture 34: Clenshaw-Curtis Quadrature"; Week 13 "Lecture 36: Integration with Weight Functions, and
  Gaussian Quadrature", "Lecture 37: Adaptive and Multidimensional Quadrature".
- SciPy User Guide "Integration (scipy.integrate)" > "General integration (quad)", "Gaussian quadrature",
  "Integrating using Samples". https://docs.scipy.org/doc/scipy/tutorial/integrate.html

## 03-05 ODE Integrators: Euler to RK4

Sources: `@ocw-18-330`, `@ocw-16-90`, `@scipy-userguide`

- 16.90 Unit 1 "Numerical Integration of ODEs": "1.2 Discretizing ODEs" (1.2.4 The Forward Euler Method,
  1.2.5 The Midpoint Method), "1.3 Order of Accuracy", "1.4 Convergence", "1.5 Zero Stability and the
  Dahlquist Equivalence Theorem", "1.9 Runge-Kutta Methods".
  https://ocw.mit.edu/courses/16-90-computational-methods-in-aerospace-engineering-spring-2014/
  Extract: local vs global error, consistency+zero-stability=convergence, RK4 Butcher table.
- 18.330 "Lecture notes 11–15: Methods for Ordinary Differential Equations".
- SciPy "Integration (scipy.integrate)" > "Ordinary differential equations (solve_ivp)".

## 03-06 Adaptive & Stiff Integrators

Sources: `@ocw-16-90`, `@scipy-userguide`, `@langtangen-fd-computing-exponential-decay`

- 16.90 Unit 1: "1.6 Systems of ODE's and Eigenvalue Stability", "1.7 Stiffness and Implicit Methods"
  (1.7.1 Stiffness, 1.7.3 Implicit Methods for Linear Systems of ODEs), "1.8 Multi-Step Methods"
  (Adams-Bashforth, Adams-Moulton, Backwards Differentiation Methods).
  Extract: stability regions, stiffness ratio, BDF; why explicit RK fails on stiff systems.
- SciPy "Ordinary differential equations (solve_ivp)" and "Solving a system with a banded Jacobian
  matrix" (methods RK45/DOP853 vs Radau/BDF/LSODA; rtol/atol).
- Langtangen, *Finite Difference Computing with Exponential Decay Models*, Chapter 3 "Generalizations":
  "3.2.9 The Odespy Software", "3.2.10 Example: Runge–Kutta Methods", "3.2.11 Example: Adaptive Runge–Kutta
  Methods" (Dormand–Prince 4(5): the difference between the 4th- and 5th-order results as the error
  indicator; tolerance sweep 1e-1 to 1e-7). https://link.springer.com/chapter/10.1007/978-3-319-29439-1_3
  Extract: embedded pair as error estimator, tolerance-driven step selection (demonstrated via Odespy).
- GAP narrowed: no verified free source derives the step-size controller itself (error-per-step vs
  error-per-unit-step, safety factor, step-rejection logic) or the embedded-pair Butcher tableaux as a
  titled section; 18.330 Chapter 5 (checked 2026-10-16: 5.2.1 Convergence, 5.2.2 Stability, 5.2.3
  Miscellaneous) does not cover it. Derive from 16.90 1.3/1.4 plus the SciPy solve_ivp docs.

## 03-07 Symplectic Integrators

Sources: `@langtangen-linge-fd-computing-pdes`, `@langtangen-fd-computing-exponential-decay`, `@ocw-18-085`

- Langtangen & Linge, *Finite Difference Computing with PDEs*, Chapter 1 "Vibration ODEs":
  "1.1.2 A Centered Finite Difference Scheme" (introduces Störmer's method / Verlet integration / Leapfrog
  in the text), "1.4 Analysis of the Numerical Scheme" (1.4.2 The Error in the Numerical Frequency, 1.4.5
  Convergence, 1.4.7 Stability), "1.5 Alternative Schemes Based on 1st-Order Equations" (1.5.1 The Forward
  Euler Scheme, 1.5.5 Runge-Kutta Methods), "1.6 Energy Considerations" (1.6.1 Derivation of the Energy
  Expression; "Growth of energy in the Forward Euler scheme"; 1.6.2 An Error Measure Based on Energy),
  "1.7 The Euler-Cromer Method" (1.7.1 Forward-Backward Discretization, 1.7.2 Equivalence with the Scheme
  for the Second-Order ODE, 1.7.4 The Störmer-Verlet Algorithm), "1.8 Staggered Mesh" (1.8.1 The
  Euler-Cromer Scheme on a Staggered Mesh); "1.9 Exercises and Problems": Problem 1.6 "Derive and
  investigate the velocity Verlet method". https://link.springer.com/chapter/10.1007/978-3-319-55456-3_1
  Extract: leapfrog/Störmer–Verlet/Euler–Cromer as the same scheme, bounded energy error vs Forward Euler
  energy growth, phase error. Caveat: the chapter uses the words "symplectic" only in the abstract; no
  heading uses "symplectic" or "Hamiltonian" (headings verified through 1.9; the rest of the chapter was
  not loaded).
- Langtangen, *Finite Difference Computing with Exponential Decay Models*, Chapter 3 "Generalizations":
  "3.2.4 Leapfrog Schemes" ("The ordinary Leapfrog scheme", "The filtered Leapfrog scheme"); Project 3.6
  "Implement and investigate the Leapfrog scheme". https://link.springer.com/chapter/10.1007/978-3-319-29439-1_3
- 18.085 "Lecture 10: Finite Differences in Time" (title confirmed 2026-10-16; the OCW page carries no
  lecture description). https://ocw.mit.edu/courses/18-085-computational-science-and-engineering-i-fall-2008/resources/lecture-10-finite-differences-in-time/
- Checked and not usable: 18.330 Chapter 5, 18.335J, 16.90 (lecture list re-verified via DSpace) and
  2.086 Spring 2013 "Lecture 15: Ordinary Differential Equations" have no symplectic/leapfrog section.
  GAP remains only for the general theory (symplectic maps, backward error analysis, modified Hamiltonian):
  the lesson writes that part itself.

## 03-08 Direct Linear Solvers

Sources: `@ocw-18-335j`, `@ocw-18-085`, `@ocw-2-092`, `@scipy-userguide`

- 18.335J Week 5 "Lecture 13: LU Factorization and Partial Pivoting", "Lecture 14: Cholesky Factorization
  and other Specialized Solvers"; Week 3 "Lecture 7: The SVD, its Applications, and Condition Numbers",
  "Lecture 8: Linear Regression and the Generalized SVD"; Week 9 "Lecture 24: Sparse-Direct Solvers".
  Extract: pivoting, growth factor, cost O(n^3), condition number and backward error.
- 18.085 "Lecture 3: Solving a Linear System"; "Lecture 7: Positive Definite Day".
- 2.092 "Physical explanation of Gauss elimination" (SES 13). https://ocw.mit.edu/courses/2-092-finite-element-analysis-of-solids-and-fluids-i-fall-2009/pages/lecture-notes/
- SciPy "Linear Algebra (scipy.linalg)".

## 03-09 Iterative & Sparse Solvers

Sources: `@ocw-18-335j`, `@scipy-userguide`

- 18.335J Week 7 "Lecture 18: Krylov Methods and the Arnoldi Algorithm", "Lecture 20: The GMRES Algorithm
  and Convergence of GMRES and Arnoldi"; Week 8 "Lecture 21: Preconditioning Techniques. The
  Conjugate-Gradient Method", "Lecture 22: Convergence of Conjugate Gradient", "Lecture 23: Biconjugate
  Gradient Algorithms". Extract: Krylov subspaces, CG for SPD, GMRES, preconditioning.
- SciPy "Sparse Arrays (scipy.sparse)"; "Sparse eigenvalue problems with ARPACK".

## 03-10 Nonlinear Systems

Sources: `@ocw-18-335j`, `@ocw-16-90`, `@scipy-userguide`

- 16.90 Unit 1 "1.7.4 Newton-Raphson Implement Implicit Methods on Nonlinear Problems", "1.7.5 Apply Newton-Rhapson" (titles as published).
- 18.335J Week 11 "Lecture 30: Quasi-Newton Methods and the BFGS Algorithm".
- SciPy "Optimization (scipy.optimize)" > "Root finding" > "Sets of equations", "Root finding for large
  problems", "Still too slow? Preconditioning." Extract: Newton–Krylov for large systems.

## 03-11 Finite Differences for PDEs

Sources: `@ocw-16-90`, `@ocw-18-085`

- 16.90 Unit 2 "Numerical Methods for PDEs": "2.2 Partial Differential Equations", "2.3 Introduction to
  Finite Difference Methods", "2.4 Analysis of Finite Difference Methods", "2.6 Upwinding and the CFL
  Condition", "2.7 Eigenvalue Stability of Finite Difference Methods".
  Extract: convection/diffusion model PDEs, consistency, CFL, von Neumann/eigenvalue stability.
- 18.085 "Lecture 23: Laplace's Equation", "Lecture 25: Fast Poisson Solver (part 1)".

## 03-12 Finite Volume Methods

Sources: `@ocw-16-90`, `@ocw-16-920j`

- 16.90 Unit 2 "2.5 Introduction to Finite Volume Methods", "2.6 Upwinding and the CFL Condition".
  Extract: conservation form, cell averages, fluxes, upwind flux.
- 16.920J "Numerical Methods for Partial Differential Equations (SMA 5212)" (Spring 2003) — topic
  "Finite Difference and Finite Volume Discretizations" (course page lists topics, not note titles).
  https://ocw.mit.edu/courses/16-920j-numerical-methods-for-partial-differential-equations-sma-5212-spring-2003/

## 03-13 Finite Element Basics

Sources: `@ocw-16-90`, `@ocw-18-085`, `@ocw-2-092`

- 16.90 Unit 2 "2.8 Method of Weighted Residuals", "2.9 Introduction to Finite Elements", "2.10 More on
  Finite Element Methods", "2.11 The Finite Element Method for Two-Dimensional Diffusion".
- 18.085 "Lecture 17: Finite Elements in 1D (part 1)", "Lecture 18: Finite Elements in 1D (part 2)",
  "Lecture 19: Quadratic/Cubic Elements", "Lecture 20: Element Matrices; 4th Order Bending Equations".
- 2.092 (Bathe) SES 4 "The principle of virtual work", SES 5 "The finite element formulation", SES 6–7
  "Finite element solution process". Extract: weak form, shape functions, assembly, stiffness matrix.

## 03-14 Unconstrained Optimization

Sources: `@ocw-18-335j`, `@scipy-userguide`

- 18.335J Week 9 "Lecture 25: Overview of Optimization Algorithms"; Week 10 "Lecture 28: Trust-Regions
  Methods and the CCSA Algorithm"; Week 11 "Lecture 30: Quasi-Newton Methods and the BFGS Algorithm",
  "Lecture 31: Derivation of the BFGS Update"; Week 12 "Lecture 32: Derivative-Free Optimization by
  Linear and Quadratic Approximations".
- SciPy "Optimization (scipy.optimize)" > "Unconstrained minimization": "Nelder-Mead Simplex algorithm",
  "Broyden-Fletcher-Goldfarb-Shanno algorithm", "Newton-Conjugate-Gradient algorithm", trust-region variants.

## 03-15 Constrained Optimization

Sources: `@ocw-18-335j`, `@ocw-16-90`, `@scipy-userguide`

- 18.335J Week 11 "Lecture 29: Lagrange Dual Problems"; Week 9 "Lecture 26: Adjoint Methods".
  Extract: KKT/Lagrangian intuition, duality, adjoint gradients.
- 16.90 Unit 3 "3.6 Introduction to Design Optimization".
- SciPy "Optimization (scipy.optimize)" > "Constrained minimization": "Trust-Region Constrained
  Algorithm (method='trust-constr')", "Sequential Least SQuares Programming (SLSQP) Algorithm".

## 03-16 Multidisciplinary Design Optimization (OpenMDAO)

Sources: `@openmdao-docs`

- OpenMDAO documentation — "Getting Started" https://openmdao.org/newdocs/versions/latest/getting_started/getting_started.html ;
  "Basic User Guide": "Single-Discipline Optimization", "Multidisciplinary Optimization", "Using the
  OpenMDAO Command Line Utility", "Recording and Reading Data".
  https://openmdao.org/newdocs/versions/latest/basic_user_guide/basic_user_guide.html
  Extract: Sellar problem, Groups/Components, N2 diagram, derivatives (total/partial).
  Caveat: docs pages show no version number and no license statement (project is Apache-2.0 on GitHub; verify).

## 03-17 Uncertainty Quantification

Sources: `@ocw-16-90`, `@jcgm-101-2008`, `@chaospy-docs`

- 16.90 Unit 3 "Probabilistic Methods and Optimization": "3.2 Review of Probability and Statistics",
  "3.3 Monte Carlo Methods" (3.3.4 Inversion Method for Sampling), "3.4 Error Estimates for the Monte Carlo
  Method" (3.4.4 Bootstrapping), "3.5 Variance Reduction Techniques for the Monte Carlo Method".
  Extract: MC estimator error ~1/sqrt(N), bootstrapping, variance reduction.
- JCGM 101:2008, "Propagation of distributions using a Monte Carlo method" (GUM Supplement 1).
  https://www.bipm.org/en/committees/jc/jcgm/publications (DOI https://doi.org/10.59161/JCGM101-2008)
- chaospy User guide (v4.3.13): "Main usage" > "Problem formulation", "Monte Carlo integration",
  "Point collocation" ("Select polynomial expansion", "Solve the linear regression problem"), "Pseudo-spectral
  projection" ("Generate quadrature nodes and weights", "Expansion of orthogonal polynomials", "Fourier
  coefficients"), "Intrusive Galerkin"; "Advanced topics" > "Generalized polynomial chaos" (sections "Monte
  Carlo integration", "Generalized polynomial chaos", "Point collocation method", "Pseudo-spectral projection",
  "Cholesky decomposition"), "Polynomial chaos Kriging".
  https://chaospy.readthedocs.io/en/latest/user_guide/main_usage/point_collocation.html ;
  https://chaospy.readthedocs.io/en/latest/user_guide/main_usage/pseudo_spectral_projection.html ;
  https://chaospy.readthedocs.io/en/latest/user_guide/advanced_topics/generalized_polynomial_chaos.html
  Extract: orthogonal polynomial basis, non-intrusive PCE by regression (collocation) and by quadrature
  (projection), moments/Sobol indices from coefficients; PCE-Kriging as a surrogate.
- GAP narrowed: software documentation, not a course; no free course or CC text with a titled polynomial
  chaos chapter was verified. Checked 2026-10-16: MIT 16.940 (uqgroup.mit.edu/teaching) lists no materials
  and is not on OCW; 16.90 Unit 3 (lectures 19-25) has no polynomial chaos lecture. The lesson derives
  the PCE theory itself.

## 03-18 Verification & Validation

Sources: `@ocw-18-335j`, `@nasa-seh-rev2`, `@nasa-std-7009b`, `@nasa-nparc-cfd-vv-tutorial`, `@langtangen-linge-fd-computing-pdes`, `@langtangen-fd-computing-exponential-decay`

- NASA-STD-7009B Change 1, "Standard for Models and Simulations" (document date 2024-03-05; "Internet
  Public", free PDF). https://standards.nasa.gov/standard/NASA/NASA-STD-7009
  Extract: V&V definitions, credibility assessment factors, reporting requirements.
- NASA Systems Engineering Handbook (NASA/SP-2016-6105 Rev 2) — "5.0 Product Realization" (verification
  and validation sections) https://www.nasa.gov/reference/systems-engineering-handbook/ ; Appendix C
  "How to Write a Good Requirement" https://www.nasa.gov/reference/appendix-c-how-to-write-a-good-requirement/
  Extract: requirement-writing checklist as the basis for test oracles.
- 18.335J Week 2 "Lecture 3: Floating-Point Summation and Backwards Stability" (code verification by
  backward error).
- NPARC Alliance CFD Verification and Validation Web Site, Tutorial (NASA Glenn; last updated 2021-02-10):
  "Verification Assessment" (steps "1. Examine the Computer Programming of the Code" ... "4. Examine Spatial
  (Grid) Convergence", "5. Examine Temporal Convergence", "6. Compare CFD Results to Highly Accurate
  Solutions") https://www.grc.nasa.gov/www/wind/valid/tutorial/verassess.html ; "Examining Spatial (Grid)
  Convergence" (sections "Order of Grid Convergence", "Asymptotic Range of Convergence", "Richardson
  Extrapolation", "Grid Convergence Index (GCI)", "Required Grid Resolution", "Effective Grid Refinement
  Ratio", "Example Grid Convergence Study", "VERIFY: A Fortran program to Perform Calculations Associated
  with a Grid Convergence Study") https://www.grc.nasa.gov/www/wind/valid/tutorial/spatconv.html ;
  also "Examining Iterative Convergence", "Examining Temporal Convergence", "Validation Assessment"
  (index https://www.grc.nasa.gov/www/wind/valid/tutorial/tutorial.html ).
  Extract: observed order p from three grids, Richardson-extrapolated value, GCI with safety factor,
  asymptotic-range check. Does not mention manufactured solutions.
- Langtangen & Linge, *Finite Difference Computing with PDEs*, Chapter 2 "Wave Equations": "2.2
  Verification" (2.2.2 Using an Analytical Solution of Physical Significance, 2.2.3 Manufactured Solution
  and Estimation of Convergence Rates, 2.2.4 Constructing an Exact Solution of the Discrete Equations),
  "2.3.4 Verification: Convergence Rates" https://link.springer.com/chapter/10.1007/978-3-319-55456-3_2 ;
  Chapter 1 "Vibration ODEs": "1.2.2 Verification" ("Testing very simple polynomial solutions", "Checking
  convergence rates"), "1.4.3 Empirical Convergence Rates and Adjusted ω".
- Langtangen, *Finite Difference Computing with Exponential Decay Models*, Chapter 3 "Generalizations":
  "3.1.4 Verifying a Constant Solution", "3.1.5 Verification via Manufactured Solutions", "3.1.6 Computing
  Convergence Rates" ("Estimating r", "Why convergence rates are important").
  https://link.springer.com/chapter/10.1007/978-3-319-29439-1_3
  Extract: MMS source-term construction, exact discrete solutions, empirical rate r from successive
  refinements as a pytest-able oracle.
