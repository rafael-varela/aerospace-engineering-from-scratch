# Reference map: Phase 02: Mathematics II: Differential Equations, Transforms & Probability

Verified 2026-10-11 (02-10 continuous FT, 02-11 to 02-14, 02-15 counting, 02-17 9.3-9.5, 02-19, 02-20 added and verified 2026-10-13) by loading the live pages; section numbers and titles are copied from the publishers' own chapter outlines or section pages. Every source here is **required-free** unless marked otherwise; sources are identified by an `@`-prefixed identifier and recorded in [sources.yaml](sources.yaml) with license, version, and caveats. **GAP** marks a lesson that no listed free source covers adequately: the lesson must carry original instruction there, and any "candidate" named is unverified until it is added to the manifest.

Conventions: OpenStax section URLs follow `https://openstax.org/books/<book-slug>/pages/<n>-<m>-<title-slug>`; a URL written as `.../pages/...` continues the last full URL above it. MIT OCW video URLs follow `https://ocw.mit.edu/courses/<course-slug>/resources/lecture-<n>-<slug>/`. OpenStax and MIT OCW are **CC BY-NC-SA 4.0**, Leishman is **CC BY-NC-ND 4.0**: the curriculum links to them and writes its own explanations; nothing is copied or adapted.


## 02-01 First-Order ODEs

Sources: `@openstax-calculus-volume-2`, `@mit-ocw-18-03-differential-equations`

- C2 4.1 Basics of Differential Equations — https://openstax.org/books/calculus-volume-2/pages/4-1-basics-of-differential-equations — order, general vs particular solution, IVPs.
- C2 4.2 Direction Fields and Numerical Methods — .../pages/4-2-direction-fields-and-numerical-methods — slope fields, Euler's method.
- C2 4.3 Separable Equations; 4.5 First-Order Linear Equations — .../pages/4-3-separable-equations, .../pages/4-5-first-order-linear-equations — separation, integrating factor.
- C2 4.4 The Logistic Equation — .../pages/4-4-the-logistic-equation — autonomous equation, equilibria (optional).
- 18.03 L1 The Geometrical View of y'= f(x,y); L3 Solving First-order Linear ODEs; L5 First-order Autonomous ODEs — https://ocw.mit.edu/courses/18-03-differential-equations-spring-2010/resources/lecture-1-the-geometrical-view-of-y-f-x-y/ etc.; notes 1-4, 8 — existence/uniqueness, phase line, stability.

## 02-02 Second-Order Linear ODEs

Sources: `@openstax-calculus-volume-3`, `@mit-ocw-18-03-differential-equations`

- C3 7.1 Second-Order Linear Equations — https://openstax.org/books/calculus-volume-3/pages/7-1-second-order-linear-equations — characteristic equation, real/complex/repeated roots.
- C3 7.3 Applications — .../pages/7-3-applications — spring-mass damping cases, RLC circuits.
- 18.03 L9 Solving Second-order Linear ODE's with Constant Coefficients; L10 Continuation: Complex Characteristic Roots; L11 Theory of General Second-order Linear Homogeneous ODEs — .../resources/lecture-9-solving-second-order-linear-odes-with-constant-coefficients/ etc.; notes 11 "Modes and the characteristic polynomial", 12 "Good vibrations, damping conditions" — under/critical/over-damping, natural frequency, damping ratio.

## 02-03 Forced Response & Resonance

Sources: `@openstax-calculus-volume-3`, `@mit-ocw-18-03-differential-equations`

- C3 7.2 Nonhomogeneous Linear Equations — https://openstax.org/books/calculus-volume-3/pages/7-2-nonhomogeneous-linear-equations — undetermined coefficients, variation of parameters.
- 18.03 L12 Continuation: General Theory for Inhomogeneous ODEs; L13 Finding Particular Solutions to Inhomogeneous ODEs; L14 Interpretation of the Exceptional Case: Resonance — .../resources/lecture-14-interpretation-of-the-exceptional-case-resonance/ etc.
- 18.03 notes 13 "Exponential response formula, spring drive", 14 "Complex gain, dashpot drive", 15 "Operators, undetermined coefficients, resonance", 16 "Frequency response", 17 "LTI systems, superposition, RLC circuits" — https://ocw.mit.edu/courses/18-03-differential-equations-spring-2010/pages/lecture-notes/ — ERF, gain and phase lag, Bode-style frequency response (flutter/vibration context).

## 02-04 Systems of ODEs & Phase Portraits

Sources: `@openstax-calculus-volume-3`, `@mit-ocw-18-03-differential-equations`

- 18.03 L24 Introduction to First-order Systems of ODEs; L25 Homogeneous Linear Systems with Constant Coefficients; L26 Continuation: Repeated Real Eigenvalues; L27 Sketching Solutions of 2x2 Homogeneous Linear System with Constant Coefficients — .../resources/lecture-27-sketching-solutions-of-2x2-homogeneous-linear-system-with-constant-coefficients/ etc.
- 18.03 notes 32 "Linear systems and matrices", 33 "Eigenvalues, eigenvectors", 34 "Complex or repeated eigenvalues", 35 "Qualitative behavior of linear systems; phase plane" — node/saddle/spiral/center classification by trace-determinant.
- No OpenStax coverage (Calculus Volume 3 Chapter 7 stops at single second-order equations).

## 02-05 Matrix Exponential & Linear Systems

Sources: `@mit-ocw-18-06-linear-algebra`, `@mit-ocw-18-03-differential-equations`

- 18.03 L29 Matrix Exponentials; L28 Matrix Methods for Inhomogeneous Systems; L30 Decoupling Linear Systems with Constant Coefficients — .../resources/lecture-29-matrix-exponentials/ etc.; notes 36 "Normal modes and the matrix exponential".
- 18.06 Lecture 23: Differential equations and exp(At) — https://ocw.mit.edu/courses/18-06-linear-algebra-spring-2010/resources/lecture-23-differential-equations-and-exp-at/ — e^{At} = S e^{Lambda t} S^{-1}; stability iff Re(lambda)<0 (state-space form x' = Ax + Bu).

## 02-06 Nonlinear ODEs, Equilibria & Linearization

Sources: `@mit-ocw-18-03-differential-equations`

- 18.03 L31 Non-linear Autonomous Systems; L32 Limit Cycles; L33 Relation Between Non-linear Systems and First-order ODEs — .../resources/lecture-31-non-linear-autonomous-systems/ etc.
- 18.03 notes 37 "Nonlinear systems", 38 "Linearization near equilibria; the nonlinear pendulum", 39 "Limitations of the linear: limit cycles and chaos" — Jacobian at equilibria, classifying with linear theory, pendulum/phugoid analogies.

## 02-07 Complex Numbers & Phasors

Sources: `@openstax-precalculus-2e`, `@openstax-calculus-volume-2`, `@mit-ocw-18-03-differential-equations`

- 18.03 L6 Complex Numbers and Complex Exponentials — .../resources/lecture-6-complex-numbers-and-complex-exponentials/; notes 5 "Complex numbers, roots of unity", 6 "Complex exponentials; sinusoidal functions", 7 "Linear system response to exponential and sinusoidal input; gain, phase lag" — Euler's formula, sinusoid as Re(A e^{i omega t}), complex gain = phasor method.
- PC 8.5 Polar Form of Complex Numbers — https://openstax.org/books/precalculus-2e/pages/8-5-polar-form-of-complex-numbers — algebraic refresher.
- C2 6.4 Working with Taylor Series — proof of Euler's formula from series.

## 02-08 Laplace Transform

Sources: `@openstax-calculus-volume-2`, `@mit-ocw-18-03-differential-equations`

- 18.03 L19 Introduction to the Laplace Transform; L20 Derivative Formulas; L21 Convolution Formula; L22 Using Laplace Transform to Solve ODEs with Discontinuous Inputs; L23 Use with Impulse Inputs — .../resources/lecture-19-introduction-to-the-laplace-transform/ etc.
- 18.03 notes 23 "Step functions and delta functions", 24 "Step response, impulse response", 25 "Convolution", 26 "Laplace transform: basic properties", 27 "Application to ODEs", 28 "Second order equations; completing the squares", 29 "The pole diagram", 30 "The transfer function and frequency response" — transfer functions and pole locations (the bridge to control theory).
- Partial fractions for inverse transforms: C2 3.4 Partial Fractions (above). No OpenStax chapter on Laplace transforms.

## 02-09 Fourier Series

Sources: `@mit-ocw-18-06-linear-algebra`, `@mit-ocw-18-03-differential-equations`

- 18.03 L15 Introduction to Fourier Series; L16 Continuation: More General Periods; L17 Finding Particular Solutions via Fourier Series — .../resources/lecture-15-introduction-to-fourier-series/ etc.; notes 20 "Fourier series", 21 "Operations on Fourier series", 22 "Periodic solutions; resonance" — coefficients via orthogonality, even/odd, Gibbs, response to periodic forcing.
- 18.06 Lecture 24: Markov matrices; fourier series — .../resources/lecture-24-markov-matrices-fourier-series/ — Fourier series as projection onto an orthogonal basis (second half of lecture).
- No OpenStax coverage.

## 02-10 Fourier Transform & the FFT

Sources: `@mit-ocw-18-06-linear-algebra`, `@mit-ocw-18-03-differential-equations`, `@mit-ocw-6-003-signals-and-systems`

- 18.06 Lecture 26: Complex matrices; fast fourier transform — https://ocw.mit.edu/courses/18-06-linear-algebra-spring-2010/resources/lecture-26-complex-matrices-fast-fourier-transform/ — DFT matrix F_n, its unitary structure, and the FFT factorization. Covers the discrete transform only.
- 6.003 Lecture 14: Fourier Representations; Lecture 15: Fourier Series; Lecture 16: Fourier Transform — https://ocw.mit.edu/courses/6-003-signals-and-systems-fall-2011/resources/lecture-16-fourier-transform/ etc. — continuous Fourier transform as the aperiodic limit of the series; relation to the Laplace transform; transform pairs.
- 6.003 Lecture 19: Relations Among Fourier Representations; Lecture 20: Applications of Fourier Transforms — .../resources/lecture-19-relations-among-fourier-representations/, .../resources/lecture-20-applications-of-fourier-transforms/ — CT/DT, series/transform relations; convolution and filtering in the frequency domain.
- 6.003 Lecture 17: Discrete-Time (DT) Frequency Representations; Lecture 18: Discrete-Time (DT) Fourier Representations — the DTFT and DFT side, bridging to the 18.06 FFT lecture.
- 6.003 Lecture 21: Sampling; Lecture 22: Sampling and Quantization — .../resources/lecture-21-sampling/, .../resources/lecture-22-sampling-and-quantization/ — sampling theorem, aliasing, reconstruction (flight-test data acquisition).
- Matching lecture-note PDFs with the same titles: https://ocw.mit.edu/courses/6-003-signals-and-systems-fall-2011/pages/lecture-notes/ .

## 02-11 PDEs: Heat, Wave & Laplace

Sources: `@mit-ocw-18-303-linear-pdes`, `@mit-ocw-18-03-differential-equations`, `@dawkins-differential-equations-notes` (supplementary-free)

- 18.303 Lecture 1 "Overview of linear PDEs and analogies with matrix algebra"; Lecture 2 "Poisson's equation and eigenfunctions in 1d: Fourier sine series" — https://ocw.mit.edu/courses/18-303-linear-partial-differential-equations-analysis-and-numerics-fall-2014/pages/lecture-summaries/ (one-line summaries per lecture with handout PDFs; no videos) — the Laplace/Poisson equation and the eigenfunction viewpoint.
- 18.303 Lecture 8 "Music and wave equations, Separation of variables, in time and space"; Lecture 9 "Separation of variables in cylindrical geometries: Bessel functions"; Lecture 10 "General Dirichlet and Neumann boundary conditions" — same page — wave equation, separation of variables, boundary conditions.
- 18.303 Lecture 22 "Von Neumann analysis and the heat equation"; Lecture 23 "Algebraic properties of wave equations and unitary time evolution, Conservation of energy in a stretched string"; Lecture 25 "Traveling waves: D'Alembert's solution"; Lecture 26 "Group-velocity derivation and dispersion" — heat equation, d'Alembert, dispersion (the acoustics/structural-vibration context).
- 18.303 syllabus — .../pages/syllabus/ — lists the topics (heat/diffusion, wave, Poisson equations; separation of variables; Green's functions) and the recommended, non-free texts (Strang CSE; Olver). Fourier-series prerequisites: 18.03 notes 20-22 (see 02-09).
- Supplementary-free, not CC: Dawkins, Differential Equations, Chapter 9 Partial Differential Equations — https://tutorial.math.lamar.edu/Classes/DE/IntroPDE.aspx — sections "The Heat Equation", "The Wave Equation", "Terminology", "Separation of Variables", "Solving the Heat Equation", "Heat Equation with Non-Zero Temperature Boundaries", "Laplace's Equation", "Vibrating String", "Summary of Separation of Variables" (slugs TheHeatEquation.aspx, TheWaveEquation.aspx, PDETerminology.aspx, SeparationofVariables.aspx, SolvingHeatEquation.aspx, HeatEqnNonZero.aspx, LaplacesEqn.aspx, VibratingString.aspx, PDESummary.aspx) — a gentler worked-example walk-through; link only.

## 02-12 Complex Variables & Conformal Maps

Sources: `@mit-ocw-18-04-complex-variables`

- 18.04 Topic 1: Complex Algebra and the Complex Plane (sections 1.1-1.6.2, 1.6.3-1.14.2); Topic 2: Analytic Functions (2.1-2.5, 2.6-2.9) — https://ocw.mit.edu/courses/18-04-complex-variables-with-applications-spring-2018/pages/lecture-notes/ (PDF notes; no videos) — Cauchy-Riemann equations, analyticity.
- 18.04 Topic 3: Cauchy's Theorem; Topic 4: Cauchy Integral Formula; Topic 5: Harmonic Functions (5.1-5.4, 5.5-5.6) — harmonic functions as the link to Laplace's equation (potential flow).
- 18.04 Topic 6: Two-Dimensional Hydrodynamics and Complex Potentials (6.1-6.4, 6.5-6.7) — complex potential, streamlines, sources/vortices: the direct aerodynamics application.
- 18.04 Topic 7: Taylor and Laurent Series; Topic 8: Residue Theorem; Topic 9: Definite Integrals Using the Residue Theorem — residues (inverse Laplace transforms, Kutta-Joukowski integrals).
- 18.04 Topic 10: Conformal Transformations (10.1-10.4, 10.5-10.8, 10.8-10.11) — conformal maps; the Joukowski airfoil transformation belongs here (verify the specific example in the PDF before citing it by name).
- 18.04 Topic 12: Laplace Transform; Topic 11: Argument Principle — optional cross-links to 02-08 and to Nyquist stability in control.

## 02-13 Calculus of Variations

Sources: `@mit-ocw-16-323-principles-of-optimal-control`, `@mit-ocw-8-033-relativity` (supplementary-free)

- 16.323 Lecture 5: Calculus of variations — https://ocw.mit.edu/courses/16-323-principles-of-optimal-control-spring-2008/resources/lec5/ — functionals, first variation, Euler-Lagrange equation (resource description: "an alternative approach to solve general optimization problems for continuous systems").
- 16.323 Lecture 6: Calculus of variations applied to optimal control — .../resources/lec6/ — the bridge to trajectory optimisation (Pontryagin-style necessary conditions). Lecture Notes index: .../pages/lecture-notes/ (topics are listed in place of formal titles).
- 16.323 Lecture 2: Nonlinear optimization: constrained nonlinear optimization, Lagrange multipliers — .../resources/lec2/ — finite-dimensional Lagrange multipliers as the warm-up.
- Supplementary-free: 8.033 Session 8 "Variational Calculus" slides — https://ocw.mit.edu/courses/8-033-relativity-fall-2006/resources/lecture8_euler/ — extremisation of path integrals, Euler-Lagrange equations and constants of motion, brachistochrone (a short, undergraduate-level introduction; the 16.323 notes are graduate-level).
- Hamilton's principle / Lagrangian mechanics for flight-vehicle dynamics is not covered by either source beyond the 8.033 optional problems: original instruction.

## 02-14 Tensor Notation

Sources: `@mit-ocw-3-11-mechanics-of-materials`, `@mit-ocw-8-07-electromagnetism-ii` (supplementary-free), `@mit-ocw-18-06-linear-algebra`

- 3.11 Appendix "Matrix and Index Notation" — https://ocw.mit.edu/courses/3-11-mechanics-of-materials-fall-1999/resources/mit3_11f99_index/ — index notation, summation convention, Kronecker delta, matrix-index correspondence (short PDF module).
- 3.11 Module "Tensor Transformations" (under "General Concepts of Stress and Strain") — .../resources/mit3_11f99_trans/; module index https://ocw.mit.edu/courses/3-11-mechanics-of-materials-fall-1999/pages/modules/ — transformation of stress/strain tensors under rotation, principal values (the stress-tensor use in structures).
- Supplementary-free: 8.07 "Electromagnetism II, Lecture Notes 2" (resource description: rotations in index notation) — https://ocw.mit.edu/courses/8-07-electromagnetism-ii-fall-2012/resources/mit8_07f12_ln2/ — rotation matrices, Levi-Civita symbol and index gymnastics in a physics setting.
- 18.06 for the eigenvalue/orthogonal-matrix background (see 01-16 to 01-18). Neither source does general curvilinear tensor calculus (covariant/contravariant components, Christoffel symbols): if the lesson needs it, that part is original instruction; 18.065 (Spring 2018) was not found to cover tensor notation and is not cited.

## 02-15 Probability Fundamentals

Sources: `@openstax-introductory-statistics-2e`, `@openstax-precalculus-2e`, `@mit-ocw-18-05-introduction-to-probability-and-statistics`

- IS 3.1 Terminology — https://openstax.org/books/introductory-statistics-2e/pages/3-1-terminology — sample space, events, P(A), conditional probability notation.
- IS 3.2 Independent and Mutually Exclusive Events — .../pages/3-2-independent-and-mutually-exclusive-events.
- IS 3.3 Two Basic Rules of Probability — .../pages/3-3-two-basic-rules-of-probability — multiplication and addition rules.
- IS 3.4 Contingency Tables; 3.5 Tree and Venn Diagrams — .../pages/3-4-contingency-tables, .../pages/3-5-tree-and-venn-diagrams — Bayes-style reasoning via trees (Bayes' theorem itself is not a named section; note this when writing the lesson).
- Counting/combinatorics: PC 11.5 Counting Principles; 11.6 Binomial Theorem; 11.7 Probability — https://openstax.org/books/precalculus-2e/pages/11-introduction-to-sequences-probability-and-counting-theory (chapter outline, confirmed 2026-10-13: 11.1 Sequences and Their Notations, 11.2 Arithmetic Sequences, 11.3 Geometric Sequences, 11.4 Series and Their Notations, 11.5 Counting Principles, 11.6 Binomial Theorem, 11.7 Probability) — permutations, combinations, binomial coefficients.
- 18.05 Reading 1b: Counting and Sets; Reading 2: Probability: Terminology and Examples; Reading 3: Conditional Probability, Independence and Bayes' Theorem — https://ocw.mit.edu/courses/18-05-introduction-to-probability-and-statistics-spring-2022/pages/readings/ — Bayes' theorem stated and named (filling the gap noted above for IS 3.4-3.5).

## 02-16 Random Variables & Distributions

Sources: `@openstax-introductory-statistics-2e`

- IS 4.1 Probability Distribution Function (PDF) for a Discrete Random Variable; 4.2 Mean or Expected Value and Standard Deviation — https://openstax.org/books/introductory-statistics-2e/pages/4-1-probability-distribution-function-pdf-for-a-discrete-random-variable, .../pages/4-2-mean-or-expected-value-and-standard-deviation — expectation and variance.
- IS 4.3 Binomial Distribution; 4.6 Poisson Distribution — .../pages/4-3-binomial-distribution, .../pages/4-6-poisson-distribution — (4.4 Geometric, 4.5 Hypergeometric optional; 4.7 is a lab.)
- IS 5.1 Continuous Probability Functions; 5.2 The Uniform Distribution; 5.3 The Exponential Distribution — .../pages/5-1-continuous-probability-functions etc. — pdf/cdf, memoryless property (failure rates, reliability).
- IS 6.1 The Standard Normal Distribution; 6.2 Using the Normal Distribution — .../pages/6-1-the-standard-normal-distribution, .../pages/6-2-using-the-normal-distribution — z-scores, empirical rule (sensor noise models).
- IS 7.1 The Central Limit Theorem for Sample Means (Averages); 7.2 The Central Limit Theorem for Sums — .../pages/7-1-the-central-limit-theorem-for-sample-means-averages, .../pages/7-2-the-central-limit-theorem-for-sums — why averaged noise is Gaussian.

## 02-17 Statistics & Hypothesis Testing

Sources: `@openstax-introductory-statistics-2e`

- IS 8.1 A Single Population Mean using the Normal Distribution; 8.2 A Single Population Mean using the Student t-Distribution — https://openstax.org/books/introductory-statistics-2e/pages/8-1-a-single-population-mean-using-the-normal-distribution, .../pages/8-2-a-single-population-mean-using-the-student-t-distribution — confidence intervals, error bound.
- IS 9.1 Null and Alternative Hypotheses; 9.2 Outcomes and the Type I and Type II Errors — .../pages/9-1-null-and-alternative-hypotheses, .../pages/9-2-outcomes-and-the-type-i-and-type-ii-errors.
- IS 9.3 Probability Distribution Needed for Hypothesis Testing; 9.4 Rare Events, the Sample, Decision and Conclusion; 9.5 Additional Information and Full Hypothesis Test Examples — .../pages/9-3-probability-distribution-needed-for-hypothesis-testing, .../pages/9-4-rare-events-the-sample-decision-and-conclusion, .../pages/9-5-additional-information-and-full-hypothesis-test-examples (headings confirmed on the live pages 2026-10-13; the 9.3 slug previously recorded as 9-3-distribution-needed-for-hypothesis-testing was wrong, which is why it errored).
- IS 9.6 Hypothesis Testing of a Single Mean and Single Proportion — .../pages/9-6-hypothesis-testing-of-a-single-mean-and-single-proportion — worked p-value procedure (acceptance testing of components).
- Descriptive statistics (mean, SD, histograms): IS Chapter 2 "Descriptive Statistics" (not fetched in this pass).

## 02-18 Regression & Curve Fitting

Sources: `@openstax-precalculus-2e`, `@openstax-introductory-statistics-2e`, `@mit-ocw-18-06-linear-algebra`

- IS 12.1 Linear Equations; 12.2 Scatter Plots — https://openstax.org/books/introductory-statistics-2e/pages/12-1-linear-equations, .../pages/12-2-scatter-plots.
- IS 12.3 The Regression Equation — .../pages/12-3-the-regression-equation — least-squares line, r and r^2, residuals.
- IS 12.4 Testing the Significance of the Correlation Coefficient — .../pages/12-4-testing-the-significance-of-the-correlation-coefficient.
- IS 12.5 Prediction; 12.6 Outliers — .../pages/12-5-prediction, .../pages/12-6-outliers — extrapolation caveats, outlier detection (wind-tunnel data cleaning).
- 18.06 Lecture 16: Projection matrices and least squares — the linear-algebra derivation of the same fit; general polynomial/multi-parameter curve fitting via normal equations or QR (Lecture 17). PC 4.8 Fitting Exponential Models to Data for nonlinear-in-parameter fits via log transform.

## 02-19 Random Processes & Power Spectra

Sources: `@mit-ocw-6-011-signals-systems-and-inference`, `@mit-ocw-6-003-signals-and-systems`

- 6.011 Lecture 12: Probabilistic Models, Random Variables; Lecture 13: Vector Picture for First- and Second-Order Statistics, MMSE and LMMSE Estimation; Lecture 14: LMMSE Estimation, Orthogonality — https://ocw.mit.edu/courses/6-011-signals-systems-and-inference-spring-2018/pages/lecture-notes/ (PDF slides; no videos; list is out of numeric order on the page) — second-order statistics, correlation, linear estimation.
- 6.011 Lecture 15: Normal Equations, Random Processes; Lecture 16: Wide-sense Stationary Processes, LTI Filtering of WSS Processes; Lecture 17: LTI Filtering of WSS Processes — autocorrelation, WSS, response of LTI systems to random inputs (gust loads, sensor noise).
- 6.011 Lecture 18: Power Spectral Density; Lecture 19: Einstein-Wiener-Khinchin Theorem, PSD Applications, Modeling Filters — PSD, Wiener-Khinchin (turbulence spectra such as Dryden/von Karman are the aerospace application; not in the source).
- 6.011 Lecture 3: Energy Spectral Density; Lecture 20: Wiener Filtering; Lecture 22: Hypothesis Testing; Lecture 24: Matched Filtering — optional extensions (estimation and detection).
- The 6.011 Readings page assigns Oppenheim & Verghese, *Signals, Systems and Inference* (Pearson, not free): do not cite it as required. Deterministic Fourier background: 6.003 Lectures 16, 19-21 (see 02-10).

## 02-20 Bayesian Inference

Sources: `@mit-ocw-18-05-introduction-to-probability-and-statistics`

- 18.05 Reading 3: Conditional Probability, Independence and Bayes' Theorem — https://ocw.mit.edu/courses/18-05-introduction-to-probability-and-statistics-spring-2022/pages/readings/ (PDF readings) — Bayes' theorem, base-rate reasoning.
- 18.05 Reading 10a: Introduction to Statistics; Reading 10b: Maximum Likelihood Estimates — likelihood as the common ground with 02-17.
- 18.05 Reading 11: Bayesian Updating with Discrete Priors; Reading 12a: Bayesian Updating: Probabilistic Prediction; Reading 12b: Bayesian Updating: Odds — prior/likelihood/posterior tables, predictive probabilities, odds form.
- 18.05 Reading 13a: Bayesian Updating with Continuous Priors; Reading 13b: Notational conventions; Reading 14: Continuous Data with Continuous Priors; Reading 15: Conjugate priors: Beta and normal — continuous updating; the normal-normal conjugate pair is the scalar Kalman-filter measurement update.
- 18.05 Reading 16a: Choosing priors; Reading 16b: Probability intervals; Reading 20: Comparison of frequentist and Bayesian inference — priors, credible intervals, frequentist vs Bayesian (reliability and fault-diagnosis framing).
- Readings numbers are not contiguous on the page (no 8, 9, 21, 25). The same materials are on MIT Open Learning Library per the course page.
