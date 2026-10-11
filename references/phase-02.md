# Reference map: Phase 02: Mathematics II: Differential Equations, Transforms & Probability

Verified 2026-10-11 by loading the live pages; section numbers and titles are copied from the publishers' own chapter outlines or section pages. Every source here is **required-free** unless marked otherwise; sources are identified by an `@`-prefixed identifier and recorded in [sources.yaml](sources.yaml) with license, version, and caveats. **GAP** marks a lesson that no listed free source covers adequately: the lesson must carry original instruction there, and any "candidate" named is unverified until it is added to the manifest.

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

Sources: `@mit-ocw-18-06-linear-algebra`, `@mit-ocw-18-03-differential-equations`

- 18.06 Lecture 26: Complex matrices; fast fourier transform — https://ocw.mit.edu/courses/18-06-linear-algebra-spring-2010/resources/lecture-26-complex-matrices-fast-fourier-transform/ — DFT matrix F_n, its unitary structure, and the FFT factorization. Covers the discrete transform only.
- Continuous Fourier transform (definition, properties, convolution theorem, transform pairs, sampling/aliasing): **GAP** across all listed sources (18.03 Spring 2010 does Fourier *series* and Laplace only). Candidate free fills to verify: MIT OCW 18.03SC (Fall 2011) has no FT either; MIT OCW 6.003 "Signals and Systems" (Fall 2011) lectures on CT/DT Fourier transforms and sampling (CC BY-NC-SA).

## 02-11 PDEs: Heat, Wave & Laplace

Sources: `@mit-ocw-18-03-differential-equations`

- **GAP.** None of the listed sources cover PDEs: the 18.03 Spring 2010 lecture-notes page lists no heat/wave sessions (confirmed), and OpenStax Calculus stops at ODEs. Candidate free fills to verify: MIT OCW 18.303 "Linear Partial Differential Equations" (CC BY-NC-SA); Paul Dawkins' "Differential Equations" notes, chapter "Partial Differential Equations" (tutorial.math.lamar.edu; free to read but not CC-licensed).

## 02-15 Probability Fundamentals

Sources: `@openstax-introductory-statistics-2e`

- IS 3.1 Terminology — https://openstax.org/books/introductory-statistics-2e/pages/3-1-terminology — sample space, events, P(A), conditional probability notation.
- IS 3.2 Independent and Mutually Exclusive Events — .../pages/3-2-independent-and-mutually-exclusive-events.
- IS 3.3 Two Basic Rules of Probability — .../pages/3-3-two-basic-rules-of-probability — multiplication and addition rules.
- IS 3.4 Contingency Tables; 3.5 Tree and Venn Diagrams — .../pages/3-4-contingency-tables, .../pages/3-5-tree-and-venn-diagrams — Bayes-style reasoning via trees (Bayes' theorem itself is not a named section; note this when writing the lesson).
- Counting/combinatorics: PC Chapter 11 "Sequences, Probability and Counting Theory" has sections on counting principles and binomial theorem (chapter outline not fetched in this pass; confirm numbers before citing).

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
- IS 9.3 Distribution Needed for Hypothesis Testing; 9.4 Rare Events, the Sample, Decision and Conclusion; 9.5 Additional Information and Full Hypothesis Test Examples — slugs 9-3-distribution-needed-for-hypothesis-testing, 9-4-rare-events-the-sample-decision-and-conclusion, 9-5-additional-information-and-full-hypothesis-test-examples (9.3 page returned an error during verification; titles for 9.3-9.5 are inferred from URL slugs seen in Previous/Next links; confirm before publishing).
- IS 9.6 Hypothesis Testing of a Single Mean and Single Proportion — .../pages/9-6-hypothesis-testing-of-a-single-mean-and-single-proportion — worked p-value procedure (acceptance testing of components).
- Descriptive statistics (mean, SD, histograms): IS Chapter 2 "Descriptive Statistics" (not fetched in this pass).

## 02-18 Regression & Curve Fitting

Sources: `@openstax-precalculus-2e`, `@openstax-introductory-statistics-2e`, `@mit-ocw-18-06-linear-algebra`

- IS 12.1 Linear Equations; 12.2 Scatter Plots — https://openstax.org/books/introductory-statistics-2e/pages/12-1-linear-equations, .../pages/12-2-scatter-plots.
- IS 12.3 The Regression Equation — .../pages/12-3-the-regression-equation — least-squares line, r and r^2, residuals.
- IS 12.4 Testing the Significance of the Correlation Coefficient — .../pages/12-4-testing-the-significance-of-the-correlation-coefficient.
- IS 12.5 Prediction; 12.6 Outliers — .../pages/12-5-prediction, .../pages/12-6-outliers — extrapolation caveats, outlier detection (wind-tunnel data cleaning).
- 18.06 Lecture 16: Projection matrices and least squares — the linear-algebra derivation of the same fit; general polynomial/multi-parameter curve fitting via normal equations or QR (Lecture 17). PC 4.8 Fitting Exponential Models to Data for nonlinear-in-parameter fits via log transform.



## Lessons without a mapped free source

- **02-12 Complex Variables & Conformal Maps**, **02-13 Calculus of Variations**, **02-14 Tensor Notation**, **02-19 Random Processes & Power Spectra**, **02-20 Bayesian Inference**: not in this verification pass. Candidates to verify: MIT OCW 18.04 Complex Variables with Applications, 18.06/18.065 for tensors, 6.011/6.041 for random processes and Bayesian inference.
