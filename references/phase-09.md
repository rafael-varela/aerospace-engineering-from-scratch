# Reference map: Phase 09: Signals, Systems & Control

**Note 2026-10-17:** first mapping of Phase 09. Pages were loaded with `curl` through the session proxy because the WebFetch tool could not resolve hosts (DNS failure). Verified 2026-10-17 by loading the live pages of MIT OCW 6.003 Fall 2011 (Lecture Notes page and lecture PDFs 1–13, 21, 22), MIT OCW 16.06 Fall 2012 (Lecture Notes and Readings pages, 25 of the 35 lecture PDFs), MIT OCW 16.30 Fall 2010 (course home, Syllabus, Lecture Notes page and 20 of the 24 lecture PDFs), MIT OCW 6.011 Spring 2018 (Lecture Notes page, lecture PDFs 1, 2, 4–11), MIT OCW 16.323 Spring 2008 (Lecture Notes page, lecture PDFs 3, 4, 13, 15, 16), MIT OCW 6.341 Fall 2005 (Lecture Notes page, lecture PDFs 7–9), MIT OCW 6.243J Fall 2003 (Lecture Notes page), the MIT OCW 18.03 and 18.06 lecture-video pages, the authors' download page and the 6th-printing PDF of Rawlings, Mayne & Diehl *Model Predictive Control* (2nd ed.), Smith's *Scientist and Engineer's Guide to DSP* (copyright page and chapters 3, 15, 16, 19, 20, 33), NTRS (NASA TN D-7843 and its papers by Deets and Szalai, records and PDFs), the SciPy 1.18.0 tutorial and API-reference pages for `scipy.signal` and `scipy.linalg`, the MicroPython PWM, Timer, ADC and RP2 quick-reference pages, and four SparkFun Learn tutorials. Seven ids were added to the manifest; `@mit-ocw-6-003-signals-and-systems`, `@mit-ocw-6-011-signals-systems-and-inference`, `@mit-ocw-16-323-principles-of-optimal-control`, `@mit-ocw-18-03-differential-equations`, `@mit-ocw-18-06-linear-algebra`, `@scipy-userguide`, `@micropython-docs` and `@sparkfun-learn-tutorials` were extended. Not cited: Åström & Murray *Feedback Systems* (fbswiki.org returned HTTP 500 "Cannot access the database" on four attempts, www.cds.caltech.edu failed TLS verification through the proxy, and a Wayback Machine copy reset the connection three times, so neither the content nor the exact personal-use terms could be checked); Michigan CTMS (ctms.engin.umich.edu returned a Cloudflare "Just a moment..." challenge with HTTP 403, so terms could not be read; it is also MATLAB-based); MIT OCW 6.302 Feedback Systems (only the Spring 2007 offering exists; it is graduate, has no Lecture Notes page — that URL returns 404 — and its topics are covered by 16.06 and 16.30); MIT OCW 16.31 as a separate course (Fall 2007 and Fall 2010 URLs return 404; the 16.30 Fall 2010 notes are headed "16.30/31"); MIT OCW 2.152 Nonlinear Control (course URL returned 404); NASA Glenn pages (no control-systems page relevant to these lessons was identified); MIT 18.335J (manifest id ocw-18-335j; nothing in it is needed beyond what 18.06 and SciPy give here).

Verified 2026-10-17 by loading the live pages; lecture titles are copied from the OCW Lecture Notes pages and section headings from the PDFs themselves. Every source here is **required-free** unless marked otherwise; sources are identified by an `@`-prefixed identifier and recorded in [sources.yaml](sources.yaml) with license, version, and caveats. **GAP** marks a lesson that no listed free source covers adequately: the lesson must carry original instruction there, and any "candidate" named is unverified until it is added to the manifest.

Conventions: MIT OCW resource pages follow `https://ocw.mit.edu/courses/<course-slug>/resources/<resource>/`, each with one PDF: 6.003 `6-003-signals-and-systems-fall-2011` / `mit6_003f11_lec<NN>`; 16.06 `16-06-principles-of-automatic-control-fall-2012` / `mit16_06f12_lecture_<n>`; 16.30 `16-30-feedback-control-systems-fall-2010` / `mit16_30f10_lec<NN>`; 6.011 `6-011-signals-systems-and-inference-spring-2018` / `mit6_011s18lec<n>`; 16.323 `16-323-principles-of-optimal-control-spring-2008` / `lec<n>`; 6.341 `6-341-discrete-time-signal-processing-fall-2005` / `lec<NN>`; 6.243J `6-243j-dynamics-of-nonlinear-systems-fall-2003` / `lec<n>_6243_2003`. A URL written `.../<resource>/` continues the prefix of the course named on that line. Where an OCW Lecture Notes page gives a topic and the PDF a different heading, both are given (page topic first). Headings in 16.30 PDFs are the page titles of the slide-style notes; 6.003 headings are slide titles. Licenses: MIT OCW is **CC BY-NC-SA 4.0**; SparkFun tutorials are CC BY-SA 4.0; MicroPython docs MIT License; SciPy BSD 3-Clause; NASA reports are US Government works ("Public Use Permitted"); Rawlings–Mayne–Diehl is "All rights reserved" (free PDF from the authors, link only) and Smith's DSP Guide is free for personal use only (link only) — both are **supplementary-free**, as are 6.341, 6.243J and the NASA F-8 report. The curriculum links to these sources and writes its own explanations; nothing is copied or adapted. SciPy is cited only for the test oracles (`scipy.signal`, `scipy.linalg`); the lessons implement every method from scratch.


## 09-01 Signals & Systems

Sources: `@mit-ocw-6-003-signals-and-systems`

- 6.003 Lecture 1 "Signals and Systems": The Signals and Systems Abstraction; Signals and Systems: Widely Applicable / Modular / Hierarchical; Example System: Leaky Tank; Analysis of the Leaky Tank; Leaky Tanks and Capacitors
  - https://ocw.mit.edu/courses/6-003-signals-and-systems-fall-2011/resources/mit6_003f11_lec01/ (loaded; PDF loaded)
- Lecture 2 "Discrete-Time (DT) Systems": Difference Equations; Operator Notation; Operator Algebra; Example: Accumulator; Finite and Infinite Impulse Responses
  - .../mit6_003f11_lec02/ (loaded; PDF loaded)
- Lecture 4 "Continuous-Time (CT) Systems": Elementary Building-Block Signals; Unit-Impulse Signal; Unit-Impulse and Unit-Step Signals
  - .../mit6_003f11_lec04/ (loaded; PDF loaded)
- Lecture 8 "Convolution": Superposition; Linearity; Time-Invariance
  - .../mit6_003f11_lec08/ (loaded; PDF loaded)
- Extract: a system as an operator on whole signals; the accumulator and delay as the canonical DT operators; linearity and time invariance as the two properties that make superposition (and therefore convolution) work.
- Partial GAP: no verified free source defines causality, memory (static vs. dynamic), BIBO stability, signal energy and power, or even/odd decomposition as a property checklist, and none tests a time-varying gain or a squarer. The property definitions, the numerical property tests on sampled operators, and the analytic verdicts are original.

## 09-02 Convolution & Impulse Response

Sources: `@mit-ocw-6-003-signals-and-systems`, `@mit-ocw-18-03-differential-equations`, `@scipy-userguide`

- 6.003 Lecture 8 "Convolution": Responses to arbitrary signals; Superposition; Structure of Superposition; Structure of Convolution; DT Convolution: Summary; CT Convolution; Applications of Convolution; Impulse Response: Summary
  - https://ocw.mit.edu/courses/6-003-signals-and-systems-fall-2011/resources/mit6_003f11_lec08/ (loaded; PDF loaded)
- 6.003 Lecture 4 "Continuous-Time (CT) Systems": Impulse Response of Acyclic CT System; Mass and Spring System
- 18.03 Lecture 21: Convolution Formula; Lecture 23: Use with Impulse Inputs
  - https://ocw.mit.edu/courses/18-03-differential-equations-spring-2010/resources/lecture-21-convolution-formula/ (listed on the loaded Video Lectures page); .../resources/lecture-23-use-with-impulse-inputs/ (pattern)
- SciPy tutorial "Signal Processing (scipy.signal)": Filtering → Convolution/Correlation (oracle only)
  - https://docs.scipy.org/doc/scipy/tutorial/signal.html (loaded)
- Extract: y = h * x as a sum of shifted, scaled impulse responses; the flip-and-slide picture; zero padding so FFT multiplication gives linear rather than circular convolution; continuous convolution as the Riemann-sum limit; step response as the running integral of h.
- Partial GAP: the RC IMU front end and first-order servo examples, and the FFT-vs-direct cost comparison, are original.

## 09-03 Transfer Functions & Block Diagrams

Sources: `@mit-ocw-16-06-principles-of-automatic-control`, `@mit-ocw-16-30-feedback-control-systems`, `@mit-ocw-6-003-signals-and-systems`, `@scipy-userguide`

- 16.06 Lecture 2 "Block diagrams, the effect of feedback" (PDF: reasons for automatic control; components of a typical control system — plant, actuator, sensor, sensor noise; open-loop vs. closed-loop control)
  - https://ocw.mit.edu/courses/16-06-principles-of-automatic-control-fall-2012/resources/mit16_06f12_lecture_2/ (loaded; PDF loaded)
- 16.06 Lecture 4 "Block diagram manipulations, Mason' rule" [sic] (PDF heading "Block Diagram Manipulations": the single-loop gain rule; "Mason's Rule")
  - .../mit16_06f12_lecture_4/ (loaded; PDF loaded)
- 16.30 Topic 2 "Basic root locus: analysis and examples": Aircraft Longitudinal Control; First Short Period Autopilot (the pitch-attitude loop with an elevator servo)
  - https://ocw.mit.edu/courses/16-30-feedback-control-systems-fall-2010/resources/mit16_30f10_lec02/ (loaded; PDF loaded)
- 6.003 Lecture 6 "Laplace Transform": Solving Differential Equations with Laplace Transforms; Properties of Laplace Transforms; Lecture 10 "Feedback and Control": Structure of a Control Problem; Analysis of wallFinder System: Block Diagram; Analysis of wallFinder System: System Function
  - https://ocw.mit.edu/courses/6-003-signals-and-systems-fall-2011/resources/mit6_003f11_lec06/ (loaded; PDF loaded); .../mit6_003f11_lec10/ (loaded; PDF loaded)
- SciPy API reference "Signal processing (scipy.signal)": "Continuous-time linear systems" (`lti`, `lsim`, `step`) as oracle
  - https://docs.scipy.org/doc/scipy/reference/signal.html (loaded)
- Extract: G(s) as the Laplace transform of the impulse response; series, parallel and feedback rules; L/(1 + L); poles, zeros and DC gain from polynomial coefficients; reading the 16.30 short-period autopilot diagram as the phase's pitch loop.
- Partial GAP: the polynomial-algebra reducer and its check by simulating the unreduced diagram are original; the pitch-loop numbers are a representative model written for this lesson.

## 09-04 Time Response & Poles

Sources: `@mit-ocw-16-06-principles-of-automatic-control`, `@mit-ocw-6-003-signals-and-systems`

- 16.06 Lecture 5 "Dynamic response of closed-loop systems" (PDF heading "Dynamic Response"); Lecture 6 "Time-domain specifications" (PDF: the step response of H(s); desired pole locations); Lecture 7 "Effect of zeros" (PDF heading "Effects of Zeros on Step Response", including "Non-Minimum Phase Zeros")
  - https://ocw.mit.edu/courses/16-06-principles-of-automatic-control-fall-2012/resources/mit16_06f12_lecture_5/ (loaded; PDF loaded); .../mit16_06f12_lecture_6/ (loaded; PDF loaded); .../mit16_06f12_lecture_7/ (loaded; PDF loaded)
- 6.003 Lecture 3 "Feedback, Poles, and Fundamental Modes": Geometric Growth: Poles; Second-Order Systems; Partial Fractions; Complex Poles; Lecture 4: CT Feedback; Convergent and Divergent Poles
  - https://ocw.mit.edu/courses/6-003-signals-and-systems-fall-2011/resources/mit6_003f11_lec03/ (loaded; PDF loaded)
- 6.003 Lecture 6 "Laplace Transform": Initial Value Theorem; Final Value Theorem
- Extract: rise time, peak time, overshoot and settling time against ζ and ω_n; the desired-pole region in the s-plane; dominant poles; how a left-half-plane zero adds overshoot and a right-half-plane zero gives undershoot; partial fractions as a sum of modes.

## 09-05 Frequency Response & Bode Plots

Sources: `@mit-ocw-6-003-signals-and-systems`, `@mit-ocw-16-06-principles-of-automatic-control`, `@mit-ocw-16-30-feedback-control-systems`, `@scipy-userguide`

- 6.003 Lecture 9 "Frequency Response": Eigenfunctions; Complex Exponentials; Rational System Functions; Vector Diagrams; Example: Mass, Spring, and Dashpot
  - https://ocw.mit.edu/courses/6-003-signals-and-systems-fall-2011/resources/mit6_003f11_lec09/ (loaded; PDF loaded)
- 6.003 Lecture 11 "Continuous-Time (CT) Frequency Response and Bode Plot": Frequency Response: H(s)|s←jω; Asymptotic Behavior: Isolated Zero; Asymptotic Behavior: Isolated Pole; Bode Plot: Adding Instead of Multiplying; Bode Plot: dB; Bode Plot: Accuracy; Frequency Response of a High-Q System
  - .../mit6_003f11_lec11/ (loaded; PDF loaded)
- 16.06 Lecture 16 "Frequency response design" (PDF: The Frequency Response; Bode Diagram); Lecture 17 "Bode plot problems" (PDF heading "Frequency Response Design": Bode Plot Construction); Lecture 18 "Bode plots (cont.)"
  - https://ocw.mit.edu/courses/16-06-principles-of-automatic-control-fall-2012/resources/mit16_06f12_lecture_17/ (loaded; PDF loaded); .../mit16_06f12_lecture_18/ (listed on the loaded Lecture Notes page)
- 16.30 Topic 3 "Frequency response methods": Frequency Response Function; Bode Example
  - https://ocw.mit.edu/courses/16-30-feedback-control-systems-fall-2010/resources/mit16_30f10_lec03/ (loaded; PDF loaded)
- SciPy API reference: "Continuous-time linear systems" (`bode`, `freqs`) as oracle
- Extract: G(jω) as the steady-state complex gain on a sinusoid; dB and log-frequency make factors add; break-frequency asymptotes and their worst-case error; resonant peak of a lightly damped pair; phase unwrapping.

## 09-06 Sampling & Aliasing

Sources: `@mit-ocw-6-003-signals-and-systems`, `@mit-ocw-6-341-discrete-time-signal-processing`, `@smith-dsp-guide`, `@scipy-userguide`

- 6.003 Lecture 21 "Sampling": Sampling and Reconstruction; The Sampling Theorem; CT Model of Sampling and Reconstruction; Aliasing; Anti-Aliasing Filter
  - https://ocw.mit.edu/courses/6-003-signals-and-systems-fall-2011/resources/mit6_003f11_lec21/ (loaded; PDF loaded)
- 6.003 Lecture 22 "Sampling and Quantization": Aliasing; Anti-Aliasing Filter; Quantization; Dithering; Discrete-Time Sampling (Resampling)
  - .../mit6_003f11_lec22/ (loaded; PDF loaded)
- 6.341 Lecture 8 "DT Filter Design: IIR Filters" (supplementary): IIR Filter Design; Butterworth Filters (maximally flat, monotonic magnitude); Chebyshev Filters; Elliptic Filters
  - https://ocw.mit.edu/courses/6-341-discrete-time-signal-processing-fall-2005/resources/lec08/ (loaded; PDF loaded)
- Smith, DSP Guide Ch. 3 "ADC and DAC" (supplementary, personal use): Quantization; The Sampling Theorem; Digital-to-Analog Conversion; Analog Filters for Data Conversion (Bessel, Butterworth, Chebyshev); Ch. 20 "Chebyshev Filters": The Chebyshev and Butterworth Responses
  - https://www.dspguide.com/ch3.htm (loaded; CH3.PDF loaded); https://www.dspguide.com/ch20.htm (loaded; CH20.PDF loaded)
- SciPy tutorial "Signal Processing (scipy.signal)": Filter Design → Analog Filter Design; API reference `butter` (analog=True) as oracle
- Extract: the folding (alias) frequency formula; why the anti-alias filter must act before the ADC; Butterworth poles equally spaced on a circle and order selection from a stopband requirement; sinc interpolation as ideal reconstruction.
- Partial GAP: no verified free source sizes an anti-alias filter for a MEMS IMU against rotor or engine vibration lines; the IMU scenario, its synthetic record and the order-selection walk-through are original.

## 09-07 Digital Filters

Sources: `@mit-ocw-6-341-discrete-time-signal-processing`, `@smith-dsp-guide`, `@mit-ocw-6-003-signals-and-systems`, `@mit-ocw-16-06-principles-of-automatic-control`, `@nasa-tn-d-7843-f8-digital-fly-by-wire`, `@scipy-userguide`

- 6.341 Lecture 7 "IIR, FIR Filter Structures": Signal Flow Graphs; IIR Filter Structures; Transposed Forms; Cascade Form; FIR Filter Structures; Direct Form; Direct Form Structure for Linear-Phase FIR Systems
  - https://ocw.mit.edu/courses/6-341-discrete-time-signal-processing-fall-2005/resources/lec07/ (loaded; PDF loaded)
- 6.341 Lecture 8 "DT Filter Design: IIR Filters": IIR Filter Design (the bilinear transformation and its tan frequency warping); Butterworth Filters
- 6.341 Lecture 9 "Filter Design: FIR Filters": FIR Filter Design by Windowing (incl. Kaiser windows); The Parks-McClellan Algorithm (not needed)
  - .../lec09/ (loaded; PDF loaded)
- Smith, DSP Guide (supplementary, personal use): Ch. 15 "Moving Average Filters" (Implementation by Convolution; Frequency Response; Recursive Implementation); Ch. 16 "Windowed-Sinc Filters" (Strategy of the Windowed-Sinc; Designing the Filter; Examples of Windowed-Sinc Filters); Ch. 19 "Recursive Filters" (The Recursive Method; Single Pole Recursive Filters; Phase Response); Ch. 33 "The z-Transform" (Analysis of Recursive Systems; Cascade and Parallel Stages; Chebyshev-Butterworth Filter Design)
  - https://www.dspguide.com/ch16.htm (loaded; CH16.PDF loaded); .../ch15.htm, .../ch19.htm, .../ch33.htm (loaded; PDFs loaded)
- 6.003 Lecture 5 "Z Transform": Delay Property; Solving Difference Equations with Z Transforms; Lecture 2: Finite and Infinite Impulse Responses
  - https://ocw.mit.edu/courses/6-003-signals-and-systems-fall-2011/resources/mit6_003f11_lec05/ (loaded; PDF loaded)
- 16.06 Lecture 31 (PDF heading "Example (continued)"): prewarping of frequencies in a Tustin design
  - https://ocw.mit.edu/courses/16-06-principles-of-automatic-control-fall-2012/resources/mit16_06f12_lecture_31/ (loaded; PDF loaded)
- Aerospace context: Szalai, "Flight Test Experience With the F-8 Digital Fly-By-Wire System" (TN D-7843 paper 6), CONTROL SYSTEM PERFORMANCE → "Stability Augmentation System Mode" (a discrete washout filter formed by the bilinear transformation)
  - https://ntrs.nasa.gov/citations/19750010179 (loaded; PDF loaded)
- SciPy API reference: `freqz`, `firwin`, `bilinear` as oracles
- Extract: FIR = finite convolution, IIR = recursion; windowed-sinc design and the window trade-off (transition width vs. ripple); the moving average as the simplest low-pass and its poor stopband; bilinear transform with prewarping at the cutoff; direct-form difference equation and its impulse response.
- Partial GAP: the gyro-record scenario (synthetic) and the comparison of filter group delay against a control-loop phase budget are original.

## 09-08 Feedback & Steady-State Error

Sources: `@mit-ocw-16-06-principles-of-automatic-control`, `@mit-ocw-16-30-feedback-control-systems`, `@mit-ocw-6-003-signals-and-systems`

- 16.06 Lecture 9 "Effect of noise, steady-state errors" (PDF heading "Unity Feedback Control With Noise": S(s) "Sensitivity transfer function" and T(s) "Complementary Sensitivity transfer function"; the tracking vs. noise trade-off; Steady-State Errors; system Type)
  - https://ocw.mit.edu/courses/16-06-principles-of-automatic-control-fall-2012/resources/mit16_06f12_lecture_9/ (loaded; PDF loaded)
- 16.06 Lecture 2 "Block diagrams, the effect of feedback" (reduce the effects of disturbances and plant variations)
- 16.30 Topic 2: Performance Issues (free integrators and steady-state error); Topic 23 "Anti-windup": Cruise control design (an automotive cruise loop with saturation)
  - https://ocw.mit.edu/courses/16-30-feedback-control-systems-fall-2010/resources/mit16_30f10_lec23/ (loaded; PDF loaded)
- 6.003 Lecture 12 "Continuous-Time (CT) Feedback and Control, Part 1": Feedback and Control; Dominant Pole; Improving Performance; Motor Controller
  - https://ocw.mit.edu/courses/6-003-signals-and-systems-fall-2011/resources/mit6_003f11_lec12/ (loaded; PDF loaded)
- 6.003 Lecture 6: Final Value Theorem
- Extract: closed-loop TF from L = KG; system type and error constants Kp, Kv, Ka; e_ss for step and ramp by the final-value theorem; disturbance-to-output TF; S + T = 1.
- Partial GAP: the airspeed-hold (aircraft cruise) model and its numbers are original.

## 09-09 Root Locus

Sources: `@mit-ocw-16-06-principles-of-automatic-control`, `@mit-ocw-16-30-feedback-control-systems`

- 16.06 Lecture 8 "The Routh criterion" (PDF: The Routh Stability Criterion; The Routh Array)
  - https://ocw.mit.edu/courses/16-06-principles-of-automatic-control-fall-2012/resources/mit16_06f12_lecture_8/ (loaded; PDF loaded)
- 16.06 Lecture 11 "The root locus method"; Lecture 12 "Root locus rules"; Lecture 13 "Root locus rules, lead compensation" (PDF heading "Root Locus Rules (cont-d)", with a Routh array for the jω crossing); Lecture 15 "Zero degree root locus"
  - .../mit16_06f12_lecture_11/, .../mit16_06f12_lecture_12/, .../mit16_06f12_lecture_13/ (loaded; PDFs loaded); .../mit16_06f12_lecture_15/ (listed on the loaded Lecture Notes page)
- 16.30 Topic 2 "Basic root locus: analysis and examples": Root Locus Basics; Root Locus Analysis (where the locus starts, behaviour as K → ∞, the real-axis rule); Root Locus Examples; Dynamic Compensation; Classic Root Locus Approach
  - https://ocw.mit.edu/courses/16-30-feedback-control-systems-fall-2010/resources/mit16_30f10_lec02/ (loaded; PDF loaded)
- Extract: angle and magnitude conditions; the construction rules (real-axis segments, asymptotes and centroid, breakaway, jω crossing by Routh); picking K for a damping-ratio target.
- Partial GAP: the altitude-hold loop is original; no verified free source treats an outer altitude loop around an inner pitch loop by root locus.

## 09-10 Nyquist & Stability Margins

Sources: `@mit-ocw-16-06-principles-of-automatic-control`, `@mit-ocw-16-30-feedback-control-systems`, `@mit-ocw-6-003-signals-and-systems`

- 16.06 Lecture 20 "The Nyquist stability criterion" (PDF heading "Bode Plots With Complex Poles": neutral stability; The Argument Principle); Lecture 21 "The Nyquist stability criterion (cont.)" (PDF: P, the number of unstable poles; Examples)
  - https://ocw.mit.edu/courses/16-06-principles-of-automatic-control-fall-2012/resources/mit16_06f12_lecture_20/ (loaded; PDF loaded); .../mit16_06f12_lecture_21/ (loaded; PDF loaded)
- 16.06 Lecture 22 "Nyquist with poles on imaginary axis" (PDF also: Nyquist Plot of Open Loop Unstable System); Lecture 23 "Stability margins, Bode gain-phase theorem" (PDF: Stability Margins; Relationship between PM and damping)
  - .../mit16_06f12_lecture_22/, .../mit16_06f12_lecture_23/ (loaded; PDFs loaded)
- 16.30 Topic 3 "Frequency response methods": Frequency Stability Tests; Gain and Phase Margins; Nyquist Stability; FR: Summary
  - https://ocw.mit.edu/courses/16-30-feedback-control-systems-fall-2010/resources/mit16_30f10_lec03/ (loaded; PDF loaded)
- 6.003 Lecture 10 "Feedback and Control": Analysis of wallFinder System: Adding Sensor Delay; Destabilizing Effect of Delay
- Extract: the argument principle and Z = N + P; indenting the contour around imaginary-axis poles; gain and phase margins from the −1 point; PM ≈ 100ζ as a rule of thumb.
- Partial GAP: none of the loaded notes defines the delay margin (PM divided by crossover frequency) or counts encirclements numerically by winding number; both are original.

## 09-11 PID Design & Tuning

Sources: `@mit-ocw-16-06-principles-of-automatic-control`, `@mit-ocw-16-30-feedback-control-systems`

- 16.06 Lecture 10 "PID control" (PDF: Proportional (P) control; Derivative (D) control; Integral (I) control; choosing gains to match a desired characteristic polynomial)
  - https://ocw.mit.edu/courses/16-06-principles-of-automatic-control-fall-2012/resources/mit16_06f12_lecture_10/ (loaded; PDF loaded)
- 16.30 Topic 2: First Short Period Autopilot; Dynamic Compensation (proportional, integral and derivative feedback); Pole Placement
- 16.30 Topic 23 "Anti-windup" (heading "Analysis of Nonlinear Systems"): Nonlinear System Analysis; Jacobian Linearization; Cruise control design; Nonlinear Simulation; Effects of the saturation; Anti-windup Logic; Anti-windup Scheme; Anti-windup Summary
  - https://ocw.mit.edu/courses/16-30-feedback-control-systems-fall-2010/resources/mit16_30f10_lec23/ (loaded; PDF loaded)
- Extract: what each of P, I and D does to the characteristic equation; PID gains by pole placement; integrator windup under actuator saturation and the anti-windup fix.
- GAP: no verified free source covers Ziegler–Nichols tuning, the filtered derivative (N or τ_f), back-calculation vs. conditional-integration (clamping) anti-windup by name, or discrete PID implementation. These are original. Candidate: Åström & Murray *Feedback Systems* Ch. 10 "PID Control", once its site and terms can be verified.

## 09-12 Lead-Lag Compensation

Sources: `@mit-ocw-16-06-principles-of-automatic-control`, `@mit-ocw-16-30-feedback-control-systems`

- 16.06 Lecture 13 "Root locus rules, lead compensation" (PDF: Lead Compensator Example); Lecture 14 "Lag compensation" (PDF heading "Lead Compensator Example (cont-d)": the angle condition; a note on closed-loop poles and zeros)
  - https://ocw.mit.edu/courses/16-06-principles-of-automatic-control-fall-2012/resources/mit16_06f12_lecture_14/ (loaded; PDF loaded)
- 16.06 Lecture 16 "Frequency response design" (PDF: Frequency Response of Lead Compensator); Lecture 24 "Bode compensation"; Lecture 25 "Lead compensation" (PDF: the maximum phase lead; a Kp constraint); Lecture 26 "Lag compensation (cont.)" (PDF: Mr predicted by PM; Kv constraint)
  - .../mit16_06f12_lecture_16/, .../mit16_06f12_lecture_25/, .../mit16_06f12_lecture_26/ (loaded; PDFs loaded); .../mit16_06f12_lecture_24/ (listed on the loaded Lecture Notes page)
- 16.30 Topic 4 "Control design using Bode plots": Bode's Gain Phase Relationship; Performance Issues; Frequency Domain Design; Lead Mechanics; Design Example; Lag Mechanics; Lag Compensation; Lag Example; Summary
  - https://ocw.mit.edu/courses/16-30-feedback-control-systems-fall-2010/resources/mit16_30f10_lec04/ (loaded; PDF loaded)
- Extract: φ_max = arcsin((α − 1)/(α + 1)) at the geometric mean of zero and pole; place the lead's maximum at the new crossover; lag raises low-frequency gain (Kv) while keeping its phase cost well below crossover; cross over with slope −1.

## 09-13 State-Space Models

Sources: `@mit-ocw-16-30-feedback-control-systems`, `@mit-ocw-6-011-signals-systems-and-inference`, `@mit-ocw-18-06-linear-algebra`, `@scipy-userguide`

- 16.30 Topic 5 "Introduction to state-space models": SS Introduction; Basic Definitions; Creating State-Space Models; Equilibrium Points; Linearization; Stability of LTI Systems
  - https://ocw.mit.edu/courses/16-30-feedback-control-systems-fall-2010/resources/mit16_30f10_lec05/ (loaded; PDF loaded)
- 16.30 Topic 6 "Developing state-space models based on transfer functions": TF's to State-Space Models; Modal Form; State-Space Models to TF's; SS to TF; State-Space Transformations
  - .../mit16_30f10_lec06/ (loaded; PDF loaded)
- 16.30 Topic 7 "State-space models: basic properties": Time Response; SS: Forced Solution; Eigenvalues and Eigenvectors; Jordan Form; Dynamic Interpretation; Diagonalization with Complex Roots; Topic 8 "System zeros and transfer function matrices"
  - .../mit16_30f10_lec07/, .../mit16_30f10_lec08/ (loaded; PDFs loaded)
- 6.011 Lecture 4 "State-Space Models"; Lecture 5 "State-Space Models, Equilibrium, Linearization"; Lecture 6 "Modal Solution of Undriven CT LTI State-Space Models"; Lecture 8 "Matrix Exponential, ZIR+ZSR, Transfer Function, Hidden Modes, Reaching Target States"
  - https://ocw.mit.edu/courses/6-011-signals-systems-and-inference-spring-2018/resources/mit6_011s18lec4/ (loaded; PDF loaded); .../mit6_011s18lec8/ (loaded; PDF loaded)
- 18.06 Lecture 23: Differential equations and exp(At); Lecture 28: Similar matrices and Jordan form
  - https://ocw.mit.edu/courses/18-06-linear-algebra-spring-2010/resources/lecture-28-similar-matrices-and-jordan-form/ (listed on the loaded Video Lectures page)
- SciPy API reference "Signal processing (scipy.signal)": "LTI representations" (`tf2ss`, `ss2tf`) and "Continuous-time linear systems" (`lti`, `StateSpace`); "Linear algebra (scipy.linalg)": Matrix Functions (`expm`) — oracles only
  - https://docs.scipy.org/doc/scipy/reference/signal.html, https://docs.scipy.org/doc/scipy/reference/linalg.html (loaded)
- Extract: companion (controller) and modal forms; G(s) = C(sI − A)⁻¹B + D; a similarity transform changes coordinates but not the TF or the eigenvalues; x(t) = e^{At}x₀ + ∫e^{A(t−τ)}Bu dτ.

## 09-14 Controllability & Observability

Sources: `@mit-ocw-16-30-feedback-control-systems`, `@mit-ocw-6-011-signals-systems-and-inference`, `@mit-ocw-18-06-linear-algebra`

- 16.30 Topic 9 "State-space model features": State-Space Model Features; Observability (unobservable states; the observability-matrix rank test)
  - https://ocw.mit.edu/courses/16-30-feedback-control-systems-fall-2010/resources/mit16_30f10_lec09/ (loaded; PDF loaded)
- 16.30 Topic 10 "Controllability": Controllability (Cayley–Hamilton argument); Further Examples; Modal Tests; Cancelation; Connection to Residue; Weaker Conditions
  - .../mit16_30f10_lec10/ (loaded; PDF loaded)
- 6.011 Lecture 7 "Full Modal Solution, Asymptotic Stability, Reachability and Observability"; Lecture 8 "... Hidden Modes, Reaching Target States"; Lecture 9 "Observers for State Estimation" (hidden modes of series, parallel and feedback connections)
  - https://ocw.mit.edu/courses/6-011-signals-systems-and-inference-spring-2018/resources/mit6_011s18lec7/ (loaded; PDF loaded)
- 18.06 Lecture 29: Singular value decomposition (numerical rank)
  - https://ocw.mit.edu/courses/18-06-linear-algebra-spring-2010/resources/lecture-29-singular-value-decomposition/ (listed on the loaded Video Lectures page)
- Extract: rank of [B AB … Aⁿ⁻¹B] and of the observability matrix; modal tests (an uncontrollable mode has a left eigenvector orthogonal to B — the PBH idea); pole-zero cancellation as a lost mode; numerical rank by SVD rather than determinant.
- Partial GAP: no loaded source names the PBH test as such, defines the controllability and observability Gramians, solves them as Lyapunov equations, or derives minimum-energy steering. Those parts, and the reaction-wheel satellite with a failed wheel or sensor, are original.

## 09-15 Pole Placement & Observers

Sources: `@mit-ocw-16-30-feedback-control-systems`, `@mit-ocw-6-011-signals-systems-and-inference`, `@scipy-userguide`

- 16.30 Topic 11 "Full-state feedback control": Full-state Feedback Controller; Ackermann's Formula; Origins of Ackermann's Formula; Reference Inputs; Pole Placement Examples; FSFB Summary
  - https://ocw.mit.edu/courses/16-30-feedback-control-systems-fall-2010/resources/mit16_30f10_lec11/ (loaded; PDF loaded)
- 16.30 Topic 14 "Open-loop and closed-loop estimators": Estimators/Observers; Open-loop Estimator; Closed-loop Estimator; Regulator/Estimator Comparison; Estimation Gain Selection; Dual Design Approach; Where Put Estimator Poles?
  - .../mit16_30f10_lec14/ (loaded; PDF loaded)
- 16.30 Topic 15 "Combined estimators and regulators": Combined Estimators and Regulators (the separation of closed-loop poles); The Compensator; Satellite Control Design; DOFB Summary
  - .../mit16_30f10_lec15/ (loaded; PDF loaded)
- 6.011 Lecture 9 "Observers for State Estimation"; Lecture 10 "Observers, State Feedback"; Lecture 11 "State Feedback, Observer-Based Feedback"
  - https://ocw.mit.edu/courses/6-011-signals-systems-and-inference-spring-2018/resources/mit6_011s18lec10/ (loaded; PDF loaded)
- SciPy API reference: `place_poles` (oracle)
- Extract: matching char. polynomials; Ackermann's formula K = [0 … 1]𝒞⁻¹φ(A); observer gain by duality (Aᵀ, Cᵀ); error dynamics ė = (A − LC)e; separation principle; place observer poles 2–5× faster than the regulator.

## 09-16 LQR

Sources: `@mit-ocw-16-30-feedback-control-systems`, `@mit-ocw-16-323-principles-of-optimal-control`, `@mit-ocw-18-06-linear-algebra`, `@scipy-userguide`

- 16.30 Topic 12 "Pole placement approach": Pole Placement Approach; Linear Quadratic Regulator; LQR Weight Matrix Selection (Bryson's Rules); Regulator Summary; the symmetric root locus; Example: B747 (B747 Lateral LQR)
  - https://ocw.mit.edu/courses/16-30-feedback-control-systems-fall-2010/resources/mit16_30f10_lec12/ (loaded; PDF loaded)
- 16.30 Topic 18 "Deterministic linear quadratic regulator (LQR)": Linear Quadratic Regulator (LQR); Constrained Optimization; LQR Optimization (the Hamiltonian matrix and the algebraic Riccati equation); LQR Stability Margins (gain margin (1/2, ∞), phase margin at least ±60°)
  - .../mit16_30f10_lec18/ (loaded; PDF loaded)
- 16.323 Lecture 4 "HJB equation: ... continuous LQR" (PDF heading "HJB Equation": the differential Riccati equation and its steady state); Lecture 3 "Dynamic programming: ... discrete LQR"
  - https://ocw.mit.edu/courses/16-323-principles-of-optimal-control-spring-2008/resources/lec4/ (loaded; PDF loaded)
- 18.06 Lecture 25: Symmetric matrices and positive definiteness (Q ⪰ 0, R ≻ 0, P ≻ 0)
- SciPy API reference "Linear algebra (scipy.linalg)": Matrix Equation Solvers (`solve_continuous_are`) — oracle
- Extract: the cost ∫(xᵀQx + uᵀRu)dt; the Hamiltonian matrix and the ARE; K = R⁻¹BᵀP; Bryson's rule for Q and R; LQR's guaranteed margins.
- Partial GAP: the loaded notes form the Hamiltonian matrix but do not derive the ARE solution from its stable invariant subspace (eigenvectors) step by step, and no verified source covers Newton–Kleinman iteration; both solvers are original.

## 09-17 Digital Control

Sources: `@mit-ocw-16-06-principles-of-automatic-control`, `@mit-ocw-16-30-feedback-control-systems`, `@mit-ocw-6-003-signals-and-systems`, `@nasa-tn-d-7843-f8-digital-fly-by-wire`, `@scipy-userguide`

- 16.06 Lecture 29 "Digital control, the z-transform" (PDF: Discrete Time Control Time line; Z-Transform)
  - https://ocw.mit.edu/courses/16-06-principles-of-automatic-control-fall-2012/resources/mit16_06f12_lecture_29/ (loaded; PDF loaded)
- 16.06 Lecture 30 "The z‐transform, design by emulation, the Tustin transform" (PDF heading "Z-transform Inversion": Design by Emulation; Tustin's approximation; matched pole-zero (MPZ))
  - .../mit16_06f12_lecture_30/ (loaded; PDF loaded)
- 16.06 Lecture 31 "Compensator design examples, time delay of ZOH" (the ZOH adds about half a sample of delay and phase lag at crossover); Lecture 32 "Discrete design"; Lecture 34 "Design examples, pre‐warping, direct design"
  - .../mit16_06f12_lecture_31/, .../mit16_06f12_lecture_32/, .../mit16_06f12_lecture_34/ (loaded; PDFs loaded)
- 16.30 Topic 20 "Digital control basics": Digital Control; Digital Control Mechanics; Sample and Hold Analysis ("Effective Delay"); "So why not just make Ts small?"; Control Law; Computer Code Layout
  - https://ocw.mit.edu/courses/16-30-feedback-control-systems-fall-2010/resources/mit16_30f10_lec20/ (loaded; PDF loaded)
- 6.003 Lecture 7 "Discrete Approximation of Continuous-Time Systems": Forward Euler: Mapping CT poles to DT poles; Backward Euler: Mapping CT poles to DT poles; Trapezoidal Rule: Mapping CT poles to DT poles; Lecture 5 "Z Transform"
  - https://ocw.mit.edu/courses/6-003-signals-and-systems-fall-2011/resources/mit6_003f11_lec07/ (loaded; PDF loaded)
- Aerospace context (supplementary): Deets, "Design and Development Experience With a Digital Fly-By-Wire Control System in an F-8C Airplane" (TN D-7843 paper 2): FLIGHT CONTROL SOFTWARE → Control Law Modes (multirate sampling, 30 ms major and 90 ms minor cycles); DIGITAL FLY-BY-WIRE DESIGN → Design Synthesis and Analysis (an acceptable range of 25 to 50 samples per second; folding effects on structural modes). Szalai (paper 6): CONTROL SYSTEM PERFORMANCE → Direct Mode (A/D quantization; effective trim quantization due to sample rate); Stability Augmentation System Mode (z-plane root locus of the pitch SAS)
  - https://ntrs.nasa.gov/citations/19750010175 (loaded; PDF loaded); https://ntrs.nasa.gov/citations/19750010179 (loaded; PDF loaded)
- SciPy API reference "Discrete-time linear systems" / `cont2discrete` (oracle)
- Extract: ZOH discretization Φ = e^{AT}, Γ = ∫e^{Aσ}dσ B (by the augmented matrix exponential); Tustin and matched pole-zero; z-plane stability |z| < 1; the ZOH and computation delay as phase lag eating phase margin; how a real flight-control computer chose its sample rate.
- Partial GAP: the augmented-matrix (Van Loan) identity for ZOH is not derived in the loaded notes; it is original here (it builds on P02.05).

## 09-18 Robustness & Uncertainty ★

Sources: `@mit-ocw-16-30-feedback-control-systems`, `@mit-ocw-16-323-principles-of-optimal-control`, `@mit-ocw-16-06-principles-of-automatic-control`

- 16.30 Topic 24 "Closed-loop system analysis": SISO Performance Objectives; SISO Design Approaches (working with S and T); Bounded Gain (the SISO Bounded Gain Theorem: sup|G(jω)| < γ tested through the eigenvalues of a Hamiltonian matrix); Typical Application; Issues; Appendix: Sketch of Proof
  - https://ocw.mit.edu/courses/16-30-feedback-control-systems-fall-2010/resources/mit16_30f10_lec24/ (loaded; PDF loaded)
- 16.323 Lecture 15 "Signals and system norms: H∞ synthesis" (PDF headings "Signals and System Norms" and "H∞ Synthesis": H∞ norm for the system; Computation — the H∞ Hamiltonian matrix and bisection on γ)
  - https://ocw.mit.edu/courses/16-323-principles-of-optimal-control-spring-2008/resources/lec15/ (loaded; PDF loaded)
- 16.323 Lecture 13 "LQG robustness" (high sensitivity when L(jω) comes close to −1)
  - .../lec13/ (loaded; PDF loaded)
- 16.06 Lecture 9: S and T and the tracking/noise trade-off; Lecture 23: Stability Margins
- Extract: Ms = max|S| as the inverse of the shortest distance from L(jω) to −1; H∞ norm by frequency sweep and by Hamiltonian bisection; robust stability as a bound on a weighted T.
- Partial GAP: no verified free source covers the disk margin and its bounds on gain and phase margin from Ms, the multiplicative-uncertainty model with the small-gain test ‖WT‖∞ < 1 as a derivation, or the Bode sensitivity integral (waterbed effect). These are original. Candidate: Åström & Murray *Feedback Systems* Ch. 13 "Robust Performance", once verifiable.

## 09-19 Nonlinear Control Basics

Sources: `@mit-ocw-16-30-feedback-control-systems`, `@mit-ocw-6-243j-dynamics-of-nonlinear-systems`, `@scipy-userguide`

- 16.30 Topic 21 "Systems with nonlinear functions": NL Example; Describing Function Analysis; Saturation Nonlinearity; Odd Nonlinearities with Memory; Limit Cycle Analysis; DF Analysis Example (saturation; relay with hysteresis); Limit Cycle Stability
  - https://ocw.mit.edu/courses/16-30-feedback-control-systems-fall-2010/resources/mit16_30f10_lec21/ (loaded; PDF loaded)
- 16.30 Topic 22 "Analysis of nonlinear systems": Lyapunov Stability Analysis; Lyapunov Stability Theorem; Example 1: Pendulum; Example 2: Linear System (AᵀP + PA = −Q, "called a Lyapunov equation"); Example 4: Local Region; Example 5: Saturation; Invariance Principle
  - .../mit16_30f10_lec22/ (loaded; PDF loaded)
- 16.30 Topic 5: Linearization; Stability of LTI Systems (Lyapunov's indirect method)
- 6.243J (supplementary, graduate) lecture titles: 5 "Lyapunov Functions and Storage Functions"; 7 "Finding Lyapunov Functions"; 8 "Local Behavior at Equilibria"; 13 "Feedback Linearization"
  - https://ocw.mit.edu/courses/6-243j-dynamics-of-nonlinear-systems-fall-2003/resources/lec13_6243_2003/ (listed on the loaded Lecture Notes page; PDF not loaded)
- SciPy API reference "Linear algebra (scipy.linalg)": Matrix Equation Solvers (`solve_continuous_lyapunov`) — oracle
- Extract: V(x) > 0, V̇ < 0; the Lyapunov equation for linear systems; a sublevel set of V as a region-of-attraction estimate; describing function N(A) of saturation and relay; limit cycle where G(jω) = −1/N(A).
- Partial GAP: the 6.243J feedback-linearization lecture was not opened and is graduate level; the input–output feedback linearization of an inverted pendulum (or nonlinear actuator) at this phase's level is original.

## 09-20 Model Predictive Control ★

Sources: `@mit-ocw-16-323-principles-of-optimal-control`, `@rawlings-mayne-diehl-mpc-2nd-edition`

- 16.323 Lecture 16 "Model predictive control" (PDF heading "Model Predictive Control": stacking predictions over a horizon into a standard QP; full-state assumption; stability with constraints; a helicopter example)
  - https://ocw.mit.edu/courses/16-323-principles-of-optimal-control-spring-2008/resources/lec16/ (loaded; PDF loaded)
- 16.323 Lecture 3 "Dynamic programming: principle of optimality, dynamic programming, discrete LQR" (the backward Riccati recursion that the unconstrained MPC must match)
  - .../lec3/ (loaded; PDF loaded)
- Rawlings, Mayne & Diehl (supplementary, all rights reserved — link only): Ch. 1 "Getting Started with Model Predictive Control": 1.2.4 Discrete Time Models; 1.2.5 Constraints; 1.3 Introductory MPC Regulator (1.3.1 Linear Quadratic Problem; 1.3.2 Optimizing Multistage Functions; 1.3.3 Dynamic Programming Solution; 1.3.4 The Infinite Horizon LQ Problem; 1.3.5 Controllability; 1.3.6 Convergence of the Linear Quadratic Regulator); Ch. 2 "Model Predictive Control—Regulation": 2.2 Model Predictive Control; 2.3 Dynamic Programming Solution; 2.4 Stability; 2.5 Examples of MPC; 2.6 Is a Terminal Constraint Set Xf Necessary?; Ch. 8 "Numerical Optimal Control": 8.7 Newton-Type Optimization with Inequalities; 8.8 Structure in Discrete Time Optimal Control
  - https://sites.engineering.ucsb.edu/~jbraw/mpc/ (loaded); .../MPC-book-2nd-edition-6th-printing.pdf (loaded; contents and copyright pages read)
- Extract: condensed prediction matrices and the QP ½UᵀHU + fᵀU; the unconstrained solution equals finite-horizon LQR from the Riccati recursion, and with terminal cost P from the ARE equals infinite-horizon LQR; receding horizon; input and rate box constraints.
- Partial GAP: no verified free source at this level walks through a from-scratch projected-gradient or small active-set QP solver for box constraints; that solver and the altitude-hold MPC example are original (P03.15 supplies the KKT background).

## 09-21 Home Lab: Servo Feedback Control

Sources: `@sparkfun-learn-tutorials`, `@micropython-docs`, `@mit-ocw-6-003-signals-and-systems`, `@mit-ocw-16-06-principles-of-automatic-control`

- SparkFun "Hobby Servo Tutorial": Servo Motor Background (Electrical Connection; Control signal; Powering Servos; Show Me The Guts — the internal potentiometer feedback loop); A Handful of Distinctions (Range Constraints; Position vs. Continuous Rotation; Analog vs. Digital); Troubleshooting
  - https://learn.sparkfun.com/tutorials/hobby-servo-tutorial/all (loaded)
- SparkFun "Pulse Width Modulation": Duty Cycle; Examples
  - https://learn.sparkfun.com/tutorials/pulse-width-modulation/all (loaded)
- SparkFun "TB6612FNG Hookup Guide" (small dual motor driver): Board Overview; Hardware Setup; Library and Example Code (Arduino); "Motors and Selecting the Right One": DC Brush Motors - The Classic (The Anatomy of a Brush Motor; Theory of Operation)
  - https://learn.sparkfun.com/tutorials/tb6612fng-hookup-guide/all (loaded); https://learn.sparkfun.com/tutorials/motors-and-selecting-the-right-one/all (loaded)
- MicroPython "class PWM – pulse width modulation": Constructors; Methods; Limitations of PWM; "class Timer – control hardware timers": Timer Types; Constructors; Methods (a fixed-rate control-loop tick); "class ADC – analog to digital conversion"; "Quick reference for the RP2": Timers; PWM (pulse width modulation); ADC (analog to digital conversion)
  - https://docs.micropython.org/en/latest/library/machine.PWM.html (loaded); .../library/machine.Timer.html (loaded); .../library/machine.ADC.html (loaded); https://docs.micropython.org/en/latest/rp2/quickref.html (loaded)
- Theory: 6.003 Lecture 12 "Continuous-Time (CT) Feedback and Control, Part 1": Motor Controller; Dominant Pole; 16.06 Lecture 6 "Time-domain specifications"
- Extract: how a hobby servo's own pot-and-amplifier loop works (and why the lab bypasses or models it); PWM duty cycle as the actuator command; a timer interrupt or fixed-period loop for a constant sample time; ADC reading of the wiper.
- Partial GAP: no verified free source gives a step-test system-identification procedure (fitting K and τ, or K/(s(τs + 1)), from measured data) or a home closed-loop servo experiment with predictions written down before the run. The identification method, the P/PI design to a stated target, the procedure, safety limits and the synthetic dataset are original (they follow p04-electronics-bench practice). The SparkFun servo and motor-driver examples use Arduino, not MicroPython.


## Optional textbooks for deeper study

Paid; never required. Not recorded in `sources.yaml`.

- G. F. Franklin, J. D. Powell and A. Emami-Naeini, *Feedback Control of Dynamic Systems* (Pearson) — the text assigned by 16.06 and recommended by 16.30; covers 09-03 … 09-17 and the Ziegler–Nichols and anti-windup gaps in 09-11.
- K. J. Åström and R. M. Murray, *Feedback Systems: An Introduction for Scientists and Engineers*, 2nd ed. (Princeton University Press) — PID, robustness (disk margin, Bode integral) and loop shaping for 09-11, 09-12, 09-18. The authors have hosted a free PDF, but that site could not be verified on 2026-10-17, so only the print edition is listed here.
- A. V. Oppenheim and A. S. Willsky, *Signals and Systems* (Pearson) — the classic text behind 6.003 for 09-01 … 09-07.
- A. V. Oppenheim and R. W. Schafer, *Discrete-Time Signal Processing* (Pearson) — the text behind 6.341 for 09-06 and 09-07.
- G. F. Franklin, J. D. Powell and M. Workman, *Digital Control of Dynamic Systems* (Ellis-Kagle) — ZOH, emulation and direct digital design for 09-17.
- S. Skogestad and I. Postlethwaite, *Multivariable Feedback Control: Analysis and Design* (Wiley) — sensitivity peaks, uncertainty and H∞ for 09-18.
- H. K. Khalil, *Nonlinear Systems* (Pearson) — Lyapunov theory and feedback linearization beyond 09-19.
- B. L. Stevens, F. L. Lewis and E. N. Johnson, *Aircraft Control and Simulation* (Wiley) — aircraft pitch, altitude and attitude loops used as this phase's examples.
