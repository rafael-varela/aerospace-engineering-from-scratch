# P10.22  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer (CL_α per degree unless stated) and say for each number whether it is computed, measured (in a cited report) or an uncertainty estimate. The measured values come from NACA TN-1669, TN-1671 and TN-1422 as transcribed in `curriculum/phase-10/data/finite_wing_cla.csv`; most were digitized from figures.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: which kind of evidence is it?)

Classify each item as code verification, solution verification, or validation, and say in one line what it can and cannot show:

(a) the span load of a symmetric wing is symmetric to 10⁻¹⁵;
(b) CL_α changes from 4.108 to 4.056 per radian when the lattice is refined from 4 × 24 to 8 × 48;
(c) the VLM predicts 0.0682 per degree for TN-1669 wing 7, measured 0.0689;
(d) a 1000-aspect-ratio unswept wing gives 6.263 per radian;
(e) the GCI of the TN-1422 prediction is 1.5 %;
(f) TN-1669 and TN-1671 report 0.0463 and 0.0449 for the same planform.

Answer form: six classifications with one line each.

## 2. (Derivation: the tools of solution verification and the section correction)

(a) For three solutions f₁ (fine), f₂, f₃ on grids refined by a constant ratio r, with error C hᵖ, derive the observed order p = ln((f₃ − f₂)/(f₂ − f₁))/ln r and the Richardson estimate f₁ + (f₁ − f₂)/(rᵖ − 1). State the condition the three values must satisfy for this to make sense.
(b) Roache's GCI on the fine grid is Fs|(f₂ − f₁)/f₁|/(rᵖ − 1) with Fs = 1.25 for a three-grid study. What is it intended to bound, and why does the lesson take u_num = GCI/2?
(c) In lifting-line theory the circulation satisfies Γ(y) = ½ V c(y) a₀ [α − α_i(y)]. Show that a wing with chord c(y) and section slope a₀ = κ 2π has the same Γ(y) as a wing with chord κc(y) and slope 2π, and hence that CL_α(AR, a₀) = κ CL_α,2π(AR/κ). Why is this only an approximation for a vortex lattice?

Answer form: two derivations with the condition; a sentence and a reason; a derivation and a sentence.

## 3. (Unfamiliar scenario: a prediction with an uncertainty for the reference trainer)

Before any test, the designer of the reference trainer (rectangular wing, AR 8.33, NACA 2412) asks for its lift-curve slope with an uncertainty.

(a) Using NACA 0012 a₀ from Ladson at Re = 2.0 × 10⁶ as a proxy for the section, compute CL_α on the three lattices of the lesson and give the fine-lattice value, the observed order, the Richardson value and u_num.
(b) Estimate u_input from the difference between Ladson's slopes at 2.0 and 5.97 × 10⁶ and the sensitivity of CL_α to κ.
(c) State the prediction with an expanded uncertainty (k = 2). Which term dominates? What would you measure or compute to shrink it?
(d) Using the validation's mean signed error and its spread, would you add a model-error term? Argue either way.

Answer form: four numbers; a number; a prediction with U and a sentence; a short argument.

## 4. (Numerical, uses your code: how the verdict depends on what you assumed)

The fraction of planforms "within uncertainty" depends on the uncertainty you assumed.

(a) Recompute the validation table with the run-to-run scatter set to 1.5 %, 3 % and 6 % (and also (3 + s/10) %). Report the fraction within k = 2 each time.
(b) Recompute it with κ = 1 (no section correction). Compare the mean signed error and the mean absolute error with the corrected case.
(c) A colleague proposes reporting "93 % of wings within uncertainty". What would make that honest, and what would make it misleading?

Answer form: four fractions; four numbers; a paragraph.

## 5. (Numerical, uses your code: does the GCI cover the error?)

For the planforms (AR 1.34, unswept), (AR 9, taper 0.4) and (AR 2.61, 60° sweep), with κ from Ladson at 2.0 × 10⁶:

(a) Compute CL_α on 2 × 12, 4 × 24, 8 × 48 and 16 × 96 lattices (chordwise × full span). Report the observed order from the first three, the GCI and the Richardson value.
(b) Treat the 16 × 96 value as "nearly exact" and compare the actual error of the 8 × 48 value with the GCI. Does the GCI bound it?
(c) Why does a uniform lattice converge at only about first order here?

Answer form: a table; three comparisons; a sentence.

## 6. (Judgment: where the linear model fails)

(a) The AR 1.34 unswept wing is the worst prediction (VLM 16 % high). The data README quotes TN-1669: for the low-aspect-ratio wings "the lift-curve slope increases" with CL. Explain what physical effect does this, why a linear lattice cannot produce it, and why it also makes the measured "slope" depend on the CL range used.
(b) The TN-1671 tapered wings (taper 0.25 and 0.5) are predicted 6–8 % high, while TN-1671 reports that taper had little effect below CL = 0.6. Name two possible causes, one in the model and one in the data.
(c) TN-1422's wings use NACA 64-210 and 65-210 sections, not NACA 0012. How does using Ladson's NACA 0012 slope as a proxy enter your uncertainty budget, and what free source would give a better a₀?
(d) Write the two sentences of a validation conclusion that a reviewer could not fault.

Answer form: four short paragraphs.
