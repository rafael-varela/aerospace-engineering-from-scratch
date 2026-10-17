# P04.14  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical. No problem here involves a real source or detector; every input is given.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: photons and doses)

(a) Compute the vacuum wavelength and photon energy of the hydrogen line n_u → 2 with n_u = 4 + (s mod 4).
(b) A green lidar (532 nm) fires (5 + s) mJ pulses. How many photons per pulse? What range corresponds to a round trip of s µs?
(c) What is the ground-state energy of He⁺ (one electron, Z = 2), ignoring the reduced mass?
(d) A fluence of 10¹⁰ protons/cm² crosses a thin silicon layer; take the mass stopping power as 34.6 MeV cm²/g (an illustrative value for protons of about 10 MeV). What dose is deposited, in Gy and in krad(Si)?

Answer form: six numbers with units.

## 2. (Derivation: the chain, its peak and its equilibrium)

(a) Solve N₂′ = λ₁N₁(0)e^{−λ₁t} − λ₂N₂ with N₂(0) = 0 by an integrating factor and obtain the Bateman solution (5).
(b) Show that the daughter activity peaks at t_max = ln(λ₂/λ₁)/(λ₂ − λ₁).
(c) Show that for λ₁ < λ₂ the ratio of daughter to parent activity tends to λ₂/(λ₂ − λ₁), and that it tends to 1 when λ₁ ≪ λ₂. Evaluate t_max and the limiting ratio for Sr-90 (T½ = 28.8 y) → Y-90 (T½ = 64.0 h).
(d) Derive the specific power λN_AQ/M and evaluate it for Pu-238 (T½ = 87.7 y, M = 238.0496 g/mol, Q from part 2.2 of the lesson).

Answer form: three derivations, two numbers, one number in W/g.

## 3. (Unfamiliar scenario: americium instead of plutonium)

Am-241 (T½ = 432.6 y, atomic mass 241.056 829 u) alpha-decays to Np-237 (237.048 173 u).

(a) Compute Q and the specific power in W/g.
(b) For a mission of (10 + s) years requiring 110 W electrical at the end, with 6 % efficiency and no degradation, compute the Pu-238 and the Am-241 inventories.
(c) Above what mission length would the two inventories be equal? (Set up the equation; solve it numerically.) Comment on whether that matters in practice.

Answer form: two numbers, two masses, one mission length, two sentences.

## 4. (Numerical, uses your code: a shielding trade)

A (MODEL, illustrative) dose-depth curve gives 20 krad/yr behind 1 mm and (2 + s/5) krad/yr behind 3 mm of aluminium.

(a) Fit D₀ and L with `fit_attenuation`. What is the half-value thickness in mm?
(b) For a 10-year mission with margin 2, compute the thickness and the mass of a 600 cm² aluminium shell for parts rated 30 krad and 100 krad.
(c) Perturb the 3 mm point by ±20 % and recompute the 30 krad thickness. Is the answer more or less sensitive than the input, and why?
(d) Name two physical reasons the real dose behind 5 mm could be higher than this model predicts.

Answer form: two numbers and a thickness; a 2 × 2 table; two thicknesses and one sentence; two sentences.

## 5. (Numerical with probability: upsets and scrubbing)

A 4 Gbit (4 × 2³⁰ bit) memory has σ = 10⁻¹⁴ (1 + s/10) cm²/bit in an (illustrative) upsetting flux of 1 cm⁻² s⁻¹. It is organised in 64-bit words protected by a code that corrects one error per word.

(a) Expected upsets per day; probability of no upset during one 95-minute orbit. Check the latter with 100 000 seeded Poisson draws.
(b) For a scrub interval T, the probability that some word collects two upsets is approximately N_words (r_w T)²/2 with r_w the per-word upset rate. Tabulate it for T = 1 s, 1 min, 1 h, 1 day.
(c) What scrub interval keeps that probability below 10⁻⁶? Explain why the approximation is good here.

Answer form: two numbers and a check, a 4-row table, one interval and one sentence.

## 6. (Judgment: matching the mitigation to the effect)

For each case, name the effect (TID, displacement damage, SEU, SEL, or charging) and the mitigation you would choose (shielding, a harder part, error-correcting code and scrubbing, current limiting with power cycling, design margin, or "accept it"), with one sentence of justification: (a) a commercial microcontroller in a 3-year LEO mission accumulating 4 krad behind 2 mm of aluminium, rated 10 krad; (b) a solar array losing output over 15 years in a proton-rich orbit; (c) a camera's frame buffer showing occasional single-pixel garbage; (d) a processor drawing 1 A suddenly and staying there; (e) an RTG-powered probe at Jupiter whose 300 krad part sits behind a shield chosen with the lesson's exponential model fitted between 1 and 3 mm, at 12 mm; (f) arcing on an exposed cable harness in geostationary orbit.

Answer form: six labelled lines.
