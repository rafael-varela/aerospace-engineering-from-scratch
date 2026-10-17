# P10.14  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact (closed form), which are numerical, and which rest on an empirical correlation. Never present a model result below Re = 2 × 10⁶ as NACA data: the authentic file stops there.

Several problems use your personal seed **s**: the day of the month you were born (1–31). Sea-level air: ν = 1.461 × 10⁻⁵ m²/s.

## 1. (Warm-up: your own small flyer)

(a) Compute the chord Reynolds number of the reference UAV at 8, 12 and 18 m/s, and the all-laminar friction-drag floor cd = 2C_f,lam(1 + 2t/c + 60(t/c)⁴) for a 12 % section at each.
(b) A model glider flies at V = 5 + s/3 m/s with chord c = 0.05 + s/200 m. Compute its Re. Is natural transition likely before laminar separation on a NACA 0012 at cl = 0.445 (use the lesson's natural-transition Re)?

Answer form: three Re and three cd; one Re and a sentence.

## 2. (Derivation: Thwaites and bubbles)

(a) For Ue = Cs^m derive θ² = 0.45νs^{1−m}/(5m + 1) and λ = 0.45m/(5m + 1) from (1). Which m gives the stagnation-point value 0.075? Which m would give separation?
(b) Prove that the λ = −0.09 location does not depend on ν, and that θ at that location scales as Re^{−1/2}.
(c) Derive Re_min of (4) from the bubble-length correlation (3), and show that the bubble length in chords is N_ℓ/(Re·Ue_s/V).
(d) Explain why the trapezoid rule for ∫Ue⁵ds is badly wrong on the first interval after a stagnation point (Ue = as), and derive the exact integral for Ue linear on an interval.

Answer form: four derivations.

## 3. (Unfamiliar scenario: a rotor blade on Mars)

A small Mars helicopter blade section has chord 0.08 m and moves at 120 m/s at 0.75R. The Martian surface air has ρ ≈ 0.017 kg/m³ and μ ≈ 1.1 × 10⁻⁵ Pa·s (assumed round values).

(a) Compute ν and the chord Reynolds number. Compare with the table of §2.1.
(b) Compute the all-laminar friction floor for a 2 % thick plate at that Re.
(c) Using §2.4 and §2.5, which kind of section would you choose for this blade and why? What does the bubble model say about closure at this Re (estimate Re_min-type reasoning qualitatively)?

Answer form: two numbers; a number; a paragraph.

## 4. (Numerical, uses your code: slow flight)

At low speed the UAV's section works at cl = 0.80 (still at Re = 1.64 × 10⁵ for simplicity).

(a) Run `bubble_analysis` for NACA 0012, 2412 and 4402: α, peak speed, separation point, bubble length, closure, ΔCp_rec and Re_min.
(b) Explain the 4402's result with the shape of its Ue(s) (plot it). What does it say about thin plates away from their design angle?
(c) Run `natural_transition_re` for the 0012 and the 4402 at cl = 0.2, 0.445 and 0.8. Interpret the trend.

Answer form: a table; a paragraph; a table and a sentence.

## 5. (Numerical, uses your code and authentic data: what grit does)

(a) From the Ladson data, pair the free-transition and grit blocks at similar Re (2, 4, 6, 9, 12 × 10⁶) and compute the percentage change of cl_max caused by the wraparound grit strip.
(b) Compare the grit strip's cd near α = 0 with the fully turbulent friction model at each Re.
(c) A model-aircraft builder adds a turbulator strip to a wing at Re = 10⁵. Using (a), (b) and §2.4, explain why the turbulator can *help* there although the grit strip *hurts* at Re = 2 × 10⁶.

Answer form: five percentages; five pairs; a paragraph.

## 6. (Judgment: how far can you trust the bubble model?)

(a) Vary N_ℓ from 3 × 10⁴ to 7 × 10⁴ and report Re_min for the 0012 and the 4402 at the cruise cl. Which conclusions of the lesson survive the whole range?
(b) Name three physical effects that the chain panel → Thwaites → Michel → (3) leaves out, and say for each whether it would make the UAV's bubble shorter or longer.
(c) Thwaites is least accurate just before separation. Using the Howarth comparison, estimate the error in x_s on the 0012 and decide whether it matters for Re_min.
(d) Which single measurement on a real UAV wing would you make first, and what would it tell you?

Answer form: a table and a sentence; three effects with directions; an estimate and a sentence; a measurement.
