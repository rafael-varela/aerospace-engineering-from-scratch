# P10.13  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which numbers are measured (Ladson's data), which are computed, and which are assumed. ISA sea level: ρ = 1.225 kg/m³, ν = 1.461 × 10⁻⁵ m²/s.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: friction and transition)

(a) The trainer's wing (chord 1.2 m) flies at V = 40 + s m/s. Compute Re and the mean skin friction for transition at x_tr = 0, 0.1, 0.3, 0.5 and 1.0 of the chord.
(b) What fraction of the fully turbulent friction is saved by keeping the flow laminar to mid-chord? Why do sailplane designers care so much?

Answer form: one Re and five C_f; a percentage and a sentence.

## 2. (Derivation: the build-up and the polar)

(a) Derive the mixed-flow formula (1) of the lesson. State the approximation it makes about the turbulent layer's origin.
(b) For C_D = C_D0 + kC_L², derive C_L* = √(C_D0/k), (L/D)max = 1/(2√(kC_D0)) and show that the induced drag equals the zero-lift drag there.
(c) Show that the speed for best L/D scales as (W/S)^{1/2} and as C_D0^{−1/4}. What happens to it if the gear is retracted?

Answer form: three derivations.

## 3. (Unfamiliar scenario: a sailplane)

A club sailplane: mass 450 kg, wing area 11.0 m², span 15.0 m, taper 0.4, 14 % thick laminar-flow sections with transition at 50 % chord, fuselage 6.5 m long and 0.6 m in diameter (assume an elliptic nose to 1.5 m and a linear taper to 0.06 m radius at the tail), tails 2.2 m² exposed with 0.6 m chord (NACA 0010, turbulent, Q = 1.05), retracted gear, 5 % miscellaneous. Fly at 30 m/s.

(a) Use `component` to build C_D0 (remember that `component` divides by the trainer's 12 m²; rescale to 11 m²).
(b) Estimate e with `lifting_line_e(AR, 0.4)`, the fuselage factor and the viscous term (4), then (L/D)max, its speed and the minimum-glide sink rate V/(L/D).
(c) Which single change to the trainer's build-up makes the biggest difference between the two aeroplanes?

Answer form: a component table and C_D0; four numbers; a sentence.

## 4. (Numerical, uses your code: the price of fixed gear)

(a) Recompute the trainer's C_D0, e and (L/D)max with retractable gear (gear drag area 0) and with the wheel fairings removed (0.25 m²).
(b) At the drag power of the baseline cruise (46.13 kW at 55 m/s), find the cruise speed of each variant (`brentq` on `cruise_drag`).
(c) A retractable gear adds about 20 kg and complexity. Is it worth it on a trainer? Argue with your numbers.

Answer form: three numbers per variant; two speeds; a paragraph.

## 5. (Data-based, uses your code: drag with lift)

(a) From Ladson's Table I, take the row nearest α = 4.2° in every Reynolds-number block. Tabulate α, cl, cd and the cd at α ≈ 0.
(b) Compute the effective transition location (`effective_transition`) for both angles. What happens to it with lift, and why?
(c) At Re = 1.89 × 10⁷ the inversion returns NaN at α ≈ 4°. Explain what that means.
(d) What does (a)–(c) say about using a zero-lift section estimate for a wing in cruise at C_L ≈ 0.3–0.5?

Answer form: a table; a table and two sentences; a sentence; a paragraph.

## 6. (Judgment: how good is C_D0 = 0.032?)

(a) List the inputs of the trainer's build-up in order of their effect on C_D0, with the change each makes over a plausible range (use the lesson's sensitivity numbers and your own runs).
(b) Give C_D0 with an uncertainty, and say how you combined the contributions.
(c) The free-transition discrepancy of the section estimate is up to 59 % at Re = 2 × 10⁶. Why does it not mean that the trainer's C_D0 is uncertain by 59 %?
(d) Propose one flight test that would measure C_D0 and e of the real aeroplane, and say what you would plot.

Answer form: a ranked list; a number with an uncertainty and a sentence; a paragraph; a test and a plot.
