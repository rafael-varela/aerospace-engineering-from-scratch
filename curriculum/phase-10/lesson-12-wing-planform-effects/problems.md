# P10.12  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. Twist in degrees in answers, radians in the code. For every optimum, report whether the optimiser succeeded and which constraints are active.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: taper and the stall)

A flat, untwisted wing of span 10 m has AR = 6 + s/5.

(a) Compute e for λ = 1, 0.5 and 0.35 with `efficiency_vs_taper`.
(b) Compute the stall-onset station for each with `stall_onset` at CL = 1.0.
(c) Explain in two sentences why e and the stall station move in opposite directions as λ falls.

Answer form: three numbers; three stations; two sentences.

## 2. (Derivation: the optimiser's mathematics)

(a) Show that CDi is the quadratic form xᵀQx of (2) and that the CL constraint is linear in x.
(b) Derive the linear-twist optimum (3) with a Lagrange multiplier, without eliminating α first, and show it agrees with (3).
(c) Show that the root bending coefficient of the elliptic load is C_B = CL/(3π).
(d) Derive f(s) of (4) and show that G₃/G₁ ≥ −1/3 is the condition for a non-negative two-term load.

Answer form: four derivations.

## 3. (Unfamiliar scenario: a sailplane wing)

A 15 m sailplane wing has S = 11.25 m² (AR 20) and taper 0.4; it cruises between thermals at CL = 0.8 and thermals at CL = 1.2.

(a) Report e and the stall-onset station of the untwisted wing (stall station at CL = 1.2).
(b) Optimise the twist at CL = 0.8 with 5 stations: unconstrained, and with a 5 % tip margin at CL = 1.2 (washout only). Report e, the twists and the success flags.
(c) Find the best linear twist. Compare all four designs and recommend one for a pilot who spends half the flight thermalling.

Answer form: two numbers; two sets of numbers; a table and a recommendation.

## 4. (Numerical, uses your code: your own bending limit)

For the rectangular trainer at the cruise CL (0.3309), limit the root bending moment to (0.85 + s/200) times the elliptic value.

(a) Run `bending_trade` and report e, the achieved C_B ratio and the twist distribution.
(b) What does the extra induced drag cost at 55 m/s, in newtons? What does the lower root moment buy in spar-cap area, if the cap area is proportional to the root moment (`P05.10`)?
(c) Is the constraint active? How do you know from the optimiser's output alone?

Answer form: three items; two numbers; a sentence.

## 5. (Numerical, uses your code: Prandtl without the floor)

(a) Evaluate f(s) of (4) at s = 1.1, 1.2247, 1.5 and 2.0, and G₃/G₁ at each.
(b) Explain why f keeps falling past √(3/2) and why that is not a better wing.
(c) Run `prandtl_bell` with n_terms = 3 + (s mod 4) and confirm that the higher harmonics stay at zero. What would you have to change in the structural model for higher harmonics to matter?

Answer form: a table; an argument; a check and a sentence.

## 6. (Judgment: the trainer's planform)

(a) A colleague proposes taper 0.36 "because it is the most efficient planform". Using the lesson's numbers, write the review comment.
(b) The recommended design has about 3° of washout. At the approach CL, what happens to its tip margin compared with the untwisted wing, and is that good?
(c) Name two things outside the lesson's model that could still make the rectangular, washed-out wing stall at the tips first, and how you would check for them.

Answer form: three short arguments.
