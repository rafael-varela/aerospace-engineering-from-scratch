# P10.06  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. Air is ISA sea level (ν = 1.461 × 10⁻⁵ m²/s) unless stated. Inside the code the chord and the free-stream speed are 1, so ν = 1/Re. Say which numbers are computed and which are measured; the only measurements here are Ladson's (NASA TM-4074) in `curriculum/phase-10/data/ladson_naca0012_M015.csv`.


## 1. (Warm-up: the boundary layer on the trainer's wing)

Treat the first part of the trainer's wing as a flat plate at 55 m/s.

(a) At x = 0.3 m compute Re_x, and θ, δ* and c_f from Blasius (use your `falkner_skan(0)` integrals) and from Thwaites (λ = 0). How large is each in millimetres? What is the percentage difference between the two methods?
(b) By Michel's criterion applied to the Blasius layer, at what Re_x and what distance does a flat-plate boundary layer at 55 m/s become turbulent?
(c) The real wing's boundary layer starts at the stagnation point, accelerates, then decelerates. Which way does each part move transition compared with (b)?

Answer form: Re_x and six numbers with a percentage; Re_x and a distance; two sentences.

## 2. (Derivation: Thwaites' method and the Squire–Young formula)

(a) Starting from the von Kármán momentum integral equation dθ/ds + (2 + H)(θ/u_e) du_e/ds = c_f/2, multiply by u_eθ/ν and write it in terms of λ = θ²u_e′/ν. Show that Thwaites' empirical linearisation, u_e d(θ²)/ds/ν ≈ 0.45 − 6λ, integrates to θ²u_e⁶ = 0.45ν∫₀^s u_e⁵ ds.
(b) For stagnation-point flow u_e = Cs, show that θ is constant and equal to √(0.075ν/C), and that λ = 0.075 everywhere. Compare with the Hiemenz value 0.2923√(ν/C) from your `falkner_skan(1)`.
(c) The Squire–Young formula c_d = 2(θ_TE/c)(u_e,TE/V)^{(H_TE + 5)/2} extrapolates the trailing-edge momentum thickness to the far wake. Show that it reduces to c_d = 2θ/c when u_e,TE = V, and explain physically why the exponent makes the drag smaller when the trailing-edge velocity is below the free stream.

Answer form: three derivations, one comparison.

## 3. (Unfamiliar scenario: the small UAV)

The reference small UAV has a 0.20 m chord and cruises at 12 m/s. Take its section as a NACA 2412.

(a) Compute the chord Reynolds number.
(b) Run `section_drag` at α = 0° and 4°. Where and by which mechanism does each surface go turbulent? What c_d does the method predict?
(c) The method assumes that a laminar separation is followed at once by a turbulent attached layer. At this Reynolds number that is a laminar separation bubble (`P10.14`), and it may not reattach. How far would you trust the c_d of (b), and in which direction is it likely to be wrong?

Answer form: a number; two lines of results; two sentences.

## 4. (Numerical, uses your code: NACA 0012 against Ladson)

(a) Read the Ladson file with `read_ladson` and extract the zero-lift drag of Table I (free transition) with `ladson_zero_lift_drag`. For each Reynolds number, predict c_d with 200 panels and tabulate the measured value, the prediction and the percentage difference.
(b) Do the same for Table VII (the NACA wraparound grit, fixed transition) against your prediction with transition forced at 5 % chord. Is the agreement better or worse, and why might a grit strip that wraps round the leading edge add drag that a "trip at 5 %" model does not?
(c) Plot c_d against Re for both, with the measured points clearly marked as measured. Describe the trends in one sentence each.

Answer form: two tables; a sentence; a plot and two sentences.

## 5. (Numerical, uses your code: the trainer's section polar)

At the trainer's cruise Reynolds number (4.52 × 10⁶):

(a) Tabulate c_l (inviscid) and c_d for α = 0°, 2°, 4°, 6° and 8°, with the transition positions on both surfaces and any turbulent separation.
(b) Run one displacement-thickness iteration at α = 0°, 2°, 4° and 6°. By how much does it reduce c_l? Why does the reduction grow with α?
(c) Estimate the wing's profile drag at cruise (55 m/s, S = 12 m²) clean and with transition fixed at 5 % chord (insects, rain or rivets), and state the difference in kilowatts of power required.

Answer form: a table; four percentages and a sentence; two forces and a power.

## 6. (Judgment: why does the prediction fall below the data at high Re?)

At Re = 1.19 × 10⁷ the model predicts c_d = 0.00567 with transition at 24.5 % chord, while Ladson measured 0.00635.

(a) Use `section_drag(..., trip_x=…)` with trips at 5, 10 and 20 % chord to find roughly where transition must be for the model to match the measurement.
(b) Give two physical reasons why transition in the tunnel could be further forward than Michel's criterion predicts, and one reason the model could be low even with the right transition point.
(c) Which single improvement to the method would you make first, and what new validation data would you need?

Answer form: a transition position; three reasons; a choice with a reason.
