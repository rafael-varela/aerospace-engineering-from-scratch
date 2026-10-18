# P10.11  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. Lift slopes per radian unless stated. Say which results are converged, and how you know.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: the classic coarse lattice)

An untapered wing of aspect ratio 5 with 45° of quarter-chord sweep is modelled with one chordwise row of 8 horseshoes across the span (4 per side), uniform spacing.

(a) Compute CL_α with `make_lattice(..., nx=1, ny=8, spacing="uniform")`.
(b) Refine to ny = 16, 32, 64 with uniform spacing and with Hough's inset. What does each converge to, and at what order?
(c) How far was the coarse answer from the converged one? What does that say about using a coarse lattice for design?

Answer form: one number; a table and two orders; a percentage and a sentence.

## 2. (Derivation: the rules behind the lattice)

(a) Derive the 1/4–3/4 rule from a single vortex on a 2-D flat plate (section 2.1).
(b) Two vortices on a flat plate (two chordwise panels, each with its own 1/4–3/4 points): set up the 2 × 2 system and show that the total Γ still gives cl = 2πα.
(c) Show that moving the moment reference by Δx changes Cm by CL cos α Δx/c_ref for the lattice's forces, and hence that (2) of the lesson defines a point independent of the reference.
(d) Derive (3), the Trefftz drag of a set of 2-D vortices, from D_i = ½ρ∫Γ w_T ds.

Answer form: four derivations.

## 3. (Unfamiliar scenario: forward sweep)

Compare an untapered AR 6 wing (span 10 m) swept −30° (forward), 0° and +30° (aft), nx = 4, ny = 40.

(a) Report CL_α, e and x_ac (in metres from the root leading edge and as a fraction of the MAC behind the MAC's leading edge).
(b) Where along the span is cl/CL largest for each? What does that mean for the stall?
(c) Forward-swept wings were long avoided for an aeroelastic reason. Find out what it is (`P05`-level reasoning is enough), and say why the lattice alone cannot predict it.

Answer form: a table; three stations and a sentence; a short argument.

## 4. (Numerical, uses your code: your own convergence study)

Use the trainer's wing (b = 10 m, c = 1.2 m), nx = 3, Hough spacing, α = 4°, and ny = n0, 2n0, 4n0 with n0 = 12 + 2 (s mod 3).

(a) Report the three CL values and the observed order.
(b) Extrapolate (Richardson) to ny → ∞ and give CL_α from the extrapolated CL. How far is your finest lattice from it?
(c) Repeat (a) with uniform spacing. Which lattice would you trust with ny = n0, and why?

Answer form: three numbers and an order; two numbers; three numbers, an order and a sentence.

## 5. (Data-based: the lattice against NACA TN-1669)

`curriculum/phase-10/data/finite_wing_cla.csv` holds measured lift slopes of untapered NACA 0012 wings of AR 1.34, 2.61 and 5.16 at sweeps of 0°, 45° and 60° (digitised from the report's figure, ±0.0005 per degree).

(a) Read the file (skip `#` lines) and compute the VLM slope (nx 6, ny 40) and the Helmbold–DATCOM value for all nine wings, per degree.
(b) Tabulate the differences from the measurements. Where does the VLM agree within the reading uncertainty, and where not?
(c) Explain the AR 1.34 unswept result. What physical effect is missing from the VLM?

Answer form: a table; a table and a sentence; an explanation.

## 6. (Judgment: lattice or line?)

(a) For the trainer, the lattice gives e = 0.970 and the lifting line 0.934. A performance engineer must choose one for the cruise polar. Which, and what difference does it make in newtons?
(b) A student's uniform-spacing lattice reports e = 1.07 for a rectangular wing. What is wrong, and what would you change?
(c) List three things the VLM cannot tell you about a wing that a designer must know, and where in this phase each comes from.

Answer form: a choice and a number; a diagnosis and a fix; three items with sources.
