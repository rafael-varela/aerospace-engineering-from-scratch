# P10.09  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. SI units throughout; ISA sea level (ρ = 1.225 kg/m³, g = 9.81 m/s²). Say which results are exact and which numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: the downwash at the tail)

Replace the trainer's elliptic wing in cruise (W = 7357.5 N, V = 55 m/s) by a single horseshoe with the same lift, its legs at the centroids of the elliptic load's trailing vorticity, a distance b′ = πb/4 apart.

(a) Find Γ of the horseshoe and the downwash at the centre of its bound vortex.
(b) Find the downwash at the tail, 4.5 m behind on the centreline, with `horseshoe_velocity`, and far downstream. Express each as an angle.
(c) Why is the downwash at the tail larger than far downstream? Compare the far-field value with 2α_i of the elliptic wing.

Answer form: two numbers; two speeds and two angles; an explanation and a comparison.

## 2. (Derivation: from Biot–Savart to the polar)

(a) Derive equation (1) of the lesson for a straight segment from dV = (Γ/4π) dl × r/|r|³. Then show the symmetric-segment result (Γ/2πh) L/√(L² + h²).
(b) Derive the semi-infinite result (2) as a limit, and show it is half the infinite line in the plane through its start.
(c) Show that Γ = Γ0 sin θ gives a uniform w = Γ0/(2b), using ∫₀^π cos nφ/(cos φ − cos θ) dφ = π sin nθ/sin θ.
(d) Derive CL* = √(πeAR CD0), (L/D)max = ½√(πeAR/CD0) and CL_mp = √3 CL*.

Answer form: four derivations.

## 3. (Unfamiliar scenario: a wake encounter)

An airliner of mass 70 t and span 60 m flies at 70 m/s. Assume an elliptic load and that each tip vortex carries Γ0 and has rolled up at ±πb/8 from the centreline.

(a) Find Γ0 and the vortex separation.
(b) Find the downwash midway between the vortices (infinite lines).
(c) A light aircraft passes 5 m from one vortex core (ignore the other). What velocity does it meet, and how does that compare with its own 55 m/s and the downwash its own wing makes (problem 1)? Why are wake-separation minima a safety rule?

Answer form: two numbers; a speed; a speed, a comparison and an argument.

## 4. (Numerical, uses your code: a non-elliptic load)

Take Γ ∝ sin φ + a₃ sin 3φ with a₃ = −s/100 on the trainer's span.

(a) Compute e with `lift_and_drag` and `span_efficiency` on 20, 40, 80 strips, and compare with 1/(1 + 3a₃²).
(b) Report the ratio of near-field to far-field drag. Why is it 1?
(c) Is the load more or less tip-loaded than elliptic for your a₃? Sketch it.

Answer form: three numbers and a closed form; a ratio and a reason; a sentence and a sketch.

## 5. (Numerical, uses your code: your own polar)

Your aircraft is the trainer with CD0 = 0.025 + s/1000 and e = 0.75 + s/200 (still assumed values).

(a) Find (L/D)max, CL*, V_md, CL_mp and V_mp.
(b) At 55 m/s find the drag, the induced share and the power.
(c) Plot drag and power against speed (`drag_and_power`) and mark the two speeds. Compare V_mp with the clean stall speed of `P10.08` (26.4 m/s). What does that mean for an endurance flight?

Answer form: five numbers; three numbers; a plot and an argument.

## 6. (Judgment: when is a span efficiency wrong?)

(a) With the control points at the y-midpoints of a cosine layout, the lesson's elliptic wing gives e = 1.031. Explain why that must be a numerical error, and what is wrong with the layout.
(b) The lesson's model keeps the wake flat. Name one quantity that a rolled-up wake would change a lot and one it would hardly change.
(c) A student uses e = 1 for the trainer's polar "because induced drag is computed exactly". Write two sentences of review.

Answer form: three short arguments.
