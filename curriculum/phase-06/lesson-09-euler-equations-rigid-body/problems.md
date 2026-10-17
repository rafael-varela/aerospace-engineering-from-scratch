# P06.09  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: a satellite with a boom)

The `P04.04` box satellite has I = diag(260.4, 208.3, 135.4) kg·m² about its centre of mass in its body axes. A (10 + s/5) kg instrument is mounted on a light boom at r = (1.5, 0.5, 0) m from the centre. Ignore the shift of the centre of mass.

(a) The instrument's contribution m(|r|²E − rrᵀ) and the new tensor. Which products of inertia appear, and why only those?
(b) The principal moments and axes with `principal_axes`.
(c) The angle between each principal axis and body axis 1. Which principal axis is the new intermediate axis?

Answer form: a matrix and a sentence; three moments and a matrix; three angles and a sentence.

## 2. (Derivation: the geometry and the exact solution)

(a) From (2) with τ = 0, show that T and H² are constant, and that H_i is constant in the inertial frame.
(b) Derive (5) and show that for H² = 2TI₂ the polhodes lie in the two planes I₁(I₂ − I₁)ω₁² = I₃(I₃ − I₂)ω₃².
(c) Show that the tip of the inertial rate vector stays in a fixed plane perpendicular to H at distance 2T/|H| from the origin.
(d) Substitute ω₁ = A₁ cn τ, ω₂ = A₂ sn τ, ω₃ = A₃ dn τ with τ = λt into (3) and, using (6), derive the three relations between A₁, A₂, A₃, λ and m. Then use (4) to obtain (7).
(e) Show that as m → 0 the period 4K/λ tends to the `P04.04` wobble period for spin about the major axis, and that as m → 1 the flip interval tends to (2/σ) ln(4/k′).

Answer form: five derivations.

## 3. (Unfamiliar scenario: approaching a tumbling satellite)

A servicing spacecraft must approach a defunct satellite with principal moments (1200, 1900, 2300) kg·m². Its tumble, measured from radar and images, is ω = (0.02, 0.03 + s/1000, 0.05) rad/s in principal axes.

(a) Its kinetic energy and |H|. Is its polhode around the major or the minor axis?
(b) The period of the body-rate motion and the parameter m. How long does one full cycle of its tumble take?
(c) Simulate one period. Which body axis stays nearest the angular-momentum direction, and between what angles from H does it swing?
(d) The servicer wants to dock along the target's angular momentum direction. Why is that direction the natural choice, and what does (c) tell it about the target's axis?

Answer form: two numbers and a word; two numbers; an axis and two angles; two sentences.

## 4. (Numerical, uses your code: your own T-handle)

Spin the lesson's T-handle at (2 + s/10) rev/s about axis 2 with a disturbance of 10^−(1 + (s mod 4)) × spin on axis 1.

(a) The polhode type, m, λ, the period and the flip interval from `elliptic_parameters`; the growth rate σ.
(b) Simulate three periods at rtol = 10⁻¹¹ (atol 10⁻¹²). The flip times, their mean interval, the largest deviation from `elliptic_solution` and the three invariant drifts.
(c) Repeat at rtol = 10⁻⁶ (atol 10⁻⁷). What degraded and what did not? Explain with (9).
(d) How many flips does the handle make in one minute?

Answer form: six numbers; a list, two numbers and three drifts; a sentence with numbers; a number.

## 5. (Numerical, uses your code: the box satellite of P04.04)

The box satellite, I = diag(260.4167, 208.3333, 135.4167) kg·m², spins at 0.5 rad/s about axis 2 with 10⁻⁶ rad/s on axis 1.

(a) σ, the flip interval and the time of the first reversal of ω₂ (half a flip interval, since the spin starts at the top of the sn curve).
(b) `P04.04` estimated the reversal from linear growth as ln(0.5/10⁻⁶)/σ. Compare, and explain the difference.
(c) Confirm (a) with `simulate` and `flip_times` over 200 s.
(d) An engineer proposes to "make the spin stable" by reducing the disturbance to 10⁻⁹ rad/s. How much later would the first reversal come?

Answer form: three numbers; two numbers and a sentence; a number; a number and a sentence.

## 6. (Judgment: free rotation in practice)

For each, give a short argued answer (one to three sentences).

(a) A colleague's tumble simulation conserves T to 10⁻¹² but the inertial H drifts by 1 % in an hour. Name two likely causes.
(b) Why can one not predict the exact time of a Dzhanibekov flip from a video of a real T-handle, while the interval between flips is well predicted?
(c) A satellite designer lists the moments (300, 301, 500) kg·m². Is spin about the 301 axis dangerous? About the 300 axis?
(d) Does a tumbling satellite's angular momentum vector move? What about its angular velocity vector?
(e) The elliptic solution gives ω(t) exactly. Why does the lesson still integrate the quaternion numerically?
(f) Which assumption of this lesson would you test first on a real spent rocket stage, and how?

Answer form: six labelled items.
