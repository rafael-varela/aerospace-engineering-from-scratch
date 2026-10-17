# P08.12  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact eigenvalues, which come from the supplied **approximate-fit** stability chart (Drela and Giles's formulas in `starter/main.py`), and which are correlations (Mack's N_crit).

Several problems use your personal seed **s**: the day of the month you were born (1–31).

Air: sea level ν = 1.4607 × 10⁻⁵ m²/s; 3000 m ν = 1.8630 × 10⁻⁵ m²/s. Water at 20 °C: ρ = 998.2 kg/m³, μ = 1.002 × 10⁻³ Pa·s.

## 1. (Warm-up: a sailplane wing)

A sailplane flies at 30 m/s near sea level; its wing chord is c = 0.6 + s/50 m. Treat the upper surface as a flat plate (H = 2.591, Re_θ = 0.664√Re_x).

(a) Compute Re_c, and with `flat_plate_transition` the instability point x₀ and the e^9 transition point x_tr. Is the layer laminar to the trailing edge?
(b) The same wing section is tested in a wind tunnel whose free-stream turbulence intensity is 0.3 %. Use Mack's correlation for N_crit and recompute x_tr/c.
(c) Which of (a) and (b) would you trust for the aircraft, and why does it matter for a sailplane in particular?

Answer form: a Reynolds number, two lengths and a yes/no; an N and a fraction; two sentences.

## 2. (Derivation: the Orr–Sommerfeld equation and what N means)

(a) Linearise the two-dimensional incompressible Navier–Stokes equations about a parallel flow U(y), with disturbance stream function ψ′ = φ(y)exp(iα(x − ct)) (so u′ = ψ′_y, v′ = −ψ′_x). Eliminate the pressure by cross-differentiation and show that φ satisfies (U − c)(φ″ − α²φ) − U″φ = (φ⁗ − 2α²φ″ + α⁴φ)/(iαRe), with φ = φ′ = 0 at a wall. State which terms you neglected and why the parallel-flow assumption is only approximate for a boundary layer.
(b) Set Re → ∞ (Rayleigh's equation). Multiply by the complex conjugate φ*/(U − c), integrate across the channel and take the imaginary part to show that an unstable inviscid mode (c_i > 0) needs U″ to change sign somewhere in the flow. Is plane Poiseuille flow inviscidly stable? Is the Blasius profile? (Use f‴ = −½ff″ from `P08.10`.) What, then, makes them unstable?
(c) A disturbance grows as e^N. How large is the amplification for N = 9 and N = 7? A typical free-stream disturbance in quiet flight has a velocity amplitude of order 10⁻⁵U at the instability point; estimate its size at N = 9 and say why that is a reasonable "transition" amplitude.

Answer form: a derivation; a proof sketch and two answers; two numbers and a sentence.

## 3. (Unfamiliar scenario: a rocket engine's cooling channel)

A small regeneratively cooled rocket engine has flat cooling channels; model one as plane Poiseuille flow of water at 20 °C between plates 2h apart, with h = (1 + s/10) mm.

(a) Using the critical Reynolds number 5772.22 (based on the centreline velocity U_c and the half-height h), find the centreline and the bulk velocity (U_b = 2U_c/3) at which linear instability first appears.
(b) Express that critical point as a bulk Reynolds number based on the full height, U_b·2h/ν, and on the hydraulic diameter of a wide channel (4h).
(c) Experiments in channels see sustained turbulence at bulk Reynolds numbers far below the linear critical value. What does that say about how transition happens in channel flow, and why would you not use the e^N method to design this cooling channel?

Answer form: two velocities; two Reynolds numbers; a short paragraph.

## 4. (Numerical, uses your code: the neutral curve)

(a) Let α₁ = 0.95 + s/500. Compute the neutral Reynolds number at α₁ with `neutral_reynolds` for N = 30, 40, 50, 60. How many digits are converged at each N?
(b) Compute the lower branch of the neutral curve at α = 0.90, 0.95, 1.00, 1.05, 1.08 (N = 60) and sketch it. Where is its minimum, and how does that compare with `critical_point(60)` and the published (5772.22, 1.02056)?
(c) Find the eigenvalue spectrum at α = 1, Re = 10 000 (N = 60): how many modes grow? Plot c_r against c_i. Which part of the spectrum is unreliable at this N, and how would you tell?

Answer form: four numbers with converged digits; five numbers, a sketch and a comparison; a count, a plot and a sentence.

## 5. (Numerical, uses your code: what the transition point depends on)

The light aeroplane's wing at 3000 m and 60 m/s, chord 1.5 m.

(a) The aircraft flies through slightly turbulent air with Tu = (0.05 + s/100) %. Compute N_crit by Mack's correlation and the flat-plate x_tr.
(b) Keep N_crit = 9 and change only the shape factor fed to the chart: H = 2.55, 2.591, 2.61, 2.65. Report Re_θ0 and x_tr. Thwaites's method (`P08.11`) gives H = 2.61 for a flat plate; Blasius gives 2.591. How much does that difference move x_tr?
(c) Using `thwaites_wedge`, compute x_tr for wedge flows with m = −0.05, 0, +0.05. Which matters more for this wing, the pressure gradient or the free-stream turbulence?

Answer form: an N and a length; a table and a percentage; three lengths and a sentence.

## 6. (Judgment: when not to trust e^N)

(a) List three routes to transition that the e^N method in this lesson does not model, and for each give an aircraft situation where it would decide the transition point.
(b) The fit used here was made for two-dimensional Falkner–Skan profiles. Name one feature of a real swept wing that would invalidate it.
(c) A designer claims laminar flow to 60 % chord on the light aeroplane on the strength of an e^9 calculation. What one flight measurement would you ask for before believing it, and what result would contradict the claim?

Answer form: three items with examples; one item; a measurement and a criterion.
