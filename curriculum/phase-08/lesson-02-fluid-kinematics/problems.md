# P08.02  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical. Sketch each line family before you compute it.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: the centreline of a wind-tunnel contraction)

Near its axis, the flow in a 2-D wind-tunnel contraction is modelled by u = U₀(1 + x/L), v = −U₀ y/L, with U₀ = 5 m/s and L = 2 + s/10 m.

(a) Write the velocity-gradient tensor. Find the divergence, the vorticity, the strain-rate tensor and the principal strain rates. Is the flow incompressible? Irrotational?
(b) Find the acceleration of the fluid particle at (1 m, 0.2 m). Which part is local and which convective? Check with `particle_acceleration`.
(c) Solve dx/dt = u for the particle on the axis that starts at x = 0, and find how long it takes to reach x = L. Show that y(1 + x/L) is constant along every pathline, and say what that means for the streamlines.

Answer form: a tensor and four quantities; an acceleration vector and a sentence; a time, a proof and a sentence.

## 2. (Derivation: lines and derivatives)

(a) For the oscillating stream u = U, v = V cos ωt, derive the pathline of the particle at (x₀, y₀) at t₀, the streamline through (x₀, y₀) at time t, and the streakline at time t from a source at (x_s, y_s). Eliminate the release time to show the streakline is a sine wave of wavelength 2πU/ω, and find its slope at the source.
(b) Prove that in a steady flow the streamline, the pathline and the streakline through a point coincide. Which property of the ODE dx/dt = u(x) do you need?
(c) Derive the material derivative Df/Dt = ∂f/∂t + u·∇f from the chain rule applied to f(x(t), y(t), t) along a pathline.
(d) Show that the rate of change of the area δA of a small material element is (∇·u) δA, and that the antisymmetric part W of the velocity gradient rotates a small element at the angular velocity ω/2, where ω = ∂v/∂x − ∂u/∂y.

Answer form: three line equations, a wavelength and a slope; a proof; a derivation; two derivations.

## 3. (Unfamiliar scenario: a smoke plume in a veering wind)

Over a few minutes on a coast, the wind keeps its speed W = 5 m/s while its direction turns steadily anticlockwise: u = W cos Ωt, v = W sin Ωt, with period 2π/Ω = 600 + 10s s. A chimney at the origin releases smoke continuously from t = 0.

(a) Find the pathline of the smoke particle released at t = 0, and the streamline through the chimney at any time. Describe their shapes.
(b) Find the streakline at t = T/4 analytically and show that it is also an arc of a circle. Where is its centre, and what is its radius?
(c) Check (a) and (b) with `pathline`, `streakline` and `streamline` (t from 0 to T/4, 300 steps). What is the largest distance between the computed points and your circles?
(d) An observer photographs the smoke at t = T/4 and draws the wind direction from its shape near the chimney. Is the drawing right there? Far from the chimney?

Answer form: two equations and two shapes; an equation, a centre and a radius; two distances; two sentences.

## 4. (Numerical, uses your code: when does smoke show streamlines?)

Your gust tunnel runs at U = 15 + s/5 m/s and swings the flow direction by ±5° at f = 1 + s/10 Hz.

(a) Compute V, ω, the streakline wavelength and the pathline amplitude V/ω.
(b) With `smoke_error`, tabulate the worst-phase smoke error over model lengths L = 0.1, 0.2, 0.3, 0.5 and 1.0 m against the reduced frequency k = ωL/U.
(c) Find the largest model length for which the smoke error stays below 5 % (use a root finder on the worst-phase error). Compare with the small-k estimate error ≈ k/2 and explain where k/2 comes from.
(d) Would you report the smoke photographs of a 0.3 m model as streamline pictures? What would you do instead?

Answer form: four numbers; a five-row table; a length, a reduced frequency and a derivation; a judgment.

## 5. (Numerical, uses your code: stretching a material loop)

(a) Place a circular material loop of radius 0.1 m centred at (0.3 m, 1.0 m) in the stagnation flow u = ax, v = −ay with a = 1 s⁻¹, represented by 800 points. Advect it with `rk4_step` to t = 2 s (dt = 0.01 s). Report the ratio of final to initial area and of final to initial perimeter, and the stretch of the loop's x-extent. Compare the stretch with exp(at).
(b) Repeat in the field u = ax, v = 0. What happens to the area, and how is that predicted by `decompose_gradient`?
(c) The area ratio in (a) is 1 to about 10⁻¹² while the perimeter ratio is near 4.7. Why does the area stay so accurate while the loop becomes a long thin ellipse? What would eventually limit the accuracy?

Answer form: three ratios and a comparison; a ratio and a sentence; a short argument.

## 6. (Judgment: which line does a flow-visualisation method show?)

For each technique say whether it shows streamlines, pathlines, streaklines or something else, and when the difference matters: (a) a smoke wand upstream of a model in a steady tunnel; (b) the same with a gust generator running; (c) a long-exposure photograph of small helium bubbles; (d) particle image velocimetry (two short exposures a few microseconds apart, giving a velocity field); (e) wool tufts on a wing in flight; (f) dye injected from a port on the nose of a model in a water tunnel while the model pitches.

Answer form: six labelled lines.
