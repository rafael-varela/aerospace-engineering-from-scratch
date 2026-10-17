# P08.16  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer, label every drag coefficient as exact, correlation (with its name and range) or assumed, and say when you use a correlation outside its range.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: a pitot mast)

The light aeroplane's pitot mast is a round tube of diameter d = 10 + s/2 mm standing 0.3 m into the stream, at 60 m/s TAS and 3000 m (ISA).

(a) Compute Re_d, c_d from White's correlation, the drag, and the power it costs.
(b) Compute the shedding frequency from Roshko's correlation.
(c) The mast's first bending frequency is 90 Hz. At what airspeed does the shedding frequency equal it? Is that speed inside the aeroplane's operating range, and what would you do about it?

Answer form: four numbers; a frequency; a speed and two sentences.

## 2. (Derivation: pressure drag, the wake model and Stokes)

(a) Derive c_d,p = ½∫₀^{2π} c_p cos θ dθ for a circular cylinder, with θ from the front stagnation point. Show that any c_p symmetric under θ → π − θ gives zero drag, and that any c_p symmetric under θ → −θ gives zero side force.
(b) Derive the closed form c_d = sin θ_s − (4/3) sin³θ_s − c_pb sin θ_s of the separated-wake model.
(c) Derive the Stokes terminal velocity from c_d = 24/Re, and show that in general the terminal Reynolds number depends only on the Archimedes number Ar = (ρ_p − ρ_f)ρ_f g d³/μ².

Answer form: three derivations.

## 3. (Numerical, uses your code: the ping-pong ball of Phase 04)

A ping-pong ball (d = 40 mm, m = 2.7 g) falls in air at 20 °C.

(a) Find its terminal velocity, Re and c_d with Morrison's and with Schiller–Naumann's correlation. How far apart are they?
(b) With `simulate_fall`, find the time and the distance to reach 95 % of the terminal velocity.
(c) In the drop experiment of `projects/p04-drag-drop` the ball falls about 2 m. What speed does it reach there, as a fraction of the terminal velocity, and how much later does it land than it would in a vacuum? Could that experiment measure c_d well? Explain.

Answer form: two sets of three numbers and a percentage; a time and a distance; a speed, a fraction, a time difference and two sentences.

## 4. (Unfamiliar scenario: a sounding-rocket payload coming down)

A 2 kg payload is recovered under a hemispherical canopy with c_d ≈ 1.4 (a typical value for such canopies, based on the canopy's projected area).

(a) Find the canopy diameter that lands it at 6 m/s at sea level. With that canopy, what is the descent speed at 3000 m?
(b) If the canopy fails, the payload (a smooth 0.12 m sphere) falls freely. Find its terminal velocity at sea level with Morrison's and with Schiller–Naumann's correlation (μ = 1.789 × 10⁻⁵ Pa·s, ρ = 1.225 kg/m³). Why do they differ by almost a factor of two, and is either one inside its range?
(c) Which number would you give to the range-safety officer, and why?

Answer form: two lengths and a speed; two speeds, an explanation and a range check; a number with a reason.

## 5. (Numerical, uses your code: how thick should a fairing be?)

(a) For the lesson's 25 mm strut at 60 m/s and 3000 m, compute the fairing's c_d on thickness and the drag saving for t/c = 0.1, 0.15, 0.25, 0.33, 0.5.
(b) Explain the optimum from the two parts of the form-factor correlation.
(c) A colleague proposes a 6:1 fairing (t/c = 0.167) "because it looks sleeker". What does it cost relative to the optimum, in newtons per metre?

Answer form: a table; two sentences; a number.

## 6. (Judgment: what the wake model can and cannot tell you)

(a) The separated-wake model with θ_s = 82° and c_pb = −1.2 gives c_d = 0.88 for a cylinder at Re = 10⁵, where White's correlation gives 1.0. List three reasons for the difference.
(b) Which of the two inputs, θ_s or c_pb, would you measure first on a wind-tunnel model, and with what instrument?
(c) A golf ball and a smooth sphere of the same size travel at 60 m/s. Which has less drag, and why? Would the same trick help a strut on the light aeroplane?

Answer form: three reasons; a choice and an instrument; a choice, a reason and a judgment.
