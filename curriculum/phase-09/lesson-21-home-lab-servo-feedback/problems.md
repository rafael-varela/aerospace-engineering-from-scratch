# P09.21  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact, which are numerical and which are SYNTHETIC. Problems that describe hardware you may build stay inside the limits of `projects/p09-servo-control/safety.md`; Problem 3 is a paper exercise only.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: a P gain for a gear motor)

A gear motor identified as K/(s(τs + 1)) has K = 400 + 10s deg/s per unit duty and τ = 30 + s ms, with a dead zone of 0.10 compensated by 0.09 in the firmware. The target is (5 + s/2) % overshoot.

(a) Find ζ, the gain Kp (duty per degree) and ω_n.
(b) Predict the 10–90 % rise time, the peak time and the 2 % settling time on the continuous model, and compare the settling time with the estimate 4/(ζω_n).
(c) What duty does the firmware ask for at the first sample of a 20° step? With a duty limit of 0.6, what is the largest step the linear prediction holds for?
(d) What steady-state error band does the under-compensated part of the dead zone leave?

Answer form: three numbers; four times; a duty and an angle; an angle.

## 2. (Derivation: the four results the lab is built on)

(a) Derive the step response (2) of K e^{−Ls}/(s(τs + 1)), show that its tail is the line Ku(t − τ − L), and that the difference between the response and that line is Kuτ e^{−(t − L)/τ} for t ≥ L.
(b) Derive ζ and ω_n of the unity-feedback loop with gain Kp around K/(s(τs + 1)), and the gain that gives an overshoot M_p.
(c) Show that PI control with Ki/Kp = 1/τ around K/(τs + 1) gives the first-order closed loop 1/(τ_cl s + 1) with τ_cl = τ/(KKp), and that its 10–90 % rise time is τ_cl ln 9 and its 2 % settling time τ_cl ln 50.
(d) Approximate a loop delay θ by 1/(θs + 1) in the loop of (c). Show that the closed loop is second order with ζ = 1/(2√(θω_c)), ω_c = KKp/τ. For θω_c = 1/4, show that the 10–90 % rise time is shorter than the first-order loop's, and by how much.

Answer form: four derivations; in (d) a percentage.

## 3. (Unfamiliar scenario, paper only: a camera-gimbal pointing loop)

A small UAV's camera gimbal has a pitch axis whose motor controller accepts a rate command; from the outside, rate command to angle is 1/(s(τs + 1)) with τ = (15 + s/2) ms (deg per deg/s). An outer pointing loop on the flight computer applies a P gain Kp (deg/s of rate command per degree of error), written one sample after it is computed. **Paper exercise: do not build a gimbal or fly anything.**

(a) Design Kp for 5 % overshoot on the continuous model; predict the rise and settling times.
(b) The outer loop runs at 100 Hz. Estimate the total loop delay, the crossover frequency and the phase margin with and without it, and predict the overshoot with the delay modelled as a lag. Check the prediction with `simulate_closed_loop`.
(c) Repeat (b) at 400 Hz.
(d) If 5 % overshoot is a hard requirement, what gain must the 100 Hz loop use and what does it cost in rise time? Which loop rate do you recommend?

Answer form: a gain and two times; a delay, a frequency, two margins, two overshoots; the same at 400 Hz; a gain, a time and a recommendation.

## 4. (Numerical, uses your code: the full analysis on your seed)

Write the SYNTHETIC dataset with `write_synthetic_csvs(folder, seed=s)` (or use your own measured data, and say which) and run `analyse_dataset(folder)`.

(a) Report the identified motor K, τ and dead zone, the three per-run fits with standard errors, and the servo K, τ, L.
(b) For each closed-loop run, give the continuous prediction, the sampled-model prediction and the synthetic values of rise time, overshoot and settling time.
(c) Repeat for seeds s to s + 9 and give the mean and standard deviation of the synthetic 15° overshoot and of the sampled model's prediction of it. Is their difference significant? What do you conclude?
(d) Label every number SYNTHETIC if it came from the generator.

Answer form: two parameter tables; a comparison table; four numbers and a sentence; the label.

## 5. (Design on the refined model: a faster servo loop)

Using the Option A model identified from the shipped SYNTHETIC dataset (or your own), keep the cancellation Ki = Kp/τ and find the Kp for which the **sampled** loop (Ts = 20 ms, one-sample delay, plant dead time L, ±60° command limit, quantisation) has a 10–90 % rise time of 0.20 s.

(a) Report Kp and Ki (rounded for the firmware) and the sampled-model rise time, overshoot and settling time.
(b) What rise time does the continuous model predict for the same gains?
(c) What would the continuous design for a 0.20 s rise have given on the sampled model? Which design would you put on the board, and what would you write in the dated prediction?

Answer form: two gains and three numbers; a time; three numbers and a decision.

## 6. (Judgment: reading a closed-loop test)

For each observation, name the effect you suspect and the one test that would decide it: (a) a PI servo loop holds the pointer within 1° but hunts slowly back and forth around the target forever; (b) the measured overshoot is twice the dated prediction and the rise time slightly shorter; (c) a 60° step rises twice as slowly as a 15° step with the same gain, with less overshoot; (d) a P loop on the gear motor stops 0.7° short, but only on steps in the negative direction; (e) the servo's identified τ and L grow with the step size; (f) at small steps the motor hums but does not move; (g) every metric drifts slowly over a 40-minute session.

Answer form: seven labelled lines.
