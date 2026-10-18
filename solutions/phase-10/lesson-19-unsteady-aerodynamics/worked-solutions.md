# P10.19  Worked solutions

Each solution states assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15. Numbers come from `solution.py` and a probe script that calls it, given to the figures shown. Model throughout: a two-dimensional thin airfoil in incompressible, attached flow with a flat wake (Theodorsen's assumptions).

## 1. How unsteady is it?

(a) Sailplane: b = 0.45 m, ω = 2π·3 = 18.85 rad/s, **k = 18.85 × 0.45/25 = 0.339**, highly unsteady; |C| = **0.670**, phase **−15.0°**.

(b) Trainer: b = 0.6 m, ω = 31.42 rad/s, **k = 0.343**, highly unsteady; |C| = **0.669**, phase **−15.0°**. (The two wings sit at almost the same k by coincidence: the ratio fb/V is nearly equal.)

(c) Rotor at 0.75 R: Ω = 210/5 = 42 rad/s, V = 157.5 m/s, b = 0.125 m, **k = 42 × 0.125/157.5 = 0.0333**, quasi-steady by the k < 0.05 rule of thumb; |C| = **0.945**, phase **−6.3°**. Total lift for pitch about c/4: **0.944 × 2πα̂, phase −3.4°**.

**Check.** |C| lies between 1/2 and 1 and the phase is negative in all three, as it must (the wake always reduces and delays the circulatory lift). Note that k = 0.033 still costs 5.5 % of the circulatory amplitude: C(k) approaches 1 only slowly (|1 − C| is still 0.049 at k = 0.01 and 0.007 at k = 0.001).

**Common mistakes.** Using the chord instead of the semi-chord (doubles k); using the tip speed with the 0.75 R chord position.

## 2. From the indicial function to C(k)

(a) With the input started at s = −∞, cl_c(s) = 2π ∫_{−∞}^{s} φ(s − σ) α₃/₄′(σ) dσ. For α₃/₄ = e^{ikσ}, α′ = ik e^{ikσ}; put τ = s − σ:
cl_c = 2π e^{iks} · ik ∫₀^∞ φ(τ) e^{−ikτ} dτ = 2π α₃/₄(s) · ik Φ(ik).
So **C(k) = ik Φ(ik)**: Theodorsen's function is the frequency response of the indicial (Wagner) system.

(b) Φ(p) = 1/p − Σ Aᵢ/(p + βᵢ), so **C(k) = 1 − Σ Aᵢ ik/(ik + βᵢ)**. As k → 0 every fraction → 0 and C → 1: slow motion has a fully developed wake and quasi-steady lift. As k → ∞ every fraction → 1 and C → 1 − ΣAᵢ = φ(0) = ½: so fast that the shed vorticity stays just behind the trailing edge and halves the circulatory lift, the same factor Wagner's function gives at the instant of a step.

(c) ik∫φe^{−ikτ}dτ = [−φe^{−ikτ}]₀^∞ + ∫φ′e^{−ikτ}dτ = φ(0) + ∫₀^∞ φ′(τ) e^{−ikτ} dτ (the boundary term at infinity vanishes in the steady-state sense). The integrand φ′ decays exponentially, so a finite numerical integral converges; φ itself tends to 1 and its transform has to be interpreted as a limit. The lesson's test computes this integral with the trapezoid rule to 600 semi-chords and matches the closed form to 2 × 10⁻⁵.

(d) 1 − ½e^{−βs} = 0.9 ⇒ s₉₀ = ln 5/β. Matching 12.7: **β = 0.1267**. A one-term fit has a single time scale, so it cannot have both Jones's quick early rise (75 % at 3.7 semi-chords) and his slow tail; its C(k) is less accurate (the two-term fit has two poles to shape the phase).

**Common mistakes.** Forgetting the factor ik (the derivative in Duhamel's integral); taking C(∞) = 0.

## 3. A rotor blade in cyclic pitch

**Assumptions.** Two-dimensional, attached flow; steady free stream 157.5 m/s; pitch about c/4; α̂ = 3° = 0.0524 rad; k = 1/30.

(a) Quasi-steady 2π α₃/₄ = 2π α̂ (1 + ik): amplitude **0.3292**, phase **+1.9°** (the α̇ term at the three-quarter chord leads). Theodorsen: amplitude **0.3105**, phase **−3.4°**. The steady value 2πα̂ is 0.3290. So the unsteady lift is 5.6 % smaller than quasi-steady and lags the pitch by **3.4° of azimuth** instead of leading it by 1.9°: a difference of 5.3° of azimuth.

(b) A rotor's flapping response already lags the aerodynamic forcing by about 90° (it is a resonant once-per-revolution system), and the pilot's cyclic input is phased to put the rotor's tilt in the right direction. A few degrees of extra aerodynamic lag rotate the response, so a forward stick input produces a little lateral tilt: cross-coupling. Control designers correct the phase; they cannot if they model the lift as quasi-steady.

(c) (1) The shed wake of a rotor does not trail away flat: it passes under the following blades (returning wake), which Loewy's extension of Theodorsen treats; it usually reduces the lift further at some frequencies. (2) Near the retreating-blade side, local angles approach stall and dynamic stall adds large overshoots and hysteresis that thin-airfoil theory cannot produce. Compressibility on the advancing side (M ≈ 0.6 at the tip) is a third.

**Common mistakes.** Treating 3.4° as negligible because the amplitude change is small; forgetting that the quasi-steady α₃/₄ already contains a phase lead.

## 4. A pitching section of your own (s = 15: k = 0.25, a = −0.35)

(a) Theodorsen: **|cl/α̂| = 4.563, phase +6.87°**; quasi-steady: **6.423, phase +12.00°**. The wake costs 29 % of the amplitude and 5° of phase.

(b) `simulate_harmonic_pitch(0.25, radians(2), -0.35)`: after 10 cycles the fitted amplitude is **4.635** at **+6.27°**, 1.9 % from exact Theodorsen. That difference is Jones's fit: against C(k) built from the same Jones approximation the time-domain result agrees to about 10⁻⁵ (discretisation). Jones's |C_J − C| is 0.0145 at most over 0 < k ≤ 5 (2.3 % of |C|).

(c) The loop area is **+0.00209** (radians × cl) from the sampled loop and the same from π Im(ᾱ* ĉl) to five figures. Positive area means the loop is traversed so that lift is higher on the upstroke than on the downstroke at the same α: the lift leads (here, because the α̇ terms with a pivot behind c/4 outweigh the wake lag). A negative area would mean the lift lags.

(d) Aerodynamic work on the pitch motion per cycle (α̂ = 1 rad): **−1.105**, negative: the air damps the pitching. Lift-loop area is not the pitch work: work is ∮ cm dα, not ∮ cl dα.

**Independent check.** Doubling the number of time steps per cycle changes the fitted amplitude in the sixth figure; the Duhamel quadrature converges at order 2.00 against a closed-form convolution (lesson §3).

**Common mistakes.** Fitting the whole simulation, transient included; confusing ∮ cl dα with the work done.

## 5. Gust alleviation (s = 15: H = 8 chords)

(a) Quasi-steady peak Δcl = 2π × 5/55 = **0.571**. With Küssner's function the peak is **0.480**, ratio **0.840**.

(b) H = 12.5 chords: ratio **0.903**. A longer gust ramps up over more chords, so the lift has time to build (Küssner's function reaches 90 % in 12.4 semi-chords, about 6 chords); the peak arrives when the response is nearly complete.

(c) W/S = 750 × 9.81/12 = 613.1 N/m². Δn = ρ V CL_α w₀/(2 W/S) = 1.225 × 55 × 4.685 × 5/(2 × 613.1) = **1.29** quasi-steady, **1.08** with the 0.840 reduction. Real aircraft also gain relief from the aircraft rising (plunging) in the gust, which reduces the effective gust angle, and from the wing's own unsteady three-dimensional lift build-up (Jones TR-681); gust-load formulas combine these into an alleviation factor that depends on the aircraft's mass ratio.

**Common mistakes.** Using 2π for the wing (a section value); measuring H in semi-chords where the problem gave chords.

## 6. Single-degree-of-freedom pitch flutter

(a) From `pitch_flutter_boundary`:

| a | pivot x/c | energy fed in below k = |
|---|---|---|
| −0.6 | 0.200 | 0.0015 |
| −0.7 | 0.150 | 0.0141 |
| −0.8 | 0.100 | 0.0268 |
| −1.0 | 0.000 (leading edge) | 0.0403 |

(b) About the pivot, cm = (apparent-mass terms) + ((a + ½)/2) × cl_circulatory. The apparent-mass terms are pure stiffness and damping of the right sign (the −(½ − a)ik term is a damping proportional to α̇). The circulatory lift acts at the quarter chord; if the pivot is at the quarter chord (a = −½) its moment is zero, and if the pivot is behind (a > −½) the lift ahead of the pivot is a nose-up moment that grows with α₃/₄ = α + (½ − a)α̇/(V/b), whose α̇ part opposes the motion: damping. Only for a pivot well ahead of the quarter chord does the lagging circulatory lift, acting behind the pivot, push in phase with the velocity at low k. A scan over 10⁻⁴ ≤ k ≤ 10 finds the work negative everywhere for −0.5 ≤ a ≤ 1.

(c) Plunge alone is always damped (lesson test), and pitch alone is damped for any practical elastic axis (between 25 % and 50 % chord). Classical flutter needs the two motions coupled with a phase between them that lets the lift from one do positive work on the other; that phase is set by the structure's mass and stiffness distribution (the distance between the elastic axis and the centre of mass, the two natural frequencies). So a flutter calculation needs the structural dynamics plus these aerodynamic coefficients, solved together (`P13.18`).

(d) Attached flow. Near stall the lift and the moment overshoot and collapse with hysteresis (dynamic stall); the moment loop can reverse direction, so the pitch damping becomes negative: stall flutter, a single-degree-of-freedom instability that Theodorsen's attached-flow theory cannot predict.

**Common mistakes.** Concluding that a forward pivot is dangerous at any k (it only feeds energy below a small k); reading lift-loop area as the pitch work.

## Check-yourself answers

1. k = ωb/V with b the semi-chord: the number of radians of oscillation per semi-chord travelled. Quasi-steady below about 0.05, highly unsteady above about 0.2. The trainer at 5 Hz: k = 0.343.
2. As k → 0, C → 1 (quasi-steady lift 2πα₃/₄); as k → ∞, C → ½ (the wake halves the circulatory lift).
3. The non-circulatory (apparent-mass) lift comes from accelerating the air around the section, proportional to ḧ and α̇, α̈; it is in phase with the acceleration, not with α, and has no wake. The circulatory lift comes from the bound circulation and is the part the wake reduces and delays through C(k).
4. One half: φ(0) = 1 − 0.165 − 0.335 = 0.5. With Jones's approximation, 90 % after 12.7 semi-chords.
5. Because the system is linear: the lift is the convolution of the indicial response with the rate of change of the input. It fails when the flow separates (dynamic stall) or the motion is large enough to change the wake geometry.
6. A lift loop (cl against α) is an ellipse whose area is π Im(ᾱ* ĉl); its direction tells whether lift leads or lags. Whether the air feeds energy into a motion is decided by the work of the generalized force on that motion, ∮ cm dα for pitch, ∮ L dh for plunge, not by the lift loop.
7. Pure plunge: always damped. Pure pitch: damped for pivots at or behind the quarter chord; for pivots ahead of about 0.2c, energy is fed in at very low k.
