# P01.15  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. Angles in degrees in your answers, radians inside your code. Say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: one attitude, by hand and by code)

An aircraft has yaw ψ = 3s°, pitch θ = s/3°, roll φ = −s/2°.

(a) Write down the three elementary matrices with numbers, multiply them in the correct order, and give C_b/n to four decimals.
(b) Verify by hand that the first row has unit length and is perpendicular to the third row, and compute the determinant by expansion along the first row.
(c) Recover (ψ, θ, φ) from your matrix with the formulas of lesson 2.5, and confirm with `euler_321_from_dcm`.
(d) Give the gravity vector in body axes and the single axis-angle form of the attitude.

Answer form: a 3×3 matrix, two checks with numbers, three angles, one vector and one axis-angle pair.

## 2. (Derivation: the group properties and gimbal lock)

(a) Prove that if C₁ and C₂ are orthogonal with determinant +1, so is C₁C₂, and so is C₁ᵀ.
(b) Prove that for any orthogonal C and any vectors u, v: (Cu)·(Cv) = u·v. Conclude that |Cv| = |v|.
(c) Starting from the full 3-2-1 matrix in lesson 2.4, show that at θ = +90° the second row reduces to (sin(φ − ψ), cos(φ − ψ), 0) and the third row to (cos(φ − ψ), −sin(φ − ψ), 0). Explain what this says about which attitudes are distinguishable there, and derive the corresponding result for θ = −90°.
(d) Show that the extraction formulas ψ = atan2(C₁₂, C₁₁), φ = atan2(C₂₃, C₃₃) return the input angles for any ψ, φ ∈ (−180°, 180°] and |θ| < 90°, including the sign quadrant.

Answer form: four short proofs.

## 3. (Unfamiliar scenario: a star tracker on a spacecraft)

A star tracker returns the DCM from the inertial frame to its own optical frame, C_t/i. It is bolted to the spacecraft so that the tracker frame is the body frame turned by 120° about the body diagonal (1, 1, 1)/√3 (this permutes the axes cyclically). The tracker reports an attitude equivalent to yaw 40°, pitch 20°, roll 10° *of the tracker frame* with respect to inertial.

(a) Write C_t/b with `dcm_from_axis_angle` and explain, from its entries, which body axis each tracker axis points along.
(b) Find C_b/i by composition. Write the product in the right order and state why the other order is wrong.
(c) A thruster fires along body +x. Give the thrust direction in the inertial frame.
(d) Give the 3-2-1 angles of the *body* attitude. Compare with the tracker's own angles and comment: are Euler angles of a composed rotation simply the sums of the angles?

Answer form: three matrices, one unit vector, three angles, two sentences.

## 4. (Numerical, uses your code: the integration experiment with your own numbers)

Repeat the computational experiment of lesson part 4 with rate ω = s °/s about axis (1, 2, 2)/3 for 20 s, N = 200, 400, 800, 1600.

(a) Tabulate orthogonality error and misalignment with and without renormalization. Fit the convergence order of each.
(b) Compare the misalignment at each N with ωT(ωΔt)²/3 and report the largest relative deviation.
(c) Replace the small-angle step by `dcm_from_axis_angle(axis, rate * dt)` and repeat. Report the orthogonality error and misalignment you now get, and explain why the misalignment is exactly zero (or as near as the arithmetic allows).
(d) With first-order steps, what step size would keep the misalignment below 0.1° over the 20 s, for your ω? Give the sample rate in Hz.

Answer form: two tables, fitted orders, one plot, one sample rate.

## 5. (Derivation: the Rodrigues formula and its inverse)

(a) Decompose a vector **v** into a part along the unit axis **n** and a part perpendicular to it. Turn the *frame* by angle *a* about **n** and write the components of **v** in the new frame in terms of **v**, **n**(**n**·**v**) and **n** × **v**. Collect terms into the matrix form of lesson 2.6 and check the sign of the skew term against R₃.
(b) Show that tr C = 1 + 2 cos *a* and C − Cᵀ = −2 sin *a* [**n**×].
(c) Explain why `axis_angle_from_dcm` uses atan2(|**w**|, tr C − 1) rather than arccos((tr C − 1)/2). Estimate the angle error of the arccos version when the trace carries a rounding error of 2 × 10⁻¹⁶ near *a* = 0.
(d) Show that for *a* = 180° the matrix is 2**nn**ᵀ − I, that it is symmetric, and that **n** and −**n** give the same matrix.

Answer form: a derivation, two identities, one estimate, one short argument.

## 6. (Judgment: when the model is the wrong model)

(a) A spinning sounding rocket rolls at 10 rev/s. Its flight computer samples the gyros at 100 Hz and updates the DCM with the first-order step I − [δθ×]. Compute the rotation per sample in radians, the orthogonality error of a single step, and the angle a single step actually turns the frame through. Is this a numerical-accuracy problem or a wrong-model problem?
(b) Using the misalignment formula, what sample rate would first-order steps need to hold 0.1° over a 10 s flight? Is that sensible? Propose the alternative.
(c) A vertical launch passes through θ = 90°. Which functions in your library stay valid there and which do not? What would you propagate instead of the three angles?
(d) Name one assumption of the lesson that fails for a 2-hour navigation run and estimate its size in degrees.

Answer form: four short paragraphs, each with at least one number.
