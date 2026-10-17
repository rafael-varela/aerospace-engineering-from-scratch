# P06.07  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. Angles in radians; `C = dcm_from_euler([Ω, i, ω], "313")`; `orthogonality_error(C)`, `np.linalg.det(C)`. 2. The orbit normal is perifocal axis 3; its inertial components are row 3 of C (because the rows of C are the body axes in reference components). Compare with the formula. 3. `euler_from_dcm(C, "321")` and `euler_from_dcm(C, "313")`.

## 2.
(a) Evaluate (2) with R₁(φ)R₂(θ)e₃, R₁(φ)e₂, e₁; solve the 3 × 3 system for (ψ̇, θ̇, φ̇) by elimination (θ̇ from q and r, then ψ̇, then φ̇). (b) det B = c₁ · (c₂ × c₃); c₂ ⊥ c₃ and c₂ × c₃ is a unit vector; for an asymmetric sequence the first axis is ⊥ to the second, and the angle between c₁ and c₂ × c₃ is a₂ (so the dot product is cos a₂); for a symmetric sequence c₁ starts parallel to c₃ and is tilted by a₂ (so sin a₂). (c) ĊCᵀ + CĊᵀ = 0 ⇒ ĊCᵀ skew = −[ω×]; d(CᵀC)/dt = ĊᵀC + CᵀĊ. (d) (2) with R₃(φ)R₁(β)e₃ ψ̇ + e₃ φ̇. (e) Write k₁…k₄ for constant A = −[ω×]: the step is (E + AΔt + (AΔt)²/2 + (AΔt)³/6 + (AΔt)⁴/24)C; the missing fifth-order term sets the orthogonality defect.

## 3.
(a) Rows: body x up = (0, 0, −1); body y east = (0, 1, 0); body z north = (1, 0, 0). (b) Loop over the twelve sequences, `euler_from_dcm`, `abs(det(euler_rate_matrix(…)))`. (c) Repeat at the four pitch angles for the eight non-singular sequences. (d) Quaternions have no singular attitude.

## 4.
(a) |ω| from `coning_body_rate(0, …)`; `integrate_dcm(coning_dcm(0, …), w, 60, n_samples=601)`; max misalignment and orthogonality error over the samples. (b) Loop over Δt; final misalignment against `coning_dcm(t[-1], …)`; log₂ of ratios. (c) `renormalize_every=1`. (d) Compare the Δt = 0.02 s error with 10⁻⁴.

## 5.
(a) r/sin δ. (b) For f in (50, 100, 200, 400): `integrate_segments_rk4(np.zeros(3), segs, 1/f, "euler")`, rebuild the DCM, `misalignment_angle` against the closed form; same with `"dcm"` from `np.eye(3)`. (c) log₂ of successive ratios; the asymptotic order needs ψ̇Δt ≪ 1. (d) Read off the table.

## 6.
(a) DCM or quaternion; angles for display only. (b) The parity sign ε. (c) Orthogonality ≠ accuracy. (d) The direction of the mapping (TM-74839 maps rotated-frame to original-frame components); transpose. (e) Mechanical gimbal lock loses a degree of freedom of the hardware; in software only the coordinates are singular. (f) Margin, rate excursions, failure cases; the cost of a singular representation in a recovery manoeuvre.
