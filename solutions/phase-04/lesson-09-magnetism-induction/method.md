# P04.09  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. v = (1 + s/10) × 10⁶ m/s. F = evB. 2. r = m v/(eB); T = 2π m/(eB). 3. Check with `gyroradius(M_PROTON, v, Q_PROTON, 30e-6)`. 4. Electron: same formulas with m_e; directions from F = q v × B with B out of the page: positive charges circle clockwise, negative anticlockwise. 5. Wire: F = ILB, perpendicular to both by L × B.

## 2.
(a) Source r′ = R(cos φ′, sin φ′, 0), dl = R(−sin φ′, cos φ′, 0) dφ′, r − r′ = (−R cos φ′, −R sin φ′, z); compute dl × (r − r′), integrate over φ′; the x and y components integrate cos φ′ and sin φ′ to zero. Expand (1 + R²/z²)^{−3/2}. (b) dB = μ₀(n dz′)IR²/(2(R² + (z − z′)²)^{3/2}); with u = arctan((z − z′)/R), integral of cos u du. Limits for z = 0 and z = L/2 as L/R → ∞. (c) F = I L × B on each side. The sides of length b: forces along ±(rotation axis), equal, opposite and collinear, so no force and no torque. The sides of length a: forces ±IaB perpendicular to B, separated perpendicular to their line of action by b sin θ: a couple IaB · b sin θ. Its direction rotates n̂ toward B, i.e. along n̂ × B.

## 3.
emf = vBL; I = emf/R; P = emf · I; F = ILB, direction from I L × B with the current direction fixed by the emf; a = F/m; Δv per orbit = a × 5600 s. Mechanical power F v = (vBL/R)LB · v = (vBL)²/R = emf²/R. Reversed current reverses the force: thrust, at the cost of battery power.

## 4.
(a) `u, B = field_along_orbit(400e3, 90.0)`; |B| min, max, mean; torque m_req|B| at the extremes. (b) Loop over N and d: `Torquer(N, 0.09, d, m_req/(N*0.09**2))`; print `I, resistance, power, inductance, time_constant`. (c) `saturation_rate(m_req, 2e4, Bmean)` against 5°/s; `detumble_time_constant`; `bdot_detumble(0.002, radians(5), Bmean, 2e4, m_req, 3000)`; first time ω < 0.1°/s; `momentum_dump_time`.

## 5.
(a) Fields of the two loops summed at (0, 0, z) for z = 0, 0.1R, 0.2R; deviation B(z)/B(0) − 1. (b) With f(z) = (R² + z²)^{−3/2}, B(z) ∝ f(z − d/2) + f(z + d/2); even in z; f″ at d/2 vanishes when d = R; the fourth derivative gives the coefficient; compare with −(144/125)(z/R)⁴. (c) B(0) = 2 × μ₀IR²/(2(R² + R²/4)^{3/2}) = (4/5)^{3/2}μ₀I/R.

## 6.
(a) |B| at r = 42 164 km from `earth_dipole_field`; compare the torque with solar-pressure torques. (b) μ₀m/(2πz³) on the axis; compare with 25–50 μT; options: duty-cycle the coil and sample the magnetometer while it is off, or subtract the known coil field. (c) τ = Σmᵢ × B is always ⊥ B. As B rotates in the body frame, the missing direction moves. (d) Orbital motion through the dipole field (once per orbit), and magnetometer noise / the torquer's own field; the noise floor sets the residual rate.
