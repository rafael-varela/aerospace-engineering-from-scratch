# P04.11  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
EIRP = 50 × 10^{1.3}; I = EIRP/(4πd²) with d = 20 200 km; E₀ = `e_from_intensity(I)`; λ = c/1575.42 MHz; A_e = λ²/(4π); P_r = I A_e; dBm = 10 log₁₀(P_r/10⁻³). N = kTB with T = 290 K, B = 2 MHz; compare in dB.

## 2.
(a) ∇ × (∇ × E) = −∂(∇ × B)/∂t = −μ₀ε₀ ∂²E/∂t²; with ∇·E = 0 this is ∇²E = μ₀ε₀ ∂²E/∂t²; c = 1/√(μ₀ε₀). For the plane wave: ∂B_y/∂t = −∂E_x/∂z = −kE₀ sin(kz − ωt) → B_y = (k/ω)E₀ cos(kz − ωt). (b) P = 2πU_max∫sin³θ dθ = 8πU_max/3; D = 4πU_max/P. Same for sin⁴θ with ∫sin⁵θ dθ. (c) I_r = P_t G_t/(4πd²); P_r = I_r A_e,r = I_r G_r λ²/(4π). Assumptions: far field, free space, matched polarisation and impedance, gains taken along the line joining the antennas. (d) Triangle O (centre), S (station), P (satellite): OS = R, OP = R + h, angle at S = 90° + e; (R + h)² = R² + d² + 2Rd sin e; solve the quadratic for d > 0.

## 3.
A = (10 + s)², m = 10 + s/3; `sail_acceleration(A, m)`; g_sun = GM/(1 au)²; ratio. The force along the normal is ∝ cos²α; its transverse component ∝ cos²α sin α; d/dα = 0 gives tan²α = 1/2. Transverse acceleration a(α*) sin α*; Δv in a year ≈ a_t × 3.156 × 10⁷ s (ignoring the change of orbit).

## 4.
λ = c/2.25 GHz; G = 0.55(πD/λ)²; A_e = Gλ²/(4π); beamwidth ≈ λ/D rad. `Link(...)` as in the hint; `budget()['margin_dB']`; rate = 10⁶ × 10^{(margin − 3)/10}. UHF: `Link(P_t=2.0, bit_rate=1e6)`. Decompose the difference: ΔFSPL = 20 log₁₀(2.25/0.437), Δ(dish − Yagi), Δ(patch − dipole), Δlosses, ΔT_sys in dB.

## 5.
L = 0.25 + s/40; `f = lambda t: dipole_pattern(t, L)`; `directivity(f)`, `directivity_quad(f)`, `half_power_beamwidth(f)`; maximum by `argmax` on a fine grid. `radiation_resistance_dipole(L)`; R_feed = R/sin²(πL). Errors of `dipole_pattern_numeric(t, L, n)` on `t = np.linspace(0.1, π − 0.1, 13)`; orders log₂(e_n/e_2n).

## 6.
(a) `to_db(dish_gain(3.0, c/8.4e9))` − 2; beamwidth λ/D. (b) The received power follows the satellite pattern sin²θ as the dipole's axis sweeps through the line of sight: fades, nulls; crossed dipoles / turnstile, circular. (c) Rank: lower-noise receiver (T_sys), bigger ground antenna, better coding, more transmit power (constrained on the satellite). (d) Compare the margin with the uncertainty of the assumed terms.
