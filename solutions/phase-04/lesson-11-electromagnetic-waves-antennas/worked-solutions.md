# P04.11  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe scripts that call it, given to the figures shown.

## 1. How weak is a navigation signal?

**Assumptions.** The transmitter power and gain are illustrative assumptions, not a system specification; free space; receiver at the sub-satellite point (d = 20 200 km).

(a) **EIRP = 50 × 10^{1.3} = 998 W**. I = 998/(4π(2.02 × 10⁷)²) = **1.95 × 10⁻¹³ W/m²**; E₀ = √(2η₀I) = **1.21 × 10⁻⁵ V/m** (12 μV/m).

(b) λ = 0.1903 m, **A_e = λ²/(4π) = 2.88 × 10⁻³ m²**; **P_r = 5.61 × 10⁻¹⁶ W = −122.5 dBm**.

(c) kTB = 1.38 × 10⁻²³ × 290 × 2 × 10⁶ = **8.0 × 10⁻¹⁵ W = −111.0 dBm**: the signal is **11.5 dB below the noise** in that bandwidth. Spread-spectrum correlation over many code chips recovers it.

**Check.** The decibel route gives the same: EIRP = 29.99 dBW = 59.99 dBm, FSPL = 20 log₁₀(4π × 2.02 × 10⁷/0.1903) = 182.5 dB, and an isotropic receiver adds 0 dBi: 59.99 − 182.5 = −122.5 dBm.

**Common mistakes.** Using the satellite's orbital radius instead of its altitude for the sub-satellite distance; forgetting that an isotropic antenna still has an aperture λ²/(4π).

## 2. Waves, directivity, Friis and geometry

(a) Curl of ∇ × E = −∂B/∂t: ∇ × (∇ × E) = −∂(∇ × B)/∂t = −μ₀ε₀ ∂²E/∂t². The left side is ∇(∇·E) − ∇²E = −∇²E since ∇·E = 0. So **∇²E = μ₀ε₀ ∂²E/∂t²**, a wave equation with speed **c = 1/√(μ₀ε₀)**. For E_x = E₀ cos(kz − ωt): (∇ × E)_y = ∂E_x/∂z = −kE₀ sin(kz − ωt) = −∂B_y/∂t, so B_y = (k/ω)E₀ cos(kz − ωt) = **E_x/c** (a static part is excluded by requiring a pure wave).

(b) U = U_max sin²θ: P = 2πU_max∫₀^π sin³θ dθ = 2πU_max × 4/3, so **D = 4πU_max/P = 3/2**. In general D = 4π/(2π∫sinⁿ⁺¹θ dθ) = 2/∫sinⁿ⁺¹θ dθ. For n = 4, ∫sin⁵θ dθ = 16/15, **D = 15/8 = 1.875** (2.73 dBi): a narrower doughnut has higher directivity.

(c) An isotropic transmitter of power P_t gives P_t/(4πd²); with gain G_t toward the receiver, I = G_tP_t/(4πd²). The receiver delivers P_r = IA_e,r = I G_rλ²/(4π). So **P_r = P_tG_tG_r(λ/(4πd))²**. Assumptions: far field of both antennas; free space (no absorption, no reflections); gains taken in the direction of the line between them; polarisations matched; receiver impedance-matched (otherwise multiply by the mismatch and polarisation factors).

(d) Triangle with Earth's centre O, station S, satellite P: OS = R, OP = R + h, SP = d; the angle at S between SO and SP is 90° + e. Law of cosines: (R + h)² = R² + d² − 2Rd cos(90° + e) = R² + d² + 2Rd sin e. Quadratic in d: **d = −R sin e + √(R² sin²e + (R + h)² − R²) = √((R + h)² − R² cos²e) − R sin e**. Checks: e = 90° gives d = h; e = 0 gives √((R + h)² − R²), the horizon distance.

**Common mistakes.** In (b), forgetting the sin θ weight of the solid angle; in (d), taking the angle at S as e rather than 90° + e.

## 3. A solar-sail cargo tug (s = 15: A = 625 m², m = 15 kg)

(a) Facing the Sun: **a = 2 × 1361 × 625/(15 × c) = 3.78 × 10⁻⁴ m/s²**, **F = 5.67 mN**. Sun's gravity at 1 au: GM/r² = **5.93 × 10⁻³ m/s²**; ratio **6.4 %**. Both the intensity and gravity fall as 1/r², so their ratio, the "lightness number", is the same at every distance; it depends only on the sail's area per mass.

(b) The force is along the sail normal with magnitude ∝ cos²α; its component perpendicular to the sunline is ∝ cos²α sin α. Setting d/dα(cos²α sin α) = cos³α − 2cos α sin²α = 0 gives **tan²α = 1/2, α = 35.26°**. Transverse acceleration there: 3.78 × 10⁻⁴ × cos²(35.26°) × sin(35.26°) = **1.46 × 10⁻⁴ m/s²** (cos² = 2/3, sin = 1/√3: a factor 0.385).

(c) Δv ≈ 1.46 × 10⁻⁴ × 3.156 × 10⁷ s ≈ **4.6 km/s per year**, a rough figure since the orbit and the Sun distance change as it accumulates. The sail spends no propellant, so the rocket equation's exhaust velocity is effectively infinite; what limits it is the area-to-mass ratio (thin films, booms, deployment), degradation of the reflective coating, the 1/r² fall-off far from the Sun, and the fact that the force can never point toward the Sun.

**Check.** At 1 au the pressure on a mirror is 9.08 μPa; × 625 m² = 5.67 mN.

**Common mistakes.** Using cos α instead of cos²α (forgetting that the projected area and the momentum component each contribute one cosine); expecting a sail to slow down an object falling toward the Sun at full strength.

## 4. An S-band downlink (s = 15: D = 1.0 m)

(a) λ = 0.1332 m; **G = 0.55(π × 1.0/0.1332)² = 24.9 dBi**; **A_e = 0.432 m²** (55 % of the 0.785 m² dish area); beamwidth ≈ λ/D = **7.6°**. A satellite at 500 km crosses the sky at up to about 1°/s near the zenith, so a 7.6° beam must track it continuously; a Yagi with 14 dBi has a beam of about 40° and forgives more.

(b) At 10° (1695 km): FSPL 164.1 dB, P_r = −104.2 dBm, N₀ = −173.8 dBm/Hz (300 K), E_b/N₀ = 9.6 dB at 1 Mbit/s, **margin −0.38 dB**: it does not close at 1 Mbit/s. **Highest rate with 3 dB margin: 459 kbit/s.**

(c) The lesson's UHF link at 2 W and 1 Mbit/s has a margin of **−5.07 dB**. The S-band link is 4.7 dB better although its path loss is **14.2 dB** higher (Friis: the isotropic aperture shrinks as λ²). It wins back 10.9 dB with the dish over the Yagi (gain ∝ 1/λ² at fixed size), 3.9 dB with the patch over the dipole, 2 dB on losses (no polarisation mismatch) and 2.2 dB on noise temperature: 10.9 + 3.9 + 2.0 + 2.2 − 14.2 = 4.7 dB. The frequency costs nothing when the receiving antenna has a fixed area; it costs 20 log₁₀(f) only for fixed gain.

**Common mistakes.** Using the UHF noise temperature at S band; quoting the dish's physical area as its aperture.

## 5. A dipole of your own length (s = 15: L = 0.625λ)

(a) **D = 1.739 (2.40 dBi)** by both the trapezoid rule and `quad` (agreeing to 10⁻¹³); **HPBW 71.5°**; maximum broadside (90°). Between the half-wave (1.641) and the full-wave (2.411) values, as expected.

(b) R_rad referred to I₀ = **131.9 Ω**; the feed current is I₀ sin(kL/2) = I₀ sin(0.625π) = 0.924I₀, so **R_feed = 154.5 Ω**. A direct 50 Ω feed would reflect a large fraction of the power (mismatch factor 1 − ((154.5 − 50)/(154.5 + 50))² = 0.74, losing 1.3 dB, ignoring the reactance, which is also non-zero off resonance); a matching network or a different length is needed.

(c)

| elements | max error | order |
|---|---|---|
| 25 | 6.71 × 10⁻⁴ | — |
| 50 | 9.83 × 10⁻⁴ | −0.55 |
| 100 | 2.46 × 10⁻⁴ | 2.00 |
| 200 | 6.14 × 10⁻⁵ | 2.00 |
| 400 | 1.54 × 10⁻⁵ | 2.00 |

The first pair is pre-asymptotic: with 25 elements the current's kink at the centre and the phase across each element are poorly resolved, and the error happens to be small by cancellation at the 13 test angles; from 50 elements on the midpoint rule's h² law holds exactly.

**Check.** D of a 0.625λ dipole must lie between the 0.5λ and 1λ values and its beam between 78° and 48°: 71.5° does.

**Common mistakes.** Using the half-wave formula for all lengths; quoting R_rad referred to I₀ as the feed resistance.

## 6. Antennas on spacecraft

(a) The high-gain dish at 8.4 GHz has 0.55(π × 3/0.0357)² = **45.8 dBi**, about **44 dB** more than a 2 dBi low-gain antenna, so it carries 25 000 times the data rate; but its beam is 0.7° wide and must be pointed at the Earth. The low-gain antenna works in any attitude: it is the link for safe mode, early operations and emergencies, when pointing cannot be trusted.

(b) The received power follows the dipole's sin²θ pattern as the tumble sweeps its axis through the line of sight, and the linear polarisation rotates against the station's antenna: deep fades, twice per rotation, to nothing at the nulls. A crossed-dipole (turnstile) antenna fed in quadrature is circularly polarised and nearly omnidirectional: fades of a few decibels instead of nulls.

(c) Doubling the rate needs 3 dB. Cheapest per decibel on the ground: a low-noise amplifier at the antenna (T_sys from 500 to 250 K is 3 dB); then better coding (a few decibels of required E_b/N₀ for the price of software); then a larger antenna (3 dB needs 1.41 times the dish diameter). Transmit power on the satellite is the most expensive: it costs power, heat and mass where all three are scarce.

(d) A 1 dB margin is less than the uncertainty of any one assumed term (T_sys, antenna gain toward the station, pointing loss), each of which is ±2 dB or worse. The budget says the link *might* close. Either reduce the uncertainty (measure the antenna pattern, the receiver's noise), add margin (lower rate, bigger dish), or adopt an adaptive rate; 3–6 dB of margin is the usual design target for the stated uncertainties.

## Check-yourself answers

1. μ₀ and ε₀: c = 1/√(μ₀ε₀); E/B = c.
2. A quarter, 340 W/m²; pressure 2I/c = 2.27 μPa.
3. The field of each element is ∝ sin θ, zero along the wire, because an accelerating charge does not radiate along its line of acceleration. D = 1.5, 1.76 dBi.
4. Its peak intensity is 1.64 times that of an isotropic antenna radiating the same power; R_rad = 73 Ω.
5. λ²/(4π) = 0.0375 m² at 437 MHz; it scales as λ², 0.0014 m² at 2.25 GHz.
6. P_r = P_t + G_t + G_r − 20 log₁₀(4πd/λ) (dB); 6.02 dB per doubling.
7. Half, 3 dB; because the received power then no longer depends on the satellite's orientation, and a fixed 3 dB is cheaper than fades to zero.
8. The margin above 3 dB is 9 dB, a factor 7.9 in power, so about 7.9 × 9600 ≈ 76 000 bit/s; E_b/N₀ falls as 1/R_b at fixed received power.
