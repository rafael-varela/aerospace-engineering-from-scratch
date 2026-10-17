# P10.03  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15. Numbers come from `solution.py` and the probe scripts that call it, given to the figures shown.

## 1. A cylinder with circulation (s = 15: g = 0.25)

**Assumptions.** Steady, inviscid, incompressible 2-D flow; α = 0; clockwise Γ.

(a) Γ = 0.25 × 4π × 0.5 × 10 = **15.708 m²/s**. sin θ_s = −g = −0.25, so the stagnation points are at **θ = −14.48° and 194.48°**, both on the lower half, symmetric about the vertical.

(b) Cp = 1 − (2 sin θ + 2g)². The minimum is at the top, θ = 90°: **Cp = 1 − 2.5² = −5.25**; the maximum, **Cp = 1**, at the two stagnation points.

(c) L′ = ρVΓ = 1.225 × 10 × 15.708 = **192.4 N/m** (upward). On the diameter, c_l = L′/(½ρV² · 2a) = 2Γ/(V · 2a) = **3.1416** (= 4πg). The pressure integration with 512 points returns the same value to 10⁻¹², and zero drag.

(d) At g = 0 the stagnation points are fore and aft (0°, 180°). As g grows they slide down and together; at **g = 1** they merge at the bottom (θ = −90°); for g > 1 a single stagnation point lifts off the body into the flow below it, inside a closed region of fluid circulating with the cylinder.

**Check.** The lift is linear in Γ and independent of a at fixed Γ, as Kutta–Joukowski says.

**Common mistakes.** Putting the stagnation points on the upper half (the sign of Γ); using the radius instead of the diameter in c_l.

## 2. The Kutta condition fixes the lift

(a) With ζ′ = ae^{iθ}, (dw/dζ)e^{iθ} = u_r − iu_θ. Multiplying, V(e^{i(θ−α)} − e^{−i(θ−α)}) + iΓ/(2πa) = 2iV sin(θ − α) + iΓ/(2πa). The real part u_r is zero (the circle is a streamline) and **u_θ = −2V sin(θ − α) − Γ/(2πa)**.

(b) The airfoil-plane velocity is w_z = w_ζ/(dz/dζ). At ζ = b, dz/dζ = 0 (the map folds the circle into a sharp edge), so w_z is infinite (the flow turns round the sharp edge) unless w_ζ(b) = 0 as well, in which case the ratio is finite. ζ = b lies on the circle at θ = −β (b − μ = ae^{−iβ}). Setting u_θ(−β) = 0: −2V sin(−β − α) = Γ/(2πa), so **Γ = 4πVa sin(α + β)**.

(c) **c_l = ρVΓ/(½ρV²c) = 8πa sin(α + β)/c.** For a thin symmetric section the circle passes through ±b, a = b, β = 0, and the chord is the distance from −2b to 2b: c = 4b, so c_l = 2π sin α ≈ 2πα. **β is the zero-lift angle** (α_L0 = −β): moving the centre up (camber) increases β. **4a/c is the thickness factor**: a thick section has a > c/4 and a lift slope above 2π (1 + 0.77 t/c for Joukowski sections, §4).

(d) On the cylinder, p − p∞ = ½ρV²[1 − (2 sin θ + Γ/(2πaV))²] (α = 0). The force is −∮(p − p∞) n a dθ with n = (cos θ, sin θ). The lift is L′ = −a∮(p − p∞) sin θ dθ. Expanding the square, only the cross term 4 sin θ · Γ/(2πaV) multiplies sin θ into a non-zero integral (∮ sin θ = ∮ sin³ θ = 0, ∮ sin² θ = π): L′ = a · ½ρV² · (4Γ/(2πaV)) · π = **ρVΓ**. The drag integrand is odd about θ = 90°, so **D′ = 0**.

**Common mistakes.** Taking the real part in (a); forgetting that β must be measured from the centre to ζ = b, not from the origin.

## 3. A rotor sail

**Assumptions.** Potential flow with a circulation equal to the surface speed times the circumference; 2-D (no end effects); sea-level air.

(a) Γ = 2πaU_s = 2π × 2.0 × 30 = **377.0 m²/s**; **g = U_s/(2V) = 1.5**. Since g > 1, there is no stagnation point on the surface: a single one sits in the flow below the rotor.

(b) L′ = ρVΓ = 1.225 × 10 × 377.0 = **4618 N/m**; over 25 m, **115 kN**; on the diameter, **c_l = 4πg = 18.85**.

(c) Real rotors fall well short of this for two reasons that potential flow leaves out: (1) the boundary layer separates, especially on the side moving against the stream, so the fluid is not carried round at the surface speed and the effective circulation is much lower, and there is a pressure drag; (2) the rotor has ends, so trailing vortices at the top and bottom reduce the circulation near them and add induced drag (`P10.09`; end plates help). The weakest assumption is (a)'s circulation: the no-slip condition drags only the boundary layer round, not the outer flow.

**Common mistakes.** Treating the potential-flow c_l as a prediction; using a for the diameter.

## 4. A section like yours (s = 15: t = 0.14)

(a) μ = **−0.06754 + 0.04473i** (b = 1), **a = 1.06848**, **β = 2.3992°**, chord **3.8419**, chord-line inclination **0.0230°** (nose up). The fitted section has t = 0.14000 and camber 0.02006.

(b) α_L0 = −β + φ = **−2.376°** (chord line); slope **1.1112 × 2π** at the cruise point; the cruise c_l = 0.3309 at **α = 0.337°** (chord line).

(c) Kutta–Joukowski 0.3309; surface pressure (n = 1024) **0.33090003**, with pressure drag −3.1 × 10⁻⁷; Blasius, momentum and the circulation line integral on the ellipse **0.3309** to 15 figures. Four routes agree; the surface integral alone carries an n-dependent error, from the trailing-edge corner (Problem 5).

**Check.** Compared with the 12 % section of §4 of the lesson (slope 1.0939 × 2π, α_L0 = −2.378°), more thickness raised the slope and left α_L0 almost unchanged: thickness sets the slope, camber the zero-lift angle.

**Common mistakes.** Mixing α from the x-axis with α from the chord line; normalising the Blasius lift by the wrong chord.

## 5. How fast does pressure integration converge?

(a) Error of c_l from pressure against the exact value, α = 4°:

| n | τ = 0 (Joukowski) | τ = 16° |
|---|---|---|
| 16 | 6.13 × 10⁻³ | 2.23 × 10⁻³ |
| 32 | 3.12 × 10⁻⁴ | 3.71 × 10⁻³ |
| 64 | 3.58 × 10⁻⁷ | 1.47 × 10⁻⁴ |
| 128 | 2.5 × 10⁻¹³ | 2.02 × 10⁻⁶ |
| 256 | 2 × 10⁻¹⁶ | 5.43 × 10⁻⁷ |

Observed orders: τ = 0: 4.3, 9.8, 20.5, then round-off; τ = 16°: −0.7, 4.7, 6.2, 1.9.

(b) The integrand is a function of the circle angle θ, periodic over 2π. The trapezoid rule on a periodic function that is analytic in a strip converges exponentially, its error falling faster than any power of 1/n: that is the Joukowski case, where dz/dζ has a simple zero at ζ = b and everything stays analytic. For τ > 0, dz/dζ ∝ r^{n−1} with the non-integer power n − 1 = 0.911: the integrand has a branch point on the contour at the trailing edge, and the convergence becomes algebraic and irregular. Panel methods (`P10.05`) see the same corner.

**Common mistakes.** Concluding that more points always help in the same way; judging convergence from one pair of n.

## 6. Is the exact solution the right answer?

(a) Least-squares slopes of Ladson's Table I c_l from −4.25° to 6.05° (six points): **5.895 /rad = 0.938 × 2π at Re = 2.0 × 10⁶**; from −4.05° to 6.20°: **6.120 /rad = 0.974 × 2π at Re = 5.97 × 10⁶**. These are measurements (authentic, NASA TM-4074, wall-corrected by the report).

(b) The potential flow assumes no boundary layer. The real boundary layer thickens towards the trailing edge, more on the suction side, and its displacement thickness effectively decambers and thins the aft section, so the circulation that satisfies the "Kutta condition" of the real flow is smaller: the exact slope 1.094 × 2π is 12–17 % above the measurements. At higher Re the boundary layer is thinner, and the measured slope rises towards the inviscid one.

(c) In an inviscid flow any circulation is allowed; it is viscosity at the sharp edge that sheds a starting vortex until the flow leaves the edge smoothly, and by Kelvin's theorem the airfoil keeps the opposite circulation. The Kutta condition writes the end state of that viscous process as a rule for the inviscid flow. It works because the boundary layer is thin at high Re, so the outer flow is nearly the inviscid one with that circulation; `P10.06` adds the boundary layer back.

**Common mistakes.** Fitting into the non-linear region near stall; treating the 2π of thin-airfoil theory as the benchmark rather than the thick-section potential flow.

## Check-yourself answers

1. Γ = ∮ u · dl round a closed curve; in this phase positive clockwise for a lifting airfoil in a left-to-right stream, so that L′ = ρVΓ is positive.
2. L′ = ρV∞Γ, per unit span, perpendicular to the free stream; there is no drag in steady 2-D inviscid flow (d'Alembert).
3. Potential flow past a body allows any circulation; the Kutta condition picks the one for which the flow leaves the sharp trailing edge smoothly (finite velocity there).
4. Γ = 4πVa sin(α + β): the trailing-edge point of the circle (θ = −β) is made a stagnation point.
5. Because the airfoil-plane velocity is the circle-plane velocity divided by dz/dζ, which is zero at the trailing edge.
6. Its circulation (the same on every enclosing contour, by Stokes' theorem in an irrotational flow); a contour that does not enclose the body gives zero.
7. Lift-slope ratio to 2π: 1.0308 at 4 %, 1.0926 at 12 %, 1.1549 at 20 % (Joukowski), close to 1 + 0.77 t/c.
