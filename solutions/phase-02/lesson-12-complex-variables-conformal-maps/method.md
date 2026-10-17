# P02.12  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. u_xx + u_yy from the two second derivatives. 2. v_y = u_x = 3x² − 3y² + s ⇒ v = 3x²y − y³ + sy + g(x); then v_x = 6xy + g′(x) = −u_y = 6xy ⇒ g′ = 0. 3. Recognise u + iv as a polynomial in z; differentiate; evaluate at 1 + i; call `complex_derivative(f, np.array([1 + 1j]))` and `cauchy_riemann_residual`. 4. Conjugating v flips the sign of one Cauchy–Riemann equation.

## 2.
(a) Put s = Re^{−iβ} in (4), multiply through by e^{−iβ}, write the bracket as a sine, solve for Γ. Clockwise-positive Γ (the + iΓ/2π log convention) gives Γ > 0 for α + β > 0. (b) On the circle s = Re^{iθ}: the same manipulation gives sin(θ − α) = −sin(α + β); two solutions in θ. (c) p from Bernoulli; F_x + iF_y = −∮ p n ds with n ds = −i dz; conjugate; replace |w|² dz̄ by w² dz using w̄ dz̄ = w dz on the surface. (d) w = Ue^{−iα} + iΓ/(2πz) + O(z⁻²); square; the residue of w² at infinity's "1/z" term; apply ∮ = 2πi Σ residues; take real and imaginary parts.

## 3.
(a) F = (iΓ/2π)[log(z − ih) − log(z + ih)]; w = (iΓ/2π)[1/(z − ih) − 1/(z + ih)]. (b) Put z = x real; combine over a common denominator; the result is real (u only), so v = 0; maximise |u| in x. (c) Velocity at ih from the image alone: w = −(iΓ/2π)/(2ih) = −Γ/(4πh): real, so horizontal; the sign gives the direction. (d) `circulation_and_flux(w, circle_path(1j*h, r))` and `circle_path(0, 2*h + something)`; the ccw integral of a clockwise vortex is −Γ.

## 4.
(a) `airfoil_circle`, `np.arctan2(y0, c - x0)`, `airfoil_coordinates`, `airfoil_geometry`. (b) `kutta_circulation`, `lift_per_span`, `lift_coefficient` with the measured chord. (c) C_l at 0 and 8° divided by the angle in radians; t/c from (a); α₀ = −β. (d) Repeat the pressure integral of the lesson's test (F = i·½ρU²∮C_p dz, lift = F_y cos α − F_x sin α).

## 5.
(a)–(b) As in problem 4 at α = s/4°; `np.argmax(np.abs(w))` for the location. (c) M = |w|_max/a; U_crit = 0.3a/(|w|_max/U). (d) Evaluate `airfoil_surface_velocity` at c(1 + 10⁻⁶); then reason about the ratio of two small numbers when Γ is perturbed by 1 %.

## 6.
(a) Viscosity (boundary layer, drag, stall), compressibility, finite span, finite trailing-edge angle; sign of each on lift. (b) Pressure recovery over the rear of the upper surface. (c) Stall: use measured section data (wind-tunnel polars) or a viscous code. (d) Agreement among three computations of the same model's lift is verification; comparison with a measurement is validation.
