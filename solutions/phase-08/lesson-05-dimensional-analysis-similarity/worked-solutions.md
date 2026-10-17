# P08.05  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe script that calls it, given to the figures shown.

## 1. Reynolds and Mach numbers you will meet (s = 15: V = 350 m/s)

(a) Sea level, μ = 1.7894 × 10⁻⁵ Pa·s: Re = 1.225 × 30 × 1.5/1.7894 × 10⁻⁵ = **3.08 × 10⁶**, M = **0.088**.

(b) At 5000 m: T = 255.65 K, ρ = 0.7361 kg/m³, μ = 1.628 × 10⁻⁵ Pa·s, a = 320.5 m/s. **M = 1.092**; Re_d = **2.37 × 10⁶**; Re_L = **4.75 × 10⁷**.

(c) The aeroplane at take-off: Re (M is far below 0.3); the wing's maximum lift, which sets the take-off speed, depends on the boundary layer and hence on Re. The rocket: M, because at M = 1.09 a bow shock and wave drag dominate the forebody drag; Re_L (4.7 × 10⁷, turbulent over almost the whole length) sets the skin-friction part, which is secondary near M = 1.

**Common mistakes.** Using the sea-level viscosity at 5 km; using the rocket's length for a nose-shock question.

## 2. The Pi theorem as linear algebra

(a) The dimension of q_j is Π_i [B_i]^{D_ij} with B the base dimensions. Then [Π] = Π_j Π_i [B_i]^{D_ij e_j} = Π_i [B_i]^{(De)_i}: **Π is dimensionless iff De = 0**. The set of such e is the nullspace of D, whose dimension by the rank–nullity theorem is **n − rank D**; any dimensionless group is a product of powers of the groups of a basis, so that is the number of independent groups.

(b) If e₁, …, e_k and e′₁, …, e′_k are two bases of the nullspace, e′_j = Σ_i C_ij e_i with C invertible; in variables, Π′_j = Π_i Π_i^{C_ij}. A relation f(Π₁, …, Π_k) = 0 becomes f′(Π′₁, …, Π′_k) = 0 by substituting the inverse. Example: St = h/(ρVc_p) = (hd/k)(μ/(ρVd))(k/(μc_p)) = Nu/(Re Pr): exponents of (Re, Nu, Pr) are (−1, 1, −1), and (Re, Nu, St) is a basis because the transformation is invertible.

(c) Repeating variables ρ, V, c give Π_F = F/(ρV²c²), Π_μ = μ/(ρVc) = 1/Re, Π_a = a/V = 1/M, so **F = ½ρV²c² C(Re, M)** (the ½ is a convention). Incompressible: F = ρV²c² g(Re). Independence of ρ: ρV²c² g(ρVc/μ) must not change with ρ, which requires g(Re) = K/Re: **F = KμVc**. Stokes's sphere: **K = 3π** with c = d; in coefficient form C_D = F/(½ρV²πd²/4) = 24/Re (test `test_stokes_drag_collapses_to_cd_re_equals_24`).

(d) The SVD's last n − r right singular vectors satisfy D v = 0 (their singular values are zero) and are orthonormal, hence independent: n − r independent vectors in a nullspace of dimension n − r span it. The SymPy basis spans the same nullspace; so the column spaces coincide, which `same_span` tests by rank(B₁) = rank(B₂) = rank([B₁ B₂]).

**Common mistakes.** Counting the groups as n minus the number of base dimensions (it is the rank, which can be smaller); believing that one basis is "the" answer.

## 3. Testing a model propeller

(a) Seven variables, rank 3: **4 groups**. With ρ, n, d repeating: **C_T = F/(ρn²d⁴)** (thrust coefficient), **J = V/(nd)** (advance ratio), **μ/(ρnd²)** (an inverse Reynolds number based on tip-speed scale), **a/(nd)** (an inverse tip Mach number, up to a factor π). SymPy returns the thrust coefficient, Fn²/(ρV⁴), Fρ/μ² and Fn²/(ρa⁴): different products of the same groups; `same_span` with the SVD basis is True.

(b) n = 45 rev/s; at 3000 m a = 328.6 m/s: **J = 60/(45 × 1.8) = 0.741**, **M_tip = π × 45 × 1.8/328.6 = 0.774** (the rotational tip speed; the helical tip Mach number is a little higher).

(c) At 288.15 K, a₀ = 340.29 m/s: n_m = 0.774 × 340.29/(π × 0.6) = **139.8 rev/s (8389 rpm)**; V_m = J n_m d_m = **62.1 m/s**. Blade Reynolds number (tip speed, chord 0.15 m full scale, 0.05 m model): 2.05 × 10⁶ full scale, 9.02 × 10⁵ model: **0.44 of full scale**.

(d) The Reynolds mismatch: J and M_tip fix the kinematics and the compressibility, but at 9 × 10⁵ the blade boundary layers may stay laminar further and separate earlier than at 2 × 10⁶, which changes the blade drag (power coefficient) more than the thrust. Trip the boundary layer, or run at two Reynolds numbers and extrapolate.

**Common mistakes.** Using rpm in the groups (n must be rev/s and consistent); matching J but not M_tip, which at 0.77 is too high to relax.

## 4. A test plan for your model (s = 15: k = 4.5)

(a) Full scale: **Re = 4.83 × 10⁶, M = 0.183**. Model chord **0.333 m**; the Re-matching speed in sea-level air is **211.7 m/s**, **M = 0.622**: more than three times the flight Mach number, in the range where compressibility changes the pressure distribution by several per cent (`P08.04` Problem 5).

(b) Pressurised at 288.15 K: **ρ = 4.173 kg/m³, p = 3.452 bar (3.41 atm)**, V = 62.1 m/s. Cryogenic at 1 atm: **T = 118.2 K**, V = 39.8 m/s; at 2 atm: **T = 192.5 K**, V = 50.8 m/s. Water tunnel: **14.55 m/s** (Mach irrelevant). The 118 K figure extrapolates Sutherland's law far below the temperatures it is normally used at, so it is indicative only; real cryogenic tunnels use nitrogen with its own property data.

(c) Full-scale flight: q = **1636 Pa**. Atmospheric Re-matched tunnel: **27 449 Pa** (16.8 times). Pressurised tunnel: **8057 Pa** (4.9 times). The atmospheric Re match puts 17 times the flight dynamic pressure on the model: its lift at a given C_L is 16.8/k² = 0.83 of the full-scale wing's although it is 4.5 times smaller, and in a geometrically similar structure the stresses scale with q alone, so the model's are about 17 times the aircraft's; a balance and a model designed for a few hundred newtons would be overloaded or deflect. Matching both numbers in a pressurised tunnel still loads the model at about five times flight q, which is why pressurised-tunnel models are made of steel.

(d) For a low-speed wing the defensible choice is the pressurised tunnel at about 3.5 bar: both numbers matched at moderate temperature, at the cost of a stiff, expensive model and a pressure vessel. An atmospheric test at the flight Mach number (Re about 1.4 × 10⁶ at 62 m/s) with transition trips and a Reynolds-number correction is cheaper and is what most low-speed projects do; it gives up the maximum-lift and drag accuracy that depend on Re.

**Common mistakes.** Believing the Mach mismatch is harmless because the full-scale flight is incompressible (the *model* is not); forgetting the model loads.

## 5. Heat transfer and three routes to the nullspace

(a) Seven variables, four base dimensions all present, **rank 4**, **3 groups**.

(b) SymPy: **Re = ρVd/μ**, **Nu = hd/k**, **St = h/(ρVc_p)**. SVD: three orthonormal vectors with max |D B| = 4.4 × 10⁻¹⁶. Repeating ρ, V, d, k: **Nu**, **μ/(ρVd) = 1/Re**, **ρVdc_p/k = Re Pr** (the Péclet number). All are dimensionless and `same_span` is True for every pair. Conversions: St = Nu/(Re Pr); Pe = Re Pr; 1/Re is Re⁻¹. So Nu = f(Re, Pr), the form in which every convective heat-transfer correlation is written.

(c) Singular values **5.639, 3.962, 0.609, 0.360**: all clearly non-zero, so rank 4 is unambiguous. A singular value of 10⁻¹⁴ would mean the matrix is rank deficient (one base dimension can be written in terms of the others for these variables, for example if no variable carried Θ independently), and there would be one more group than the count from the number of base dimensions.

(d) The SVD vectors are irrational combinations of all seven variables: correct, but no one can measure, tabulate or interpret "ρ^−0.2041 V^0.2659 …". The repeating-variable basis gives integer groups each of which isolates one physical ratio (one dependent variable per group, the rest made from a fixed reference set), which is what correlations and test reports need.

**Common mistakes.** Choosing repeating variables that do not include Θ (then k or c_p cannot be made dimensionless); thinking the three bases are three different answers.

## 6. Which number do you match?

(a) Landing wing, 1:5 atmospheric: Re matters (maximum lift) and M ≈ 0.09 can be relaxed up to about 0.3; Re cannot be matched (it would need 5 times the speed); use a pressurised or larger tunnel, or test at the highest Re available and extrapolate, with trips to fix transition.
(b) Transonic airliner wing: M first (shock position), then Re (shock–boundary-layer interaction); cryogenic pressurised tunnels exist for this; otherwise match M and correct for Re.
(c) Planing hull: Froude number V/√(gL) (wave making) must match; Re cannot be matched at the same time in water, so match Fr and correct the friction drag with a skin-friction line (the Froude method).
(d) Parachute opening: a reduced time scale (opening time × V/D) and the fabric's porosity and elastic properties matter more than Re, which is high in both; the mass ratio of canopy to displaced air must also scale.
(e) Base heating: Mach and Reynolds numbers of the external flow and the plume's stagnation temperature, pressure ratio and Prandtl number; no tunnel matches them all, so tests use hot-gas plumes and analysis.
(f) Launch tower: Re (vortex shedding off members), and the Strouhal number for vortex-induced vibration; usually tested at much lower Re with roughness, and checked against full-scale anemometry.

## Check-yourself answers

1. The number of independent dimensionless groups is the number of variables minus the rank of the dimension matrix.
2. Its columns are the base-dimension exponents of the variables; a dimensionless group's exponent vector is in its nullspace.
3. The SymPy basis is exact and rational (integers after scaling); the SVD basis is orthonormal floating-point; both span the same space, and either can be converted to the other by an invertible matrix.
4. Re = ρVL/μ, inertial to viscous forces; M = V/a, flow speed to sound speed (compressibility).
5. For a 1:k model in the same air, matching Re needs k times the speed, which raises M by k; both cannot be matched unless ρa/μ is raised by k.
6. Raising the pressure raises ρ at constant a and μ; lowering the temperature lowers μ (by Sutherland) and a, raising Re at fixed M; using a different fluid (water) changes ν by a large factor.
7. Dimensionally homogeneous equations keep the same form whatever the units, so a group computed in SI, cgs or imperial units has the same value; if it does not, a unit has been mixed up.
8. Partial similarity: matching the most important group (for example Mach for a transonic wing) and correcting for the others (Reynolds number by trips or extrapolation).
