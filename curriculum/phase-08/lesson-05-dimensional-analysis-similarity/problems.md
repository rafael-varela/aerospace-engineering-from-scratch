# P08.05  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical. For every dimensional analysis, list the variables with their dimensions in a table before you form any group.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

Shared numbers: ISA from `P08.01`; Sutherland's law for air; water at 20 °C ρ = 998.2 kg/m³, μ = 1.002 × 10⁻³ Pa·s.

## 1. (Warm-up: Reynolds and Mach numbers you will meet)

(a) The light aeroplane at take-off: chord 1.5 m, 30 m/s, ISA sea level. Compute Re and M.
(b) A sounding rocket of diameter 0.15 m and length 3 m flies at V = 200 + 10s m/s at 5000 m (ISA). Compute M and the Reynolds numbers based on diameter and on length.
(c) For each case, say which number you would expect to control the drag coefficient most, and why.

Answer form: two numbers; three numbers; three sentences.

## 2. (Derivation: the Pi theorem as linear algebra)

(a) Show that a product Π = Π_j q_j^{e_j} of n variables is dimensionless if and only if the exponent vector e is in the nullspace of the dimension matrix D. Deduce that the number of independent groups is n − rank D (`P01.14`).
(b) Show that if Π₁, …, Π_k is one basis of groups, any other basis consists of products of powers of them, and that a relation f(Π₁, …, Π_k) = 0 can be rewritten in the new basis. Illustrate with (Re, Nu, Pr) and (Re, Nu, St), St = Nu/(Re Pr).
(c) Apply the theorem to the force on a wing, F = f(ρ, V, c, μ, a), to obtain F = ½ρV²c² C(Re, M). Then show that if F must not depend on ρ (creeping flow, Re → 0) and the flow is incompressible (no a), F = const × μVc. Compare with Stokes's drag on a sphere, 3πμVd.
(d) The SVD gives the nullspace as orthonormal vectors that are not integers. Prove that they span the same space as the SymPy basis, using rank arguments only.

Answer form: two derivations; a derivation with an example; two derivations and a constant; a proof.

## 3. (Unfamiliar scenario: testing a model propeller)

A propeller's thrust depends on the air density ρ, the rotational speed n (rev/s), the diameter d, the flight speed V, the viscosity μ and the speed of sound a.

(a) Using the code (add `n` with dimensions T⁻¹ to a copy of `DIMENSIONS`), find the number of groups and form them with ρ, n, d as repeating variables. Name the groups (thrust coefficient, advance ratio, a Reynolds number, a tip Mach number). Check that the SymPy basis spans the same space even though it looks different.
(b) The light aeroplane's propeller has d = 1.8 m and turns at 2700 rpm in cruise at 60 m/s at 3000 m. Compute the advance ratio J = V/(nd) and the tip Mach number πnd/a.
(c) A 0.6 m model is tested in a sea-level tunnel at 288.15 K. Find the rotational speed and tunnel speed that match J and the tip Mach number. By what factor does the blade Reynolds number (based on the tip speed and a blade chord of 0.15 m full scale, scaled with the model) fall short?
(d) Which mismatch would worry you most for a propeller test, and why?

Answer form: a count, four groups and a check; two numbers; two speeds and a factor; a sentence.

## 4. (Numerical, uses your code: a test plan for your model)

Scale the light aeroplane's wing (chord 1.5 m, cruise at 60 m/s at 3000 m) by 1:k with k = 3 + s/10.

(a) Compute the full-scale Re and M, the model chord, the speed that matches Re in a sea-level atmospheric tunnel, and the Mach number it implies.
(b) Find the pressure that matches both at 288.15 K, the temperature that matches both at 1 atm and at 2 atm (`cryogenic_tunnel`), and the speed for a Re match in a water tunnel.
(c) Compute the dynamic pressure on the model in the full-scale flight, in the atmospheric Re-matched tunnel and in the pressurised tunnel. Which of the options would overload a model balance or the model?
(d) Choose one option and defend it in three sentences, including what you give up.

Answer form: four numbers; four numbers; three dynamic pressures and a sentence; a short argument.

## 5. (Numerical, uses your code: heat transfer and three routes to the nullspace)

The heat-transfer coefficient h (W/(m²·K), dimensions M T⁻³ Θ⁻¹) from a duct wall to the air in it depends on ρ, V, d, μ, the conductivity k and the specific heat c_p.

(a) Build the dimension matrix (add `h` to a copy of `DIMENSIONS`), find its rank and the number of groups.
(b) Find the groups by SymPy, by SVD and by repeating variables ρ, V, d, k. Show that all three are dimensionless and that they span the same space. Express each SymPy and each repeating-variable group in terms of Re, Nu = hd/k and Pr = μc_p/k.
(c) Give the singular values of D. What would a singular value of 10⁻¹⁴ have told you?
(d) Explain why the SVD basis is mathematically correct but useless for reporting results, and why the repeating-variable basis is the one engineers use.

Answer form: a rank and a count; three bases and three conversions; four numbers and a sentence; two sentences.

## 6. (Judgment: which number do you match?)

For each model test, say which dimensionless number or numbers must be matched, which can be relaxed, and what you would do about the mismatch: (a) the light aeroplane's wing at landing (30 m/s, flaps down) in a small atmospheric tunnel at 1:5; (b) an airliner wing at M = 0.82 in a transonic tunnel; (c) a seaplane hull planing on water at 1:10 in a towing tank; (d) a parachute canopy opening, at 1:4 in a vertical tunnel; (e) a launch vehicle's base heating; (f) wind loads on a launch tower.

Answer form: six labelled lines.
