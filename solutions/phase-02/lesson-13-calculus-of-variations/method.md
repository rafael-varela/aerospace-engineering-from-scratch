# P02.13  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. r = x_f/y_f; f(θ) = (θ − sin θ)/(1 − cos θ) − r; f(π) = π/2 − r tells you the side; bisect on (π, 2π) or (0, π). 2. k = y_f/(1 − cos θ_f); T = θ_f√(k/g); T_ramp = √(x_f² + y_f²)·√(2/(g y_f)). 3. Dip iff θ_f > π; lowest point y = 2k at x = πk. 4. `travel_time(*brachistochrone(x_f, y_f)[:2], g)` minus (7).

## 2.
(a) J[y + εη]; d/dε at 0; integrate the F_{y′}η′ term by parts; boundary term vanishes by η(a) = η(b) = 0; fundamental lemma. (b) d/dx(F − y′F_{y′}) = F_x + F_y y′ + F_{y′}y″ − y″F_{y′} − y′ d/dx F_{y′} = F_x + y′(F_y − d/dx F_{y′}); the bracket is zero on extremals; F_x remains if present. (c) F_{y′} = y′/[√(2gy)√(1 + y′²)]; form F − y′F_{y′}; square; rename the constant; then y′ = √((2k − y)/y) = cot(θ/2) with y = k(1 − cos θ); dx = dy/y′ = k(1 − cos θ) dθ; integrate. (d) sin θ_path = dx/ds = (1 − cos θ)/√(2(1 − cos θ)) = sin(θ/2); v = √(2gk(1 − cos θ)) = 2√(gk) sin(θ/2); divide.

## 3.
(a) Apex: θ = π/2, c = c₀ + ah; ground: c = c₀; equate sin θ/c. Circle geometry: with centre at height −c₀/a and radius R, sin θ = (y + c₀/a)/R, so sin θ/c = 1/(aR), constant; at the apex y = h, R = (c₀ + ah)/a. (b) Symmetric arc: centre at (D/2, −c₀/a), R² = (D/2)² + (c₀/a)²; h = R − c₀/a; elevation = 90° − θ₀. (c) Parametrise the arc by angle from the centre between the two ground points; `path_time(x, y, lambda x, y: c0 + a*y)`; compare with `path_time` on y = 0. (d) `snell_invariant` on the arc; on the ground path c = c₀ and θ = 90° so sin θ/c is constant *along it*, but an extremal must also satisfy the Euler–Lagrange equation, which contains c′(y) = a ≠ 0: the ray curves wherever the gradient is non-zero; the ground path is a boundary case (a path along the edge of the medium), not an interior extremal.

## 4.
(a) `ritz_minimum_time(x_f, y_f, n, g)` for n = 1, 2, 4, 8; excess = T_r/T − 1. (b) Monotone because each family contains the previous (set the new coefficient to zero); upper bound because the minimum over a subset is ≥ the minimum over the whole set. (c) `(y[1] − y[0])/(x[1] − x[0])` of the Ritz path against the cycloid's infinite slope. (d) y = y_f x/x_f + a₁x^{2/3}(1 − x/x_f) + a₂ sin(πx/x_f); Nelder–Mead on `travel_time`.

## 5.
(a) `perturbed_time_increase(x, y, eps, g)` for both paths; exponent = log(ΔT₁/ΔT₂)/log 2. (b) εδJ dominates unless δJ = 0; then ½ε²δ²J. (c) x_f/y_f against π/2; `cycloid_parameter`; depth 2k, position πk. (d) Speed is bought by depth; a deeper middle gives a faster average speed even though some of it is repaid at the end.

## 6.
(a) Friction: work μmg × (horizontal run) against mg y_f; point mass: rescales g only; rest: shifts the cusp; planar: nothing for a guided part. (b) Replace the vertical start by a straight drop at the friction angle, then join the cycloid. (c) Two tracks, same end points, a ball, light gates; predict the ratio from (7) and the ramp formula; a ratio below about 1.15 for the lesson's geometry means the model's assumptions have failed. (d) Stationary, not minimal; second variation, competing paths and the Ritz bound were the rest.
