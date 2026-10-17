# P01.09  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer, and for every gradient say which variables it is taken with respect to and in what units.

Several problems use your personal seed **s**: the day of the month you were born (1–31). $\rho = 1.225$ kg/m³ unless stated.

## 1. (Warm-up: minimum drag at fixed weight)

A fixed-weight aircraft has $W = 12\,000 + 100s$ N, $S = 14$ m², $C_{D0} = 0.028$, $e = 0.8$ and a wing of aspect ratio 8 that is not going to change. Model it as `Aircraft(S=14, cd0=0.028, e=0.8, W0=W, c_w=0.0)`.

Find: (a) $K$, the minimum-drag speed $V_{md}$, $D_{\min}$ and $(L/D)_{\max}$; (b) the lift coefficient at $V_{md}$ and a check that it equals $\sqrt{C_{D0}/K}$; (c) at $V = 1.2V_{md}$, the drag, its parasite and induced parts, and the gradient $\nabla D = (\partial D/\partial V, \partial D/\partial\mathrm{AR})$ with units; (d) the directional derivative of $D$ at that point along the direction $(1, 0.1)$ in the $(V, \mathrm{AR})$ plane, and the largest directional derivative available there. Check everything against `minimum_drag_speed`, `drag`, `drag_gradient` and `directional_derivative`. Answer form: numbers with units, four significant figures.

## 2. (Derivation: the drag landscape in closed form)

For $D(V, \mathrm{AR})$ of equation (1) with $K = 1/(\pi e\,\mathrm{AR})$ and $W = W_0 + c_w\mathrm{AR}^2$:

(a) Derive both components of $\nabla D$, equations (2) and (3), naming the variable held fixed in each.
(b) Show that at $\partial D/\partial V = 0$ the parasite and induced drags are equal and $D_{\min}(\mathrm{AR}) = 2W\sqrt{KC_{D0}}$.
(c) Minimise $D_{\min}(\mathrm{AR})$ over AR and show $\mathrm{AR}^* = \sqrt{W_0/3c_w}$; show the wing penalty $c_w\mathrm{AR}^{*2}$ is exactly $W_0/3$.
(d) Show that at the optimum $D_{VV} = 4\rho SC_{D0}$, $D_{V,\mathrm{AR}} = 0$ and $D_{\mathrm{AR,AR}} = 6c_w\sqrt{K^*C_{D0}}$, and evaluate all three for the lesson's aircraft. Apply the second-derivative test.
(e) In the $(V, C_L)$ plane, set up $\nabla D = \lambda\nabla L$ with $L = W$ and derive $C_L^* = \sqrt{C_{D0}/K}$ and $\lambda = 2\sqrt{KC_{D0}} = 1/(L/D)_{\max}$. Then show directly, by differentiating $D_{\min} = 2W\sqrt{KC_{D0}}$ with respect to $W$, that $\lambda = \mathrm{d}D_{\min}/\mathrm{d}W$.

Answer form: derivations; (d) with numbers.

## 3. (Unfamiliar scenario: a glider hunting a thermal)

A thermal's updraft near the ground is modelled as $w(x, y) = w_0\exp[-(x^2 + y^2)/R^2]$ m/s with $w_0 = 3$ m/s and $R = 400$ m, $x$ east and $y$ north of the core. A glider is at $(x, y) = (300, 200)$ m flying on a heading of north-east.

(a) Compute $w$, $\nabla w$ and $|\nabla w|$ at the glider, with units. In which compass direction does the updraft increase fastest, and at what rate per metre?
(b) Compute the directional derivative of $w$ along the glider's current heading. Is it climbing into stronger or weaker lift?
(c) Which way, and by how many degrees, should the pilot turn to follow the steepest ascent? Why is this not the same as "turn toward the core"? (It is, here; explain why, and give a field for which it would not be.)
(d) Compute the Hessian of $w$ at the core and at the glider. Classify both points. Why is the glider's position neither a maximum nor a critical point and yet has a Hessian with eigenvalues of opposite sign?
(e) A real pilot cannot measure $\nabla w$; the variometer reads $w$ alone. Describe, in one paragraph, how circling lets the pilot estimate the direction of the gradient from a scalar reading, and what the directional derivative has to do with it.

Answer form: numbers with units; (c) and (e) in words.

## 4. (Uses your code: step size, scaling and conditioning)

Using `gradient_descent` from $(V, \mathrm{AR}) = (35, 6)$ on the lesson's aircraft with `tol=1e-6`:

(a) Record the iteration count for fixed steps $s = 10^{-4}, 10^{-3}, 10^{-2}, 0.03, 0.1, 0.25, 0.5$ in raw variables. Plot count against $s$ on log–log axes.
(b) Repeat with `scaled_drag` (variables $V/50$, $\mathrm{AR}/10$, drag in kN) for $s = 0.05, 0.1, 0.2, 0.3, 0.4, 0.5$ and `tol=1e-9`.
(c) Compute the Hessian at the optimum in both scalings with `hessian_fd`, its eigenvalues, the condition number $\lambda_{\max}/\lambda_{\min}$ and the stability limit $2/\lambda_{\max}$. Mark the limit on each plot. Explain why the count first falls as $1/s$, and what happens at and above the limit given that the code halves the step when the function fails to decrease.
(d) The scaled problem has the *worse* condition number. Explain why scaling was still worth doing, and propose a scaling (different reference values) that makes the scaled Hessian's condition number exactly 1. Verify it.
(e) Run `newton_minimize` from the same start and report the sequence of $|\nabla f|$ values. Fit the convergence order from the last three.

Answer form: two tables, two plots, the Hessian numbers, and a paragraph for (c)–(d).

## 5. (Judgment: is the optimum worth having?)

The lesson's optimum is $\mathrm{AR}^* = 10.41$, $V^* = 46.77$ m/s, $D^* = 1071.7$ N. Suppose the wing-weight coefficient $c_w$ is uncertain by ±20 % and $e$ by ±10 %.

(a) Propagate both uncertainties into $\mathrm{AR}^*$ and $D^*$, using the closed forms (the log-derivative trick of `P01.08` Problem 3 is enough).
(b) Evaluate $D_{\min}$ at $\mathrm{AR} = 8$ and $\mathrm{AR} = 13$ and compare with $D^*$. Given (a), is there any engineering content in the second decimal place of $\mathrm{AR}^*$?
(c) The structures group says each unit of AR above 9 adds €4000 to the wing's cost. Each newton of cruise drag costs €15 over the aircraft's life in fuel. Where is the cost optimum? Set it up as a one-variable problem and solve it by Newton's method or `scipy.optimize.minimize_scalar`.
(d) Write five sentences for a design review: the recommended AR, the two numbers that drive it, the one uncertainty that matters most, and the measurement that would reduce it.

Answer form: numbers for (a)–(c), prose for (d).

## 6. (Unfamiliar scenario: a hangar door as a constraint)

The aircraft of the lesson must fit a hangar with a 10 m door, so the span is limited to $b \le 10$ m; with $S = 16$ m², $\mathrm{AR} = b^2/S \le 6.25$.

(a) Compute the constrained optimum: $V_{md}$ and $D_{\min}$ at $\mathrm{AR} = 6.25$, and the weight. By how many newtons does the door cost the aircraft in drag, compared with $D^*$?
(b) Compute $\mathrm{d}D_{\min}/\mathrm{d}\mathrm{AR}$ at the limit, by finite difference and from the closed form of Problem 2(c), and convert it to newtons per metre of span.
(c) Set the problem up as $\min D_{\min}(\mathrm{AR})$ subject to $g(\mathrm{AR}) = \mathrm{AR} - 6.25 \le 0$, write the Lagrange condition at the active constraint, and identify the multiplier with (b). What is the sign of the multiplier, and what does that sign mean?
(d) A folding wing-tip adds 0.5 m of span per side (so $b = 11$ m) at a cost of 150 N of weight ($W_0 \to 13\,150$ N). Using only the two multipliers you now have ($\lambda_W = \mathrm{d}D/\mathrm{d}W$ and $\lambda_b = \mathrm{d}D/\mathrm{d}b$), estimate whether it pays. Then compute the exact answer and compare.

Answer form: numbers with units; (c) a short derivation.
