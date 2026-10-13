# P01.09  Gradients, Directional Derivatives & Optimization

> The gradient points uphill. Everything else in optimization is deciding how far to walk the other way.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the gradient of the drag of an aircraft in level flight with respect to speed and aspect ratio, the directional derivative, the Hessian and the second-derivative test, the unconstrained optimum in closed form, and the Lagrange-multiplier solution of minimum drag subject to lift equals weight. 2. Implement and test steepest descent and Newton's method with analytic and finite-difference gradients and Hessians, and compare them with `scipy.optimize.minimize`. 3. Measure, by computational experiment, how the step size and the scaling of the variables change the number of iterations, and relate it to the Hessian's eigenvalues. 4. Decide, from the Lagrange multiplier, what a newton of extra weight costs in drag, and from the span constraint what a metre of span is worth. |
| **Prerequisites** | `P01.08` (partial derivatives, tangent plane), `P01.07` (dot product, projection), `P01.02` (Newton's method, optimization in one variable), `P00.09` (SymPy `solve`). Skill check: you can find the minimum of $AV^2 + B/V^2$ by setting a derivative to zero, and you know what the eigenvalues of a symmetric $2 \times 2$ matrix are. |
| **Engineering question** | At what cruise speed and with what wing aspect ratio does a light aircraft have the least drag, when a higher aspect ratio reduces induced drag but makes the wing heavier; and how much drag does each extra newton of weight, or each metre of span, buy or cost? |
| **Assumptions and limits** | Steady level flight (lift equals weight); the parabolic drag polar $C_D = C_{D0} + KC_L^2$ with $K = 1/(\pi e\,\mathrm{AR})$ and constant $C_{D0}$, $e$; incompressible flow; a wing-weight surrogate $W = W_0 + c_w\mathrm{AR}^2$ that stands in for a structural analysis; fixed wing area. The surrogate is the weakest link: real wing weight depends on span, loading, taper and materials, and $c_w$ here is chosen so the optimum lands in a plausible place, not derived. The drag polar fails above Mach 0.3, near stall, and for $C_{D0}$ that varies with Reynolds number. Minimum drag is the right objective for maximum glide range and for the maximum range of a jet only after a factor of $V$ is included; for a piston aircraft maximum range is at a different speed again (Phase 10). |
| **Analytical work** | $\nabla D$ by hand; directional derivatives and the steepest-descent direction; $V_{md}$ at fixed weight; $D_{\min}(\mathrm{AR})$ and $\mathrm{AR}^* = \sqrt{W_0/3c_w}$; the Hessian at the optimum in closed form ($D_{VV} = 4\rho SC_{D0}$, $D_{\mathrm{AR,AR}} = 6c_w\sqrt{K^*C_{D0}}$, zero cross term) and the second-derivative test; Lagrange multipliers for $\min D$ subject to $L = W$ in the $(V, C_L)$ plane, giving $C_L^* = \sqrt{C_{D0}/K}$ and $\lambda = 2\sqrt{KC_{D0}} = 1/(L/D)_{\max}$. |
| **Computational deliverable** | `starter/main.py`: `induced_factor`, `weight`, `drag`, `drag_gradient`, `gradient_fd`, `hessian_fd`, `directional_derivative`, `classify_critical_point`, `gradient_descent`, `newton_minimize`, `scaled_drag`, `minimum_drag_speed`, `optimum_analytic`, `lagrange_min_drag`, `lagrange_sympy`. Tests prove: parasite and induced parts scale as $V^2$ and $V^{-2}$; the analytic gradient equals central differences at second order; the gradient vanishes at $V_{md}$ where the two drag parts are equal; the directional derivative is maximal along the gradient and zero along the level curve; the analytic optimum is a critical point with a positive-definite Hessian and lower drag than fifty random neighbours; steepest descent and Newton converge to it, Newton quadratically; the Lagrange conditions hold with $\lambda = \mathrm{d}D/\mathrm{d}W$ by finite difference and the SymPy solution matches the closed forms. |
| **Practical activity** | **computational experiment**: iteration count of steepest descent against step size in raw and scaled variables, against the stability limit $2/\lambda_{\max}$ from the Hessian; and the drag cost of a span limit through the constrained optimum. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `drag_landscape.png`, the hand calculation of $\nabla D$ at (35 m/s, AR 6) and of the optimum with its comparison to the code, the step-size table with the Hessian eigenvalues, and the Lagrange interpretation in one paragraph. |
| **Assessment** | `assessments/phase-01/grading-rubric.md`; competencies 1–4. Pass: the gradient, the optimum and the Lagrange conditions derived with every held-fixed variable named; all public tests pass; the hand calculation agrees with the code to four figures; the step-size experiment reports the stability limit with a number and compares it to $2/\lambda_{\max}$. |
| **References** | `@openstax-calculus-volume-3` §4.6 "Directional Derivatives and the Gradient" (required-free); §4.7 "Maxima/Minima Problems" (required-free); §4.8 "Lagrange Multipliers" (required-free); `@openstax-calculus-volume-3` §4.4 "Tangent Planes and Linear Approximations" for the Hessian as the second-order term (supplementary-free). |
| **Extensions** | ★ Replace the fixed-step descent by a line search that minimises $f$ along $-\nabla f$ at each step and compare iteration counts. ★ Minimise $D/V$ (jet range) and $DV$ (power, piston endurance) instead of $D$ and show the three optimum speeds are $1.32\,V_{md}$, $V_{md}$ and $0.76\,V_{md}$. |
| **Time** | unmeasured |

## 1. The engineering problem

A light aircraft's drag in level flight has two parts. Parasite drag, from skin friction and pressure, grows with the square of the speed. Induced drag, the price of making lift with a finite wing, falls with the square of the speed and is smaller for a wing of higher **aspect ratio** (span squared over area), because a long slender wing sheds weaker tip vortices. There is therefore a speed of least drag, and the designer can push the drag lower still by stretching the wing. But a longer wing is heavier, and weight raises the lift required and with it the induced drag. Somewhere in the plane of (cruise speed, aspect ratio) the drag is least, and the design team wants that point, how flat the drag surface is around it, and two prices: how many newtons of drag an extra newton of weight costs, and how many a metre of span saves. Those prices decide whether a lighter battery or a longer wing is the better investment.

The mathematics is a function of two variables and its **gradient**, the vector of partial derivatives from `P01.08`. The gradient points in the direction of steepest ascent; walking against it is the simplest optimization algorithm there is, and the Hessian, the matrix of second partials, tells us whether a point where the gradient vanishes is a valley bottom, a hilltop or a saddle. Constraints such as "lift must equal weight" are handled by **Lagrange multipliers**, which not only find the constrained optimum but hand back the price of relaxing the constraint. This lesson is the one-variable optimization of `P01.02` grown up, and it is also a preview of Phase 03, where the same algorithms run on functions that cost hours to evaluate.

Success criteria, written before we start: the analytic gradient agrees with central differences to $10^{-7}$ and converges at second order; the optimum found by steepest descent, by Newton's method and by `scipy.optimize.minimize` agrees with the closed form to $10^{-5}$; the Hessian at the optimum is positive definite with the entries predicted in §2.5; and the Lagrange multiplier equals $\mathrm{d}D/\mathrm{d}W$ by finite difference to $10^{-6}$.

## 2. Understand and derive

### 2.1 The drag landscape

**Assumption 1 (parabolic drag polar):** $C_D = C_{D0} + KC_L^2$ with $K = 1/(\pi e\,\mathrm{AR})$, $e \approx 0.8$ the span-efficiency factor (Phase 10 derives $K$ from lifting-line theory; here it is given). **Assumption 2 (level flight):** $L = W$, so $C_L = 2W/(\rho V^2S)$. Then, with $q = \tfrac12\rho V^2$,

$$D(V, \mathrm{AR}) = qSC_{D0} + \frac{KW^2}{qS} = \underbrace{\tfrac12\rho V^2SC_{D0}}_{\text{parasite}} + \underbrace{\frac{2K(\mathrm{AR})\,W(\mathrm{AR})^2}{\rho V^2S}}_{\text{induced}}.$$   (1)

**Assumption 3 (weight surrogate):** $W(\mathrm{AR}) = W_0 + c_w\mathrm{AR}^2$. Our aircraft: $S = 16$ m², $C_{D0} = 0.025$, $e = 0.8$, $W_0 = 13\,000$ N, $c_w = 40$ N, $\rho = 1.225$ kg/m³. At $\mathrm{AR} = 10$ the wing penalty is 4000 N, a quarter of the total; a crude but not absurd share for a long wing on a small aircraft.

Units of (1): Pa · m² = N for the first term; N²/(Pa · m²) = N²/N = N for the second ✓. A **level curve** $D = $ const in the $(V, \mathrm{AR})$ plane is a contour of the drag landscape; `main()` draws them and they look like a tilted valley.

### 2.2 The gradient

The gradient is the vector of partials (`@openstax-calculus-volume-3` §4.6):

$$\nabla D = \left(\frac{\partial D}{\partial V},\ \frac{\partial D}{\partial\mathrm{AR}}\right).$$

Differentiate (1) with AR held fixed (so $K$ and $W$ are constants):

$$\frac{\partial D}{\partial V} = \rho VSC_{D0} - \frac{4KW^2}{\rho V^3S} \quad [\mathrm{N/(m/s)}].$$   (2)

Now with $V$ fixed; $K' = \mathrm{d}K/\mathrm{d}\mathrm{AR} = -K/\mathrm{AR}$ and $W' = 2c_w\mathrm{AR}$:

$$\frac{\partial D}{\partial\mathrm{AR}} = \frac{2}{\rho V^2S}\left[K'W^2 + 2KWW'\right] = \frac{K'W^2 + 2KWW'}{qS} \quad [\mathrm{N}].$$   (3)

The first bracketed term is negative (more AR, less induced drag per unit of lift squared), the second positive (more AR, more weight). **Hand calculation at $V = 35$ m/s, AR = 6:** $K = 1/(\pi \times 0.8 \times 6) = 0.06631$; $W = 13\,000 + 1440 = 14\,440$ N; $q = 750.3$ Pa; parasite $= 750.3 \times 16 \times 0.025 = 300.1$ N; induced $= 0.06631 \times 14\,440^2/(750.3 \times 16) = 1151.8$ N; $D = 1451.9$ N. $\partial D/\partial V = 1.225 \times 35 \times 16 \times 0.025 - 4 \times 0.06631 \times 14\,440^2/(1.225 \times 35^3 \times 16) = 17.15 - 65.82 = -48.67$ N/(m/s). $\partial D/\partial\mathrm{AR}$: $K' = -0.011052$, $W' = 480$ N; $[-0.011052 \times 14\,440^2 + 2 \times 0.06631 \times 14\,440 \times 480]/(750.3 \times 16) = [-2.3045 \times 10^6 + 9.193 \times 10^5]/12\,005 = -115.4$ N. So $\nabla D = (-48.67, -115.4)$: at this point both more speed and more aspect ratio reduce drag, and the drag falls fastest in the direction $-\nabla D/|\nabla D| = (0.389, 0.921)$, mostly toward higher AR.

The gradient has two geometric properties that the tests check (§4.6): it is **perpendicular to the level curve** through the point (moving along a contour changes nothing, so the rate of change in that direction is zero), and its length $|\nabla D|$ is the largest rate of change in any direction.

### 2.3 Directional derivatives

The rate of change of $D$ per unit distance in the direction of a unit vector $\hat{\mathbf u}$ is the **directional derivative**

$$D_{\hat{\mathbf u}}f = \nabla f \cdot \hat{\mathbf u} = |\nabla f|\cos\theta,$$   (4)

the projection of the gradient on the direction (`P01.07`). It is largest ($= |\nabla f|$) along the gradient, smallest ($= -|\nabla f|$) against it, and zero perpendicular to it. "Per unit distance" in the $(V, \mathrm{AR})$ plane mixes m/s with a dimensionless number, which is a warning sign we return to in §2.6: $|\nabla D|$ depends on the units you chose for the axes, and so does "steepest". At (35, 6): along $+V$ alone the slope is $-48.67$; along $+\mathrm{AR}$ alone $-115.4$; along the steepest-descent direction $-|\nabla D| = -125.2$.

### 2.4 Critical points and the one-variable optimum

A minimum in the interior has $\nabla D = \mathbf 0$ (`@openstax-calculus-volume-3` §4.7). Set (2) to zero at fixed AR:

$$\rho VSC_{D0} = \frac{4KW^2}{\rho V^3S} \iff qSC_{D0} = \frac{KW^2}{qS} \iff V_{md} = \sqrt{\frac{2W}{\rho S}}\left(\frac{K}{C_{D0}}\right)^{1/4}.$$   (5)

At the minimum-drag speed the parasite and induced drags are *equal*, and $D_{\min}(\mathrm{AR}) = 2W\sqrt{KC_{D0}}$, so the best lift-to-drag ratio is $(L/D)_{\max} = 1/(2\sqrt{KC_{D0}})$. Now minimise over AR: $D_{\min} \propto W(\mathrm{AR})\sqrt{K(\mathrm{AR})} \propto (W_0 + c_w\mathrm{AR}^2)\mathrm{AR}^{-1/2}$. Its derivative vanishes when $2c_w\mathrm{AR}^{1/2} = \tfrac12(W_0 + c_w\mathrm{AR}^2)\mathrm{AR}^{-3/2}$, i.e. $4c_w\mathrm{AR}^2 = W_0 + c_w\mathrm{AR}^2$:

$$\mathrm{AR}^* = \sqrt{\frac{W_0}{3c_w}}, \qquad V^* = V_{md}(\mathrm{AR}^*).$$   (6)

For our numbers $\mathrm{AR}^* = \sqrt{13\,000/120} = 10.41$, $W^* = 17\,333$ N (the wing penalty is exactly $W_0/3$ at the optimum, a quarter of the total), $K^* = 0.03823$, $V^* = \sqrt{2 \times 17\,333/(1.225 \times 16)}(0.03823/0.025)^{1/4} = 42.05 \times 1.1122 = 46.77$ m/s, $D^* = 2 \times 17\,333\sqrt{0.03823 \times 0.025} = 1071.7$ N, $(L/D)^* = 16.2$. A 12.9 m span on 16 m² at 1.8 tonnes: a motor-glider.

### 2.5 The Hessian and the second-derivative test

Whether a critical point is a minimum, a maximum or a saddle is decided by the **Hessian**, the symmetric matrix of second partials (§4.7):

$$H = \begin{pmatrix} D_{VV} & D_{V,\mathrm{AR}} \\ D_{V,\mathrm{AR}} & D_{\mathrm{AR,AR}} \end{pmatrix}.$$

The second-order Taylor term is $\tfrac12\Delta\mathbf x^{\mathsf T}H\Delta\mathbf x$ (`P01.08`, equation (8)): the point is a minimum if that quadratic form is positive for every $\Delta\mathbf x$, i.e. if both eigenvalues of $H$ are positive; a maximum if both are negative; a saddle if they differ in sign; undecided if one is zero. (The textbook's $D = f_{xx}f_{yy} - f_{xy}^2 > 0$ test is the same statement: the determinant is the product of the eigenvalues.) At our optimum the entries have closed forms. Write $D = AV^2 + B/V^2$ with $A = \tfrac12\rho SC_{D0}$; $D_{VV} = 2A + 6B/V^4$, and at $V_{md}$, $B = AV^4$, so $D_{VV} = 8A = 4\rho SC_{D0} = 1.96$ N/(m/s)². Along AR, $D_{\min}(\mathrm{AR}) = 2\sqrt{C_{D0}/\pi e}\,(W_0\mathrm{AR}^{-1/2} + c_w\mathrm{AR}^{3/2})$, whose second derivative at $\mathrm{AR}^*$ (use $W_0 = 3c_w\mathrm{AR}^{*2}$) is $3c_w/\sqrt{\mathrm{AR}^*}$ times the prefactor, giving $D_{\mathrm{AR,AR}} = 6c_w\sqrt{K^*C_{D0}} = 7.42$ N. The cross term $D_{V,\mathrm{AR}} = -4(KW^2)'/(\rho V^3S)$ vanishes at the optimum because $(KW^2)' = 0$ is exactly the condition (6). So $H^* = \mathrm{diag}(1.96, 7.42)$, both positive: a minimum, and a very shallow one in $V$: moving 5 m/s off $V^*$ costs $\tfrac12 \times 1.96 \times 25 = 25$ N, 2 % of the drag.

### 2.6 Walking downhill: steepest descent and Newton

If you cannot solve $\nabla f = \mathbf 0$ by hand, walk: from $\mathbf x_k$, step against the gradient,

$$\mathbf x_{k+1} = \mathbf x_k - s\,\nabla f(\mathbf x_k),$$   (7)

with a step $s$ small enough that $f$ decreases (our code halves $s$ until it does). For a quadratic bowl with Hessian eigenvalues $\lambda_{\min} \le \lambda_{\max}$, (7) converges if $s < 2/\lambda_{\max}$ and the error shrinks each step by a factor no better than $(\lambda_{\max} - \lambda_{\min})/(\lambda_{\max} + \lambda_{\min})$: a long narrow valley (large $\lambda_{\max}/\lambda_{\min}$) means zig-zagging and thousands of steps. In raw units $\lambda_{\max} = 7.42$ N, so $s$ must be below 0.27 (m/s)²/N and the best fixed step is near 0.1; choosing $s = 10^{-3}$ because "small is safe" costs 6943 iterations where 0.1 costs 74. **Scaling** the variables, $z = (V/50, \mathrm{AR}/10)$ and $f = D/1000$, makes the step dimensionless and the gradient components comparable; it does not by itself improve the conditioning (here the scaled eigenvalue ratio is 6.6 against 3.8 raw), but it lets you choose $s$ without thinking about units, and in a real design code with ten variables of wildly different magnitudes it is the difference between converging and not.

**Newton's method** uses the Hessian to jump straight to the minimum of the local quadratic model:

$$\mathbf x_{k+1} = \mathbf x_k - H(\mathbf x_k)^{-1}\nabla f(\mathbf x_k),$$   (8)

the multivariable version of `P01.02`'s Newton step. It converges quadratically near a minimum (six iterations here against 97 for scaled descent) and in one step on an exact quadratic, but each step needs the Hessian, which `hessian_fd` builds from $2n^2 + 1$ function evaluations, and it will happily converge to a saddle or a maximum if started badly, since it only seeks $\nabla f = \mathbf 0$. `scipy.optimize.minimize` with BFGS builds an approximate Hessian from the gradients it has already seen: 14 iterations, no second derivatives.

### 2.7 Lagrange multipliers: minimum drag subject to lift equals weight

In §2.1 we eliminated the constraint $L = W$ by substituting $C_L = 2W/(\rho V^2S)$. The method of **Lagrange multipliers** (`@openstax-calculus-volume-3` §4.8) handles the constraint without eliminating it, in the $(V, C_L)$ plane where

$$D = qS(C_{D0} + KC_L^2), \qquad L = qSC_L = W.$$

At the constrained minimum the level curve of $D$ is tangent to the constraint curve, so their gradients are parallel: $\nabla D = \lambda\nabla L$. Components:

$$\frac{\partial D}{\partial V} = \rho VS(C_{D0} + KC_L^2) = \lambda\,\rho VSC_L, \qquad \frac{\partial D}{\partial C_L} = 2qSKC_L = \lambda\,qS.$$   (9)

The second gives $\lambda = 2KC_L$; substituting in the first, $C_{D0} + KC_L^2 = 2KC_L^2$, so

$$C_L^* = \sqrt{C_{D0}/K}, \qquad V^* = \sqrt{\frac{2W}{\rho SC_L^*}}, \qquad \lambda^* = 2\sqrt{KC_{D0}}.$$   (10)

$V^*$ is exactly $V_{md}$ of (5), as it must be. The multiplier is the useful extra: $\lambda = \mathrm{d}D_{\min}/\mathrm{d}W$, the drag cost of one more newton of weight, and $2\sqrt{KC_{D0}} = 1/(L/D)_{\max} = 0.0618$ for our aircraft. Each newton of weight costs 0.062 N of drag at the optimum; a 100 N heavier battery costs 6.2 N of drag, 0.6 %. SymPy solves the three equations (9) plus the constraint symbolically and returns the same three expressions, which is a test.

### 2.8 Where the model stops being valid

- **The weight surrogate.** $W = W_0 + c_w\mathrm{AR}^2$ encodes "longer wings are heavier" and nothing else. Change $c_w$ by 20 % and $\mathrm{AR}^*$ moves 10 %. The optimum is only as good as the structural model behind $c_w$, which is Phase 06's job.
- **Constant $C_{D0}$.** A higher-AR wing at fixed area has a smaller chord, a lower Reynolds number and slightly higher skin-friction coefficient; a sweep of $C_{D0}$ with AR would shift the optimum a little toward lower AR.
- **Objective.** Minimum drag is the maximum-glide-range condition and the maximum-range condition for a piston aircraft; a jet's range is best at minimum $D/V$ and a piston aircraft's endurance at minimum power $DV$ (Phase 10). The *algorithm* is identical; only the function changes.
- **Discrete reality.** Aspect ratio is not continuously adjustable: spars come in sizes, spans must fit hangars (Problem 6), and certification fixes much of $W_0$.

## 3. Build the computational model

Open `starter/main.py`. The aircraft lives in a frozen `Aircraft` dataclass; `drag` broadcasts over grids for the contour plot.

### Step 1: the function and its gradient

```python
def drag(V, AR, ac=DEFAULT_AIRCRAFT):
    q = 0.5 * ac.rho * V**2
    return q * ac.S * ac.cd0 + induced_factor(AR, ac) * weight(AR, ac) ** 2 / (q * ac.S)     # (1)

def drag_gradient(V, AR, ac=DEFAULT_AIRCRAFT):
    K, W = induced_factor(AR, ac), weight(AR, ac)
    dK, dW = -K / AR, 2.0 * ac.c_w * AR
    q = 0.5 * ac.rho * V**2
    dD_dV = ac.rho * V * ac.S * ac.cd0 - 4.0 * K * W**2 / (ac.rho * V**3 * ac.S)               # (2)
    dD_dAR = (dK * W**2 + 2.0 * K * W * dW) / (q * ac.S)                                        # (3)
    return np.array([dD_dV, dD_dAR])
```

### Step 2: gradient and Hessian of anything, by finite differences

```python
def gradient_fd(f, x, h):
    g = np.empty_like(x)
    for i in range(x.size):
        e = np.zeros_like(x); e[i] = h[i]
        g[i] = (f(x + e) - f(x - e)) / (2.0 * h[i])
    return g

def hessian_fd(f, x, h):
    ...  # diagonal: (f(x+e_i) - 2 f(x) + f(x-e_i)) / h_i**2
    ...  # off-diagonal: (f(++) - f(+-) - f(-+) + f(--)) / (4 h_i h_j), then symmetrise
```

`h` may be one step per coordinate, because 1 m/s and 1 unit of AR are not the same size; the tests use `[1e-3, 1e-4]` for the drag.

### Step 3: directions and classification

```python
def directional_derivative(grad, direction):
    d = np.asarray(direction, float)
    return float(np.dot(grad, d / np.linalg.norm(d)))                                          # (4)

def classify_critical_point(H):
    w = np.linalg.eigvalsh(0.5 * (H + H.T))
    ...  # all > 0: "minimum"; all < 0: "maximum"; any ~0: "degenerate"; else "saddle"
```

### Step 4: the two searches

```python
def gradient_descent(f, grad, x0, step=0.1, tol=1e-8, max_iter=10000):
    x = x0.copy(); path = [x.copy()]; fx = f(x)
    for _ in range(max_iter):
        g = grad(x)
        if np.linalg.norm(g) < tol: break
        s = step
        while True:                                   # backtracking: halve s until f decreases
            x_new = x - s * g
            f_new = f(x_new)                          # (7); a ValueError (V <= 0) counts as "no good"
            if np.isfinite(f_new) and f_new < fx: break
            s *= 0.5
        x, fx = x_new, f_new; path.append(x.copy())
    return x, np.array(path)

def newton_minimize(f, grad, hess, x0, tol=1e-10, max_iter=50):
    ...  # x <- x - solve(hess(x), grad(x))                                                    # (8)
```

`scaled_drag` wraps `drag` and `drag_gradient` in the dimensionless variables of §2.6 (the chain rule multiplies each gradient component by its reference value). `minimum_drag_speed`, `optimum_analytic`, `lagrange_min_drag` are (5), (6), (10); `lagrange_sympy` hands the three equations to `sp.solve` and keeps the positive root.

### Verify

Run `python -m pytest curriculum/phase-01/lesson-09-gradients-and-optimization` from the repo root. The tests, in order of how much they would hurt if wrong:

1. **The function.** Parasite part $\propto V^2$, induced part $\propto V^{-2}$ and $\propto 1/\mathrm{AR}$; zero weight kills the induced part; broadcasting; $V \le 0$ and $\mathrm{AR} \le 0$ raise.
2. **Gradient three ways.** Analytic against central differences at three points, rates above 1.9; `gradient_fd` and `hessian_fd` on $x^2y + \sin z$ against the hand Hessian, symmetric.
3. **Critical points.** $\partial D/\partial V = 0$ at $V_{md}$ for three aspect ratios with parasite = induced and $D = 2W\sqrt{KC_{D0}}$; negative slope below, positive above; $V_{md} \propto \sqrt W$; the analytic optimum has zero gradient, a positive-definite Hessian classified "minimum", lower drag than fifty random neighbours, $\mathrm{AR}^* = \sqrt{W_0/3c_w}$, and halves when $c_w$ quadruples; known Hessians classify correctly.
4. **Directions.** $\nabla f \cdot \hat{\mathbf u}$ is $|\nabla f|$ along the gradient, $-|\nabla f|$ against, 0 along the contour, no larger in fifty random directions, and matches a finite difference along four lines through the point.
5. **Searches.** Scaled descent reaches the optimum to $10^{-5}$ with monotone decrease; descent and Newton solve a $2 \times 2$ quadratic bowl exactly (Newton in one step); Newton on the scaled drag converges in at most 12 iterations with a final error ratio below $10^{-2}$.
6. **Lagrange.** $C_L^*$, $V^*$ and the stationarity conditions (9); $\lambda = \mathrm{d}D/\mathrm{d}W$ by central difference to $10^{-6}$ and $= D^*/W$; the SymPy solution equals the closed forms symbolically.
7. **Interface.** `main()` writes the figure.

### Compare

`main()` runs the same minimisation five ways:

```
Analytic optimum: V* = 46.767 m/s, AR* = 10.4083, D* = 1071.696 N (W = 17333 N, L/D = 16.17)
  raw-variable descent, step 0.001: 6943 iterations -> V = 46.767, AR = 10.4083, |grad| = 1.9e-05
  raw-variable descent, step 0.1: 74 iterations -> V = 46.767, AR = 10.4083, |grad| = 1.0e-06
  scaled descent: 97 iterations -> V = 46.767, AR = 10.4083
  Newton: 6 iterations -> V = 46.767, AR = 10.4083
  Hessian at optimum: [[1.96, 0.0], [0.0, 7.4194]], eigenvalues [1.96 7.41943] -> minimum
  scipy.optimize.minimize (BFGS): V = 46.767, AR = 10.4083, 14 iterations
```

All five agree to the printed digits. `scipy.optimize.minimize` is the established tool and what you will call in Phase 03 and beyond; BFGS needs no Hessian and takes 14 iterations because it learns the curvature from the gradients it has seen. The point of writing your own descent and Newton steps once is to know what those 14 iterations are doing and why a badly scaled problem can turn 14 into 14 000. The Hessian entries 1.96 and 7.419 match the closed forms $4\rho SC_{D0}$ and $6c_w\sqrt{K^*C_{D0}}$ of §2.5 to four figures, and the off-diagonal is zero to $10^{-6}$, as predicted.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.**

**Experiment A: step size and scaling.** From $(35, 6)$ run `gradient_descent` on the raw drag with fixed steps $s = 10^{-4}, 10^{-3}, 10^{-2}, 0.03, 0.1, 0.25$ and record the iteration count to $|\nabla D| < 10^{-6}$. Prediction, written before running: the count falls roughly as $1/s$ until $s$ approaches $2/\lambda_{\max} = 2/7.42 = 0.27$, where the backtracking starts halving steps and the count rises again. Repeat in the scaled variables with $s = 0.05, 0.1, 0.2, 0.3, 0.4$ (limit $2/4.90 = 0.41$). Tabulate, and compute $\lambda_{\max}/\lambda_{\min}$ in both scalings from `hessian_fd` at the optimum. Then answer: which scaling has the better conditioning, and is "scaled" the same as "well-conditioned"?

**Experiment B: the price of span.** A hangar limits the span to $b_{\max} = 10$ m, so $\mathrm{AR} \le b_{\max}^2/S = 6.25$. Compute $D_{\min}$ at $\mathrm{AR} = 6.25$ and at $\mathrm{AR}^*$, and the slope $\mathrm{d}D_{\min}/\mathrm{d}\mathrm{AR}$ at the limit (by finite difference of `optimum`-style evaluation at fixed AR, or from the closed form of §2.5). Convert to newtons of drag per metre of span using $\mathrm{d}\mathrm{AR}/\mathrm{d}b = 2b/S$. Prediction: the constrained optimum sits *on* the constraint, the drag penalty is of order 100 N, and each extra metre of span would save of order 60 N. Problem 6 turns this into a second Lagrange multiplier.

Raw data format: two tables, $(s, \text{iterations}, \text{scaling})$ and $(\mathrm{AR}, V_{md}, D_{\min}, \mathrm{d}D/\mathrm{d}\mathrm{AR})$, kept in your lesson notes.

## 5. Validate and improve

Hand (§2), code and the library side by side:

| quantity | hand | code | independent check | difference |
|---|---|---|---|---|
| $D(35, 6)$ (N) | 1451.9 | 1451.94 | | 0 |
| $\nabla D(35, 6)$ | $(-48.67, -115.4)$ | $(-48.668, -115.394)$ | FD $(-48.668, -115.394)$ | $10^{-9}$ relative |
| $\mathrm{AR}^*$ | 10.41 | 10.4083 | descent 10.4083, Newton 10.4083, BFGS 10.4083 | $< 10^{-5}$ |
| $V^*$ (m/s) | 46.77 | 46.767 | all searches 46.767 | $< 10^{-5}$ |
| $D^*$ (N) | 1071.7 | 1071.696 | | 0 |
| $H^*$ | diag(1.96, 7.42) | diag(1.96, 7.419), off-diag $10^{-6}$ | closed forms §2.5 | $10^{-5}$ |
| $\lambda^*$ | 0.0618 | 0.061829 | FD $\mathrm{d}D/\mathrm{d}W$ 0.061829; $1/(L/D) = 1/16.17$ | $< 10^{-6}$ |

Every residual has a cause. The $10^{-9}$ on the gradient is the $h^2$ truncation of the central difference with $h = 10^{-3}$ m/s and $10^{-4}$ in AR. The $10^{-5}$ on the searches is the stopping tolerance on $|\nabla f|$, not an error in the optimum: a gradient of $10^{-9}$ in scaled units corresponds to a position error of about $10^{-9}/\lambda_{\min}$, and the test asks for $10^{-5}$ relative. The off-diagonal Hessian entry is $10^{-6}$ rather than 0 because the finite-difference optimum is itself off by $10^{-6}$ and the cross term grows linearly away from it.

What is still unverified when everything agrees: the model. Five algorithms finding the same minimum of the same function proves the arithmetic, not that an aspect ratio of 10.4 is the right wing for this aircraft. That depends on $c_w$, which we invented, on $e$ and $C_{D0}$, which Phase 10 computes, and on what the aircraft is for: minimum drag is the glide-range optimum and not the fuel-burn optimum. A sensitivity sweep of $\mathrm{AR}^*$ against $c_w$ (it goes as $c_w^{-1/2}$) is the honest companion to the optimum, and it costs one line.

Improvements worth making: a line search along $-\nabla f$ (extension), which removes the step-size guess entirely; a trust region around the Newton step so it cannot leap to a saddle. Improvements worth *not* making: tightening the tolerance below $10^{-8}$ on a function whose inputs are known to two figures.

## 6. Defend your work

Submit: the completed `starter/main.py` with all tests passing; `drag_landscape.png`; the hand calculation of $\nabla D$ at (35, 6), of the optimum and of $H^*$ with the comparison table; the two experiment tables with one paragraph each; and a one-page note (sections 1, 2, 5, 6 and 7 of `projects/_template/report-template.md`) stating what the drag model and the weight surrogate assume and what $\lambda$ means for the design.

### Adversarial self-checks

- **What if the weight doubles?** At fixed AR, $V_{md} \propto \sqrt W$ rises by 41 % and $D_{\min} \propto W$ doubles; $(L/D)_{\max}$ is unchanged. In the surrogate, doubling $W_0$ raises $\mathrm{AR}^*$ by 41 %. Check both with `optimum_analytic(Aircraft(W0=26000))`.
- **What if $c_w \to 0$? What if $C_{D0} \to 0$?** With no weight penalty $\mathrm{AR}^* \to \infty$: the surrogate is what makes the problem well posed. With no parasite drag, $V_{md} \to \infty$ and $D_{\min} \to 0$: induced drag alone has no interior minimum in $V$. Both are limits the code must handle without `nan`, and `minimum_drag_speed` raises rather than dividing by zero.
- **Which assumption fails first in the real world?** The weight surrogate, always: a real wing's weight scales roughly with span cubed divided by chord, with load factor and with material, and $c_w$ hides all of that.
- **What measurement could falsify the model?** A glide test gives $(L/D)_{\max}$ and $V_{md}$ directly; if the measured $V_{md}$ differs from (5) by more than 5 %, either $K$ or $C_{D0}$ is wrong, and the two are separable because $D_{\min} = 2W\sqrt{KC_{D0}}$ and $V_{md}^4 \propto K/C_{D0}$.
- **Does the result survive a change in solver?** Steepest descent, Newton and BFGS must agree to $10^{-5}$; if they do not, one of them has stopped at a saddle (check the Hessian) or has hit the domain boundary ($V \le 0$, caught by the `ValueError`).
- **Which decision changes if the uncertainty doubled?** Suppose $c_w$ is uncertain by ±40 % instead of ±20 %. $\mathrm{AR}^*$ is then $10.4^{+2.4}_{-1.5}$ instead of $10.4^{+1.2}_{-0.9}$, but the drag surface is so flat near the optimum ($D_{\mathrm{AR,AR}} = 7.4$ N) that an AR anywhere between 8 and 13 costs under 30 N (3 %). The decision that changes is not the wing but the structures programme: with ±40 % on $c_w$ the first thing to fund is a better weight model, because that is where the 3 % lives.

## Check yourself

1. Write $\nabla D$ for $D = AV^2 + B/V^2$ with $B = B(\mathrm{AR})$ and say in words what each component measures.
2. At the minimum-drag speed, what is the ratio of parasite to induced drag?
3. The directional derivative of $f$ along $\hat{\mathbf u}$ is $\nabla f \cdot \hat{\mathbf u}$. In which direction is it zero, and what does that say about the gradient and the level curves?
4. A critical point has Hessian eigenvalues $3$ and $-1$. What kind of point is it, and what does the determinant test give?
5. Why does steepest descent with a fixed step diverge when the step exceeds $2/\lambda_{\max}$? (Think of a one-variable parabola $f = \tfrac12\lambda x^2$.)
6. State the Lagrange condition for minimising $D(V, C_L)$ subject to $L(V, C_L) = W$ and the meaning of $\lambda$.
7. At the lesson's optimum, $\lambda = 0.0618$. How much drag does a 50 N heavier instrument panel cost, and what is $(L/D)_{\max}$?
8. Newton's method converged in 6 iterations and steepest descent in 97. Name one reason to prefer steepest descent anyway.
