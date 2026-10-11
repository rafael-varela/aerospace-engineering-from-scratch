# P00.09  Symbolic Math with SymPy

> Let the computer do the algebra; keep the physics, and the checking, for yourself.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the equation of motion of a point mass falling under gravity with quadratic drag, state where it is valid, and find its terminal velocity and exact solution. 2. Encode the derivation in SymPy: symbols with assumptions, `Eq`, `solve`, `diff`, `integrate`, `limit`, `series`, `simplify`. 3. Turn a symbolic expression into a fast NumPy function with `lambdify` and check it against hand values, limiting cases and a SciPy integration. 4. Check a symbolic derivative against a finite difference and explain the error curve (truncation versus rounding). |
| **Prerequisites** | `P00.03`, `P00.06`, `P00.07`; Bridge B1 (algebra, trig); skill check: can differentiate $\tanh$ and integrate $1/(a^2 - v^2)$ with a table; knows what an initial condition is. |
| **Engineering question** | How fast does a 1.2 kg instrument package fall after release, when does drag stop it accelerating, and how far has it fallen by then? |
| **Assumptions and limits** | Point mass; one-dimensional vertical motion with the object never moving upward ($v \ge 0$); constant gravity $g$; constant air density and drag coefficient, so $k = \tfrac12\rho C_D A$ is constant; drag exactly quadratic in speed. Fails for a 10 km drop (density changes), near Mach 1 ($C_D$ changes), at very low speed where drag becomes linear in $v$ (Stokes regime), and for a body that tumbles ($A$ and $C_D$ change). |
| **Analytical work** | Separation of variables to the $\tanh$ solution; terminal velocity; distance fallen; small-time series; no-drag limit. |
| **Computational deliverable** | `starter/main.py`: `acceleration_expr`, `terminal_velocity_expr`, `velocity_expr`, `position_expr`, `ode_residual`, `no_drag_limit`, `small_time_series`, `acceleration_derivative_expr`, four `lambdify_*` functions, `central_difference`, `derivative_check`, `velocity_derivative_check`. Tests prove: the symbolic solution satisfies the ODE identically, $x' = v$, the $k \to 0$ limit is $gt$, the series coefficients are right, lambdified functions match hand values and limits, and the finite-difference check converges at second order. |
| **Practical activity** | **computational experiment**: measure the error of the finite-difference derivative against the symbolic one over eight decades of step size and locate the truncation/rounding crossover. |
| **Equipment** | None |
| **Safety** | None beyond normal computer use. The lesson does not ask you to drop anything; P00.15 does, with its own safety box. |
| **Evidence** | Completed `main.py`, passing tests, `fall_with_drag.png`, the error-versus-$h$ table and log-log plot, and a short note explaining the shape of that curve. |
| **Assessment** | `assessments/phase-00/grading-rubric.md`; competencies 1 (reasoning), 2 (analytical), 3 (computational); pass when the derivation is complete with assumptions named, all public tests pass, and the error-curve note identifies both regimes with numbers. |
| **References** | `@sympy-tutorial` "Introductory Tutorial": "Gotchas" (Symbols, Equals signs), "Basic Operations", "Simplification", "Calculus", "Solvers" (required-free); `@numpy-userguide` "NumPy quickstart" for the arrays `lambdify` returns (supplementary-free). |
| **Extensions** | ★ Solve the ODE directly with `sp.dsolve` and show the result is the same `tanh`. ★ Add the upward-motion branch ($v < 0$) with `sp.Piecewise` or `sp.sign` and verify against `solve_ivp`. ★ Let $\rho$ follow the ISA from P00.03 and show the symbolic approach stops helping. |
| **Time** | unmeasured |

## 1. The engineering problem

An instrument package of mass 1.2 kg, frontal area 0.01 m² and drag coefficient about 1.0 is released from a balloon or a slow aircraft. The parachute team needs three numbers before the flight: the speed it will reach if the parachute fails to open (terminal velocity), how long it takes to get there, and how far it falls meanwhile. Those numbers set the opening altitude and the size of the recovery area.

The physics is Newton's second law with two forces. The hard part is not the physics; it is the algebra of solving a non-linear differential equation and then *trusting* the result. Done by hand, that algebra is where sign errors and lost factors of two live. Done by a computer algebra system, it is fast and exact, but a computer algebra system will happily simplify a wrong equation. So the lesson has two halves: use **SymPy** (the Python library for symbolic mathematics) to derive and solve, and then use the tools of P00.07 to check the symbolic result against limiting cases, a finite-difference derivative and a numerical integration.

Success criterion: an exact expression for $v(t)$ and $x(t)$ that satisfies the equation of motion *identically* (not just numerically), a NumPy function made from it that agrees with `scipy.integrate.solve_ivp` to better than $10^{-7}$ m/s, and a plot of the three predictions.

## 2. Understand and derive

### 2.1 Forces and the equation of motion

Take **downward as positive** and let $v$ be the speed (so $v \ge 0$ while falling). Two forces act: weight $mg$ downward, and aerodynamic drag $D = \tfrac12\rho C_D A v^2$ opposing the motion, i.e. upward.

```
            ▲  D = k v²   (drag, opposes the motion)
            │
         ┌──┴──┐
         │  m  │   package: mass m, frontal area A, drag coefficient C_D
         └──┬──┘
            │
            ▼  m g       (weight)
            ▼  v > 0     (downward positive; the body is falling)
```

Define the **drag constant**

$$k = \tfrac12 \rho C_D A \qquad [\mathrm{kg/m}],$$

so $D = kv^2$. Units: $\mathrm{kg\,m^{-3}\cdot m^2 = kg\,m^{-1}}$, and $kv^2$ has units $\mathrm{kg\,m^{-1}\,m^2\,s^{-2} = N}$. Newton's second law:

$$m\frac{\mathrm{d}v}{\mathrm{d}t} = mg - kv^2 \quad\Longrightarrow\quad \frac{\mathrm{d}v}{\mathrm{d}t} = g - \frac{k}{m}v^2.$$

Assumptions used so far: point mass (no rotation), constant $g$ (altitude change small compared with Earth's radius), constant $\rho$, $C_D$ and $A$ (no tumbling, low Mach, altitude change small compared with the atmosphere's scale height of about 8 km), and $v \ge 0$ (if the body could move upward the drag would have to be written $-k\,v|v|$). For the example: $k = 0.5 \times 1.225 \times 1.0 \times 0.01 = 0.006125$ kg/m.

### 2.2 Terminal velocity

Acceleration vanishes when drag balances weight:

$$g - \frac{k}{m}v_t^2 = 0 \;\Rightarrow\; v_t = \sqrt{\frac{mg}{k}}.$$

Units: $\sqrt{\mathrm{kg \cdot m\,s^{-2} / (kg\,m^{-1})}} = \mathrm{m/s}$. For the example, $v_t = \sqrt{1.2 \times 9.80665/0.006125} = 43.8$ m/s (158 km/h). Note the scaling: doubling the mass raises $v_t$ by $\sqrt2$; doubling the area lowers it by $\sqrt2$. A 1 kg payload under a 2 m² parachute ($C_D \approx 1.5$) has $v_t \approx 2.3$ m/s. This is the first **limiting case** a test should check: $a(v_t) = 0$ exactly.

### 2.3 Exact solution by separation of variables

Write the equation as $\dfrac{\mathrm{d}v}{\mathrm{d}t} = \dfrac{g}{v_t^2}\,(v_t^2 - v^2)$, using $k/m = g/v_t^2$. Separate:

$$\int_0^{v} \frac{\mathrm{d}v'}{v_t^2 - v'^2} = \frac{g}{v_t^2}\,t.$$

The left side is a standard integral: $\dfrac{1}{v_t}\operatorname{artanh}(v/v_t)$ for $0 \le v < v_t$. Hence

$$v(t) = v_t \tanh\!\left(\frac{g\,t}{v_t}\right), \qquad v(0) = 0.$$

Check it: $\mathrm{d}v/\mathrm{d}t = g\,(1 - \tanh^2) = g - (g/v_t^2)\,v^2 = g - (k/m)v^2$. ✓ The **time constant** is $\tau = v_t/g$ (4.5 s for the example): after $t = 2.65\tau$ the speed is 99 % of terminal.

Integrate once more for the distance fallen with $x(0) = 0$:

$$x(t) = \int_0^t v_t \tanh\frac{g t'}{v_t}\,\mathrm{d}t' = \frac{v_t^2}{g}\,\ln\cosh\frac{g\,t}{v_t}.$$

Check: $\mathrm{d}x/\mathrm{d}t = (v_t^2/g)\cdot(g/v_t)\tanh(\cdot) = v_t\tanh(\cdot) = v$. ✓

### 2.4 Two more limiting cases and a series

- **No drag:** as $k \to 0$, $v_t \to \infty$ and $\tanh(gt/v_t) \approx gt/v_t$, so $v \to gt$: free fall. A test asserts `limit(v, k, 0) == g*t`.
- **Long time:** $\tanh \to 1$, so $v \to v_t$ and $x \approx v_t t - (v_t^2/g)\ln 2$: the body falls as if at terminal speed, minus a fixed distance lost while accelerating.
- **Small time:** $\tanh u = u - u^3/3 + 2u^5/15 - \dots$ with $u = gt/v_t$ gives

$$v(t) = gt - \frac{g^2 k}{3m}\,t^3 + \frac{2g^3k^2}{15m^2}\,t^5 - \dots$$

The $t^3$ term is drag's first bite. Its relative size is $gkt^2/(3m)$, so drag changes the speed by 1 % after $t = \sqrt{0.03\,m/(gk)} = 0.77$ s for the example. Series like this are how you decide when a simpler model (free fall) is good enough.

### 2.5 Checking a derivative numerically

If $f$ is smooth, the **central difference**

$$f'(x) \approx \frac{f(x+h) - f(x-h)}{2h}$$

has truncation error $\tfrac{h^2}{6}f'''(\xi)$, and rounding error about $\varepsilon |f(x)|/h$ with $\varepsilon \approx 2.2 \times 10^{-16}$. The total is smallest where the two balance, $h_{\text{opt}} \approx (3\varepsilon |f|/|f'''|)^{1/3}$, which is around $10^{-5}$ to $10^{-4}$ for quantities of order one. For $a(v) = g - kv^2/m$, a quadratic, $f''' = 0$: the central difference is *exact* up to rounding, so the check in `derivative_check` measures pure rounding (about $10^{-13}$). For $v(t)$, which has non-zero $v'''$, the error in `velocity_derivative_check` falls as $h^2$ until rounding takes over. Seeing both behaviours is the point of part 4.

### 2.6 How SymPy thinks

Three facts from the `@sympy-tutorial` "Gotchas" page that save hours. First, a SymPy symbol is not a Python variable: `t = sp.symbols("t")` makes a mathematical object named $t$; assigning to the Python name later does not change expressions already built. Second, `==` tests *structural* equality, not mathematical equality: `(x+1)**2 == x**2 + 2*x + 1` is `False`; use `sp.simplify(a - b) == 0` or `sp.expand`. Third, assumptions matter: declaring `positive=True` lets SymPy simplify $\sqrt{gm/k}^2$ to $gm/k$ and pick the positive root of $v_t^2 = gm/k$ without a case split. Finally, `lambdify((t, m, g, k), expr, modules="numpy")` compiles an expression to a plain Python function that accepts NumPy arrays: that is the bridge from derivation to computation.

Every step of the hand derivation in 2.1–2.5 has a one-line SymPy counterpart, and the code in part 3 is nothing more than this table written out:

| Derivation step | SymPy call | What to check |
|---|---|---|
| declare the quantities | `sp.symbols("t v m g k", positive=True)` | assumptions match the physics ($m, k, g > 0$) |
| write Newton's law | `g - k*v**2/m` | units by hand (SymPy does not track them) |
| find $v_t$ from $a = 0$ | `sp.solve(sp.Eq(a, 0), v)` | exactly one positive root |
| verify a candidate $v(t)$ | `sp.simplify(sp.diff(vc, t) - a.subs(v, vc))` | returns `0`, and a wrong candidate does not |
| distance fallen | `sp.integrate(v_expr, t)` | fix the constant so $x(0) = 0$; differentiate back |
| no-drag limit | `sp.limit(v_expr, k, 0)` | equals `g*t` |
| small-time behaviour | `sp.series(v_expr, t, 0, 4).removeO()` | leading term `g*t`, no $t^2$ term |
| sensitivity $\partial a/\partial v$ | `sp.diff(a, v)` | matches a central difference |
| make it fast | `sp.lambdify((t, m, g, k), v_expr, "numpy")` | hand values, limits, `solve_ivp` |

## 3. Build the computational model

Open `starter/main.py`. The symbols are declared once at the top so every function shares them.

### Step 1: the equation of motion and terminal velocity

```python
import sympy as sp
t, v, m, g, k = sp.symbols("t v m g k", positive=True)

def acceleration_expr():
    return g - k * v**2 / m

def terminal_velocity_expr():
    roots = sp.solve(sp.Eq(acceleration_expr(), 0), v)
    (root,) = [r for r in roots if r.is_positive]      # exactly one, thanks to positive=True
    return root
```

### Step 2: the exact solution, and the identity check

```python
def velocity_expr():
    vt = terminal_velocity_expr()
    return vt * sp.tanh(g * t / vt)

def ode_residual(v_candidate):
    return sp.simplify(sp.diff(v_candidate, t) - acceleration_expr().subs(v, v_candidate))
```

`ode_residual(velocity_expr())` must return exactly `0`. Try `ode_residual(g*t)` too: it returns $-g^2kt^2/m$, which is drag's missing term, and the test `test_ode_residual_is_nonzero_for_a_wrong_candidate` makes sure your residual function is not accidentally always zero. `position_expr` is $v_t^2/g \cdot \log\cosh(gt/v_t)$; the test differentiates it and compares with `velocity_expr`.

### Step 3: limits and series

`no_drag_limit` is `sp.limit(velocity_expr(), k, 0)`. `small_time_series(order)` is `sp.series(velocity_expr(), t, 0, order).removeO()`; the test reads off the $t$, $t^2$ and $t^3$ coefficients with `sp.Poly`. `acceleration_derivative_expr` is `sp.diff(acceleration_expr(), v)`.

### Step 4: lambdify

```python
def lambdify_velocity():
    return sp.lambdify((t, m, g, k), velocity_expr(), modules="numpy")
```

Same pattern for acceleration, position and $\mathrm{d}a/\mathrm{d}v$. Argument order is a contract: the tests call `vel(ts, M, G, K)`.

### Step 5: finite differences

`central_difference(f, x, h)` is three lines and must reject $h \le 0$. `derivative_check` compares the lambdified $\mathrm{d}a/\mathrm{d}v$ with the central difference of the lambdified $a$ at speed $v_0$. `velocity_derivative_check` compares $a(v(t_0))$ with the central difference of $v(t)$ at $t_0$: this is the ODE verified *numerically*, which is what you would do if no closed form existed.

### Verify

- **Analytic case:** `test_velocity_solves_the_ode_with_zero_initial_speed` (residual exactly 0), `test_position_is_the_integral_of_velocity`.
- **Limiting cases:** $v(0) = 0$; $a(v_t) = 0$; $k \to 0$ gives $gt$; $t \to \infty$ gives $v_t$; scaling of $v_t$ with $m$ and $k$.
- **Conservation-like property:** the trapezoid integral of the lambdified $v$ reproduces the lambdified $x$ to $10^{-4}$ m over 10 s; $v$ is monotone and bounded by $v_t$.
- **Convergence:** `test_velocity_derivative_check_converges_at_second_order` requires a $\log_2$ rate above 1.8 as $h$ halves.

### Compare

`main()` integrates the same ODE with `scipy.integrate.solve_ivp` (`rtol=1e-10, atol=1e-12`) and reports the maximum difference from the exact solution over 15 s: about $9 \times 10^{-9}$ m/s. The difference is the integrator's error, not the formula's: tighten the tolerances and it shrinks; loosen them to the defaults (`rtol=1e-3`) and it grows to $10^{-2}$ m/s. The exact solution is the reference here; in Phase 03 the roles reverse and the integrator is tested on problems like this one precisely because the answer is known.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Map the error of the finite-difference derivative against the symbolic one.

1. Fix the example parameters ($m = 1.2$ kg, $g = 9.80665$ m/s², $k = 0.006125$ kg/m) and $t_0 = 3$ s.
2. For $h \in \{10^{-1}, 10^{-2}, \dots, 10^{-8}\}$ record `velocity_derivative_check(3.0, m, g, k, h)` and `derivative_check(20.0, m, g, k, h)`.
3. Before running, write down your prediction: at which $h$ is the first check smallest, and what slope do you expect on a log-log plot on each side of that minimum?
4. Plot both errors against $h$ on log-log axes, add the lines $C h^2$ and $\varepsilon |v|/h$, and save as `fd_error_vs_h.png`.

Raw data format: CSV with columns `h, err_velocity_ode, err_dadv`.

## 5. Validate and improve

Put the three error estimates side by side: the analytic truncation term $\tfrac{h^2}{6}|v'''(t_0)|$, the rounding estimate $\varepsilon|v(t_0)|/h$, and the measurement. You should see the velocity check fall by a factor of 100 per decade from $h = 10^{-1}$ ($3.1 \times 10^{-5}$ m/s²) to $10^{-3}$ ($3.2 \times 10^{-9}$), reach a minimum near $h = 10^{-4}$ ($5 \times 10^{-11}$), and then *rise* by a factor of 10 per decade: $2.8 \times 10^{-9}$ at $10^{-6}$ and $4 \times 10^{-7}$ at $10^{-8}$. The slope +2 on the left is truncation; the slope −1 on the right is rounding, because the subtraction $v(t+h) - v(t-h)$ loses leading digits as $h$ shrinks. The optimum sits near $(3\varepsilon|v|/|v'''|)^{1/3} \approx 1.3 \times 10^{-4}$ with $|v| \approx 25$ m/s and $|v'''| \approx 0.02$ m/s⁴ for this case, not at the textbook $\varepsilon^{1/3} \approx 6 \times 10^{-6}$, which assumes $|v| \sim |v'''| \sim 1$. Explain that discrepancy in your note: the rule of thumb has hidden scale factors.

The acceleration check behaves differently: flat at $10^{-15}$ to $10^{-13}$ for every $h$ above $10^{-4}$, because a quadratic has no $h^2$ truncation term; only rounding is visible, and it too grows as $h$ shrinks. If you had seen a smooth $h^2$ decay there you would have to suspect your `acceleration_expr` is not quadratic.

Improvement: for the velocity check, use the five-point stencil $(-f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h))/(12h)$ and show the left slope becomes +4 while the right slope stays −1, moving the optimum to a larger $h$ and a smaller minimum error.

## 6. Defend your work

Submit `main.py`, the test log, `fall_with_drag.png`, `fd_error_vs_h.png` with its CSV, and a one-page note (template in `projects/_template/report-template.md`) containing the derivation with assumptions, the three predictions for the package ($v_t$, time to 99 % of $v_t$, distance fallen by then) and the error-curve explanation.

### Adversarial self-checks

- What if the mass doubles? $v_t$ rises by $\sqrt2$ and the time constant $v_t/g$ with it; the distance fallen in the first second barely changes. Check all three with the lambdified functions.
- What if $v \to 0$? Drag vanishes quadratically, so the early motion is pure free fall; is there a speed below which quadratic drag is the *wrong* model (hint: Reynolds number below about 1)?
- Which assumption fails first in the real world? For a 1 km drop, none seriously; for a 30 km stratospheric drop, constant density: $\rho$ at 30 km is 1.5 % of sea level, so the early "terminal velocity" is eight times larger and the body decelerates as it descends. The symbolic solution no longer exists in closed form.
- What measurement could falsify the model? A video of the drop (P00.15) giving $v(t)$: a fit of $v_t\tanh(gt/v_t)$ with $v_t$ free should recover $k$ within the uncertainty of $\rho C_D A$.
- Does the result survive a change of tool? `solve_ivp` agrees to $10^{-8}$ m/s; `sp.dsolve` returns the same $\tanh$ (extension). Which decision changes if $C_D$ is uncertain by ±20 %? $v_t$ by ∓10 %: the parachute opening altitude must carry that margin.

## Check yourself

1. Write the equation of motion for a falling body with quadratic drag and name the two assumptions that make $k$ a constant.
2. Derive the units of $k$ and of $v_t = \sqrt{mg/k}$.
3. Why is `sp.simplify(a - b) == 0` the right test of equality and `a == b` the wrong one?
4. State the no-drag limit of $v(t)$ and the long-time limit, and say which SymPy function finds each.
5. The central difference of a quadratic is exact. Why, then, does `derivative_check` not return exactly zero?
6. For the example package, how long until the speed is within 1 % of terminal, and how far has it fallen by then? (Use $\operatorname{artanh}(0.99) = 2.647$.)
7. What does `positive=True` buy you when solving $v_t^2 = gm/k$?
8. On a log-log plot of finite-difference error against $h$, what slopes do you expect on each side of the minimum, and which physical constant sets the right-hand side?
