# P01.08  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and name the variable held fixed in every partial derivative you write.

Several problems use your personal seed **s**: the day of the month you were born (1–31). $\rho = 1.225$ kg/m³ and $g = 9.80665$ m/s² unless stated.

## 1. (Warm-up: a linear wing)

A wing has $S = 10 + 0.2s$ m², lift-curve slope $a = 4.5 + 0.02s$ per radian, zero-lift angle $\alpha_0 = -2°$ and no stall roll-off ($k_3 = 0$). It flies at $V = 50$ m/s and $\alpha = 5°$.

Find: (a) $C_L$ and $L$; (b) $\partial L/\partial\alpha$ in N/rad and N/deg, and $\partial L/\partial V$; (c) the lift increments for $+1°$ alone, for $+1$ m/s alone, and for both, first from the total differential and then exactly; (d) the percentage error of the total differential in each case. Check (a)–(c) with `lift`, `lift_partials` and `lift_linearized` using `Wing(S=..., a=..., k3=0.0)`. Answer form: numbers with units, four significant figures.

## 2. (Derivation: the calculus of the lift function)

For $L = \tfrac12\rho V^2 S\,C_L(\alpha)$ with $C_L = a(\alpha - \alpha_0) - k_3(\alpha - \alpha_0)^3$:

(a) Derive $L_\alpha$, $L_V$, $L_{\alpha\alpha}$, $L_{VV}$, $L_{\alpha V}$ and $L_{V\alpha}$ and confirm the last two are equal.
(b) Show that along the level curve $L = W$, $\mathrm{d}V/\mathrm{d}\alpha = -\tfrac12 V\,C_L'(\alpha)/C_L(\alpha)$ in two ways: by implicit differentiation of $L(\alpha, V(\alpha)) = W$, and by differentiating the explicit $V(\alpha) = \sqrt{2W/(\rho S C_L)}$.
(c) Show that for a pure speed change $\Delta V$ at fixed $\alpha$ the tangent-plane error is *exactly* $\tfrac12\rho S C_L\,\Delta V^2$ (no higher-order terms), and that for a pure $\alpha$ change it is $\tfrac12 L_{\alpha\alpha}\Delta\alpha^2 + \tfrac16 L_{\alpha\alpha\alpha}\Delta\alpha^3$ exactly.
(d) Write the second-order two-variable Taylor remainder (8) and use (a) to show that the cross term $L_{\alpha V}\Delta\alpha\Delta V$ equals $2\,\mathrm{d}L_\alpha\,\Delta V/V$, where $\mathrm{d}L_\alpha = L_\alpha\Delta\alpha$ is the linear $\alpha$-increment. Interpret the factor $2\Delta V/V$.

Answer form: derivations with each held-fixed variable named.

## 3. (Unfamiliar scenario: airspeed from a Pitot tube, with uncertainty)

An airspeed indicator computes $V = \sqrt{2q/\rho}$ from the measured dynamic pressure $q$ and an assumed density $\rho$. On one flight $q = 1500 \pm 15$ Pa and $\rho = 1.225 \pm 0.020$ kg/m³ (independent, one standard deviation each).

(a) Compute $V$ and the partials $\partial V/\partial q$ and $\partial V/\partial\rho$ with units.
(b) Propagate the two uncertainties with the total differential in quadrature: $\delta V^2 = (V_q\delta q)^2 + (V_\rho\delta\rho)^2$. Which input dominates?
(c) Show that the relative sensitivities are $\partial\ln V/\partial\ln q = \tfrac12$ and $\partial\ln V/\partial\ln\rho = -\tfrac12$, and redo (b) in percentages.
(d) On a hot day at a high airfield the density is 15 % below the assumed value and the pilot flies by indicated airspeed (the $V$ the instrument shows for the assumed $\rho$). Is the true airspeed higher or lower than indicated, by what percentage, and does the *lift* at a given indicated airspeed change? (Use $L = \tfrac12\rho V^2 S C_L = qSC_L$.)

Answer form: one speed, two partials, one uncertainty with the dominant term named, two percentages, one short argument.

## 4. (Numerical, uses your code: the step-size study)

Run Experiment A of lesson part 4 for the reference point ($\alpha = 4°$, $V = 55$ m/s, default `Wing`).

(a) Tabulate the errors of `partial_derivative` in $\alpha$ and in $V$ for $h = 10^{-1}, \dots, 10^{-8}$ and plot them on log–log axes.
(b) Fit the slope of the $\alpha$ error for $h \ge 10^{-4}$ and for $h \le 10^{-7}$; compare the large-$h$ branch with the predicted $\tfrac16 h^2|L_{\alpha\alpha\alpha}|$ using $L_{\alpha\alpha\alpha} = -6k_3 qS$.
(c) Explain the $V$ column: why is there no sloped branch on the left, and at what $h$ does rounding become visible?
(d) Repeat (a) for the mixed partial with `mixed_partial` in both orders and report the largest $h$ at which they agree to $10^{-6}$ relative.

Answer form: a table, a plot, two slopes, one comparison, two short explanations.

## 5. (Judgment: where the linear sensitivity stops)

For each situation, say whether the total differential (7) is adequate (error below 2 % of $L_0$), and if not what you would use instead; give the number that decides:

(a) a $2°$ vertical gust at the reference point;
(b) a $6°$ vertical gust with a simultaneous $12$ m/s headwind gust;
(c) a $4°$ gust at $\alpha = 14°$, close to the model's stall;
(d) the same $2°$ gust during a descent at 95 m/s at sea level (hint: Mach number).

Answer form: four verdicts, each with a number and one sentence.

## 6. (Hard: lift through altitude, the chain rule with three variables)

Let density depend on altitude through the ISA troposphere, $\rho(h) = \rho_0(1 - 2.2558 \times 10^{-5}h)^{4.2559}$ with $h$ in metres and $\rho_0 = 1.225$ kg/m³ (`P00.03`), so that $L = L(\alpha, V, \rho(h))$.

(a) Write $\partial L/\partial\rho$ and $\mathrm{d}\rho/\mathrm{d}h$, and evaluate both at $h = 2000$ m.
(b) The aircraft climbs at $\dot h = 5$ m/s with $\alpha = 4°$ and $V = 55$ m/s held constant. Compute $\mathrm{d}L/\mathrm{d}t$ by the chain rule and check it against a central difference of $L(\alpha, V, \rho(h(t)))$ in $t$.
(c) To hold the lift constant during the climb the pilot must raise $\alpha$. Find the required $\dot\alpha$ in degrees per second from $L_\alpha\dot\alpha + L_\rho\dot\rho = 0$, and the altitude at which $\alpha$ would reach the model's stall angle if $V$ stays at 55 m/s.
(d) Add the time derivative of $V$ as a third term. If the climb also bleeds speed at $\dot V = -0.1$ m/s², which of the two terms dominates $\mathrm{d}L/\mathrm{d}t$, by what factor, and what $\dot\alpha$ is now needed to hold the lift? State in one sentence what the chain rule is doing in this problem that a single-variable derivative could not.

Answer form: two expressions with values, one rate with its finite-difference check, two $\dot\alpha$ values and one altitude, one factor, one sentence.
