# P00.14  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units and uncertainty on every numerical answer; report expanded uncertainties with k = 2 and say so.

## 1. (Warm-up: one point of a budget)

A 0.1 g pocket scale reads a 20.00 g reference coin eight times: 19.9, 20.1, 20.0, 19.8, 20.0, 20.1, 19.9, 20.0 g. The mint's tolerance for the coin is ±0.05 g.

Find: the mean, the Type A standard uncertainty of the mean, the Type B uncertainty from resolution, the Type B uncertainty from the tolerance (rectangular), the combined standard uncertainty, and the expanded uncertainty. Report the mean correctly. Then state whether the scale's offset at 20 g is consistent with zero. Answer form: a budget table and one sentence. (If you have your own scale and a reference, use your own eight readings instead.)

## 2. (Derivation: the rules from the law)

Starting from the first-order expansion Δy ≈ Σ (∂f/∂xᵢ) Δxᵢ and the assumption of independent inputs:

(a) Derive u(y) for y = a + b and y = a − b.
(b) Derive the relative-uncertainty rule for y = a b and y = a / b.
(c) Derive the rule for y = aᵖ.
(d) Apply (b) to the ISA sea-level density ρ = p / (R T) with p = 101 325 ± 100 Pa, T = 288.15 ± 0.5 K and R = 287.05 J/(kg·K) taken as exact. Give ρ with k = 2 and name the dominant input.

## 3. (Unfamiliar scenario: a wind-tunnel drag coefficient)

A home-built balance measures the drag on a model: D = 0.120 ± 0.004 N. The airspeed is v = 12.0 ± 0.3 m/s, the air density ρ = 1.20 ± 0.01 kg/m³, and the reference area S = 0.0120 ± 0.0002 m² (all standard uncertainties).

Find C_D = D / (½ ρ v² S) with its expanded uncertainty. Tabulate the relative contribution of each input and say which measurement you would improve first and why. Then: by how much would the uncertainty in C_D fall if the airspeed uncertainty were halved, everything else unchanged?

## 4. (Numerical: uses your code)

Use your own `data/calibration.csv` if you have it. Otherwise generate a stand-in with your seed (for example the day of your birth), clearly labelled as synthetic in your write-up:

```python
rng = np.random.default_rng(SEED)
reference = np.arange(0.0, 61.0, 10.0)                   # g
reading = 0.985 * reference + 0.3 + rng.normal(0.0, 0.08, reference.size)
```

(a) Run `calibration_fit` and report slope and offset with their standard uncertainties. Is the slope consistent with 1.000 and the offset with 0.0? Quantify with |difference| / u.
(b) Plot the residuals with error bars of ±u_c where u_c combines `type_b_resolution(0.1)` and the fit's `s_residual`. Describe what you see.
(c) An unknown object reads 23.4 g; take its reading uncertainty as 0.05 g. Use `corrected_value` to report its mass with k = 2 and list the three contributions.
(d) Compare slope, offset and their uncertainties with `np.polyfit(reference, reading, 1, cov=True)`. Explain any factor between the two sets of uncertainties.

## 5. (Judgment: when the straight line is wrong)

Under what conditions does the linear calibration model y = b x + a stop being valid for a scale? Name at least three distinct physical causes, and for each give one measurable sign in the data (a shape in the residual plot, a trend against time or order, a difference between loading and unloading). Then say what you would do in each case: change the model, restrict the range, or fix the procedure.

## 6. (Correlation, and the difference between an uncertainty and a mistake)

(a) Two readings on the same scale, 31.2 g and 18.7 g, each with u = 0.05 g from resolution and repeatability, share the scale's unknown zero offset, which you estimate as 0.0 ± 0.3 g. Compute the uncertainty of their difference (i) treating the offset as an independent input of each reading, and (ii) recognising it as fully correlated. Use `propagate` with a correlation matrix to confirm by computation.
(b) The Mars Climate Orbiter was lost in 1999 because one piece of ground software produced thruster impulse in pound-force-seconds while the receiving software expected newton-seconds, a factor of 4.45. Explain, in a paragraph, why this was not a measurement uncertainty, why no amount of uncertainty analysis of the thruster data would have found it, and which check from `P00.05` would have.

## 7. (Coverage: when k = 2 is not 95 %)

A budget is dominated by a Type A term from n = 3 readings.

(a) Explain why k = 2 then understates the 95 % interval. Look up (any free table) the Student's t factor for 95 % with 2 degrees of freedom and compare.
(b) How many readings would you need before k = 2 is within about 10 % of the correct 95 % factor?
(c) ★ Optional: for t = √(2h/g) with h = 0.30 ± 0.10 m (a deliberately large uncertainty), compare the linear-law u_t with a Monte Carlo estimate (sample h from a normal distribution 10⁵ times, compute t, take the standard deviation). Report the ratio and explain the sign of the difference.
