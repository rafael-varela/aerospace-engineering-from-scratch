# P00.14  Measurement, Calibration & Uncertainty 🔬

> A number without an uncertainty is not a measurement. It is a guess with decimals.

| | |
|---|---|
| **Learning outcomes** | 1. Evaluate a Type A standard uncertainty from repeated readings and a Type B standard uncertainty from a resolution or a stated tolerance. 2. Derive and apply the law of propagation of uncertainty to a sum, a product, a quotient and a power law, by hand and with a finite-difference Jacobian. 3. Calibrate an instrument against a known reference, fit slope and offset, and read the residuals. 4. Report a result as value ± expanded uncertainty with coverage factor k = 2 and the right number of significant figures. 5. Decide which input dominates an uncertainty budget and what to improve first. |
| **Prerequisites** | `P00.03`, `P00.04`, `P00.07`, `P00.13`; skill check: can compute a mean and a standard deviation by hand for five numbers; can differentiate x², √x and 1/x. |
| **Engineering question** | What does your kitchen or pocket scale actually read for a known mass, how far can you trust it, and what is the mass of an unknown object with a defensible ± afterwards? |
| **Assumptions and limits** | Inputs are independent unless a correlation is stated; the measurand depends smoothly on its inputs so a first-order (linear) propagation holds; the instrument responds linearly over the calibrated range; k = 2 gives about 95 % coverage only when the combined distribution is close to normal and based on enough readings. The method fails for a hysteretic or drifting instrument, for strongly non-linear functions with large input uncertainties, and for a budget dominated by a Type A term from very few readings. |
| **Analytical work** | Derivation of the propagation law from a first-order Taylor expansion; the sum, product, quotient and power rules; the uncertainty of a stack of n coins; dynamic pressure q = ½ρv² and drag coefficient C_D = D/(qS) as aerospace cases. |
| **Computational deliverable** | `starter/main.py`: `type_a`, `type_b_resolution`, `type_b_half_width`, `combine`, `jacobian_fd`, `propagate`, `expanded`, `format_measurement`, `calibration_fit`, `apply_calibration`, `corrected_value`. Tests check the propagation against analytic sum, product, quotient and power-law cases, the correlated-difference limit, second-order convergence of the Jacobian, exact recovery of a known calibration line, and the reporting format. |
| **Practical activity** | **physical experiment**: calibrate a 0.1 g pocket scale or a 1 g kitchen scale with coins of known nominal mass or with measured volumes of water, then measure an unknown with a full uncertainty budget; optionally check a ruler against a standard card and a phone accelerometer against g. See `projects/p00-calibration/`. |
| **Equipment** | Level 1: a scale, coins or water with a syringe or measuring cup, a ruler; optional: a phone with a sensor app. See `projects/p00-calibration/bom.yaml`. |
| **Safety** | See `projects/p00-calibration/safety.md`. Low energies; the only hazards are spilled water near electronics and small parts near children. |
| **Evidence** | `data/calibration.csv`, the residual plot with error bars, a complete uncertainty budget table, the unknown's mass as value ± U (k = 2), and the problem set. |
| **Assessment** | `assessments/phase-00/grading-rubric.md`; competencies 1 (derivation), 2 (budget by hand), 3 (toolkit and tests), 4 (which term dominates and why). Threshold: the reported unknown must carry a k = 2 uncertainty whose budget accounts for the scale's resolution, its repeatability, the reference's tolerance, and the calibration fit. |
| **References** | `@jcgm-100-2008` §4.2 "Type A evaluation of standard uncertainty", §4.3 "Type B evaluation of standard uncertainty", §5.1 "Uncorrelated input quantities", §6 "Determining expanded uncertainty", §7 "Reporting uncertainty" (required-free); `@nist-tn-1297` §2 "Classification of components of uncertainty", §4 "Type B evaluation of standard uncertainty", §5 "Combined standard uncertainty", §6 "Expanded uncertainty" (required-free). |
| **Extensions** | ★ Weighted least squares with per-point uncertainties and a χ² goodness-of-fit. ★ Monte Carlo propagation (sample the inputs, run the function, read the output spread) and compare with the linear law for t = √(2h/g) at large u_h. |
| **Time** | unmeasured |

## 1. The engineering problem

Every number this course will ever ask you to defend comes with the same question: how far can it be wrong? In `P00.15` you will predict the fall time of a ball to the millisecond and then measure it. Whether the two agree is a question that cannot be answered until both carry an uncertainty. In Phase 10 you will put a model on a home-made wind-tunnel balance and read a drag force of 0.12 N. A drag coefficient computed from that force is worthless unless you know whether the balance reads 0.12 ± 0.01 N or 0.12 ± 0.05 N. In Phase 14 a rocket's mass depletion m(t) = m₀ − ṁt feeds every performance number; the uncertainty in ṁ, measured by weighing before and after a static test, sets how much you can claim.

The instrument you will calibrate here is a kitchen or pocket scale, because it is cheap, you probably own one, and it has every feature of a real instrument: a finite resolution (the last digit), a repeatability (it does not read the same twice), an offset (it may not read zero at zero), a gain error (it may read 0.99 g per real gram), and a range outside which it is meaningless. The procedure you learn transfers unchanged to a load cell, a pressure sensor, a thermocouple, or the accelerometer in your phone.

The success criteria: you produce a calibration line (slope and offset with their uncertainties), a residual plot showing no structure larger than the uncertainty bars, and the mass of an unknown object as a value ± U (k = 2) whose budget lists every contribution and names the largest. When you can do that, "the scale says 17.3 g" becomes "the object's mass is 17.45 ± 0.16 g (k = 2), limited by the coins' tolerance". The first is a reading; the second is a measurement.

## 2. Understand and derive

Vocabulary first, following the GUM (`@jcgm-100-2008`). The **measurand** is the quantity you want (the mass of the object). An **input quantity** is anything the measurand is computed from (the raw reading, the calibration slope and offset). The **standard uncertainty** u(x) of an input is one standard deviation of where its true value could be. Combining the standard uncertainties of the inputs gives the **combined standard uncertainty** u_c(y) of the measurand. Multiplying by a **coverage factor** k gives the **expanded uncertainty** U = k·u_c, the half-width of the interval you report. With k = 2 that interval contains the true value with roughly 95 % probability when the distribution is close to normal.

### 2.1 Type A: uncertainty from repeated readings

Put the same 10 g reference on the scale five times, lifting it off and re-centring it each time: 9.9, 10.0, 9.9, 10.0, 9.8 g. The mean is x̄ = 9.92 g. The sample standard deviation,

s = √( Σ(xᵢ − x̄)² / (n − 1) ) = 0.084 g,

describes the scatter of a *single* reading. But you report the mean, and the mean of n readings scatters less than one reading does:

u_A = s / √n = 0.084 / √5 = 0.037 g.

This is the **Type A** evaluation (`@jcgm-100-2008` §4.2): a statistical analysis of a series of observations. Two things to notice. The n − 1 is there because x̄ was computed from the same data and used up one degree of freedom. And the 1/√n means that quadrupling the readings halves the uncertainty: a slow way to improve, which is why the next section matters.

### 2.2 Type B: uncertainty from everything else

The scale displays 0.1 g steps. When it shows 10.0, the true input was anywhere between 9.95 and 10.05 with, as far as you know, equal probability. That is a **rectangular distribution** of half-width a = 0.05 g, whose standard deviation is a/√3:

u_B(resolution) = (0.1 / 2) / √3 = 0.1 / √12 = 0.029 g.

This is a **Type B** evaluation (`@jcgm-100-2008` §4.3): an uncertainty from anything other than statistics on your own readings: a specification, a tolerance, a resolution, a calibration certificate, or physics. Three distributions cover almost every case you will meet:

| You know | Distribution | u |
|---|---|---|
| the value is somewhere within ± a, nothing more | rectangular | a / √3 |
| within ± a, and the middle is more likely | triangular | a / √6 |
| a quoted "± a at 95 %" (k = 2) | normal | a / 2 |

A coin is the reference you are most likely to use. A mint publishes a nominal mass and a tolerance for each coin; look up yours. For example the United States Mint specifies a nickel at 5.000 g with a tolerance of ±0.194 g. Treat that as rectangular: u(coin) = 0.194 / √3 = 0.112 g. That is three times the scale's resolution uncertainty, and it is the honest number: you do not know where in the tolerance band your particular coin sits. (You can do better by using a coin fresh from a bank roll, whose wear is zero, or by using water: see §2.5.)

Resolution, repeatability and reference tolerance are independent sources, so they add in quadrature:

u_c = √( u_A² + u_B,res² + u_B,ref² ).

The quadrature sum is not a convention; it is what the variance of a sum of independent random variables does. Which brings us to the general rule.

### 2.3 The law of propagation of uncertainty

Suppose the measurand is y = f(x₁, x₂, …, x_N) and each input has a standard uncertainty u(xᵢ). Expand f to first order about the estimates:

Δy ≈ Σᵢ (∂f/∂xᵢ) Δxᵢ.

The partial derivative cᵢ = ∂f/∂xᵢ is the **sensitivity coefficient**: how much y moves per unit move of xᵢ. Square both sides and take expectations. If the Δxᵢ are independent, every cross term averages to zero, and

u_c²(y) = Σᵢ cᵢ² u²(xᵢ).      (GUM eq. 10; `@jcgm-100-2008` §5.1)

Each product |cᵢ| u(xᵢ) is that input's **contribution** to the budget, in units of y. Comparing them tells you where to spend effort. Three special cases you should be able to write without thinking:

- **Sum or difference**, y = a ± b: c = ±1, so u(y) = √(u_a² + u_b²). Note: a difference has the *same* uncertainty as the sum, not a smaller one.
- **Product or quotient**, y = ab or a/b: divide by y, and (u(y)/y)² = (u_a/a)² + (u_b/b)². Relative uncertainties add in quadrature.
- **Power**, y = aᵖ: u(y)/y = |p| u_a/a. For y = √a the relative uncertainty halves; for y = 1/a it is unchanged; for y = a² it doubles.

**Worked example: fall time.** t = √(2h/g), with h = 1.200 m, u_h = 0.005 m, g = 9.81 m/s², u_g = 0.03 m/s². The power rule twice: t ∝ h^½ g^−½, so

u_t / t = ½ √( (u_h/h)² + (u_g/g)² ) = ½ √( (0.0042)² + (0.0031)² ) = 0.0026.

With t = 0.4946 s, u_t = 0.0013 s. The height term (0.0042) dominates the g term (0.0031), but not by much; measuring h to 2 mm would make g the limit. That is exactly the kind of decision the budget exists to inform.

**Worked example: dynamic pressure and drag coefficient.** In a wind tunnel, q = ½ ρ v². Relative uncertainty: (u_q/q)² = (u_ρ/ρ)² + (2 u_v/v)². The factor 2 on the velocity term is the power rule, and it is why airspeed is the quantity worth measuring well. Then C_D = D/(q S) is a quotient of three things, so (u_C/C)² = (u_D/D)² + (u_q/q)² + (u_S/S)². A 3 % velocity uncertainty alone puts 6 % on C_D before the balance has been touched.

**When the law fails.** It is first order. If f curves noticeably across the range ±u of an input (think 1/x near x = 0, or √h with u_h comparable to h), the linear estimate is wrong, and the GUM's Supplement 1 (Monte Carlo propagation) is the remedy. For the measurements in this course, with relative uncertainties of a few percent, the linear law is accurate to a fraction of a percent of itself.

### 2.4 Correlated inputs

If two inputs move together, the cross terms do not vanish. The general form is u_c²(y) = Σᵢ Σⱼ cᵢ cⱼ u(xᵢ) u(xⱼ) r(xᵢ, xⱼ), with r the correlation coefficient. The important home case: two readings taken on the same scale share its offset error. Their *difference* (the tare-and-weigh procedure) cancels that shared error completely: with r = 1 and equal u, the difference has u = 0 from that source, while the sum has u = 2u, not √2·u. The toolkit's `propagate` takes an optional correlation matrix so you can see this happen.

### 2.5 Calibration: the straight-line fit

Place known references x₁ … x_n (0, 5, 10, … g) and record the readings y₁ … y_n. If the scale is linear, y = b x + a, where b is the **slope** (gain: 1.000 for a perfect scale) and a the **offset** (zero error). Ordinary least squares gives

b = Σ(xᵢ − x̄)(yᵢ − ȳ) / Σ(xᵢ − x̄)²,   a = ȳ − b x̄.

The **residuals** rᵢ = yᵢ − (b xᵢ + a) are what the line failed to explain. Their standard deviation s = √(Σrᵢ²/(n − 2)) (two parameters were fitted, so n − 2) is the scale's scatter about the line, and from it the uncertainties of the fitted parameters follow: u(b) = s/√S_xx and u(a) = s √(1/n + x̄²/S_xx), where S_xx = Σ(xᵢ − x̄)².

To use the calibration, invert it: an unknown reading y_u corresponds to a mass x_u = (y_u − a)/b. That is a function of three inputs (y_u, a, b), so §2.3 applies and the toolkit's `corrected_value` does it.

**Reading the residual plot.** Plot rᵢ against xᵢ with error bars of ± u_c(reading). Residuals that scatter randomly inside the bars mean the line is adequate. A bow (negative at the ends, positive in the middle, or the reverse) means the scale is non-linear and a straight line is the wrong model. A step means something changed mid-session. A trend in the *order the points were taken* means drift (warm-up, battery). The plot is the single most informative picture in a calibration; always make it.

**Water as a reference.** If no coin specification is available, water works: at 20 °C its density is 998 kg/m³, so 10.0 mL has a mass of 9.98 g. The limit is now your volume measurement. A 10 mL kitchen syringe read to ±0.25 mL (rectangular) gives u = 0.14 g, worse than a nickel's tolerance; a laboratory-style 10 mL measuring cylinder is somewhat better. State which reference you used and its uncertainty: that is half of what a calibration report *is*.

### 2.6 Reporting: coverage factor and significant figures

Report U = k u_c with k = 2 and say so: "m = 17.45 g, U = 0.16 g (k = 2)" or "17.45 ± 0.16 g (k = 2)". Two rules for digits (`@jcgm-100-2008` §7; `@nist-tn-1297` §7): give the uncertainty to **two significant figures**, and round the value to the **same decimal place** as the uncertainty. 17.4512 ± 0.1587 g becomes 17.45 ± 0.16 g. Writing 17.4512 ± 0.16 claims digits the uncertainty says you do not have; writing 17.5 ± 0.2 throws away information you do have. Keep all digits in intermediate calculations; round once, at the end.

A warning about k = 2 and "95 %". The 95 % holds when the combined uncertainty is close to normally distributed and is not dominated by a Type A term from very few readings. With n = 3 readings the right factor for 95 % is about 4.3, not 2 (Student's t with 2 degrees of freedom). The practical rule for this course: take at least 5 repeated readings, prefer 10, and if a Type A term from n < 6 dominates your budget, say so in the report.

## 3. Build the computational model

The toolkit is a small set of functions you will reuse in every later lesson that touches data. Build it in this order; each step has a test.

### Step 1: Type A and Type B

```python
def type_a(readings):
    x = np.asarray(readings, dtype=float).ravel()
    n = x.size
    if n < 2:
        raise ValueError("Type A evaluation needs at least two readings")
    return float(x.mean()), float(x.std(ddof=1) / math.sqrt(n))

def type_b_half_width(a, distribution="rectangular"):
    divisors = {"rectangular": math.sqrt(3), "triangular": math.sqrt(6), "normal": 2.0}
    return a / divisors[distribution]

def type_b_resolution(resolution):
    return type_b_half_width(resolution / 2, "rectangular")
```

`ddof=1` gives the n − 1 denominator. `type_b_resolution(0.1)` must equal `type_b_half_width(0.05)`: a test checks it.

### Step 2: The Jacobian by central differences

You could differentiate every function by hand. Instead, compute each ∂f/∂xᵢ numerically: nudge input i up and down by a small step h and take the slope.

```python
def jacobian_fd(f, x, rel_step=1e-6):
    x0 = np.asarray(x, dtype=float).ravel()
    jac = np.empty_like(x0)
    for i in range(x0.size):
        h = rel_step * max(abs(x0[i]), 1.0)
        xp, xm = x0.copy(), x0.copy()
        xp[i] += h; xm[i] -= h
        jac[i] = (f(xp) - f(xm)) / (2 * h)
    return jac
```

The central difference has an error proportional to h², so shrinking the step tenfold shrinks the error a hundredfold until round-off takes over near h ≈ 10⁻⁸ relative. A test checks that convergence. The step scales with |xᵢ| so that inputs of very different size (a height in metres, a g in m/s²) are treated alike.

### Step 3: The propagation law

```python
def propagate(f, x, u, corr=None):
    x0, u0 = np.asarray(x, float).ravel(), np.asarray(u, float).ravel()
    c = jacobian_fd(f, x0)
    cu = c * u0                          # each input's contribution, in units of y
    var = cu @ cu if corr is None else cu @ np.asarray(corr) @ cu
    return Propagated(value=float(f(x0)), u=math.sqrt(max(var, 0.0)),
                      sensitivities=c, contributions=np.abs(cu))
```

`cu @ cu` is Σ(cᵢuᵢ)²: GUM equation 10 in one line. With a correlation matrix the same line becomes the full double sum. Return the contributions too: the budget table in your report is just those numbers.

### Step 4: Reporting

`round_to_uncertainty(value, u, sig=2)` finds the decimal place of the second significant figure of u (`decimals = sig - 1 - floor(log10(u))`), rounds both numbers there, and `format_measurement` prints them with the unit. Watch the case where rounding u carries it into the next decade (0.0996 → 0.10): recompute the place once.

### Step 5: Calibration fit

Implement the formulas of §2.5 directly with NumPy (`np.sum`, no `polyfit`, so that you own every number), return a `CalibrationFit` with slope, offset, their uncertainties, the residuals and s, and write `apply_calibration` and `corrected_value` on top. `corrected_value` is nothing more than `propagate(lambda p: (p[0] - p[2]) / p[1], [reading, slope, offset], [u_reading, u_slope, u_offset])`: the toolkit calling itself.

### Verify

`tests/test_main.py` checks, against closed forms: `type_a([1, 2, 3])` gives (2, 1/√3); u_A falls as 1/√n; the resolution rule; the three half-width divisors; the Jacobian against analytic partials and its second-order convergence; propagation through a sum (quadrature), a product and a quotient (relative quadrature), a constant multiple, a power law, zero inputs, and a fully correlated difference (cancels); exact recovery of a known line; residuals summing to zero on noisy data; `apply_calibration` inverting the line; and the reporting format on four cases. Run `pytest curriculum/phase-00/lesson-14-measurement-calibration-uncertainty`.

### Compare

NumPy has `np.polyfit(x, y, 1, cov=True)`. Run it on your calibration data and compare slope, offset and the square roots of the covariance diagonal with your own `calibration_fit`. They should agree to round-off. If `polyfit`'s uncertainties differ by a constant factor, read its documentation on how it scales the covariance (it rescales by the residual variance unless told otherwise) and explain the factor. Then compare `propagate` with the hand-derived power rule for t = √(2h/g) from §2.3; the two should agree to five figures. Any difference beyond that is a bug in one of them.

## 4. Experiment at home, or investigate real data

**Mode: physical experiment.** The full bill of materials, safety box and step-by-step procedure are in `projects/p00-calibration/` (`bom.yaml`, `safety.md`, `procedure.md`). Read the safety box first, every session.

> **Safety (summary; the authoritative box is `projects/p00-calibration/safety.md`).** Hazards: water near a battery-powered scale and a phone; small coins near children or pets. Energy limits: masses under 0.5 kg, no height above the bench, battery voltage ≤ 5 V, no mains. Controls: scale on a dry tray, water measured away from electronics, coins counted back at the end. Supervision: none needed for an adult. Stop criteria: liquid on any electronics; a reading that keeps changing with nothing on the pan (drift) until the cause is found. Disposal: none. Not permitted: heating water, stacking anything above 0.5 kg, mains-powered scales near water.

The procedure in outline:

1. **Predict first.** Before any reading, write down what you expect the slope and offset to be (a new scale: 1.00 and 0.0), and the uncertainty you expect for a single reading from its resolution alone (0.029 g for a 0.1 g scale). Record it in `data/prediction.md`.
2. **Warm up and zero.** Switch on, wait two minutes, tare with the empty pan. Read the empty pan ten times, lifting nothing: that is your zero repeatability.
3. **Repeatability.** Place one reference coin ten times, removing and re-centring it each time. Compute x̄ and u_A.
4. **Calibration.** Stack 0, 1, 2, … 6 coins (or 0, 5, 10, … mL of water). Read each twice, in random order, and record: `reference_g, u_reference_g, reading_g, order` in `data/calibration.csv`. Remember that a stack of n coins each with tolerance u_coin has u = √n · u_coin if the coins' errors are independent, not n · u_coin.
5. **The unknown.** Weigh an object of unknown mass (a key, a small bolt) five times.
6. **Optional, a second instrument.** Check a ruler against a standard-size card (ISO/IEC 7810 ID-1 cards are 85.60 mm wide; treat ±0.2 mm as a conservative half-width) at several positions along the ruler. Or lay your phone flat and record its accelerometer for 30 s with a free sensor app: the mean should be close to 9.81 m/s² and the scatter is its Type A uncertainty; the difference from your local g is its offset.

Raw data are never edited. Analysis happens in a copy, with the toolkit.

## 5. Validate and improve

Now put the three sources of knowledge side by side.

**Analysis against measurement.** Your prediction said slope 1.00 and offset 0. Your fit says, for example, slope 0.9914 ± 0.0026 and offset 0.000 ± 0.046 g. Is 0.9914 consistent with 1.0? The difference is 0.0086, which is 3.3 times u(slope): the scale under-reads by almost 1 %, and that is a real finding, not noise. Is the offset consistent with zero? 0.000 ± 0.046: yes. Write both sentences in the report.

**Computation against hand calculation.** For one point of the budget (say the 10 g reference), compute u_c by hand from §2.1 and §2.2 and compare with `combine(type_a(...), type_b_resolution(...))`. They must agree exactly. Then compare `propagate` for the corrected unknown with the quotient rule applied by hand.

**The budget table.** Every report in this course carries one, in this form (numbers from the illustrative example in `main()`, for a raw reading of 17.3 g):

| Source | Type | Distribution | Half-width or s | u (g) | Sensitivity c | Contribution \|c\|·u (g) |
|---|---|---|---|---|---|---|
| reading repeatability (n = 5) | A | normal | s = 0.084 | 0.037 | 1/b = 1.009 | 0.038 |
| display resolution 0.1 g | B | rectangular | a = 0.05 | 0.029 | 1/b = 1.009 | 0.029 |
| calibration slope b | A (fit) | normal | u(b) = 0.0026 | 0.0026 | −(y−a)/b² = −17.6 | 0.045 |
| calibration offset a | A (fit) | normal | u(a) = 0.046 | 0.046 | −1/b = −1.009 | 0.046 |
| reference tolerance (enters through fit) | B | rectangular | mint value | per reference | through b and a | (state separately if large) |
| **combined u_c** | | | | | | **0.080** |
| **expanded U, k = 2** | | | | | | **0.16** |

Reported: m = 17.45 ± 0.16 g (k = 2). The table is the proof behind the ±: a reader can check every row and see at a glance which one to attack.

**The budget.** Tabulate the contributions returned by `corrected_value`: reading repeatability, scale resolution, slope, offset, reference tolerance (this last one enters through the fit; if you used coins with a large tolerance, propagate it explicitly as an extra input). One of them will be largest. In the example above the reading, slope and offset terms are all near 0.045 g, so no single improvement helps much; with nickels as references the ±0.194 g tolerance, entering as 0.112/√n per stack, would dominate everything and the honest conclusion is "this calibration is limited by the reference, not the scale".

**Discrepancies.** If a residual lies outside its bar, find out why; "experimental error" is not an explanation. Common real causes: the pan was tared with a coin's residue of water on it; the scale auto-zeroed during a slow stack; a coin from circulation was worn 0.1 g light; the reading was taken before the display settled. Each one has a measurable test (re-tare and repeat; weigh that coin alone; wait and re-read). Do the test and report what it showed.

**Improve and re-run.** Pick the dominant term and reduce it: more repeats if Type A dominates (slow: 1/√n), a better reference if tolerance dominates (fast), a scale with finer resolution if resolution dominates (costs money). Re-run the fit and show the new budget next to the old one. The report's redesign section is this paragraph with your numbers in it.

## 6. Defend your work

Submit: `data/prediction.md`, `data/calibration.csv` (raw), the analysis script, the residual plot with error bars, the budget table, the unknown's mass as value ± U (k = 2), the problem set, and a report using `projects/_template/report-template.md`.

### Adversarial self-checks

- What if the unknown's mass doubles? Which terms of the budget scale with it (slope) and which do not (offset, resolution)? Does the *relative* uncertainty improve or worsen?
- What if the reference tolerance were half as large? Recompute. What if it were zero: what then limits you?
- Which assumption fails first in the real world: independence of the coins' errors, linearity of the scale, or the normality behind k = 2?
- What measurement could falsify the straight-line model? (Hint: it is a shape in the residual plot.)
- Does your slope survive a change in the order the points were taken? Re-fit using only the first half of the readings and only the second half.
- Which decision changes if the measurement uncertainty were twice as large: would you still conclude the scale under-reads?

## Check yourself

1. Define Type A and Type B evaluation in one sentence each. Which one is "statistical"?
2. A display reads in steps of 0.5 g. What is the standard uncertainty from resolution alone?
3. Five readings have sample standard deviation 0.20 g. What is the standard uncertainty of their mean? How many readings would halve it?
4. y = a − b with u_a = 0.3 and u_b = 0.4. What is u_y? Does the minus sign matter?
5. q = ½ρv². If v has 2 % relative uncertainty and ρ 1 %, what is the relative uncertainty of q?
6. Report 3.14159 with u = 0.0996 correctly.
7. Your residual plot bows upward in the middle. What does that tell you about the scale, and what should you do?
8. Two masses weighed on the same un-tared scale are subtracted. Does the scale's zero offset enter the difference? Why?
