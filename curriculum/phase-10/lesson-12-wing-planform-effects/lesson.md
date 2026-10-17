# P10.12  Wing Planform Effects

> Taper and twist shape the span load; the span load sets induced drag, root bending and where the stall begins. Optimise one of them alone and you will pay in the others, so put all three in the problem.

| | |
|---|---|
| **Learning outcomes** | 1. Compute (not quote) the span efficiency of straight-tapered wings and locate the minimum of the induced-drag factor near taper 0.35. 2. Compute the local lift coefficient along the span, the station where a wing of uniform section first stalls, and how taper, washout and sweep move it. 3. Formulate and solve the twist-optimisation problem (minimum induced drag at a design CL) with SLSQP: unconstrained (it recovers the elliptic load, e → 1), with a linear twist (checked against the Lagrange-multiplier closed form), and with constraints on the root bending moment and the tip's stall margin. 4. Reproduce Prandtl's 1933 bell-shaped span load by optimisation, and explain which constraint makes it optimal. 5. Recommend a planform for the reference trainer and defend it. |
| **Prerequisites** | `P10.10` (lifting line, Glauert coefficients, δ), `P10.11` (vortex lattice, for sweep), `P10.09` (induced drag, span efficiency), `P03.14` (constrained optimisation, SLSQP, Lagrange multipliers), `P03.16` (multidisciplinary trade-offs: aerodynamics against structure), `P05.12` (root bending moment of a span load). Skill check: you can minimise a quadratic form subject to one linear equality by Lagrange multipliers. |
| **Engineering question** | Should the reference trainer's rectangular, untwisted wing be tapered or twisted, and what would each choice buy in induced drag, cost in stall behaviour and do to the root bending moment? Answer (lifting line, same span and area): taper 0.36 maximises e (0.987 against 0.934, 5.4 N less induced drag at cruise), but moves the first-stalling station from the root to 66 % of the semi-span; restoring an inboard stall with a 5 % margin at the approach CL needs 4.1° of washout and drops e to 0.766 (121 N, 22 N more than the rectangular wing). The rectangular wing with 3.0° of tip washout instead reaches e = 0.9998 (93 N) at cruise, keeps a root stall and lowers the root bending moment by 7.6 %. Recommendation: keep the rectangular planform and add about 3° of washout. |
| **Assumptions and limits** | Lifting line (`P10.10`) for straight wings, flat sections with a0 = 2π and camber only as an angle shift, linear section lift (the "stall station" is where cl/CL is largest; real stall depends on the section's cl_max along the span); sweep effects by the vortex lattice (`P10.11`). Twist is piecewise linear between 5 stations (root twist 0). Root bending moment from the aerodynamic load only (no inertia relief). Prandtl's problem uses ∫L′y² dy as the structural-weight proxy and requires the lift to be non-negative. One parameter of the trainer is varied at a time, and the text says which. |
| **Analytical work** | The induced-drag quadratic form CDi = xᵀQx in the design variables; the Lagrange solution for a linear twist; the root bending coefficient C_B = M_root/(qSb) and its elliptic value CL/(3π); Prandtl's problem reduced to f(s) = 1/s² + 3(1/s³ − 1/s)² with the non-negative-lift bound s ≤ √(3/2). |
| **Computational deliverable** | `starter/main.py`: `wing_from_coefficients`, `efficiency_vs_taper`, `cl_ratio`, `stall_onset`, `sweep_stall_station`, `hat_twists`, `twist_basis`, `optimize_twist`, `linear_twist_optimum`, `prandtl_bell`, `bending_trade` (a compact `P10.10` lifting line `ll_solve`, `tapered_chord`, `ll_gamma`, a compact `P10.11` lattice `vlm_strip_cl` and `main` are given). Tests prove: an elliptic load has e = 1 and C_B = CL/(3π); e peaks at taper 0.3–0.42 for three AR; the local lift integrates to CL; taper moves the stall outboard, washout inboard, aft sweep outboard; the twist basis superposes; the unconstrained optimum tends to e = 1 as twist stations are added and gives the elliptic root moment; the linear-twist SLSQP optimum equals the Lagrange closed form and is a minimum; an elliptic planform needs no twist; constraints hold and cost drag monotonically; Prandtl's bell (span × √(3/2), drag × 8/9) is recovered. |
| **Practical activity** | **computational experiment**: e against taper for four aspect ratios; stall-onset station against taper, washout and sweep; twist optimisation with three constraint sets; the bending-moment trade; Prandtl's bell. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `planform.png`, the derivations of the Lagrange optimum and of Prandtl's f(s), the experiment tables, and a one-page planform recommendation for the trainer naming the active constraints and their price. |
| **Assessment** | `assessments/phase-10/grading-rubric.md`; competencies 1–4. Pass: the Lagrange optimum is derived and matched by the optimiser; all public tests pass; the e(λ) curve is computed, not quoted; the recommendation states the induced-drag, stall and bending consequences with numbers. |
| **References** | `@ntrs-naca-tr-921-span-loading` (DeYoung & Harper, https://ntrs.nasa.gov/citations/19930091985): Effect of Plan-Form Variation on the Basic Loading; Effect of Plan-Form Variation on the Additional Loading (supplementary-free): span loads of tapered, swept and twisted wings. `@mit-ocw-16-unified-fluids` Fluids (spring) Lecture 8 "Wing Geometry, Wing Design Problem" (https://ocw.mit.edu/courses/16-01-unified-engineering-i-ii-iii-iv-fall-2005-spring-2006/resources/f08_sp/) and Lecture 9 "General Wings" (required-free). `@leishman-intro-aerospace-flight-vehicles` "Wing Shapes & Nomenclature" (https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/wing-shapes-and-nomenclature/): Geometric Definition of a Wing; Calculation of Mean Chords; and "Aerodynamics of Finite Wings": Effects of Compressibility & Wing Sweep (required-free). `@ntrs-naca-tn-1422-sivells-three-wings` (https://ntrs.nasa.gov/citations/19930082118) (required-free): an AR-9, taper-0.4 wing tested with 0° and 2° washout. `@ntrs-naca-tn-1671-swept-wings-taper-ratio` (https://ntrs.nasa.gov/citations/19930082298) (required-free): measured effects of taper on swept wings. `@scipy-userguide` "Optimization (scipy.optimize)" > "Constrained minimization": "Sequential Least SQuares Programming (SLSQP) Algorithm" (required-free): the driver. Optional textbooks for deeper study (optional-paid, not required; library or used copy): Raymer, *Aircraft Design: A Conceptual Approach*, for planform selection in design practice; Anderson, *Fundamentals of Aerodynamics*, for e against taper in lifting-line form. |
| **Extensions** | ★ Add the wing's weight as inertia relief to the bending moment and redo the bending trade (`P05.12`). ★ Make the chord distribution a design variable too (taper and twist together) with a constraint on the tip chord for manufacture. ★ Solve Jones's 1950 variant: fixed lift and fixed *root* bending moment with the span free, and compare its span and drag with Prandtl's. |
| **Time** | unmeasured |

## 1. The engineering problem

`P10.10` and `P10.11` analysed a given wing. A designer chooses the wing, and the main levers after span and area are **taper** (how the chord shrinks towards the tip), **twist** (washout: the tips set at a lower angle than the root) and **sweep**. Each changes the span load. The span load sets three things that pull in different directions:

- **Induced drag**, least for an elliptic load (`P10.09`);
- **Root bending moment**, least when the load is pushed inboard, which costs induced drag;
- **Where the stall begins**: if the local lift coefficient cl(y) peaks near the tips, the tips stall first; the ailerons lose effectiveness and the aircraft rolls off into a spin entry. For a trainer, a root-first stall is a safety requirement, not a preference.

Textbooks often say that a taper ratio near 0.35 gives a nearly elliptic load and that washout prevents tip stall. This lesson computes both statements for the trainer, and then poses the design question properly, as an optimisation with constraints (`P03.14`, `P03.16`): minimise induced drag at the design lift coefficient by choosing the twist, subject to a bending-moment limit and a stall-margin requirement. It ends with Prandtl's 1933 discovery that, once structure is accounted for, the elliptic load is not optimal at all.

Success criteria, written first: e(λ) computed for four aspect ratios with the location of its maximum; the unconstrained twist optimum reaching e > 0.9999 with enough twist stations; the optimiser agreeing with the Lagrange closed form for a linear twist to 10⁻⁶ rad; every constrained optimum feasible and its cost reported; Prandtl's numbers (span × 1.2247, drag × 0.8889) recovered.

## 2. Understand and derive

### 2.1 What taper does

In lifting-line theory (`P10.10`) the load follows the chord and the angle: a rectangular wing is "too full" near the tips compared with an ellipse (δ = 0.0707 at AR 8.33), a pointed wing is too full inboard (δ = 0.16). Somewhere between, a straight taper approximates the ellipse best; the code finds the minimum of δ, i.e. the maximum of e, by scanning λ. The local lift coefficient is a different matter:

cl(y)/CL = 2Γ(y)/(V∞ c(y) CL),

and dividing by a chord that shrinks towards the tip *raises* cl outboard. A tapered wing can carry a nearly elliptic *load* while its *lift coefficient* peaks well outboard. With the same section cl_max everywhere, the wing stalls first where cl/CL is largest.

### 2.2 Twist as a design variable

Twist adds a **basic load** that does not depend on α (`@ntrs-naca-tr-921-span-loading`: basic and additional loading). The lifting line is linear, so with α and the twists t_k at stations η_k (piecewise linear in |η| = |2y/b|, root twist zero) collected in x = [α, t₁, …, t_K],

**A = B x,**   (1)

where column 0 of B holds the Glauert coefficients for α = 1 rad, and column k those for a unit "hat" of twist at η_k. Then

CL = πAR (B x)₁,  **CDi = πAR Σ n (Bx)ₙ² = xᵀ Q x**,  Q = πAR Bᵀ diag(n) B.   (2)

The problem *minimise xᵀQx subject to CL = CL_design* is a quadratic programme with one linear equality. With enough twist stations its solution is the elliptic load, whatever the planform: twist can make any planform elliptic at *one* CL.

### 2.3 The linear-twist optimum by Lagrange (worked example)

With one twist variable (tip twist t, linear in |η|), A = αa + tr, where a and r are the coefficient vectors for α = 1 and t = 1. Fix A₁ = CL/(πAR); eliminate α = (A₁ − t r₁)/a₁; then Aₙ = pₙ + t qₙ with p = A₁a/a₁ and q = r − r₁a/a₁ (q₁ = 0). Minimise Σ n(pₙ + tqₙ)²:

**t* = −Σ n pₙqₙ / Σ n qₙ².**   (3)

For the rectangular trainer at CL = 0.3309: t* = −2.614°, e = 0.99064. The SLSQP run with one twist station must return the same.

### 2.4 Constraints

- **Root bending moment.** For one wing, M_root = ρV∞∫₀^{b/2} Γ y dy; the coefficient C_B = M_root/(qSb). An elliptic load puts each half's lift at 4(b/2)/(3π) from the root: **C_B,ell = CL/(3π)**. Pushing load inboard lowers C_B and raises CDi.
- **Stall margin at the tips.** At a high-lift CL (the approach CL = 1.112 here, with the same twist and α re-trimmed), require cl(η) ≤ (1 − m) CL_high for η ≥ 0.7: the outboard 30 % runs below the wing's average, so the stall begins inboard. Optionally the twist may only decrease outboard ("washout only"), as a builder would make it.

### 2.5 Prandtl's bell (1933)

Prandtl asked a different question: for a given lift and a given *structural weight*, which span and load give the least induced drag? A spar sized to its bending moment weighs roughly in proportion to ∫M dy, which for a cantilever equals ½∫L′(y)y² dy. With Γ = Σ Gₙ sin nθ (odd n) and y = −(b/2) cos θ:

L ∝ bG₁,  I = ∫L′y² dy ∝ b³(G₁ + G₃),  D_i ∝ Σ nGₙ² (independent of b).

In units of the elliptic wing of span b₀ (G₁ = 1, G₃ = 0): s = b/b₀, G₁ = 1/s, G₃ = 1/s³ − 1/s, and higher harmonics only add drag, so

**D_i/D_i,ell = f(s) = 1/s² + 3(1/s³ − 1/s)².**   (4)

f falls as s grows, but G₃/G₁ = 1/s² − 1 becomes more negative, and the load near the tips, proportional to (1 + 3G₃/G₁)θ, turns negative beyond **G₃/G₁ = −1/3, s = √(3/2) = 1.2247**. Download at the tips does not lighten a spar sized to |M|, so non-negative lift is a constraint, and the optimum sits on it: **22.5 % more span, 11.1 % less induced drag (8/9)**, with the load Γ ∝ sin θ − ⅓ sin 3θ = (4/3) sin³θ = (4/3)(1 − η²)^{3/2}, the **bell**. There f′(s) = 0 as well: the bell is also the stationary point of the Lagrangian, which is why it appears in the classical derivation. Its root bending moment is 2 % below the elliptic wing's of span b₀.

## 3. Build the computational model

Open `starter/main.py`. Given: `tapered_chord`, `ll_solve` (compact `P10.10` lifting line returning the odd coefficients), `ll_odd_n`, `ll_gamma`, `vlm_strip_cl` (compact `P10.11` lattice with Hough spacing), `main`.

### Step 1: loads and the stall station

`wing_from_coefficients` returns CL, CDi, e, C_B (by the midpoint rule over the right wing) and cl(y). `efficiency_vs_taper` scans λ; `cl_ratio` trims α to a CL by linearity; `stall_onset` returns the η of the largest cl/CL; `sweep_stall_station` does the same with the lattice.

### Step 2: the optimiser

`hat_twists` and `twist_basis` build B of (1). `optimize_twist` forms Q of (2), the CL equality with its gradient, and the optional constraints, and calls `scipy.optimize.minimize(method="SLSQP")` (`@scipy-userguide`). The stall constraint is linear in x: re-trimming to CL_high only adds a multiple of column 0.

```python
Q = (np.pi * AR) * (Bm.T * n) @ Bm                      # CDi = x' Q x
cons = [dict(type="eq", fun=lambda x: np.pi * AR * (Bm[0] @ x) - CL_design, jac=lambda x: np.pi * AR * Bm[0])]
res = minimize(lambda x: x @ Q @ x, x0, jac=lambda x: 2 * Q @ x, constraints=cons, method="SLSQP")
```

`linear_twist_optimum` is (3); `bending_trade` loops over C_B limits; `prandtl_bell` minimises (4) on [1, √(3/2)] and checks it with SLSQP over (s, G₁, …, G₁₁) with the two equalities and Γ ≥ 0.

Running `main()` prints (abridged):

```
flat untwisted wings, AR 8.333: e at taper 0, 0.2, 0.35, 0.5, 1.0 = 0.8624, 0.9746, 0.9868, 0.9822, 0.9339; maximum e = 0.98689 at taper 0.36
twist optimisation at the cruise CL = 0.3309 (5 twist stations, root twist 0):
  rectangular trainer: untwisted e = 0.9339; optimum e = 0.99978, twist (deg) at eta [0.   0.25 0.5  0.75 1.  ] = [ 0.    -0.112 -0.513 -1.218 -2.998], alpha 4.583 deg, C_B 0.03511 (elliptic 0.03511)
  linear twist, Lagrange closed form: tip twist -2.6142 deg, e = 0.99064
  linear twist, SLSQP: tip twist -2.6142 deg, e = 0.99064
Prandtl 1933 (lift and int L y^2 dy fixed, span free, lift >= 0): span x 1.2247 (sqrt(3/2) = 1.2247), df/ds there -3.3e-10, induced drag x 0.8889 (-11.11 %), G3/G1 = -0.3333 (bell: -1/3), root moment x 0.9798
```

and writes `planform.png`: e against taper, cl/CL along the span at the approach CL for three tapers and the optimised twist, and the elliptic and bell loads of equal lift.

### Verify

`python -m pytest curriculum/phase-10/lesson-12-wing-planform-effects` runs 13 tests, listed in the metadata. The optimiser tests are of two kinds: against mathematics (the Lagrange closed form, stationarity of the optimum under perturbation, the elliptic planform needing no twist, Prandtl's 8/9) and against feasibility (every constraint holds at the returned point and its cost rises monotonically as it tightens).

### Compare

Three independent checks agree. The optimiser and the Lagrange formula (3) give the same tip twist to 10⁻⁶ rad. The unconstrained optimum's root bending coefficient equals the elliptic value CL/(3π) = 0.03511, as it must if the load is elliptic. And Prandtl's problem solved two ways (a one-variable bounded minimisation of (4), and SLSQP over six Fourier coefficients and the span) gives s = 1.2247 and 0.8889 both times, with G₅…G₁₁ at zero.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction before running it.

1. **e against taper** (flat, untwisted). *Finding:* AR 8.33: e = 0.862, 0.975, 0.987, 0.982, 0.934 at λ = 0, 0.2, 0.35, 0.5, 1; maximum 0.98689 at λ = 0.36. AR 6: maximum 0.991 at 0.37 (rectangular 0.954); AR 12: 0.980 at 0.35 (0.907); AR 20: 0.967 at 0.34 (0.861). The dip of δ (peak of e) sits at λ ≈ 0.34–0.37 for every AR, and the penalty of a rectangular wing grows with AR.
2. **Where the stall begins** (AR 8.33, approach CL, largest cl/CL): λ = 1: root (η = 0), peak 1.123; λ = 0.7: η = 0.36; 0.5: η = 0.53; 0.35: η = 0.66 (peak 1.069, and cl/CL = 0.955 at η = 0.9); 0.2: η = 0.79 (and 1.076 at η = 0.9). Linear washout on λ = 0.35: −2° moves it to η = 0.52, −4° to 0.34. Sweep (lattice, λ = 0.35, AR 8.33): −20° (forward) η = 0.42; 0° 0.62; +20° 0.72; +35° 0.77, with the peak cl/CL rising from 1.048 to 1.157: aft sweep and taper both push the stall towards the ailerons.
3. **Twist optimisation for the rectangular trainer** at the cruise CL. *Finding:* 5 stations give e = 0.99978 with twist 0, −0.11, −0.51, −1.22, −3.00° (α = 4.58° at the root), and the elliptic root moment. More stations: e = 0.998484 (2), 0.999783 (4), 0.999970 (8), 0.999996 (16). Linear twist: −2.614°, e = 0.99064, from both the optimiser and (3).
4. **The tapered variant (λ = 0.35, same b and S) with a stall constraint** at CL_high = 1.112, washout only:

| constraint | e at cruise | CDi vs its unconstrained optimum | twist at η = 0, 0.25, 0.5, 0.75, 1 (°) | stall onset η |
|---|---|---|---|---|
| none | 0.99913 | — | 0, +0.42, +0.73, +0.93, −0.78 (wash-in allowed) | 0.664 |
| tips ≤ CL_high | 0.90031 | +10.98 % | 0, 0, −0.09, −2.37, −2.37 | 0.470 |
| tips ≤ 0.95 CL_high | 0.76615 | +30.41 % | 0, 0, −0.65, −4.11, −4.11 | 0.383 |
| tips ≤ 0.90 CL_high | 0.62477 | +59.92 % | 0, 0, −1.22, −5.85, −5.85 | 0.308 |

Predict the price of the 5 % margin before running. *Finding:* at the cruise CL the induced drag rises from 93 N (unconstrained) to 121 N; the tapered wing made safe is worse than the untwisted rectangular wing (99.6 N), which stalls at the root anyway.
5. **Bending trade** (rectangular trainer, cruise CL, fixed span; the untwisted wing's C_B is 0.0380, 8.3 % above elliptic):

| C_B limit / elliptic | e | tip twist |
|---|---|---|
| 1.00 | 0.99978 | −3.00° |
| 0.97 | 0.99270 | −3.77° |
| 0.94 | 0.97191 | −4.55° |
| 0.91 | 0.93908 | −5.33° |
| 0.88 | 0.89666 | −6.11° |

A 9 % cut in root moment costs as much e as the untwisted rectangular wing already pays; a 12 % cut costs 0.10 of e.
6. **Prandtl's bell.** *Finding:* span × 1.2247, induced drag × 0.8889, G₃/G₁ = −0.3333, root moment × 0.98; SLSQP agrees with all higher harmonics at zero.

Raw data format: CSV tables `AR, taper, e` (1); `taper, washout_deg, sweep_deg, eta_stall, max_cl_over_CL` (2); `n_stations, e, twist_deg...` (3); `margin, e, twist_deg..., eta_stall` (4); `CB_ratio, e, tip_twist_deg` (5); plus `planform.png`.

## 5. Validate and improve

| quantity | theory | code | difference | cause |
|---|---|---|---|---|
| linear twist optimum | Lagrange (3): −2.6142° | SLSQP −2.6142° | < 10⁻⁶ rad | — |
| unconstrained optimum, 16 stations | e → 1 | 0.999996 | 4 × 10⁻⁶ | finite twist freedom |
| C_B of the elliptic optimum | CL/(3π) = 0.03511 | 0.03511 | < 2 × 10⁻³ rel. | — |
| Prandtl's span and drag | √(3/2), 8/9 | 1.2247, 0.8889 | < 10⁻⁵ | flat f′ at the bound |
| e(λ) maximum | near 0.35 (classical) | 0.36 (AR 8.33) | — | computed |

The optimiser is verified against closed forms. The model's validity is that of the lifting line (`P10.10`: about 5 % on the lift slope against the lattice, more for e near the tips) and of a stall criterion that assumes the same section cl_max everywhere: a real wing with a thinner tip section or a lower tip Reynolds number stalls at the tips sooner than this predicts, which makes the washout requirement *larger*, not smaller.

The engineering answer: **for the trainer, taper buys at most 5 N of induced drag at cruise but moves the stall outboard, and buying the stall back with washout costs 22 N; the rectangular wing with about 3° of tip washout gives the least induced drag (93 N), a root-first stall and a 7.6 % smaller root bending moment.** That is why simple trainers keep the "Hershey-bar" rectangular wing; a long-span sailplane, whose induced drag matters far more, tapers and washes out and accepts the complexity. To improve: re-run with the lattice of `P10.11` as the analysis engine, with section cl_max varying along the span, and with the bending moment and the twist in one multidisciplinary problem (extension 2, `P03.16`).

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (2), (3) and (4); `planform.png`; the six experiment tables; and a one-page planform recommendation (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) naming the active constraints and their price in newtons of drag.

### Adversarial self-checks

- **What if the design CL doubles (a climb)?** The twist optimal at cruise is no longer optimal: the basic (twist) load is fixed while the additional (α) load doubles, so the load is closer to the untwisted shape. Optimise at a mission-weighted set of CLs, not one.
- **What if the stall margin is demanded at a higher CL?** At CL_high → ∞ the additional load dominates and twist hardly helps; tip stall must then be fixed by the planform or the section (a more cambered tip section, stall strips, slots).
- **Which assumption fails first?** The uniform section cl_max of the stall criterion, then the lifting line near swept or tapered tips.
- **What measurement could falsify the model?** Tufts on a wing in a stall test show where separation starts; if a "root-stalling" design shows tip separation first, the uniform-cl_max assumption is wrong.
- **Does the result survive a change of discretisation?** More twist stations change e by about 2 × 10⁻⁴ beyond 4 stations; the optimum's character (washout concentrated outboard) does not change.
- **Which decision changes if the uncertainty doubles?** If the required stall margin were uncertain between 0 and 10 %, the tapered wing's induced drag would span 103–149 N; the recommendation (rectangular with washout) stands across that whole range.

## Check yourself

1. Near which taper ratio is the induced-drag factor δ of a straight-tapered wing smallest, and why does that not make it the best trainer wing?
2. Why does taper raise the local lift coefficient at the tips even when the load is nearly elliptic?
3. Why is induced drag a quadratic function of the twist variables?
4. What does the unconstrained twist optimum converge to as twist freedom grows?
5. What is the root bending coefficient of an elliptic load, C_B = M_root/(qSb)?
6. Which constraint makes Prandtl's bell optimal, and what span and drag does it give relative to the elliptic wing?
7. Which way does aft sweep move the first-stalling station, and why is that dangerous?
8. Why might the twist that is optimal at cruise be wrong for a climb?
