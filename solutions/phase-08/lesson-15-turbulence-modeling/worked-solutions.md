# P08.15  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe scripts that call it, given to the figures shown. "Model" means the mixing-length closure of the lesson with κ = 0.41 and A⁺ = 26; Dean's and Blasius's formulas are correlations.

## 1. Wall units of a duct (s = 15: gap 55 mm, h = 27.5 mm)

**Assumptions.** Fully developed, smooth walls, infinite span, air at 20 °C (ν = 1.503 × 10⁻⁵ m²/s).

(a) Re_b = 4 × 0.055/1.503 × 10⁻⁵ = **14 634**.

| closure | Re_τ | u_τ (m/s) | τ_w (Pa) | dP/dx (Pa/m) | Δp over 2 m (Pa) |
|---|---|---|---|---|---|
| Prandtl outer | 437.4 | 0.2391 | 0.0688 | −2.503 | 5.01 |
| Nikuradse outer | 405.0 | 0.2214 | 0.0590 | −2.146 | 4.29 |

(b) Dean: c_f = 0.073 × 14 634^(−1/4) = 0.006637, τ_w = 0.006637 × ½ × 1.204 × 16 = **0.0639 Pa** (Δp over 2 m 4.65 Pa). The model is **+7.7 %** (Prandtl) and **−7.7 %** (Nikuradse): it brackets the correlation.

(c) y⁺ = 1: ν/u_τ = **63 µm** (Prandtl) or 68 µm; y⁺ = 30: **1.9 mm** or 2.0 mm.

**Check.** The momentum balance h|dP/dx| = τ_w holds to round-off (a test checks it): 0.0275 × 2.503 = 0.0688.

**Common mistakes.** Using h instead of 2h in Re_b (Dean's correlation uses the full height); quoting one closure's number as "the" answer.

## 2. Derivations

(a) Averaged x-momentum, fully developed: 0 = −dP/dx + dτ/dy with τ = μ dU/dy − ρ⟨u′v′⟩. Integrate from the wall: τ(y) = τ_w + y dP/dx (at the wall ⟨u′v′⟩ = 0 by no slip). Symmetry: τ(h) = 0 ⇒ dP/dx = −τ_w/h, so τ/τ_w = 1 − y/h = **1 − y⁺/Re_τ** (since y/h = y⁺/Re_τ). Control volume of length L and unit span between the walls: pressure force on the ends 2h(P₁ − P₂) = 2hL(−dP/dx); shear on two walls 2Lτ_w; momentum flux in = out (fully developed). Hence −dP/dx = τ_w/h, the same result.

(b) With −⟨u′v′⟩⁺ = l⁺²g² (g ≥ 0): g + l⁺²g² = τ⁺. Roots g = [−1 ± √(1 + 4l⁺²τ⁺)]/(2l⁺²); take +. Multiply numerator and denominator by 1 + √(1 + 4l⁺²τ⁺): numerator (1 + 4l⁺²τ⁺) − 1 = 4l⁺²τ⁺, so **g = 2τ⁺/(1 + √(1 + 4l⁺²τ⁺))**. Limits: l⁺ → 0 gives g = τ⁺ (laminar); l⁺²τ⁺ ≫ 1 gives g → √τ⁺/l⁺ (fully turbulent). The textbook form subtracts two numbers close to 1 when l⁺ is small (in the sublayer l⁺ ~ 10⁻³), losing most significant digits, and it is 0/0 at l⁺ = 0; the rationalised form has no cancellation.

(c) Constant stress: ν_t dU/dy = u_τ². Equilibrium: P = ν_t(dU/dy)² = ε. With dU/dy = u_τ/(κy): ν_t = κu_τy; ε = u_τ² · u_τ/(κy) = **u_τ³/(κy)**; ν_t = C_μk²/ε ⇒ k² = κu_τy · u_τ³/(κy)/C_μ = u_τ⁴/C_μ ⇒ **k = u_τ²/√C_μ**. ε equation: diffusion d/dy[(κu_τy/σ_ε)(−u_τ³/(κy²))] = d/dy[−u_τ⁴/(σ_ε y)] = u_τ⁴/(σ_ε y²). Source: (C_ε1 − C_ε2)ε²/k = (C_ε1 − C_ε2)(u_τ⁶/(κ²y²))(√C_μ/u_τ²) = (C_ε1 − C_ε2)√C_μ u_τ⁴/(κ²y²). Sum zero ⇒ 1/σ_ε = (C_ε2 − C_ε1)√C_μ/κ² ⇒ **κ² = σ_ε√C_μ(C_ε2 − C_ε1)**.

**Common mistakes.** Dropping the absolute value in l²|dU/dy| (harmless here, since dU/dy ≥ 0 in the lower half); forgetting that k is constant in the log layer, which makes the k equation trivially satisfied.

## 3. The log-law friction law

(a) (1/Re_τ)∫₀^{Re_τ}[(1/κ)ln y⁺ + B]dy⁺ = (1/(κRe_τ))[Re_τ ln Re_τ − Re_τ] + B = **(1/κ)ln Re_τ + B − 1/κ**. Ignored: the viscous and buffer layers, where the real u⁺ lies below the log law (the integral counts the log law's large negative values near y⁺ → 0, which roughly mimics this, so the error is small); and the wake, where the real u⁺ rises above the log law towards the centreline (raises U_b⁺, lowers c_f).

(b) Re_b = 2 × 10⁴: Re_τ = 556.2, **c_f = 0.006188**, Dean 0.006139 (**+0.8 %**). Re_b = 10⁵: Re_τ = 2329, **c_f = 0.004338**, Dean 0.004105 (**+5.7 %**). The missing wake matters more at high Re.

(c) κ = 0.4327: **c_f = 0.006624** and **0.004671**, 7.0 % and 7.7 % higher than with κ = 0.41. A smaller log-law slope means a slower velocity rise across the layer, a lower U_b⁺ and more friction. A k–ε code with wall functions that impose κ = 0.41 at the first cell while its interior equations produce κ = 0.433 has a small inconsistency in the log layer; its friction then depends a little on where the first cell sits, which is one reason wall-function results drift with y⁺ of the first cell.

**Check.** At Re_τ = 500 the log-law friction law is within 0.7 % of Dean (the lesson's sweep).

## 4. Which constant matters? (Re_τ = 1000)

(a)

| case | c_f | against Dean | B (fit) |
|---|---|---|---|
| κ 0.41, A⁺ 26, Prandtl (baseline) | 0.00560 | +7.0 % | 4.90 |
| κ = 0.38 | 0.00512 | −1.1 % | 4.69 |
| κ = 0.42 | 0.00577 | +9.8 % | 4.96 |
| A⁺ = 22 | 0.00611 | +15.5 % | 4.39 |
| A⁺ = 30 | 0.00518 | −0.1 % | 5.31 |
| Nikuradse outer | 0.00477 | −7.0 % | (wake) |

(b) A⁺ (±8 % for ±4 in A⁺: it sets the buffer layer and hence B, the offset of the whole log profile) ≈ the outer closure (−15 % from Prandtl to Nikuradse: it sets the wake and the core velocity) > κ (about 3 % per 0.01: the slope of the log layer).

(c) Larger: the model's spread over plausible constants is about ±8 %, against a few per cent for the correlation. A duct designed from the model alone needs that ±8 % margin; a designer with Dean's correlation (for a straight channel) has a better number, and the model's value is in cases the correlation does not cover.

**Common mistakes.** Fitting the log law over the wake (B then means nothing; the Nikuradse fit gives κ = 0.36).

## 5. A propellant feed line is a pipe

(a) A cylinder of radius r and length L: pressure force πr²(P₁ − P₂) = πr²L(−dP/dx) balances the shear 2πrLτ(r) (fully developed, no momentum change), so τ(r) = (r/2)(−dP/dx), linear in r, with τ_w = (R/2)(−dP/dx). Then τ/τ_w = r/R = 1 − y/R. Flow rate Q = ∫₀^R U 2πr dr; with r = R − y, U_b = Q/(πR²) = (2/R²)∫₀^R U(R − y)dy, and in wall units **U_b⁺ = (2/R⁺²)∫₀^{R⁺} u⁺(R⁺ − y⁺)dy⁺**.

(b)

| R⁺ | closure | Re_D | c_f (model) | Blasius | difference |
|---|---|---|---|---|---|
| 500 | Prandtl | 15 995 | 0.00782 | 0.00703 | +11.1 % |
| 500 | Nikuradse | 17 125 | 0.00682 | 0.00691 | −1.4 % |
| 2000 | Prandtl | 78 533 | 0.00519 | 0.00473 | +9.8 % |
| 2000 | Nikuradse | 83 312 | 0.00461 | 0.00466 | −1.0 % |

(c) Nikuradse's outer length, within 1.5 % of Blasius. It was fitted to Nikuradse's own pipe measurements, so it should do well in a pipe; in a channel (part 4 of the lesson) it underpredicts by about 7 %. A closure is only as general as its calibration.

(d) u_cl/U_b = **1.250** (R⁺ = 500) and **1.198** (2000) with the Nikuradse length (1.177 and 1.135 with Prandtl's). A meter on the axis reads 20–25 % above the bulk velocity, and the correction depends on Re and on the profile model: calibrate it in place or use the model's ratio with its uncertainty.

## 6. Meshing the duct for CFD

(a) y⁺ = 1 at **63 µm** (Problem 1). With growth ratio 1.2, the cell count to the centreline (27.5 mm) is n = ln(1 + 0.2 × 27.5/0.063)/ln 1.2 = **25 cells** on each half, about 50 across the gap, more in practice to resolve the bend.

(b) The first cell centre in 30 ≤ y⁺ ≤ 0.1–0.2Re_τ, here 1.9 mm to about 4–9 mm (Re_τ ≈ 440). At Re_τ = 180 the log layer spans y⁺ = 30 to 18–36: there is almost no room, and a small change in speed moves the first cell out of it. Wall functions are for high Re; at low Re use a low-Re model.

(c) Separation on the inside of the bend and its reattachment; the secondary (Dean) vortices from curvature, whose strength the eddy-viscosity model under-predicts because it ignores stress anisotropy and curvature effects on the turbulence; the slow recovery downstream; any relaminarisation in a strong acceleration.

(d) The model's own uncertainty in straight-duct friction is about ±8 % (Problem 4), larger than 2 %. A 2 % answer needs a measured pressure drop on a prototype duct or a validated correlation, with the CFD used for the differences between designs, not their absolute value.

## Check-yourself answers

1. The Reynolds stress −ρ⟨u′_i u′_j⟩, the average of products of fluctuations; it is a new unknown, and averaging produces no equation for it.
2. Because the mean momentum balance, integrated from the wall, gives τ = τ_w + y dP/dx with dP/dx constant; symmetry fixes dP/dx = −τ_w/h. No closure enters.
3. ν_t = l²|dU/dy|; near a wall the largest eddies that can mix fluid across a distance y are of size proportional to y.
4. It reduces the mixing length to about κy⁺²/A⁺ in the viscous sublayer, giving u⁺ = y⁺ there and a buffer layer; it sets the log-law intercept B.
5. It avoids the subtraction of nearly equal numbers when l⁺ is small and is well defined at l⁺ = 0.
6. k = u_τ²/√C_μ, ε = u_τ³/(κy), ν_t = κu_τy.
7. κ = √(1.3 × 0.3 × 0.48) = 0.4327.
8. Low-Re model: y⁺ ≈ 1 at the first cell. Wall functions: first cell centre at 30 < y⁺ < about 0.1–0.2 Re_τ.
