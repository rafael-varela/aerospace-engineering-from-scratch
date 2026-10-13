# P01.04  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. s = b/2, τ = 0.5. 2. mass = σ c_r s (1 + τ). 3. ȳ = s(1 + 2τ)/(3(1 + τ)). 4. I_xx = 2σ c_r s³[1/3 − (1 − τ)/4]. 5. `roll_inertia_numeric(lambda y: chord_linear(y, c_r, 0.5, s), sigma, s, 2)`; the integrand y² c(y) is a cubic and Simpson is exact for cubics. 6. Rectangle of equal mass: c̄ = mass/(2σ s), I = 2σ c̄ s³/3; radius of gyration √(I_xx/mass).

## 2.
(a) u = y², dv = e^(−y/ℓ)dy → v = −ℓe^(−y/ℓ); boundary term −ℓ y² e^(−y/ℓ) plus 2ℓ∫y e^(−y/ℓ); repeat with u = y; evaluate from 0 to s (every boundary term vanishes at 0 except the constant 2ℓ³). (b) y = s sin θ, dy = s cos θ dθ, √(1 − sin²θ) = cos θ ≥ 0 on [0, π/2]; sin²θ cos²θ = ¼sin²2θ = ⅛(1 − cos 4θ); integrate. (c) T − k v² = (√T − √k v)(√T + √k v); write 1/(…) = A/(√T − √k v) + B/(√T + √k v), find A = B = 1/(2√T); integrate each to a log; combine. Expansion: x = v/v_t, ln((1 + x)/(1 − x)) = 2x + 2x³/3 + …; multiply by m/(2√(Tk)) and substitute v_t. (d) S_n = (2M_n + T_n)/3 with M_n = I − (h²/24)Δf′ + O(h⁴), T_n = I + (h²/12)Δf′ + O(h⁴); the h² terms cancel; the surviving O(h⁴) term is proportional to f⁽⁴⁾ (state from §3.6), which vanishes for cubics.

## 3.
(a) ∫₀^∞ e^(−st)dt = 1/s (s > 0). ∫ t e^(−st)dt by parts = 1/s² (s > 0). ∫ e^(−(s+a)t)dt = 1/(s + a) (s > −a). ∫ e^(−st) sin ωt dt by parts twice gives I = ω/s² − (ω²/s²)I, so I = ω/(s² + ω²) (s > 0). (b) `improper_integral(lambda t: np.exp(-s*t)*f(t), 0.0, np.array([2, 5, 10, 20]), n=4000)` for each f. (c) f = eᵗ, s = 2: integrand e^(−t), converges to 1; s = 0.5: integrand e^(0.5t), the sequence grows without bound. (d) The transform exists only for s larger than the fastest growth rate of f; the region of convergence is where the system's response is defined.

## 4.
(a) `roll_inertia_numeric` for each rule and n; `np.polyfit(np.log(n), np.log(err), 1)`; the singular term ∝ h^(3/2) at the tip is the same for every rule and dominates every rule's own h² or h⁴ term. (b) Write the θ-integrand and call `midpoint`, `trapezoid`, `simpson` with n = 4; compare. The midpoint rule is also a periodic equal-weight rule (shifted by h/2) and is exact too; Simpson's 1-4-2-4-1 weights are not translation-invariant on the period. (c) `trapezoid_error_bound(2 * 2 * sigma * c_r, 0, s, n)` and the measured error; the ratio is ∫f″/((b − a) max|f″|) × (constant), a property of f alone. (d) Solve (b − a)³ M₂/(12n²) ≤ 10⁻⁴ I and (b − a)⁵ M₄/(180n⁴) ≤ 10⁻⁴ I for n (round up, even for Simpson); run and compare.

## 5.
(a) ∫₀ˢ √(1 − y²/s²) dy = s π/4 (y = s sin θ, ∫cos²θ dθ = π/4), so L₀ = 2W/(πs). (b) M₀ = L₀ s² ∫₀^(π/2) sin θ cos²θ dθ = L₀ s²[−cos³θ/3]₀^(π/2) = L₀ s²/3; substitute L₀. (c) Tapered: (W/2) × s(1 + 2τ)/(3(1 + τ)); uniform: (W/2)(s/2). (d) `simpson(lambda y: L0*np.sqrt(1-(y/s)**2), 2.0, s, n)` and with (y − 2); error estimate from the n^(−1.5) law or from halving n.

## 6.
(a) Add 2 × m_fuel × y² to I_xx; compare with the wing's own value. (b) Replace the inner section by a rectangle (polynomial correction) and add a fuselage estimate (a cylinder of radius r: m r²/2 about its axis is small; its mass is close to the axis). (c) A zero tip chord with √ behaviour between the last two points: the trapezoid rule on the last panel underestimates by a known fraction; estimate it. (d) m dv/dt = T(v) − k v² − μ(W − L): the new v_t solves T(v) = k v² + μ(W − L(v)); with T linear in v the denominator is still a quadratic and partial fractions still apply, with different roots.

## 7.
(a) Linear: 2σ′c_r² s³ ∫₀¹ u²[1 − (1 − τ)u]² du, expand and integrate term by term. Elliptic: 2σ′c_r² ∫₀ˢ (y² − y⁴/s²) dy = 2σ′c_r² s³(1/3 − 1/5). (b) Substitute; masses 2σ′∫c² dy by the same expansions; k = √(I/m). (c) `simpson` with n = 2 vs the closed form; the error is exactly (s⁵/(180 × 2⁴)) × f⁽⁴⁾ with f = 2σ′ y² c² (f⁽⁴⁾ = 2σ′ × 24 × c_r²(1 − τ)²/s² is constant); find n from the bound. (d) Match the law to where the mass is: uniform foam → λ ∝ c × thickness ∝ c²; a single spar sized for root bending → heavier than c² at the root; ballast → add point/line masses to whichever law.
