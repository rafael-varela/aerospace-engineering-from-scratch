# P02.06  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Every quantity in (a) is a one-line function call; build the `case` dictionary with your numbers and pass it to `swing_test_error`. For (d), `max_amplitude_for_inertia_error` does the bisection for you.
2. For (b), the characteristic polynomial at θ = π is λ² + 2ζω₀λ − ω₀² = 0: the constant term is negative, so the product of the roots is negative. For (c), start from ½θ̇² = ω₀²(cos θ − cos θ₀), use cos θ = 1 − 2 sin²(θ/2), and let sin(θ/2) = k sin φ. For (d), you need k² to order θ₀⁴ only if you want the θ₀⁴ term; for the θ₀² term k² ≈ θ₀²/4 is enough.
3. sin θ cos θ = ½ sin 2θ. The equilibria are where sin 2θ = 0. Write the system as (θ, θ̇) and the Jacobian has −3n²k cos 2θ in the lower-left corner.
4. A peak is a sample larger than both its neighbours. The log decrement assumes the linear model; the amplitude term in the period does not care about damping and the damping term in the period is O(ζ²): compute both sizes before deciding which dominates.
5. I_cg = m g d T²/(4π²) − m d²; differentiate with respect to T and to d separately. For (c), add the three relative errors in quadrature to see the total, then look at which one is biggest.
6. For (c), with x = A cos t, ẋ = −A sin t, the integrand is μ(1 − A² cos² t) A² sin² t. You need ∫₀^{2π} sin² t dt = π and ∫₀^{2π} sin² t cos² t dt = π/4.
