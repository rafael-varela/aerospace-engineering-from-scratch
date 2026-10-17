# P02.07  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. |z|⁵ and 5 arg z; subtract 2π from the argument until it lands in (−π, π]. The product of the n-th roots of z is (−1)^{n−1} z for n odd... check: the roots are the zeros of wⁿ − z, so their product is −(−z) = z when n = 3.
2. For (c), write cos ωt = Re(e^{iωt}) and move Re outside each term; that step needs real coefficients. For (d), minimise the denominator D(r) = (1 − r²)² + 4ζ²r² by setting dD/dr² = 0.
3. 0.9 cos(2πf₂t) = Re(0.9 e^{i2πΔf t} e^{i2πf₁t}): the second phasor rotates slowly at Δf = f₂ − f₁. Then the two-source formula (3) of the lesson applies with δ = 2πΔf t.
4. Floor displacement amplitude is acceleration/ω². For the base-excited ODE, the right-hand side c ẏ + k y is a known function of t; write the state as (x, ẋ) and pass y and ẏ through the lambda.
5. Minimising a sum of two cosine-rule amplitudes has no closed form; a sweep of δ in 1° steps is a perfectly good method. Think about what a 80 Hz line and a broadband background would do to a *measured* minimum.
6. For (d), Σ_{k=0}^{N−1} q^k = (1 − q^N)/(1 − q) with q = e^{iδ}; then factor e^{iNδ/2} from the numerator and e^{iδ/2} from the denominator and use e^{iα} − e^{−iα} = 2i sin α.
