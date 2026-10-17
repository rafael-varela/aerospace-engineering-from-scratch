# P01.03  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. The two sums Σ aₖ/(eₖ + 1) = 0.0685083 and Σ aₖ/(eₖ + 2) = 0.0288033 do not depend on t or c; only the prefactors 10 t c² and c do. For (d) scale the lesson's 0.022 by t/0.12 (the error is linear in t) and by c² (it is an area).
2. For (b), write f(x) = f(x_m) + f′(x_m)(x − x_m) + ½f″(x_m)(x − x_m)² + … and integrate each term over [x_m − h/2, x_m + h/2]; odd powers of (x − x_m) vanish. For (c), every (x/c)^e becomes u^(2e) and dx becomes 2cu du; collect powers of u.
3. The thrust curve is linear between the given points, so any rule whose panel edges include the breakpoints integrates each piece exactly if the rule is exact for linear functions; check which of the three rules is, and why the left and right sums agree anyway here (look at F(0) and F(1.5)).
4. Signed errors: the first-order term is ∓(h/2)[f(c) − f(0)], the singular term is −K h^(3/2) for every rule. Averaging the left and right sums cancels the first term exactly. For (c), the accumulation function is least smooth where the integrand's derivative is largest.
5. ȳ = ∫ y_c · 2y_t dx / ∫ 2y_t dx is a weighted average of y_c with weights 2y_t; y_c is proportional to m, and a weighted average of a function that is zero at both ends and m at one point is less than m. The horizontal shift y_t sin θ moves the upper and lower surfaces in opposite directions, so the first-order area change cancels.
6. For each: is f bounded? Is it continuous? Does it vary on a scale smaller than h? Is the sum's error larger than the data's? One of them has no finite answer at all.
7. Near the leading edge y_t³ ∝ x^(3/2), so the singular error term is h^(5/2); the smooth term is (h²/24)[f′(c) − f′(0)] with f = (2/3)y_t³, and both f′(0) and f′(c) are small. Which term wins decides the order you observe.
