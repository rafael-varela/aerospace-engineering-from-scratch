# P01.17  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. AᵀA for a line fit is [[m, Σα], [Σα, Σα²]] and Aᵀ**b** = (ΣC_L, Σα C_L). With α = 0, 4, 8 the sums are tiny to do by hand.
2. For (a), expand the quadratic and use ∇(**x**ᵀM**x**) = 2M**x** for symmetric M. For (d), the diagonal of AᵀA is 1 + ε² and 1 + 10⁻¹⁶ is not representable in double precision.
3. Complete the square: c₀ + c₁C_L + c₂C_L² = (c₀ − c₁²/4c₂) + c₂(C_L + c₁/2c₂)². A residual sign pattern of +,+,+,−,−,…,−,+,+,+ is the signature of a missing linear term.
4. The standard error of e is e × (standard error of K)/K. Standard errors go as 1/√m, so the required m scales as (current SE/target SE)².
5. For (a) think about which direction the stalled points pull the parabola and what sign the low-α residuals then have. For (c), `condition_squaring` gives cond(AᵀA); the number of digits lost is about log₁₀ of it.
6. CGS2's loss is at rounding level. Lesson 2.5 blamed the explicit-Qᵀ**b** failure on the loss of orthogonality of Q divided by the tiny R diagonal; with an orthogonal Q that first factor is gone, so predict what is left before you run it.
