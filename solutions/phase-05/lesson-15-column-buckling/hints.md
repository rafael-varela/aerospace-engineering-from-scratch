# P05.15  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. r² = I/A; for a tube r = √(d_o² + d_i²)/4. Compare KL/r with λ_c = √(2π²E/F_cy) before choosing a formula.
2. Write w = A sin kx + B cos kx + Cx + D, apply the four boundary conditions, and set the determinant of the 4 × 4 system to zero (or eliminate step by step). For (c), equate the two curves and their derivatives with respect to the slenderness.
3. Wrap the column curve in a function of the diameter that returns P_cr − 40 kN, with `tube_properties(d, 0.001)`, and bracket the root between 10 mm and 200 mm.
4. Richardson extrapolation removes the h² term: if P(h) = P + Ch², then (4P(h/2) − P(h))/3 = P + O(h⁴). For the mode shape, the coefficients follow from w(0) = w'(0) = 0 and tan β = β.
5. Assume the loaded shape is also a half sine of amplitude a and substitute; the sines cancel. For (b) solve σ(P) = F_cy with `brentq` between 0 and P_E.
6. Look up where L/r = 60 sits against λ_c for 4130, and think about what a weld does to straightness and to where the load line passes through the joint.
