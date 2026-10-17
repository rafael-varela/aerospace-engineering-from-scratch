# P05.18  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Work in SI (metres, newtons) and convert at the end. The stress formula is 6Pa/(bh²); the safe load inverts it. The self-weight load per metre is ρ g b h.
2. (a) Each integration brings a constant; the clamp conditions make both zero. Beyond the load M = 0, so w″ = 0. (b) Take logarithms before differentiating. (c) The clamp's rotation adds a rigid line θ₀x to the deflection; evaluate it at x = a and factor out a³/(3EI). (d) Superposition: the total deflection is the sum of the deflections under each load acting alone.
3. Convert the mass to a force first; E I for the rule is about 0.2 N·m², twenty times smaller than the aluminium strip's. For the 1 % mass, `scipy.optimize.brentq` on the ratio linear/elastica − 1.01. For the shape, compare the elastica and linear deflections at the top load against the reading noise of one point.
4. Copy the reduction lines from `main()`, replacing 518 with your seed. For (d), build w from `deflection_point_load` with no noise; `fit_profile` returns the uncertainty of θ₀ from the design of the fit and u_w, not from the data.
5. The terms add in quadrature, so reducing a term that is already small gains nothing. Find the term that is largest after the dial gauge and what halves it.
6. All three candidates scale the apparent E by a factor. Set each factor equal to 56/69 and solve for the size of the error. Then ask which error is easiest to make with the measuring method they describe.
