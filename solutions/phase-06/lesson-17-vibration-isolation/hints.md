# P06.17  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Four mounts in parallel add their stiffnesses. The sag under n g is n g/ωₙ², and it does not depend on the mass once fₙ is fixed.
2. (a) Divide the complex amplitude equation by k and write r = ω/ωₙ. (b) Square T and compare numerator and denominator as functions of r². (c) Differentiate T² with respect to u = r²; the result is a quadratic in u. (d) The damper force cω grows with frequency; a hysteretic force does not.
3. Two requirements pull fₙ in opposite directions. Compute each bound before you look for a compromise.
4. A binding requirement is one that is met with (almost) no margin at the optimum and would be broken by moving towards a better Grms.
5. Damping helps when the input's energy is near fₙ (where it lowers the resonance) and hurts when the energy is far above fₙ (where it raises T ≈ 2ζ/r).
6. A viscous model forced through hysteretic data will look acceptable at resonance and fail at high frequency; the χ² and the high-frequency slope show where.
