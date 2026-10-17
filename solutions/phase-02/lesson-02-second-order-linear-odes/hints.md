# P02.02  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. ζ = c/(2√(km)): compute √(km) first. For release from rest, v₀ = 0 and the sine coefficient in (5) is ζω_n x₀/ω_d, not zero.
2. In (c) the first-derivative terms of u cancel exactly because r = −ω_n is a *double* root; write out ẍ + 2ω_n ẋ + ω_n²x with x = u e^{−ω_n t} and watch them go. In (d) differentiate W, then replace ẍ₁ and ẍ₂ using the equation; the x₁x₂ terms cancel because both solutions have the same k/m. In (e) peaks of e^{−σt} cos(ω_d t − φ) are spaced by exactly one period.
3. The correspondence is L ↔ m, R ↔ c, 1/C ↔ k; then everything is the same formula. Critical: R² = 4L/C.
4. Build the table with a loop over `touchdown_response`, `max_stroke`, `strut_force`, `loses_contact`, `settling_time`; the band for the settling time is 2 % of the static deflection, about the static value (subtract x_s first).
5. A damper force proportional to ẋ|ẋ| is zero at zero velocity but the linear one is cv₀ at the first instant; the real peak load therefore comes later, near full stroke. For the bounce, think about what force the tyre can and cannot exert on the runway.
6. Differentiate x(t) and set ẋ = 0; collect the cos and sin terms. For (b), Newton's law for the mass is m Ẍ = mg − F, nothing more. For (c), the touchdown spike cv₀ is proportional to ζ and to v₀; the force at maximum stroke falls with ζ. Where two competing terms cross is where the minimum sits.
