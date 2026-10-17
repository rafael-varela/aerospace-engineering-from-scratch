# P05.12  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. I of a tube is π(D⁴ − d⁴)/64; the tube's weight per metre is ρAg, multiplied by the load factor like the probe's. The two deflections add because the beam equation is linear.
2. (a) From EI w⁗ = q the four conditions are w(0) = w′(0) = 0 and w″(L) = w‴(L) = 0; no reaction is needed. (b) Integrating with R unknown gives two constants plus R: three unknowns for three kinematic conditions. (c) With u = 1 − t², t³√(1 − t²) dt becomes −½(1 − u)√u du.
3. Count the unknown reactions (root force, root moment, strut force) against the two equilibrium equations. The strut condition is w(a) = 0: the cantilever's deflection at a under the lift must be cancelled by the strut force acting at a.
4. Deflection scales as 1/EI when the shape of EI along the span is kept, so the factor is the ratio of the computed deflection to the allowed one.
5. For the clamped–clamped closed form, integrate EI w⁗ = q four times and use w = w′ = 0 at both ends. A point load between nodes is moved by up to h/2; how does that error scale with h?
6. ∫₀ˢ V dy = M(0) for a cantilever; shear strain V/(GA) integrated along the span gives the shear deflection. Compare each effect's sign with the sign of the discrepancy.
