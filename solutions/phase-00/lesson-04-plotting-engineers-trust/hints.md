# P00.04  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Ask of each label: which quantity, which unit, and could the reader confuse this with another quantity that shares the symbol?
2. Substitute `ρ = p / (R T)` into the hydrostatic equation and separate variables; a decade is a factor of 10, so divide a natural-log slope by `ln 10`.
3. Take logarithms of `D = ½ ρ V² S CD` with `CD` constant; what is the exponent of `V`?
4. `check_figure` only sees what the Matplotlib objects record. A flaw it cannot see is one that is not stored in any object property; think about what "joined by a line" means to a `Line2D`.
5. The bars are ±1 σ and the noise is Gaussian with exactly that σ. What fraction of a normal distribution lies outside ±1 σ?
6. Every rule in the style exists to protect a particular reading. A rule should break when that reading is not the one that matters, or when it is impossible.
7. The rule is really "the unit is unambiguously marked". Think about accepting `[unit]` *anywhere* in the label, and about what a false positive costs compared with a missed defect.
