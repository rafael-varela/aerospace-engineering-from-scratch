# P06.14  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. det(K − ω²M) is a quadratic in ω²: write its coefficients from the four entries. Each shape comes from either row of (K − ω²M)φ = 0.
2. (a) Pre-multiply each eigen-equation by the other mode and subtract. (b) r = ΦΦᵀMr because ΦΦᵀM is the identity on the space spanned by the modes. (c) Each modal coordinate is an SDOF system driven by −Γᵢÿ; the absolute acceleration is ÿ r + Φ η̈.
3. A free body can translate without stretching any spring: that is a zero-frequency mode. Removing k₁ changes every mode, not only the first.
4. A step input to an SDOF system gives a peak of twice the static response when undamped; in a chain, the modal peaks can coincide in time in some springs.
5. Look at how much the service and payload modules move *relative to each other* in each mode: a damper only dissipates energy where it is stretched.
6. A model is a hypothesis; a low-level sweep is the measurement that tests it before the full-level test can damage anything.
