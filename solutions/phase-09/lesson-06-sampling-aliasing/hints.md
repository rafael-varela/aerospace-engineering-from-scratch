# P09.06  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Find the multiple of 400 Hz nearest to each line; the alias is the distance to it.
2. (a) Put s = jω, so (ω/ω_c)² = −s²/ω_c², and solve (−s²/ω_c²)ⁿ = −1 for s. (b) Take logarithms of both inequalities and eliminate f_c. (c) Expand log H(jω) to first order in ω; for the Butterworth, a(s) = Π(s − p_k) so a₁/a₀ is a sum over the poles.
3. The stop edge is f_s − f_band in both cases; only the ratio f_st/f_p in (5) changes.
4. The prediction is the line amplitude times \|H\| at the *true* frequency, not at the alias. Then compare the 780 Hz prediction with the size of the noise in one DFT bin.
5. Each missing sample outside the record removes a sinc tail that decays only as 1/distance; a tone near Nyquist changes sign almost every sample, so its tails do not cancel.
6. Ask in each case: has the information already been destroyed at the sampler, and where does the delay come from?
