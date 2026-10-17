# P09.07  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. The zeros of an N-point average are at multiples of f_s/N; 180 Hz is nulled only if 180 N/f_s is a whole number.
2. (a) Factor e^{−jωN/2} out of the numerator and e^{−jω/2} out of the denominator of the summed series. (b) Pair the tap k with the tap N − 1 − k. (c) The trapezoidal rule is an integrator; compare \|K + s\| with \|K − s\| when Re s < 0.
3. A notch has zeros on the unit circle at e^{±j2π·180/1000} and poles just inside them; Q sets how close.
4. The RMS error against the true rate contains the delay: a delay T applied to a motion of amplitude A and frequency f costs up to 2πfAT.
5. Invert the warping law: the analog corner f_c lands at (f_s/π) arctan(πf_c/f_s).
6. (a) Poles near z = 1 make A(1) = Σa_k a small difference of large numbers. (c) A filter that may look into the future can cancel its own delay.
