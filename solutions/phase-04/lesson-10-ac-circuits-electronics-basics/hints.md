# P04.10  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Z = R + jωL; S = V²/Z*. A parallel capacitor supplies vars: it must absorb −Q, so ωCV² = Q.
2. (a) Divide numerator and denominator of Z_C/(R + Z_C) by Z_C. (b) Kirchhoff's current law at each capacitor node with admittances 1/R and jωC; eliminate node 1. For −3 dB set |1 + 3x + x²|² = 2 with x = jy and solve the quadratic in y². (c) Rule 2 makes v₋ = v_in; rule 1 makes R_f and R_g a divider.
3. Near f₀ the inductor and capacitor reactances cancel and only R limits the current that charges C; the voltage across C is Q times the input. Two reactive elements above f₀ give −40 dB/decade.
4. Use `aliased_frequency`, `adc_lsb(3.3, 10)` and the ratio 20 mV / LSB. A filter only has to reject what is actually there; a corner above Nyquist protects nothing between Nyquist and the corner.
5. Without noise every reading is the same count: the average of identical numbers is that number. Noise spreads the readings over neighbouring counts in proportion to where v₀ lies between them.
6. Ask in each case: is the information still there after this step, and how many LSB is the signal or the error?
