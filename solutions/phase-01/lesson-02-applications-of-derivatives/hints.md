# P01.02  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Everything follows from K and √(C_D0/K). Keep C_D0 and K to five figures: (L/D)max is sensitive to their product. The glide angle is arctan(1/(L/D)max), and the sink rate is V* divided by (L/D)max in the small-angle picture.
2. For (a), write the sink as a constant times C_D/C_L^(3/2) and differentiate with the quotient rule; you only need the numerator. For (c), (L/D)max = ½(C_D0 K)^(−1/2) is a power of C_D0: differentiate twice as a power law, then use Taylor's second-order term ½f″Δx².
3. Write h(C_L) = C_L^(1/2)/(C_D0 + K C_L²); its derivative's numerator is ½C_L^(−1/2)(C_D0 + K C_L²) − C_L^(1/2)·2K C_L. Setting it to zero gives K C_L² in terms of C_D0 with a factor of 3 in a different place than before.
4. In a tailwind the tangent from (V_w, 0) with V_w < 0 touches the polar left of V*, but the polar's slope is zero at the minimum-sink speed and a tangent from a point on the axis must have negative slope there; think about where the tangent point can be. For (c), compare the failures with where (L/D)″ changes sign.
5. R = d/h, both decreasing: dR/dt = (ḋh − dḣ)/h². Substitute the two rates and factor out w/h; the sign is set by d/h against (V_b − V_w)/w.
6. Ask for each: which of "parabolic polar", "steady", "uniform wind", "fixed weight" is broken, and which quantity the pilot can see (airspeed, vertical speed, groundspeed) would disagree with the model.
7. The condition is w + u = V w′(V): the tangent to the shifted polar from the origin. Adding u to w changes F(V) by a constant and leaves F′ unchanged, so the Newton code changes by one term.
