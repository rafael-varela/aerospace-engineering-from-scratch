# P00.12  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Volume × density for (a); count × per-person mass for (b); $\tfrac12 m V^2$ for (c), then divide by 43 MJ. Your error factor is roughly the product of the uncertainty factors of your inputs.
2. Differentiate $T/W$ with respect to $C_L$ once. Then invert $C_L^* = W/(\tfrac12 \rho V^{*2} S)$ for $V^*$; the scaling with each variable is read off the square root.
3. Write the units of $V/c$ and of $V/(cg)$; one of them is a length. The limiting case: with no fuel burned the range must be zero.
4. The flags that bound the passing interval are Mach at the top and $C_L$ at both ends; find where each crosses its threshold. For (b), convert $C_L^*$ to a speed with the ISA density at 11 km.
5. Every flag is a true statement about the drone; none is an error in its data. The checker's ranges are the problem, not the aircraft.
6. Concorde's Mach and the turboprop's altitude and Mach will trip. Each option has a cost: (i) lets the "speed in km/h" mistake through for airliners; (ii) needs the user to know the class; (iii) needs a human to read every flag.
7. The expected effect of one thruster firing on the trajectory can be estimated from the thruster's rated thrust, the firing duration, and the spacecraft mass; compare that with the trajectory change the navigators actually inferred. A factor of 4.45 is far outside any reasonable uncertainty.
