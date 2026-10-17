# P08.06  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. With no pressure gradient along the film, the Couette shear stress μU/h is the same everywhere in the gap; the torque is that stress times the shaft's surface area times the radius.
2. (a) Fully developed means ∂u/∂x = 0; continuity then gives ∂v/∂y = 0, and v = 0 on a wall makes v = 0 everywhere. (b) In cylindrical coordinates the viscous term for w(r) is μ(1/r)(r w′)′; the constant of integration multiplying ln r must vanish for a finite speed on the axis. (c) Write u_t and u_yy with the chain rule through η; the variables t and y must drop out of the resulting equation.
3. Put x down the slope: gravity has a component g sin θ along x (driving the film) and −g cos θ along y (setting a hydrostatic pressure across it). The free surface condition is μ du/dy = 0 at y = H.
4. The pressure-driven part of the force, Gh/2 per unit area times the land area, does not depend on the piston speed; the Couette part μU/h does. The leakage at fixed Δp scales as h³.
5. Write e^{−ky} cos(ωt − ky) as the real part of e^{iωt − (1 + i)ky} and differentiate once in t and twice in y. The start-up transient is a Stokes-first-problem-like wave that leaves through the far boundary.
6. Viscosity of hydraulic fluid falls as it warms; leakage is proportional to 1/μ and to h³. Compare the sizes of the factors before ranking them.
