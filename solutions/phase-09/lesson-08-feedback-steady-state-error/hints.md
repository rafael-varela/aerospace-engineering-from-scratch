# P09.08  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. L(0) = K_p/c; the error to a step of A is A/(1 + L(0)); for PI, K_v = K_i/c.
2. (b) Multiply E(s) = R/(1 + L) by s and take s → 0 with R = A/s^{q+1}; the factor sⁿ in L's denominator decides which limit is finite. (d) Try a loop whose closed-loop poles are on the imaginary axis or in the right half plane.
3. The plant 1/(Js²) already has two integrators; a disturbance torque enters *before* them, after the controller, so the controller's own integrators are what count.
4. The step error falls as 1/(1 + K_p/c), but the damping falls too; the critical gain bounds what you can try.
5. mg sin γ with γ in radians; the open-loop loss is F/c with time constant m/c.
6. \|S\| > 1 means feedback amplifies a disturbance at that frequency; a saturated actuator breaks the loop.
