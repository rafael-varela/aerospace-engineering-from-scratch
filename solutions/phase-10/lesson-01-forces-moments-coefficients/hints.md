# P10.01  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Re = Vc/ν and the scale factor multiply: to keep Vc the same with a chord five times smaller, the speed must be five times larger. Then compare that speed with the speed of sound.
2. (a) Write the force as one vector and project it on the stream direction and on its normal. (b) The moment of a force F at r about P is the cross product (r − r_P) × F; pick the component that is nose-up for x aft and z up. (c) Write the moment about x_ref as the moment about the a.c. plus N times the lever arm, and divide by qSc. (d) Hold N and A fixed: the balance measures them, the angle is what is uncertain.
3. The rotation uses the angle between the *balance* axes and the stream, not the chord's. The moment transfer is done in the balance axes, where the offsets are given; a moment is the same number in any axes.
4. `np.polyfit(x, y, 1, cov=True)` returns the slope first. The slope of c_m against c_n is (x_ref − x_ac)/c.
5. ∂c_d/∂α ≈ c_l, so the α contribution grows with lift; the q error scales c_d in proportion and c_d is small.
6. Ask which quantities depend on the boundary layer (and so on Re) and which are set mainly by the potential flow around the shape.
