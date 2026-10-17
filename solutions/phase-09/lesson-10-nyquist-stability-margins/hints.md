# P09.10  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. The phase crossover of this servo does not depend on K; only \|L\| there does.
2. (a) The zeros of 1 + L are the closed-loop poles and its poles are the open-loop poles; encircling 0 by 1 + L is encircling −1 by L. (b) Near s = 0, L ≈ c/s^k; put s = εe^{jφ}. (c) Set the imaginary part of the denominator at s = jω to zero.
3. With P = 1 the stable loop needs one counterclockwise encirclement; a delay rotates every point of the curve clockwise by ωT.
4. The delay margin is PM divided by ω_gc, and both numbers move against you as K rises.
5. The phase of the lead network arctan(T_dω) − arctan(τω) has a maximum; differentiate it.
6. A delay adds phase lag at every frequency, so it moves the phase crossover too.
