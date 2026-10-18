# P08.07  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. ε/D is the only roughness number Colebrook needs. The fan power is the volume flow times the pressure drop it overcomes.
2. (a) Substitute V = 4Q/(πD²) and Re = ρVD/μ into f(L/D)ρV²/2 and watch ρ cancel. (b) Differentiate g twice; a concave increasing function crosses zero once, and a Newton step from either side lands on the left of the root (where g ≤ 0), so the later iterates increase. (c) At fixed ṁ, V ∝ D⁻² and Re ∝ D⁻¹.
3. Liquid oxygen is about ten times less viscous than RP-1 and a little denser. At very high Re the friction factor stops depending on Re: the roughness takes over.
4. The loss falls roughly as D⁻⁵, so a small change in bore makes a large change in loss; the inside diameter is the outside diameter minus two walls.
5. Fixed-point iteration converges linearly with rate |g′(x*)| of the map; Newton squares the error. In the hydrazine line, look at where Re falls relative to 2300 and 4000.
6. Change one input at a time with `feed_line_budget` and compare the changes in p_tank with each other, not with the tank pressure itself.
