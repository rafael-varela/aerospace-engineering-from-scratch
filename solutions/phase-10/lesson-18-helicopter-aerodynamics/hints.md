# P10.18  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. `hover_power(Rotor(mass=...), rho=...)` takes the density; the induced power goes as T^{3/2}/√ρ while the profile power goes as ρ.
2. (a) The local lift per unit span is ½ρ(Ωr)²c a(θ − λ/r); sum B blades and divide by ρA(ΩR)². (c) Average sin ψ and sin²ψ over a revolution: 0 and ½. (d) For μ ≫ λ, √(μ² + λ²) ≈ μ.
3. With the mass fixed, C_T ∝ 1/(ΩR)² but C_T^{3/2}(ΩR)³ is independent of ΩR; the profile power is σC_d0/8 × ρA(ΩR)³.
4. The rate of climb is the excess power divided by the weight; the maximum level speed is where the required power equals the available power (`brentq`).
5. In the reverse-flow region u_T < 0, so u_T³ is negative there but the drag still absorbs power, ∝ |u_T|³.
6. At V_mp the profile power is about half the total, at hover about a quarter: the inputs matter in different places.
