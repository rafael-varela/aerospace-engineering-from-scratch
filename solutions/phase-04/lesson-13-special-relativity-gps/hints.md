# P04.13  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Both terms are fractional rates; multiply each by 36 000 s. For (b), expand 1/R − 1/(R + h) for h ≪ R and you get gh/(c²) with g = GM/R².
2. (a) Substitute x = ct, x′ = ct′ into both equations, divide one by the other, and solve for A². (b) Expand γ²(cΔt − βΔx)² − γ²(Δx − βcΔt)²; the cross terms cancel. (c) Write each boost as a matrix in (ct, x) and multiply; compare with the single-boost matrix. For |w| < c look at c − w and c + w. (d) cos 0 = 1, cos π = −1, cos π/2 = 0; γ(1 − β)γ(1 + β) = 1.
3. `gps_clock_budget(R + h, mu=GM_MARS, R=R_MARS, omega=0.0, J2=0.0)`. The crossover condition GM/2a = GM(1/R − 1/a) does not contain GM. Kepler's third law gives the areostationary radius.
4. `eccentricity_correction(e, pi/2)` gives the amplitude; n = √(GM/a³). For (b), the satellite position is a(cos θ, sin θ sin p, 0) with cos θ = R/a.
5. ∫(1/r)dt over one period = ∫ dE/(an) from 0 to 2π. Vis-viva: v²/2 = GM/r − GM/2a, so v²/2 + GM/r = 2GM/r − GM/2a.
6. Divide the requirement by each term's range rate (fractional rate times c) for the secular terms; for the periodic ones compare the bounded amplitude with the requirement.
