# P05.06  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Area first: πd²/4 in mm² with the force in N gives MPa directly. The diameter changes by the lateral strain −νε times d. True stress is σ(1 + ε); at strains of a few thousandths, ask whether the difference is bigger than the uncertainty of E.
2. (b) Write F = σ_t A, set dF = 0 and use dA/A = −dε_t. For the power law, differentiate Kε_tᵐ and set it equal to itself. (d) Put the principal stresses along the diagonals; each diagonal's strain has a direct part and a Poisson part, and the shear strain of the square equals the sum of the diagonal strains.
3. Cut the sphere through its centre: the pressure on the projected area πR² is carried by the wall stress over the ring 2πRt. Both tangential stresses are equal, the radial stress is negligible next to them, so `hooke_strain(np.diag([s, s, 0]), E, nu)`.
4. `considere_point` gives the ultimate the generator was built around; the yield of the record is lower than 880 MPa for the same reason as in the lesson (engineering stress is σ_t/(1 + ε)). A sharper knee (larger n) means the plastic strain in the modulus window is tiny.
5. Equation (7) gives c from two moduli. For (c), write the recorded strain as σ/E_app + ε_p and put it into the offset-line equation. For (d), call `offset_yield(x, y, 73.1e3, r.eps0)`.
6. (a) Compare 52 GPa with the modulus of every aluminium alloy you know, then look at the compliance table of part 4. (b) Mass = energy / (toughness per unit mass); divide the toughness by the density. (c) Think about what "typical" says about the next batch. (d) Forming a sheet into a curved skin takes it far beyond 1 % strain.
