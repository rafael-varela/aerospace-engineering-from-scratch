# P05.10  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. The support moment is the load times the arm. For a rectangle I = bh³/12 and c = h/2, so S = bh²/6: stress goes as 1/h², which is why turning the bar on its side hurts so much. A downward tip load on a cantilever puts the top edge in tension.
2. (a) Write ∫σ dA = 0 with σ = −Eκy and ask what it says about ∫y dA. (b) With x, y from the centroid the terms ∫x dA and ∫y dA vanish; what remains is a 2 × 2 system in b and c whose matrix is [[I_yy, I_xy], [I_xy, I_xx]]. (c) Each cap's own I about its centroid is bt³/12 = At²/12; add the parallel-axis term A(h_c/2)².
3. The tube's wall lies at distance R from the axis everywhere, so I = ∫y² dA with y = R sin φ and dA = Rt dφ. The axial stress is the same all round; the bending stress is ±MR/I at the two points in the plane of bending.
4. `brentq` needs a function that changes sign: 1.5 M (h/2)/I(b_f) − F_tu. When MS_ultimate = 0 the limit stress is exactly F_tu/1.5, whatever the load.
5. In (c), b = 0 means the strain plane is ε = cy. Write the moment conditions with that strain and solve for c and M_h; one of them is no longer free.
6. Each effect changes the root moment or the section. Root moments: elliptic 4Ls/(3π), uniform Ls/2. A uniform weight q gives qs²/2 at the root. Holes move the centroid up and remove I; recompute with the net-section rectangles.
