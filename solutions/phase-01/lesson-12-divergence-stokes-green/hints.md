# P01.12  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Divergence adds the three "own-coordinate" partials; curl mixes the others. A divergence-free field has zero net flux through *every* closed surface, whatever its shape. For **G**, each component depends on its own coordinate at most quadratically: think about what the midpoint rule does to a linear integrand.
2. For (a), integrate ∂Q/∂x in x first with the fundamental theorem, then recognise the two boundary terms as the line integrals up the right side and down the left side. For (b), write the rate of change of mass as minus the outward mass flux and convert the surface integral to a volume integral. For (c), the side walls of a stream tube carry no flux because **v** · **n** = 0 there.
3. **g** is linear in **r**, so its divergence is a constant times 3. A constant divergence integrates to (constant) × (volume), and ρ × volume is the mass. At the surface of a uniform sphere g = (4/3)πGρR; invert for ρ.
4. Fit slopes with `np.polyfit(np.log(n), np.log(res), 1)`. For (b), write the midpoint sum over the column of cells in x for a fixed (y_j, z_k) and compare with F_x(x₁) − F_x(x₀) at the same (y_j, z_k): when is Σ f′(x_i) Δx exactly f(x₁) − f(x₀)? For (c), the loop never touches the vorticity step; the area integral does.
5. The flux of a curl through a closed surface is zero because ∇ · (∇ × **F**) = 0 and the divergence theorem. On the hemisphere **n** = **r**/R, so (0, 0, 2ω) · **n** = 2ωz/R. Swapping parameters reverses the cross product.
6. Ask separately: is ∂ρ/∂t zero? Is ρ constant? Is the surface closed and impermeable? Each failure maps to one of the four cases.
7. Apply the divergence theorem to **F** = **r**/3 (divergence 1) and to **F** = (x²/2, 0, 0) (divergence x). For a flat triangle, **n** dA = ½(**p**₂ − **p**₁) × (**p**₃ − **p**₁) and **r** · **n** is the same at every point of the plane.
