# P10.11  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. b = 1, c = 0.2, sweep 45°, `lift_slope_and_ac(make_lattice(..., nx=1, ny=8, spacing="uniform"))`. 2. Loop ny and both spacings; `observed_order` on the last three. 3. Compare with the finest Hough value.

## 2.
(a) Γ/(2πd) = Vα; cl = 2Γ/(Vc); match 2πα and the centre of pressure. (b) Vortices at c/8 and 5c/8, control points at 3c/8 and 7c/8; 2 × 2 system; sum Γ. (c) My(x₂) = My(x₁) + (x₂ − x₁)ΣF_z; F = ΓV × l gives ΣF_z = L cos α. (d) A 2-D vortex field in the y–z plane, normal component on each strip, D = −½ρΣGv·nΔs.

## 3.
1. `make_lattice(10, 10/6, sweep_c4=±30°, nx=4, ny=40)`; `lift_slope_and_ac`; `analyse` at 5°. 2. The strip of largest cl/CL on the right wing. 3. Bending of a forward-swept wing twists its tips nose-up (divergence).

## 4.
1. Three lattices, `observed_order`. 2. Richardson; CL_α = CL/sin(4°). 3. Uniform spacing; compare orders and errors.

## 5.
1. `np.genfromtxt(..., delimiter=",", comments="#", names=True, dtype=None, encoding=None)` or the csv module; loop over rows. 2. Differences in per degree and per cent; compare with ±0.0005. 3. Side-edge vortex lift (nonlinear) and the linear-slope reading.

## 6.
(a) Induced drag with each e at the cruise CL. (b) e > 1 is impossible for a planar wing; tip inset or finer lattice. (c) Profile drag (`P10.06`, `P10.13`), cl_max and stall (`P10.15`), surface pressures with thickness (`P10.05`).
