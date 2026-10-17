# P06.15  Continuous Systems

> A beam has infinitely many modes, but a transcendental equation and a handful of cubic elements will give you the first few to six digits, and tell you how far to trust the rest.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the Euler–Bernoulli beam vibration equation EI w⁗ + ρA ẅ = 0, separate variables, and obtain the frequency equations and mode shapes of cantilever, pinned–pinned and clamped–clamped beams, including a cantilever with a tip mass. 2. Solve the frequency equations numerically with `brentq`, and check the mode shapes against the differential equation, their boundary conditions and orthogonality. 3. Implement a finite-difference eigenproblem (pinned–pinned) and a finite-element eigenproblem with Hermite cubic elements and consistent mass, and show that they converge to the closed form at second and fourth order respectively. 4. Use the finite-element model where no closed form exists (a tank part-way along the span, tapered stiffness and mass), and judge the single-degree-of-freedom estimate of `P06.13` against the continuous answer. 5. Decide where to place a mass on a wing to move its first two bending frequencies, and say which modes a model of n elements can be trusted for. |
| **Prerequisites** | `P05.12` (beam deflection, EI w⁗ = q, boundary conditions), `P02.11` (separation of variables, the wave equation), `P06.13` (SDOF model of a cantilever, Rayleigh's effective mass), `P06.14` (the generalized eigenproblem K φ = ω² M φ, orthogonality, Rayleigh quotient), `P03.13` (finite elements: weak form, assembly, the axial bar's eigenproblem), `P03.11` (finite differences, observed order), `P03.02` (root finding, `brentq`). Skill check: you can write the general solution of y⁗ = β⁴y and apply four boundary conditions to a 4 × 4 determinant. |
| **Engineering question** | The light-aeroplane wing of `P05.12` (EI = 6.958 × 10⁵ N·m² at the root, a 5 m semi-span) carries 45 kg of wing structure spread along the span (illustrative). What are its first bending frequencies, how good is a one-degree-of-freedom estimate, how many finite elements does a model need to give them to 10⁻⁵, and how do a 40 kg tank at 2 m, the same mass at the tip, and a tapered spar change them? The answers: 6.224, 39.00 and 109.2 Hz (spar-alone model); the SDOF estimate of `P06.13` is 1.5 % high; 8 cubic elements give mode 1 to 2 × 10⁻⁶ and 16 give mode 3 to 4 × 10⁻⁵; the tank at 2 m lowers f₁ to 5.69 Hz and f₂ to 24.8 Hz, at the tip it lowers f₁ to 2.89 Hz; tapering the stiffness to 30 % and the mass 3 : 1 towards the tip raises f₁ to 6.91 Hz. |
| **Assumptions and limits** | Euler–Bernoulli beam: plane sections stay plane, shear deformation and rotary inertia are neglected (good when the wavelength of the mode is many times the beam's depth); linear elastic, small deflections; a clamped root that allows no rotation at all; mass per unit length ρA and stiffness EI that vary only slowly; bending in one plane, uncoupled from torsion (a real wing's bending and torsion modes couple through the offset of its centres of mass and stiffness); point masses rigidly attached and without rotary inertia; no aerodynamic forces (in still air). |
| **Analytical work** | The beam equation from a slice's equilibrium with d'Alembert's inertia force; separation of variables; the general spatial solution in cosh, sinh, cos and sin; the four boundary-condition sets and the frequency equations 1 + cos βL cosh βL = 0, sin βL = 0, cos βL cosh βL = 1; the tip-mass equation; asymptotes; mode shapes; orthogonality from the self-adjointness of d⁴/dx⁴; the dispersion relation ω = β²√(EI/ρA); the Hermite element stiffness and consistent mass matrices from the weak form; the exact eigenvalues of the finite-difference stencil. |
| **Computational deliverable** | `starter/main.py`: `frequency_equation`, `beta_L`, `natural_frequencies`, `mode_shape`, `fd_pinned_frequencies`, `beam_element`, `fem_matrices`, `free_dofs`, `fem_frequencies`, `rayleigh_tip_mass_estimate` (`wing_spar`, `fd_pinned_exact`, `observed_order` and `main` are given). Tests prove: the roots solve their equations, 1.875104 for the cantilever, nπ for the pinned beam, and approach their asymptotes; frequencies scale as √EI/L²; a tip mass lowers f₁, the SDOF estimate is an upper bound and becomes exact for a heavy tip mass; the mode shapes satisfy w⁗ = β⁴w, their boundary conditions and orthogonality; the FD eigenvalues equal their discrete closed form and converge at order 2.0 ± 0.1, from below; the element matrices are symmetric with rigid-body null vectors and the right total mass; the FE frequencies are upper bounds converging at order 4 ± 0.5 for three supports, the FE mode shapes match the closed form, the tip-mass FE model matches its frequency equation to 10⁻⁵; the FE stiffness reproduces the static cantilever of `P05.12` exactly at the nodes; added mass lowers every frequency; variable properties converge. |
| **Practical activity** | **computational experiment**: convergence of consistent-mass and lumped-mass finite elements and of finite differences against the frequency equations; a tip-mass sweep against the SDOF estimate; the effect of a 40 kg tank's spanwise position on the first three frequencies; and the effect of stiffness and mass taper. |
| **Equipment** | None. A computer with the environment from `P00.01`. (The `P05.18` cantilever and the `p05-cantilever` kit can be plucked and its first frequency timed with a phone, as an optional check; nothing in the lesson depends on it.) |
| **Safety** | None beyond normal computer use. If you pluck a clamped ruler or strip, keep your face clear of its tip. |
| **Evidence** | Completed `main.py`, passing tests, `beam_modes.png`, the derivations of the beam equation, the cantilever frequency equation and the element matrices, the experiment tables with your predictions, and a half-page note answering the engineering question. |
| **Assessment** | `assessments/phase-06/grading-rubric.md`; competencies 1–4. Pass: the cantilever frequency equation is derived (not quoted); all public tests pass; the observed orders of convergence of the FE and FD eigenproblems are reported with their reasons; the note gives the wing's first frequencies, the effect of the tank's position, and the assumption that most limits them. |
| **References** | `@mit-ocw-2-003sc-engineering-dynamics` Lecture 27: Vibration of Continuous Structures: Strings, Beams, Rods, etc. (required-free): the wave equation for strings and rods and the vibration of beams, with dispersion. `@libretexts-cline-variational-principles` 14.10 Discrete Lattice Chain; 16.2 The Continuous Uniform Linear Chain (required-free): the passage from a discrete chain to a continuum. `@ocw-18-085` "Lecture 17: Finite Elements in 1D (part 1)", "Lecture 18: Finite Elements in 1D (part 2)", "Lecture 19: Quadratic/Cubic Elements", "Lecture 20: Element Matrices; 4th Order Bending Equations" (required-free): cubic elements and the bending element matrices. `@ocw-16-90` Unit 2 "2.9 Introduction to Finite Elements" (supplementary-free). `@mit-ocw-3-11-mechanics-of-materials` Roylance module "Beam Displacements": Multiple integration (supplementary-free): the static equation behind the dynamic one. `@scipy-userguide` "Linear Algebra (scipy.linalg)" (required-free). |
| **Extensions** | ★ Add rotary inertia and shear deformation (a Timoshenko element) and find the mode number at which the wing spar's Euler–Bernoulli frequency is 5 % too high. ★ Give the tank a rotary inertia about the beam axis (a 4 × 4 point-mass matrix with an entry on w′) and show which mode it affects. ★ Add a torsion degree of freedom with an offset between the elastic axis and the centre of mass, and watch the first bending and torsion modes couple (the start of flutter analysis). |
| **Time** | unmeasured |

## 1. The engineering problem

A wing is a beam, and a beam has bending modes. The first of them matters for almost everything that moves the wing: gusts and turbulence excite it, the pilot's control inputs and the landing impact excite it, an engine or propeller excites it at its rotation frequencies, and in flight the aerodynamic forces that depend on the wing's motion couple it to the torsion modes, which is how flutter starts. Certification requires the wing's frequencies to be known, and every external store or tank added to a wing requires them to be known again, because a mass hung on a wing moves its modes. The engineer's tools are a closed form for uniform beams, which says how frequencies scale and gives an exact benchmark, and a finite-element model for the real, tapered wing with masses on it.

`P06.13` treated a cantilever with a tip mass as one mass on one spring, with an effective mass for the beam itself; `P06.14` treated a structure as a few lumped masses. Both are approximations to a **continuous system**, one with mass and stiffness everywhere along it and, therefore, infinitely many modes. This lesson writes the beam's partial differential equation, solves it exactly for uniform beams with the standard supports, and then builds two discretisations, finite differences and finite elements, and measures how fast they converge to the exact frequencies. That convergence is the lesson's independent check, and it is what tells you how many elements a real model needs.

The wing (illustrative): the spar of `P05.12`, EI = 6.958 × 10⁵ N·m², on a 5 m semi-span, carrying 45 kg of wing structure (9 kg/m). The spar-alone stiffness makes this a soft wing; the real wing's skins would raise its frequencies. Success criteria, written first: the roots of the frequency equations solve them to 10⁻⁹ relative; the mode shapes satisfy w⁗ = β⁴w and their boundary conditions; the FD frequencies equal their own discrete closed form and converge at order 2; the FE frequencies are upper bounds converging at order 4, and match the tip-mass frequency equation to 10⁻⁵.

## 2. Understand and derive

### 2.1 The beam equation

`P05.12` derived the static beam equation from the equilibrium of a slice, V′ = q, M′ = V, with EI w″ = M: (EI w″)″ = q. In vibration the load on the slice includes d'Alembert's inertia force −ρA ẅ per unit length (ρA is the mass per unit length, kg/m). With no external load:

**(EI w″)″ + ρA ẅ = 0,  and for a uniform beam  EI w⁗ + ρA ẅ = 0.**   (1)

Units: N·m² × m⁻⁴ = N/m and kg/m × m/s² = N/m. The boundary conditions are those of `P05.12`: clamped w = w′ = 0; pinned w = w″ = 0; free w″ = w‴ = 0.

**Separation of variables** (`P02.11`). Try w(x, t) = φ(x) cos ωt. Then EI φ⁗ = ρA ω² φ, i.e.

**φ⁗ = β⁴ φ,  β⁴ = ρA ω²/EI,**   (2)

whose general solution is φ = C₁ cosh βx + C₂ sinh βx + C₃ cos βx + C₄ sin βx. Four boundary conditions give four homogeneous equations in C₁…C₄, and a non-trivial solution exists only when their determinant vanishes: the **frequency equation**, a transcendental equation for βL. Each root βᵢL gives

**ωᵢ = (βᵢL)² √(EI/(ρA L⁴)).**   (3)

The frequency scales as √(EI/ρA)/L²: doubling the span divides every frequency by four.

### 2.2 The frequency equations

**Cantilever** (clamped at x = 0, free at x = L). φ(0) = 0 gives C₃ = −C₁; φ′(0) = 0 gives C₄ = −C₂. Then φ = C₁(cosh − cos) + C₂(sinh − sin), and the two free-end conditions φ″(L) = φ‴(L) = 0 give

C₁(cosh βL + cos βL) + C₂(sinh βL + sin βL) = 0,
C₁(sinh βL − sin βL) + C₂(cosh βL + cos βL) = 0.

The determinant, (cosh + cos)² − (sinh² − sin²), simplifies with cosh² − sinh² = 1 and cos² + sin² = 1 to 2 + 2 cos βL cosh βL:

**1 + cos βL cosh βL = 0:  βL = 1.875104, 4.694091, 7.854757, 10.995541, …**   (4)

with the mode shape **φ = cosh βx − cos βx − σ(sinh βx − sin βx), σ = (cosh βL + cos βL)/(sinh βL + sin βL)**. Because cosh βL grows exponentially, the roots approach those of cos βL = 0, βL → (2n − 1)π/2, quickly: the fifth root is 9π/2 to 10⁻⁶ (a test). Dividing (4) by cosh βL gives cos βL + 1/cosh βL = 0, a function of order one that `brentq` handles without overflow.

**Pinned–pinned.** φ(0) = φ″(0) = 0 removes the cosh and cos terms; φ(L) = φ″(L) = 0 then requires sinh βL · C₂ = 0 and sin βL · C₄ = 0, so **sin βL = 0, βL = nπ, φ = sin(nπx/L)**: the exact sine modes.

**Clamped–clamped.** The same steps give **cos βL cosh βL = 1: βL = 4.730041, 7.853205, 10.995608, 14.137165**, approaching (2n + 1)π/2.

**Cantilever with a tip mass.** A mass m_t at the free end changes the shear condition: the beam must accelerate the mass, EI φ‴(L) = −m_t ω² φ(L). With μ = m_t/(ρAL), the frequency equation becomes

**1 + cos βL cosh βL + μβL (cos βL sinh βL − sin βL cosh βL) = 0.**   (5)

For μ → ∞ the first root falls to zero as the mass dominates, and ω₁ → √(3EI/(L³m_t)): the tip mass on the beam's static stiffness, the model of `P06.13` (a test).

For the wing, √(EI/(ρAL⁴)) = 11.122 s⁻¹, so (3) and (4) give **f₁ = 6.224 Hz, f₂ = 39.00 Hz, f₃ = 109.2 Hz**: ratios 1 : 6.27 : 17.5, as (βᵢL)² dictates. The SDOF estimate of `P06.13` (stiffness 3EI/L³ on 33/140 of the beam's mass) gives 6.315 Hz, **1.5 % high**: Rayleigh's method always overestimates, because the static deflection shape is not quite the mode shape.

### 2.3 Orthogonality and dispersion

Multiply φᵢ⁗ = βᵢ⁴φᵢ by φⱼ, integrate over the span twice by parts, and use the boundary conditions (each makes the boundary terms vanish for clamped, pinned or free ends): ∫φⱼφᵢ⁗ dx = ∫φⱼ″φᵢ″ dx, which is symmetric in i and j. Subtracting the same with i and j swapped:

**(βᵢ⁴ − βⱼ⁴) ∫₀ᴸ φᵢ φⱼ dx = 0:  modes are orthogonal**   (6)

(with weight ρA, which is constant here; with a tip mass the weight includes the point mass). This is the continuous version of `P06.14`'s φᵢᵀMφⱼ = 0, and it is why a beam's response can be expanded in its modes exactly as a lumped structure's can.

**Dispersion.** A bending wave w = cos(βx − ωt) on an infinite beam satisfies (1) when ω = β²√(EI/ρA): its phase speed ω/β = √ω (EI/ρA)^{1/4} grows with frequency (2.003SC Lecture 27). Unlike the string of `P04.06`, whose waves all travel at one speed, a beam spreads a pulse out. At the wing's first frequency the phase speed is 104 m/s. Dispersion is also why a beam's frequencies are not harmonics of the first one.

### 2.4 Finite differences: the pinned beam

On n + 1 equal intervals of length h, the five-point fourth difference (wᵢ₋₂ − 4wᵢ₋₁ + 6wᵢ − 4wᵢ₊₁ + wᵢ₊₂)/h⁴ approximates w⁗ with error (h²/6)w⁽⁶⁾. At a pinned end w₀ = 0, and w″ = 0 is enforced with a ghost node, w₋₁ = −w₁. The discrete problem D w = λ w is then exactly diagonalised by the sine vectors sin(jπxᵢ/L), with eigenvalues

**λⱼ = (2 − 2 cos(jπh/L))²/h⁴,  so  ω_h = √(EI/ρA) (2 − 2 cos(jπh/L))/h² ≈ ωⱼ (1 − (jπh/L)²/12),**   (7)

from below and at second order (a test checks both the eigenvalues against (7) and the order). The error grows as j²: a grid that resolves mode 1 to 0.2 % resolves mode 3 only to 1.8 % (part 4).

### 2.5 Finite elements: Hermite cubic beam elements

`P03.13` built linear elements for the axial bar, whose weak form has first derivatives. The beam's weak form, ∫EI w″v″ dx = ω²∫ρA w v dx for all admissible v, has second derivatives, so the element shape functions must have continuous slopes: **Hermite cubics**, with the deflection and the slope at each end as the four degrees of freedom (18.085 Lecture 20). Integrating the products of their second derivatives and of the functions themselves over an element of length h gives

**k = (EI/h³) [[12, 6h, −12, 6h], [6h, 4h², −6h, 2h²], [−12, −6h, 12, −6h], [6h, 2h², −6h, 4h²]],**
**m = (ρAh/420) [[156, 22h, 54, −13h], [22h, 4h², 13h, −3h²], [54, 13h, 156, −22h], [−13h, −3h², −22h, 4h²]].**   (8)

k has the rigid translation (1, 0, 1, 0) and rotation (0, 1, h, 1) as null vectors, and m gives the element's total mass ρAh for a translation (tests). m is the **consistent** mass matrix, derived from the same shape functions as k. Assembled and restricted to the free degrees of freedom, they form the generalized eigenproblem of `P06.14`, K φ = ω² M φ. Two theorems predict the results: because the FE solution minimises the Rayleigh quotient over a subspace of the admissible functions, **every FE frequency is an upper bound** on the exact one; and for cubic elements the eigenvalue error falls as h^{2(p+1−m)} = h⁴ (p = 3 the polynomial degree, m = 2 the derivative order of the energy), so **the frequency error is fourth order**. A static load test closes the loop with `P05.12`: under a uniform load with consistent nodal loads, cubic elements give the exact cantilever deflection at the nodes.

### 2.6 Where the model stops being valid

- **Shear and rotary inertia.** Euler–Bernoulli treats the beam as infinitely stiff in shear. For mode n the error grows roughly as (n × depth/span)²; for the 0.18 m spar on 5 m it is small for the first few modes, but higher modes, and short deep beams, need Timoshenko theory (Extension 1).
- **The clamped root.** A real wing root rotates a little under moment (`P05.12` §2.5); any root flexibility lowers the frequencies, the first one most.
- **Bending alone.** A wing's bending and torsion couple through the chordwise offset of its mass and stiffness centres; the modes of a real wing are mixtures. This uncoupled model is the starting point of flutter analysis, not its result.
- **The spar alone.** Skins and stringers raise EI; fuel, systems and control surfaces add mass. The numbers here are for the illustrative spar-and-mass model.
- **Large meshes.** The FE eigen-solution's round-off grows with the stiffness matrix's condition number (∝ n⁴): at 32 elements the first frequency's error has reached the 10⁻⁸ floor and the observed order falls from 4 to 3.2 (part 3's printout).

## 3. Build the computational model

Open `starter/main.py`. `wing_spar`, `fd_pinned_exact`, `observed_order` and `main` are given.

### Step 1: the closed form

`frequency_equation(bL, bc, mu)` returns the left side of (4), sin βL, or cos βL cosh βL − 1, each divided by cosh βL, and (5) for a tip mass. `beta_L(bc, n, mu)` scans (0.01, (n + 2)π] for sign changes and refines each with `brentq` (xtol 10⁻¹⁴). `natural_frequencies` is (3). `mode_shape` evaluates the closed forms of 2.2 and scales to max |φ| = 1, positive at the free end.

### Step 2: finite differences

`fd_pinned_frequencies` assembles the pentadiagonal stencil on the interior nodes, subtracts 1 from the first and last diagonal entries (the ghost node w₋₁ = −w₁ adds −w₁ to row 1), and takes `eigvalsh`.

### Step 3: finite elements

`beam_element` is (8). `fem_matrices` assembles n_el elements, evaluates EI and ρA at element midpoints when they are functions, and adds point masses to the nearest node's deflection degree of freedom. `free_dofs` removes (w, w′) at a clamp and w at a pin. `fem_frequencies` calls `scipy.linalg.eigh(K, M, subset_by_index=[0, n − 1])` on the free degrees of freedom and returns the frequencies and the nodal deflections of each mode. `rayleigh_tip_mass_estimate` is `P06.13`'s SDOF model.

```python
for e in range(n_el):
    k, m = beam_element(EI_e, rhoA_e, h)
    K[2*e:2*e+4, 2*e:2*e+4] += k          # dofs (w_e, w'_e, w_e+1, w'_e+1)
    M[2*e:2*e+4, 2*e:2*e+4] += m
```

Running `main()` prints

```
clamped-free     beta L = [ 1.875104  4.694091  7.854757 10.995541]
pinned-pinned    beta L = [ 3.141593  6.283185  9.424778 12.566371]
clamped-clamped  beta L = [ 4.730041  7.853205 10.995608 14.137165]

wing spar (EI 6.958e+05 N m^2, 9.0 kg/m, 5.0 m): bending f = [  6.2237  39.0035 109.2111] Hz
   SDOF estimate (P06.13, no tip mass) 6.3149 Hz

FEM convergence, cantilever, mode 1 and 3 relative errors
   n_el   2: mode 1 4.834e-04   mode 3 2.182e-01
   n_el   4: mode 1 3.271e-05   mode 3 7.742e-03
   n_el   8: mode 1 2.083e-06   mode 3 6.082e-04
   n_el  16: mode 1 1.312e-07   mode 3 3.975e-05
   n_el  32: mode 1 1.480e-08   mode 3 2.512e-06
   observed orders mode 1 [3.89 3.97 3.99 3.15], mode 3 [3.67 3.94 3.98]
FD convergence, pinned-pinned, mode 1 relative error
     10 intervals: 8.198e-03
     20 intervals: 2.054e-03
     40 intervals: 5.139e-04
     80 intervals: 1.285e-04
    160 intervals: 3.212e-05
   observed orders [2. 2. 2. 2.]

tip mass 10 kg: closed form f = [ 4.5131 31.9284] Hz, FEM (40 el) [ 4.5131 31.9284] Hz; SDOF estimate f1 4.5306 Hz

tip mass 40 kg: closed form f = [ 2.8892 28.9238] Hz, FEM (40 el) [ 2.8892 28.9238] Hz; SDOF estimate f1 2.8911 Hz
40 kg tank at 2 m (FEM): f = [ 5.6918 24.8199 94.4728] Hz
tapered spar (EI to 30 %, mass per metre 3:1 root to tip): f = [ 6.9142 34.9256 91.5181] Hz (80 el: [ 6.9138 34.9237 91.5128])
```

and writes `beam_modes.png`: the first three cantilever modes in closed form with the 8-element FE nodal values on top; the log–log convergence of the FE (slope 4) and FD (slope 2) frequency errors; and the cantilever frequency equation with and without a 40 kg tip mass.

### Verify

Run `python -m pytest curriculum/phase-06/lesson-15-continuous-systems`. The tests:

1. **Roots**: (4) and the clamped–clamped equation satisfied to 10⁻⁹ relative; 1.875104; nπ; the asymptotes; an unknown support is refused.
2. **Scaling**: (3); four times the EI and twice the span halves every frequency.
3. **Tip mass**: no mass leaves (4); heavier masses lower f₁; the SDOF estimate is an upper bound, within 2 % with no tip mass and exact (10⁻³) as μ → ∞.
4. **Mode shapes**, three supports: w⁗ = β⁴w by fourth differences, the boundary conditions, and orthogonality to 10⁻⁶.
5. **Finite differences**: eigenvalues equal (7) to 10⁻¹⁰; order 2 ± 0.1; from below.
6. **Elements**: symmetric k and m; rigid-body null vectors; total mass; positive definite m; the assembled mass with a point mass.
7. **Finite elements**, three supports: upper bounds; order 4 ± 0.5 on 8–32 elements; 32 elements within 10⁻⁵; the 16-element mode shapes match the closed form to 2 × 10⁻³; the tip-mass model matches (5) to 10⁻⁵; the static cantilever under a uniform load equals `P05.12`'s qx²(6L² − 4Lx + x²)/(24EI) at the nodes to 10⁻⁹; a mid-span tank lowers every frequency, mode 2 most; tapered properties converge; doubling EI multiplies the frequencies by √2.
8. `main` writes the figure.

### Compare

The independent check of this lesson is the discretised eigenproblem against the closed-form frequency equation. Finite differences converge at order 2.00 for the pinned beam, exactly as (7) predicts (it is not an estimate: the discrete eigenvalues are known in closed form). Finite elements converge at order 3.9–4.0 for the cantilever until round-off: 8 elements give mode 1 to 2.1 × 10⁻⁶ and mode 3 to 6.1 × 10⁻⁴; 16 elements give mode 3 to 4.0 × 10⁻⁵. The FE model with a 40 kg tip mass reproduces the tip-mass frequency equation (5), 2.8892 and 28.9238 Hz, to all printed digits. The established tool, `scipy.linalg.eigh` for the generalized symmetric problem, is used inside both discretisations; what is being checked is the modelling, not the solver.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** All results are outputs of the model, not measurements of a wing. Write each prediction before running it.

1. **Consistent against lumped mass.** Repeat the FE convergence study with a lumped (diagonal) mass matrix: half an element's mass at each of its nodes on the deflection degrees of freedom, none on the slopes (the slopes are then condensed out statically). Predict the sign and order of its error. *Finding:* relative errors of modes 1 and 3:

   | elements | consistent, mode 1 | consistent, mode 3 | lumped, mode 1 | lumped, mode 3 |
   |---|---|---|---|---|
   | 4 | +3.3 × 10⁻⁵ | +7.7 × 10⁻³ | −2.8 × 10⁻² | −0.138 |
   | 8 | +2.1 × 10⁻⁶ | +6.1 × 10⁻⁴ | −7.1 × 10⁻³ | −0.039 |
   | 16 | +1.3 × 10⁻⁷ | +4.0 × 10⁻⁵ | −1.8 × 10⁻³ | −0.010 |
   | 32 | +1.0 × 10⁻⁸ | +2.5 × 10⁻⁶ | −4.5 × 10⁻⁴ | −0.0025 |

   Consistent mass: from above, order 4. Lumped mass: from below, order 2. The lumped model is cheaper (a diagonal M) and still useful, but at equal accuracy it needs far more elements.
2. **Tip-mass sweep.** μ = m_t/(ρAL) = 0, 0.1, 0.25, 0.5, 1, 2 and 10 (m_t = 0 to 450 kg). Compare f₁ from (5) with the SDOF estimate. *Finding:* f₁ = 6.2237, 5.2534, 4.3840, 3.5691, 2.7566, 2.0501, 0.9583 Hz; SDOF errors +1.465, +0.725, +0.346, +0.150, +0.053, +0.016, +0.001 %. f₂ falls from 39.00 Hz to 27.46 Hz at μ = 10, towards the clamped–pinned beam's first frequency, (3.9266)² × 11.122/(2π) = 27.29 Hz: a heavy tip mass acts as a pin for the second mode. The SDOF model is excellent once the tip mass exceeds a quarter of the beam's mass, which justifies `P06.13`'s camera bracket.
3. **Where to put the tank.** A 40 kg tank at x = 0, 1, 2, 3, 4, 5 m (FEM, 40 elements). Predict which positions move f₂ least. *Finding:*

   | tank at | f₁ (Hz) | f₂ (Hz) | f₃ (Hz) |
   |---|---|---|---|
   | 0 m | 6.224 | 39.00 | 109.2 |
   | 1 m | 6.178 | 33.02 | 73.81 |
   | 2 m | 5.692 | 24.82 | 94.47 |
   | 3 m | 4.676 | 29.69 | 96.87 |
   | 4 m | 3.671 | 38.88 | 98.25 |
   | 5 m | 2.889 | 28.92 | 90.28 |

   A point mass lowers a mode in proportion to the square of the mode's deflection at that point (the Rayleigh quotient, `P06.14`). The cantilever's second mode has its largest inboard deflection near 2 m (φ₂ = −0.68 there) and a node near 4 m (φ₂ = 0.07), so a tank at 2 m drops f₂ by 36 % and a tank at 4 m barely touches it while dropping f₁ by 41 %. Mass placement is a frequency-placement tool, which flutter engineers use to keep bending and torsion frequencies apart.
4. **Taper.** EI falling linearly to (1 − taper) of the root value, taper = 0, 0.3, 0.5, 0.7, with uniform mass and with mass per metre falling 3 : 1 to the tip (same 45 kg). *Finding:* uniform mass: f₁ = 6.224, 6.032, 5.890, 5.731 Hz; tapered mass: 7.469, 7.252, 7.092, 6.914 Hz. Removing stiffness outboard costs little frequency (8 % at 70 % taper) because the outboard wing bends little; moving mass inboard gains 20 %.
5. **Finite differences, mode by mode.** The pinned beam with 10, 20 and 40 intervals: errors of modes 1–3 against −(jπh/L)²/12 from (7). *Finding:* 10 intervals: −0.0082, −0.0325, −0.0719 (estimate −0.0082, −0.0329, −0.0740); 20: −0.0021, −0.0082, −0.0184; 40: −0.00051, −0.0021, −0.0046. Mode 3 needs four times mode 1's grid for the same accuracy: count points per wavelength, not points per beam.

Raw data format: CSV tables `n_el, consistent_err_1, consistent_err_3, lumped_err_1, lumped_err_3` (1); `mu, m_tip_kg, f1_Hz, f2_Hz, f1_sdof_Hz, sdof_err_pct` (2); `x_tank_m, f1_Hz, f2_Hz, f3_Hz` (3); `taper, f1_uniform_mass_Hz, f1_tapered_mass_Hz` (4); `intervals, err_mode1, err_mode2, err_mode3` (5); plus `beam_modes.png`.

## 5. Validate and improve

| quantity | closed form (§2) | discrete | difference | cause |
|---|---|---|---|---|
| cantilever f₁ (Hz) | (4): 6.22375 | FE 8 elements | +2.1 × 10⁻⁶ | discretisation, order 4 |
| cantilever f₃ (Hz) | (4): 109.211 | FE 16 elements | +4.0 × 10⁻⁵ | discretisation, order 4 |
| pinned f₁ | sin βL = 0 | FD 160 intervals | −3.2 × 10⁻⁵ | discretisation, order 2, (7) |
| f₁, f₂ with 40 kg tip mass (Hz) | (5): 2.8892, 28.9238 | FE 40 elements: 2.8892, 28.9238 | < 10⁻⁵ | — |
| f₁, SDOF estimate (Hz) | 6.2237 | 6.3149 (`P06.13`) | +1.5 % | Rayleigh upper bound |
| static tip deflection, uniform load | `P05.12` qL⁴/(8EI) | FE nodal (test) | < 10⁻⁹ | cubic elements are nodally exact |
| tapered spar f₁ (Hz) | none | 6.9138 (80 el), 6.9142 (160 el) | 6 × 10⁻⁵ | midpoint property sampling |

The required independent check, finite-difference and finite-element eigenproblems converging to the closed-form frequency equations, holds at the predicted orders for every support. The engineering answer: **the spar-and-mass wing has its first three bending modes at 6.22, 39.0 and 109.2 Hz; an 8-element cubic model gives the first to two parts in a million and 16 elements the third to four parts in 10⁵; a 40 kg tank at 2 m lowers f₁ to 5.69 Hz and f₂ to 24.8 Hz, at 4 m it lowers f₁ to 3.67 Hz but leaves f₂ at 38.9 Hz, and at the tip it lowers f₁ to 2.89 Hz.** The tapered model has no closed form; its 80- and 160-element answers differ by 6 × 10⁻⁵, which is its own check (sampling varying properties at element midpoints converges more slowly than the uniform case).

What to improve first: the physics, not the mesh. Coupling with torsion (Extension 3) and the skins' contribution to EI change the answers by tens of per cent; refining beyond 16 elements changes them by parts per million.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (1), (4), (5) (at least the modified shear condition), (6) and (8); `beam_modes.png`; the five experiment tables with your predictions; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) giving the wing's first frequencies, the SDOF estimate's error, the mesh needed, and the effect of the tank's position.

### Adversarial self-checks

- **What if the span doubles?** Every frequency falls by four (3); the wing's f₁ would be 1.56 Hz, close to the short-period and gust frequencies of the aircraft, where structure and flight dynamics couple.
- **What if the tip mass goes to infinity, or to zero?** f₁ → √(3EI/(L³m_t)) (a test) or the bare cantilever; f₂ tends to the clamped–pinned value. Both limits are checked.
- **Which assumption fails first in the real world?** The uncoupled bending model: a real wing's first bending mode carries some twist. Then the rigid root, which lowers f₁.
- **What measurement could falsify the model?** A ground vibration test (a shaker on the wing, accelerometers along the span): the measured first bending frequency and shape. A shape with a kink at the root means the root is flexible; a mode with twist means coupling; a frequency well above 6.2 Hz means the skins stiffen the wing.
- **Does the result survive a change of method?** FE, FD and the closed form agree at the orders theory predicts; consistent and lumped mass bracket the answer from above and below (experiment 1).
- **Which decision changes if the uncertainty doubles?** If EI is uncertain by ±20 % instead of ±10 %, f₁ spans 5.57–6.82 Hz instead of 5.90–6.53 Hz (frequency ∝ √EI); a requirement to keep f₁ above, say, 6 Hz would no longer be met with confidence, and a ground vibration test would become necessary rather than confirmatory.

## Check yourself

1. Write the Euler–Bernoulli vibration equation and the units of each term.
2. What are the boundary conditions at a clamped, a pinned and a free end, and which one changes when a mass is attached to a free end?
3. Write the cantilever frequency equation. Why do its roots approach (2n − 1)π/2?
4. How does a beam's natural frequency scale with its length, with EI and with its mass per unit length?
5. Why are the mode shapes orthogonal, and with what weight?
6. Why must beam finite elements have continuous slopes, and what are the four degrees of freedom of a Hermite element?
7. A consistent-mass FE model gives frequencies from above and a lumped-mass one from below. Why, and what are their orders of convergence?
8. Where on a cantilever would you put a mass to lower the second bending frequency without moving the first much, and why?
