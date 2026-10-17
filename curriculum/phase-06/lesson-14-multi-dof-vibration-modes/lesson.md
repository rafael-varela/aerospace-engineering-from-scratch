# P06.14  Multi-DOF Vibration & Modes

> A structure with n masses is n single-degree-of-freedom systems in disguise. Find the disguise (the modes) and every response is a sum of lessons-13 answers.

| | |
|---|---|
| **Learning outcomes** | 1. Assemble the mass and stiffness matrices of a lumped chain and solve the generalized eigenproblem K φ = ω² M φ by Cholesky reduction, with mass-normalised modes. 2. Prove and check M- and K-orthogonality, use the Rayleigh quotient as an upper bound, and compute modal participation factors and effective masses for base excitation. 3. Fit proportional (Rayleigh) damping to two modal damping ratios and test whether a damping matrix is diagonalised by the modes. 4. Compute free, forced and frequency responses by modal superposition and show that they equal direct integration of the coupled equations and the direct inverse (K − ω²M + iωC)⁻¹; quantify modal truncation. 5. Predict the accelerations and the interface force of a spacecraft stack in a base sine test, and decide how many modes, and which damping assumption, the prediction needs. |
| **Prerequisites** | `P06.13` (SDOF free and forced response, base excitation, transmissibility), `P01.16` (eigenvalues of symmetric matrices, the Rayleigh quotient), `P01.13` (linear systems), `P02.05` (the matrix exponential; systems of ODEs), `P06.11` (Lagrange's equations; normal modes of the double pendulum), `P03.06` (`solve_ivp`), `P03.08` (Cholesky factorisation). Skill check: you can find the eigenvalues of a 2 × 2 symmetric matrix by hand and explain why its eigenvectors are orthogonal. |
| **Engineering question** | A 1100 kg spacecraft (an adapter ring of 150 kg, a service module of 600 kg and a payload module of 350 kg on axial stiffnesses of 6, 4 and 2.5 × 10⁸ N/m: illustrative numbers) sits on the launch-vehicle interface. What are its axial modes and how much mass does each carry into the interface; in a 1 g axial sine test at the first resonance, with 2 % damping, how hard is the payload shaken and what force passes through the interface; and how many modes does the prediction need? The answers: 73.5, 180.0 and 420.2 Hz, with effective masses 1008, 45 and 47 kg (91.6 %, 4.1 %, 4.3 %); 31.0 g on the payload and 247.4 kN through the interface, 23 times the 10.8 kN of the rigid stack; one mode plus the residual mass gives the interface force to 3 × 10⁻⁵, but the answer scales as 1/ζ, so the damping, not the mode count, sets the uncertainty. |
| **Assumptions and limits** | Lumped masses joined by massless linear springs in one direction (axial); the interface is a prescribed, rigid base motion; viscous damping, proportional (Rayleigh) unless stated; small motions; steady sinusoidal input for the test prediction. A real spacecraft has lateral, torsional and local modes, distributed mass and joints that are not linear; the three-mass model gives only the first few axial modes, and only roughly. |
| **Analytical work** | M and K of a fixed–free chain; the generalized eigenproblem and its reduction by Cholesky; orthogonality from the symmetry of M and K; mass normalisation; modal coordinates and the uncoupled equations; modal initial conditions; the Rayleigh quotient and its bound; participation factors Γ = Φᵀ M r and effective masses Γᵢ², summing to the total mass; Rayleigh damping ζ(ω) = a/(2ω) + bω/2; modal superposition of the FRF; the apparent mass Σ m_eff,i T_i(ω) with a residual mass for truncated modes; the uniform chain's closed-form frequencies. |
| **Computational deliverable** | `starter/main.py`: `chain_matrices`, `modes`, `rayleigh_quotient`, `participation`, `rayleigh_coefficients`, `modal_damping`, `modal_free_response`, `direct_integration`, `modal_integration`, `frf_direct`, `frf_modal`, `apparent_mass`, `base_sine_accelerations` (`spacecraft_stack` and `main` are given). Tests prove: M and K symmetric, K positive definite, the static flexibility of the chain; the Cholesky modes equal `scipy.linalg.eigh(K, M)`, are mass-normalised, M- and K-orthogonal, and the first has no node; a two-mass closed form and the n-mass uniform chain's closed form with j − 1 nodes in mode j; the Rayleigh quotient is exact at a mode and bounded by ω₁² and ω_n²; effective masses sum to rᵀ M r; Rayleigh damping hits its targets and has zero modal coupling while a single damper does not; modal superposition equals direct integration for the free response (to 10⁻⁹) and a forced response (all modes, 10⁻⁸; truncation errors falling with each mode); the modal FRF equals the direct inverse at five frequencies, K⁻¹ at zero and is symmetric; the undamped modal response conserves energy; the apparent mass tends to the total mass and equals the direct inertia force; at low frequency every mass moves with the base. |
| **Practical activity** | **computational experiment**: modal truncation for accelerations and interface force across the test band; four damping models and the damping they assign to the third mode; non-proportional damping from a single discrete damper against its diagonal (modal) approximation and the exact complex modes; sensitivity of the modes to payload mass and to each stiffness; and the payload's peak response in a 5–200 Hz sweep at three damping levels. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `modal_analysis.png`, the derivations of orthogonality, the uncoupled modal equations, the effective-mass sum and the apparent mass, the experiment tables with your predictions, and a half-page note answering the engineering question. |
| **Assessment** | `assessments/phase-06/grading-rubric.md`; competencies 1–4. Pass: K φ = ω² M φ and orthogonality are derived (not quoted); all public tests pass; modal superposition is shown to agree with direct integration with numbers; the note gives the modes, effective masses, the test response and the interface force, and says which assumption most limits the prediction. |
| **References** | `@mit-ocw-2-003sc-engineering-dynamics` Lecture 22: Finding Natural Frequencies and Mode Shapes of a 2 DOF System; Lecture 23: Vibration by Mode Superposition; Lecture 24: Modal Analysis: Orthogonality, Mass Stiffness, Damping Matrix; Lecture 25: Modal Analysis: Response to IC's and to Harmonic Forces; Lecture 26: Response of 2-DOF Systems by the Use of Transfer Functions; Recitation 12 Notes: Modal Analysis (required-free). `@mit-ocw-16-07-dynamics` L19 "Vibration, normal modes, natural frequencies, instability": Vibration of Multi-Degree of Freedom Systems; Normal Modes of Multi-Degree of Freedom Systems; Characteristic Value Problem; Solution for Unequal Masses (required-free). `@libretexts-cline-variational-principles` Ch. 14 "Coupled Linear Oscillators": 14.3 Normal Modes; 14.6 General Analytic Theory for Coupled Linear Oscillators; 14.8 Three-body coupled linear oscillator systems; 14.11 Damped Coupled Linear Oscillators (required-free). `@mit-ocw-18-06-linear-algebra` Lecture 25: Symmetric matrices and positive definiteness (supplementary-free). `@mit-ocw-18-03-differential-equations` L30 Decoupling Linear Systems with Constant Coefficients (supplementary-free). `@scipy-userguide` "Linear Algebra (scipy.linalg)" (required-free). |
| **Extensions** | ★ Let the interface itself be a mass on a spring (a flexible launch-vehicle upper stage) and show how the coupled system's first mode falls below both the stack's and the vehicle's ("coupled loads"). ★ Add a tuned mass damper to the payload module, tuned to mode 1, and find the mass ratio that halves the payload's resonant acceleration. ★ Replace the axial chain by a lateral one with bending springs and rotational inertia, and compare its first mode with a cantilever beam of `P06.15`. |
| **Time** | unmeasured |

## 1. The engineering problem

A launch vehicle shakes its payload. At the low-frequency end of the spectrum, below about 100 Hz, the vehicle's own structural modes, engine thrust oscillations and transients at ignition and stage separation drive the spacecraft through its interface ring, and the spacecraft answers with its own modes. Two numbers dominate the conversation between a spacecraft team and a launch-vehicle team: the spacecraft's **fundamental frequencies** (axial and lateral), which the launcher's user's guide requires to stay above stated minimums so that the two structures do not couple, and the **interface forces** that the spacecraft's modes put back into the vehicle. Both come from a modal analysis: the eigenproblem of the spacecraft's mass and stiffness, and the participation of each mode in motion of the interface.

`P06.13` treated a single mass on a spring. A spacecraft is a stack of masses on springs, and its equations are coupled: pushing one mass moves the others. The central fact of linear vibration is that this coupling is only apparent. There are special shapes, the **normal modes**, in which every mass moves in proportion at a single frequency, and in the coordinates that measure "how much of each mode", the equations separate into n independent single-degree-of-freedom systems. `P06.11` found the two modes of a double pendulum; this lesson builds the general machinery, proves the orthogonality that makes it work, and uses it to predict a sine test.

The model (illustrative): three lumped masses on axial springs, fixed at the launch-vehicle interface. Success criteria, written first: the modes equal `scipy.linalg.eigh(K, M)` to 10⁻¹² and satisfy ΦᵀMΦ = I; the effective masses sum to the total mass to round-off; modal superposition with all modes equals direct integration to 10⁻⁸ of the response for both a free and a forced motion, and the modal FRF equals the direct inverse to 10⁻⁹.

## 2. Understand and derive

### 2.1 The model and its matrices

Number the masses from the interface up: the adapter ring m₁, the service module m₂, the payload module m₃. Spring k₁ joins the interface to m₁, k₂ joins m₁ to m₂, k₃ joins m₂ to m₃. With xᵢ the displacement of mass i, Newton's law for each mass (or Lagrange's equations with T = ½Σmᵢẋᵢ² and V = ½Σkᵢ(xᵢ − xᵢ₋₁)², x₀ = 0) gives

**M ẍ + K x = f(t),  M = diag(m₁, m₂, m₃),  K = [[k₁ + k₂, −k₂, 0], [−k₂, k₂ + k₃, −k₃], [0, −k₃, k₃]].**   (1)

K is symmetric (Maxwell reciprocity, `P05.12`) and positive definite because the chain is fixed at the interface: every non-zero x stretches some spring. Its inverse is the flexibility matrix; a unit force on the payload stretches every spring by 1/kᵢ, so the payload moves Σ1/kᵢ (a test).

```
     payload  m3 = 350 kg   ── x3
        |  k3 = 2.5e8 N/m
     service  m2 = 600 kg   ── x2
        |  k2 = 4.0e8 N/m
     adapter  m1 = 150 kg   ── x1
        |  k1 = 6.0e8 N/m
  ======= LV interface (prescribed motion y)
```

### 2.2 The generalized eigenproblem

Look for motions in which every mass moves in proportion, x = φ cos ωt. Substituting into (1) with f = 0:

**K φ = ω² M φ.**   (2)

This is a **generalized** eigenproblem because of M. It has non-trivial solutions only when det(K − ω²M) = 0, a cubic in ω² for three masses (16.07 L19 works the characteristic value problem; 2.003SC Lecture 22 the two-mass case). The two-mass version, a test, is m₁m₂ω⁴ − (m₁k₂ + m₂(k₁ + k₂))ω² + k₁k₂ = 0.

**Reduction to a standard symmetric problem.** M is symmetric positive definite, so it has a Cholesky factor M = LLᵀ (`P03.08`). With φ = L⁻ᵀu, (2) becomes

L⁻¹ K L⁻ᵀ u = ω² u,   (3)

a standard eigenproblem with a symmetric matrix A = L⁻¹KL⁻ᵀ. Its eigenvalues ω² are real and positive (A is positive definite) and its eigenvectors u can be chosen orthonormal (`P01.16`, 18.06 Lecture 25). This is what `modes` does, and what `scipy.linalg.eigh(K, M)` does internally.

**Orthogonality.** Take two modes, Kφᵢ = ωᵢ²Mφᵢ and Kφⱼ = ωⱼ²Mφⱼ. Multiply the first by φⱼᵀ, the second by φᵢᵀ, and subtract; because K and M are symmetric, φⱼᵀKφᵢ = φᵢᵀKφⱼ and the same for M, so

(ωᵢ² − ωⱼ²) φᵢᵀ M φⱼ = 0  ⇒  **φᵢᵀ M φⱼ = 0 and φᵢᵀ K φⱼ = 0 for ωᵢ ≠ ωⱼ.**   (4)

Modes are not orthogonal in the ordinary sense, but in the "mass-weighted" one. Scaling each mode so that φᵢᵀMφᵢ = 1 (**mass normalisation**) gives, with Φ the matrix of modes as columns,

**ΦᵀMΦ = I,  ΦᵀKΦ = diag(ωᵢ²).**   (5)

**The stack.** For the illustrative numbers:

| mode | f (Hz) | shape (adapter, service, payload), mass-normalised |
|---|---|---|
| 1 | 73.542 | (+0.01130, +0.02734, +0.03900) |
| 2 | 180.042 | (−0.01432, −0.02892, +0.03654) |
| 3 | 420.213 | (+0.07959, −0.00908, +0.00104) |

Mode 1 moves everything one way, more at the top: the whole stack bouncing on its springs. Mode 2 has the payload against the rest. Mode 3 is the light adapter ring rattling between its two stiff springs while the rest stays almost still. For a uniform chain of n equal masses and springs the frequencies have a closed form, ω_j = 2√(k/m) sin((2j − 1)π/(2(2n + 1))), and mode j has j − 1 nodes (a test for n = 8).

**The Rayleigh quotient.** For any trial vector x, R(x) = xᵀKx/xᵀMx lies between ω₁² and ω_n², and equals ωᵢ² when x is mode i. Because R is stationary at a mode, a reasonable guess gives a frequency with a second-order error: the static deflection of the stack under its own weight gives f = 73.821 Hz, 0.4 % above the true 73.542 Hz. (`P06.13`'s effective mass was the same idea for a beam.)

### 2.3 Modal coordinates: n uncoupled systems

Write any motion as a combination of modes, x = Φη. Substituting into M ẍ + C ẋ + K x = f and multiplying by Φᵀ, (5) gives

η̈ + ΦᵀCΦ η̇ + diag(ωᵢ²) η = Φᵀ f.   (6)

If **ΦᵀCΦ is diagonal**, with entries 2ζᵢωᵢ, the equations separate:

**η̈ᵢ + 2ζᵢωᵢ η̇ᵢ + ωᵢ² ηᵢ = φᵢᵀ f(t),  i = 1 … n,**   (7)

n independent single-degree-of-freedom oscillators, each with its own frequency and damping, each forced by the projection of the load on its mode. The initial conditions follow from (5): **η(0) = ΦᵀM x(0), η̇(0) = ΦᵀM ẋ(0)** (2.003SC Lecture 25). Solve each (7) by `P06.13` and superpose: that is **modal superposition** (2.003SC Lecture 23).

**When is ΦᵀCΦ diagonal?** Always if C = aM + bK (**Rayleigh**, or proportional, damping), since then ΦᵀCΦ = aI + b diag(ωᵢ²) and

**ζᵢ = a/(2ωᵢ) + bωᵢ/2.**   (8)

Two targets fix a and b. For the stack with ζ = 0.02 in modes 1 and 2: a = 13.1229 s⁻¹, b = 2.5105 × 10⁻⁵ s, and (8) then assigns ζ₃ = 0.03563 to mode 3 whether you like it or not. A damper between two particular masses is not proportional: its modal matrix has off-diagonal terms that couple the modes (`modal_damping` measures them), and the exact solution has complex modes (Cline 14.11). Part 4 measures how much the diagonal approximation then costs.

### 2.4 Frequency response by modes

For a harmonic force f e^{iωt}, each (7) has the SDOF complex gain, and superposing,

**H(ω) = (K − ω²M + iωC)⁻¹ = Σᵢ φᵢφᵢᵀ/(ωᵢ² − ω² + 2iζᵢωᵢω),**   (9)

exact for proportional damping with all modes kept (a test at five frequencies). Each term is a rank-one matrix: near ωᵢ one term dominates and H is nearly φᵢφᵢᵀ/(2iζᵢωᵢ²), which is how a vibration test measures mode shapes. H is symmetric: the response at j to a force at i equals the response at i to a force at j.

### 2.5 Base excitation: participation and effective mass

In the sine test the interface moves, y(t). With relative coordinates z = x − r y, where **r** = (1, 1, 1) is the rigid-body motion of every mass with the interface (springs carry no force when z = 0), (1) becomes M z̈ + C ż + K z = −M r ÿ, the matrix version of `P06.13`'s (6). In modal coordinates z = Φη:

**η̈ᵢ + 2ζᵢωᵢη̇ᵢ + ωᵢ²ηᵢ = −Γᵢ ÿ,  Γᵢ = φᵢᵀ M r,**   (10)

the **modal participation factor**: how strongly a rigid motion of the base drives mode i. Its square is the **effective mass** m_eff,i = Γᵢ² (in kg for mass-normalised modes), and expanding r in modes, r = ΦΦᵀMr = ΦΓ, gives

**Σᵢ m_eff,i = rᵀ M r = total mass.**   (11)

For the stack: Γ = 31.749, −6.713, 6.850 kg^½, effective masses **1008.0, 45.1 and 46.9 kg (91.6 %, 4.1 %, 4.3 %)**, summing to 1100.000 kg. A launch-vehicle analyst reads this as: 92 % of the spacecraft's mass "lives" in the 73.5 Hz mode; the higher modes carry little into the interface.

**The interface force.** The force the interface must supply is the total inertia force, Σmᵢẍᵢ. Mode by mode, each modal mass behaves like an SDOF system on the base, so the interface force per unit base acceleration, the **apparent mass**, is

**F/ÿ = Σᵢ m_eff,i (1 + 2iζᵢrᵢ)/(1 − rᵢ² + 2iζᵢrᵢ),  rᵢ = ω/ωᵢ,**   (12)

each term an effective mass times its absolute transmissibility (`P06.13` (8)). At ω → 0 it is the total mass (the stack moves rigidly); at a resonance one term is multiplied by about 1/(2ζᵢ). If modes above some cut-off are dropped, their effective masses still move with the interface at frequencies well below their own, so they are added back as **residual mass**; that keeps (11) and the low-frequency limit exact.

**The 1 g sine test at f₁.** With ζ = 0.02, the direct solve gives absolute accelerations of 9.03 g (adapter), 21.73 g (service module) and **30.97 g (payload)**, and an **interface force of 247.4 kN**, 23 times the rigid-body 10.8 kN. In mode-1 terms: the payload's response is about Γ₁φ₃₁/(2ζ) = 1.238/0.04 = 30.95 g per g, and the interface force about m_eff,1/(2ζ) × g = 1008 × 25 × 9.81 N.

### 2.6 Where the model stops being valid

- **Lumped masses.** A real stack has distributed mass and many local modes (panels, propellant, instruments); a three-mass model gets the first one or two axial modes roughly and nothing above them.
- **One direction.** Lateral modes, usually much lower than axial (bending of the whole stack), and torsion are separate analyses, coupled in reality by off-axis mass.
- **Damping.** ζ is assumed, not derived; at test level it is typically 1–5 % per mode and the resonant response scales as 1/ζ (part 4).
- **A rigid, prescribed interface.** The launch vehicle is itself flexible. At resonance the spacecraft's 247 kN reaction would change the interface motion, which is why tests "force-limit" the input and why coupled-loads analyses combine both models (Extension 1).
- **Proportional damping.** Discrete dampers, isolators and joints give non-proportional damping; modal superposition with diagonal ζ is then an approximation (part 4, experiment 3).

## 3. Build the computational model

Open `starter/main.py`. `spacecraft_stack` and `main` are given.

### Step 1: matrices and modes

`chain_matrices(m, k)` builds (1) for any length. `modes(K, M)` is (3): `L = np.linalg.cholesky(M)`, A = L⁻¹KL⁻ᵀ (symmetrise it before `eigh` to remove round-off asymmetry), Φ = L⁻ᵀU, signs fixed so that each mode's largest entry is positive. `rayleigh_quotient` is one line.

```python
L = np.linalg.cholesky(M)
Li = np.linalg.inv(L)
lam, U = np.linalg.eigh(0.5 * (A + A.T))     # A = Li @ K @ Li.T
Phi = Li.T @ U                                # mass-normalised: Phi.T @ M @ Phi = I
```

### Step 2: participation and damping

`participation(Phi, M, r)` returns Γ and Γ². `rayleigh_coefficients` solves (8) at two frequencies; `modal_damping(Phi, C, w)` returns the diagonal of ΦᵀCΦ divided by 2ωᵢ and the largest normalised off-diagonal term.

### Step 3: the two routes

`modal_free_response` projects the initial conditions, evaluates each modal SDOF closed form and returns Φη. `modal_integration` integrates each (7) separately with `solve_ivp` (the first n_modes only) and superposes. `direct_integration` integrates the coupled (1) as a 2n-state first-order system. `frf_direct` inverts the dynamic stiffness; `frf_modal` is (9). `apparent_mass` is (12) with the residual mass of omitted modes; `base_sine_accelerations` solves the coupled base-driven problem directly and returns the absolute accelerations and the interface force Σmᵢẍᵢ.

Running `main()` prints

```
stack: masses [150. 600. 350.] kg, springs [6.0e+08 4.0e+08 2.5e+08] N/m, total 1100 kg
 mode   f (Hz)    shape (adapter, service, payload), mass-normalised      Gamma    m_eff (kg)   share
  1     73.542   [+0.01130 +0.02734 +0.03900]     31.7492     1008.01    91.6 %
  2    180.042   [-0.01432 -0.02892 +0.03654]     -6.7131       45.07     4.1 %
  3    420.213   [+0.07959 -0.00908 +0.00104]      6.8500       46.92     4.3 %
 sum of effective masses 1100.000000 kg; max|Phi^T M Phi - I| 4.4e-16; max off-diagonal |Phi^T K Phi| 2.8e-09
 Rayleigh quotient of the static 1 g deflection: f = 73.821 Hz (f1 = 73.542 Hz)

Rayleigh damping fitted to 2 % in modes 1 and 2: a = 13.1229 1/s, b = 2.5105e-05 s; modal zeta = [0.02    0.02    0.03563]; coupling 1.9e-16
free response from a 5 g static deflection (max 0.275 mm): modal vs direct max diff 1.2e-14 m
   sine force on the payload at f1, 1 mode(s): max |modal - direct| 5.72e-06 m (response peak 0.5126 mm)
   sine force on the payload at f1, 2 mode(s): max |modal - direct| 4.49e-08 m (response peak 0.5126 mm)
   sine force on the payload at f1, 3 mode(s): max |modal - direct| 1.38e-13 m (response peak 0.5126 mm)

1 g interface sine at f1 = 73.54 Hz: |acceleration| (g) adapter 9.03, service 21.73, payload 30.97; interface force 247.372 kN (rigid-body 10.79 kN)
   apparent mass with 1 mode(s) + residual: interface force 247.366 kN
   apparent mass with 2 mode(s) + residual: interface force 247.371 kN
   apparent mass with 3 mode(s) + residual: interface force 247.372 kN
```

and writes `modal_analysis.png`: the three mode shapes up the stack; the payload's driving-point receptance |H₃₃| by the direct inverse and by modes (all three, and mode 1 alone); and the interface apparent mass with three modes and with one mode plus residual mass.

### Verify

Run `python -m pytest curriculum/phase-06/lesson-14-multi-dof-vibration-modes`. The tests:

1. **Matrices**: symmetric M and K, K positive definite, the static flexibility Σ1/kᵢ.
2. **Modes**: equal to `eigh(K, M)` to 10⁻¹²; ascending; ΦᵀMΦ = I and ΦᵀKΦ = diag(ω²); KΦ = MΦΩ²; mode 1 without a node; the two-mass characteristic polynomial; the uniform 8-mass chain's closed form and node count.
3. **Rayleigh quotient**: exact at each mode, between ω₁² and ω₃² for 20 random vectors, within 1 % of ω₁² for the static deflection.
4. **Effective masses**: sum to 1100 kg (and to rᵀMr for another r), all non-negative, mode 1 above half.
5. **Damping**: Rayleigh targets met, ζ₃ by (8), coupling below 10⁻¹⁰; a single discrete damper has coupling above 0.05.
6. **Superposition**: free response from a 5 g static deflection and a velocity, modal against direct to 10⁻⁹ of the response; a sine force at f₁ with 1, 2 and 3 modes, the error falling with each mode and below 10⁻⁸ with all three; energy conserved by the undamped modal response.
7. **FRF**: modal (9) equals the direct inverse at five frequencies to 10⁻⁹; K⁻¹ at zero; symmetric.
8. **Base excitation**: the apparent mass tends to 1100 kg; it equals the direct inertia force at four frequencies to 10⁻⁹; one mode plus residual mass is within 1 % at f₁ and exact at zero frequency; at low frequency every mass moves with the base.
9. `main` writes the figure.

### Compare

The established tool is `scipy.linalg.eigh(K, M)`, which solves (2) with LAPACK's generalized symmetric driver; the lesson's Cholesky route reproduces its frequencies to 10⁻¹². The independent check required of this lesson is modal superposition against direct integration: for the free response from a 5 g static deflection (the stack's state at a sudden thrust cut-off) the two agree to 1.2 × 10⁻¹⁴ m on a 0.275 mm motion; for a sine force on the payload at f₁ they agree to 1.4 × 10⁻¹³ m with all three modes, 4.5 × 10⁻⁸ m with two and 5.7 × 10⁻⁶ m (1.1 % of the 0.513 mm peak) with mode 1 alone. The FRF route agrees with the direct inverse to 10⁻⁹, and the apparent mass (12) agrees with the direct inertia force, 247.372 kN, to six figures.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** All results are outputs of the model, not measurements of a spacecraft. Write each prediction before running it.

1. **Truncation across the band.** 1 g base sine at 73.5 Hz (f₁), 100 Hz, 180.0 Hz (f₂) and 420.2 Hz (f₃): payload acceleration and interface force, direct against modal with 1 and 2 modes (accelerations without residual correction; force with residual mass). Predict where one mode fails. *Finding:*

   | frequency | direct: payload (g), force (kN) | 1 mode | 2 modes |
   |---|---|---|---|
   | 73.5 Hz | 30.968, 247.372 | 30.971, 247.366 | 30.968, 247.371 |
   | 100 Hz | 1.802, 10.520 | 1.694, 10.742 | 1.803, 10.548 |
   | 180.0 Hz | 6.122, 11.327 | 0.486, 1.098 | 6.123, 11.332 |
   | 420.2 Hz | 0.105, 6.543 | 0.277, 0.595 | 0.023, 0.098 |

   A truncated model is good below its last retained mode and wrong at the resonance of any mode it drops. Keep every mode up to about 1.5 times the highest test frequency.
2. **Which damping model?** Four ways to give mode 1 a 2 % damping ratio, and what each implies for modes 2 and 3. *Finding:* Rayleigh fitted at modes 1 and 2: (0.02, 0.02, 0.0356); Rayleigh at 1 and 3: (0.02, 0.0143, 0.02); mass-proportional only: (0.02, 0.0082, 0.0035); stiffness-proportional only: (0.02, 0.049, 0.114). A choice made for convenience changes the predicted resonant response of mode 3 by a factor of 30. Fit the two coefficients at the two modes that matter most and state the others as consequences.
3. **Non-proportional damping.** A single dashpot c_d between the service and payload modules, c_d = 2 × 10⁴ and 10⁵ N·s/m. Compare the modal damping ratios from the diagonal of ΦᵀCΦ with the exact ones from the eigenvalues of the 6 × 6 state matrix, and the payload's receptance peaks with and without the off-diagonal terms. *Finding:* at 2 × 10⁴: diagonal ζ (0.00294, 0.03788, 0.00039) against exact (0.00294, 0.03789, 0.00039); peaks within 0.3 %. At 10⁵: diagonal (0.01471, 0.18941, 0.00194) against exact (0.01426, 0.1902, 0.00186); the diagonal approximation's peaks are 2.1 % low, 5.6 % high and 0.4 % high. The approximation is fine for light damping and degrades as the damper grows; the normalised coupling measure (1.0 here) overstates the effect because the modes are well separated.
4. **Sensitivity.** Add 50 kg to the payload; multiply k₁, k₂ or k₃ by 1.5. Predict which change raises f₁ most. *Finding:* baseline (73.542, 180.042, 420.213) Hz; +50 kg payload (70.865, 174.781, 420.204); k₃ × 1.5 (75.455, 214.697, 420.640); k₁ × 1.5 (78.304, 183.334, 474.680); k₂ × 1.5 (80.191, 182.463, 465.722). Stiffening the middle spring raises f₁ most: it carries the large service-module mass. The cheapest place to add stiffness is where the mode's strain energy is largest, which you can read from the mode shape.
5. **Damping and the test.** Sweep the 1 g base sine from 5 to 200 Hz with ζ = 0.01, 0.02 and 0.05 (Rayleigh at modes 1 and 2). *Finding:* payload peaks 61.91 g, 30.97 g and 12.42 g near 73.5 Hz; the one-mode estimate Γ₁φ₃₁/(2ζ) = 61.91, 30.95, 12.38 g. The damping assumption changes the answer by a factor of five; the mode count does not.

Raw data format: CSV tables `f_Hz, payload_direct_g, force_direct_kN, payload_1mode_g, force_1mode_kN, payload_2mode_g, force_2mode_kN` (1); `model, a, b, zeta1, zeta2, zeta3` (2); `cd, zeta_diag_1..3, zeta_exact_1..3, peak_ratio_1..3` (3); `case, f1_Hz, f2_Hz, f3_Hz, meff1_kg` (4); `zeta, payload_peak_g, f_peak_Hz, one_mode_estimate_g` (5); plus `modal_analysis.png`.

## 5. Validate and improve

| quantity | analysis (§2) | computation | difference | cause |
|---|---|---|---|---|
| frequencies (Hz) | eigh(K, M): 73.542, 180.042, 420.213 | Cholesky `modes` | < 10⁻¹² | — |
| ΦᵀMΦ − I | 0 (5) | 4.4 × 10⁻¹⁶ | — | round-off |
| sum of effective masses (kg) | 1100 (11) | 1100.000000 | < 10⁻¹² | — |
| f₁ from the static deflection | Rayleigh bound ≥ f₁ | 73.821 Hz vs 73.542 Hz | +0.4 % | trial vector not exact |
| free response, 5 g release | modal (7) | direct `solve_ivp` | 1.2 × 10⁻¹⁴ m | solver tolerance |
| sine force at f₁, 3 modes | modal | direct | 1.4 × 10⁻¹³ m | solver tolerance |
| same, 1 mode | modal | direct | 5.7 × 10⁻⁶ m (1.1 %) | truncation |
| interface force, 1 g at f₁ | apparent mass (12): 247.372 kN | direct inertia force 247.372 kN | < 10⁻⁹ | — |
| payload peak, 1 g at f₁ | Γ₁φ₃₁/(2ζ) = 30.95 g | direct 30.97 g | 0.07 % | other modes, √(1 − ζ²) |

The required independent check, modal superposition against direct integration, holds to the solver tolerance with all modes, and orthogonality holds to round-off. The engineering answer: **the stack's axial modes are at 73.5, 180.0 and 420.2 Hz; 92 % of its mass participates in the first; a 1 g sine at 73.5 Hz with 2 % damping shakes the payload at 31 g and puts 247 kN through the interface, 23 times the static weight.** One mode plus the residual mass predicts the interface force to 3 × 10⁻⁵ at f₁, so the mode count is not the issue; damping is. Between ζ = 0.01 and 0.05 the payload response spans 12–62 g. In a real test that is why the input is notched at the measured resonance and the interface force is limited during the run.

What to improve: a measured damping per mode (from a low-level sweep, `P06.13`), then a lateral model, then coupling with the launcher (Extension 1). What not to improve: the eigen-solver or the integration, already at round-off.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (2)–(5), (7), (8), (10)–(12); `modal_analysis.png`; the five experiment tables with your predictions; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) giving the modes, effective masses, the test accelerations and interface force, and the assumption that most limits the prediction.

### Adversarial self-checks

- **What if the payload mass doubles?** f₁ falls (to below 71 Hz for only +50 kg, experiment 4), the effective mass of mode 1 grows, and the interface force at resonance grows with it. The sum of effective masses must still equal the total mass: check it.
- **What if a spring becomes infinitely stiff?** Its two masses merge into one and one mode goes to infinite frequency; with k₁ → ∞ the adapter is locked to the interface and the stack becomes a two-mass chain. The code handles 10⁶ times stiffer springs with a warning-free Cholesky; try it.
- **Which assumption fails first in the real world?** The prescribed rigid interface: at 247 kN the launch vehicle's adapter would flex and the interface motion would not stay at 1 g. Then the assumed damping.
- **What measurement could falsify the model?** A low-level sine sweep with accelerometers at the three masses: the measured resonances, mode shapes (amplitude ratios at each peak) and the interface force (from a force-measuring ring) should match within the model's few-per-cent accuracy. A resonance the model does not have is a mode it leaves out (lateral, local).
- **Does the result survive a change of method?** Modal and direct routes agree to round-off; the FRF and the apparent mass agree with direct solves; the frequencies agree with `eigh`.
- **Which decision changes if the uncertainty doubles?** If ζ is known to ±50 % instead of ±25 %, the predicted payload response at f₁ spans 21–62 g instead of 25–41 g: a notch set by analysis would need a margin twice as large, and a low-level pre-test becomes mandatory rather than prudent.

## Check yourself

1. Why is K φ = ω² M φ called a generalized eigenproblem, and how does a Cholesky factor of M turn it into a standard symmetric one?
2. State the orthogonality relations of the modes and the property of M and K that guarantees them.
3. What does mass normalisation do, and what are ΦᵀMΦ and ΦᵀKΦ afterwards?
4. Write the uncoupled modal equations and the modal initial conditions.
5. What condition on C lets the modes uncouple the damped equations, and what damping does Rayleigh damping give a mode far above the two fitted ones?
6. Define the participation factor and the effective mass of a mode for base excitation. What do the effective masses sum to?
7. Why can a truncated modal model keep the omitted modes' effective mass as a rigid residual mass at low frequency?
8. At a resonance, roughly how large is the interface force compared with the rigid-body force, in terms of ζ and the effective mass?
