"""P08.18 Incompressible CFD from Scratch: starter code.

A two-dimensional lid-driven cavity solver in the stream function-vorticity formulation, written with
NumPy and SciPy only. The unit square cavity is closed on three sides; its lid (y = 1) slides to the
right at speed U = 1. In dimensionless form the steady incompressible Navier-Stokes equations become

    laplacian(psi) = -omega                                    (kinematics: defines psi from omega)
    u d(omega)/dx + v d(omega)/dy = (1/Re) laplacian(omega)     (vorticity transport)
    u = d(psi)/dy,  v = -d(psi)/dx                              (continuity holds identically)

with psi = 0 on every wall (no flow through the walls) and the wall vorticity from Thom's formula
(no slip). Second-order central differences on a uniform N x N grid of nodes, h = 1/(N - 1).

Two independent routes to the same discrete steady state:

* ``solve_cavity_newton``: Newton's method on the coupled sparse nonlinear system (P03.10), starting
  from the Stokes (Re -> 0) solution; a few iterations, each one sparse direct solve (P03.08);
* ``march_cavity``: explicit pseudo-time marching of the vorticity equation with a Poisson solve for
  psi at every step, run until the time derivative vanishes.

Verification and validation:

* the Poisson operator against a manufactured solution (second order), Thom's formula against a known
  stream function (first order at the wall), the convection operator exact on quadratics;
* the discrete divergence of (u, v) is zero to round-off (mass conservation);
* grid convergence with observed order on three grids, and Richardson extrapolation;
* the centreline velocities and the primary-vortex centre against Ghia, Ghia & Shin (1982), Tables
  I, II and V, Re = 100, typed in below and labelled as published data.

Array convention: a[i, j] is the value at (x_i, y_j), x_i = i h, y_j = j h; the lid is j = N - 1.

Run:  python main.py [out_dir]     (prints the lesson's numbers, writes cavity.png)
Test: pytest curriculum/phase-08/lesson-18-incompressible-cfd-lid-driven-cavity
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.interpolate import CubicSpline

# ---------------------------------------------------------------- published data (Ghia, Ghia & Shin 1982)
# PUBLISHED DATA, not computed here. U. Ghia, K. N. Ghia and C. T. Shin, "High-Re solutions for
# incompressible flow using the Navier-Stokes equations and a multigrid method", Journal of
# Computational Physics 48 (1982) 387-411. Their solution used a 129 x 129 grid at Re = 100.
# Table I: u-velocity along the vertical line through the geometric centre (x = 0.5), Re = 100 column.
GHIA_RE100_U = np.array([
    [1.0000, 1.00000], [0.9766, 0.84123], [0.9688, 0.78871], [0.9609, 0.73722], [0.9531, 0.68717],
    [0.8516, 0.23151], [0.7344, 0.00332], [0.6172, -0.13641], [0.5000, -0.20581], [0.4531, -0.21090],
    [0.2813, -0.15662], [0.1719, -0.10150], [0.1016, -0.06434], [0.0703, -0.04775], [0.0625, -0.04192],
    [0.0547, -0.03717], [0.0000, 0.00000]])          # columns: y, u
# Table II: v-velocity along the horizontal line through the geometric centre (y = 0.5), Re = 100 column.
GHIA_RE100_V = np.array([
    [1.0000, 0.00000], [0.9688, -0.05906], [0.9609, -0.07391], [0.9531, -0.08864], [0.9453, -0.10313],
    [0.9063, -0.16914], [0.8594, -0.22445], [0.8047, -0.24533], [0.5000, 0.05454], [0.2344, 0.17527],
    [0.2266, 0.17507], [0.1563, 0.16077], [0.0938, 0.12317], [0.0781, 0.10890], [0.0703, 0.10091],
    [0.0625, 0.09233], [0.0000, 0.00000]])           # columns: x, v
# Table V: primary vortex at Re = 100: stream function and location (x, y).
GHIA_RE100_PSI_MIN = -0.103423
GHIA_RE100_VORTEX_XY = (0.6172, 0.7344)


# ---------------------------------------------------------------- grid and operators (given)

def grid(N: int) -> tuple[np.ndarray, float]:
    """Node coordinates 0, h, ..., 1 (N nodes in each direction) and the spacing h = 1/(N - 1)."""
    x = np.linspace(0.0, 1.0, N)
    return x, 1.0 / (N - 1)


def observed_order(f_coarse: float, f_mid: float, f_fine: float, r: float = 2.0) -> tuple[float, float]:
    """Observed order p = ln((f1 - f2)/(f2 - f3))/ln r from three grids refined by the ratio r, and the
    Richardson extrapolate f3 + (f3 - f2)/(r^p - 1) (P03.11, P03.18)."""
    p = np.log(abs((f_coarse - f_mid) / (f_mid - f_fine))) / np.log(r)
    return float(p), float(f_fine + (f_fine - f_mid) / (r ** p - 1.0))


def _interior_index(N: int) -> np.ndarray:
    """Flat (C-order) indices of the interior nodes 1..N-2 in both directions."""
    i, j = np.meshgrid(np.arange(1, N - 1), np.arange(1, N - 1), indexing="ij")
    return (i * N + j).ravel()


# ---------------------------------------------------------------- the building blocks

def laplacian_matrix(N: int) -> sp.csr_matrix:
    """Five-point Laplacian on the (N-2)^2 interior nodes with zero Dirichlet values on the walls,
    as a sparse CSR matrix: (phi_{i-1,j} + phi_{i+1,j} + phi_{i,j-1} + phi_{i,j+1} - 4 phi_ij)/h^2.
    Unknowns ordered k = (i-1)(N-2) + (j-1). Symmetric negative definite."""
    raise NotImplementedError


def solve_poisson(f: np.ndarray, boundary: np.ndarray | None = None) -> np.ndarray:
    """Solve laplacian(phi) = f on the unit square with Dirichlet values: f is an N x N array (only its
    interior is used); ``boundary`` an N x N array whose edge values are imposed (zero if None).
    Returns the full N x N array phi, second-order accurate."""
    raise NotImplementedError


def velocities(psi: np.ndarray, U: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
    """u = d(psi)/dy and v = -d(psi)/dx by central differences at the interior nodes; on the walls
    the no-slip values (u = U on the lid y = 1, zero elsewhere; v = 0 everywhere)."""
    raise NotImplementedError


def wall_vorticity(psi: np.ndarray, U: float = 1.0) -> np.ndarray:
    """Thom's formula: an N x N array that is zero inside and holds the wall vorticity on the edges
    (corners left at zero; the five-point stencil never uses them). From the Taylor series of psi
    normal to a wall with psi = 0 and the wall velocity U_w: omega_w = -2 psi_adjacent/h^2 - 2 U_w/h
    on the lid (U_w = U) and -2 psi_adjacent/h^2 on the fixed walls. First-order accurate at the wall."""
    raise NotImplementedError


def divergence(u: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Discrete divergence du/dx + dv/dy by central differences at the interior nodes ((N-2)^2 array).
    For u, v from ``velocities`` it vanishes to round-off: the two difference operators commute."""
    raise NotImplementedError


def convection(psi: np.ndarray, omega: np.ndarray) -> np.ndarray:
    """The convection term u d(omega)/dx + v d(omega)/dy at the interior nodes, central differences,
    with u = d(psi)/dy and v = -d(psi)/dx also centred. Exact when psi and omega are quadratics."""
    raise NotImplementedError


def _lap_interior(a: np.ndarray, h: float) -> np.ndarray:
    return (a[2:, 1:-1] + a[:-2, 1:-1] + a[1:-1, 2:] + a[1:-1, :-2] - 4 * a[1:-1, 1:-1]) / h ** 2


def cavity_residual(psi: np.ndarray, omega: np.ndarray, Re: float, U: float = 1.0) -> tuple[float, float]:
    """Max-norm residuals of the two discrete steady equations at the interior nodes, with the wall
    vorticity recomputed from psi by Thom's formula: |lap(psi) + omega| and
    |convection - lap(omega)/Re|. Both are ~1e-10 or less at a converged solution."""
    raise NotImplementedError


def difference_operators(N: int) -> dict:
    """Sparse full-grid operators (N^2 x N^2, C-order flattening k = i N + j), given:
    Dx, Dy (central first differences), L (five-point Laplacian), all with zero rows on the walls;
    S (N^2 x n^2) places interior values into a full vector with zeros on the walls; P = S^T picks the
    interior; W (N^2 x n^2) and b (N^2) give Thom's wall vorticity, omega_wall = W psi_int + b."""
    _, h = grid(N)
    e = np.ones(N)
    D1 = sp.lil_matrix(sp.diags([-e[:-1], e[:-1]], [-1, 1]) / (2 * h))
    D2 = sp.lil_matrix(sp.diags([e[:-1], -2 * e, e[:-1]], [-1, 0, 1]) / h ** 2)
    for M in (D1, D2):
        M[0, :] = 0.0
        M[N - 1, :] = 0.0
    D1, D2, I = D1.tocsr(), D2.tocsr(), sp.identity(N, format="csr")
    Dx, Dy = sp.kron(D1, I).tocsr(), sp.kron(I, D1).tocsr()
    L = (sp.kron(D2, I) + sp.kron(I, D2)).tocsr()
    inner = _interior_index(N)
    n2 = inner.size
    S = sp.csr_matrix((np.ones(n2), (inner, np.arange(n2))), shape=(N * N, n2))
    col = {k: c for c, k in enumerate(inner)}
    rows, cols = [], []
    b = np.zeros(N * N)
    for m in range(1, N - 1):
        for wall, adj in (((m, 0), (m, 1)), ((m, N - 1), (m, N - 2)), ((0, m), (1, m)), ((N - 1, m), (N - 2, m))):
            rows.append(wall[0] * N + wall[1])
            cols.append(col[adj[0] * N + adj[1]])
        b[m * N + N - 1] = -2.0 / h                    # lid term, per unit lid speed
    W = sp.csr_matrix((np.full(len(rows), -2.0 / h ** 2), (rows, cols)), shape=(N * N, n2))
    return {"Dx": Dx, "Dy": Dy, "L": L, "S": S, "P": S.T.tocsr(), "W": W, "b": b, "h": h}


# ---------------------------------------------------------------- route 1: Newton on the steady equations

def solve_cavity_newton(N: int, Re: float, U: float = 1.0, tol: float = 1e-10,
                        max_iter: int = 30, initial: dict | None = None) -> dict:
    """Steady lid-driven cavity by Newton's method. Unknowns z = [psi_int, omega_int] (2 (N-2)^2).
    Residuals F1 = P L S psi + omega and F2 = P[(Dy psi_f)(Dx omega_f) - (Dx psi_f)(Dy omega_f)]
    - P L omega_f / Re, with psi_f = S psi and omega_f = S omega + W psi + U b. The Jacobian is
    assembled exactly from the sparse operators; the start is the Stokes solution (convection
    dropped, one linear solve), or the interior of ``initial`` (a previous solution on the same grid,
    for continuation in Re). Returns a dict with x, h, psi, omega (N x N, walls filled),
    iterations, residual history and wall-clock time."""
    raise NotImplementedError


# ---------------------------------------------------------------- route 2: pseudo-time marching

def march_cavity(N: int, Re: float, U: float = 1.0, dt: float | None = None, tol: float = 1e-9,
                 max_steps: int = 200_000) -> dict:
    """Explicit pseudo-time marching to the steady state: each step solves laplacian(psi) = -omega
    (one sparse LU factorisation, reused), sets the wall vorticity by Thom's formula and advances the
    interior vorticity by forward Euler, d(omega)/dt = -(u omega_x + v omega_y) + lap(omega)/Re.
    Stability of forward Euler with central differences needs dt <= h^2 Re/4 and dt <= 2/(Re U^2);
    the default takes 0.8 of the smaller. Stops when max|d(omega)/dt| < tol. Returns the same dict
    as ``solve_cavity_newton`` with ``steps`` in place of iterations."""
    raise NotImplementedError


# ---------------------------------------------------------------- post-processing

def centreline_profiles(sol: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """u(0.5, y) along the vertical centreline and v(x, 0.5) along the horizontal one (N odd, so both
    lines are grid lines). Returns y, u_c, x, v_c."""
    psi, x = sol["psi"], sol["x"]
    N = psi.shape[0]
    if N % 2 == 0:
        raise ValueError("use an odd N so that the centrelines are grid lines")
    u, v = velocities(psi, sol["U"])
    m = N // 2
    return x, u[m, :], x, v[:, m]


def compare_with_ghia(sol: dict) -> dict:
    """Interpolate the computed centreline profiles (cubic spline through the nodes) to the y and x
    stations of Ghia et al. Tables I and II (published data, Re = 100) and return the differences
    and their maximum magnitudes."""
    raise NotImplementedError


def vortex_centre(sol: dict) -> tuple[float, float, float]:
    """Primary-vortex centre: the minimum of psi, refined by fitting a quadratic surface to the 3 x 3
    nodes around the smallest nodal value and taking its stationary point. Returns (psi_min, x, y)."""
    raise NotImplementedError


def centre_values(sol: dict) -> tuple[float, float]:
    """psi and u at the geometric centre (0.5, 0.5), a node on every odd grid: the quantities used for
    the grid-convergence study (no interpolation error)."""
    m = sol["psi"].shape[0] // 2
    u, _ = velocities(sol["psi"], sol["U"])
    return float(sol["psi"][m, m]), float(u[m, m])


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".", grids: tuple[int, ...] = (17, 33, 65, 129)) -> None:
    """Print the lesson's numbers (grid sweep at Re = 100, Ghia comparison, observed order, Newton
    against marching, mass conservation, the Stokes limit) and save cavity.png. ``grids`` are odd N,
    each one refining the last by 2; the default ends at Ghia's 129 x 129."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    Re = 100.0
    print("lid-driven cavity, Re = 100, stream function-vorticity, central differences")
    print("    N   unknowns  Newton its  time (s)   psi_min     x_c     y_c    max|u-Ghia|  max|v-Ghia|  max|div|")
    sols = {}
    for N in grids:
        s = solve_cavity_newton(N, Re)
        sols[N] = s
        pm, xc, yc = vortex_centre(s)
        g = compare_with_ghia(s)
        u, v = velocities(s["psi"])
        div = np.max(np.abs(divergence(u, v)))
        print(f"  {N:3d}   {2 * (N - 2) ** 2:6d}     {s['iterations']:3d}     {s['seconds']:7.3f}   {pm:.6f}  {xc:.4f}"
              f"  {yc:.4f}    {g['max_du']:.5f}      {g['max_dv']:.5f}    {div:.1e}")
    print(f"Ghia et al. (1982), published data, 129 x 129: psi_min = {GHIA_RE100_PSI_MIN}, "
          f"centre ({GHIA_RE100_VORTEX_XY[0]}, {GHIA_RE100_VORTEX_XY[1]})")
    Nm, Nf = grids[-2], grids[-1]
    gm, gf = compare_with_ghia(sols[Nm]), compare_with_ghia(sols[Nf])
    du_ex = gf["du"] + (gf["du"] - gm["du"]) / 3.0          # Richardson, order 2, r = 2
    dv_ex = gf["dv"] + (gf["dv"] - gm["dv"]) / 3.0
    print(f"Richardson-extrapolated profiles ({Nm}, {Nf}) minus Ghia: max|du| = {np.max(np.abs(du_ex)):.5f}, "
          f"max|dv| = {np.max(np.abs(dv_ex)):.5f}")

    print("\ngrid convergence at the centre node (0.5, 0.5):")
    for k in range(len(grids) - 2):
        trio = grids[k:k + 3]
        ps = [centre_values(sols[n])[0] for n in trio]
        us = [centre_values(sols[n])[1] for n in trio]
        pp, pe = observed_order(*ps)
        pu, ue = observed_order(*us)
        print(f"  grids {trio}: psi {ps[0]:.6f} {ps[1]:.6f} {ps[2]:.6f}  order {pp:.2f}  extrapolated {pe:.6f};"
              f"  u {us[0]:.5f} {us[1]:.5f} {us[2]:.5f}  order {pu:.2f}  extrapolated {ue:.5f}")
    trio = grids[-3:]
    pms = [vortex_centre(sols[n])[0] for n in trio]
    pp, pe = observed_order(*pms)
    print(f"  psi_min on {trio}: order {pp:.2f}, Richardson extrapolate {pe:.6f} "
          f"(Ghia {GHIA_RE100_PSI_MIN}, difference {pe - GHIA_RE100_PSI_MIN:+.6f})")

    m = march_cavity(17, Re)
    nw = sols[17] if 17 in sols else solve_cavity_newton(17, Re)
    print(f"\nN = 17: marching ({m['steps']} steps of dt = {m['dt']:.4f}, {m['seconds']:.2f} s) against Newton "
          f"({nw['iterations']} iterations, {nw['seconds']:.3f} s): max|psi difference| = "
          f"{np.max(np.abs(m['psi'] - nw['psi'])):.1e}")
    r1, r2 = cavity_residual(sols[Nf]["psi"], sols[Nf]["omega"], Re)
    print(f"N = {Nf} residuals: kinematics {r1:.1e}, vorticity transport {r2:.1e}")
    st = solve_cavity_newton(Nm, 1e-3)
    pm, xc, yc = vortex_centre(st)
    print(f"Stokes limit (Re = 0.001, N = {Nm}): psi_min {pm:.6f} at ({xc:.4f}, {yc:.4f}); "
          f"symmetry error max|psi(x) - psi(1-x)| = {np.max(np.abs(st['psi'] - st['psi'][::-1, :])):.1e}")
    pm, xc, yc = vortex_centre(sols[Nf])
    c = Nf // 2
    print(f"Re = 100 (N = {Nf}): vortex centre moved downstream to x = {xc:.4f}; |omega| at the centre node "
          f"{abs(sols[Nf]['omega'][c, c]):.4f}")

    s = sols[Nf]
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    X, Y = np.meshgrid(s["x"], s["x"], indexing="ij")
    levels = np.sort(np.concatenate([np.linspace(-0.1, -0.01, 10), [-1e-3, -1e-4, 0, 1e-6, 1e-5]]))
    ax[0].contour(X, Y, s["psi"], levels=levels, colors="k", linewidths=0.7)
    ax[0].set_aspect("equal")
    ax[0].set_title(f"streamlines, Re = 100, {Nf} x {Nf} (computed)")
    ax[0].set_xlabel("x")
    ax[0].set_ylabel("y")
    styles = [":", "--", "-.", "-"]
    for k, N in enumerate(grids):
        y, uc, x, vc = centreline_profiles(sols[N])
        ax[1].plot(uc, y, styles[k % 4], label=f"computed, {N} x {N}")
        ax[2].plot(x, vc, styles[k % 4], label=f"computed, {N} x {N}")
    ax[1].plot(GHIA_RE100_U[:, 1], GHIA_RE100_U[:, 0], "ko", mfc="none", label="Ghia et al. 1982 (published data)")
    ax[2].plot(GHIA_RE100_V[:, 0], GHIA_RE100_V[:, 1], "ko", mfc="none", label="Ghia et al. 1982 (published data)")
    ax[1].set_xlabel("u(0.5, y)")
    ax[1].set_ylabel("y")
    ax[1].set_title("vertical centreline")
    ax[2].set_xlabel("x")
    ax[2].set_ylabel("v(x, 0.5)")
    ax[2].set_title("horizontal centreline")
    ax[1].legend(fontsize=8)
    ax[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "cavity.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'cavity.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
