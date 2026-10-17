"""P06.15 Continuous Systems: starter code.

Bending vibration of an Euler-Bernoulli beam, EI w'''' + rho A w_tt = 0, three ways:

* closed form: separation of variables, the frequency equations of the cantilever
  (1 + cos bL cosh bL = 0), the pinned-pinned beam (sin bL = 0), the clamped-clamped beam
  (cos bL cosh bL = 1) and the cantilever with a tip mass, solved with brentq; mode shapes and their
  orthogonality;
* a finite-difference eigenproblem for the pinned-pinned beam (the five-point fourth difference), whose
  discrete frequencies are known exactly and converge at second order;
* a finite-element eigenproblem with Hermite cubic beam elements and consistent mass matrices, for any
  of the supports, variable EI and rho A, and point masses, converging at fourth order in frequency;
* a wing-like cantilever: the P05.12 spar with the wing's mass, a tip mass, a mid-span tank and a taper.

The wing numbers are illustrative: the light aeroplane's spar of P05.12 (EI = 6.958e5 N m^2 at the
root, 5 m semi-span) carrying a 45 kg wing spread uniformly.

Conventions: x from the root (clamped end of a cantilever); w transverse deflection; FE degrees of
freedom (w, w') at each node; SI units; frequencies in rad/s unless named f (Hz).

Run:  python main.py [out_dir]     (prints the tables, writes beam_modes.png)
Test: pytest curriculum/phase-06/lesson-15-continuous-systems
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Callable

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq

G0 = 9.80665
BCS = ("clamped-free", "pinned-pinned", "clamped-clamped")


def wing_spar() -> dict[str, float]:
    """Illustrative wing: the P05.12 spar (EI = 6.958e5 N m^2), semi-span 5 m, 45 kg of wing spread
    uniformly (9 kg/m)."""
    return dict(EI=6.958e5, L=5.0, rhoA=9.0)


# ---------------------------------------------------------------- closed form

def frequency_equation(bL: float | np.ndarray, bc: str, mu: float = 0.0) -> np.ndarray:
    """Characteristic function whose roots are the dimensionless eigenvalues beta L, written divided by
    cosh(bL) so that it stays of order one:
    'clamped-free': cos + 1/cosh + mu bL (cos tanh - sin) with mu = m_tip/(rho A L) (mu = 0: 1 + cos cosh = 0);
    'pinned-pinned': sin bL;  'clamped-clamped': cos - 1/cosh (cos cosh = 1)."""
    raise NotImplementedError


def beta_L(bc: str, n: int, mu: float = 0.0) -> np.ndarray:
    """First n positive roots beta_i L of the frequency equation, by brentq on a fine sign-change scan
    of (0.01, (n + 2) pi]."""
    raise NotImplementedError


def natural_frequencies(EI: float, rhoA: float, L: float, bc: str, n: int, m_tip: float = 0.0) -> np.ndarray:
    """omega_i = (beta_i L)^2 sqrt(EI / (rho A L^4)) (rad/s); a tip mass only for 'clamped-free'."""
    raise NotImplementedError


def mode_shape(x: np.ndarray, L: float, bL: float, bc: str) -> np.ndarray:
    """Mode shape of the uniform beam for the root bL (no tip mass), scaled so that max |phi| = 1 with a
    positive value at the largest-|phi| point nearest the free end:
    'clamped-free': cosh bx - cos bx - s (sinh bx - sin bx), s = (cosh bL + cos bL)/(sinh bL + sin bL);
    'pinned-pinned': sin bx;  'clamped-clamped': cosh - cos - s (sinh - sin), s = (cosh - cos)/(sinh - sin)."""
    raise NotImplementedError


# ---------------------------------------------------------------- finite differences (pinned-pinned)

def fd_pinned_frequencies(EI: float, rhoA: float, L: float, n_int: int, n_modes: int) -> np.ndarray:
    """Finite-difference frequencies of a pinned-pinned beam: n_int + 1 intervals, unknowns at the
    n_int interior nodes, w'''' by the five-point stencil (1, -4, 6, -4, 1)/h^4 with w = 0 at the ends
    and w'' = 0 there by the ghost-node rule w_{-1} = -w_1. Returns the lowest n_modes omega (rad/s)."""
    raise NotImplementedError


def fd_pinned_exact(EI: float, rhoA: float, L: float, n_int: int, n_modes: int) -> np.ndarray:
    """The discrete frequencies of ``fd_pinned_frequencies`` in closed form: the sine modes are exact
    eigenvectors of the stencil, with eigenvalues (2 - 2 cos(j pi h / L))^2 / h^4 (given)."""
    h = L / (n_int + 1)
    j = np.arange(1, n_modes + 1)
    return np.sqrt(EI / rhoA) * (2 - 2 * np.cos(j * np.pi * h / L)) / h ** 2


# ---------------------------------------------------------------- finite elements

def beam_element(EI: float, rhoA: float, h: float) -> tuple[np.ndarray, np.ndarray]:
    """Hermite cubic Euler-Bernoulli element, dofs (w1, w1', w2, w2'): stiffness
    EI/h^3 [[12, 6h, -12, 6h], [6h, 4h^2, -6h, 2h^2], [-12, -6h, 12, -6h], [6h, 2h^2, -6h, 4h^2]] and
    consistent mass rho A h/420 [[156, 22h, 54, -13h], [22h, 4h^2, 13h, -3h^2], [54, 13h, 156, -22h],
    [-13h, -3h^2, -22h, 4h^2]]."""
    raise NotImplementedError


def fem_matrices(L: float, EI: float | Callable[[np.ndarray], np.ndarray], rhoA: float | Callable[[np.ndarray], np.ndarray],
                 n_el: int, point_masses: list[tuple[float, float]] = ()) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Assemble global K and M (2 (n_el + 1) square) on n_el equal elements. EI and rho A may be
    functions of x, evaluated at each element's midpoint. Point masses (x, m) are added to the
    translational dof of the nearest node. Returns (nodes, K, M) before supports are applied."""
    raise NotImplementedError


def free_dofs(n_el: int, bc: str) -> np.ndarray:
    """Indices of the unconstrained dofs: a clamped end fixes w and w', a pinned end fixes w only."""
    raise NotImplementedError


def fem_frequencies(L: float, EI, rhoA, n_el: int, bc: str, n_modes: int,
                    point_masses: list[tuple[float, float]] = ()) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Lowest n_modes natural frequencies (rad/s) of the FE beam, with the nodal deflections of each mode
    (columns, scaled to max |w| = 1, positive at the largest point nearest the free end). Uses
    scipy.linalg.eigh(K, M) on the free dofs. Returns (omega, nodes, W)."""
    raise NotImplementedError


def rayleigh_tip_mass_estimate(EI: float, rhoA: float, L: float, m_tip: float) -> float:
    """SDOF estimate of a cantilever's first frequency from P06.13: k = 3 EI / L^3 on an effective
    mass m_tip + (33/140) rho A L. An upper bound (Rayleigh), exact as m_tip / (rho A L) -> infinity."""
    raise NotImplementedError


def observed_order(errors: np.ndarray, ratio: float = 2.0) -> np.ndarray:
    """Observed convergence orders of errors on grids refined by ``ratio`` (given)."""
    e = np.asarray(errors, dtype=float)
    return np.log(e[:-1] / e[1:]) / np.log(ratio)


def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save beam_modes.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    for bc in BCS:
        print(f"{bc:16s} beta L = {np.array2string(beta_L(bc, 4), precision=6)}")
    w = wing_spar()
    EI, L, rA = w["EI"], w["L"], w["rhoA"]
    om = natural_frequencies(EI, rA, L, "clamped-free", 3)
    print(f"\nwing spar (EI {EI:.4g} N m^2, {rA} kg/m, {L} m): bending f = {np.array2string(om / 2 / np.pi, precision=4)} Hz")
    print(f"   SDOF estimate (P06.13, no tip mass) {rayleigh_tip_mass_estimate(EI, rA, L, 0.0) / 2 / np.pi:.4f} Hz")

    print("\nFEM convergence, cantilever, mode 1 and 3 relative errors")
    ns = [2, 4, 8, 16, 32]
    e1, e3 = [], []
    for n in ns:
        wf, _, _ = fem_frequencies(L, EI, rA, n, "clamped-free", 3)
        e1.append(wf[0] / om[0] - 1)
        e3.append(wf[2] / om[2] - 1)
        print(f"   n_el {n:3d}: mode 1 {e1[-1]:.3e}   mode 3 {e3[-1]:.3e}")
    print(f"   observed orders mode 1 {np.array2string(observed_order(e1), precision=2)}, mode 3 {np.array2string(observed_order(e3[1:]), precision=2)}")
    print("FD convergence, pinned-pinned, mode 1 relative error")
    ep = natural_frequencies(EI, rA, L, "pinned-pinned", 1)[0]
    efd = []
    for n in (9, 19, 39, 79, 159):
        efd.append(abs(fd_pinned_frequencies(EI, rA, L, n, 1)[0] / ep - 1))
        print(f"   {n + 1:4d} intervals: {efd[-1]:.3e}")
    print(f"   observed orders {np.array2string(observed_order(efd), precision=2)}")

    # tip mass and tank
    for mt in (10.0, 40.0):
        wt = natural_frequencies(EI, rA, L, "clamped-free", 2, m_tip=mt)
        wf, _, _ = fem_frequencies(L, EI, rA, 40, "clamped-free", 2, point_masses=[(L, mt)])
        print(f"\ntip mass {mt:.0f} kg: closed form f = {np.array2string(wt / 2 / np.pi, precision=4)} Hz, FEM (40 el) "
              f"{np.array2string(wf / 2 / np.pi, precision=4)} Hz; SDOF estimate f1 {rayleigh_tip_mass_estimate(EI, rA, L, mt) / 2 / np.pi:.4f} Hz")
    wf, _, _ = fem_frequencies(L, EI, rA, 40, "clamped-free", 3, point_masses=[(2.0, 40.0)])
    print(f"40 kg tank at 2 m (FEM): f = {np.array2string(wf / 2 / np.pi, precision=4)} Hz")
    EIt = lambda x: EI * (1 - 0.7 * x / L)
    rAt = lambda x: rA * 1.5 * (1 - 2.0 / 3.0 * x / L)          # 45 kg, three times as heavy per metre at the root as at the tip
    wt80, _, _ = fem_frequencies(L, EIt, rAt, 80, "clamped-free", 3)
    wt160, _, _ = fem_frequencies(L, EIt, rAt, 160, "clamped-free", 3)
    print(f"tapered spar (EI to 30 %, mass per metre 3:1 root to tip): f = {np.array2string(wt160 / 2 / np.pi, precision=4)} Hz "
          f"(80 el: {np.array2string(wt80 / 2 / np.pi, precision=4)})")

    x = np.linspace(0, L, 401)
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    bls = beta_L("clamped-free", 3)
    wf, nodes, W = fem_frequencies(L, EI, rA, 8, "clamped-free", 3)
    for i in range(3):
        ax[0].plot(x, mode_shape(x, L, bls[i], "clamped-free"), label=f"mode {i + 1}: {om[i] / 2 / np.pi:.2f} Hz")
        ax[0].plot(nodes, W[:, i], "o", ms=4, color=f"C{i}")
    ax[0].set_xlabel("span x (m)")
    ax[0].set_ylabel("w / max|w|")
    ax[0].set_title("cantilever modes: closed form and FEM (8 elements)")
    ax[0].legend(fontsize=8)
    h = L / np.array(ns)
    ax[1].loglog(h, np.abs(e1), "o-", label="FEM mode 1")
    ax[1].loglog(h, np.abs(e3), "s-", label="FEM mode 3")
    hf = L / np.array([10, 20, 40, 80, 160])
    ax[1].loglog(hf, efd, "^-", label="FD pinned mode 1")
    ax[1].loglog(h, abs(e1[0]) * (h / h[0]) ** 4, "k--", lw=0.8, label="slope 4")
    ax[1].loglog(hf, efd[0] * (hf / hf[0]) ** 2, "k:", lw=0.8, label="slope 2")
    ax[1].set_xlabel("element / grid size h (m)")
    ax[1].set_ylabel("relative frequency error")
    ax[1].set_title("convergence to the frequency equation")
    ax[1].legend(fontsize=8)
    bx = np.linspace(0.05, 12, 2000)
    ax[2].plot(bx, frequency_equation(bx, "clamped-free"), label="1 + cos cosh (/cosh)")
    ax[2].plot(bx, frequency_equation(bx, "clamped-free", mu=40.0 / 45.0), label="tip mass mu = 0.89")
    ax[2].axhline(0, color="0.5", lw=0.8)
    ax[2].set_ylim(-1.5, 1.5)
    ax[2].set_xlabel("beta L")
    ax[2].set_title("cantilever frequency equation")
    ax[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "beam_modes.png", dpi=110)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'beam_modes.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
