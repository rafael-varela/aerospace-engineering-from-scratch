"""P05.14 Failure Theories: starter code.

A structural part fails when some measure of its stress state reaches what the material can carry.
For a ductile metal under combined stress the question "has it yielded?" is answered by a yield
criterion that turns the six components of the stress tensor into one equivalent stress that can be
compared with the uniaxial yield strength from a tensile test:

* von Mises (distortion energy): sigma_vm = sqrt(3 J2), the stress that stores the same distortion
  energy as uniaxial tension;
* Tresca (maximum shear): sigma_T = sigma_1 - sigma_3, twice the largest shear stress;
* for a brittle material that is much stronger in compression than in tension, Coulomb-Mohr:
  sigma_1 / S_t - sigma_3 / S_c = 1 / n, with a tension cut-off.

Aerospace sizing reports margins of safety, MS = allowable / (FS x applied) - 1, with a factor of
safety of 1.0 on yield at limit load and 1.5 on ultimate. MS >= 0 passes.

The worked case is a 4130 steel landing-gear tube at the fuselage attachment carrying bending,
torsion and axial compression at limit load.

Units: SI throughout (Pa, m, N, N m). Stress tensors are 3 x 3 symmetric arrays; plane stress is
(sigma_x, sigma_y, tau_xy) with sigma_z = tau_xz = tau_yz = 0. Tension positive.

Run:  python main.py [out_dir]     (prints the margin tables, writes failure_theories.png)
Test: pytest curriculum/phase-05/lesson-14-failure-theories
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import brentq


# ---------------------------------------------------------------- materials

@dataclass(frozen=True)
class Material:
    """Typical room-temperature properties (not design allowables). Strengths in Pa, E in Pa.

    Values are typical figures of the kind published on manufacturers' and ASM data sheets; a design
    uses the statistically based A- or B-basis allowables of MMPDS, which are lower."""
    name: str
    E: float
    nu: float
    F_ty: float
    F_tu: float
    rho: float


MATERIALS = {
    "2024-T3": Material("2024-T3", 73.1e9, 0.33, 345e6, 483e6, 2780.0),
    "7075-T6": Material("7075-T6", 71.7e9, 0.33, 503e6, 572e6, 2810.0),
    "Ti-6Al-4V": Material("Ti-6Al-4V", 113.8e9, 0.342, 880e6, 950e6, 4430.0),
    "4130": Material("4130", 205e9, 0.29, 435e6, 670e6, 7850.0),
}

FS_YIELD = 1.0       # limit load must not cause yield (no detrimental permanent set)
FS_ULTIMATE = 1.5    # ultimate load = 1.5 x limit must not cause rupture


# ---------------------------------------------------------------- stress states

def plane_stress_tensor(sx: float, sy: float, txy: float) -> np.ndarray:
    """The 3 x 3 tensor of a plane stress state: sigma_z = tau_xz = tau_yz = 0."""
    return np.array([[sx, txy, 0.0], [txy, sy, 0.0], [0.0, 0.0, 0.0]], dtype=float)


def principal_stresses(sigma: np.ndarray) -> np.ndarray:
    """Principal stresses sigma_1 >= sigma_2 >= sigma_3 of a symmetric 3 x 3 tensor (eigenvalues)."""
    raise NotImplementedError


def principal_stresses_plane(sx: float, sy: float, txy: float) -> tuple[float, float]:
    """In-plane principal stresses from Mohr's circle: centre (sx + sy)/2, radius
    sqrt(((sx - sy)/2)^2 + txy^2). Returns (sigma_a, sigma_b) with sigma_a >= sigma_b."""
    raise NotImplementedError


def von_mises(sigma: np.ndarray) -> float:
    """von Mises equivalent stress from the six components of a 3 x 3 tensor:
    sqrt(((sx - sy)^2 + (sy - sz)^2 + (sz - sx)^2 + 6 (txy^2 + tyz^2 + tzx^2)) / 2)."""
    raise NotImplementedError


def von_mises_principal(p: np.ndarray) -> float:
    """von Mises stress from the three principal stresses: sqrt(((s1-s2)^2 + (s2-s3)^2 + (s3-s1)^2)/2)."""
    raise NotImplementedError


def tresca(sigma: np.ndarray) -> float:
    """Tresca equivalent stress sigma_1 - sigma_3 (twice the absolute maximum shear stress), using all
    three principal stresses, so that plane stress counts the zero out-of-plane principal stress."""
    raise NotImplementedError


def max_shear_stress(sigma: np.ndarray) -> float:
    """Absolute maximum shear stress (sigma_1 - sigma_3) / 2."""
    raise NotImplementedError


def strain_energy_density(sigma: np.ndarray, E: float, nu: float) -> float:
    """Total elastic strain energy per unit volume, u = (1/2) sigma : epsilon, with epsilon from
    isotropic Hooke's law epsilon = ((1 + nu) sigma - nu tr(sigma) I) / E."""
    raise NotImplementedError


def volumetric_energy_density(sigma: np.ndarray, E: float, nu: float) -> float:
    """Energy of the volume change alone: the mean stress p = tr(sigma)/3 acting through the
    volumetric strain, u_v = (1 - 2 nu) / (6 E) (tr sigma)^2."""
    raise NotImplementedError


def distortion_energy_density(sigma: np.ndarray, E: float, nu: float) -> float:
    """Distortion (shape-change) energy per unit volume: total minus volumetric. The von Mises
    criterion says yield starts when this equals its value at uniaxial yield, (1 + nu) F_ty^2 / (3 E)."""
    raise NotImplementedError


# ---------------------------------------------------------------- margins

def margin_of_safety(allowable: float, applied: float, factor: float = 1.0) -> float:
    """MS = allowable / (factor * applied) - 1. Returns +inf for a zero applied stress."""
    raise NotImplementedError


def margins(sigma_limit: np.ndarray, mat: Material, fs_yield: float = FS_YIELD,
            fs_ult: float = FS_ULTIMATE) -> dict[str, float]:
    """Yield and ultimate margins of a stress state at LIMIT load, by von Mises and by Tresca.

    Yield: equivalent stress at limit load against F_ty with fs_yield. Ultimate: the stresses are taken
    as linear in load, so the stress at ultimate load is fs_ult times the limit value, against F_tu."""
    raise NotImplementedError


def load_factor_to_yield(sigma_unit: np.ndarray, F_y: float, criterion: str = "von_mises") -> float:
    """Load multiplier lambda at which lambda * sigma_unit first reaches F_y, found by root finding
    (brentq) on the criterion. For a proportional load it must equal F_y / sigma_eq(sigma_unit):
    the independent check that both criteria are homogeneous of degree one."""
    raise NotImplementedError


# ---------------------------------------------------------------- brittle materials

def coulomb_mohr_factor(sigma: np.ndarray, S_t: float, S_c: float) -> float:
    """Factor of safety n of a brittle material by Coulomb-Mohr with a tension cut-off.

    With principal stresses sigma_1 >= sigma_2 >= sigma_3 (all three, including any zero):
    1/n = max(sigma_1/S_t - sigma_3/S_c, sigma_1/S_t). S_t and S_c are the (positive) ultimate
    strengths in tension and compression. Uniaxial tension gives S_t/sigma, uniaxial compression
    S_c/|sigma|, pure shear tau gives S_t S_c / ((S_t + S_c) tau). With S_t = S_c = S it is Tresca."""
    raise NotImplementedError


def rankine_factor(sigma: np.ndarray, S_t: float, S_c: float) -> float:
    """Factor of safety by the maximum-normal-stress (Rankine) criterion: the largest principal
    stress against S_t, the most compressive against S_c."""
    raise NotImplementedError


# ---------------------------------------------------------------- the landing-gear tube

def tube_section(d_o: float, t: float) -> dict[str, float]:
    """Area, second moment I, polar moment J and outer radius c of a circular tube of outside
    diameter d_o and wall t."""
    raise NotImplementedError


def tube_stress_at(theta: float, r: float, M: float, T: float, P: float, d_o: float, t: float) -> np.ndarray:
    """Stress tensor at polar angle theta and radius r of a tube under bending moment M (about the
    section's y axis, so the fibre at theta = 0, z = r cos theta is the most tensile for M > 0), torque T
    and axial force P (tension positive). The tube axis is x; the returned tensor is in local
    (axial, hoop, radial) coordinates: sigma_axial = P/A + M z / I, tau_axial-hoop = T r / J."""
    raise NotImplementedError


def critical_point(M: float, T: float, P: float, d_o: float, t: float, n: int = 721) -> tuple[float, float]:
    """Scan the outer surface and return (theta_critical, sigma_vm_max). The outer fibre carries the
    largest bending and shear stress, so the critical point is on it."""
    raise NotImplementedError


def tube_resultants(M: float, T: float, P: float, d_o: float, t: float,
                    n_r: int = 40, n_th: int = 360) -> tuple[float, float, float]:
    """Integrate the stress distribution back over the annulus (midpoint rule in r and theta) to recover
    the axial force, the bending moment and the torque: the equilibrium check of the section stresses."""
    raise NotImplementedError


# ---------------------------------------------------------------- envelopes

def sigma_tau_envelope(criterion: str, F: float, n: int = 361, S_c: float | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Failure envelope in the (sigma, tau) plane of a state sigma_x = sigma, tau_xy = tau (all else 0),
    the state of a shaft in bending and torsion. von Mises: sigma^2 + 3 tau^2 = F^2; Tresca:
    sigma^2 + 4 tau^2 = F^2; Coulomb-Mohr with S_t = F and S_c: solved numerically on rays.
    Returns (sigma, tau) arrays tracing the closed curve."""
    raise NotImplementedError


# ---------------------------------------------------------------- driver

def main(out_dir: str | Path = ".") -> None:
    """Print the worked case and the experiment tables; save failure_theories.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    MPa = 1e6

    # ---- the landing-gear tube at limit load
    mat = MATERIALS["4130"]
    d_o, t = 0.050, 0.003
    M, T, P = 1.2e3, 0.8e3, -20e3
    s = tube_section(d_o, t)
    print(f"4130 tube {d_o*1e3:.0f} x {t*1e3:.0f} mm: A = {s['A']*1e6:.1f} mm^2, I = {s['I']*1e12:.0f} mm^4, "
          f"J = {s['J']*1e12:.0f} mm^4")
    sb = M * s["c"] / s["I"]
    sa = P / s["A"]
    tq = T * s["c"] / s["J"]
    print(f"limit load M = {M/1e3:.1f} kN m, T = {T/1e3:.1f} kN m, P = {P/1e3:.0f} kN: bending {sb/MPa:.1f} MPa, "
          f"axial {sa/MPa:.1f} MPa, torsional shear {tq/MPa:.1f} MPa")
    th_c, vm_c = critical_point(M, T, P, d_o, t)
    sig_c = tube_stress_at(th_c, d_o / 2, M, T, P, d_o, t)
    print(f"critical point at theta = {np.degrees(th_c):.0f} deg: sigma = {sig_c[0,0]/MPa:.1f} MPa, tau = {sig_c[0,1]/MPa:.1f} MPa")
    sig_t = tube_stress_at(0.0, d_o / 2, M, T, P, d_o, t)
    print(f"tension side (theta = 0): sigma = {sig_t[0,0]/MPa:.1f} MPa, von Mises {von_mises(sig_t)/MPa:.1f} MPa")
    p = principal_stresses(sig_c)
    print(f"principal stresses at the critical point: {p[0]/MPa:.1f}, {p[1]/MPa:.1f}, {p[2]/MPa:.1f} MPa")
    m = margins(sig_c, mat)
    print(f"von Mises {m['von_mises']/MPa:.1f} MPa, Tresca {m['tresca']/MPa:.1f} MPa (ratio {m['tresca']/m['von_mises']:.3f})")
    print(f"  MS yield (FS 1.0, F_ty {mat.F_ty/MPa:.0f}):   von Mises {m['MS_yield_vm']:+.3f}   Tresca {m['MS_yield_tresca']:+.3f}")
    print(f"  MS ultimate (FS 1.5, F_tu {mat.F_tu/MPa:.0f}): von Mises {m['MS_ult_vm']:+.3f}   Tresca {m['MS_ult_tresca']:+.3f}")
    N, Mb, Tq = tube_resultants(M, T, P, d_o, t)
    print(f"resultants from the integrated stresses: N = {N/1e3:.4f} kN, M = {Mb/1e3:.4f} kN m, T = {Tq/1e3:.4f} kN m")
    ud = distortion_energy_density(sig_c, mat.E, mat.nu)
    print(f"distortion energy {ud/1e3:.2f} kJ/m^3 = (1+nu) sigma_vm^2 / 3E = {(1+mat.nu)*vm_c**2/(3*mat.E)/1e3:.2f} kJ/m^3")

    # ---- experiment 1: Tresca / von Mises across the stress ratio
    print("\ntau/sigma   von Mises/sigma   Tresca/sigma   ratio")
    for r in (0.0, 0.25, 0.5, 1.0, 2.0, 1e6):
        st = plane_stress_tensor(1.0, 0.0, r) if r < 1e5 else plane_stress_tensor(0.0, 0.0, 1.0)
        vm, tr = von_mises(st), tresca(st)
        lab = "inf" if r >= 1e5 else f"{r:.2f}"
        print(f"{lab:>9}   {vm:15.4f}   {tr:12.4f}   {tr/vm:.4f}")

    # ---- experiment 2: material comparison for the same tube and loads
    print("\nmaterial     F_ty  F_tu   MS_y(vM)  MS_u(vM)  MS_y(Tr)  MS_u(Tr)  governs")
    for k, mt in MATERIALS.items():
        mm = margins(sig_c, mt)
        gov = "ultimate" if mm["MS_ult_vm"] < mm["MS_yield_vm"] else "yield"
        print(f"{k:10s}  {mt.F_ty/MPa:5.0f} {mt.F_tu/MPa:5.0f}  {mm['MS_yield_vm']:+8.3f}  {mm['MS_ult_vm']:+8.3f}  "
              f"{mm['MS_yield_tresca']:+8.3f}  {mm['MS_ult_tresca']:+8.3f}  {gov}")
    print("ultimate governs when F_tu / F_ty < 1.5:",
          ", ".join(f"{k} {mt.F_tu/mt.F_ty:.2f}" for k, mt in MATERIALS.items()))

    # ---- experiment 3: a brittle cast part (typical gray cast iron, S_t = 214, S_c = 752 MPa)
    St, Sc = 214e6, 752e6
    print("\nbrittle check, S_t = 214 MPa, S_c = 752 MPa; state at 100 MPa scale")
    for lab, st in [("uniaxial tension", plane_stress_tensor(100e6, 0, 0)),
                    ("uniaxial compression", plane_stress_tensor(-100e6, 0, 0)),
                    ("pure shear", plane_stress_tensor(0, 0, 100e6)),
                    ("biaxial tension", plane_stress_tensor(100e6, 100e6, 0)),
                    ("tension + compression", plane_stress_tensor(100e6, -100e6, 0))]:
        print(f"  {lab:22s} Coulomb-Mohr n = {coulomb_mohr_factor(st, St, Sc):.3f}   Rankine n = {rankine_factor(st, St, Sc):.3f}")

    # ---- figure
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
    for crit, sty in (("von_mises", "-"), ("tresca", "--")):
        sg, ta = sigma_tau_envelope(crit, mat.F_ty)
        ax[0].plot(sg / MPa, ta / MPa, sty, label=crit.replace("_", " "))
    sg, ta = sigma_tau_envelope("coulomb_mohr", St, S_c=Sc)
    ax[0].plot(sg / MPa, ta / MPa, ":", label="Coulomb-Mohr (cast iron)")
    ax[0].plot(sig_c[0, 0] / MPa, sig_c[0, 1] / MPa, "ko", label="tube, limit load")
    ax[0].plot(1.5 * sig_c[0, 0] / MPa, 1.5 * sig_c[0, 1] / MPa, "k^", label="tube, ultimate load")
    ax[0].set_xlabel("sigma [MPa]"); ax[0].set_ylabel("tau [MPa]"); ax[0].set_aspect("equal")
    ax[0].set_title("envelopes for sigma + tau (4130 F_ty)"); ax[0].legend(fontsize=7); ax[0].grid(alpha=0.3)
    th = np.linspace(0, 2 * np.pi, 361)
    vm = [von_mises(tube_stress_at(a, d_o / 2, M, T, P, d_o, t)) / MPa for a in th]
    trs = [tresca(tube_stress_at(a, d_o / 2, M, T, P, d_o, t)) / MPa for a in th]
    ax[1].plot(np.degrees(th), vm, label="von Mises"); ax[1].plot(np.degrees(th), trs, "--", label="Tresca")
    ax[1].axhline(mat.F_ty / MPa, color="k", lw=0.8, label="F_ty")
    ax[1].set_xlabel("position around the tube [deg]"); ax[1].set_ylabel("equivalent stress [MPa]")
    ax[1].set_title("outer surface of the gear tube, limit load"); ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / "failure_theories.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'failure_theories.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
