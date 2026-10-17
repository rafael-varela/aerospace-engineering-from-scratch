"""P05.07 Axial Loading & Thermal Stress: reference solution.

A bar under axial load stretches by delta = N L / (E A); a bar that is not free to stretch carries
a stress set by compatibility, and a bar that is heated but restrained carries -E alpha dT. The
bolted joint is the aerospace device in which all three meet: a preloaded bolt and the clamped
members are two springs in parallel, an external load is shared between them in proportion to
their stiffnesses, and a temperature change alters the preload when bolt and members expand by
different amounts. This module builds:

* axial deformation of uniform, stepped and tapered bars, and a 1-D bar finite-element solver
  (direct stiffness, any connectivity, initial strains for thermal expansion or misfit) used as the
  independent check of every closed form;
* the fixed-fixed two-segment bar with a load at the junction and a temperature change, in closed
  form from compatibility;
* the bolt: tensile stress area, stiffness of the shank and threaded part in series;
* the clamped members: a 30-degree pressure-cone (frustum) stiffness, integrated in closed form
  for any stack of layers;
* the joint diagram: joint constant, bolt and member loads with an external separating load,
  separation load, thermal change of preload, torque-preload relation.

Units: SI-consistent mm, N, MPa (N/mm^2), N/mm for stiffness, N mm for torque, kelvin for
temperature differences. Material values are typical, not design allowables.

Run:  python solution.py [out_dir]     (prints the lesson's tables, writes bolted_joint.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-05/lesson-07-axial-loading-thermal-stress
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

# Typical room-temperature properties: E (MPa), alpha (1/K, about 20-100 C), Fty, Ftu (MPa).
MATERIALS: dict[str, dict[str, float]] = {
    "2024-T3":   {"E": 73.1e3, "alpha": 23.2e-6, "Fty": 345.0, "Ftu": 483.0},
    "7075-T6":   {"E": 71.7e3, "alpha": 23.6e-6, "Fty": 503.0, "Ftu": 572.0},
    "Ti-6Al-4V": {"E": 113.8e3, "alpha": 8.6e-6, "Fty": 880.0, "Ftu": 950.0},
    "4130":      {"E": 205.0e3, "alpha": 12.2e-6, "Fty": 435.0, "Ftu": 670.0},
}


# ---------------------------------------------------------------- bars

def axial_deformation(N: float, L: float, A: float, E: float) -> float:
    """Elongation of a uniform bar carrying internal force N: delta = N L / (E A)."""
    return N * L / (E * A)


def stepped_bar_deformation(N: np.ndarray, L: np.ndarray, A: np.ndarray, E: np.ndarray) -> float:
    """Total elongation of segments in series, each with its own internal force, length, area and
    modulus: delta = sum N_i L_i / (E_i A_i)."""
    N, L, A, E = (np.asarray(v, dtype=float) for v in (N, L, A, E))
    return float(np.sum(N * L / (E * A)))


def tapered_bar_deformation(P: float, L: float, A0: float, A1: float, E: float) -> float:
    """Elongation of a bar whose area varies linearly from A0 to A1, carrying P: integrating
    P dx / (E A(x)) gives delta = P L ln(A0 / A1) / (E (A0 - A1)) (and P L / (E A0) when A1 = A0)."""
    if np.isclose(A0, A1):
        return P * L / (E * A0)
    return P * L * np.log(A0 / A1) / (E * (A0 - A1))


def constrained_thermal_stress(E: float, alpha: float, dT: float) -> float:
    """Stress in a bar held at fixed length while its temperature changes by dT:
    sigma = -E alpha dT (compressive when heated)."""
    return -E * alpha * dT


def bar_fe_solve(x: np.ndarray, elements: np.ndarray | None, EA: np.ndarray, loads: np.ndarray,
                 fixed: list[int], eps0: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """1-D bar finite elements by the direct stiffness method.

    x: node coordinates (n,). elements: (m, 2) node pairs, or None for consecutive nodes. EA: axial
    rigidity of each element (m,). loads: external nodal forces (n,). fixed: indices of nodes with
    zero displacement. eps0: initial (stress-free) strain of each element, e.g. alpha dT for
    thermal expansion or -misfit / L for a member made short.

    Element stiffness (EA / L) [[1, -1], [-1, 1]]; the initial strain adds equivalent nodal forces
    EA eps0 (-1, +1). Returns (u, N, R): nodal displacements, element axial forces
    N = EA ((u_j - u_i) / L - eps0) (tension positive), and nodal reactions R = K u - f_equiv - loads
    (zero at free nodes). Equilibrium: R + loads sums to zero.
    """
    x = np.asarray(x, dtype=float)
    n = x.size
    el = np.column_stack([np.arange(n - 1), np.arange(1, n)]) if elements is None else np.asarray(elements, int)
    EA = np.broadcast_to(np.asarray(EA, dtype=float), (len(el),))
    e0 = np.zeros(len(el)) if eps0 is None else np.broadcast_to(np.asarray(eps0, dtype=float), (len(el),))
    loads = np.asarray(loads, dtype=float)
    K = np.zeros((n, n))
    f0 = np.zeros(n)
    Ls = np.abs(x[el[:, 1]] - x[el[:, 0]])
    for (i, j), ea, L, e in zip(el, EA, Ls, e0):
        k = ea / L
        K[np.ix_([i, j], [i, j])] += k * np.array([[1.0, -1.0], [-1.0, 1.0]])
        f0[i] -= ea * e
        f0[j] += ea * e
    free = np.setdiff1d(np.arange(n), fixed)
    u = np.zeros(n)
    u[free] = np.linalg.solve(K[np.ix_(free, free)], (loads + f0)[free])
    N = EA * ((u[el[:, 1]] - u[el[:, 0]]) / Ls - e0)
    R = K @ u - f0 - loads
    R[free] = 0.0
    return u, N, R


def fixed_fixed_two_segment(P: float, L1: float, EA1: float, alpha1: float, L2: float, EA2: float,
                            alpha2: float, dT: float = 0.0) -> tuple[float, float, float]:
    """A bar of two segments in series between two rigid walls, axial force P (+x) applied at the
    junction, uniform temperature change dT. Compatibility (total length fixed) with
    k_i = EA_i / L_i gives the junction displacement

        u = (P + k1 alpha1 L1 dT - k2 alpha2 L2 dT) / (k1 + k2),

    and the internal forces N1 = k1 (u - alpha1 L1 dT), N2 = k2 (-u - alpha2 L2 dT).
    Returns (u, N1, N2), tension positive."""
    k1, k2 = EA1 / L1, EA2 / L2
    u = (P + k1 * alpha1 * L1 * dT - k2 * alpha2 * L2 * dT) / (k1 + k2)
    return u, k1 * (u - alpha1 * L1 * dT), k2 * (-u - alpha2 * L2 * dT)


# ---------------------------------------------------------------- the bolt

def tensile_stress_area(d: float, pitch: float) -> float:
    """ISO metric thread tensile stress area A_t = (pi / 4)(d - 0.9382 p)^2 [mm^2]."""
    return np.pi / 4.0 * (d - 0.9382 * pitch) ** 2


def bolt_stiffness(d: float, pitch: float, l_shank: float, l_thread: float, E: float) -> float:
    """Axial stiffness of the bolt within the grip: the unthreaded shank (area pi d^2 / 4, length
    l_shank) and the threaded part in the grip (tensile stress area, length l_thread) in series."""
    Ad = np.pi * d ** 2 / 4.0
    At = tensile_stress_area(d, pitch)
    return 1.0 / (l_shank / (E * Ad) + l_thread / (E * At))


# ---------------------------------------------------------------- the members: pressure cones

def _cone_compliance(s_a: float, s_b: float, E: float, d: float, D_w: float, tan_a: float) -> float:
    """Compliance of the part of a frustum between distances s_a < s_b from its small end:
    integral of ds / (E pi ((D(s)/2)^2 - (d/2)^2)) with D(s) = D_w + 2 s tan(alpha), in closed form
    (partial fractions): [ln((D - d)/(D + d))] / (pi E d tan(alpha))."""
    g = lambda s: np.log((D_w + 2 * s * tan_a - d) / (D_w + 2 * s * tan_a + d))
    return (g(s_b) - g(s_a)) / (np.pi * E * d * tan_a)


def frustum_stiffness(E: float, d: float, D_w: float, t: float, alpha_deg: float = 30.0) -> float:
    """Stiffness of one frustum (cone of half-angle alpha) of height t, hole diameter d, small-end
    diameter D_w:  k = pi E d tan(a) / ln[(2 t tan(a) + D_w - d)(D_w + d) / ((2 t tan(a) + D_w + d)(D_w - d))]."""
    ta = np.tan(np.radians(alpha_deg))
    return 1.0 / _cone_compliance(0.0, t, E, d, D_w, ta)


def member_stiffness(d: float, D_w: float, layers: list[tuple[float, float]], alpha_deg: float = 30.0) -> float:
    """Clamped-member stiffness by the pressure-cone model for a stack of layers (thickness, E),
    listed from the head side. Two cones grow from the two bearing faces (diameter D_w) and meet at
    the mid-plane of the grip; every piece of every layer on each side of the mid-plane is a slice
    of a frustum, and all slices are springs in series."""
    ta = np.tan(np.radians(alpha_deg))
    l = sum(t for t, _ in layers)
    c, z = 0.0, 0.0
    for t, E in layers:
        a, b = z, z + t
        if a < l / 2:                                         # head-side cone: s = x
            c += _cone_compliance(a, min(b, l / 2), E, d, D_w, ta)
        if b > l / 2:                                         # nut-side cone: s = l - x
            c += _cone_compliance(l - b, l - max(a, l / 2), E, d, D_w, ta)
        z = b
    return 1.0 / c


# ---------------------------------------------------------------- the joint diagram

def joint_constant(k_b: float, k_m: float) -> float:
    """Fraction of an external separating load that the bolt carries: C = k_b / (k_b + k_m)."""
    return k_b / (k_b + k_m)


def separation_load(F_i: float, k_b: float, k_m: float) -> float:
    """External load at which the member clamp force falls to zero: P_sep = F_i / (1 - C)."""
    return F_i / (1.0 - joint_constant(k_b, k_m))


def joint_loads(F_i: float, P: np.ndarray | float, k_b: float, k_m: float) -> tuple[np.ndarray, np.ndarray]:
    """Bolt tension F_b and member clamp force F_m (compression, as a positive number) for preload
    F_i and external separating load P (applied at the bolt head and nut faces).

    Before separation, compatibility (the bolt stretches as much as the members unload) and
    equilibrium (F_b - F_m = P) give F_b = F_i + C P and F_m = F_i - (1 - C) P. After separation
    the bolt carries all of P and F_m = 0."""
    P = np.asarray(P, dtype=float)
    C = joint_constant(k_b, k_m)
    Fb = F_i + C * P
    Fm = F_i - (1.0 - C) * P
    sep = Fm <= 0.0
    return np.where(sep, P, Fb), np.where(sep, 0.0, Fm)


def thermal_preload_change(dT: float, alpha_b: float, layers_alpha: list[tuple[float, float]],
                           k_b: float, k_m: float) -> float:
    """Change of preload when the joint changes temperature by dT. The members (thickness, alpha)
    try to grow by sum(alpha_i t_i) dT, the bolt by alpha_b l dT; the difference is taken up by the
    bolt and members in series:  dF = (sum(alpha_i t_i) - alpha_b l) dT k_b k_m / (k_b + k_m).
    Positive means the preload rises. Stiffness changes with temperature are ignored."""
    l = sum(t for t, _ in layers_alpha)
    mismatch = (sum(a * t for t, a in layers_alpha) - alpha_b * l) * dT
    return mismatch * k_b * k_m / (k_b + k_m)


def preload_from_torque(T: float, K: float, d: float) -> float:
    """Short-form torque-tension relation F_i = T / (K d), nut factor K (about 0.2 dry steel,
    with a scatter of tens of per cent)."""
    return T / (K * d)


@dataclass
class Joint:
    """A single-bolt joint: bolt (d, pitch, material, shank and thread lengths in the grip),
    members as (thickness, material) from the head side, bearing-face diameter D_w, preload F_i."""
    d: float
    pitch: float
    bolt_material: str
    l_shank: float
    l_thread: float
    layers: list[tuple[float, str]]
    D_w: float
    F_i: float

    @property
    def grip(self) -> float:
        return sum(t for t, _ in self.layers)

    @property
    def k_b(self) -> float:
        return bolt_stiffness(self.d, self.pitch, self.l_shank, self.l_thread, MATERIALS[self.bolt_material]["E"])

    @property
    def k_m(self) -> float:
        return member_stiffness(self.d, self.D_w, [(t, MATERIALS[m]["E"]) for t, m in self.layers])

    @property
    def A_t(self) -> float:
        return tensile_stress_area(self.d, self.pitch)

    def preload_at(self, dT: float) -> float:
        """Preload after a uniform temperature change dT from assembly."""
        la = [(t, MATERIALS[m]["alpha"]) for t, m in self.layers]
        return self.F_i + thermal_preload_change(dT, MATERIALS[self.bolt_material]["alpha"], la, self.k_b, self.k_m)


LESSON_JOINT = Joint(d=8.0, pitch=1.25, bolt_material="Ti-6Al-4V", l_shank=14.0, l_thread=6.0,
                     layers=[(10.0, "7075-T6"), (10.0, "7075-T6")], D_w=12.0, F_i=20e3)
P_LIMIT = 12e3                    # N, limit separating load per bolt (the lesson's load case)
FS_ULT = 1.5                      # ultimate factor of safety
DT_COLD, DT_HOT = -75.0, 50.0     # K: assembled at 20 C, cold soak -55 C, hot soak +70 C


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save bolted_joint.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)

    # ---- bars
    al, st = MATERIALS["2024-T3"], MATERIALS["4130"]
    print(f"2024-T3 rod 1.2 m, 100 mm^2, 15 kN: delta = {axial_deformation(15e3, 1200.0, 100.0, al['E']):.3f} mm")
    print(f"constrained heating by 50 K: 2024-T3 {constrained_thermal_stress(al['E'], al['alpha'], 50.0):.1f} MPa, "
          f"4130 {constrained_thermal_stress(st['E'], st['alpha'], 50.0):.1f} MPa, "
          f"Ti-6Al-4V {constrained_thermal_stress(113.8e3, 8.6e-6, 50.0):.1f} MPa")
    # fixed-fixed: 4130 steel 300 mm x 80 mm^2 joined to 2024-T3 400 mm x 150 mm^2, 10 kN at the junction, +40 K
    args = (10e3, 300.0, st["E"] * 80.0, st["alpha"], 400.0, al["E"] * 150.0, al["alpha"])
    for dT in (0.0, 40.0):
        u, N1, N2 = fixed_fixed_two_segment(*args, dT=dT)
        x = np.array([0.0, 300.0, 700.0])
        uf, Nf, Rf = bar_fe_solve(x, None, [st["E"] * 80.0, al["E"] * 150.0], [0.0, 10e3, 0.0], [0, 2],
                                  eps0=[st["alpha"] * dT, al["alpha"] * dT])
        print(f"steel-aluminium bar, 10 kN at junction, dT = {dT:+.0f} K: u = {u:.4f} mm, N1 = {N1:.0f} N "
              f"({N1 / 80:.1f} MPa), N2 = {N2:.0f} N ({N2 / 150:.1f} MPa); FE |diff| = "
              f"{max(abs(uf[1] - u), abs(Nf[0] - N1), abs(Nf[1] - N2)):.1e}, reactions + load = {Rf.sum() + 10e3:.1e} N")
    # tapered bar convergence
    P, L, A0, A1 = 20e3, 500.0, 200.0, 50.0
    exact = tapered_bar_deformation(P, L, A0, A1, al["E"])
    print(f"\ntapered 2024-T3 bar 200 -> 50 mm^2, 500 mm, 20 kN: exact delta = {exact:.5f} mm")
    print("elements   FE delta (mm)   rel. error")
    for m in (2, 4, 8, 16, 32):
        x = np.linspace(0, L, m + 1)
        xm = 0.5 * (x[1:] + x[:-1])
        u, _, _ = bar_fe_solve(x, None, al["E"] * (A0 + (A1 - A0) * xm / L), np.r_[np.zeros(m), P], [0])
        print(f"{m:8d}   {u[-1]:.5f}        {abs(u[-1] - exact) / exact:.2e}")

    # ---- the lesson joint
    J = LESSON_JOINT
    kb, km = J.k_b, J.k_m
    C = joint_constant(kb, km)
    print(f"\njoint: M{J.d:.0f}x{J.pitch} {J.bolt_material} bolt, grip {J.grip:.0f} mm of 7075-T6, D_w = {J.D_w:.0f} mm")
    print(f"  A_t = {J.A_t:.2f} mm^2, shank area {np.pi * J.d ** 2 / 4:.2f} mm^2")
    print(f"  k_b = {kb / 1e3:.1f} kN/mm, k_m = {km / 1e3:.1f} kN/mm, C = {C:.3f}")
    print(f"  preload {J.F_i / 1e3:.1f} kN = {J.F_i / J.A_t:.0f} MPa in the thread ({J.F_i / J.A_t / 880:.2f} Fty); "
          f"torque at K = 0.2: {0.2 * J.F_i * J.d / 1e3:.1f} N m")
    print("  case             preload  P_sep   at limit (12 kN):           at ultimate (18 kN):")
    print("                   (kN)     (kN)    F_b (kN)  stress  clamp    F_b (kN)  stress  clamp (kN)")
    for label, dT in (("assembly 20 C", 0.0), ("cold -55 C", DT_COLD), ("hot +70 C", DT_HOT)):
        Fi = J.preload_at(dT)
        row = f"  {label:15s}  {Fi / 1e3:6.2f}  {separation_load(Fi, kb, km) / 1e3:6.2f}  "
        for P in (P_LIMIT, FS_ULT * P_LIMIT):
            Fb, Fm = joint_loads(Fi, P, kb, km)
            row += f"  {float(Fb) / 1e3:7.2f}  {float(Fb) / J.A_t:6.0f}  {float(Fm) / 1e3:6.2f}  "
        print(row)
    # independent check: the joint as two parallel bar elements with a misfit
    l = J.grip
    misfit = J.F_i * (1 / kb + 1 / km)
    u, N, R = bar_fe_solve(np.array([0.0, l]), [[0, 1], [0, 1]], [kb * l, km * l], [0.0, P_LIMIT], [0],
                           eps0=[-misfit / l, 0.0])
    Fb, Fm = joint_loads(J.F_i, P_LIMIT, kb, km)
    print(f"  FE check (two parallel elements, misfit {misfit * 1e3:.2f} um): F_b = {N[0]:.1f} N vs {float(Fb):.1f} N, "
          f"F_m = {-N[1]:.1f} N vs {float(Fm):.1f} N")
    # steel bolt for comparison
    Js = Joint(**{**J.__dict__, "bolt_material": "4130"})
    print(f"  same joint, 4130 bolt: k_b = {Js.k_b / 1e3:.1f} kN/mm, C = {joint_constant(Js.k_b, Js.k_m):.3f}, "
          f"cold preload {Js.preload_at(DT_COLD) / 1e3:.2f} kN, hot {Js.preload_at(DT_HOT) / 1e3:.2f} kN")

    print("  cone half-angle sensitivity: " + ", ".join(
        f"{a:.0f} deg: k_m {member_stiffness(J.d, J.D_w, [(t, MATERIALS[mm]['E']) for t, mm in J.layers], a) / 1e3:.1f} kN/mm, "
        f"C {joint_constant(kb, member_stiffness(J.d, J.D_w, [(t, MATERIALS[mm]['E']) for t, mm in J.layers], a)):.3f}"
        for a in (25.0, 30.0, 35.0)))

    # ---- experiment: torque scatter Monte Carlo
    rng = np.random.default_rng(507)
    T = 0.2 * J.F_i * J.d
    K = rng.normal(0.20, 0.025, 100_000)
    Fi = preload_from_torque(T, K, J.d)
    dF = J.preload_at(DT_COLD) - J.F_i
    Psep = separation_load(Fi + dF, kb, km)
    print(f"\nMonte Carlo, nut factor K ~ N(0.20, 0.025), 100000 joints torqued to {T / 1e3:.1f} N m (seed 507):")
    print(f"  preload 1st-99th percentile {np.percentile(Fi, 1) / 1e3:.1f}-{np.percentile(Fi, 99) / 1e3:.1f} kN; "
          f"cold P_sep below limit in {np.mean(Psep < P_LIMIT) * 100:.2f} %, "
          f"below ultimate in {np.mean(Psep < FS_ULT * P_LIMIT) * 100:.2f} %")
    Fb_hot, _ = joint_loads(Fi + J.preload_at(DT_HOT) - J.F_i, P_LIMIT, kb, km)
    Fb_hot_u, _ = joint_loads(Fi + J.preload_at(DT_HOT) - J.F_i, FS_ULT * P_LIMIT, kb, km)
    print(f"  hot: bolt stress at limit above Fty (880 MPa) in {np.mean(Fb_hot / J.A_t > 880) * 100:.2f} %, "
          f"at ultimate above Ftu (950 MPa) in {np.mean(Fb_hot_u / J.A_t > 950) * 100:.2f} %")
    for sd in (0.0125, 0.05):
        Fi_sd = preload_from_torque(T, rng.normal(0.20, sd, 100_000), J.d)
        Fb_sd, _ = joint_loads(Fi_sd + J.preload_at(DT_HOT) - J.F_i, P_LIMIT, kb, km)
        print(f"  nut-factor sd {sd}: cold P_sep below ultimate in "
              f"{np.mean(separation_load(Fi_sd + dF, kb, km) < FS_ULT * P_LIMIT) * 100:.2f} %, "
              f"hot bolt stress at limit above Fty in {np.mean(Fb_sd / J.A_t > 880) * 100:.2f} %")

    # ---- figure: joint diagram
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    db = np.linspace(0, J.F_i / kb * 1.6, 50)
    dm = np.linspace(0, J.F_i / km, 50)
    ax[0].plot(db * 1e3, kb * db / 1e3, label="bolt: F = k_b x stretch")
    ax[0].plot((J.F_i / kb + J.F_i / km - dm) * 1e3, km * dm / 1e3, label="members: F = k_m x squeeze")
    Fb, Fm = joint_loads(J.F_i, P_LIMIT, kb, km)
    x_l = (J.F_i + C * P_LIMIT) / kb
    ax[0].vlines(x_l * 1e3, float(Fm) / 1e3, float(Fb) / 1e3, color="k", label=f"P = {P_LIMIT / 1e3:.0f} kN")
    ax[0].set_xlabel("deformation (um)"); ax[0].set_ylabel("force (kN)")
    ax[0].set_title("joint diagram (assembly temperature)"); ax[0].legend(fontsize=8)
    P = np.linspace(0, 40e3, 400)
    for label, dT in (("20 C", 0.0), ("-55 C", DT_COLD), ("+70 C", DT_HOT)):
        Fb, Fm = joint_loads(J.preload_at(dT), P, kb, km)
        line, = ax[1].plot(P / 1e3, Fb / 1e3, label=f"bolt {label}")
        ax[1].plot(P / 1e3, Fm / 1e3, "--", color=line.get_color(), label=f"clamp {label}")
    ax[1].axvline(P_LIMIT / 1e3, color="k", lw=0.8)
    ax[1].set_xlabel("external load P (kN)"); ax[1].set_ylabel("force (kN)"); ax[1].legend(fontsize=7)
    ax[1].set_title("load sharing and separation")
    fig.tight_layout()
    fig.savefig(out_dir / "bolted_joint.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'bolted_joint.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
