"""Public tests for P10.11 Vortex Lattice Method. They assert the method, not magic numbers: the
lattice tiles the planform with upward unit normals and its points where the 1/4-3/4 rule puts them;
a strip in the middle of a very long wing recovers the 2-D slope 2 pi; the span load is symmetric;
lift is linear in alpha and induced drag quadratic in lift; Hough's quarter-panel tip inset gives
second-order lattice convergence and the plain uniform lattice first order; the VLM agrees with the
lifting line (a compact P10.10 solver) for unswept wings with a difference that shrinks like 1/AR;
swept-wing slopes follow the Helmbold-DATCOM trend; moments transfer between reference points as
statics requires and the aerodynamic centre is where dCm/dCL = 0; dihedral and taper act the right way.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

B, C = 10.0, 1.2
A5 = np.radians(5.0)


def test_lattice_geometry():
    lat = impl.make_lattice(B, 2.0, taper=0.5, sweep_c4=np.radians(20.0), dihedral=np.radians(6.0), nx=3, ny=8,
                            spacing="uniform")
    assert len(lat["A"]) == 24
    assert np.sum(lat["area"]) == pytest.approx(lat["S"] / np.cos(np.radians(6.0)), rel=1e-12)   # tilted panels
    npt.assert_allclose(np.linalg.norm(lat["normal"], axis=1), 1.0)
    assert np.all(lat["normal"][:, 2] > 0.9)
    # root strip of a rectangular, unswept lattice: bound legs at 1/4 and control points at 3/4 of each panel
    r = impl.make_lattice(B, C, nx=2, ny=4, spacing="uniform")
    npt.assert_allclose(r["A"][:2, 0], [0.25 * 0.6, 0.6 + 0.25 * 0.6])
    npt.assert_allclose(r["cp"][:2, 0], [0.75 * 0.6, 0.6 + 0.75 * 0.6])
    h = impl.make_lattice(B, C, nx=1, ny=10, spacing="hough")
    d = B / 10.5
    assert h["y_edges"][0] == pytest.approx(-B / 2 + d / 4) and np.allclose(np.diff(h["y_edges"]), d)


def test_middle_of_a_very_long_wing_is_two_dimensional():
    lat = impl.make_lattice(2000.0, 1.0, nx=1, ny=200, spacing="hough")
    r = impl.analyse(lat, np.radians(2.0))
    mid = np.argmin(np.abs(r["y_strip"]))
    assert r["cl_strip"][mid] / (2 * np.pi * np.sin(np.radians(2.0))) == pytest.approx(1.0, abs=0.01)


def test_unit_horseshoe_far_wake_is_two_dimensional():
    s = 3.0
    A, Bp = np.array([[0.0, -s / 2, 0.0]]), np.array([[0.0, s / 2, 0.0]])
    v = impl.horseshoe_unit_velocity(np.array([[1e5, 0.0, 0.0]]), A, Bp)[0, 0]
    assert -v[2] == pytest.approx(2 / (np.pi * s), rel=1e-8)                 # two infinite trailing vortices
    v0 = impl.horseshoe_unit_velocity(np.array([[0.0, 0.0, 0.0]]), A, Bp)[0, 0]
    assert -v0[2] == pytest.approx(1 / (np.pi * s), rel=1e-12)                # half of it at the bound leg


def test_two_dimensional_limit_puts_the_aerodynamic_centre_at_the_quarter_chord():
    r = impl.lift_slope_and_ac(impl.make_lattice(400.0, 1.0, nx=4, ny=80))
    assert r["x_ac_mac"] == pytest.approx(0.25, abs=0.003)
    assert r["CL_alpha"] / (2 * np.pi) > 0.95


def test_symmetric_load_linear_lift_quadratic_drag():
    lat = impl.make_lattice(B, 1.6, taper=0.5, sweep_c4=np.radians(25.0), nx=3, ny=20)
    r1, r2 = impl.analyse(lat, np.radians(2.0)), impl.analyse(lat, np.radians(4.0))
    npt.assert_allclose(r1["G_strip"], r1["G_strip"][::-1], rtol=1e-10)
    assert np.all(r1["G_strip"] > 0)
    assert r2["CL"] / r1["CL"] == pytest.approx(np.sin(np.radians(4.0)) / np.sin(np.radians(2.0)), rel=1e-9)
    assert r2["CDi"] / r1["CDi"] == pytest.approx((r2["CL"] / r1["CL"]) ** 2, rel=1e-9)
    assert 0.9 < r1["e"] <= 1.0
    assert impl.analyse(lat, 0.0)["CL"] == pytest.approx(0.0, abs=1e-14)


def test_zero_lift_angle_shift():
    lat = impl.make_lattice(B, C, nx=2, ny=16)
    aL0 = np.radians(-2.0)
    assert impl.analyse(lat, aL0, alpha_L0=aL0)["CL"] == pytest.approx(0.0, abs=1e-14)
    assert impl.analyse(lat, 0.0, alpha_L0=aL0)["CL"] == pytest.approx(impl.analyse(lat, -aL0)["CL"], rel=1e-12)


def test_observed_order_helper():
    f = [1 + 0.3 * h ** 2 for h in (0.4, 0.2, 0.1)]
    assert impl.observed_order(f) == pytest.approx(2.0)


def test_lattice_convergence_orders():
    cl_h = [impl.analyse(impl.make_lattice(B, C, nx=2, ny=n, spacing="hough"), A5)["CL"] for n in (10, 20, 40)]
    cl_u = [impl.analyse(impl.make_lattice(B, C, nx=2, ny=n, spacing="uniform"), A5)["CL"] for n in (10, 20, 40)]
    assert 1.7 < impl.observed_order(cl_h) < 2.3                              # quarter-panel inset: order 2
    assert 0.6 < impl.observed_order(cl_u) < 1.2                              # tip vortex at the tip: order 1
    assert abs(cl_h[-1] - cl_h[-2]) < abs(cl_u[-1] - cl_u[-2]) / 5
    cl_x = [impl.analyse(impl.make_lattice(B, C, nx=n, ny=20), A5)["CL"] for n in (2, 4, 8)]
    assert 1.7 < impl.observed_order(cl_x) < 2.5                              # chordwise refinement


def test_vlm_agrees_with_lifting_line_within_a_one_over_AR_difference():
    diffs = []
    for AR in (8.33, 20.0):
        c = B / AR
        vlm = impl.lift_slope_and_ac(impl.make_lattice(B, c, nx=2, ny=60))["CL_alpha"]
        ll = impl.lifting_line_CL_e(B, lambda y: np.full_like(y, c), 0.1, N=100)[0] / 0.1
        diffs.append((vlm / ll - 1) * AR)
        assert vlm < ll                                                       # chordwise extent at the tips
    assert all(-0.5 < d < -0.35 for d in diffs)                               # difference ~ -0.43/AR
    assert abs(diffs[1] / 20.0) < 0.025


def test_swept_wings_follow_the_helmbold_datcom_trend():
    c = B / 6.0
    v = {L: impl.lift_slope_and_ac(impl.make_lattice(B, c, sweep_c4=np.radians(L), nx=4, ny=30))["CL_alpha"]
         for L in (0.0, 30.0, 45.0)}
    h = {L: impl.helmbold_datcom(6.0, np.radians(L)) for L in v}
    assert v[0.0] > v[30.0] > v[45.0]
    for L in (30.0, 45.0):
        assert v[L] / v[0.0] == pytest.approx(h[L] / h[0.0], rel=0.03)
        assert 0.90 < v[L] / h[L] < 0.97


def test_moment_transfer_and_aerodynamic_centre():
    lat = impl.make_lattice(B, 1.8, taper=0.45, sweep_c4=np.radians(25.0), nx=3, ny=20)
    cref = impl.mac(1.8, 0.45)
    g = impl.solve_vlm(lat, A5)
    r1 = impl.vlm_forces(lat, g, A5, 0.0, cref)
    r2 = impl.vlm_forces(lat, g, A5, 1.0, cref)
    assert r2["Cm"] - r1["Cm"] == pytest.approx(r1["CL"] * np.cos(A5) * 1.0 / cref, rel=1e-10)   # M2 = M1 + dx Fz
    ac = impl.lift_slope_and_ac(lat)
    m1 = impl.analyse(lat, np.radians(2.0), x_ref=ac["x_ac"], c_ref=cref)["Cm"]
    m2 = impl.analyse(lat, np.radians(6.0), x_ref=ac["x_ac"], c_ref=cref)["Cm"]
    assert m1 == pytest.approx(m2, abs=2e-4)                                  # dCm/dCL = 0 about the a.c.
    rect = impl.lift_slope_and_ac(impl.make_lattice(B, C, nx=4, ny=30))
    assert 0.22 < rect["x_ac_mac"] < 0.26                                     # near the quarter chord
    sw = impl.lift_slope_and_ac(impl.make_lattice(B, C, sweep_c4=np.radians(30.0), nx=4, ny=30))
    assert sw["x_ac"] > rect["x_ac"] + 0.5                                    # sweep moves it aft


def test_dihedral_and_taper_act_the_right_way():
    flat = impl.lift_slope_and_ac(impl.make_lattice(B, C, nx=2, ny=30))["CL_alpha"]
    dih = impl.lift_slope_and_ac(impl.make_lattice(B, C, dihedral=np.radians(10.0), nx=2, ny=30))["CL_alpha"]
    assert np.cos(np.radians(10.0)) < dih / flat < 1.0                       # strip theory: cos(dihedral)
    e_rect = impl.analyse(impl.make_lattice(B, C, nx=2, ny=30), A5)["e"]
    e_tap = impl.analyse(impl.make_lattice(B, 2 * 12.0 / (B * 1.4), taper=0.4, nx=2, ny=30), A5)["e"]
    assert e_rect < e_tap <= 1.0


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "vlm.png").stat().st_size > 10_000
