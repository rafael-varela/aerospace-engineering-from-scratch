"""Public tests for P10.16 Fuselage & Body Aerodynamics. They assert physics, not magic numbers:
slender-body loading integrates to N = 2 q alpha S_base (C_N_alpha = 2 on the base area) and to zero
for a closed body; a closed body feels a pure couple 2 q alpha Vol independent of the moment
reference (Munk); the paraboloid's uniform loading puts the centre of pressure at L/2; the exact
spheroid moment (Lamb's k2 - k1) vanishes for a sphere and tends to the slender-body value as the
fineness ratio grows; volume and wetted area match closed forms; drag area has an interior
optimum in fineness; the fuselage moves the trainer's neutral point forward.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import pytest

from conftest import load_impl

impl = load_impl(__file__)

Q, ALPHA = 1500.0, np.radians(3.0)


def cosine_grid(L, n=4001):
    t = np.linspace(0, np.pi, n)
    return L * (1 - np.cos(t)) / 2


# ---------------------------------------------------------------- geometry

def test_volume_and_wetted_area_of_a_sphere_and_spheroid():
    R = 0.8
    x = cosine_grid(2 * R)
    r = impl.ellipsoid_radius(x, 2 * R, R)
    assert impl.body_volume(x, r) == pytest.approx(4 / 3 * np.pi * R ** 3, rel=1e-5)
    assert impl.wetted_area(x, r) == pytest.approx(4 * np.pi * R ** 2, rel=1e-5)
    # prolate spheroid surface area: 2 pi b^2 (1 + (a / (b e)) arcsin e)
    a, b = 3.0, 0.5
    e = np.sqrt(1 - (b / a) ** 2)
    x = cosine_grid(2 * a, 8001)
    r = impl.ellipsoid_radius(x, 2 * a, b)
    assert impl.wetted_area(x, r) == pytest.approx(2 * np.pi * b ** 2 * (1 + a / (b * e) * np.arcsin(e)), rel=1e-5)
    assert impl.section_area(np.array([b]))[0] == pytest.approx(np.pi * b ** 2)


# ---------------------------------------------------------------- slender-body theory

def test_closed_body_has_zero_normal_force_and_a_pure_munk_couple():
    L, R = 6.0, 0.5
    x = cosine_grid(L)
    r = impl.ellipsoid_radius(x, L, R)
    vol = 4 / 3 * np.pi * (L / 2) * R ** 2
    moments = []
    for xr in (-2.0, 0.0, 2.5, 6.0, 9.0):
        N, M = impl.normal_force_and_moment(x, r, ALPHA, Q, xr)
        assert abs(N) < 1e-8 * Q * ALPHA * np.pi * R ** 2
        moments.append(M)
    np.testing.assert_allclose(moments, 2 * Q * ALPHA * vol, rtol=1e-5)          # reference-independent
    assert moments[0] > 0                                                         # nose-up: destabilising


def test_munk_moment_is_linear_in_alpha_and_q():
    L, R = 6.0, 0.5
    x = cosine_grid(L)
    r = impl.ellipsoid_radius(x, L, R)
    _, M1 = impl.normal_force_and_moment(x, r, ALPHA, Q, 1.0)
    _, M2 = impl.normal_force_and_moment(x, r, 2 * ALPHA, Q, 1.0)
    _, M3 = impl.normal_force_and_moment(x, r, ALPHA, 3 * Q, 1.0)
    assert M2 == pytest.approx(2 * M1, rel=1e-12) and M3 == pytest.approx(3 * M1, rel=1e-12)


def test_paraboloid_normal_force_and_centre_of_pressure():
    L, R = 2.0, 0.4
    x = np.linspace(0, L, 4001)
    r = impl.paraboloid_radius(x, L, R)
    dN = impl.slender_body_loading(x, r, ALPHA, Q)
    np.testing.assert_allclose(dN[5:-5], 2 * Q * ALPHA * np.pi * R ** 2 / L, rtol=1e-8)   # uniform loading
    N, M = impl.normal_force_and_moment(x, r, ALPHA, Q, 0.0)
    assert N / (Q * ALPHA * np.pi * R ** 2) == pytest.approx(2.0, rel=1e-6)                 # C_N_alpha = 2
    assert -M / N == pytest.approx(L / 2, rel=1e-6)


def test_numerical_loading_matches_closed_form_for_open_body():
    x = np.linspace(0, 7.0, 7001)
    r = impl.trainer_fuselage_radius(x)
    S_b = np.pi * r[-1] ** 2
    vol = impl.body_volume(x, r)
    for xr in (0.0, 2.3, 7.0):
        N, M = impl.normal_force_and_moment(x, r, ALPHA, Q, xr)
        Nc, Mc = impl.slender_body_closed_form(S_b, vol, 7.0, ALPHA, Q, xr)
        assert N == pytest.approx(Nc, rel=1e-3)
        assert M == pytest.approx(Mc, rel=1e-4)


# ---------------------------------------------------------------- independent check: exact spheroid

def test_lamb_coefficients_sphere_and_limits():
    k1, k2 = impl.lamb_coefficients(1.0)
    assert k1 == pytest.approx(0.5) and k2 == pytest.approx(0.5)
    k1, k2 = impl.lamb_coefficients(1.0001)
    assert k1 == pytest.approx(0.5, abs=1e-3) and k2 == pytest.approx(0.5, abs=1e-3)   # continuous at the sphere
    fs = [1.5, 3.0, 6.0, 12.0, 50.0]
    k = np.array([impl.lamb_coefficients(f) for f in fs])
    assert np.all(np.diff(k[:, 0]) < 0) and np.all(np.diff(k[:, 1]) > 0)
    assert k[-1, 0] < 0.002 and k[-1, 1] > 0.995


def test_exact_spheroid_moment_tends_to_slender_body():
    gaps = []
    for f in (4.0, 8.0, 16.0, 32.0, 64.0):
        L = 6.0
        R = L / (2 * f)
        vol = 4 / 3 * np.pi * (L / 2) * R ** 2
        x = cosine_grid(L)
        _, M_sb = impl.normal_force_and_moment(x, impl.ellipsoid_radius(x, L, R), 1e-4, Q, 0.0)
        M_ex = impl.spheroid_munk_moment(f, vol, 1e-4, Q)
        assert M_ex < M_sb                                                            # slender body over-predicts
        gaps.append(1 - M_ex / M_sb)
    assert np.all(np.diff(gaps) < 0) and gaps[-1] < 0.005
    assert impl.spheroid_munk_moment(1.0, 1.0, 0.3, Q) == pytest.approx(0.0, abs=1e-12)   # sphere: no moment


def test_allen_perkins_reduces_to_potential_at_small_alpha():
    S_b, A_p, S_ref = 0.5, 4.0, 10.0
    a = 1e-3
    assert impl.allen_perkins_cn(a, S_b, A_p, S_ref) == pytest.approx(2 * a * S_b / S_ref, rel=1e-2)
    assert impl.allen_perkins_cn(a, S_b, A_p, S_ref, eta=0.0) == pytest.approx(
        S_b / S_ref * np.sin(2 * a) * np.cos(a / 2), rel=1e-12)
    # the viscous term grows as alpha^2 and dominates at large alpha
    c10 = impl.allen_perkins_cn(np.radians(10), S_b, A_p, S_ref) - impl.allen_perkins_cn(np.radians(10), S_b, A_p, S_ref, eta=0)
    c20 = impl.allen_perkins_cn(np.radians(20), S_b, A_p, S_ref) - impl.allen_perkins_cn(np.radians(20), S_b, A_p, S_ref, eta=0)
    assert c20 / c10 == pytest.approx((np.sin(np.radians(20)) / np.sin(np.radians(10))) ** 2, rel=1e-10)


# ---------------------------------------------------------------- drag

def test_form_factor_and_best_fineness():
    assert impl.body_form_factor(3.0) > impl.body_form_factor(6.0) > 1.0
    assert impl.body_form_factor(100.0) > impl.body_form_factor(30.0)        # long bodies: friction length
    f_opt = impl.best_fineness(4.0, 55.0)
    d_opt = impl.spheroid_drag_area(f_opt, 4.0, 55.0)["D_over_q"]
    for f in (f_opt * 0.6, f_opt * 1.6):
        assert impl.spheroid_drag_area(f, 4.0, 55.0)["D_over_q"] > d_opt
    d = impl.spheroid_drag_area(5.0, 4.0, 55.0)
    assert np.pi / 6 * d["L"] * d["D"] ** 2 == pytest.approx(4.0, rel=1e-10)       # same volume
    assert d["D_over_q"] == pytest.approx(d["Cf"] * d["FF"] * d["S_wet"])


# ---------------------------------------------------------------- trainer

def test_fuselage_is_destabilising_and_moves_neutral_point_forward():
    x = np.linspace(0, 7.0, 7001)
    r = impl.trainer_fuselage_radius(x)
    full = impl.fuselage_pitch_derivatives(x, r, 2.3)
    cut = impl.fuselage_pitch_derivatives(x, r, 2.3, viscous_cut=True)
    vol = impl.body_volume(x, r)
    assert full["Cm_alpha"] > 0 and cut["Cm_alpha"] > 0
    assert full["Cm_alpha"] == pytest.approx(2 * vol / (12.0 * 1.2), rel=0.05)      # nearly closed: Munk
    assert cut["CN_alpha"] > full["CN_alpha"]                                        # cut keeps the nose lift
    base = impl.trainer_neutral_point()
    with_f = impl.trainer_neutral_point(full["Cm_alpha"], full["CN_alpha"])
    assert 0.35 < base["h_n"] < 0.7
    assert with_f["h_n"] < base["h_n"]


def test_neutral_point_bookkeeping():
    base = impl.trainer_neutral_point()
    # a pure moment derivative dCm moves h_n by -dCm / CL_alpha
    m = impl.trainer_neutral_point(0.2, 0.0)
    assert m["h_n"] - base["h_n"] == pytest.approx(-0.2 / base["CL_alpha"], rel=1e-12)
    # without a tail the neutral point is the wing's aerodynamic centre
    no_tail = impl.trainer_neutral_point(eta_t=0.0)
    assert no_tail["h_n"] == pytest.approx(0.25)


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "fuselage_slender_body.png").stat().st_size > 10_000
