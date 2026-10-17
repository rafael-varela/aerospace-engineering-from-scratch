"""Public tests for P05.10 Bending Stress in Beams. They assert mechanics, not magic numbers: the
section properties obey the parallel-axis theorem and agree with independent sums; the flexure stress
is linear, zero on the centroidal axis and has zero axial resultant and the right moment resultant;
the unsymmetric formula reduces to the flexure formula when I_xy = 0 and agrees with an independent
fibre (plane-strain) model when it does not; the elliptic-lift shear and moment satisfy equilibrium
with the load and match numerical integration; an independent 1-D finite-element cantilever gives the
same bending stress through its curvature, converging at second order; the two-flange idealisation
and the margins of safety behave as the lesson derives.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)


def spar_rects():
    return impl.i_section(0.180, 0.060, 0.010, 0.003)


# ---------------------------------------------------------------- section properties

def test_rectangle_properties_match_closed_forms():
    b, h = 0.04, 0.10
    p = impl.section_properties(np.array([[b, h, 0.3, -0.2]]))
    assert p["A"] == pytest.approx(b * h)
    assert p["xc"] == pytest.approx(0.3) and p["yc"] == pytest.approx(-0.2)
    assert p["Ixx"] == pytest.approx(b * h ** 3 / 12)
    assert p["Iyy"] == pytest.approx(h * b ** 3 / 12)
    assert abs(p["Ixy"]) < 1e-20


def test_i_section_is_symmetric_and_obeys_parallel_axis_theorem():
    h, bf, tf, tw = 0.180, 0.060, 0.010, 0.003
    p = impl.section_properties(impl.i_section(h, bf, tf, tw))
    assert abs(p["xc"]) < 1e-15 and abs(p["yc"]) < 1e-15
    assert abs(p["Ixy"]) < 1e-18
    # independent: outer box minus the two side cut-outs
    I_box = bf * h ** 3 / 12 - (bf - tw) * (h - 2 * tf) ** 3 / 12
    assert p["Ixx"] == pytest.approx(I_box, rel=1e-12)
    assert p["A"] == pytest.approx(2 * bf * tf + tw * (h - 2 * tf))
    # moving the whole section does not change centroidal properties
    r = impl.i_section(h, bf, tf, tw).copy()
    r[:, 2] += 0.7
    r[:, 3] -= 0.4
    q = impl.section_properties(r)
    for k in ("A", "Ixx", "Iyy", "Ixy"):
        assert q[k] == pytest.approx(p[k], rel=1e-9, abs=1e-18)
    assert q["xc"] == pytest.approx(0.7) and q["yc"] == pytest.approx(-0.4)


def test_z_section_has_centroid_at_origin_and_a_product_of_inertia():
    p = impl.section_properties(impl.z_section(0.040, 0.020, 0.002))
    assert abs(p["xc"]) < 1e-15 and abs(p["yc"]) < 1e-15
    assert p["Ixy"] > 0.1 * np.sqrt(p["Ixx"] * p["Iyy"])        # strongly coupled
    assert p["Ixx"] * p["Iyy"] - p["Ixy"] ** 2 > 0                 # positive definite
    # independent: fibre sums
    x, y, dA = impl.fiber_mesh(impl.z_section(0.040, 0.020, 0.002), 60)
    assert np.sum(x * y * dA) == pytest.approx(p["Ixy"], rel=1e-3)
    assert np.sum(y * y * dA) == pytest.approx(p["Ixx"], rel=1e-3)


# ---------------------------------------------------------------- flexure

def test_flexure_stress_is_linear_with_zero_axial_resultant_and_correct_moment():
    rects = spar_rects()
    p = impl.section_properties(rects)
    M = 3.0e4
    x, y, dA = impl.fiber_mesh(rects, 60)
    sig = impl.flexure_stress(M, y - p["yc"], p["Ixx"])
    assert abs(np.sum(sig * dA)) < 1e-9 * M / 0.09                 # N = 0: neutral axis at centroid
    assert -np.sum(sig * y * dA) == pytest.approx(M, rel=1e-3)     # moment resultant
    assert impl.flexure_stress(M, 0.0, p["Ixx"]) == 0.0
    assert impl.flexure_stress(M, 0.09, p["Ixx"]) < 0               # positive M compresses the top
    npt.assert_allclose(impl.flexure_stress(2 * M, 0.05, p["Ixx"]), 2 * impl.flexure_stress(M, 0.05, p["Ixx"]))


def test_general_formula_reduces_to_flexure_formula_without_product_of_inertia():
    rng = np.random.default_rng(510)
    x, y = rng.uniform(-0.05, 0.05, 20), rng.uniform(-0.09, 0.09, 20)
    Ixx, Iyy = 9.7e-6, 3.6e-7
    npt.assert_allclose(impl.bending_stress_general(2.5e4, 0.0, x, y, Ixx, Iyy, 0.0),
                        impl.flexure_stress(2.5e4, y, Ixx), rtol=1e-12)
    # pure lateral moment: stress varies with x only
    s = impl.bending_stress_general(0.0, 500.0, x, y, Ixx, Iyy, 0.0)
    npt.assert_allclose(s, -500.0 * x / Iyy, rtol=1e-12)


def test_general_formula_matches_independent_fibre_model_on_z_section():
    rects = impl.z_section(0.040, 0.020, 0.002)
    p = impl.section_properties(rects)
    pts_x = np.array([0.019, -0.019, 0.0, 0.009, -0.001])
    pts_y = np.array([0.019, -0.019, 0.0, 0.019, 0.010])
    for M, Mh in ((100.0, 0.0), (0.0, 40.0), (80.0, -30.0)):
        s_gen = impl.bending_stress_general(M, Mh, pts_x, pts_y, p["Ixx"], p["Iyy"], p["Ixy"])
        s_fib = impl.fiber_stress(rects, M, Mh, pts_x, pts_y, 70e9, 40)
        scale = np.max(np.abs(s_gen))
        npt.assert_allclose(s_fib, s_gen, atol=2e-3 * scale)       # fibre midpoint sums: O(1/n^2)


def test_z_section_stress_is_zero_on_the_neutral_axis_and_differs_from_naive():
    p = impl.section_properties(impl.z_section(0.040, 0.020, 0.002))
    M = 100.0
    a = impl.neutral_axis_angle(M, 0.0, p["Ixx"], p["Iyy"], p["Ixy"])
    r = np.linspace(-0.02, 0.02, 9)
    s_na = impl.bending_stress_general(M, 0.0, r * np.cos(a), r * np.sin(a), p["Ixx"], p["Iyy"], p["Ixy"])
    assert np.max(np.abs(s_na)) < 1e-6 * M / p["Ixx"] * 0.02
    assert abs(np.degrees(a)) > 10                                  # the neutral axis is not horizontal
    # at the top-flange tip the naive formula gets even the sign wrong
    s_gen = impl.bending_stress_general(M, 0.0, 0.019, 0.019, p["Ixx"], p["Iyy"], p["Ixy"])
    assert s_gen > 0 > impl.flexure_stress(M, 0.019, p["Ixx"])
    # symmetric section: neutral axis horizontal
    assert impl.neutral_axis_angle(M, 0.0, 9.7e-6, 3.6e-7, 0.0) == pytest.approx(0.0, abs=1e-12)


def test_fibre_model_reproduces_flexure_and_finds_the_centroid():
    rects = spar_rects()
    E = 71.7e9
    a, b, c = impl.fiber_strain_plane(rects, 3.0e4, 0.0, E, 0.0, 40)
    assert abs(b) < 1e-12 * abs(c)                                   # no lateral curvature
    assert abs(a / c) < 1e-9                                         # strain zero at y = 0, the centroid
    I = impl.section_properties(rects)["Ixx"]
    assert -E * c == pytest.approx(3.0e4 / I, rel=1e-3)              # curvature = M / (E I)
    # an axial force alone gives uniform stress N / A
    A = impl.section_properties(rects)["A"]
    a, b, c = impl.fiber_strain_plane(rects, 0.0, 0.0, E, 1.0e4, 20)
    assert E * a == pytest.approx(1.0e4 / A, rel=1e-9) and abs(c) < 1e-12 * abs(a)


# ---------------------------------------------------------------- wing load

def test_elliptic_load_integrates_to_the_semispan_lift():
    L, s = 13974.5, 5.0
    y = np.linspace(0, s, 20001)
    w = impl.elliptic_lift(y, L, s)
    assert np.trapezoid(w, y) == pytest.approx(L, rel=1e-4)
    assert w[-1] == 0.0 and w[0] == pytest.approx(4 * L / (np.pi * s))


def test_elliptic_shear_and_moment_satisfy_equilibrium():
    L, s = 13974.5, 5.0
    V0, M0 = impl.elliptic_shear_moment(0.0, L, s)
    assert float(V0) == pytest.approx(L)                              # root shear = all the lift
    assert float(M0) == pytest.approx(L * 4 * s / (3 * np.pi))        # lift acts at 4s/(3 pi)
    Vt, Mt = impl.elliptic_shear_moment(s, L, s)
    assert abs(float(Vt)) < 1e-9 and abs(float(Mt)) < 1e-9            # free tip
    y = np.linspace(0.2, 4.8, 7)
    h = 1e-5
    V, M = impl.elliptic_shear_moment(y, L, s)
    Vp, Mp = impl.elliptic_shear_moment(y + h, L, s)
    Vm, Mm = impl.elliptic_shear_moment(y - h, L, s)
    npt.assert_allclose((Mp - Mm) / (2 * h), -V, rtol=1e-6)           # dM/dy = -V
    npt.assert_allclose((Vp - Vm) / (2 * h), -impl.elliptic_lift(y, L, s), rtol=1e-5)   # dV/dy = -w


def test_numerical_shear_and_moment_match_closed_form_and_converge():
    L, s = 13974.5, 5.0
    errs = []
    for n in (201, 401, 801):
        y = np.linspace(0, s, n)
        V, M = impl.shear_moment_numeric(y, impl.elliptic_lift(y, L, s))
        Vc, Mc = impl.elliptic_shear_moment(y, L, s)
        assert V[-1] == 0.0 and M[-1] == 0.0
        errs.append(abs(M[0] - Mc[0]) / Mc[0])
    assert errs[-1] < 1e-4 and errs[1] < errs[0] and errs[2] < errs[1]
    # uniform load: exact for the trapezoid rule on V, M = w (s-y)^2 / 2
    y = np.linspace(0, 2.0, 11)
    V, M = impl.shear_moment_numeric(y, np.full_like(y, 100.0))
    npt.assert_allclose(V, 100.0 * (2.0 - y), atol=1e-9)
    npt.assert_allclose(M, 50.0 * (2.0 - y) ** 2, rtol=1e-2, atol=1e-9)


# ---------------------------------------------------------------- independent FEM check

def test_fem_is_exact_for_no_load_and_uniform_load():
    L, EI = 2.0, 1.0e5
    r = impl.cantilever_fem(L, EI, 4, lambda x: np.zeros_like(x))
    assert np.allclose(r["w"], 0.0)                                   # no load: no deflection
    # uniform load: Hermite elements reproduce the exact tip deflection w L^4 / (8 EI)
    q = 250.0
    r = impl.cantilever_fem(L, EI, 4, lambda x: np.full_like(x, q))
    assert r["w"][-1] == pytest.approx(q * L ** 4 / (8 * EI), rel=1e-10)
    npt.assert_allclose(r["M_nodes"], q * (L - r["x"]) ** 2 / 2, rtol=1e-9, atol=1e-8)
    assert r["w"][0] == 0.0 and r["theta"][0] == 0.0                  # clamped root


def test_fem_curvature_stress_matches_flexure_formula_at_second_order():
    L, s = 13974.5, 5.0
    E, rects = 71.7e9, spar_rects()
    I = impl.section_properties(rects)["Ixx"]
    c = 0.09
    errs = []
    for n_el in (10, 20, 40):
        r = impl.cantilever_fem(s, E * I, n_el, lambda x: impl.elliptic_lift(x, L, s))
        sig_fem = E * r["kappa_mid"] * c                              # bottom fibre: sigma = -E kappa (-c)
        _, M = impl.elliptic_shear_moment(r["x_mid"], L, s)
        sig_flex = -impl.flexure_stress(M, c, I)
        errs.append(np.max(np.abs(sig_fem - sig_flex)) / np.max(sig_flex))
    orders = np.log(np.array(errs[:-1]) / np.array(errs[1:])) / np.log(2)
    assert np.all(orders > 1.7)
    assert errs[-1] < 5e-4
    assert r["w"][-1] > 0                                             # lift bends the wing up


# ---------------------------------------------------------------- caps and margins

def test_two_flange_idealisation_brackets_the_flexure_formula():
    sp = impl.Spar()
    M = 2.9655e4
    s_two = impl.two_flange_stress(M, sp.A_cap, sp.h_c)
    s_cap_centroid = M * (sp.h_c / 2) / sp.I
    s_outer = M * sp.c / sp.I
    assert s_cap_centroid < s_two                                      # ignores the web: conservative
    assert s_two == pytest.approx(s_outer, rel=0.1)                    # but within 10 % of the extreme fibre
    # thin caps far apart, no web: the idealisation becomes exact at the cap centroid
    thin = impl.i_section(1.0, 0.05, 0.0005, 1e-12)
    I = impl.section_properties(thin)["Ixx"]
    assert impl.two_flange_stress(M, 0.05 * 0.0005, 0.9995) == pytest.approx(M * 0.49975 / I, rel=1e-6)


def test_cap_area_gives_zero_margin_and_margin_sign():
    M, h_c, Ftu = 2.9655e4, 0.17, 572e6
    A = impl.cap_area_required(M, h_c, Ftu, 1.5)
    assert impl.margin_of_safety(Ftu, impl.two_flange_stress(M, A, h_c), 1.5) == pytest.approx(0.0, abs=1e-12)
    assert impl.cap_area_required(2 * M, h_c, Ftu, 1.5) == pytest.approx(2 * A)     # linear in load
    assert impl.cap_area_required(M, 2 * h_c, Ftu, 1.5) == pytest.approx(A / 2)     # deeper spar, half the cap
    assert impl.margin_of_safety(500e6, 250e6, 1.0) == pytest.approx(1.0)
    assert impl.margin_of_safety(500e6, -400e6, 1.5) < 0                             # compression counts too


def test_lesson_spar_has_positive_margins_at_limit_and_ultimate():
    sp = impl.Spar()
    L = impl.wing_semispan_lift(750.0, 3.8)
    _, M0 = impl.elliptic_shear_moment(0.0, L, 5.0)
    sig = impl.flexure_stress(float(M0), -sp.c, sp.I)
    assert 100e6 < sig < sp.material.F_ty                               # realistic, below yield
    assert impl.margin_of_safety(sp.material.F_tu, sig, 1.5) > 0
    assert impl.wing_semispan_lift(750.0, 1.0) == pytest.approx(750 * 9.80665 / 2)


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "wing_spar_stress.png").stat().st_size > 10_000
