"""Public tests for P10.07 Airfoil Design & Selection. They assert the inverse problem, not magic
numbers: the sine series of a load round-trips; the ideal angle closes the camber line; a uniform
load recovers the NACA a = 1.0 mean line at the observed order of the series truncation; the forward
thin-airfoil analysis (the given P10.04 solver) of every recovered camber line reproduces its target
load, cl and cm; the inversion is linear in the load; thickness is wrapped correctly; the panel
method approaches the thin-airfoil design load as the thickness goes to zero and departs from it
as it grows; the decision matrix is normalised, scale-invariant and never picks a dominated section;
the weight-sensitivity tools agree with brute force.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

UNIFORM = lambda x: np.full_like(x, 0.4)


def roof(x, a=0.6, c=0.5):
    """Roof-top load: uniform to x = a, linear to zero at the trailing edge, total cl = c."""
    return c * np.where(x < a, 1.0, (1 - x) / (1 - a)) / (1 - 0.5 * (1 - a))


# ---------------------------------------------------------------- forward load and its series

def test_glauert_load_integrates_to_thin_airfoil_cl_and_cm():
    rng = np.random.default_rng(3)
    for _ in range(3):
        A0 = rng.uniform(-0.05, 0.05)
        An = rng.uniform(-0.1, 0.1, size=6)
        cl, cm = impl.load_cl_cm(lambda x: impl.glauert_load(A0, An, x))
        assert cl == pytest.approx(2 * np.pi * (A0 + An[0] / 2), abs=2e-4)
        assert cm == pytest.approx(np.pi / 4 * (An[1] - An[0]), abs=2e-4)


def test_load_coefficients_round_trip():
    An = np.array([0.12, -0.03, 0.02, 0.0, -0.01])
    back = impl.load_coefficients(lambda x: impl.glauert_load(0.0, An, x), n_terms=8)
    npt.assert_allclose(back[:5], An, atol=1e-8)
    npt.assert_allclose(back[5:], 0.0, atol=1e-8)


def test_ideal_angle_closes_the_camber_line():
    # a load symmetric about mid-chord (odd n only) needs no angle; two terms give A2/3
    assert impl.ideal_angle(np.array([0.1, 0.0, 0.05])) == pytest.approx(0.0, abs=1e-15)
    assert impl.ideal_angle(np.array([0.1, 0.06])) == pytest.approx(0.06 / 3)
    An = np.array([0.1, 0.04, -0.02, 0.01])
    x = np.array([0.0, 1.0])
    z, _, _ = impl.camber_from_coefficients(An, x)
    npt.assert_allclose(z, 0.0, atol=1e-14)


def test_camber_slope_is_the_derivative_of_the_camber():
    An = np.array([0.1, 0.04, -0.02, 0.01])
    x = np.linspace(0.05, 0.95, 19)
    h = 1e-6
    z, s, _ = impl.camber_from_coefficients(An, x)
    zp, _, _ = impl.camber_from_coefficients(An, x + h)
    zm, _, _ = impl.camber_from_coefficients(An, x - h)
    npt.assert_allclose((zp - zm) / (2 * h), s, atol=1e-7)


# ---------------------------------------------------------------- the classical result: a = 1.0

def test_naca_a1_closed_form():
    x = np.linspace(0.01, 0.99, 99)
    z, s = impl.naca_a1_camber(x, 0.4)
    npt.assert_allclose(z, z[::-1], atol=1e-15)                               # symmetric about mid-chord
    assert impl.naca_a1_camber(np.array([0.5]), 0.4)[0][0] == pytest.approx(0.4 * np.log(2) / (4 * np.pi))
    h = 1e-6
    d = (impl.naca_a1_camber(x + h, 0.4)[0] - impl.naca_a1_camber(x - h, 0.4)[0]) / (2 * h)
    npt.assert_allclose(d, s, rtol=1e-6, atol=1e-9)
    npt.assert_allclose(impl.naca_a1_camber(np.array([0.0, 1.0]), 0.4)[0], 0.0)


def test_uniform_load_recovers_the_a1_mean_line_at_second_order():
    x = 0.5 * (1 - np.cos(np.linspace(0, np.pi, 201)))
    z_ex, _ = impl.naca_a1_camber(x, 0.4)
    errs = []
    for N in (50, 100, 200):
        z, _, ai = impl.camber_from_coefficients(impl.load_coefficients(UNIFORM, N), x)
        errs.append(np.max(np.abs(z - z_ex)))
        assert abs(ai) < 1e-12
    orders = np.log(np.array(errs[:-1]) / np.array(errs[1:])) / np.log(2)
    assert np.all(orders > 1.8) and np.all(orders < 2.2)                       # tail sum ~ 1/N^2
    assert errs[-1] < 2e-6


def test_forward_analysis_of_recovered_camber_reproduces_the_uniform_load():
    d = impl.inverse_design(UNIFORM, n_terms=200)
    f = impl.thin_airfoil(d["slope"], d["alpha_i"], n_terms=60, m=20000)
    assert abs(f["A0"]) < 1e-10                                              # ideal angle: no LE singularity
    assert f["cl"] == pytest.approx(0.4, rel=1e-3)
    assert f["cm_c4"] == pytest.approx(-0.1, rel=1e-3)                        # centre of pressure at mid-chord
    xs = np.linspace(0.05, 0.95, 19)
    npt.assert_allclose(impl.glauert_load(f["A0"], f["An"], xs), 0.4, atol=0.01)


def test_roof_top_load_round_trip():
    d = impl.inverse_design(roof, n_terms=200)
    assert d["cl"] == pytest.approx(0.5, rel=1e-6)
    f = impl.thin_airfoil(d["slope"], d["alpha_i"], n_terms=200, m=20000)
    assert f["cl"] == pytest.approx(0.5, rel=1e-3) and f["cm_c4"] == pytest.approx(d["cm_c4"], abs=1e-3)
    xs = np.linspace(0.02, 0.98, 49)
    npt.assert_allclose(impl.glauert_load(f["A0"], f["An"], xs), roof(xs), atol=0.02)
    assert d["alpha_i"] > 0                                                   # aft-unloaded: needs positive alpha


def test_inversion_is_linear_in_the_load():
    d1 = impl.inverse_design(roof, n_terms=100)
    d2 = impl.inverse_design(lambda x: 2 * roof(x), n_terms=100)
    npt.assert_allclose(d2["z"], 2 * d1["z"], rtol=1e-10, atol=1e-14)
    assert d2["alpha_i"] == pytest.approx(2 * d1["alpha_i"])


def test_camber_from_cl_cm_hits_its_targets():
    for cl, cm in ((0.3, -0.053), (0.6, -0.15), (0.4, 0.0)):
        An, ai = impl.camber_from_cl_cm(cl, cm)
        f = impl.thin_airfoil(lambda x: impl.camber_from_coefficients(An, x)[1], ai)
        assert f["cl"] == pytest.approx(cl, abs=1e-6)
        assert f["cm_c4"] == pytest.approx(cm, abs=1e-6)
        assert abs(f["A0"]) < 1e-9


# ---------------------------------------------------------------- thickness and the panel method

def test_thicken_wraps_thickness_around_the_camber_line():
    cam = lambda x: impl.naca4_camber("2412")(x)
    slo = lambda x: impl.naca4_slope("2412")(x)
    X, Y = impl.thicken(cam, slo, 0.12, 60)
    X2, Y2 = impl.naca4("2412", 60)
    npt.assert_allclose(X, X2, atol=1e-12)                                    # same construction as P10.02
    npt.assert_allclose(Y, Y2, atol=1e-12)
    assert X[0] == pytest.approx(X[-1], abs=1e-12) and Y[0] == pytest.approx(Y[-1], abs=1e-12)   # closed TE
    area2 = np.sum(X[:-1] * Y[1:] - X[1:] * Y[:-1])
    assert area2 < 0                                                           # clockwise ordering
    X0, Y0 = impl.thicken(cam, slo, 0.0, 60)
    npt.assert_allclose(Y0[:61][::-1], Y0[60:], atol=1e-15)                    # zero thickness: one line


def test_panel_approaches_thin_airfoil_as_thickness_vanishes():
    cam = lambda x: impl.naca_a1_camber(x, 0.4)[0]
    slo = lambda x: impl.naca_a1_camber(np.clip(x, 1e-9, 1 - 1e-9), 0.4)[1]
    xe = np.linspace(0.1, 0.9, 41)
    cls, devs = [], []
    for t in (0.02, 0.06, 0.12):
        X, Y = impl.thicken(cam, slo, t, 100)
        cls.append(impl.panel_solve([(X, Y)], 0.0)["cl"])
        devs.append(np.mean(np.abs(impl.panel_load(X, Y, 0.0, xe) - 0.4)))
    assert cls[0] == pytest.approx(0.4, rel=0.02)                              # thin limit
    assert cls[0] < cls[1] < cls[2]                                            # thickness adds lift
    assert devs[0] < 0.015 and devs[0] < devs[1] < devs[2]                     # and distorts the load


def test_panel_load_integrates_to_panel_cl():
    X, Y = impl.naca4("2412", 100)
    xe = 0.5 * (1 - np.cos(np.linspace(0, np.pi, 801)))
    d = impl.panel_load(X, Y, 4.0, xe)
    cl_int = np.trapezoid(d, xe)
    assert cl_int == pytest.approx(impl.panel_solve([(X, Y)], 4.0)["cl"], rel=0.02)


# ---------------------------------------------------------------- selection

def test_section_metrics_symmetry_and_camber_trends():
    m0 = impl.section_metrics(*impl.naca4("0012", 60), impl.naca4_slope("0012"), 0.3, 1.1)
    assert abs(m0["alpha_L0_deg"]) < 1e-12 and abs(m0["cm0"]) < 1e-12 and abs(m0["cl_ideal"]) < 1e-12
    assert m0["t_c"] == pytest.approx(0.12, abs=1e-3)
    m2 = impl.section_metrics(*impl.naca4("2412", 60), impl.naca4_slope("2412"), 0.3, 1.1)
    m4 = impl.section_metrics(*impl.naca4("4412", 60), impl.naca4_slope("4412"), 0.3, 1.1)
    assert m4["alpha_L0_deg"] == pytest.approx(2 * m2["alpha_L0_deg"], rel=1e-6)   # linear in camber
    assert m4["cm0"] == pytest.approx(2 * m2["cm0"], rel=1e-6)
    assert m0["peak_suction"] > m2["peak_suction"] > m4["peak_suction"]       # camber softens the LE peak


def test_decision_matrix_properties():
    rng = np.random.default_rng(11)
    raw = rng.uniform(0, 1, size=(6, 4))
    lib = [True, False, True, False]
    w = np.array([0.4, 0.3, 0.2, 0.1])
    s, tot, rank = impl.decision_matrix(raw, lib, w)
    assert s.min() >= 0 and s.max() <= 1
    npt.assert_allclose(s.max(axis=0), 1.0)
    npt.assert_allclose(s.min(axis=0), 0.0)
    _, tot2, rank2 = impl.decision_matrix(raw, lib, 10 * w)                   # weights' scale irrelevant
    npt.assert_allclose(tot2, tot)
    _, tot3, _ = impl.decision_matrix(3 * raw + 1, lib, w)                    # units of a criterion irrelevant
    npt.assert_allclose(tot3, tot)
    # a candidate dominated on every criterion never wins
    raw[0] = [raw[:, 0].min() - 1, raw[:, 1].max() + 1, raw[:, 2].min() - 1, raw[:, 3].max() + 1]
    for _ in range(20):
        _, t, r = impl.decision_matrix(raw, lib, rng.dirichlet(np.ones(4)))
        assert r[0] != 0 and t[0] == pytest.approx(0.0)


def test_weight_monte_carlo():
    norm = np.array([[1.0, 0.0], [0.0, 1.0], [0.4, 0.4]])
    win = impl.weight_monte_carlo(norm, np.array([0.5, 0.5]), n_samples=20000, seed=1)
    assert win.sum() == pytest.approx(1.0)
    # symmetric problem: the two specialists win equally often (binomial s.d. 0.0035 at n = 20000)
    assert abs(win[0] - win[1]) < 0.02 and win[2] == 0.0
    win2 = impl.weight_monte_carlo(norm, np.array([0.8, 0.2]), concentration=1e6, seed=1)
    assert win2[0] == pytest.approx(1.0)                                      # tight weights: nominal winner


def test_weight_flip_points_match_brute_force():
    rng = np.random.default_rng(5)
    norm = rng.uniform(0, 1, size=(5, 4))
    w0 = np.array([0.25, 0.35, 0.15, 0.25])
    for k in range(4):
        flips = impl.weight_flip_points(norm, w0, k)
        ws = np.linspace(0, 1, 20001)
        rest = np.delete(np.arange(4), k)
        winners = []
        for w in ws:
            wv = np.zeros(4); wv[k] = w; wv[rest] = (1 - w) * w0[rest] / w0[rest].sum()
            winners.append(int(np.argmax(norm @ wv)))
        brute = [(ws[i], winners[i]) for i in range(1, len(ws)) if winners[i] != winners[i - 1]]
        assert len(brute) == len(flips)
        for (wb, jb), (wf, jf) in zip(brute, flips):
            assert jb == jf and abs(wb - wf) < 1e-4


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "airfoil_design.png").stat().st_size > 10_000
