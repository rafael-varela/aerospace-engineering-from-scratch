"""Public tests for P10.12 Wing Planform Effects. They assert planform physics and optimiser
behaviour, not magic numbers: an elliptic load has e = 1 and root bending coefficient CL/(3 pi);
the span efficiency of flat tapered wings peaks at a moderate taper; the local lift integrates to
the wing lift; taper moves the stall outboard and washout moves it back, aft sweep moves it
outboard and forward sweep inboard (VLM); the twist basis superposes; the unconstrained twist
optimum recovers the elliptic load (e -> 1) and improves with more twist stations; the linear-twist
optimum equals the Lagrange closed form and is stationary; bending-moment and tip-stall constraints
hold and cost induced drag monotonically; Prandtl's bell is recovered (span x sqrt(3/2), drag x 8/9).

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

B, S = 10.0, 12.0
AR = B ** 2 / S
ZERO = lambda y: np.zeros_like(np.asarray(y, dtype=float))
RECT = lambda y: np.full_like(np.asarray(y, dtype=float), 1.2)
ETA5 = np.array([0.0, 0.25, 0.5, 0.75, 1.0])


def test_elliptic_coefficients_give_e_one_and_the_elliptic_root_moment():
    A = np.zeros(20); A[0] = 0.4 / (np.pi * AR)
    w = impl.wing_from_coefficients(B, S, RECT, A)
    assert w["CL"] == pytest.approx(0.4) and w["e"] == pytest.approx(1.0)
    assert w["CB"] == pytest.approx(0.4 / (3 * np.pi), rel=1e-6)              # half-lift x 4(b/2)/(3 pi)
    A[1] = -A[0] / 3                                                          # bell: tips unloaded
    wb = impl.wing_from_coefficients(B, S, RECT, A)
    assert wb["e"] == pytest.approx(0.75) and wb["CB"] < w["CB"]


def test_span_efficiency_peaks_at_moderate_taper():
    t = np.linspace(0.0, 1.0, 41)
    for ar in (6.0, 8.333, 12.0):
        e = impl.efficiency_vs_taper(ar, t, N=40)
        assert np.all(e < 1.0)
        assert 0.3 <= t[int(np.argmax(e))] <= 0.42
        assert e[-1] < e.max() - 0.03 and e[0] < e[-1]


def test_local_lift_integrates_to_the_wing_lift():
    ch = impl.tapered_chord(B, S, 0.4)
    tw = lambda y: np.radians(-3.0) * np.abs(2 * np.asarray(y) / B)
    th = (np.arange(3000) + 0.5) * np.pi / 3000
    y = -0.5 * B * np.cos(th)
    r = impl.cl_ratio(B, S, ch, tw, 0.6, y)
    assert np.sum(r * ch(y) * 0.5 * B * np.sin(th)) * np.pi / 3000 / S == pytest.approx(1.0, rel=1e-6)


def test_taper_moves_stall_outboard_and_washout_moves_it_back():
    etas = [impl.stall_onset(B, S, impl.tapered_chord(B, S, t), ZERO, 1.0) for t in (1.0, 0.6, 0.35, 0.2)]
    assert etas[0] == pytest.approx(0.0, abs=1e-12)                           # rectangular: root first
    assert np.all(np.diff(etas) > 0)
    ch = impl.tapered_chord(B, S, 0.35)
    wash = [impl.stall_onset(B, S, ch, lambda y, t=t: np.radians(t) * np.abs(2 * np.asarray(y) / B), 1.0)
            for t in (0.0, -2.0, -4.0)]
    assert wash[0] > wash[1] > wash[2]


def test_sweep_moves_the_stall_station():
    eta = impl.sweep_stall_station(B, S, 0.35, [-20.0, 0.0, 30.0])
    assert eta[0] < eta[1] < eta[2]                                           # aft sweep loads the tips


def test_twist_basis_superposes():
    ch = impl.tapered_chord(B, S, 0.5)
    Bm = impl.twist_basis(B, ch, ETA5, N=30)
    x = np.array([0.07, -0.01, 0.02, -0.03, -0.05])
    tw = lambda y: np.interp(np.abs(2 * np.asarray(y) / B), ETA5, np.concatenate([[0.0], x[1:]]))
    npt.assert_allclose(Bm @ x, impl.ll_solve(B, ch, tw, x[0], 30), atol=1e-14)


def test_unconstrained_twist_recovers_the_elliptic_load():
    es = [impl.optimize_twist(B, S, RECT, 0.33, np.linspace(0, 1, K + 1))["e"] for K in (2, 4, 8)]
    assert es[0] < es[1] < es[2] and es[2] > 0.9999
    r = impl.optimize_twist(B, S, RECT, 0.33, ETA5)
    assert r["CL"] == pytest.approx(0.33, rel=1e-8) and r["success"]
    assert np.all(np.diff(r["twist"]) < 0)                                    # washout, growing outboard
    assert r["CB"] == pytest.approx(0.33 / (3 * np.pi), rel=2e-3)             # elliptic root moment


def test_linear_twist_matches_lagrange_and_is_stationary():
    for ch in (RECT, impl.tapered_chord(B, S, 0.6)):
        lag = impl.linear_twist_optimum(B, S, ch, 0.5)
        num = impl.optimize_twist(B, S, ch, 0.5, np.array([0.0, 1.0]))
        assert num["twist"][1] == pytest.approx(lag["tip_twist"], abs=1e-6)
        assert num["e"] == pytest.approx(lag["e"], rel=1e-8)
    # perturb the tip twist, re-trim alpha: induced drag rises on both sides (a minimum)
    a = impl.ll_solve(B, RECT, ZERO, 1.0, 40)
    r = impl.ll_solve(B, RECT, lambda y: np.abs(2 * np.asarray(y) / B), 0.0, 40)
    lag = impl.linear_twist_optimum(B, S, RECT, 0.5)
    def cdi(t):
        al = (0.5 / (np.pi * AR) - t * r[0]) / a[0]
        return impl.wing_from_coefficients(B, S, RECT, al * a + t * r)["CDi"]
    t0 = lag["tip_twist"]
    assert cdi(t0) == pytest.approx(lag["CDi"], rel=1e-10)
    assert cdi(t0 + 0.01) > cdi(t0) and cdi(t0 - 0.01) > cdi(t0)


def test_elliptic_planform_needs_no_twist():
    c0 = 4 * S / (np.pi * B)
    ell = lambda y: c0 * np.sqrt(np.clip(1 - (2 * np.asarray(y) / B) ** 2, 1e-12, None))
    r = impl.optimize_twist(B, S, ell, 0.5, ETA5)
    assert np.max(np.abs(r["twist"])) < 1e-4 and r["e"] > 0.99999


def test_bending_constraint_holds_and_costs_drag_monotonically():
    res = impl.bending_trade(B, S, RECT, 0.33, ETA5, [1.0, 0.95, 0.9])
    cb_ell = 0.33 / (3 * np.pi)
    for r, ratio in zip(res, [1.0, 0.95, 0.9]):
        assert r["CB"] <= ratio * cb_ell * (1 + 1e-7)
        assert r["CL"] == pytest.approx(0.33, rel=1e-8)
    assert res[0]["e"] > res[1]["e"] > res[2]["e"]
    assert res[2]["twist"][-1] < res[0]["twist"][-1]                         # more washout unloads the tips


def test_tip_stall_constraint_moves_stall_inboard_at_a_drag_price():
    ch = impl.tapered_chord(B, S, 0.35)
    free = impl.optimize_twist(B, S, ch, 0.33, ETA5)
    con = impl.optimize_twist(B, S, ch, 0.33, ETA5, tip_margin=0.05, CL_high=1.1, washout_only=True)
    assert con["success"] and con["CDi"] > free["CDi"]
    tw = lambda y: np.interp(np.abs(2 * np.asarray(y) / B), ETA5, con["twist"])
    y = 0.5 * B * np.linspace(0.7, 0.98, 15)
    assert np.all(impl.cl_ratio(B, S, ch, tw, 1.1, y) <= 0.95 + 1e-6)
    assert impl.stall_onset(B, S, ch, tw, 1.1) < 0.7
    assert np.all(np.diff(con["twist"]) <= 1e-9)                             # washout only


def test_prandtl_bell():
    p = impl.prandtl_bell()
    assert p["span_ratio"] == pytest.approx(np.sqrt(1.5), rel=1e-5)              # f is flat there (f' = 0)
    assert p["drag_ratio"] == pytest.approx(8 / 9, rel=1e-6)
    assert p["G3_over_G1"] == pytest.approx(-1 / 3, rel=1e-4)
    assert abs(p["dfds"]) < 1e-4                                              # Lagrange-stationary as well
    assert p["slsqp_span"] == pytest.approx(p["span_ratio"], rel=1e-4)
    assert p["slsqp_drag"] == pytest.approx(p["drag_ratio"], rel=1e-6)
    npt.assert_allclose(p["G"][2:], 0.0, atol=1e-4)                           # no higher harmonics


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "planform.png").stat().st_size > 10_000
