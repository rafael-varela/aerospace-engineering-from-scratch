"""Public tests for P10.09 Finite Wings & Induced Drag. They assert vortex physics, not magic
numbers: the finite-segment Biot-Savart law reproduces the exact symmetric-segment result and tends
to Gamma/(2 pi h); the right-hand rule; a semi-infinite line is the limit of a long segment and gives
half the infinite line at its end; Stokes and Kelvin/Helmholtz (circulation round each leg equals its
strength, round the whole wake zero, a divergence-free field); a single horseshoe's downwash doubles
from the wing to the Trefftz plane; an elliptic load built from horseshoes has uniform downwash
Gamma_0/(2b), e = 1, near-field drag equal to Trefftz-plane drag, and converges at second order; a
two-term load matches Glauert's e = 1/(1 + delta); the parabolic polar's optima are where they must be.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)


# ---------------------------------------------------------------- Biot-Savart

def test_symmetric_segment_matches_exact_and_infinite_line():
    h = 0.7
    for L in (0.5, 3.0, 50.0, 5000.0):
        v = impl.segment_velocity(np.array([[0.0, 0.0, h]]), np.array([[0.0, -L, 0.0]]), np.array([[0.0, L, 0.0]]), 2.0)
        exact = 2.0 / (2 * np.pi * h) * L / np.sqrt(L ** 2 + h ** 2)
        assert np.linalg.norm(v) == pytest.approx(exact, rel=1e-12)
    assert np.linalg.norm(v) == pytest.approx(2.0 / (2 * np.pi * h), rel=1e-7)   # long segment -> infinite line


def test_right_hand_rule_linearity_and_reversal():
    P = np.array([[0.3, 0.2, 0.5]])
    A, B = np.array([[0.0, -1.0, 0.0]]), np.array([[0.0, 1.0, 0.0]])
    v = impl.segment_velocity(np.array([[0.0, 0.0, 1.0]]), A, B)[0, 0]
    assert v[0] > 0 and abs(v[1]) < 1e-15 and abs(v[2]) < 1e-15             # +y vortex, point above: flow +x
    v1 = impl.segment_velocity(P, A, B, 1.0)
    npt.assert_allclose(impl.segment_velocity(P, A, B, 3.0), 3 * v1)
    npt.assert_allclose(impl.segment_velocity(P, B, A, 1.0), -v1)
    on_line = impl.segment_velocity(np.array([[0.0, 0.5, 0.0]]), A, B)
    assert np.all(np.isfinite(on_line)) and np.allclose(on_line, 0.0)


def test_semi_infinite_line_is_the_limit_of_a_long_segment():
    rng = np.random.default_rng(2)
    P = rng.uniform(-2, 2, size=(6, 3))
    A = np.array([[0.1, 0.4, -0.2]])
    d = np.array([1.0, 0.0, 0.0])
    v_si = impl.semi_infinite_velocity(P, A, d, 1.5)
    v_long = impl.segment_velocity(P, A, A + 1e6 * d, 1.5)
    npt.assert_allclose(v_si, v_long, rtol=1e-6, atol=1e-12)
    # in the plane through its end, half the infinite line
    v = impl.semi_infinite_velocity(np.array([[0.0, 0.0, 0.8]]), np.array([[0.0, 0.0, 0.0]]), d, 1.0)
    assert np.linalg.norm(v) == pytest.approx(0.5 / (2 * np.pi * 0.8), rel=1e-12)


# ---------------------------------------------------------------- horseshoe: Stokes, Kelvin, Helmholtz

HS = lambda P: impl.horseshoe_velocity(P, np.array([-1.0]), np.array([1.0]), np.array([2.0]))


def test_circulation_round_each_leg_and_the_whole_wake():
    ex, ey = np.array([1.0, 0, 0]), np.array([0, 1.0, 0])
    assert impl.loop_circulation(HS, np.array([4.0, 1.0, 0.0]), ex, 0.3) == pytest.approx(2.0, rel=1e-9)
    assert impl.loop_circulation(HS, np.array([4.0, -1.0, 0.0]), ex, 0.3) == pytest.approx(-2.0, rel=1e-9)
    assert impl.loop_circulation(HS, np.array([0.0, 0.2, 0.0]), ey, 0.3) == pytest.approx(2.0, rel=1e-9)
    assert abs(impl.loop_circulation(HS, np.array([6.0, 0.0, 0.0]), ex, 3.0)) < 1e-9   # no net wake circulation
    assert abs(impl.loop_circulation(HS, np.array([0.0, 3.0, 2.0]), ex, 0.4)) < 1e-12  # empty loop


def test_horseshoe_field_is_divergence_free():
    rng = np.random.default_rng(4)
    P = rng.uniform([-1, -2, 0.2], [3, 2, 1.0], size=(8, 3))
    h = 1e-4
    div = np.zeros(len(P))
    for k in range(3):
        e = np.zeros(3); e[k] = h
        div += (HS(P + e)[:, k] - HS(P - e)[:, k]) / (2 * h)
    scale = np.linalg.norm(HS(P), axis=1) / 1.0
    assert np.all(np.abs(div) < 1e-6 * np.maximum(scale, 1e-3))


def test_single_horseshoe_downwash_doubles_downstream():
    b, g = 4.0, 3.0
    hs = lambda P: impl.horseshoe_velocity(P, np.array([-b / 2]), np.array([b / 2]), np.array([g]))
    w0 = -hs(np.array([[0.0, 0.0, 0.0]]))[0, 2]
    wT = -hs(np.array([[1e5 * b, 0.0, 0.0]]))[0, 2]
    assert w0 == pytest.approx(g / (np.pi * b), rel=1e-12)
    assert wT == pytest.approx(2 * w0, rel=1e-8)
    up = -hs(np.array([[-1e5 * b, 0.0, 0.0]]))[0, 2]
    assert abs(up) < 1e-9                                                     # nothing far upstream


# ---------------------------------------------------------------- elliptic wing from horseshoes

def test_elliptic_load_has_uniform_downwash_and_e_equal_one():
    b, g0 = 10.0, 13.9
    ed, yc = impl.cosine_layout(b, 40)
    r = impl.lift_and_drag(ed, impl.elliptic_gamma(yc, b, g0), yc, 55.0)
    npt.assert_allclose(r["w"], g0 / (2 * b), rtol=1e-3)
    npt.assert_allclose(r["wT"], 2 * r["w"], rtol=1e-12)                      # Trefftz = 2 x lifting line
    assert r["Di_near"] == pytest.approx(r["Di_trefftz"], rel=1e-12)
    assert impl.span_efficiency(r["L"], r["Di_trefftz"], b, 55.0) == pytest.approx(1.0, abs=1e-10)


def test_elliptic_drag_converges_at_second_order_to_the_closed_form():
    b, S, V, CL = 10.0, 12.0, 55.0, 0.33
    g0 = impl.elliptic_gamma0(CL, V, S, b)
    q = 0.5 * 1.225 * V ** 2
    errs = []
    for n in (10, 20, 40, 80):
        ed, yc = impl.cosine_layout(b, n)
        r = impl.lift_and_drag(ed, impl.elliptic_gamma(yc, b, g0), yc, V)
        errs.append(abs(r["Di_trefftz"] / (q * S) - CL ** 2 / (np.pi * b ** 2 / S)))
    orders = np.log(np.array(errs[:-1]) / np.array(errs[1:])) / np.log(2)
    assert np.all(orders > 1.9) and np.all(orders < 2.1)
    assert r["L"] == pytest.approx(q * S * CL, rel=1e-3)


def test_two_term_load_matches_glauert_span_efficiency():
    b, V = 8.0, 30.0
    for a3 in (0.1, 0.3, -0.1):
        ed, yc = impl.cosine_layout(b, 120)
        phi = np.arccos(-2 * yc / b)
        G = np.sin(phi) + a3 * np.sin(3 * phi)
        r = impl.lift_and_drag(ed, G, yc, V)
        e = impl.span_efficiency(r["L"], r["Di_trefftz"], b, V)
        assert e == pytest.approx(1 / (1 + 3 * a3 ** 2), rel=2e-4)            # e = 1/(1 + sum n (A_n/A_1)^2)
        npt.assert_allclose(r["w"], r["w"][::-1], rtol=1e-9)                 # symmetric load, symmetric downwash


def test_uniform_load_drag_grows_without_bound_as_the_tip_vortices_concentrate():
    # a constant Gamma (one big horseshoe) sheds point vortices at the tips: the Trefftz drag
    # diverges like ln(n) as the layout resolves them, so e falls towards zero
    b, V = 10.0, 50.0
    es = []
    for n in (10, 40, 160):
        ed, yc = impl.cosine_layout(b, n)
        r = impl.lift_and_drag(ed, np.ones(n), yc, V)
        es.append(impl.span_efficiency(r["L"], r["Di_trefftz"], b, V))
    assert es[0] < 1 and es[0] > es[1] > es[2]


# ---------------------------------------------------------------- drag polar

def test_best_glide_and_minimum_power_are_optima():
    CD0, e, AR = 0.03, 0.8, 8.333
    bg = impl.best_glide(CD0, e, AR)
    CL = np.linspace(0.05, 2.5, 200001)
    CD = impl.drag_polar(CL, CD0, e, AR)
    assert CL[np.argmax(CL / CD)] == pytest.approx(bg["CL_md"], abs=1e-4)
    assert np.max(CL / CD) == pytest.approx(bg["LD_max"], rel=1e-8)
    assert CL[np.argmax(CL ** 1.5 / CD)] == pytest.approx(bg["CL_mp"], abs=1e-4)
    assert bg["CL_mp"] == pytest.approx(np.sqrt(3) * bg["CL_md"])
    # induced drag equals parasite drag at best L/D, and three times it at minimum power
    assert impl.drag_polar(bg["CL_md"], 0.0, e, AR) == pytest.approx(CD0)
    assert impl.drag_polar(bg["CL_mp"], 0.0, e, AR) == pytest.approx(3 * CD0)


def test_speeds_and_power_curve():
    m, S, CD0, e, AR = 750.0, 12.0, 0.03, 0.8, 8.333
    bg = impl.best_glide(CD0, e, AR)
    Vmd, Vmp = impl.speed_for_cl(m, S, bg["CL_md"]), impl.speed_for_cl(m, S, bg["CL_mp"])
    assert Vmp == pytest.approx(Vmd / 3 ** 0.25)
    V = np.linspace(15, 80, 65001)
    d = impl.drag_and_power(V, m, S, CD0, e, AR)
    assert V[np.argmin(d["D"])] == pytest.approx(Vmd, abs=2e-3)
    assert V[np.argmin(d["P"])] == pytest.approx(Vmp, abs=2e-3)
    dm = impl.drag_and_power(Vmd, m, S, CD0, e, AR)
    assert float(dm["D0"]) == pytest.approx(float(dm["Di"]))
    assert float(dm["D"]) == pytest.approx(m * 9.81 / bg["LD_max"])
    assert impl.speed_for_cl(4 * m, S, 0.5) == pytest.approx(2 * impl.speed_for_cl(m, S, 0.5))


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "induced_drag.png").stat().st_size > 10_000
