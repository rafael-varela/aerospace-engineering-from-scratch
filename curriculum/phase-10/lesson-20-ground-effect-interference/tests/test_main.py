"""Public tests for P10.20 Ground Effect & Interference. They assert properties, not magic numbers:
the image system leaves no normal velocity on the ground plane; far from the ground the free-air
results return; the span load stays symmetric; lift is linear in angle of attack; lift rises and
induced drag falls monotonically as the wing approaches the ground; the lattice reproduces the
closed-form downwash reduction of a single horseshoe and agrees with the classical empirical
induced-drag ratio within its scatter; the Trefftz-plane downwash equals the 3-D Biot-Savart sum
far downstream; the downwash at a tail is consistent with the wake; main() writes its figure.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

B, C = 10.0, 1.2
ALPHA = np.radians(5.0)


@pytest.fixture(scope="module")
def free():
    return impl.wing_in_ground_effect(B, C, None, ALPHA, ns=24, nc=3)


def test_mirror_reflects_in_the_ground_plane():
    p = np.array([[1.0, 2.0, 3.0], [0.0, -1.0, 0.5]])
    q = impl.mirror(p)
    npt.assert_allclose(q[:, :2], p[:, :2])
    npt.assert_allclose(q[:, 2], -p[:, 2])
    npt.assert_allclose(impl.mirror(q), p)
    assert impl.mirror(p, 1.0)[0, 2] == pytest.approx(-1.0)


def test_images_cancel_the_normal_velocity_on_the_ground():
    g = impl.wing_in_ground_effect(B, C, 1.0, ALPHA, ns=24, nc=3)
    rng = np.random.default_rng(7)
    pts = np.column_stack([rng.uniform(-5, 20, 100), rng.uniform(-8, 8, 100)])
    assert impl.ground_plane_residual(g["lat"], g["gamma"], pts) < 1e-12


def test_far_from_the_ground_the_free_air_result_returns(free):
    g = impl.wing_in_ground_effect(B, C, 1000.0, ALPHA, ns=24, nc=3)
    assert g["CL"] == pytest.approx(free["CL"], rel=1e-5)
    assert g["CDi"] == pytest.approx(free["CDi"], rel=1e-4)
    npt.assert_allclose(g["G"], free["G"], rtol=1e-5)


def test_span_load_is_symmetric_in_and_out_of_ground_effect(free):
    npt.assert_allclose(free["G"], free["G"][::-1], rtol=1e-10)
    g = impl.wing_in_ground_effect(B, C, 0.5, ALPHA, ns=24, nc=3)
    npt.assert_allclose(g["G"], g["G"][::-1], rtol=1e-10)
    assert np.all(g["G"] > 0)


def test_lift_is_linear_in_sin_alpha():
    lat = impl.rectangular_lattice(B, C, 16, 2, z=1.0)
    g1 = impl.solve_lattice(lat, np.radians(2.0), ground=True)
    g2 = impl.solve_lattice(lat, np.radians(6.0), ground=True)
    ratio = impl.lift_coefficient(lat, g2) / impl.lift_coefficient(lat, g1)
    assert ratio == pytest.approx(np.sin(np.radians(6.0)) / np.sin(np.radians(2.0)), rel=1e-12)


def test_large_aspect_ratio_tends_to_two_pi():
    lat = impl.rectangular_lattice(2000.0, 1.0, 30, 2, spacing="uniform")
    g = impl.solve_lattice(lat, 0.05)
    assert impl.lift_coefficient(lat, g) / 0.05 == pytest.approx(2 * np.pi, rel=0.01)


def test_lift_rises_and_induced_drag_factor_falls_near_the_ground(free):
    K0 = free["CDi"] / free["CL"] ** 2
    CL, K = [], []
    for h in (5.0, 2.0, 1.0, 0.5):
        g = impl.wing_in_ground_effect(B, C, h, ALPHA, ns=24, nc=3)
        CL.append(g["CL"]); K.append(g["CDi"] / g["CL"] ** 2 / K0)
    assert np.all(np.diff(CL) > 0) and CL[0] > free["CL"]
    assert np.all(np.diff(K) < 0) and K[0] < 1


def test_single_horseshoe_matches_its_closed_form():
    # one strip, one chordwise panel: the Trefftz downwash at the centre falls by (4h/b)^2/(1+(4h/b)^2)
    y_edges = np.array([-0.5, 0.5])
    G = np.array([1.0])
    w0 = impl.trefftz_downwash(y_edges, G, np.array([0.0]))[0]
    assert w0 == pytest.approx(2 / np.pi)                      # two tip vortices at distance 1/2
    for h in (0.05, 0.2, 1.0):
        wh = impl.trefftz_downwash(y_edges, G, np.array([0.0]), h)[0]
        assert wh / w0 == pytest.approx(float(impl.single_horseshoe_trefftz_ratio(h / 1.0)), rel=1e-12)


def test_trefftz_downwash_equals_biot_savart_far_downstream():
    lat = impl.rectangular_lattice(B, C, 16, 2, z=0.8)
    gam = impl.solve_lattice(lat, ALPHA, ground=True)
    G = impl.strip_circulation(lat, gam)
    yc = 0.5 * (lat.y_edges[1:] + lat.y_edges[:-1])
    wT = impl.trefftz_downwash(lat.y_edges, G, yc, h=0.8)
    P = np.column_stack([np.full(len(yc), 1e6), yc, np.full(len(yc), 0.8)])
    V = impl.horseshoe_velocity(P, lat.A, lat.B) - impl.horseshoe_velocity(P, impl.mirror(lat.A), impl.mirror(lat.B))
    npt.assert_allclose(wT, -(V[..., 2] @ gam), rtol=1e-6, atol=1e-12)


def test_vlm_agrees_with_the_classical_empirical_ratio(free):
    # Wieselsberger's fit is an empirical correlation (stated range 0.033 < h/b < 0.25);
    # a 0.05 tolerance is about its typical scatter against other estimates
    K0 = free["CDi"] / free["CL"] ** 2
    for hb in (0.05, 0.1, 0.2):
        g = impl.wing_in_ground_effect(B, C, hb * B, ALPHA, ns=24, nc=3)
        assert g["CDi"] / g["CL"] ** 2 / K0 == pytest.approx(float(impl.wieselsberger_ratio(hb)), abs=0.05)


def test_far_field_downwash_gradient_formula():
    assert impl.downwash_gradient_far_field(2 * np.pi, 1e9) == pytest.approx(0.0, abs=1e-8)
    assert impl.downwash_gradient_far_field(5.0, 8.0) == pytest.approx(2 * 5.0 / (np.pi * 8.0))


def test_tail_downwash_is_consistent_with_the_wake(free):
    lat, gam = free["lat"], free["gamma"]
    G = impl.strip_circulation(lat, gam)
    bt = 2.83
    far = impl.tail_downwash(lat, gam, 1e5, 0.0, bt, n=20)
    yt = np.linspace(-bt / 2, bt / 2, 20)
    assert far == pytest.approx(np.mean(impl.trefftz_downwash(lat.y_edges, G, yt)), rel=1e-6)
    e1 = impl.tail_downwash(lat, gam, 4.8, 0.0, bt) / ALPHA
    e2 = impl.tail_downwash(lat, gam, 4.8, 1.0, bt) / ALPHA
    CLa = free["CL"] / ALPHA
    assert e2 < e1 < 1.0
    assert e1 == pytest.approx(impl.downwash_gradient_far_field(CLa, free["AR"]), rel=0.25)
    g = impl.wing_in_ground_effect(B, C, 1.0, ALPHA, ns=24, nc=3)
    eg = impl.tail_downwash(g["lat"], g["gamma"], 4.8, 1.0, bt, ground=True) / ALPHA
    assert eg < e1                                             # the ground straightens the flow at the tail


def test_lattice_convergence_is_monotone_with_a_positive_order():
    CL = [impl.wing_in_ground_effect(B, C, 1.0, ALPHA, ns=ns, nc=2)["CL"] for ns in (10, 20, 40, 80)]
    p = impl.observed_order(np.abs(np.diff(CL)))
    assert np.all(np.diff(CL) < 0) and np.all((p > 0.5) & (p < 2.5))


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path, quick=True)
    assert (tmp_path / "ground_effect.png").stat().st_size > 10_000
