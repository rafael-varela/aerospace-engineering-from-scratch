"""Public tests for P07.12 Heat Transfer: Conduction. They assert physics and numerics, not magic numbers:
the finite-volume solver reproduces the Fourier-series solution of a rectangle with one heated edge and
converges at second order; a layered 1D wall matches its series-resistance closed form to round-off;
a 2D sheet losing heat from its faces matches the 1D fin closed form at second order; uniform
generation between two clamped edges gives the exact parabola; the global energy balance closes to
round-off; the transient matches the semi-infinite erfc solutions (fixed surface temperature and
convection) with second-order error, conserves energy exactly, and the explicit scheme is stable only
below Fourier number 1/2; and the heat-sink wall behaves as the semi-infinite, lumped and finite-wall
limits say it should.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import pytest

from conftest import load_impl

impl = load_impl(__file__)

CU = impl.MATERIALS["copper"]
DIRICHLET_TOP = {"left": ("T", 0.0), "right": ("T", 0.0), "bottom": ("T", 0.0), "top": ("T", 1.0)}


# ---------------------------------------------------------------- closed forms

def test_wall_series_resistances():
    layers = [(0.002, 10.0), (0.005, 400.0)]
    q, Ts = impl.wall_heat_flux(layers, 1000.0, 300.0)                      # surfaces at the fluid T
    assert q == pytest.approx(700.0 / (0.002 / 10 + 0.005 / 400), rel=1e-14)
    assert Ts[0] == pytest.approx(1000.0) and Ts[-1] == pytest.approx(300.0)
    q2, Ts2 = impl.wall_heat_flux(layers, 1000.0, 300.0, 2000.0, 5000.0)
    assert q2 < q and 300.0 < Ts2[-1] < Ts2[1] < Ts2[0] < 1000.0
    assert (1000.0 - Ts2[0]) * 2000.0 == pytest.approx(q2, rel=1e-12)       # same flux through every resistance


def test_fourier_series_satisfies_its_boundary_conditions():
    x = np.linspace(0.05, 0.95, 7)
    assert np.allclose(impl.rectangle_fourier(x, 0.0, 1.0, 1.0), 0.0, atol=1e-15)
    assert np.allclose(impl.rectangle_fourier(0.0, x, 1.0, 1.0), 0.0, atol=1e-12)
    # by symmetry (superposing the four rotations gives T = 1) the centre of a square is exactly 1/4
    assert float(impl.rectangle_fourier(0.5, 0.5, 1.0, 1.0)) == pytest.approx(0.25, abs=1e-12)


def test_semi_infinite_solutions_limits():
    a, k = 1e-4, 400.0
    x = np.array([0.0, 0.01, 0.1])
    T = impl.semi_infinite_step(x, 1.0, a, 300.0, 1000.0)
    assert T[0] == 1000.0 and T[-1] == pytest.approx(300.0, abs=1e-6)
    # convection with h -> infinity recovers the fixed-surface solution; h -> 0 leaves the solid at T_i
    Th = impl.semi_infinite_convection(x, 1.0, a, k, 1e12, 300.0, 1000.0)   # surface lag ~ 700/(beta sqrt(pi)) K
    np.testing.assert_allclose(Th, T, atol=1e-4)
    assert np.allclose(impl.semi_infinite_convection(x, 1.0, a, k, 1e-12, 300.0, 1000.0), 300.0, atol=1e-9)


# ---------------------------------------------------------------- steady finite volumes

def test_rectangle_with_one_heated_edge_matches_fourier_series():
    errs = []
    for n in (21, 41):
        s = impl.conduction_system_2d(1.0, 1.0, n, n, 1.0, 0.0, DIRICHLET_TOP)
        T = impl.solve_steady(s)
        X, Y = np.meshgrid(s["x"], s["y"])
        away = Y <= 0.75                                                        # the corners at y = 1 are singular
        errs.append(np.max(np.abs(T - impl.rectangle_fourier(X, Y, 1.0, 1.0))[away]))
    assert errs[-1] < 1e-3
    assert errs[0] / errs[1] == pytest.approx(4.0, rel=0.1)                    # second order


def test_second_order_convergence_on_a_smooth_problem():
    errs = []
    for n in (11, 21, 41, 81):
        s = impl.conduction_system_2d(2.0, 1.0, 2 * n - 1, n, 1.0, 0.0, {"left": ("T", 0.0), "right": ("T", 0.0),
                                      "bottom": ("T", 0.0), "top": ("T", lambda x: np.sin(np.pi * x / 2))})
        T = impl.solve_steady(s)
        X, Y = np.meshgrid(s["x"], s["y"])
        errs.append(np.max(np.abs(T - np.sin(np.pi * X / 2) * np.sinh(np.pi * Y / 2) / np.sinh(np.pi / 2))))
    orders = impl.observed_order(errs)
    assert np.all(np.abs(orders - 2.0) < 0.1)


def test_uniform_generation_between_clamped_edges_is_the_exact_parabola():
    kt, q, H, Tw = 0.05, 900.0, 0.1, 343.15
    s = impl.conduction_system_2d(0.16, H, 17, 11, kt, q, {"bottom": ("T", Tw), "top": ("T", Tw)})
    T = impl.solve_steady(s)
    exact = Tw + q * s["y"] * (H - s["y"]) / (2 * kt)
    np.testing.assert_allclose(T, np.tile(exact[:, None], (1, 17)), rtol=1e-12)   # quadratic: no truncation error


def test_global_energy_balance_closes_to_round_off():
    kt = impl.board_sheet_conductance()
    T, s, hb = impl.solve_card(kt, 41)
    assert hb["generated"] == pytest.approx(15.0, rel=1e-12)                  # chip overlap areas exact on any grid
    assert abs(hb["imbalance"]) < 1e-9 * hb["generated"]
    assert hb["out_fixed"] == pytest.approx(15.0, rel=1e-10)
    # flux and convection edges too
    s = impl.conduction_system_2d(0.1, 0.05, 21, 11, 0.2, 300.0, {"left": ("flux", 40.0), "right": ("h", 15.0, 290.0),
                                                                   "bottom": ("T", 300.0)}, h_face=8.0, T_face=280.0)
    T = impl.solve_steady(s)
    hb = impl.heat_balance(s, T)
    assert hb["residual_free"] < 1e-9 and abs(hb["imbalance"]) < 1e-10 * (abs(hb["generated"]) + 1.0)


def test_layered_wall_matches_closed_form_to_round_off():
    layers = [(0.001, 11.6), (0.004, 398.0)]
    q, Ts = impl.wall_heat_flux(layers, 3000.0, 300.0, 3000.0, 2.0e4)
    n = 51
    k = np.where(np.arange(n - 1) < 10, 11.6, 398.0)                           # interface on node 10
    s = impl.conduction_system_1d(0.005, n, k, 1.0, ("h", 3000.0, 3000.0), ("h", 2.0e4, 300.0))
    T = impl.solve_steady(s)
    assert T[0] == pytest.approx(Ts[0], rel=1e-12) and T[10] == pytest.approx(Ts[1], rel=1e-12)
    assert T[-1] == pytest.approx(Ts[2], rel=1e-12)
    hb = impl.heat_balance(s, T)
    assert abs(hb["imbalance"]) < 1e-6 and hb["in_edges"] == pytest.approx(0.0, abs=1e-6)   # in = out


def test_sheet_with_face_losses_matches_fin_closed_form():
    kt, L, h, W = 167.0e-3, 0.04, 50.0, 0.02
    m = np.sqrt(h / kt)
    errs = []
    for n in (11, 21, 41):
        s = impl.conduction_system_2d(L, W, n, 5, kt, 0.0, {"left": ("T", 40.0)}, h_face=h, T_face=0.0)
        T = impl.solve_steady(s)
        assert np.allclose(T, T[0], rtol=1e-12)                               # uniform across the width
        errs.append(np.max(np.abs(T[0] - impl.fin_closed_form(s["x"], L, m, 40.0))))
    assert np.all(impl.observed_order(errs) > 1.9)
    hb = impl.heat_balance(s, T)
    assert hb["out_faces"] == pytest.approx(impl.fin_heat_rate(kt * W, L, m, 40.0), rel=1e-3)
    assert hb["out_fixed"] == pytest.approx(-hb["out_faces"], rel=1e-10)      # the base supplies what the faces lose


def test_card_heat_frame_lowers_the_hot_spot():
    T0 = impl.solve_card(impl.board_sheet_conductance(), 41)[0].max()
    T1 = impl.solve_card(impl.board_sheet_conductance(t_frame=1.5e-3), 41)[0].max()
    assert 343.15 < T1 < T0
    # temperature rise above the wall scales as 1 / (k t): linear problem
    kt = impl.board_sheet_conductance()
    T2 = impl.solve_card(2 * kt, 41)[0].max()
    assert (T2 - 343.15) == pytest.approx((T0 - 343.15) / 2, rel=1e-10)


# ---------------------------------------------------------------- transients

@pytest.mark.parametrize("theta", [1.0, 0.5])
def test_transient_matches_erfc_at_second_order(theta):
    a = impl.diffusivity("copper")
    errs = []
    for n in (26, 51, 101):
        s = impl.conduction_system_1d(0.05, n, CU["k"], CU["rho"] * CU["cp"], ("T", 1000.0), ("adiabatic",))
        dx = 0.05 / (n - 1)
        steps = int(round(0.5 / (0.5 * dx ** 2 / a)))                         # dt proportional to dx^2
        out = impl.march(s, 300.0, 0.5 / steps, steps, theta=theta)
        errs.append(np.max(np.abs(out["T"] - impl.semi_infinite_step(s["x"], 0.5, a, 300.0, 1000.0))))
    assert errs[-1] < 0.25                                                    # K, out of a 700 K step
    assert np.all(impl.observed_order(errs) > 1.8)


def test_transient_energy_conservation_and_absorbed_heat():
    a = impl.diffusivity("copper")
    s = impl.conduction_system_1d(0.05, 101, CU["k"], CU["rho"] * CU["cp"], ("T", 1000.0), ("adiabatic",))
    out = impl.march(s, 300.0, 1e-3, 500, theta=1.0)
    assert out["E_stored"] == pytest.approx(out["Q_in"], rel=1e-10)          # exact for the scheme
    E = out["E_stored"] + s["C"][0] * 700.0                                   # plus the surface half-cell's jump
    assert E == pytest.approx(impl.semi_infinite_step_energy(0.5, CU["k"], a, 300.0, 1000.0), rel=2e-3)
    s = impl.conduction_system_1d(0.05, 101, CU["k"], CU["rho"] * CU["cp"], ("h", 3000.0, 3000.0), ("adiabatic",))
    out = impl.march(s, 300.0, 1e-3, 1000, theta=0.5)
    assert out["E_stored"] == pytest.approx(out["Q_in"], rel=1e-10)
    np.testing.assert_allclose(out["T"], impl.semi_infinite_convection(s["x"], 1.0, a, CU["k"], 3000.0, 300.0, 3000.0), atol=0.3)


def test_explicit_scheme_stability_limit():
    a = impl.diffusivity("copper")
    s = impl.conduction_system_1d(0.05, 51, CU["k"], CU["rho"] * CU["cp"], ("T", 1000.0), ("adiabatic",))
    dx = 0.05 / 50
    stable = impl.march(s, 300.0, 0.45 * dx ** 2 / a, 400, theta=0.0)
    unstable = impl.march(s, 300.0, 0.55 * dx ** 2 / a, 400, theta=0.0)
    assert 300.0 - 1e-9 <= stable["T"].min() and stable["T"].max() <= 1000.0 + 1e-9   # maximum principle
    assert np.max(np.abs(unstable["T"])) > 1e6
    implicit = impl.march(s, 300.0, 5.0 * dx ** 2 / a, 100, theta=1.0)        # Fo = 5: still bounded
    assert 300.0 - 1e-9 <= implicit["T"].min() and implicit["T"].max() <= 1000.0 + 1e-9


def test_steady_limit_of_the_transient():
    s = impl.conduction_system_2d(0.16, 0.1, 17, 11, 0.06, 900.0, {"bottom": ("T", 343.15), "top": ("T", 343.15)},
                                  rho_c_t=1800 * 1600 * 1.6e-3)
    out = impl.march(s, 343.15, 50.0, 400, theta=1.0)
    np.testing.assert_allclose(out["T"].reshape(11, 17), impl.solve_steady(s), rtol=1e-8)


def test_heat_sink_wall_limits():
    # a thick copper wall at short times is semi-infinite; a thin one heats faster
    t_si = impl.semi_infinite_time_to_limit("copper", 3000.0, 3000.0, 300.0, 800.0)
    t_thick, res = impl.time_to_surface_limit("copper", 0.10, 3000.0, 3000.0, 300.0, 800.0)
    t_thin, _ = impl.time_to_surface_limit("copper", 0.010, 3000.0, 3000.0, 300.0, 800.0)
    assert t_thick == pytest.approx(t_si, rel=0.01)
    assert t_thin < t_thick
    assert res["E_stored"] == pytest.approx(res["Q_in"], rel=1e-10)
    # a thin wall with Bi = hL/k = 0.015 behaves almost as a lumped mass; its hot face runs slightly ahead
    # of the mean temperature (by a fraction of order Bi of the gas-to-wall difference), so the
    # lumped time is a little long: within 3 %
    L, h = 0.002, 3000.0
    t_l = CU["rho"] * CU["cp"] * L / h * np.log(2700.0 / 2200.0)
    t_fd = impl.time_to_surface_limit("copper", L, h, 3000.0, 300.0, 800.0, dt=1e-3)[0]
    assert t_fd == pytest.approx(t_l, rel=0.03) and t_fd < t_l
    # a low-conductivity wall reaches the limit far sooner than copper
    assert impl.time_to_surface_limit("Inconel X-750", 0.025, h, 3000.0, 300.0, 800.0, dt=1e-3)[0] < 0.1 * t_si


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "conduction.png").stat().st_size > 10_000
