"""Public tests for P03.12 Finite Volume Methods. They assert properties, with SciPy as the oracle.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-03/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy.integrate import quad
from scipy.linalg import circulant

from conftest import load_impl

impl = load_impl(__file__)

RNG = np.random.default_rng(312)
A, L = 12.0, 2.4
GAUSS = lambda x: impl.gaussian_pulse(x, 0.6, 0.1)                           # noqa: E731
SQUARE = lambda x: impl.square_pulse(x, 0.6, 0.3)                            # noqa: E731
LINEAR = ("upwind", "lax_friedrichs", "lax_wendroff")
ALL = LINEAR + ("minmod",)


# ---------------------------------------------------------------- grid, averages, exact solution

def test_cell_centres_and_cell_averages_match_scipy_quad():
    x, dx = impl.cell_centres(L, 8)
    assert dx == pytest.approx(0.3) and x[0] == pytest.approx(0.15) and x[-1] == pytest.approx(2.25)
    with pytest.raises(ValueError):
        impl.cell_centres(L, 0)
    n = 24
    avg = impl.cell_averages(GAUSS, L, n)
    x, dx = impl.cell_centres(L, n)
    for i in (3, 6, 7, 20):
        exact, _ = quad(GAUSS, x[i] - dx / 2, x[i] + dx / 2)
        assert avg[i] == pytest.approx(exact / dx, rel=1e-6)                 # 4 Gauss points on a cell of width sigma
    # Averaging a cubic is exact with 2 Gauss points or more; averaging differs from sampling at O(dx^2).
    cubic = lambda xx: xx ** 3 - xx                                           # noqa: E731
    npt.assert_allclose(impl.cell_averages(cubic, L, n, n_gauss=2), ((x + dx / 2) ** 4 - (x - dx / 2) ** 4) / (4 * dx) - x, rtol=1e-12)
    assert np.max(np.abs(avg - GAUSS(x))) > 1e-4
    # Total mass of the averaged Gaussian equals its integral.
    assert np.sum(avg) * dx == pytest.approx(0.1 * np.sqrt(2 * np.pi), rel=1e-8)


def test_exact_advection_shifts_periodically_and_pulses_have_the_right_support():
    x, dx = impl.cell_centres(L, 240)
    u = impl.exact_advection(GAUSS, x, A, 0.1, L)
    assert x[np.argmax(u)] == pytest.approx(0.6 + 1.2, abs=dx)
    npt.assert_allclose(impl.exact_advection(GAUSS, x, A, L / A, L), GAUSS(x), atol=1e-12)   # one lap
    npt.assert_allclose(impl.exact_advection(GAUSS, x, -A, 0.1, L), impl.exact_advection(GAUSS, x, A, L / A - 0.1, L), atol=1e-12)
    sq = SQUARE(x)
    assert np.sum(sq) * dx == pytest.approx(0.3, abs=dx)
    assert set(np.unique(sq)) == {0.0, 1.0}
    assert impl.courant_number(A, 0.001, 0.02) == pytest.approx(0.6)
    with pytest.raises(ValueError):
        impl.courant_number(A, -0.1, 0.02)


# ---------------------------------------------------------------- fluxes and the update

def test_fluxes_reduce_to_the_textbook_forms_and_are_consistent():
    u = RNG.uniform(0, 1, 10)
    nu = 0.4
    npt.assert_allclose(impl.numerical_flux(u, A, nu, "upwind"), A * u)
    npt.assert_allclose(impl.numerical_flux(u, A, nu, "lax_wendroff"), 0.5 * A * ((1 + nu) * u + (1 - nu) * np.roll(u, -1)))
    F_lf = impl.numerical_flux(u, A, nu, "lax_friedrichs")
    npt.assert_allclose(F_lf, 0.5 * A * (u + np.roll(u, -1)) - 0.5 * (A / nu) * (np.roll(u, -1) - u))
    # Consistency: a uniform state has flux a u through every interface, for every scheme.
    for scheme in ALL:
        npt.assert_allclose(impl.numerical_flux(np.full(10, 0.7), A, nu, scheme), A * 0.7)
    # At nu = 1 Lax-Wendroff and Lax-Friedrichs reduce to upwind.
    npt.assert_allclose(impl.numerical_flux(u, A, 1.0, "lax_wendroff"), A * u)
    npt.assert_allclose(impl.numerical_flux(u, A, 1.0, "lax_friedrichs"), A * u)
    # minmod on a linear profile has the full (unlimited) slope away from the periodic wrap.
    lin = np.arange(10.0)
    F = impl.numerical_flux(lin, A, nu, "minmod")
    npt.assert_allclose(F[1:-1], A * (lin[1:-1] + 0.5 * (1 - nu)))
    npt.assert_allclose(impl.minmod(np.array([1.0, -2.0, 3.0, 0.0]), np.array([2.0, -1.0, -3.0, 5.0])), [1.0, -1.0, 0.0, 0.0])
    with pytest.raises(ValueError):
        impl.numerical_flux(u, A, nu, "godunov-roe")
    with pytest.raises(ValueError):
        impl.numerical_flux(u, A, 0.0, "lax_friedrichs")


def test_update_conserves_mass_refuses_cfl_violations_and_handles_negative_speed():
    n = 60
    x, dx = impl.cell_centres(L, n)
    u0 = impl.cell_averages(SQUARE, L, n)
    dt = 0.5 * dx / A
    for scheme in ALL:
        u = impl.advect(u0, A, dx, dt, 37, scheme)
        assert np.sum(u) * dx == pytest.approx(np.sum(u0) * dx, abs=1e-13)      # conservation to rounding
        # Reversing the flow direction mirrors the solution.
        u_neg = impl.advect(u0[::-1], -A, dx, dt, 37, scheme)
        npt.assert_allclose(u_neg[::-1], u, atol=1e-12)
    with pytest.raises(ValueError):
        impl.fv_step(u0, A, 1.1 * dx / A, dx, "upwind")
    # Just past the CFL limit the solution grows; within it, it does not.
    for scheme in LINEAR:
        with pytest.raises(ValueError):
            impl.advect(u0, A, dx, 1.2 * dx / A, 10, scheme)
    grown = u0.copy()
    for _ in range(60):
        grown = grown - (1.2 * dx / A / dx) * (impl.numerical_flux(grown, A, 1.2, "upwind") - np.roll(impl.numerical_flux(grown, A, 1.2, "upwind"), 1))
    assert np.max(np.abs(grown)) > 100.0


def test_nu_equal_one_shifts_the_data_exactly_by_one_cell():
    n = 48
    x, dx = impl.cell_centres(L, n)
    u0 = impl.cell_averages(GAUSS, L, n)
    for scheme in ALL:
        u = impl.advect(u0, A, dx, dx / A, 7, scheme)
        npt.assert_allclose(u, np.roll(u0, 7), atol=1e-12)


def test_upwind_is_monotone_and_lax_wendroff_is_not():
    n = 120
    x, dx = impl.cell_centres(L, n)
    u0 = impl.cell_averages(SQUARE, L, n)
    dt = 0.5 * dx / A
    tv0 = impl.total_variation(u0)
    assert tv0 == pytest.approx(2.0, abs=1e-12)
    u_up = impl.advect(u0, A, dx, dt, 120, "upwind")
    u_mm = impl.advect(u0, A, dx, dt, 120, "minmod")
    u_lw = impl.advect(u0, A, dx, dt, 120, "lax_wendroff")
    for u in (u_up, u_mm):
        assert np.min(u) >= -1e-12 and np.max(u) <= 1.0 + 1e-12            # no new extrema
        assert impl.total_variation(u) <= tv0 + 1e-12                         # TVD
    assert np.min(u_lw) < -0.1 and np.max(u_lw) > 1.1                          # overshoot and undershoot
    assert impl.total_variation(u_lw) > tv0 + 0.5
    # The limited scheme is sharper than upwind: smaller L1 error on the jump.
    ref = impl.cell_averages(lambda xx: impl.exact_advection(SQUARE, xx, A, 120 * dt, L), L, n)
    assert impl.l1_error(u_mm, ref, dx) < 0.5 * impl.l1_error(u_up, ref, dx)


# ---------------------------------------------------------------- von Neumann analysis

def test_amplification_factors_match_circulant_eigenvalues_and_the_cfl_limit():
    n = 40
    th = 2 * np.pi * np.arange(n) / n
    e = np.zeros(n)
    e[0] = 1.0
    x, dx = impl.cell_centres(L, n)
    for scheme in LINEAR:
        for nu in (0.3, 0.8, 1.0):
            dt = nu * dx / A
            # One step applied to the unit vectors gives the update matrix, which is circulant.
            G = np.column_stack([impl.fv_step(np.roll(e, j), A, dt, dx, scheme) for j in range(n)])
            npt.assert_allclose(G, circulant(G[:, 0]), atol=1e-12)
            eig = np.linalg.eigvals(G)
            g = impl.amplification_factor(scheme, nu, th)
            npt.assert_allclose(np.sort_complex(np.round(eig, 10)), np.sort_complex(np.round(g, 10)), atol=1e-8)
            assert np.max(np.abs(g)) <= 1.0 + 1e-12
        assert np.max(np.abs(impl.amplification_factor(scheme, 1.2, th))) > 1.0 + 1e-6   # unstable past CFL
        assert impl.amplification_factor(scheme, 0.5, 0.0) == pytest.approx(1.0)          # constant mode untouched
    # Closed forms for |g|^2.
    nu, t = 0.4, np.linspace(0, np.pi, 50)
    npt.assert_allclose(np.abs(impl.amplification_factor("upwind", nu, t)) ** 2, 1 - 2 * nu * (1 - nu) * (1 - np.cos(t)), atol=1e-12)
    npt.assert_allclose(np.abs(impl.amplification_factor("lax_wendroff", nu, t)) ** 2, 1 - nu ** 2 * (1 - nu ** 2) * (1 - np.cos(t)) ** 2, atol=1e-12)
    npt.assert_allclose(np.abs(impl.amplification_factor("lax_friedrichs", nu, t)) ** 2, np.cos(t) ** 2 + nu ** 2 * np.sin(t) ** 2, atol=1e-12)
    # Lax-Wendroff lags the exact phase for resolved modes; at nu = 1 every scheme is exact.
    assert np.all(impl.phase_speed_ratio("lax_wendroff", 0.5, t[1:]) < 1.0 + 1e-12)
    npt.assert_allclose(impl.amplification_factor("upwind", 1.0, t), np.exp(-1j * t), atol=1e-12)
    with pytest.raises(ValueError):
        impl.amplification_factor("minmod", 0.5, t)


# ---------------------------------------------------------------- the modified equation

def test_upwind_diffuses_a_gaussian_exactly_as_its_modified_equation_predicts():
    n = 240
    x, dx = impl.cell_centres(L, n)
    u0 = impl.cell_averages(GAUSS, L, n)
    var0 = impl.pulse_variance(x, u0, dx, L)
    assert var0 == pytest.approx(0.1 ** 2 + dx ** 2 / 12, rel=1e-6)          # averaging adds dx^2/12
    for nu in (0.25, 0.5, 0.8):
        dt = nu * dx / A
        steps = 120
        u = impl.advect(u0, A, dx, dt, steps, "upwind")
        D = impl.numerical_diffusivity(A, dx, nu)
        assert D == pytest.approx(0.5 * A * dx * (1 - nu))
        growth = impl.pulse_variance(x, u, dx, L) - var0
        assert growth == pytest.approx(2 * D * steps * dt, rel=2e-3)        # the modified equation holds to 0.2 %
    assert impl.numerical_diffusivity(A, dx, 1.0) == 0.0
    with pytest.raises(ValueError):
        impl.pulse_variance(x, np.zeros(n), dx, L)


def test_small_courant_upwind_tends_to_the_semi_discrete_scipy_solution():
    n = 60
    x, dx = impl.cell_centres(L, n)
    u0 = impl.cell_averages(GAUSS, L, n)
    t_final = 0.02
    ref = impl.semi_discrete_upwind_oracle(u0, A, dx, t_final)
    errs = []
    for nu in (0.2, 0.1, 0.05):
        dt = nu * dx / A
        u = impl.advect(u0, A, dx, dt, int(round(t_final / dt)), "upwind")
        errs.append(np.max(np.abs(u - ref)))
    npt.assert_allclose(impl.observed_order(np.array(errs)), 1.0, atol=0.1)   # O(nu) gap to the method of lines
    assert errs[-1] < 2e-2


# ---------------------------------------------------------------- convergence (observed order)

def test_observed_orders_on_a_smooth_slug_are_one_for_upwind_and_two_for_lax_wendroff():
    ns = [240, 480, 960, 1920]
    e_up = impl.convergence_study("upwind", ns, GAUSS, A, L, 0.1, 0.5)
    e_lw = impl.convergence_study("lax_wendroff", ns, GAUSS, A, L, 0.1, 0.5)
    e_mm = impl.convergence_study("minmod", ns, GAUSS, A, L, 0.1, 0.5)
    e_lf = impl.convergence_study("lax_friedrichs", ns, GAUSS, A, L, 0.1, 0.5)
    for e in (e_up, e_lw, e_mm, e_lf):
        assert np.all(np.diff(e) < 0)
    o_up, o_lw, o_mm = impl.observed_order(e_up), impl.observed_order(e_lw), impl.observed_order(e_mm)
    assert 0.8 < o_up[-1] < 1.05 and np.all(np.diff(o_up) > 0)              # first order, approached from below
    npt.assert_allclose(o_lw, 2.0, atol=0.05)                                 # second order
    assert 1.5 < o_mm[-1] < 2.0                                               # limiter clips the extremum: between 1 and 2
    assert np.all(e_lf > e_up) and np.all(e_lw < e_up)
    npt.assert_allclose(impl.observed_order(np.array([8.0, 4.0, 2.0])), [1.0, 1.0])


def test_observed_order_on_a_discontinuity_drops_to_one_half_for_upwind():
    ns = [240, 480, 960, 1920]
    e_up = impl.convergence_study("upwind", ns, SQUARE, A, L, 0.1, 0.5)
    e_mm = impl.convergence_study("minmod", ns, SQUARE, A, L, 0.1, 0.5)
    npt.assert_allclose(impl.observed_order(e_up), 0.5, atol=0.05)           # L1 order 1/2 at a jump
    assert np.all(impl.observed_order(e_mm) > 0.6)                            # the TVD scheme does better (about 2/3)
    assert np.all(e_mm < e_up)


# ---------------------------------------------------------------- the sensor

def test_sensor_history_reports_arrival_and_the_schemes_distort_the_reading_as_predicted():
    c = impl.feed_line()
    n = 120
    x, dx = impl.cell_centres(c["L"], n)
    dt = 0.5 * dx / c["a"]
    u0 = impl.cell_averages(lambda xx: impl.square_pulse(xx, c["x_pulse"], c["width"]), c["L"], n)
    i_s = int(c["x_sensor"] / dx)
    t_front = (c["x_sensor"] - c["x_pulse"] - 0.5 * c["width"]) / c["a"]
    reads = {}
    for scheme in ALL:
        h = impl.sensor_history(u0, c["a"], dx, dt, 240, scheme, i_s)
        assert h.shape == (240,)
        t50 = impl.arrival_time(h, dt, 0.5)
        assert t50 == pytest.approx(t_front, abs=0.004)                      # 50 % crossing within a few ms
        reads[scheme] = (impl.arrival_time(h, dt, 0.1), np.max(h), np.min(h))
    assert reads["upwind"][1] < 0.9 and reads["lax_friedrichs"][1] < reads["upwind"][1]   # smeared peaks
    assert reads["lax_wendroff"][1] > 1.1 and reads["lax_wendroff"][2] < -0.1              # overshoot
    assert reads["minmod"][1] > 0.95 and reads["minmod"][2] > -1e-12                       # sharp and bounded
    assert reads["upwind"][0] < reads["minmod"][0] < t_front                                # diffusion arrives early
    assert np.isnan(impl.arrival_time(np.zeros(5), dt, 0.5))


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "fv_advection.png").stat().st_size > 1000
