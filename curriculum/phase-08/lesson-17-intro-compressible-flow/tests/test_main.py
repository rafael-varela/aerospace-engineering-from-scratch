"""Public tests for P08.17 Intro to Compressible Flow. They assert gas dynamics, not magic numbers:
the ISA reproduces its published anchor values and the hydrostatic equation; a^2 = gamma p / rho;
the isentropic ratios are consistent with each other and with the ideal-gas law at stagnation, are
inverted exactly, and agree with an independent numerical integration of Euler's equation along a
streamline; the stagnation pressure rise tends to the Bernoulli dynamic pressure as M -> 0, with the
known M^2/4 correction; the compressible pitot formula returns the speed it was given; the normal
shock relations conserve mass, momentum and energy, are the identity at M = 1 and lose total
pressure above it; and the Rayleigh pitot formula agrees with shock-then-isentropic compression,
matches the subsonic formula at M = 1, and is inverted exactly.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy.integrate import solve_ivp

from conftest import load_impl

impl = load_impl(__file__)

g = impl.GAMMA


# ---------------------------------------------------------------- atmosphere and sound speed

def test_isa_reproduces_published_anchor_values():
    # ISA anchor values (published data, ICAO / U.S. Standard Atmosphere 1976)
    T, p, rho = impl.isa(0.0)
    assert (T, p) == (288.15, 101325.0) and rho == pytest.approx(1.225, abs=1e-4)
    for h, p_ref in ((11000.0, 22632.06), (20000.0, 5474.889), (32000.0, 868.019), (47000.0, 110.906)):
        assert impl.isa(h)[1] == pytest.approx(p_ref, rel=5e-6)            # tabulated to 6-7 figures
    assert impl.isa(11000.0)[2] == pytest.approx(0.36392, rel=1e-4)
    assert impl.isa(15000.0)[0] == pytest.approx(216.65)                    # isothermal layer
    with pytest.raises(ValueError):
        impl.isa(50000.0)


def test_isa_is_hydrostatic_and_continuous():
    dh = 0.5
    for h in (1000.0, 11000.0 - 3.0, 15000.0, 25000.0, 40000.0):
        dpdh = (impl.isa(h + dh)[1] - impl.isa(h - dh)[1]) / (2 * dh)
        assert dpdh == pytest.approx(-impl.isa(h)[2] * impl.G0, rel=1e-6)
    for hb in (11000.0, 20000.0, 32000.0):
        a, b = impl.isa(hb - 1e-6), impl.isa(hb + 1e-6)
        npt.assert_allclose(a, b, rtol=1e-8)


def test_speed_of_sound():
    T, p, rho = impl.isa(0.0)
    a = float(impl.speed_of_sound(T))
    assert a == pytest.approx(340.294, abs=1e-3)                             # published ISA sea-level value
    assert a ** 2 == pytest.approx(g * p / rho, rel=1e-12)                   # a^2 = gamma p / rho
    assert float(impl.speed_of_sound(4 * T)) == pytest.approx(2 * a)        # ~ sqrt(T)
    assert float(impl.speed_of_sound(300.0, 5 / 3, 2077.1)) == pytest.approx(np.sqrt(5 / 3 * 2077.1 * 300.0))  # helium


# ---------------------------------------------------------------- isentropic relations

def test_isentropic_ratios_are_consistent_and_invert():
    M = np.array([0.0, 0.1, 0.5, 1.0, 2.0, 3.5])
    T0, p0, r0 = impl.T0_over_T(M), impl.p0_over_p(M), impl.rho0_over_rho(M)
    npt.assert_allclose(p0, r0 * T0, rtol=1e-13)                            # p = rho R T at stagnation too
    npt.assert_allclose(p0, T0 ** (g / (g - 1)), rtol=1e-13)
    npt.assert_allclose(impl.mach_from_p0_over_p(p0), M, atol=1e-12)        # inversion
    assert float(impl.p0_over_p(1.0)) == pytest.approx(((g + 1) / 2) ** (g / (g - 1)))


def test_isentropic_pressure_agrees_with_integrating_euler_along_a_streamline():
    # dp/dV = -rho V with p / rho^gamma constant, from the free stream (V, p) down to rest
    T, p, rho = 250.0, 50000.0, 50000.0 / (impl.R_AIR * 250.0)
    for M in (0.3, 0.8, 1.5):
        V = M * float(impl.speed_of_sound(T))
        K = p / rho ** g
        sol = solve_ivp(lambda v, y: [-(y[0] / K) ** (1 / g) * v], (V, 0.0), [p], rtol=1e-11, atol=1e-6)
        assert sol.y[0, -1] / p == pytest.approx(float(impl.p0_over_p(M)), rel=1e-8)


def test_low_mach_limit_recovers_bernoulli():
    for M in (1e-3, 1e-2, 0.05):
        cr = float(impl.pressure_coefficient_ratio(M))
        assert cr == pytest.approx(1 + M ** 2 / 4 + (2 - g) * M ** 4 / 24, rel=1e-8)
    assert float(impl.pressure_coefficient_ratio(1e-4)) == pytest.approx(1.0, abs=1e-8)
    # the airspeed error is about M^2 / 8 and grows monotonically
    M = np.linspace(0.01, 0.9, 50)
    e = impl.bernoulli_error(M)
    assert np.all(np.diff(e) > 0) and e[0] == pytest.approx(1e-4 / 8, rel=1e-3)


def test_one_percent_mach_number_and_energy_equation():
    Me = impl.mach_for_bernoulli_error(0.01)
    assert float(impl.bernoulli_error(Me)) == pytest.approx(0.01, rel=1e-10)
    assert 0.25 < Me < 0.32                                                # the "M < 0.3" rule of thumb
    T, V = 268.65, 250.0
    M = V / float(impl.speed_of_sound(T))
    assert impl.stagnation_temperature(T, V) == pytest.approx(T * float(impl.T0_over_T(M)), rel=1e-12)


def test_compressible_pitot_returns_the_true_airspeed():
    for h, V in ((0.0, 50.0), (3000.0, 60.0), (9000.0, 250.0)):
        T, p, rho = impl.isa(h)
        M = V / float(impl.speed_of_sound(T))
        qc = p * (float(impl.p0_over_p(M)) - 1)
        Vc, Mc = impl.subsonic_pitot_airspeed(qc, p, T)
        assert Vc == pytest.approx(V, rel=1e-12) and Mc == pytest.approx(M, rel=1e-12)
        assert impl.bernoulli_airspeed(qc, rho) > V                         # Bernoulli always reads high
    assert impl.bernoulli_airspeed(0.0, 1.2) == 0.0


# ---------------------------------------------------------------- normal shock and Rayleigh pitot

def test_normal_shock_conserves_mass_momentum_and_energy():
    T1, p1 = 216.65, 22632.0
    rho1 = p1 / (impl.R_AIR * T1)
    for M1 in (1.2, 2.0, 4.0):
        s = impl.normal_shock(M1)
        u1 = M1 * float(impl.speed_of_sound(T1))
        p2, rho2, T2 = p1 * s["p2_p1"], rho1 * s["rho2_rho1"], T1 * s["T2_T1"]
        u2 = s["M2"] * float(impl.speed_of_sound(T2))
        assert rho1 * u1 == pytest.approx(rho2 * u2, rel=1e-12)                       # mass
        assert p1 + rho1 * u1 ** 2 == pytest.approx(p2 + rho2 * u2 ** 2, rel=1e-12)   # momentum
        assert impl.CP_AIR * T1 + u1 ** 2 / 2 == pytest.approx(impl.CP_AIR * T2 + u2 ** 2 / 2, rel=1e-12)  # energy
        assert p2 == pytest.approx(rho2 * impl.R_AIR * T2, rel=1e-12)                  # equation of state
        assert s["M2"] < 1 and s["p02_p01"] < 1                                         # subsonic behind; entropy rises


def test_normal_shock_limits():
    s = impl.normal_shock(1.0)
    for k in ("M2", "p2_p1", "rho2_rho1", "T2_T1", "p02_p01"):
        assert s[k] == pytest.approx(1.0, abs=1e-12)
    strong = impl.normal_shock(1e3)
    assert strong["rho2_rho1"] == pytest.approx((g + 1) / (g - 1), rel=1e-5)          # density ratio limit 6
    assert strong["M2"] == pytest.approx(np.sqrt((g - 1) / (2 * g)), rel=1e-5)
    losses = [impl.normal_shock(M)["p02_p01"] for M in (1.1, 1.5, 2.0, 3.0)]
    assert np.all(np.diff(losses) < 0)
    with pytest.raises(ValueError):
        impl.normal_shock(0.9)


def test_rayleigh_pitot_is_shock_then_isentropic_and_continuous_at_mach_one():
    for M1 in (1.0, 1.3, 2.0, 5.0):
        s = impl.normal_shock(M1)
        assert float(impl.rayleigh_pitot(M1)) == pytest.approx(s["p2_p1"] * float(impl.p0_over_p(s["M2"])), rel=1e-12)
    assert float(impl.rayleigh_pitot(1.0)) == pytest.approx(float(impl.p0_over_p(1.0)), rel=1e-14)
    # the probe reads less than isentropic compression would give (the shock loses total pressure)
    assert float(impl.rayleigh_pitot(2.0)) < float(impl.p0_over_p(2.0))


def test_mach_from_pitot_inverts_both_branches():
    for M in (0.2, 0.7, 0.999, 1.0, 1.001, 1.6, 2.0, 6.0):
        r = float(impl.rayleigh_pitot(M)) if M > 1 else float(impl.p0_over_p(M))
        assert impl.mach_from_pitot(r) == pytest.approx(M, rel=1e-10)
    rs = [float(impl.rayleigh_pitot(M)) for M in np.linspace(1.0, 4.0, 20)]
    assert np.all(np.diff(rs) > 0)                                         # monotone, so the root is unique
    with pytest.raises(ValueError):
        impl.mach_from_pitot(0.5)


def test_results_do_not_depend_on_pressure_units():
    # ratios only: the same reading in kPa and in Pa gives the same Mach number and airspeed
    T, p = 268.65, 70108.5
    qc = 1650.0
    V_pa, _ = impl.subsonic_pitot_airspeed(qc, p, T)
    V_kpa, _ = impl.subsonic_pitot_airspeed(qc / 1e3, p / 1e3, T)
    assert V_pa == pytest.approx(V_kpa, rel=1e-13)


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "compressible.png").stat().st_size > 10_000
