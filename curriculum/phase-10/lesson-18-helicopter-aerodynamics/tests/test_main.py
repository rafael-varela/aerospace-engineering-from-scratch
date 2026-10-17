"""Public tests for P10.18 Helicopter Aerodynamics. They assert physics, not magic numbers:
momentum theory conserves momentum and energy in hover; blade-element thrust with uniform inflow
matches its closed form and strip integration, and its induced power equals momentum theory;
the figure of merit is below one and rises with disk loading; Glauert's inflow tends to the hover
value as mu -> 0 and to C_T/(2 mu) at high speed; the profile-power factor is exactly 1 + 3 mu^2
without radial flow and near 4.65 with it; the forward-flight power satisfies the energy balance
C_T lambda = induced + parasite; the power curve has its bucket at the minimum-power speed and the
maximum-range speed is faster; autorotation follows from W V_d = P.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import pytest

from conftest import load_impl

impl = load_impl(__file__)

RHO = 1.225


@pytest.fixture(scope="module")
def rot():
    return impl.Rotor()


# ---------------------------------------------------------------- hover

def test_hover_momentum_conserves_momentum_and_energy(rot):
    T = rot.W
    h = impl.hover_momentum(T, rot.A)
    mdot = RHO * rot.A * h["v_h"]
    assert mdot * h["v_wake"] == pytest.approx(T, rel=1e-12)                     # momentum
    assert h["P_ideal"] == pytest.approx(0.5 * mdot * h["v_wake"] ** 2, rel=1e-12)  # energy
    assert RHO * h["A_wake"] * h["v_wake"] == pytest.approx(mdot)                 # continuity: wake contracts to A/2


def test_coefficients_and_lambda_h(rot):
    CT = impl.thrust_coefficient(rot.W, rot)
    v_h = impl.hover_momentum(rot.W, rot.A)["v_h"]
    assert np.sqrt(CT / 2) == pytest.approx(v_h / rot.tip_speed, rel=1e-12)
    assert float(impl.power_from_coefficient(CT ** 1.5 / np.sqrt(2), rot)) == pytest.approx(
        impl.hover_momentum(rot.W, rot.A)["P_ideal"], rel=1e-12)


def test_figure_of_merit_bounds_and_trend(rot):
    T = rot.W
    P_ideal = impl.hover_momentum(T, rot.A)["P_ideal"]
    assert impl.figure_of_merit(T, P_ideal, rot.A) == pytest.approx(1.0)
    fms = [impl.hover_power(rot, m * 9.81)["FM"] for m in (600, 900, 1200, 1500)]
    assert all(0 < f < 1 for f in fms)
    assert np.all(np.diff(fms) > 0)                    # profile power is a smaller share at higher loading
    hp = impl.hover_power(impl.Rotor(cd0=0.0, kappa=1.0))
    assert hp["FM"] == pytest.approx(1.0, rel=1e-12)  # no profile drag, uniform inflow: ideal


def test_blade_element_closed_form_vs_strips():
    for th0, tw, lam, sig in [(0.2, 0.0, 0.03, 0.05), (0.25, -0.14, 0.04, 0.07), (0.1, -0.3, 0.0, 0.03)]:
        cf = impl.be_thrust_coefficient(th0, lam, sig, 5.7, tw)
        st = impl.be_thrust_strips(th0, lam, sig, 5.7, tw, 400)
        assert st == pytest.approx(cf, rel=1e-5, abs=1e-12)
    # linear in collective, decreasing in inflow
    a = impl.be_thrust_coefficient(0.2, 0.03, 0.05, 5.7)
    b = impl.be_thrust_coefficient(0.3, 0.03, 0.05, 5.7)
    c = impl.be_thrust_coefficient(0.4, 0.03, 0.05, 5.7)
    assert b - a == pytest.approx(c - b)
    assert impl.be_thrust_coefficient(0.2, 0.05, 0.05, 5.7) < a


def test_collective_reproduces_thrust_and_induced_power_matches_momentum(rot):
    CT = impl.thrust_coefficient(rot.W, rot)
    th0, th75 = impl.hover_collective(CT, rot)
    lam = np.sqrt(CT / 2)
    assert impl.be_thrust_coefficient(th0, lam, rot.sigma, rot.a, rot.theta_tw) == pytest.approx(CT, rel=1e-12)
    assert th75 == pytest.approx(th0 + 0.75 * rot.theta_tw)
    # uniform inflow: blade-element induced power = lambda C_T = momentum C_T^(3/2)/sqrt(2)
    r = (np.arange(2000) + 0.5) / 2000
    dCT = 0.5 * rot.sigma * rot.a * ((th0 + rot.theta_tw * r) * r ** 2 - lam * r)
    CPi = np.sum(lam * dCT) / 2000
    assert CPi == pytest.approx(CT ** 1.5 / np.sqrt(2), rel=1e-5)


def test_hover_power_parts(rot):
    hp = impl.hover_power(rot)
    assert hp["P"] == pytest.approx(hp["P_induced"] + hp["P_profile"])
    assert hp["CP0"] == pytest.approx(rot.sigma * rot.cd0 / 8)
    assert hp["P_induced"] == pytest.approx(rot.kappa * impl.hover_momentum(rot.W, rot.A)["P_ideal"], rel=1e-12)
    assert 0.5 < hp["FM"] < 0.8                                                   # a light helicopter's range


# ---------------------------------------------------------------- forward flight inflow

def test_glauert_inflow_hover_limit_and_residual(rot):
    CT = impl.thrust_coefficient(rot.W, rot)
    lam, lam_i, _ = impl.glauert_inflow(0.0, CT, 0.0)
    assert lam_i == pytest.approx(np.sqrt(CT / 2), rel=1e-12)
    for mu, al in [(1e-4, 0.0), (0.05, 0.02), (0.2, 0.05), (0.35, -0.03)]:
        lam, lam_i, it = impl.glauert_inflow(mu, CT, al)
        assert abs(lam - mu * np.tan(al) - CT / (2 * np.sqrt(mu ** 2 + lam ** 2))) < 1e-13
        assert it < 30
    assert impl.glauert_inflow(1e-4, CT, 0.0)[1] == pytest.approx(np.sqrt(CT / 2), rel=1e-4)


def test_glauert_inflow_high_speed_limit(rot):
    CT = impl.thrust_coefficient(rot.W, rot)
    errs = []
    for mu in (0.1, 0.2, 0.4):
        lam_i = impl.glauert_inflow(mu, CT, 0.0)[1]
        errs.append(abs(lam_i / (CT / (2 * mu)) - 1))
    assert errs[-1] < 2e-3 and np.all(np.diff(errs) < 0)                         # v -> T / (2 rho A V)
    li = [impl.glauert_inflow(m, CT, 0.0)[1] for m in np.linspace(0, 0.3, 31)]
    assert np.all(np.diff(li) < 0)                                               # inflow falls with speed


# ---------------------------------------------------------------- profile power

def test_profile_power_factor_three_without_radial_flow():
    for mu in (0.05, 0.15, 0.3):
        k = (impl.profile_power_integral(mu, radial_flow=False, n=200) - 1) / mu ** 2
        assert 2.95 < k < 3.1                                                    # 1 + 3 mu^2 (|u_T| adds reverse flow)
    assert impl.profile_power_integral(0.0) == pytest.approx(1.0, rel=1e-4)


def test_profile_power_factor_with_radial_flow_near_four_point_six():
    ks = [(impl.profile_power_integral(mu, radial_flow=True, n=200) - 1) / mu ** 2 for mu in (0.1, 0.2, 0.3)]
    assert all(4.4 < k < 4.9 for k in ks)


# ---------------------------------------------------------------- power curve

def test_hover_and_forward_flight_power_agree_at_zero_speed(rot):
    p = impl.forward_flight_power(0.0, rot)
    assert p["P"][0] == pytest.approx(impl.hover_power(rot)["P"], rel=1e-10)


def test_energy_balance_in_forward_flight():
    rot = impl.Rotor(kappa=1.0)
    for V in (10.0, 30.0, 50.0, 70.0):
        p = impl.forward_flight_power(V, rot)
        # with kappa = 1, C_T lambda = (induced) + (parasite): T (V sin a + v) = T v + D V
        assert p["CT"][0] * p["lam"][0] == pytest.approx(p["CPi"][0] + p["CPp"][0], rel=1e-10)
        D = 0.5 * RHO * V ** 2 * rot.f
        assert p["P_parasite"][0] == pytest.approx(D * V, rel=1e-10)


def test_power_curve_bucket_and_speeds(rot):
    bs = impl.best_speeds(rot)
    assert 10 < bs["V_min_power"] < bs["V_max_range"] < 70
    P = lambda v: impl.forward_flight_power(v, rot)["P"][0]
    assert P(bs["V_min_power"]) < P(0.0) and P(bs["V_min_power"]) < P(60.0)
    # tangency: at the max-range speed, dP/dV = P/V
    V = bs["V_max_range"]
    dPdV = (P(V + 0.01) - P(V - 0.01)) / 0.02
    assert dPdV == pytest.approx(P(V) / V, rel=1e-3)
    # a cleaner airframe (smaller f) flies faster for range
    assert impl.best_speeds(impl.Rotor(f=0.6))["V_max_range"] > V


def test_heavier_helicopter_needs_more_power_and_climb_adds_W_Vc(rot):
    heavy = impl.Rotor(mass=1300.0)
    for V in (0.0, 30.0, 60.0):
        assert impl.forward_flight_power(V, heavy)["P"][0] > impl.forward_flight_power(V, rot)["P"][0]
    p0 = impl.forward_flight_power(30.0, rot)["P"][0]
    p1 = impl.forward_flight_power(30.0, rot, climb=2.0)["P"][0]
    T = np.hypot(rot.W, 0.5 * RHO * 30.0 ** 2 * rot.f)
    assert p1 - p0 == pytest.approx(T * 2.0, rel=1e-10)


def test_autorotation_energy_balance(rot):
    au = impl.autorotation(rot)
    bs = impl.best_speeds(rot)
    assert au["Vd_min"] * rot.W == pytest.approx(bs["P_min"], rel=1e-10)
    assert au["Vd_best_glide"] > au["Vd_min"]
    assert au["glide_ratio"] == pytest.approx(au["V_best_glide"] / au["Vd_best_glide"])
    assert 3.0 < au["Vd_min"] < 12.0                                             # metres per second, not tens


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "helicopter_power.png").stat().st_size > 10_000
