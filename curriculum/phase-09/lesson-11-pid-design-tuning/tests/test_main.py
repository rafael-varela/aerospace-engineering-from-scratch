"""Public tests for P09.11 PID Design & Tuning.

They assert control properties, not magic numbers: the plant transfer function equals C(sI - A)^-1 B of the
state-space model; the PID and closed-loop transfer functions equal their defining expressions pointwise;
the sampled linear loop without saturation converges (first order in dt) to the transfer-function
prediction computed independently by ``scipy.signal.step``; the ultimate gain puts roots on the j-omega
axis; the three-pole placement makes those poles exact roots; the discrete PID does what its equations
say (integration, derivative kick, clamping); anti-windup changes nothing when the actuator never
saturates and reduces overshoot when it does.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy import signal

from conftest import load_impl

impl = load_impl(__file__)

P = impl.pitch_plant()
TF_ = 0.02
GAINS = (2.4, 5.8, 0.46, TF_)                      # close to the lesson's placed design


def _ss_eval(A, B, C, s):
    return (C @ np.linalg.solve(s * np.eye(len(A)) - A, B))[0, 0]


# ---------------------------------------------------------------- transfer functions

def test_plant_tf_equals_state_space_pointwise():
    num, den = impl.pitch_plant_tf(P)
    A, B, C = impl.pitch_state_space(P)
    for s in (0.3 + 1j, 2j, -0.7 + 4j, 15.0 + 30j):
        assert np.polyval(num, s) / np.polyval(den, s) == pytest.approx(_ss_eval(A, B, C, s), rel=1e-10)
    # poles of the TF are the eigenvalues of A (no cancellation in this model)
    npt.assert_allclose(np.sort_complex(np.roots(den)), np.sort_complex(np.linalg.eigvals(A)), atol=1e-9)


@pytest.mark.filterwarnings("ignore::scipy.signal.BadCoefficients")
def test_plant_tf_matches_scipy_realisation():
    A, B, C = impl.pitch_state_space(P)
    ref = signal.lti(A, B, C, np.zeros((1, 1))).to_tf()
    num, den = impl.pitch_plant_tf(P)
    w = np.logspace(-1, 2, 30)
    _, H_ref = signal.freqs(ref.num, ref.den, worN=w)
    _, H = signal.freqs(num, den, worN=w)
    npt.assert_allclose(H, H_ref, rtol=1e-9)


def test_pid_tf_is_kp_plus_ki_over_s_plus_filtered_derivative():
    Kp, Ki, Kd, Tf = 1.7, 3.2, 0.4, 0.05
    num, den = impl.pid_tf(Kp, Ki, Kd, Tf)
    for s in (0.5j, 1 + 2j, 40j):
        direct = Kp + Ki / s + Kd * s / (Tf * s + 1)
        assert np.polyval(num, s) / np.polyval(den, s) == pytest.approx(direct, rel=1e-12)
    # Tf -> 0: the ideal PID Kp + Ki/s + Kd s at moderate frequency
    num0, den0 = impl.pid_tf(Kp, Ki, Kd, 1e-9)
    s = 3j
    assert np.polyval(num0, s) / np.polyval(den0, s) == pytest.approx(Kp + Ki / s + Kd * s, rel=1e-7)


@pytest.mark.parametrize("b,dm", [(1.0, False), (1.0, True), (0.0, True), (0.5, True)])
def test_closed_loop_tf_is_P_Cr_over_1_plus_PC(b, dm):
    num_p, den_p = impl.pitch_plant_tf(P)
    Kp, Ki, Kd, Tf = GAINS
    num, den = impl.closed_loop_tf(num_p, den_p, GAINS, b=b, d_on_measurement=dm)
    for s in (0.4 + 0.8j, 3j, -1 + 5j):
        Ps = np.polyval(num_p, s) / np.polyval(den_p, s)
        C = Kp + Ki / s + Kd * s / (Tf * s + 1)
        Cr = b * Kp + Ki / s + (0.0 if dm else Kd * s / (Tf * s + 1))
        assert np.polyval(num, s) / np.polyval(den, s) == pytest.approx(Ps * Cr / (1 + Ps * C), rel=1e-9)
    # integral action: unit DC gain from reference to pitch whatever b and the derivative input
    assert np.polyval(num, 0.0) / np.polyval(den, 0.0) == pytest.approx(1.0, rel=1e-12)


# ---------------------------------------------------------------- the independent check: linear loop vs TF

@pytest.mark.parametrize("b,dm", [(1.0, True), (0.0, True), (1.0, False)])
def test_linear_loop_without_saturation_matches_the_transfer_function(b, dm):
    num_p, den_p = impl.pitch_plant_tf(P)
    ref = np.radians(10.0)
    sim = impl.simulate_pitch_loop(P, GAINS, ref, 4.0, 1e-3, b=b, d_on_measurement=dm)
    num, den = impl.closed_loop_tf(num_p, den_p, GAINS, b=b, d_on_measurement=dm)
    _, y_ref = signal.step(signal.lti(num, den), T=sim["t"])          # independent oracle
    # dt = 1 ms against loop poles of a few rad/s: the sample-and-hold delay is ~dt/2 and the backward-Euler
    # derivative filter has Tf = 20 dt, so the error is a few 1e-3 of the step (largest with the derivative kick)
    assert np.abs(sim["theta"] - ref * y_ref).max() < 6e-3 * ref
    assert sim["theta"][-1] == pytest.approx(ref, rel=2e-3)


def test_sampling_error_is_first_order_in_dt():
    num_p, den_p = impl.pitch_plant_tf(P)
    num, den = impl.closed_loop_tf(num_p, den_p, GAINS)
    errs = []
    for dt in (4e-3, 2e-3, 1e-3):
        sim = impl.simulate_pitch_loop(P, GAINS, 1.0, 3.0, dt)
        _, y_ref = signal.step(signal.lti(num, den), T=sim["t"])
        errs.append(np.abs(sim["theta"] - y_ref).max())
    assert errs[0] / errs[1] == pytest.approx(2.0, rel=0.15)
    assert errs[1] / errs[2] == pytest.approx(2.0, rel=0.15)


def test_linear_loop_is_superposable_and_zero_input_gives_zero():
    a = impl.simulate_pitch_loop(P, GAINS, np.radians(3.0), 2.0, 2e-3)
    b = impl.simulate_pitch_loop(P, GAINS, np.radians(9.0), 2.0, 2e-3)
    npt.assert_allclose(b["theta"], 3.0 * a["theta"], rtol=1e-9, atol=1e-14)
    z = impl.simulate_pitch_loop(P, GAINS, 0.0, 1.0, 2e-3)
    npt.assert_allclose(z["theta"], 0.0, atol=0.0)


# ---------------------------------------------------------------- tuning

def test_ultimate_gain_puts_the_loop_on_the_imaginary_axis():
    num_p, den_p = impl.pitch_plant_tf(P)
    Ku, wu, Pu = impl.ultimate_gain(num_p, den_p)
    assert Pu == pytest.approx(2 * np.pi / wu)
    roots = np.roots(np.polyadd(den_p, Ku * num_p))
    assert np.min(np.abs(roots - 1j * wu)) < 1e-6 * wu
    # the plant phase is -180 deg there and |K P| = 1
    Pjw = np.polyval(num_p, 1j * wu) / np.polyval(den_p, 1j * wu)
    assert abs(Ku * Pjw) == pytest.approx(1.0, rel=1e-9)
    assert np.angle(Pjw) == pytest.approx(np.pi, abs=1e-7) or np.angle(Pjw) == pytest.approx(-np.pi, abs=1e-7)
    # just below Ku the loop is stable, just above it is not
    assert np.roots(np.polyadd(den_p, 0.98 * Ku * num_p)).real.max() < 0
    assert np.roots(np.polyadd(den_p, 1.02 * Ku * num_p)).real.max() > 0


def test_ultimate_gain_on_a_textbook_plant():
    # 1/(s (s+1)(s+2)): crossing at w = sqrt(2), Ku = 6 (Routh)
    Ku, wu, _ = impl.ultimate_gain(np.array([1.0]), np.polymul([1.0, 0.0], np.polymul([1.0, 1.0], [1.0, 2.0])))
    assert Ku == pytest.approx(6.0, rel=1e-9)
    assert wu == pytest.approx(np.sqrt(2.0), rel=1e-9)
    with pytest.raises(ValueError):
        impl.ultimate_gain(np.array([1.0]), np.array([1.0, 1.0]))     # first order: never reaches the axis


def test_ziegler_nichols_rules():
    zn = impl.ziegler_nichols(5.0, 0.8)
    assert zn["P"] == pytest.approx((2.5, 0.0, 0.0))
    Kp, Ki, Kd = zn["PID"]
    assert Kp == pytest.approx(3.0) and Kp / Ki == pytest.approx(0.4) and Kd / Kp == pytest.approx(0.1)
    Kp, Ki, Kd = zn["PI"]
    assert Kp == pytest.approx(2.25) and Kp / Ki == pytest.approx(0.8 / 1.2) and Kd == 0.0


@pytest.mark.parametrize("poles", [[-4 + 4j, -4 - 4j, -1.5], [-3 + 3j, -3 - 3j, -2.0], [-2.0, -3.0, -5.0]])
def test_three_pole_placement_makes_the_chosen_poles_exact(poles):
    num_p, den_p = impl.pitch_plant_tf(P)
    Kp, Ki, Kd = impl.pid_place_three_poles(num_p, den_p, poles, TF_)
    _, den = impl.closed_loop_tf(num_p, den_p, (Kp, Ki, Kd, TF_))
    roots = np.roots(den)
    for pole in poles:
        assert np.min(np.abs(roots - pole)) < 1e-6 * max(1.0, abs(pole))
        assert abs(np.polyval(den, pole)) < 1e-8 * np.abs(den).max()


def test_three_pole_placement_on_a_design_model_places_every_pole():
    # b/(s(s+a)) with an ideal-ish PID (tiny Tf): the closed loop is third order plus the filter pole
    a, bb = 2.5, 12.0
    poles = [-3 + 3j, -3 - 3j, -3.0]
    Kp, Ki, Kd = impl.pid_place_three_poles(np.array([bb]), np.array([1.0, a, 0.0]), poles, 1e-6)
    c = np.real(np.poly(poles))
    assert Kd == pytest.approx((c[1] - a) / bb, rel=1e-4)
    assert Kp == pytest.approx(c[2] / bb, rel=1e-4)
    assert Ki == pytest.approx(c[3] / bb, rel=1e-4)


# ---------------------------------------------------------------- the discrete controller

def test_pid_step_integrates_and_kicks_as_its_equations_say():
    Kp, Ki, Kd, Tf = 2.0, 4.0, 0.5, 0.05
    dt = 0.01
    state = np.zeros(3)
    u, v, state = impl.pid_step(state, 1.0, 0.0, (Kp, Ki, Kd, Tf), dt, d_on_measurement=False)
    # first sample after a unit step in the error: P term Kp, no integral yet, derivative kick Kd/(Tf + dt)
    assert v == pytest.approx(Kp + Kd / (Tf + dt))
    assert state[0] == pytest.approx(Ki * dt)
    u2, _, _ = impl.pid_step(np.zeros(3), 1.0, 0.0, (Kp, Ki, Kd, Tf), dt, d_on_measurement=True)
    assert u2 == pytest.approx(Kp)                                    # derivative on measurement: no kick
    u3, _, s3 = impl.pid_step(np.zeros(3), 0.0, 0.0, (Kp, Ki, Kd, Tf), dt)
    assert u3 == 0.0 and np.all(s3 == 0.0)


def test_output_respects_the_limit_and_clamping_freezes_the_integrator():
    lim = 0.1
    u, v, s = impl.pid_step(np.array([0.05, 0.0, 0.0]), 1.0, 0.0, GAINS, 0.01, u_limit=lim, anti_windup="clamp")
    assert u == pytest.approx(lim) and v > lim
    assert s[0] == pytest.approx(0.05)                                # saturated and pushing: integrator frozen
    u, v, s = impl.pid_step(np.array([0.05, 0.0, 0.0]), 1.0, 0.0, GAINS, 0.01, u_limit=lim, anti_windup="none")
    assert s[0] == pytest.approx(0.05 + GAINS[1] * 0.01)
    sim = impl.simulate_pitch_loop(P, GAINS, np.radians(20.0), 4.0, 2e-3, u_limit=P["delta_max"], anti_windup="none")
    assert np.abs(sim["u"]).max() <= P["delta_max"] + 1e-15


@pytest.mark.parametrize("scheme", ["clamp", "backcalc"])
def test_anti_windup_is_inactive_when_the_actuator_never_saturates(scheme):
    a = impl.simulate_pitch_loop(P, GAINS, np.radians(2.0), 3.0, 2e-3, u_limit=P["delta_max"], anti_windup="none")
    assert np.abs(a["v"]).max() < P["delta_max"]                      # small step: never saturates
    b = impl.simulate_pitch_loop(P, GAINS, np.radians(2.0), 3.0, 2e-3, u_limit=P["delta_max"], anti_windup=scheme)
    npt.assert_allclose(b["theta"], a["theta"], atol=1e-15)


@pytest.mark.parametrize("scheme", ["clamp", "backcalc"])
def test_anti_windup_reduces_overshoot_under_saturation(scheme):
    ref = np.radians(20.0)
    wound = impl.simulate_pitch_loop(P, GAINS, ref, 8.0, 2e-3, u_limit=P["delta_max"], anti_windup="none")
    fixed = impl.simulate_pitch_loop(P, GAINS, ref, 8.0, 2e-3, u_limit=P["delta_max"], anti_windup=scheme)
    assert np.abs(wound["v"]).max() > 2 * P["delta_max"]              # the case really saturates
    os_wound = impl.step_metrics(wound["t"], wound["theta"], ref)["overshoot_pct"]
    os_fixed = impl.step_metrics(fixed["t"], fixed["theta"], ref)["overshoot_pct"]
    assert os_fixed < 0.6 * os_wound
    assert np.abs(fixed["I"]).max() < np.abs(wound["I"]).max()        # the integrator wound up less
    assert fixed["theta"][-1] == pytest.approx(ref, rel=1e-2)         # and still reaches the target


def test_integral_action_removes_a_trim_disturbance():
    pd = (GAINS[0], 0.0, GAINS[2], TF_)
    ref = np.radians(5.0)
    d = np.radians(-2.0)
    with_i = impl.simulate_pitch_loop(P, GAINS, ref, 15.0, 2e-3, disturbance=d, t_disturbance=3.0)
    without = impl.simulate_pitch_loop(P, pd, ref, 15.0, 2e-3, disturbance=d, t_disturbance=3.0)
    assert abs(with_i["theta"][-1] - ref) < 1e-4 * ref
    # PD only: steady error = the disturbance moment balanced by Kp e: e = -d * (dc gain)/... must be clearly non-zero
    assert abs(without["theta"][-1] - ref) > 0.05 * ref


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "pid_pitch.png").stat().st_size > 1000
