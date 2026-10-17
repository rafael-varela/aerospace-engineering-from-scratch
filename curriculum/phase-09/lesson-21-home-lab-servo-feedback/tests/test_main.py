"""Public tests for P09.21 Home Lab: Servo Feedback Control.

They assert: the step-response models' limits; identification recovering known parameters from noiseless
and from seeded noisy data (log-linear and Levenberg-Marquardt), and the dead zone from several steps;
the sampled-data metrics against the closed-form second-order and first-order results; the P and PI
designs meeting their targets; the closed-loop prediction against scipy.signal.step (the oracle); the
sampled-loop simulator reducing to the prediction when sampling is fast and ideal; the duty (command)
limit never exceeded; quantisation, noise and the one-sample delay behaving as stated; and the lesson's
independent check (the identified model's predicted closed-loop rise time, overshoot and settling against
the SYNTHETIC closed-loop runs). The SYNTHETIC CSV files are regenerated identically from the seed (only
into tmp_path) and the shipped files in projects/p09-servo-control/data/synthetic/ equal the generator's.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
from pathlib import Path

import numpy as np
import pytest
from scipy import signal

from conftest import load_impl

impl = load_impl(__file__)
ROOT = Path(__file__).resolve().parents[4]
SHIPPED = ROOT / "projects" / "p09-servo-control" / "data" / "synthetic"


def sampled(fn, t_end, dt, *args, t_pre=0.05):
    t = np.arange(-t_pre, t_end + dt / 2, dt)
    return t, fn(t, *args)


# ---------------------------------------------------------------- models and identification

def test_step_response_models_limits():
    K, tau, L, u = 0.9, 0.06, 0.012, 20.0
    t = np.array([-0.1, 0.0, L, L + tau, L + 40 * tau])
    y = impl.first_order_step_response(t, K, tau, L, u)
    assert y[0] == 0 and y[1] == 0 and y[2] == 0                       # nothing before the dead time
    assert y[3] == pytest.approx(K * u * (1 - np.exp(-1)))              # 63.2 % one time constant later
    assert y[4] == pytest.approx(K * u, rel=1e-12)                      # final value K u
    K, tau, L, u = 540.0, 0.045, 0.005, 0.3
    tt = np.array([0.0, L, 2.0, 3.0])
    y = impl.integrator_step_response(tt, K, tau, L, u)
    assert y[0] == 0 and y[1] == 0
    assert (y[3] - y[2]) == pytest.approx(K * u, rel=1e-12)             # tail slope K u
    assert 3.0 - y[3] / (K * u) == pytest.approx(tau + L, rel=1e-9)     # tail lags the ramp by tau + L
    # both closed forms against the oracle's step of the same transfer function
    tg = np.linspace(0, 0.5, 501)
    _, ys = signal.step(signal.lti([K], [tau, 1.0, 0.0]), T=tg)
    np.testing.assert_allclose(impl.integrator_step_response(tg, K, tau, 0.0, 1.0), ys, rtol=1e-9, atol=1e-9)
    _, ys = signal.step(signal.lti([0.93], [0.055, 1.0]), T=tg)
    np.testing.assert_allclose(impl.first_order_step_response(tg, 0.93, 0.055, 0.0, 1.0), ys, atol=1e-12)


@pytest.mark.parametrize("K,tau,L,u,dt", [(0.93, 0.055, 0.0126, 20.0, 0.002), (2.0, 0.3, 0.0, 1.0, 0.01)])
def test_first_order_loglinear_recovers_noiseless(K, tau, L, u, dt):
    t, y = sampled(impl.first_order_step_response, 20 * tau + L, dt, K, tau, L, u)
    f = impl.fit_first_order_loglinear(t, y + 5.0, u)                    # an offset must not matter
    assert f["y0"] == pytest.approx(5.0)
    assert f["K"] == pytest.approx(K, rel=1e-6)
    assert f["tau"] == pytest.approx(tau, rel=1e-6)
    assert f["delay"] == pytest.approx(L, abs=1e-6)


@pytest.mark.parametrize("K,tau,L,u,dt", [(540.0, 0.045, 0.0, 0.6, 0.002), (120.0, 0.2, 0.01, 1.0, 0.01)])
def test_integrator_loglinear_recovers_noiseless(K, tau, L, u, dt):
    t, y = sampled(impl.integrator_step_response, 20 * tau + L, dt, K, tau, L, u)
    f = impl.fit_integrator_loglinear(t, y, u)
    assert f["K"] == pytest.approx(K, rel=1e-6)
    assert f["tau"] == pytest.approx(tau, rel=1e-4)
    assert f["delay"] == pytest.approx(L, abs=1e-4 * tau)


def test_least_squares_converges_from_a_poor_start():
    K, tau, L, u = 540.0, 0.045, 0.004, 0.45
    t, y = sampled(impl.integrator_step_response, 0.5, 0.01, K, tau, L, u)
    model = lambda tt, p: impl.integrator_step_response(tt, p[0], p[1], p[2], u)  # noqa: E731
    r = impl.refine_least_squares(model, t, y, np.array([0.7 * K, 1.4 * tau, 0.0]), lower=np.zeros(3))
    np.testing.assert_allclose(r["p"], [K, tau, L], rtol=1e-6, atol=1e-8)
    assert r["rms"] < 1e-8


def test_identification_on_seeded_noisy_data():
    # integrator plant with wiper noise 0.08 deg and 12-bit quantisation, as in the synthetic bench
    num, den = np.array([540.0]), np.array([0.045, 1.0, 0.0])
    sim = impl.simulate_open_loop(num, den, 0.6, 0.01, 0.5, quant=300 / 4096, noise_sd=0.08, t_pre=0.1, seed=7)
    g = impl.fit_integrator_loglinear(sim["t"], sim["y"], 0.6)
    model = lambda tt, p: g["y0"] + impl.integrator_step_response(tt, p[0], p[1], p[2], 0.6)  # noqa: E731
    r = impl.refine_least_squares(model, sim["t"], sim["y"], np.array([g["K"], g["tau"], g["delay"]]),
                                  lower=np.zeros(3))
    # tolerance: 4 linearised standard errors (a 4-sigma band fails by chance about 1 time in 16 000)
    assert abs(r["p"][0] - 540.0) < 4 * r["se"][0] + 1e-9
    assert abs(r["p"][1] - 0.045) < 4 * r["se"][1]
    assert r["p"][2] < 4 * r["se"][2] + 1e-12
    assert r["se"][1] < 0.05 * 0.045                                     # tau to better than 5 % (1 sigma)
    assert r["rms"] == pytest.approx(np.sqrt(0.08 ** 2 + (300 / 4096) ** 2 / 12), rel=0.25)


def test_dead_zone_from_several_steps():
    K, udz = 540.0, 0.1
    duties = np.array([0.3, 0.45, 0.6])
    Kf, dzf = impl.fit_dead_zone(duties, K * (duties - udz))
    assert (Kf, dzf) == pytest.approx((K, udz), rel=1e-12)
    with pytest.raises(ValueError):
        impl.fit_dead_zone(np.array([0.3]), np.array([100.0]))


# ---------------------------------------------------------------- metrics

@pytest.mark.parametrize("zeta", [0.2, 0.6])
def test_metrics_match_second_order_closed_form(zeta):
    wn = 18.0
    t, y = sampled(impl.second_order_step, 3.0, 1e-4, zeta, wn)
    m = impl.step_metrics(t, y, y_initial=0.0, y_final=1.0)
    ref = impl.second_order_metrics(zeta, wn)
    assert m["overshoot_pct"] == pytest.approx(100 * np.exp(-np.pi * zeta / np.sqrt(1 - zeta ** 2)), abs=1e-4)
    assert ref["overshoot_pct"] == pytest.approx(m["overshoot_pct"], abs=1e-4)
    assert m["peak_time"] == pytest.approx(np.pi / (wn * np.sqrt(1 - zeta ** 2)), abs=1e-4)
    assert m["rise_time"] == pytest.approx(ref["rise_time"], abs=1e-6)
    assert m["settling_time"] == pytest.approx(ref["settling_time"], abs=1e-6)
    # independent check of the closed form against the oracle
    _, ys = signal.step(signal.lti([wn ** 2], [1, 2 * zeta * wn, wn ** 2]), T=t[t >= 0])
    np.testing.assert_allclose(y[t >= 0], ys, atol=1e-9)


def test_first_order_metrics_and_inverses():
    tau = 0.2
    t, y = sampled(impl.first_order_step_response, 3.0, 1e-4, 1.0, tau, 0.0, 1.0)
    m = impl.step_metrics(t, y, y_initial=0.0, y_final=1.0)
    assert m["rise_time"] == pytest.approx(tau * np.log(9), abs=1e-6)
    assert m["settling_time"] == pytest.approx(tau * np.log(50), abs=1e-6)
    assert m["overshoot_pct"] == 0.0
    for os_ in (1.0, 9.0, 25.0, 60.0):
        assert impl.overshoot_from_zeta(impl.zeta_from_overshoot(os_)) == pytest.approx(os_, rel=1e-12)
    assert impl.overshoot_from_zeta(1.0) == 0.0
    assert impl.zeta_from_overshoot(10.0) == pytest.approx(0.5912, abs=1e-4)


# ---------------------------------------------------------------- design and prediction

@pytest.mark.parametrize("target", [5.0, 20.0])
def test_p_design_meets_overshoot(target):
    K, tau = 540.0, 0.045
    d = impl.design_p_position(K, tau, target)
    p = impl.predict_closed_loop(np.array([K]), np.array([tau, 1.0, 0.0]), d["kp"], 0.0, 15.0, 1.5)
    assert p["metrics"]["overshoot_pct"] == pytest.approx(target, abs=1e-3)
    assert p["metrics"]["y_final"] == pytest.approx(15.0)                  # type 1: no steady-state error
    poles = np.sort_complex(p["poles"])
    assert np.abs(poles[0]) == pytest.approx(d["wn"], rel=1e-9)            # |pole| = wn
    assert -poles[0].real / np.abs(poles[0]) == pytest.approx(d["zeta"], rel=1e-9)


def test_pi_cancellation_meets_rise_time():
    K, tau, tr = 0.93, 0.055, 0.4
    d = impl.design_pi_first_order(K, tau, tr)
    assert d["ki"] / d["kp"] == pytest.approx(1 / tau)                     # the zero sits on the plant pole
    p = impl.predict_closed_loop(np.array([K]), np.array([tau, 1.0]), d["kp"], d["ki"], 20.0, 2.0)
    assert p["metrics"]["rise_time"] == pytest.approx(tr, rel=1e-4)
    assert p["metrics"]["overshoot_pct"] < 1e-9
    assert p["metrics"]["settling_time"] == pytest.approx(d["tau_cl"] * np.log(50), rel=1e-4)


@pytest.mark.parametrize("plant,kp,ki", [((np.array([540.0]), np.array([0.045, 1.0, 0.0])), 0.029, 0.0),
                                         ((np.array([0.93]), np.array([0.055, 1.0])), 0.33, 5.9),
                                         ((np.array([0.93]), np.array([0.055, 1.0])), 2.0, 0.0)])
def test_prediction_against_scipy_oracle(plant, kp, ki):
    num, den = plant
    p = impl.predict_closed_loop(num, den, kp, ki, 1.0, 1.0, dt=1e-4)
    _, ys = signal.step(signal.lti(p["num_cl"], p["den_cl"]), T=p["t"])
    np.testing.assert_allclose(p["y"], ys, atol=1e-9)
    # closed-loop algebra against the unreduced loop evaluated pointwise: T = CG/(1 + CG)
    s = 1j * np.array([0.5, 3.0, 40.0])
    C = kp + ki / s
    G = np.polyval(num, s) / np.polyval(den, s)
    np.testing.assert_allclose(np.polyval(p["num_cl"], s) / np.polyval(p["den_cl"], s), C * G / (1 + C * G),
                               rtol=1e-12)
    # P on a type-0 plant leaves the error 1/(1 + K kp); PI or a type-1 plant leaves none
    dc = p["num_cl"][-1] / p["den_cl"][-1]
    expect = 1.0 if (ki or den[-1] == 0) else num[-1] * kp / (den[-1] + num[-1] * kp)
    assert dc == pytest.approx(expect)


# ---------------------------------------------------------------- the sampled loop

@pytest.mark.parametrize("plant,kp,ki,r", [((np.array([540.0]), np.array([0.045, 1.0, 0.0])), 0.029, 0.0, 15.0),
                                           ((np.array([0.93]), np.array([0.055, 1.0])), 0.33, 5.9, 20.0)])
def test_fast_ideal_simulation_matches_prediction(plant, kp, ki, r):
    num, den = plant
    pred = impl.predict_closed_loop(num, den, kp, ki, r, 1.5)["metrics"]
    sim = impl.simulate_closed_loop(num, den, kp, ki, r, 0.0005, 1.5, delay_samples=0, substeps=4)
    m = impl.step_metrics(sim["t"], sim["y"], ref=r)
    # a zero-order hold at Ts = 0.5 ms acts as a 0.25 ms delay: a few per cent at most on these loops
    assert m["rise_time"] == pytest.approx(pred["rise_time"], rel=0.02)
    assert m["settling_time"] == pytest.approx(pred["settling_time"], rel=0.03)
    assert m["overshoot_pct"] == pytest.approx(pred["overshoot_pct"], abs=0.3)


def test_open_loop_simulation_matches_closed_form():
    num, den = np.array([540.0]), np.array([0.045, 1.0, 0.0])
    sim = impl.simulate_open_loop(num, den, 0.45, 0.01, 0.5, plant_delay=0.005, substeps=20, t_pre=0.05)
    exact = impl.integrator_step_response(sim["t"], 540.0, 0.045, 0.005, 0.45)
    np.testing.assert_allclose(sim["y"], exact, atol=1e-6)             # RK4 at 0.5 ms on a 45 ms lag
    dz = impl.simulate_open_loop(num, den, 0.3, 0.01, 0.5, dz=0.1)
    assert dz["y"][-1] - dz["y"][-2] == pytest.approx(540.0 * 0.2 * 0.01, rel=1e-4)   # slope K (u - u_dz); e^(-0.5/0.045) left
    nothing = impl.simulate_open_loop(num, den, 0.08, 0.01, 0.5, dz=0.1)
    assert np.all(nothing["y"] == 0.0)                                   # inside the dead zone: no motion


@pytest.mark.parametrize("kp,ki,r,lim", [(0.2, 0.0, -30.0, 0.4), (0.05, 2.0, 90.0, 0.5)])
def test_duty_limit_never_exceeded(kp, ki, r, lim):
    num, den = np.array([540.0]), np.array([0.045, 1.0, 0.0])
    sim = impl.simulate_closed_loop(num, den, kp, ki, r, 0.01, 1.5, u_min=-lim, u_max=lim, dz=0.1, dz_comp=0.09,
                                    quant=300 / 4096, noise_sd=0.08, seed=3)
    assert np.max(np.abs(sim["u"])) <= lim + 1e-15
    assert np.max(np.abs(sim["u"])) == pytest.approx(lim)                # the limit is actually reached


def test_quantisation_noise_and_zero_gain():
    num, den = np.array([0.93]), np.array([0.055, 1.0])
    q = 0.066
    a = impl.simulate_closed_loop(num, den, 0.33, 5.9, 20.0, 0.02, 1.0, quant=q, noise_sd=0.06, seed=11)
    b = impl.simulate_closed_loop(num, den, 0.33, 5.9, 20.0, 0.02, 1.0, quant=q, noise_sd=0.06, seed=11)
    c = impl.simulate_closed_loop(num, den, 0.33, 5.9, 20.0, 0.02, 1.0, quant=q, noise_sd=0.06, seed=12)
    np.testing.assert_allclose(a["y"] / q, np.round(a["y"] / q), atol=1e-9)   # multiples of one count
    np.testing.assert_array_equal(a["y"], b["y"])                              # seeded: reproducible
    assert not np.array_equal(a["y"], c["y"])
    z = impl.simulate_closed_loop(num, den, 0.0, 0.0, 20.0, 0.02, 1.0)
    assert np.all(z["y"] == 0.0) and np.all(z["u"] == 0.0)                     # no gain, no action


def test_one_sample_delay_shifts_and_destabilises():
    num, den = np.array([540.0]), np.array([0.045, 1.0, 0.0])
    s0 = impl.simulate_closed_loop(num, den, 0.029, 0.0, 15.0, 0.01, 1.0, delay_samples=0)
    s1 = impl.simulate_closed_loop(num, den, 0.029, 0.0, 15.0, 0.01, 1.0, delay_samples=1)
    k0 = int(np.flatnonzero(s0["u"] != 0)[0])
    k1 = int(np.flatnonzero(s1["u"] != 0)[0])
    assert k1 == k0 + 1                                                        # applied one sample later
    m0 = impl.step_metrics(s0["t"], s0["y"])
    m1 = impl.step_metrics(s1["t"], s1["y"])
    assert m1["overshoot_pct"] > m0["overshoot_pct"] + 5                       # the delay costs damping


def test_sampled_redesign_hits_target():
    num, den = np.array([540.0]), np.array([0.045, 1.0, 0.0])
    kw = dict(u_min=-0.6, u_max=0.6, delay_samples=1, dz=0.1, dz_comp=0.1)
    kp = impl.design_p_sampled(num, den, 15.0, 0.01, 10.0, 0.005, 0.05, t_end=1.0, **kw)
    sim = impl.simulate_closed_loop(num, den, kp, 0.0, 15.0, 0.01, 1.0, **kw)
    assert impl.step_metrics(sim["t"], sim["y"])["overshoot_pct"] == pytest.approx(10.0, abs=0.05)
    assert kp < impl.design_p_position(540.0, 0.045, 10.0)["kp"]              # delay forces a lower gain


# ---------------------------------------------------------------- the lesson's independent check

@pytest.fixture(scope="module")
def lab(tmp_path_factory):
    out = tmp_path_factory.mktemp("p09_l21")
    return out, impl.main(out)


def test_identification_recovers_the_synthetic_bench(lab):
    _, res = lab
    m, s = res["motor"], res["servo"]
    assert m["K"] == pytest.approx(impl.MOTOR_TRUE["K"], rel=0.01)
    assert m["tau"] == pytest.approx(impl.MOTOR_TRUE["tau"], rel=0.03)
    assert m["dead_zone"] == pytest.approx(impl.MOTOR_TRUE["dz_pos"], abs=0.01)
    assert m["delay"] < 0.002                                                   # the motor bench has none
    # the servo's dead band biases K and tau low; the identified delay is the frame latency
    assert s["K"] == pytest.approx(impl.SERVO_TRUE["K"], rel=0.04)
    assert s["tau"] == pytest.approx(impl.SERVO_TRUE["tau"], rel=0.10)
    assert s["delay"] == pytest.approx(impl.SERVO_TRUE["delay"], abs=0.002)
    # the firmware gains are the designs, rounded (the motor gains rounded down)
    assert impl.MOTOR_KP == pytest.approx(m["design"]["kp"], rel=0.05)
    assert impl.MOTOR_KP_REDESIGN == pytest.approx(m["kp_sampled"], rel=0.05)
    assert impl.SERVO_KP == pytest.approx(s["design"]["kp"], rel=0.05)
    assert impl.SERVO_KI == pytest.approx(s["design"]["ki"], rel=0.05)


def test_model_predicted_metrics_against_synthetic_closed_loop(lab):
    _, res = lab
    runs = res["motor"]["runs"] + res["servo"]["runs"]
    for run in runs:
        assert np.max(np.abs(run["sim"]["u"])) <= (impl.MOTOR_DUTY_LIMIT if run["ki"] == 0 else impl.SERVO_CMD_LIMIT)
        sim, meas = run["sim_metrics"], run["meas"]
        # the sampled model (identified plant + Ts + delay + limit + dead zone + quantisation) predicts the
        # synthetic runs; the residual is noise, the asymmetric dead zone and the servo's slew and dead band
        assert meas["rise_time"] == pytest.approx(sim["rise_time"], rel=0.10)
        assert meas["overshoot_pct"] == pytest.approx(sim["overshoot_pct"], abs=3.5)
        assert meas["settling_time"] == pytest.approx(sim["settling_time"], rel=0.15)
    # and the continuous model alone misses the delay: the motor's measured overshoot is about double
    r1 = res["motor"]["runs"][0]
    assert r1["meas"]["overshoot_pct"] > 1.5 * r1["pred"]["metrics"]["overshoot_pct"]


def test_main_writes_labelled_outputs(lab):
    out, _ = lab
    assert (out / "servo_feedback_lab.png").stat().st_size > 10_000
    for name in impl.SYNTHETIC_FILES:
        assert "synthetic" in name
        lines = (out / name).read_text().splitlines()
        assert lines[0].startswith("# SYNTHETIC DATA") and "not a measurement" in lines[0]
        assert all(line.endswith(",synthetic") for line in lines[2:])
    assert len(impl.read_runs_csv(out / impl.SYNTHETIC_FILES[1])) == len(impl.MOTOR_CLOSED_RUNS)


# ---------------------------------------------------------------- the SYNTHETIC dataset

def test_synthetic_csvs_regenerate_identically(tmp_path):
    a = impl.write_synthetic_csvs(tmp_path / "a", seed=impl.SEED)
    b = impl.write_synthetic_csvs(tmp_path / "b", seed=impl.SEED)
    c = impl.write_synthetic_csvs(tmp_path / "c", seed=impl.SEED + 1)
    for pa, pb, pc in zip(a, b, c):
        assert pa.read_bytes() == pb.read_bytes()
        assert pa.read_bytes() != pc.read_bytes()


def test_shipped_synthetic_files_equal_the_generator(tmp_path):
    fresh = impl.write_synthetic_csvs(tmp_path, seed=impl.SEED)
    for p in fresh:
        shipped = SHIPPED / p.name
        assert shipped.exists(), f"missing shipped file {shipped}"
        assert shipped.read_bytes() == p.read_bytes()
