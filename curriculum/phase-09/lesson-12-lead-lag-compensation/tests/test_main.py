"""Public tests for P09.12 Lead-Lag Compensation.

They assert loop-shaping properties, not magic numbers: the frequency response equals
``scipy.signal.freqs``; a lead reaches its maximum phase asin((1 - a)/(1 + a)) at the geometric mean of its
zero and pole; a lag has DC gain beta and costs little phase a decade above its zero; the margin finder
reproduces textbook margins and the margins read off ``scipy.signal.bode``; the margins measured on the
compensated loop equal the design targets (the lesson's independent check), Kv meets its target, the
closed-loop step agrees with ``scipy.signal.step`` and the ramp error with 1/Kv from ``scipy.signal.lsim``.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy import signal

from conftest import load_impl

impl = load_impl(__file__)

NG, DG = impl.servo_tf(impl.servo_plant())


def bode_margins(num, den):
    """Margins read off scipy.signal.bode on a dense grid by linear interpolation (independent of impl)."""
    w = np.logspace(-2, 4, 200001)
    w, mag, ph = signal.bode((num, den), w=w)
    i = np.where(np.diff(np.sign(mag)) != 0)[0][0]
    f = mag[i] / (mag[i] - mag[i + 1])
    wc = w[i] + f * (w[i + 1] - w[i])
    pm = 180.0 + ph[i] + f * (ph[i + 1] - ph[i])
    return wc, pm


# ---------------------------------------------------------------- building blocks

def test_freq_response_matches_scipy_freqs():
    w = np.logspace(-1, 3, 50)
    nl, dl = impl.lead_compensator(2.0, 0.3, 0.05)
    n, d = impl.series(nl, dl, NG, DG)
    _, H = signal.freqs(n, d, worN=w)
    npt.assert_allclose(impl.freq_response(n, d, w), H, rtol=1e-12)


@pytest.mark.parametrize("alpha", [0.1, 0.3, 0.6])
def test_lead_reaches_its_maximum_phase_at_the_geometric_mean(alpha):
    T = 0.04
    n, d = impl.lead_compensator(1.0, alpha, T)
    w = np.logspace(0, 3, 200001)
    ph = np.degrees(np.angle(impl.freq_response(n, d, w)))
    phi, gain = impl.lead_max_phase(alpha)
    assert ph.max() == pytest.approx(phi, abs=1e-6)
    assert w[np.argmax(ph)] == pytest.approx(1.0 / (T * np.sqrt(alpha)), rel=1e-3)
    assert abs(impl.freq_response(n, d, np.array([1.0 / (T * np.sqrt(alpha))]))[0]) == pytest.approx(gain, rel=1e-12)
    assert impl.lead_alpha(phi) == pytest.approx(alpha, rel=1e-12)       # the two formulas are inverses


def test_lead_limits():
    assert impl.lead_max_phase(1.0)[0] == pytest.approx(0.0, abs=1e-12)   # zero on the pole: no lead
    assert impl.lead_max_phase(1e-12)[0] == pytest.approx(90.0, abs=1e-3)  # pole at infinity: a pure zero
    n, d = impl.lead_compensator(3.0, 0.25, 0.1)
    assert abs(impl.freq_response(n, d, np.array([1e-6]))[0]) == pytest.approx(3.0, rel=1e-9)
    assert abs(impl.freq_response(n, d, np.array([1e6]))[0]) == pytest.approx(3.0 / 0.25, rel=1e-4)


def test_lag_has_dc_gain_beta_and_little_phase_a_decade_above_its_zero():
    beta, z = 6.0, 3.0
    n, d = impl.lag_compensator(beta, z)
    assert abs(impl.freq_response(n, d, np.array([1e-7]))[0]) == pytest.approx(beta, rel=1e-6)
    assert abs(impl.freq_response(n, d, np.array([1e5]))[0]) == pytest.approx(1.0, rel=1e-6)
    ph = np.degrees(np.angle(impl.freq_response(n, d, np.array([10 * z]))[0]))
    # exact: atan(10) - atan(10 beta) for the lag (s + z)/(s + z/beta) ... at w = 10 z: atan(10) - atan(10 beta)
    assert ph == pytest.approx(np.degrees(np.arctan(10.0) - np.arctan(10.0 * beta)), abs=1e-9)
    assert -6.0 < ph < 0.0


def test_velocity_constant_of_a_type_one_loop():
    assert impl.velocity_constant(NG, DG) == pytest.approx(20.0)
    with pytest.raises(ValueError):
        impl.velocity_constant(np.array([1.0]), np.array([1.0, 1.0]))   # type 0


# ---------------------------------------------------------------- margins

def test_margins_of_a_textbook_loop():
    # K/(s(s+1)(s+2)) with K = 1: phase crossover at sqrt(2) rad/s, gain margin 6 (Routh)
    den = np.polymul([1.0, 0.0], np.polymul([1.0, 1.0], [1.0, 2.0]))
    m = impl.margins(np.array([1.0]), den)
    assert m["w180"] == pytest.approx(np.sqrt(2.0), rel=1e-8)
    assert m["gm"] == pytest.approx(6.0, rel=1e-8)
    # at the gain crossover |L| = 1 and the phase margin is 180 + angle
    L = impl.freq_response(np.array([1.0]), den, np.array([m["wc"]]))[0]
    assert abs(L) == pytest.approx(1.0, rel=1e-9)
    assert m["pm"] == pytest.approx(180.0 + np.degrees(np.angle(L)), abs=1e-7)
    assert m["delay_margin"] == pytest.approx(np.radians(m["pm"]) / m["wc"])


def test_margins_agree_with_scipy_bode_and_second_order_has_infinite_gm():
    m = impl.margins(10.0 * NG, DG)
    wc, pm = bode_margins(10.0 * NG, DG)
    assert m["wc"] == pytest.approx(wc, rel=1e-5)
    assert m["pm"] == pytest.approx(pm, abs=1e-3)
    assert np.isinf(m["gm"])                                            # phase -> -180 only as w -> inf


# ---------------------------------------------------------------- the independent check: compensated margins

@pytest.mark.parametrize("wc,pm", [(30.0, 55.0), (20.0, 45.0), (50.0, 60.0)])
def test_lead_design_hits_the_target_margins(wc, pm):
    Kc, alpha, T = impl.design_lead(NG, DG, wc, pm)
    nl, dl = impl.lead_compensator(Kc, alpha, T)
    nL, dL = impl.series(nl, dl, NG, DG)
    m = impl.margins(nL, dL)
    assert m["wc"] == pytest.approx(wc, rel=1e-7)
    assert m["pm"] == pytest.approx(pm, abs=1e-5)
    # measured independently on scipy.signal.bode
    wcb, pmb = bode_margins(nL, dL)
    assert wcb == pytest.approx(wc, rel=1e-5)
    assert pmb == pytest.approx(pm, abs=1e-3)
    # the lead's maximum sits at the crossover
    assert 1.0 / (T * np.sqrt(alpha)) == pytest.approx(wc, rel=1e-12)


@pytest.mark.parametrize("ratio", [5.0, 10.0, 20.0])
def test_lead_lag_design_hits_margins_and_velocity_constant(ratio):
    d = impl.design_lead_lag(NG, DG, 30.0, 55.0, 200.0, lag_ratio=ratio)
    nL, dL = d["loop"]
    m = impl.margins(nL, dL)
    assert m["wc"] == pytest.approx(30.0, rel=1e-7)
    assert m["pm"] == pytest.approx(55.0, abs=1e-5)
    assert impl.velocity_constant(nL, dL) == pytest.approx(200.0, rel=1e-8)
    assert d["beta"] > 1.0
    assert d["z"] == pytest.approx(30.0 / ratio)
    wcb, pmb = bode_margins(nL, dL)
    assert pmb == pytest.approx(55.0, abs=1e-3)


def test_lead_lag_with_lax_kv_needs_no_lag_and_impossible_lead_raises():
    d = impl.design_lead_lag(NG, DG, 30.0, 55.0, 10.0)
    assert d["beta"] == pytest.approx(1.0)
    with pytest.raises(ValueError):
        impl.design_lead(NG, DG, 300.0, 80.0)                           # ~ 160 deg of lead needed


# ---------------------------------------------------------------- closed-loop responses

def test_closed_loop_step_matches_scipy_and_overshoot_falls_with_phase_margin():
    t = np.linspace(0.0, 1.0, 2001)
    overs = []
    for pm in (40.0, 50.0, 60.0):
        Kc, alpha, T = impl.design_lead(NG, DG, 30.0, pm)
        n, d = impl.unity_feedback(*impl.series(*impl.lead_compensator(Kc, alpha, T), NG, DG))
        y = impl.tf_response(n, d, t, np.ones_like(t))
        _, y_ref = signal.step(signal.lti(n, d), T=t)
        npt.assert_allclose(y, y_ref, atol=1e-9)
        assert y[-1] == pytest.approx(1.0, abs=1e-4)                    # type 1: no step error (transient e^-10 left at 1 s)
        overs.append(impl.step_metrics(t, y)["overshoot_pct"])
    assert overs[0] > overs[1] > overs[2]


def test_ramp_error_is_one_over_kv():
    d = impl.design_lead_lag(NG, DG, 30.0, 55.0, 200.0)
    nL, dL = d["loop"]
    t = np.linspace(0.0, 15.0, 15001)
    _, e, _ = signal.lsim((dL, np.polyadd(dL, nL)), t, t)              # E/R = 1/(1 + L), unit ramp
    assert e[-1] == pytest.approx(1.0 / 200.0, rel=1e-4)
    # the lead alone has the same crossover but five times the ramp error
    Kc, alpha, T = impl.design_lead(NG, DG, 30.0, 55.0)
    nL2, dL2 = impl.series(*impl.lead_compensator(Kc, alpha, T), NG, DG)
    _, e2, _ = signal.lsim((dL2, np.polyadd(dL2, nL2)), t, t)
    assert e2[-1] == pytest.approx(1.0 / impl.velocity_constant(nL2, dL2), rel=1e-6)
    assert e2[-1] > 4.0 * e[-1]


def test_gain_alone_cannot_meet_both_specs():
    nL, dL = impl.series(np.array([10.0]), np.array([1.0]), NG, DG)   # Kv = 200 by gain alone
    m = impl.margins(nL, dL)
    assert impl.velocity_constant(nL, dL) == pytest.approx(200.0)
    assert m["pm"] < 25.0 and m["wc"] > 50.0


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "lead_lag.png").stat().st_size > 1000
