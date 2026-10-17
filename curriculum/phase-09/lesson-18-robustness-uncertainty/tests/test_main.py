"""Tests for P09.18 Robustness & Uncertainty.

Run against the starter (default) or the reference solution (AEFS_IMPL=solution).
The independent checks: the H-infinity norm by frequency sweep against Hamiltonian bisection; the disk-margin
bounds implied by Ms against the gain and phase margins measured on the loop.
"""
from __future__ import annotations

import numpy as np
import numpy.testing as npt
import pytest
from conftest import load_impl
from scipy.signal import freqs

impl = load_impl(__file__)

RNG = np.random.default_rng(20250918)
W = np.logspace(-1, 4, 40000)


def random_stable_tf(n: int) -> tuple[np.ndarray, np.ndarray]:
    poles = []
    while len(poles) < n:
        if n - len(poles) >= 2 and RNG.random() < 0.6:
            wn, z = RNG.uniform(0.5, 20.0), RNG.uniform(0.05, 0.9)
            poles += [complex(-z * wn, wn * np.sqrt(1 - z * z)), complex(-z * wn, -wn * np.sqrt(1 - z * z))]
        else:
            poles.append(-RNG.uniform(0.2, 30.0))
    den = np.real(np.poly(poles))
    num = RNG.standard_normal(RNG.integers(1, n + 1)) * den[-1]
    return num, den


# ---------------------------------------------------------------- realisation and closed-loop maps

@pytest.mark.parametrize("n", [1, 2, 4])
def test_tf_to_ss_reproduces_the_frequency_response(n):
    num, den = random_stable_tf(n)
    A, B, C, D = impl.tf_to_ss(num, den)
    w = np.logspace(-2, 3, 50)
    G = np.array([(C @ np.linalg.solve(1j * wi * np.eye(n) - A, B) + D)[0, 0] for wi in w])
    _, ref = freqs(num, den, worN=w)
    npt.assert_allclose(G, ref, rtol=1e-9, atol=1e-12)
    npt.assert_allclose(np.sort_complex(np.linalg.eigvals(A)), np.sort_complex(np.roots(den)), rtol=1e-8, atol=1e-10)
    with pytest.raises(ValueError):
        impl.tf_to_ss([1.0, 0.0, 0.0], [1.0, 1.0])


def test_sensitivity_identities_and_limits():
    nG, dG = impl.servo_plant()
    nC, dC = impl.lead(0.8, 20.0, 80.0)
    t = impl.closed_loop_tfs(nC, dC, nG, dG)
    L = impl.freqresp(*t["L"], W)
    S, T = impl.freqresp(*t["S"], W), impl.freqresp(*t["T"], W)
    npt.assert_allclose(S + T, 1.0, atol=1e-12)                                # S + T = 1
    npt.assert_allclose(S, 1.0 / (1.0 + L), rtol=1e-10)
    npt.assert_allclose(impl.freqresp(*t["KS"], W), impl.freqresp(nC, dC, W) * S, rtol=1e-10)
    npt.assert_allclose(impl.freqresp(*t["GS"], W), impl.freqresp(nG, dG, W) * S, rtol=1e-10)
    # type-1 loop: S(0) = 0, T(0) = 1; at high frequency S -> 1, T -> 0, KS -> C(inf) = Kc p / z
    assert abs(impl.freqresp(*t["S"], [1e-6])[0]) < 1e-6
    npt.assert_allclose(impl.freqresp(*t["T"], [1e-6])[0], 1.0, atol=1e-6)
    npt.assert_allclose(impl.freqresp(*t["KS"], [1e7])[0], 0.8 * 80.0 / 20.0, rtol=1e-3)
    assert impl.tf_hinf_norm(*t["KS"]) >= 3.2 * (1 - 1e-6)


def test_closed_loop_stability_test():
    nG, dG = impl.servo_plant()
    assert impl.is_closed_loop_stable(*impl.lead(0.8, 20.0, 80.0), nG, dG)
    nU, dU = impl.unstable_airframe()
    assert impl.is_closed_loop_stable(*impl.lead(0.6, 1.0, 40.0), nU, dU)
    assert not impl.is_closed_loop_stable(*impl.lead(0.05, 1.0, 40.0), nU, dU)    # too little gain: diverges
    assert not impl.is_closed_loop_stable(*impl.lead(3.2, 25.0, 200.0), *impl.true_servo_plant(4))


# ---------------------------------------------------------------- H-infinity norm: the independent check

@pytest.mark.parametrize("n", [1, 2, 3, 5])
def test_hinf_sweep_agrees_with_hamiltonian_bisection(n):
    num, den = random_stable_tf(n)
    A, B, C, D = impl.tf_to_ss(num, den)
    sweep, wpk = impl.hinf_norm_sweep(A, B, C, D, np.logspace(-3, 4, 3000))
    bis = impl.hinf_norm_bisection(A, B, C, D)
    npt.assert_allclose(sweep, bis, rtol=1e-6)


@pytest.mark.parametrize("zeta", [0.05, 0.2, 0.5])
def test_hinf_of_a_resonance_has_the_closed_form_peak(zeta):
    wn = 3.0
    num, den = np.array([wn ** 2]), np.array([1.0, 2 * zeta * wn, wn ** 2])
    peak = 1.0 / (2 * zeta * np.sqrt(1 - zeta ** 2))
    npt.assert_allclose(impl.tf_hinf_norm(num, den, "bisection"), peak, rtol=1e-7)
    npt.assert_allclose(impl.tf_hinf_norm(num, den, "sweep"), peak, rtol=1e-7)
    # first order: ||1/(s + a)|| = 1/a, at w = 0
    npt.assert_allclose(impl.tf_hinf_norm([1.0], [1.0, 4.0]), 0.25, rtol=1e-7)


def test_hamiltonian_test_switches_at_the_norm():
    num, den = np.array([9.0]), np.array([1.0, 0.6, 9.0])
    A, B, C, D = impl.tf_to_ss(num, den)
    g = impl.hinf_norm_bisection(A, B, C, D)
    assert impl.hamiltonian_has_imaginary_eigs(A, B, C, D, 0.99 * g)
    assert not impl.hamiltonian_has_imaginary_eigs(A, B, C, D, 1.01 * g)
    with pytest.raises(ValueError):
        impl.hinf_norm_bisection(np.array([[0.5]]), np.array([[1.0]]), np.array([[1.0]]), np.array([[0.0]]))


# ---------------------------------------------------------------- margins: the disk-margin check

def test_measured_margins_on_a_known_loop():
    """L = 2/(s(s + 1)(s + 2)): phase crossover at sqrt(2) rad/s where |L| = 1/3, so GM = 3."""
    w = np.logspace(-2, 2, 200001)
    L = impl.freqresp([2.0], np.poly([0.0, -1.0, -2.0]), w)
    m = impl.margins(w, L)
    npt.assert_allclose(m["gm_high"], 3.0, rtol=1e-6)
    assert m["gm_low"] == 0.0
    npt.assert_allclose(abs(impl.freqresp([2.0], np.poly([0.0, -1.0, -2.0]), [m["wgc"]])[0]), 1.0, rtol=1e-4)
    assert 0.0 < m["pm"] < 90.0


def test_disk_margin_formula():
    gl, gh, pm = impl.disk_margin_bounds(2.0)
    npt.assert_allclose([gl, gh, pm], [2.0 / 3.0, 2.0, 2.0 * np.degrees(np.arcsin(0.25))], rtol=1e-12)
    gl, gh, pm = impl.disk_margin_bounds(1.0 + 1e-9)               # Ms -> 1: GM -> (1/2, inf), PM -> 60 deg
    assert gh > 1e8 and abs(gl - 0.5) < 1e-8 and abs(pm - 60.0) < 1e-6


@pytest.mark.parametrize("which", ["servo-0.5", "servo-1.2", "true-servo-0.8", "true-servo-2.5", "unstable-airframe"])
def test_disk_margin_bound_from_ms_holds_for_measured_margins(which):
    w = np.logspace(-4, 4, 240000)
    if which.startswith("unstable"):
        L = impl.freqresp(*impl.lead(0.6, 1.0, 40.0), w) * impl.freqresp(*impl.unstable_airframe(), w)
    else:
        Kc = float(which.split("-")[-1])
        z, p = (20.0, 80.0) if Kc < 2 else (25.0, 150.0)
        L = impl.freqresp(*impl.lead(Kc, z, p), w) * impl.freqresp(*impl.servo_plant(), w)
        if which.startswith("true"):
            L = L * (1.0 + impl.true_servo_relative_error(w))
    Ms = np.max(np.abs(1.0 / (1.0 + L)))
    gl, gh, pmb = impl.disk_margin_bounds(Ms)
    m = impl.margins(w, L)
    assert m["pm"] >= pmb - 0.05
    assert m["gm_high"] >= gh * (1 - 1e-3)
    assert m["gm_low"] <= gl * (1 + 1e-3)
    if which.startswith("unstable"):
        npt.assert_allclose(m["gm_low"], 1.0 / 6.0, rtol=1e-3)   # L(0) = -6: gain can fall to 1/6 before divergence


# ---------------------------------------------------------------- uncertainty and small gain

def test_weight_covers_the_true_actuator_error():
    w = np.logspace(-1, 4, 40000)
    assert impl.weight_covers(w, impl.true_servo_relative_error(w), *impl.uncertainty_weight()) <= 1.0
    # and the Pade model of the true plant is the same plant at low and mid frequency
    G0 = impl.freqresp(*impl.servo_plant(), w[:20000])
    Gt = impl.freqresp(*impl.true_servo_plant(4), w[:20000])
    npt.assert_allclose(Gt / G0 - 1.0, impl.true_servo_relative_error(w[:20000]), atol=1e-6)


def test_small_gain_passed_means_every_covered_perturbation_is_stable():
    nG, dG = impl.servo_plant()
    nC, dC = impl.lead(0.8, 20.0, 80.0)
    Wn, Wd = impl.uncertainty_weight()
    assert impl.robust_stability_index(nC, dC, nG, dG, Wn, Wd) < 1.0
    assert impl.is_closed_loop_stable(nC, dC, *impl.true_servo_plant(4))
    for _ in range(20):                                # random all-pass and constant Delta with |Delta| <= 0.99
        c = RNG.uniform(-0.99, 0.99)
        a = RNG.uniform(1.0, 1000.0)
        dn, dd = (np.array([-c, c * a]), np.array([1.0, a])) if RNG.random() < 0.7 else (np.array([c]), np.array([1.0]))
        nGp = np.polymul(nG, np.polyadd(np.polymul(Wd, dd), np.polymul(Wn, dn)))
        dGp = np.polymul(dG, np.polymul(Wd, dd))
        assert impl.is_closed_loop_stable(nC, dC, nGp, dGp)


def test_small_gain_failed_means_a_destabilising_perturbation_exists():
    nG, dG = impl.servo_plant()
    nC, dC = impl.lead(2.5, 25.0, 150.0)
    Wn, Wd = impl.uncertainty_weight()
    nT, dT = impl.closed_loop_tfs(nC, dC, nG, dG)["T"]
    nWT, dWT = np.polymul(Wn, nT), np.polymul(Wd, dT)
    peak, wst = impl.hinf_norm_sweep(*impl.tf_to_ss(nWT, dWT))
    assert peak > 1.0
    target = -1.0 / impl.freqresp(nWT, dWT, [wst])[0]
    dn, dd = impl.destabilising_perturbation(wst, target)
    npt.assert_allclose(impl.freqresp(dn, dd, [wst])[0], target, rtol=1e-10)
    npt.assert_allclose(np.abs(impl.freqresp(dn, dd, np.logspace(-2, 5, 50))), abs(target), rtol=1e-10)   # all-pass
    assert abs(target) < 1.0
    nGp = np.polymul(nG, np.polyadd(np.polymul(Wd, dd), np.polymul(Wn, dn)))
    dGp = np.polymul(dG, np.polymul(Wd, dd))
    r = np.roots(np.polyadd(np.polymul(dC, dGp), np.polymul(nC, nGp)))
    assert np.min(np.abs(r - 1j * wst)) < 1e-6 * wst


def test_design_table_orders_robustness_with_gain():
    rows = impl.design_table([(0.5, 20.0, 80.0), (1.2, 20.0, 80.0), (3.2, 25.0, 200.0)])
    assert rows[0]["WT"] < rows[1]["WT"] < 1.0 < rows[2]["WT"]
    assert rows[0]["true_stable"] and rows[1]["true_stable"] and not rows[2]["true_stable"]
    assert rows[0]["pm_true"] > rows[1]["pm_true"]


# ---------------------------------------------------------------- Bode sensitivity integral

def test_bode_integral_is_zero_for_a_stable_loop_and_pi_p_for_an_unstable_one():
    nG, dG = impl.servo_plant()
    t = impl.closed_loop_tfs(*impl.lead(0.8, 20.0, 80.0), nG, dG)
    assert abs(impl.bode_sensitivity_integral(*t["S"])) < 1e-3
    tu = impl.closed_loop_tfs(*impl.lead(0.6, 1.0, 40.0), *impl.unstable_airframe())
    npt.assert_allclose(impl.bode_sensitivity_integral(*tu["S"]), 2.0 * np.pi, rtol=1e-4)
    # a different stabilising controller changes |S| but not the area
    tu2 = impl.closed_loop_tfs(*impl.lead(1.0, 1.5, 60.0), *impl.unstable_airframe())
    assert impl.is_closed_loop_stable(*impl.lead(1.0, 1.5, 60.0), *impl.unstable_airframe())
    npt.assert_allclose(impl.bode_sensitivity_integral(*tu2["S"]), 2.0 * np.pi, rtol=1e-4)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "robustness.png").stat().st_size > 1000
