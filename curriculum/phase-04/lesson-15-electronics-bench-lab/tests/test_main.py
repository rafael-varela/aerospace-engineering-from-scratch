"""Public tests for P04.15 Electronics Bench Lab. They assert the divider and LED design rules, standard
values, meter loading against the Thevenin formula, the ADC LSB and the LSB/sqrt(12) quantization noise
(checked statistically with a seeded ramp), the calibration fit, the beta and Steinhart-Hart models and
their fits (round trips and recovery of a known truth), the ratiometric divider, the uncertainty
propagation against an analytic case, and the SYNTHETIC dataset's labelling and reproducibility.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-04/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)
K = 273.15


# ---------------------------------------------------------------- design

def test_divider_limits_and_symmetry():
    assert impl.divider_output(5.0, 1e3, 1e3) == pytest.approx(2.5)
    assert impl.divider_output(5.0, 0.0, 1e3) == pytest.approx(5.0)          # no top resistor: full supply
    assert impl.divider_output(5.0, 1e3, 0.0) == 0.0
    assert impl.divider_output(3.3, 4.7e3, 10e3) + impl.divider_output(3.3, 10e3, 4.7e3) == pytest.approx(3.3)
    assert impl.divider_output(5.0, 1e3, 1e3, r_load=np.inf) == impl.divider_output(5.0, 1e3, 1e3)
    # dimensionless in the resistances: scaling both by 1000 changes nothing
    assert impl.divider_output(5.0, 8.2, 12.0) == pytest.approx(impl.divider_output(5.0, 8.2e3, 12e3))


def test_divider_design_round_trip_and_standard_values():
    rt, rb = impl.divider_design(5.0, 3.0, 20e3)
    assert (rt, rb) == pytest.approx((8e3, 12e3))
    assert impl.divider_output(5.0, rt, rb) == pytest.approx(3.0)
    with pytest.raises(ValueError):
        impl.divider_design(3.3, 5.0, 10e3)
    assert impl.nearest_e_series(8000.0) == pytest.approx(8200.0)
    assert impl.nearest_e_series(300.0, direction="up") == pytest.approx(330.0)
    assert impl.nearest_e_series(300.0, direction="down") == pytest.approx(270.0)
    assert impl.nearest_e_series(10e3) == pytest.approx(10e3)
    assert impl.nearest_e_series(0.47) == pytest.approx(0.47)
    with pytest.raises(ValueError):
        impl.nearest_e_series(-1.0)


def test_led_resistor_current_and_power():
    r = impl.led_resistor(5.0, 2.0, 10e-3)
    assert r == pytest.approx(300.0)
    assert impl.led_current(5.0, 2.0, r) == pytest.approx(10e-3)
    ru = impl.nearest_e_series(r, direction="up")
    assert impl.led_current(5.0, 2.0, ru) < 10e-3                            # rounding up never overdrives
    assert impl.resistor_power(impl.led_current(5.0, 2.0, 330.0), 330.0) == pytest.approx(27.3e-3, rel=2e-3)
    assert impl.led_current(1.5, 2.0, 100.0) == 0.0                          # below forward voltage: off
    with pytest.raises(ValueError):
        impl.led_resistor(1.8, 2.0, 10e-3)


# ---------------------------------------------------------------- instruments

def test_meter_loading_matches_thevenin_formula():
    for rt, rb in ((1e3, 1e3), (100e3, 47e3), (1e6, 1e6), (2.2e6, 330e3)):
        rth = impl.thevenin_resistance(rt, rb)
        assert impl.loading_error(rt, rb, 10e6) == pytest.approx(-rth / (rth + 10e6), rel=1e-12)
    assert impl.loading_error(1e6, 1e6, 10e6) == pytest.approx(-1 / 21)       # -4.76 %
    assert impl.loading_error(1e3, 1e3, np.inf) == 0.0                        # ideal voltmeter
    errs = [impl.loading_error(r, r, 10e6) for r in (1e3, 1e4, 1e5, 1e6)]
    assert np.all(np.diff(errs) < 0)                                          # worse for larger resistors


def test_ammeter_shunt_reduces_the_current():
    i0 = impl.led_current(5.0, 2.0, 330.0)
    i1 = impl.led_current(5.0, 2.0, 330.0, r_extra=10.0)
    assert i1 / i0 - 1 == pytest.approx(-10.0 / 340.0, rel=1e-12)
    assert impl.led_current(5.0, 2.0, 330.0, r_extra=0.0) == i0


# ---------------------------------------------------------------- the ADC

def test_lsb_and_quantization_noise():
    assert impl.lsb(3.3, 12) == pytest.approx(3.3 / 4096)
    assert impl.lsb(3.3, 10) == pytest.approx(4 * impl.lsb(3.3, 12))          # 2 bits = factor 4
    assert impl.quantization_noise_rms(3.3, 12) == pytest.approx(3.3 / 4096 / np.sqrt(12))
    # statistical check: a dense ramp quantized by the ideal quantizer has rms error LSB/sqrt(12);
    # 200 001 uniformly spaced inputs make the sample rms exact to ~1e-4 relative, 1 % tolerance is ample
    v = np.linspace(0.5, 2.5, 200_001)
    err = impl.counts_to_volts(impl.volts_to_counts(v, 3.3, 12), 3.3, 12) - v
    assert np.sqrt(np.mean(err ** 2)) == pytest.approx(impl.quantization_noise_rms(3.3, 12), rel=0.01)
    assert np.max(np.abs(err)) <= impl.lsb(3.3, 12) / 2 + 1e-12


def test_quantizer_clips_and_inverts():
    c = impl.volts_to_counts(np.array([-0.1, 0.0, 1.65, 3.3, 5.0]), 3.3, 12)
    npt.assert_array_equal(c[[0, 1]], [0, 0])
    assert c[2] == 2048
    npt.assert_array_equal(c[[3, 4]], [4095, 4095])
    assert impl.u16_to_volts(0) == 0.0
    assert impl.u16_to_volts(65535, 3.3) == pytest.approx(3.3)
    assert impl.u16_to_volts(32768, 5.0) == pytest.approx(2.5, rel=1e-4)


def test_adc_calibration_fit_recovers_a_known_line():
    raw = np.array([5000.0, 15000.0, 30000.0, 45000.0, 60000.0])
    volts = 5.02e-5 * raw - 0.008
    g, o, s = impl.adc_calibration_fit(raw, volts)
    assert g == pytest.approx(5.02e-5, rel=1e-10) and o == pytest.approx(-0.008, abs=1e-12)
    assert s == pytest.approx(0.0, abs=1e-12)


# ---------------------------------------------------------------- the thermistor

def test_beta_model_round_trip_and_reference_point():
    T = np.linspace(260.0, 340.0, 9)
    R = impl.thermistor_beta_r(T, 10e3, 3950.0)
    npt.assert_allclose(impl.thermistor_beta_t(R, 10e3, 3950.0), T, rtol=1e-12)
    assert impl.thermistor_beta_r(298.15, 10e3, 3950.0) == pytest.approx(10e3)
    assert np.all(np.diff(R) < 0)                                             # NTC: R falls with T
    r0, beta = impl.fit_beta(R, T)
    assert (r0, beta) == pytest.approx((10e3, 3950.0), rel=1e-10)


def test_steinhart_hart_inverse_and_fit_recover_the_truth():
    a, b, c = impl.SYNTH_SH
    T = np.array([270.0, 285.0, 300.0, 315.0, 330.0])
    R = impl.steinhart_hart_r(T, a, b, c)
    npt.assert_allclose(impl.steinhart_hart_t(R, a, b, c), T, rtol=1e-12)
    fit = impl.fit_steinhart_hart(R, T)
    npt.assert_allclose(fit, (a, b, c), rtol=1e-6)
    assert impl.steinhart_hart_r(298.15, a, b, c) == pytest.approx(10e3, rel=0.01)  # a 10 kOhm part
    with pytest.raises(ValueError):
        impl.fit_steinhart_hart(R[:2], T[:2])


def test_beta_is_steinhart_hart_with_c_zero():
    r0, beta, T0 = 10e3, 3950.0, 298.15
    a = 1 / T0 - np.log(r0) / beta
    R = np.array([30e3, 10e3, 4e3])
    npt.assert_allclose(impl.steinhart_hart_t(R, a, 1 / beta, 0.0), impl.thermistor_beta_t(R, r0, beta), rtol=1e-12)


def test_ratiometric_divider_inverse_and_fixed_resistor_choice():
    rf = 10e3
    Rth = np.array([3e3, 10e3, 33e3])
    ratio = rf / (Rth + rf)
    npt.assert_allclose(impl.thermistor_r_from_ratio(ratio, rf), Rth, rtol=1e-12)
    assert impl.thermistor_r_from_ratio(0.5, rf) == pytest.approx(rf)
    with pytest.raises(ValueError):
        impl.thermistor_r_from_ratio(1.0, rf)
    # the geometric mean maximises the swing between two temperatures
    Rlo, Rhi = impl.thermistor_beta_r(273.15, 10e3, 3950.0), impl.thermistor_beta_r(323.15, 10e3, 3950.0)
    rbest = impl.best_fixed_resistor(273.15, 323.15, 10e3, 3950.0)
    swing = lambda r: r / (Rhi + r) - r / (Rlo + r)
    assert swing(rbest) > swing(0.8 * rbest) and swing(rbest) > swing(1.25 * rbest)


def test_temperature_from_u16_round_trip():
    sh = impl.SYNTH_SH
    for Tc in (5.0, 25.0, 45.0):
        R = impl.steinhart_hart_r(Tc + K, *sh)
        raw = 65535 * 10e3 / (R + 10e3)
        assert impl.temperature_from_u16(raw, 10e3, sh) == pytest.approx(Tc, abs=1e-9)


# ---------------------------------------------------------------- uncertainty

def test_dmm_uncertainty_and_propagation_against_analytic_cases():
    assert impl.dmm_standard_uncertainty(3.000, 0.5, 2, 0.001) == pytest.approx((0.015 + 0.002) / np.sqrt(3))
    assert impl.dmm_standard_uncertainty(0.0, 0.5, 0, 0.001) == 0.0
    # product: relative uncertainties add in quadrature
    v, u = impl.propagate(lambda p: p[0] * p[1], np.array([2.0, 5.0]), np.array([0.02, 0.1]))
    assert v == pytest.approx(10.0)
    assert u == pytest.approx(10.0 * np.hypot(0.01, 0.02), rel=1e-6)
    # zero input uncertainty: zero output uncertainty
    assert impl.propagate(lambda p: np.exp(p[0]), np.array([1.0]), np.array([0.0]))[1] == 0.0
    # power law: LED current through R has u_I/I = u_R/R
    v, u = impl.propagate(lambda p: 3.0 / p[0], np.array([330.0]), np.array([3.3]))
    assert u / v == pytest.approx(0.01, rel=1e-6)


# ---------------------------------------------------------------- SYNTHETIC data

def test_synthetic_dataset_is_reproducible_and_physical():
    d1, d2 = impl.make_synthetic_bench_data(7), impl.make_synthetic_bench_data(7)
    for k in d1:
        npt.assert_array_equal(d1[k], d2[k])
    assert not np.array_equal(impl.make_synthetic_bench_data(8)["u16_mean"], d1["u16_mean"])
    assert np.all(np.diff(d1["u16_mean"]) > 0)                               # hotter: lower R, higher voltage
    assert np.all((d1["u16_mean"] > 0) & (d1["u16_mean"] < 65535))


def test_analysis_of_synthetic_data_recovers_the_thermistor():
    # the reference temperatures carry sigma = 0.25 C noise and R_fixed is 0.5 % off, so the fits
    # should land within about 1 C everywhere and the beta value within a few percent of the truth's 3973 K
    d = impl.make_synthetic_bench_data(2026)
    a = impl.analyse_thermistor(d)
    assert np.max(np.abs(a["res_sh_C"])) < 1.0
    assert np.max(np.abs(a["res_beta_C"])) < 1.0
    r0, beta = a["beta"]
    assert r0 == pytest.approx(10e3, rel=0.03) and beta == pytest.approx(3973.0, rel=0.03)
    g, o, s = impl.adc_calibration_fit(d["pot_u16"], d["pot_dmm_V"])
    assert g == pytest.approx(3.3 / 65535 / 1.008, rel=0.005) and abs(o) < 0.02 and s < 0.005


def test_synthetic_csv_is_labelled_and_reads_back(tmp_path):
    p1, p2 = impl.write_synthetic_csv(tmp_path, seed=3)
    for p in (p1, p2):
        assert "synthetic" in p.name
        first = p.read_text().splitlines()[0]
        assert first.startswith("#") and "SYNTHETIC" in first
    back = impl.read_bench_csv(p1)
    d = impl.make_synthetic_bench_data(3)
    npt.assert_allclose(back["u16_mean"], d["u16_mean"], atol=0.006)
    assert set(back["note"]) == {"synthetic"}


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "bench_lab.png").exists()
