"""Public tests for P04.14 Atomic & Radiation Physics. They assert the Bohr levels and the Rydberg
formula, photon bookkeeping, binding energies from masses against the liquid-drop model, the Pu-238
Q-value, the half-life law and its consequences (activity, decay heat, RTG decline), the Bateman chain
against solve_ivp and its secular-equilibrium limit, dose-unit relations, the Poisson upset probability,
and the exponential shielding MODEL (zero thickness, monotonicity, fit round-trip).

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-04/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)
YEAR = impl.YEAR_S
M = impl.ATOMIC_MASS_U


# ---------------------------------------------------------------- atoms and photons

def test_bohr_levels_and_reduced_mass():
    assert impl.energy_level(1, reduced_mass=False) == pytest.approx(-13.6057, abs=1e-4)
    assert impl.energy_level(1) == pytest.approx(-13.5983, abs=1e-4)
    n = np.arange(1, 8)
    E = impl.energy_level(n)
    npt.assert_allclose(E * n ** 2, E[0])                                    # 1/n^2 law
    assert np.all(np.diff(E) > 0) and np.all(E < 0)
    assert impl.energy_level(1, Z=2, reduced_mass=False) == pytest.approx(4 * impl.energy_level(1, reduced_mass=False))
    with pytest.raises(ValueError):
        impl.energy_level(0)


def test_hydrogen_lines_match_rydberg_formula_and_measured_values():
    # Rydberg: 1/lambda = R_H (1/nl^2 - 1/nu^2), R_H = R_inf / (1 + me/mp); vacuum wavelengths
    R_H = 10973731.568 / (1 + impl.M_E_U / impl.M_P_U)
    for nu, nl in ((2, 1), (3, 2), (4, 2), (5, 3)):
        lam = impl.transition_wavelength(nu, nl)
        assert 1 / lam == pytest.approx(R_H * (1 / nl ** 2 - 1 / nu ** 2), rel=1e-8)
    assert impl.transition_wavelength(2, 1) == pytest.approx(121.57e-9, rel=1e-4)    # Lyman alpha
    assert impl.transition_wavelength(3, 2) == pytest.approx(656.47e-9, rel=1e-4)    # H-alpha (vacuum)
    with pytest.raises(ValueError):
        impl.transition_wavelength(2, 3)


def test_photon_energy_and_lidar_bookkeeping():
    assert impl.photon_energy_ev(1239.84198e-9) == pytest.approx(1.0, rel=1e-8)  # hc = 1239.84 eV nm
    lam = np.array([400e-9, 800e-9])
    E = impl.photon_energy_ev(lam)
    assert E[0] == pytest.approx(2 * E[1])
    n = impl.photons_per_pulse(10e-3, 1064e-9)
    assert n == pytest.approx(5.356e16, rel=1e-3)
    assert n * impl.photon_energy_ev(1064e-9) * impl.E_CHARGE == pytest.approx(10e-3)
    assert impl.lidar_range(2 * 1000.0 / impl.C) == pytest.approx(1000.0)


# ---------------------------------------------------------------- nuclei

def test_binding_energy_from_masses():
    assert impl.binding_energy_mev(1, 1, impl.M_H_U) == pytest.approx(0.0, abs=1e-9)  # a lone proton is not bound
    assert impl.binding_energy_mev(2, 4, M["He-4"]) == pytest.approx(28.30, abs=0.01)
    bpa = {k: impl.binding_energy_mev(Z, A, M[k]) / A for k, Z, A in
           (("He-4", 2, 4), ("C-12", 6, 12), ("O-16", 8, 16), ("Fe-56", 26, 56), ("Ni-62", 28, 62), ("U-238", 92, 238))}
    assert bpa["Fe-56"] == pytest.approx(8.790, abs=0.002)
    assert max(bpa, key=bpa.get) in ("Fe-56", "Ni-62")                       # the peak of the curve
    assert bpa["U-238"] < bpa["Fe-56"] and bpa["He-4"] < bpa["C-12"]


def test_semf_is_a_good_model_for_heavy_nuclei_and_a_poor_one_for_light():
    for k, Z, A in (("Fe-56", 26, 56), ("U-238", 92, 238), ("Pu-238", 94, 238)):
        assert impl.semf_binding_energy_mev(Z, A) == pytest.approx(impl.binding_energy_mev(Z, A, M[k]), rel=0.01)
    assert abs(impl.semf_binding_energy_mev(2, 4) / impl.binding_energy_mev(2, 4, M["He-4"]) - 1) > 0.1


def test_pu238_alpha_decay_q_value():
    q = impl.q_value_mev(M["Pu-238"], M["U-234"], M["He-4"])
    assert q == pytest.approx(5.593, abs=0.003)
    # Q equals the gain in binding energy, an independent route through the same masses
    dB = (impl.binding_energy_mev(92, 234, M["U-234"]) + impl.binding_energy_mev(2, 4, M["He-4"])
          - impl.binding_energy_mev(94, 238, M["Pu-238"]))
    assert q == pytest.approx(dB, rel=1e-9)
    assert impl.q_value_mev(M["Fe-56"], M["He-4"], 51.940505) < 0              # Fe-56 cannot alpha-decay (Cr-52)


# ---------------------------------------------------------------- decay law

def test_half_life_law():
    T = 87.7
    assert impl.remaining_fraction(0.0, T) == 1.0
    assert impl.remaining_fraction(T, T) == pytest.approx(0.5)
    assert impl.remaining_fraction(3 * T, T) == pytest.approx(0.125)
    t = np.array([5.0, 12.0])
    assert impl.remaining_fraction(t.sum(), T) == pytest.approx(np.prod(impl.remaining_fraction(t, T)))  # memoryless
    assert impl.decay_constant(T) * T == pytest.approx(np.log(2))
    with pytest.raises(ValueError):
        impl.decay_constant(0.0)


def test_activity_specific_activity_and_decay_heat():
    T = impl.PU238_HALF_LIFE_Y * YEAR
    assert impl.activity_bq(1e20, T) == pytest.approx(np.log(2) / T * 1e20)
    sa = impl.specific_activity_bq_per_g(T, impl.PU238_MOLAR_MASS)
    assert sa == pytest.approx(6.336e11, rel=2e-3)
    assert sa / 3.7e10 == pytest.approx(17.1, abs=0.1)                       # Ci/g
    q = impl.q_value_mev(M["Pu-238"], M["U-234"], M["He-4"])
    w = impl.decay_heat_w_per_g(T, impl.PU238_MOLAR_MASS, q)
    assert w == pytest.approx(0.568, abs=0.002)
    # longer half-life, lower specific power: U-238 (4.468e9 y) is ~5e7 times weaker per gram
    assert impl.specific_activity_bq_per_g(4.468e9 * YEAR, 238.05) < sa / 1e7


def test_bateman_chain_matches_solve_ivp_and_reaches_secular_equilibrium():
    lam1, lam2 = impl.decay_constant(2.0), impl.decay_constant(0.3)
    t = np.linspace(0.0, 6.0, 61)
    N1, N2 = impl.bateman_two(t, 1e6, lam1, lam2)
    n1, n2 = impl.chain_numeric(t, 1e6, lam1, lam2)
    npt.assert_allclose(N1, n1, rtol=1e-7, atol=1e-3)
    npt.assert_allclose(N2, n2, rtol=1e-7, atol=1e-3)
    assert N2[0] == 0.0
    # long-lived parent: the activity ratio tends to l2/(l2 - l1), which -> 1 (secular equilibrium)
    l1, l2 = impl.decay_constant(1e4), impl.decay_constant(1.0)
    M1, M2 = impl.bateman_two(np.array([40.0]), 1.0, l1, l2)
    assert l2 * M2[0] / (l1 * M1[0]) == pytest.approx(l2 / (l2 - l1), rel=1e-9)
    assert l2 * M2[0] == pytest.approx(l1 * M1[0], rel=2e-4)
    # conservation: nothing is lost before the stable end of the chain is counted
    assert np.all(N1 + N2 <= 1e6 * (1 + 1e-12))
    with pytest.raises(ValueError):
        impl.bateman_two(t, 1.0, 0.5, 0.5)


# ---------------------------------------------------------------- dose

def test_dose_units_and_weighting():
    assert impl.absorbed_dose_gy(1.0, 1.0) == 1.0
    assert impl.absorbed_dose_gy(0.5, 70.0) == pytest.approx(0.5 / 70.0)
    assert impl.gy_to_krad(10.0) == pytest.approx(1.0)
    assert impl.equivalent_dose_sv({"photon": 1e-3}) == pytest.approx(1e-3)
    assert impl.equivalent_dose_sv({"alpha": 1e-3}) == pytest.approx(20e-3)
    assert impl.equivalent_dose_sv({"photon": 1e-3, "proton": 1e-3}) == pytest.approx(3e-3)
    with pytest.raises(ValueError):
        impl.absorbed_dose_gy(1.0, 0.0)


def test_dose_from_fluence_is_linear_and_dimensionally_right():
    # 1 MeV/g = 1.602e-13 J / 1e-3 kg = 1.602e-10 Gy
    assert impl.dose_from_fluence_gy(1.0, 1.0) == pytest.approx(1.602176634e-10)
    assert impl.dose_from_fluence_gy(2e9, 1.6) == pytest.approx(2 * impl.dose_from_fluence_gy(1e9, 1.6))
    assert impl.dose_from_fluence_gy(0.0, 1.6) == 0.0


def test_tid_and_seu_bookkeeping():
    assert impl.total_ionizing_dose(3.0, 10.0) == pytest.approx(60.0)
    assert impl.total_ionizing_dose(3.0, 10.0, margin=1.0) == pytest.approx(30.0)
    r = impl.seu_rate_per_day(1e-14, 1.0, 8 * 2 ** 30)
    assert r == pytest.approx(7.42, rel=1e-3)
    assert impl.seu_rate_per_day(2e-14, 1.0, 8 * 2 ** 30) == pytest.approx(2 * r)
    assert impl.prob_at_least_one(0.0, 10.0) == 0.0
    assert impl.prob_at_least_one(1.0, 1.0) == pytest.approx(1 - np.exp(-1))
    assert impl.prob_at_least_one(1e-12, 1.0) == pytest.approx(1e-12, rel=1e-6)  # small-rate limit, no cancellation
    rng = np.random.default_rng(1414)
    # 20 000 Poisson draws at mean 0.31: the fraction with >= 1 event has a standard error of 0.003,
    # so a 0.015 tolerance is five standard errors
    assert np.mean(rng.poisson(0.31, 20_000) >= 1) == pytest.approx(impl.prob_at_least_one(0.31, 1.0), abs=0.015)


# ---------------------------------------------------------------- RTG

def test_rtg_power_follows_the_half_life():
    P = impl.rtg_electrical_power(np.array([0.0, 87.7, 175.4]), 2000.0, 0.06)
    npt.assert_allclose(P, [120.0, 60.0, 30.0])
    Pd = impl.rtg_electrical_power(14.0, 2000.0, 0.06, degradation_per_year=0.015)
    assert Pd < impl.rtg_electrical_power(14.0, 2000.0, 0.06)
    assert impl.rtg_electrical_power(14.0, 2000.0, 0.06) / 120.0 == pytest.approx(0.8953, abs=1e-4)


def test_rtg_fuel_mass_round_trip():
    q = impl.q_value_mev(M["Pu-238"], M["U-234"], M["He-4"])
    m = impl.rtg_fuel_mass_g(110.0, 14.0, 0.06, q)
    assert m == pytest.approx(3610.0, rel=2e-3)
    w = impl.decay_heat_w_per_g(impl.PU238_HALF_LIFE_Y * YEAR, impl.PU238_MOLAR_MASS, q)
    assert impl.rtg_electrical_power(14.0, m * w, 0.06) == pytest.approx(110.0)
    assert impl.rtg_fuel_mass_g(110.0, 14.0, 0.06, q, 0.015) > m


# ---------------------------------------------------------------- shielding model

def test_shielding_zero_thickness_and_monotonicity():
    assert impl.shielded_dose(44.7, 0.0, 0.3355) == 44.7                    # zero thickness: no attenuation
    x = np.linspace(0, 3, 31)
    D = impl.shielded_dose(10.0, x, 0.5)
    assert np.all(np.diff(D) < 0)
    npt.assert_allclose(D[10] * D[10], D[0] * D[20])                          # equal steps, equal factors
    assert impl.areal_density_g_cm2(1.0) == pytest.approx(0.27)


def test_attenuation_fit_round_trip_and_thickness_for_target():
    d0, L = impl.fit_attenuation(0.27, 20.0, 0.81, 4.0)
    assert impl.shielded_dose(d0, 0.27, L) == pytest.approx(20.0)
    assert impl.shielded_dose(d0, 0.81, L) == pytest.approx(4.0)
    assert L == pytest.approx(0.3355, abs=1e-4)
    x = impl.areal_density_for_dose(1.5, d0, L)
    assert impl.shielded_dose(d0, x, L) == pytest.approx(1.5)
    assert x / 0.27 == pytest.approx(4.22, abs=0.01)                          # mm of aluminium
    assert impl.areal_density_for_dose(2 * d0, d0, L) == 0.0
    # narrow-beam photons: L = 1/(mu/rho); 1 MeV in Al, mu/rho = 0.0615 cm^2/g -> half-value layer 11.3 g/cm^2
    hvl = impl.areal_density_for_dose(0.5, 1.0, 1 / 0.0615)
    assert hvl == pytest.approx(np.log(2) / 0.0615)
    with pytest.raises(ValueError):
        impl.fit_attenuation(0.27, 20.0, 0.27, 4.0)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "radiation_dose.png").exists()
