"""Public tests for P07.01 Systems, Properties & State. They assert thermodynamic identities and
independent data, not magic numbers: cp = dh/dT and ds deg/dT = cp/T for the NASA polynomials; h(298.15 K)
equals the tabulated enthalpy of formation; the two NASA intervals agree at 1000 K; the equipartition
limits (5/2 R for helium, 7/2 R for N2 near room temperature); dry-air cp against Lienhard Table A.6;
p v = R T against the kinetic-theory pressure of P04.07; the state postulate (any independent pair fixes
the state, a dependent pair is refused); and the size of the constant-cp error at turbine temperatures.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. The only file written goes to tmp_path.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

R = impl.R_U
SPECIES = ["N2", "O2", "Ar", "CO2", "H2O", "H2"]


# ---------------------------------------------------------------- single species

@pytest.mark.parametrize("sp", SPECIES)
def test_cp_is_the_temperature_derivative_of_h_and_t_ds_is_cp_dt(sp):
    for T in (350.0, 800.0, 1500.0, 4000.0):
        dT = 1e-3 * T
        dhdT = (impl.h_molar(sp, T + dT) - impl.h_molar(sp, T - dT)) / (2 * dT)
        dsdT = (impl.s0_molar(sp, T + dT) - impl.s0_molar(sp, T - dT)) / (2 * dT)
        cp = impl.cp_molar(sp, T)
        # central difference error ~ (dT^2/6) h''' : relative 1e-6 is generous for smooth polynomials
        assert dhdT == pytest.approx(cp, rel=1e-6)
        assert T * dsdT == pytest.approx(cp, rel=1e-6)


def test_enthalpy_at_298_15_k_is_the_enthalpy_of_formation():
    # the database was fitted with R = 8.314510; with CODATA R the reproduction is good to about 2 J/mol
    for sp in SPECIES:
        assert impl.h_molar(sp, 298.15) == pytest.approx(impl.NASA9[sp]["Hf298"], abs=3.0)
    assert impl.h_molar("CO2", 298.15) < -3.9e5 and impl.h_molar("H2O", 298.15) < -2.4e5


@pytest.mark.parametrize("sp", SPECIES)
def test_the_two_fits_are_continuous_at_1000_k(sp):
    dcp, dh, ds = impl.fit_jump(sp)
    assert abs(dcp) < 1e-5 and abs(dh) < 1e-5 and abs(ds) < 1e-5
    # and fit_jump really compares the two intervals (forcing an interval changes the coefficients)
    assert not np.allclose(impl.coefficients(sp, 1000.0, 0)[0], impl.coefficients(sp, 1000.0, 1)[0]) or sp == "He"


def test_equipartition_limits():
    npt.assert_allclose(impl.cp_molar("He", np.array([300.0, 1000.0, 5000.0])) / R, 2.5, rtol=1e-12)
    assert impl.cp_molar("Ar", 500.0) / R == pytest.approx(2.5, rel=1e-12)
    # N2 near room temperature: 3 translational + 2 rotational, vibration frozen
    assert impl.cp_molar("N2", 300.0) / R == pytest.approx(3.5, rel=2e-3)
    # vibration wakes up: cp rises monotonically towards 9/2 R
    T = np.linspace(300.0, 3000.0, 28)
    cp = impl.cp_molar("N2", T) / R
    assert np.all(np.diff(cp) > 0) and 4.3 < cp[-1] < 4.5


def test_works_on_arrays_and_rejects_out_of_range_temperatures():
    T = np.array([300.0, 900.0, 1100.0, 2500.0])
    npt.assert_allclose(impl.cp_molar("O2", T), [impl.cp_molar("O2", t) for t in T])
    with pytest.raises(ValueError):
        impl.cp_molar("N2", 150.0)


# ---------------------------------------------------------------- dry air

def test_dry_air_molar_mass_and_gas_constant():
    assert sum(impl.AIR.values()) == pytest.approx(1.0, abs=1e-12)
    assert impl.M_AIR == pytest.approx(28.965e-3, rel=1e-4)
    assert impl.R_AIR == pytest.approx(287.05, abs=0.05)


def test_air_cp_against_lienhard_table_a6():
    cp = impl.air_cp(impl.LIENHARD_AIR_T)
    dev = (cp - impl.LIENHARD_AIR_CP) / impl.LIENHARD_AIR_CP
    # Lienhard quotes 0.2-2 % uncertainty for gas cp; ideal-gas NASA against 1-atm real air: within 0.5 %
    assert np.max(np.abs(dev)) < 5e-3
    assert np.max(np.abs(dev[impl.LIENHARD_AIR_T >= 400])) < 2e-3


def test_air_mixture_identities():
    T = np.array([300.0, 700.0, 1400.0])
    assert np.all(impl.air_gamma(T) < 1.4005) and np.all(np.diff(impl.air_gamma(T)) < 0)
    npt.assert_allclose(impl.air_h(T) - impl.air_u(T), impl.R_AIR * T, rtol=1e-12)        # h = u + p v = u + R T
    # entropy falls by R ln(p2/p1) at constant T; mixing entropy is positive
    assert impl.air_s(500.0, 2e5) - impl.air_s(500.0, 1e5) == pytest.approx(-impl.R_AIR * np.log(2), rel=1e-10)
    s_mix = impl.air_s(500.0, 1e5) - sum(x * impl.s0_molar(sp, 500.0) for sp, x in impl.AIR.items()) / impl.M_AIR
    assert s_mix > 0
    # a single-species "mixture" reduces to the species itself
    assert impl.mixture_cp({"N2": 1.0}, 600.0) == pytest.approx(impl.cp_molar("N2", 600.0) / impl.molar_mass("N2"))


# ---------------------------------------------------------------- equation of state

def test_ideal_gas_pressure_equals_the_kinetic_theory_pressure():
    for rho, T, M in ((1.225, 288.15, impl.M_AIR), (0.1, 216.65, impl.M_AIR), (40.0, 300.0, impl.molar_mass("He"))):
        p = impl.ideal_gas_pressure(rho, T, M)
        assert impl.kinetic_pressure(rho, T, M) == pytest.approx(p, rel=1e-12)
        assert p * impl.specific_volume(p, T, M) == pytest.approx(impl.R_U / M * T, rel=1e-12)    # p v = R T
    assert impl.ideal_gas_pressure(1.225, 288.15, impl.M_AIR) == pytest.approx(101325.0, rel=1e-3)


# ---------------------------------------------------------------- state postulate

def test_any_independent_pair_fixes_the_state():
    ref = impl.air_state(p=2.5e6, T=820.0)
    for pair in (("p", "v"), ("p", "h"), ("p", "s"), ("v", "T"), ("v", "s"), ("h", "s"), ("u", "v"), ("T", "s")):
        st = impl.air_state(**{k: ref[k] for k in pair})
        for k in ref:
            assert st[k] == pytest.approx(ref[k], rel=1e-8, abs=1e-6), (pair, k)


def test_dependent_or_wrong_count_pairs_are_refused():
    ref = impl.air_state(p=1e5, T=400.0)
    for pair in (("T", "h"), ("T", "u"), ("h", "u")):
        with pytest.raises(ValueError):
            impl.air_state(**{k: ref[k] for k in pair})
    with pytest.raises(ValueError):
        impl.air_state(p=1e5)
    with pytest.raises(ValueError):
        impl.air_state(p=1e5, T=400.0, v=1.0)


def test_temperature_inversions_round_trip():
    for T in (250.0, 999.0, 1001.0, 2200.0):
        assert impl.temperature_from_h(impl.air_h(T)) == pytest.approx(T, abs=1e-7)
        assert impl.temperature_from_s(impl.air_s(T, 3e5), 3e5) == pytest.approx(T, abs=1e-7)


# ---------------------------------------------------------------- the engineering comparison

def test_constant_cp_is_good_near_room_temperature_and_poor_at_turbine_inlet():
    # near 300 K the constant 1005 is within 0.5 % of the variable-cp enthalpy change
    dv = impl.delta_h_variable(280.0, 320.0)
    assert impl.delta_h_constant(280.0, 320.0) == pytest.approx(dv, rel=5e-3)
    # between compressor delivery and turbine inlet it underestimates by more than 10 %
    dv = impl.delta_h_variable(850.0, 1700.0)
    err = (impl.delta_h_constant(850.0, 1700.0) - dv) / dv
    assert -0.20 < err < -0.10
    # the variable-cp model itself is consistent: delta h = integral of cp dT (trapezoid on a fine grid)
    T = np.linspace(850.0, 1700.0, 2001)
    assert np.trapezoid(impl.air_cp(T), T) == pytest.approx(dv, rel=1e-7)


def test_exit_temperature_inverts_the_heat_addition():
    q = 6.0e5
    T2 = impl.exit_temperature(800.0, q)
    assert impl.delta_h_variable(800.0, T2) == pytest.approx(q, rel=1e-10)
    assert impl.exit_temperature(800.0, q, constant_cp=True) == pytest.approx(800.0 + q / 1005.0)
    assert impl.exit_temperature(800.0, q, constant_cp=True) > T2          # constant cp overpredicts the exit T
    assert impl.exit_temperature(800.0, 0.0) == pytest.approx(800.0, abs=1e-8)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "ideal_gas_properties.png").exists()
