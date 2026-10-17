"""Public tests for P04.07 Kinetic Theory of Gases. They assert the momentum-transfer pressure against
n k T, the Maxwell-Boltzmann density (normalisation, the three characteristic speeds, scaling with T and
M), equipartition and gamma, the speed of sound and its ratio to v_rms, the mean free path, the ISA
layers against the hydrostatic equation and standard values, and two seeded SIMULATIONS (box pressure,
relaxation by elastic collisions) with conservation checks and statistical tolerances.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. The only file written goes to tmp_path.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy.integrate import quad

from conftest import load_impl

impl = load_impl(__file__)

KB, NA, R = impl.K_B, impl.N_A, impl.R_GAS


# ---------------------------------------------------------------- pressure and the characteristic speeds

def test_momentum_transfer_pressure_equals_ideal_gas_pressure():
    T, p = 300.0, 101325.0
    n = impl.number_density(p, T)
    m = impl.M_N2 / NA
    assert impl.ideal_gas_pressure(n, T) == pytest.approx(p, rel=1e-12)
    assert impl.pressure_from_momentum_transfer(n, m, impl.rms_speed(T, impl.M_N2)) == pytest.approx(p, rel=1e-10)
    # doubling the rms speed quadruples the pressure (momentum per hit x hit rate)
    assert impl.pressure_from_momentum_transfer(n, m, 2 * 400.0) == pytest.approx(4 * impl.pressure_from_momentum_transfer(n, m, 400.0))


def test_characteristic_speeds_order_and_scaling():
    T, M = 300.0, impl.M_N2
    vp, vm, vr = impl.most_probable_speed(T, M), impl.mean_speed(T, M), impl.rms_speed(T, M)
    assert vp < vm < vr
    assert vm / vp == pytest.approx(np.sqrt(4 / np.pi))
    assert vr / vp == pytest.approx(np.sqrt(1.5))
    assert impl.rms_speed(4 * T, M) == pytest.approx(2 * vr)                    # ~ sqrt(T)
    assert impl.rms_speed(T, 4 * M) == pytest.approx(vr / 2)                    # ~ 1/sqrt(M)
    assert impl.rms_speed(T, impl.M_H2) > impl.rms_speed(T, impl.M_O2)
    assert 480 < impl.rms_speed(288.15, impl.M_AIR) < 520                       # about 500 m/s for air


def test_maxwell_boltzmann_is_normalised_with_the_right_moments():
    T, M = 350.0, impl.M_O2
    f = lambda v: impl.maxwell_boltzmann_pdf(v, T, M)
    assert quad(f, 0, np.inf)[0] == pytest.approx(1.0, rel=1e-8)
    assert quad(lambda v: v * f(v), 0, np.inf)[0] == pytest.approx(impl.mean_speed(T, M), rel=1e-8)
    assert quad(lambda v: v ** 2 * f(v), 0, np.inf)[0] == pytest.approx(impl.rms_speed(T, M) ** 2, rel=1e-8)
    v = np.linspace(1.0, 3000.0, 30000)
    assert v[np.argmax(f(v))] == pytest.approx(impl.most_probable_speed(T, M), rel=1e-3)
    assert f(0.0) == 0.0
    # equipartition: the mean kinetic energy is 3 k T / 2
    m = M / NA
    assert 0.5 * m * impl.rms_speed(T, M) ** 2 == pytest.approx(1.5 * KB * T)


def test_fraction_above_limits_and_the_escape_argument():
    T, M = 300.0, impl.M_N2
    assert impl.fraction_above(0.0, T, M) == pytest.approx(1.0, rel=1e-8)
    assert impl.fraction_above(1e5, T, M) < 1e-100
    assert impl.fraction_above(impl.mean_speed(T, M), T, M) == pytest.approx(1 - quad(lambda v: impl.maxwell_boltzmann_pdf(v, T, M), 0, impl.mean_speed(T, M))[0], rel=1e-8)
    # at 1000 K a hydrogen molecule exceeds the escape speed 10^80 times more often than nitrogen
    fh = impl.fraction_above(10800.0, 1000.0, impl.M_H2)
    fn = impl.fraction_above(10800.0, 1000.0, impl.M_N2)
    assert 1e-7 < fh < 1e-5
    assert fn < 1e-60


def test_sampled_speeds_match_the_distribution():
    rng = np.random.default_rng(71)
    T, M = 250.0, impl.M_HE
    s = impl.sample_speeds(200_000, T, M, rng)
    # with 2e5 samples the standard error of the mean is ~ 0.1 %; 0.5 % is a 5-sigma tolerance
    assert s.mean() == pytest.approx(impl.mean_speed(T, M), rel=5e-3)
    assert np.sqrt(np.mean(s ** 2)) == pytest.approx(impl.rms_speed(T, M), rel=5e-3)
    assert np.all(s >= 0.0)
    # fraction above v_rms matches the integral to ~1 %
    frac = np.mean(s > impl.rms_speed(T, M))
    assert frac == pytest.approx(impl.fraction_above(impl.rms_speed(T, M), T, M), abs=5e-3)


# ---------------------------------------------------------------- equipartition, gamma, sound

def test_gamma_and_heat_capacities_from_degrees_of_freedom():
    assert impl.gamma_from_dof(3) == pytest.approx(5 / 3)
    assert impl.gamma_from_dof(5) == pytest.approx(1.4)
    assert impl.gamma_from_dof(7) == pytest.approx(9 / 7)
    cv, cp = impl.molar_heat_capacities(5)
    assert cp - cv == pytest.approx(R)
    assert cp / cv == pytest.approx(impl.gamma_from_dof(5))
    assert cv == pytest.approx(2.5 * R)


def test_speed_of_sound_and_its_ratio_to_the_rms_speed():
    assert impl.speed_of_sound(288.15) == pytest.approx(340.3, abs=0.2)
    assert impl.speed_of_sound(216.65) == pytest.approx(295.1, abs=0.2)
    assert impl.speed_of_sound(4 * 288.15) == pytest.approx(2 * impl.speed_of_sound(288.15))
    T = 300.0
    assert impl.speed_of_sound(T) / impl.rms_speed(T, impl.M_AIR) == pytest.approx(impl.sound_to_rms_ratio(impl.GAMMA_AIR))
    assert impl.sound_to_rms_ratio(1.4) == pytest.approx(np.sqrt(1.4 / 3))
    assert impl.speed_of_sound(T, 5 / 3, impl.M_HE) > impl.speed_of_sound(T)   # helium: light and monatomic


# ---------------------------------------------------------------- mean free path

def test_mean_free_path_scaling_and_sea_level_value():
    T, p = 288.15, 101325.0
    lam = impl.mean_free_path(T, p)
    assert 5e-8 < lam < 8e-8                                                    # tens of nanometres
    assert impl.mean_free_path(T, p / 10) == pytest.approx(10 * lam)             # ~ 1/p
    assert impl.mean_free_path(2 * T, p) == pytest.approx(2 * lam)               # ~ T at fixed p
    assert impl.mean_free_path(T, p, 2 * impl.D_N2) == pytest.approx(lam / 4)    # ~ 1/d^2
    assert impl.collision_frequency(T, p) == pytest.approx(impl.mean_speed(T, impl.M_AIR) / lam)
    assert impl.knudsen_number(T, p, 1.0) == pytest.approx(lam)
    assert impl.knudsen_number(T, p, 1e-9) > 10.0                                # a nanometre-scale body is in free-molecule flow even at sea level


# ---------------------------------------------------------------- the ISA

def test_isa_standard_values_and_continuity():
    assert impl.isa_temperature(0.0) == pytest.approx(288.15)
    assert impl.isa_pressure(0.0) == pytest.approx(101325.0)
    assert impl.isa_density(0.0) == pytest.approx(1.225, rel=1e-3)
    assert impl.isa_temperature(11000.0) == pytest.approx(216.65)
    assert impl.isa_pressure(11000.0) == pytest.approx(22632.0, rel=1e-3)
    assert impl.isa_pressure(20000.0) == pytest.approx(5474.9, rel=1e-3)
    assert impl.isa_temperature(32000.0) == pytest.approx(228.65)
    assert impl.isa_pressure(32000.0) == pytest.approx(868.0, rel=1e-3)
    assert impl.isa_temperature(47000.0) == pytest.approx(270.65)
    assert impl.isa_temperature(71000.0) == pytest.approx(214.65)
    assert impl.isa_temperature(86000.0) == pytest.approx(184.65)
    h = np.linspace(0.0, 86000.0, 861)
    assert np.all(np.diff(impl.isa_pressure(h)) < 0)                            # pressure falls monotonically
    assert np.all(np.diff(impl.isa_density(h)) < 0)
    for hb in impl.ISA_BASES[1:-1]:                                              # continuous across every layer boundary
        assert impl.isa_pressure(hb - 1e-3) == pytest.approx(impl.isa_pressure(hb + 1e-3), rel=1e-6)
        assert impl.isa_temperature(hb - 1e-3) == pytest.approx(impl.isa_temperature(hb + 1e-3), abs=1e-4)
    assert np.isnan(impl.isa_pressure(90000.0))                                   # outside the model


def test_isa_satisfies_the_hydrostatic_equation_and_the_sound_speed_follows_t():
    h = np.array([500.0, 5000.0, 15000.0, 25000.0, 40000.0, 60000.0, 80000.0])
    npt.assert_allclose(impl.hydrostatic_residual(h), 0.0, atol=1e-5)
    npt.assert_allclose(impl.isa_speed_of_sound(h), np.sqrt(1.4 * impl.R_AIR * impl.isa_temperature(h)))
    assert impl.isa_speed_of_sound(11000.0) == pytest.approx(295.1, abs=0.2)
    # the isothermal layer 11-20 km: pressure falls exponentially with the scale height R T / g
    H = impl.R_AIR * 216.65 / impl.G0
    assert impl.isa_pressure(20000.0) / impl.isa_pressure(11000.0) == pytest.approx(np.exp(-9000.0 / H), rel=1e-9)


# ---------------------------------------------------------------- simulations

def test_box_simulation_reproduces_the_ideal_gas_pressure():
    # SIMULATION with a fixed seed; with 10^5 molecules the run-to-run relative scatter is ~0.5 % (measured over 40 seeds),
    # so 2 % is a 4-sigma tolerance.
    rng = np.random.default_rng(7)
    p_sim, p_ideal = impl.simulate_box_pressure(100_000, 300.0, impl.M_N2, 1e-6, 1e-7, rng)
    assert p_sim == pytest.approx(p_ideal, rel=0.02)
    assert p_ideal == pytest.approx(100_000 / 1e-18 * KB * 300.0)
    # the pressure is a statistical limit: the fluctuation with 100 molecules is much larger
    rng = np.random.default_rng(7)
    ratios = [impl.simulate_box_pressure(100, 300.0, impl.M_N2, 1e-6, 1e-7, rng)[0] / (100 / 1e-18 * KB * 300.0) for _ in range(8)]
    assert np.std(ratios) > 0.02                                                  # 100-molecule runs scatter by several percent


def test_collision_simulation_conserves_energy_and_momentum_and_relaxes_to_maxwell_boltzmann():
    rng = np.random.default_rng(11)
    N, T, M = 20_000, 300.0, impl.M_N2
    dirs = rng.normal(size=(N, 3))
    dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
    v0 = impl.rms_speed(T, M) * dirs                                              # all at one speed: far from equilibrium
    v1 = impl.relax_to_maxwellian(v0, 8 * N, rng)
    assert v1.shape == v0.shape
    assert np.sum(v1 ** 2) == pytest.approx(np.sum(v0 ** 2), rel=1e-12)          # kinetic energy exact
    npt.assert_allclose(v1.sum(axis=0), v0.sum(axis=0), atol=1e-6 * np.abs(v0).sum())   # momentum exact
    assert impl.temperature_from_velocities(v1, M) == pytest.approx(T, rel=1e-6)
    s = np.linalg.norm(v1, axis=1)
    # Maxwell-Boltzmann moments: <v> / v_rms = sqrt(8 / 3 pi) = 0.921, <v^4>/<v^2>^2 = 5/3; the single-speed
    # start has 1.0 and 1.0. With N = 2e4 the standard error of these ratios is ~0.3 %; 2 % is wide.
    assert s.mean() / impl.rms_speed(T, M) == pytest.approx(np.sqrt(8 / (3 * np.pi)), rel=0.02)
    assert np.mean(s ** 4) / np.mean(s ** 2) ** 2 == pytest.approx(5 / 3, rel=0.02)
    # a histogram against the density: chi-square-like check on 20 bins
    counts, edges = np.histogram(s, bins=20, range=(0.0, 1500.0))
    centres = 0.5 * (edges[1:] + edges[:-1])
    expected = N * impl.maxwell_boltzmann_pdf(centres, T, M) * (edges[1] - edges[0])
    big = expected > 50
    assert np.max(np.abs(counts[big] - expected[big]) / np.sqrt(expected[big])) < 5.0
    # the identity: zero collisions leaves the velocities unchanged
    npt.assert_array_equal(impl.relax_to_maxwellian(v0, 0, rng), v0)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "kinetic_theory.png").exists()
