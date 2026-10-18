"""Public tests for P10.01 Forces, Moments & Coefficients. They assert properties, not magic numbers:
the axis rotation is a rotation (preserves magnitude, inverts, reduces to identity at alpha = 0);
moment transfer is consistent (round trip, composition, a pure couple is unchanged); coefficients
are dimensionless and scale correctly; reducing noise-free SYNTHETIC readings recovers the generator's
truth; the aerodynamic centre does not depend on the reference point used to find it; the cm-vs-cn
slope changes by exactly the reference shift; first-order uncertainty propagation agrees with Monte
Carlo; and a small alpha error moves cd by about cn * d(alpha).

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`.
"""
from pathlib import Path

import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

DATA = Path(__file__).resolve().parents[1] / "data" / "SYNTHETIC_seed042_balance.csv"


# ---------------------------------------------------------------- similarity parameters

def test_dynamic_pressure_reynolds_mach_scale_correctly():
    assert impl.dynamic_pressure(1.225, 55.0) == pytest.approx(0.5 * 1.225 * 55.0 ** 2)
    assert impl.dynamic_pressure(1.0, 20.0) == pytest.approx(4 * impl.dynamic_pressure(1.0, 10.0))
    Re = impl.reynolds_number(55.0, 1.2, 1.461e-5)
    assert Re == pytest.approx(55.0 * 1.2 / 1.461e-5)
    assert impl.reynolds_number(110.0, 1.2, 1.461e-5) == pytest.approx(2 * Re)       # linear in V
    assert impl.mach_number(340.3, 340.3) == pytest.approx(1.0)


# ---------------------------------------------------------------- axes

def test_body_to_wind_is_a_rotation():
    rng = np.random.default_rng(3)
    N, A, a = rng.normal(50, 20, 40), rng.normal(0, 3, 40), rng.uniform(-0.3, 0.3, 40)
    L, D = impl.body_to_wind(N, A, a)
    npt.assert_allclose(np.hypot(L, D), np.hypot(N, A), rtol=1e-12)                # magnitude kept
    N2, A2 = impl.wind_to_body(L, D, a)
    npt.assert_allclose(N2, N, atol=1e-10)
    npt.assert_allclose(A2, A, atol=1e-10)
    L0, D0 = impl.body_to_wind(N, A, np.zeros_like(a))
    npt.assert_allclose(L0, N)
    npt.assert_allclose(D0, A)


def test_rotation_signs_follow_the_conventions():
    # pure normal force, nose-up alpha: lift positive, drag positive (the force tilts back)
    L, D = impl.body_to_wind(100.0, 0.0, np.radians(10.0))
    assert L > 0 and D > 0
    # leading-edge suction: negative axial force at alpha > 0 with zero drag
    a = np.radians(6.0)
    N, A = impl.wind_to_body(1.0, 0.0, a)
    assert A < 0 and A == pytest.approx(-np.sin(a))


# ---------------------------------------------------------------- moments

def test_moment_transfer_round_trip_and_composition():
    M, N, A = 2.0, 60.0, -3.0
    M1 = impl.transfer_moment(M, N, A, 0.12, -0.02, 0.075, 0.0)
    back = impl.transfer_moment(M1, N, A, 0.075, 0.0, 0.12, -0.02)
    assert back == pytest.approx(M)
    via = impl.transfer_moment(impl.transfer_moment(M, N, A, 0.12, -0.02, 0.3, 0.05), N, A, 0.3, 0.05, 0.075, 0.0)
    assert via == pytest.approx(M1)
    # a pure couple (no force) is the same about every point
    assert impl.transfer_moment(M, 0.0, 0.0, 0.0, 0.0, 5.0, 3.0) == pytest.approx(M)


def test_moment_transfer_signs():
    # an upward normal force at the quarter chord gives a nose-down moment about the leading edge
    assert impl.transfer_moment(0.0, 10.0, 0.0, 0.25, 0.0, 0.0, 0.0) == pytest.approx(-2.5)
    # an aft axial force acting below the reference point gives a nose-down moment
    assert impl.transfer_moment(0.0, 0.0, 4.0, 0.0, -0.5, 0.0, 0.0) == pytest.approx(-2.0)


def test_coefficients_are_dimensionless():
    q, S, c = 737.0, 0.18, 0.30
    cl, cd, cm = impl.coefficients(100.0, 2.0, -1.5, q, S, c)
    # doubling every force and q leaves the coefficients unchanged; doubling S and c (lengths) scales
    cl2, cd2, cm2 = impl.coefficients(200.0, 4.0, -3.0, 2 * q, S, c)
    npt.assert_allclose([cl2, cd2, cm2], [cl, cd, cm])
    assert cm == pytest.approx(-1.5 / (q * S * c))


# ---------------------------------------------------------------- reduction and the aerodynamic centre

@pytest.fixture(scope="module")
def clean():
    """Noise-free readings: the generator's true loads as the balance would read them."""
    d = impl.synthetic_balance_data(seed=7)
    return dict(N=d["N_true"], A=d["A_true"], M_bal=d["M_bal_true"], alpha_deg=d["alpha_true_deg"], q=d["q_true"])


def test_reduction_of_noise_free_readings_recovers_the_truth(clean):
    m = impl.MODEL
    S, c = m["c"] * m["span"], m["c"]
    r = impl.reduce_balance(clean["N"], clean["A"], clean["M_bal"], clean["alpha_deg"], clean["q"], S, c,
                            m["x_bal"], m["z_bal"])
    cl, cd, cm = impl.true_section(clean["alpha_deg"], 0.25)
    npt.assert_allclose(r["cl"], cl, atol=1e-12)
    npt.assert_allclose(r["cd"], cd, atol=1e-12)
    npt.assert_allclose(r["cm"], cm, atol=1e-12)


def test_aerodynamic_centre_is_independent_of_the_reference_point(clean):
    m = impl.MODEL
    S, c = m["c"] * m["span"], m["c"]
    xs = []
    for xr in (0.0, 0.25, 0.6, 1.0):
        r = impl.reduce_balance(clean["N"], clean["A"], clean["M_bal"], clean["alpha_deg"], clean["q"], S, c,
                                m["x_bal"], m["z_bal"], x_ref=xr * c)
        x_ac, cm_ac, k = impl.fit_aerodynamic_centre(r["cn"], r["cm"], xr * c, c)
        assert k == pytest.approx(xr - impl.TRUE_SECTION["x_ac"], abs=1e-10)        # slope = (x_ref - x_ac)/c
        assert cm_ac == pytest.approx(impl.TRUE_SECTION["cm_ac"], abs=1e-10)
        xs.append(x_ac)
    npt.assert_allclose(np.array(xs) / c, impl.TRUE_SECTION["x_ac"], atol=1e-10)


def test_ignoring_the_balance_offset_biases_the_moment_by_A_dz(clean):
    m = impl.MODEL
    S, c = m["c"] * m["span"], m["c"]
    args = (clean["N"], clean["A"], clean["M_bal"], clean["alpha_deg"], clean["q"], S, c, m["x_bal"])
    right = impl.reduce_balance(*args, m["z_bal"])
    wrong = impl.reduce_balance(*args, 0.0)
    npt.assert_allclose(wrong["cm"] - right["cm"], -m["z_bal"] * right["ca"] / c, atol=1e-12)
    npt.assert_allclose(wrong["cl"], right["cl"])                                      # forces unaffected


def test_noisy_reduction_is_unbiased_within_statistics():
    # 30 points per seed, 20 seeds: the mean cl error over 600 points has sigma ~ 0.0035/sqrt(600) = 1.4e-4;
    # 5 sigma = 7e-4. The x_ac fit scatter is ~1e-4 c per seed; the 20-seed mean is held to 3e-4 c.
    m = impl.MODEL
    S, c = m["c"] * m["span"], m["c"]
    errs, xacs = [], []
    for seed in range(20):
        d = impl.synthetic_balance_data(seed)
        r = impl.reduce_balance(d["N"], d["A"], d["M_bal"], d["alpha_deg"], d["q"], S, c, m["x_bal"], m["z_bal"])
        errs.append(r["cl"] - impl.true_section(d["alpha_true_deg"])[0])
        xacs.append(impl.fit_aerodynamic_centre(r["cn"], r["cm"], 0.25 * c, c)[0] / c)
    assert abs(np.mean(errs)) < 7e-4
    assert abs(np.mean(xacs) - impl.TRUE_SECTION["x_ac"]) < 3e-4


# ---------------------------------------------------------------- uncertainty

def test_linear_propagation_of_a_known_function():
    f = lambda x: np.array([x[0] * x[1], x[0] / x[1]])
    x, u = np.array([3.0, 2.0]), np.array([0.03, 0.01])
    uy, contrib = impl.propagate_uncertainty(f, x, u)
    npt.assert_allclose(uy[0], np.hypot(2.0 * 0.03, 3.0 * 0.01), rtol=1e-6)        # d(xy) = y dx + x dy
    npt.assert_allclose(contrib.sum(axis=1), uy ** 2, rtol=1e-12)
    assert contrib.shape == (2, 2)


def test_linear_propagation_agrees_with_monte_carlo():
    m = impl.MODEL
    S, c = m["c"] * m["span"], m["c"]
    d = impl.synthetic_balance_data(42)
    i = 10
    x = np.array([d["N"][i], d["A"][i], d["M_bal"][i], d["alpha_deg"][i], d["q"][i]])
    u = np.array([impl.NOISE["N"], impl.NOISE["A"], impl.NOISE["M"], impl.NOISE["alpha_deg"],
                  impl.NOISE["q_rel"] * d["q"][i]])
    f = lambda xx: impl.point_coefficients(xx, S, c, m["x_bal"], m["z_bal"], 0.25 * c)
    uy, _ = impl.propagate_uncertainty(f, x, u)
    mean, std = impl.monte_carlo_uncertainty(f, x, u, n=4000, seed=5)
    # with n = 4000 the sample standard deviation has relative sigma 1/sqrt(2n) = 1.1 %; allow 5 %
    npt.assert_allclose(std, uy, rtol=0.05)
    assert np.all(np.abs(mean - f(x)) < 4 * uy / np.sqrt(4000) + 1e-6)              # the model is nearly linear


def test_alpha_error_moves_cd_by_cn_dalpha():
    m = impl.MODEL
    S, c = m["c"] * m["span"], m["c"]
    d = impl.synthetic_balance_data(42)
    f = lambda xx: impl.point_coefficients(xx, S, c, m["x_bal"], m["z_bal"], 0.25 * c)
    for i in (4, 16, 26):
        x = np.array([d["N"][i], d["A"][i], d["M_bal"][i], d["alpha_deg"][i], d["q"][i]])
        da = 0.1
        dcd = f(x + np.array([0, 0, 0, da, 0]))[1] - f(x)[1]
        cl = f(x)[0]
        assert dcd == pytest.approx(cl * np.radians(da), rel=0.02)                  # dD/dalpha = L


# ---------------------------------------------------------------- data file and main

def test_synthetic_generator_is_seeded():
    a, b = impl.synthetic_balance_data(11), impl.synthetic_balance_data(11)
    npt.assert_array_equal(a["N"], b["N"])
    assert not np.allclose(impl.synthetic_balance_data(12)["N"], a["N"])


@pytest.mark.skipif(not DATA.exists(), reason="data file not shipped")
def test_shipped_csv_is_labelled_synthetic_and_reproducible():
    first = DATA.read_text().splitlines()[0]
    assert first.startswith("# SYNTHETIC: generated by synthetic_balance_data(seed=42)")
    d = impl.read_balance_csv(DATA)
    g = impl.synthetic_balance_data(42)
    for k in ("alpha_deg", "q", "N", "A", "M_bal"):
        npt.assert_allclose(d[k], g[k], atol=1e-6)


def test_main_writes_outputs(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "balance_coefficients.png").stat().st_size > 1000
    assert (tmp_path / "SYNTHETIC_seed042_balance.csv").read_text().startswith("# SYNTHETIC")
