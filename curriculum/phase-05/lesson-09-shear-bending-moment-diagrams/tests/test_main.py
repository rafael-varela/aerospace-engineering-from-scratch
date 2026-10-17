"""Public tests for P05.09 Shear & Bending Moment Diagrams. They assert physics, not magic numbers:
reactions satisfy sum F = 0 and sum M = 0 to round-off; V and M vanish beyond the right end
(the whole beam is in equilibrium); V jumps by P at a point force and M by -M0 at a couple;
dM/dx = V and dV/dx = w; the singularity-function diagram equals an independent method of
sections with quadrature at random cuts on random beams; the classic closed forms (w L^2 / 8,
P a b / L, P L, w L^2 / 2) and their locations come out of the extreme-moment search; symmetric
loads give symmetric M and antisymmetric V; superposition holds; indeterminate supports are
refused; the panelled elliptical lift converges to the closed-form root shear and moment.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

PF, PM, DL = impl.PointForce, impl.PointMoment, impl.DistributedLoad


def random_beam(rng, L=4.0):
    loads = []
    for _ in range(3):
        loads.append(PF(rng.uniform(0, L), rng.uniform(-5e3, 5e3)))
    a, b = np.sort(rng.uniform(0, L, 2))
    loads.append(DL(a, b, rng.uniform(-2e3, 2e3), rng.uniform(-2e3, 2e3)))
    loads.append(PM(rng.uniform(0, L), rng.uniform(-3e3, 3e3)))
    return loads


def resultant(loads):
    return np.sum([impl.load_resultant(l) for l in loads], axis=0)


# ---------------------------------------------------------------- singularity functions and resultants

def test_macaulay_brackets():
    x = np.array([0.0, 1.0, 2.0, 3.0])
    npt.assert_allclose(impl.macaulay(x, 1.5, 0), [0, 0, 1, 1])
    npt.assert_allclose(impl.macaulay(x, 1.5, 1), [0, 0, 0.5, 1.5])
    npt.assert_allclose(impl.macaulay(x, 1.5, 2), [0, 0, 0.25, 2.25])
    assert impl.macaulay(1.5, 1.5, 0) == 1.0                       # right-continuous step
    npt.assert_allclose(impl.macaulay(x, 1.0, -1), 0.0)             # delta: zero away from a


def test_distributed_load_resultant_and_centroid():
    l = DL(1.0, 4.0, 1000.0, 4000.0)                                 # trapezoid
    F, M = impl.load_resultant(l)
    assert F == pytest.approx(0.5 * (1000 + 4000) * 3)
    xbar = 1.0 + 3.0 * (1000 + 2 * 4000) / (3 * (1000 + 4000))       # trapezoid centroid
    assert M / F == pytest.approx(xbar)
    Ff, Mf = impl.load_resultant(impl.DistributedFunction(1.0, 4.0, lambda s: 1000 + 1000 * (s - 1)))
    assert Ff == pytest.approx(F) and Mf == pytest.approx(M)


# ---------------------------------------------------------------- reactions and equilibrium

def test_reactions_satisfy_equilibrium_on_random_beams():
    rng = np.random.default_rng(509)
    for supports in ([("pin", 0.0), ("roller", 4.0)], [("pin", 0.8), ("roller", 3.1)], [("fixed", 0.0)],
                     [("fixed", 4.0)]):
        for _ in range(5):
            loads = random_beam(rng)
            R = impl.solve_reactions(supports, loads)
            F, M = resultant(loads + R)
            assert abs(F) < 1e-9 * 1e4 and abs(M) < 1e-9 * 1e4
            # moments about any other point vanish too
            F_all = [impl.load_resultant(l) for l in loads + R]
            p = rng.uniform(-3, 7)
            assert abs(sum(m - p * f for f, m in F_all)) < 1e-7


def test_indeterminate_and_unstable_supports_are_refused():
    loads = [PF(1.0, -1e3)]
    with pytest.raises(ValueError):
        impl.solve_reactions([("pin", 0.0), ("roller", 2.0), ("roller", 4.0)], loads)
    with pytest.raises(ValueError):
        impl.solve_reactions([("fixed", 0.0), ("roller", 4.0)], loads)
    with pytest.raises(ValueError):
        impl.solve_reactions([("roller", 0.0)], loads)
    with pytest.raises(ValueError):
        impl.solve_reactions([("pin", 2.0), ("roller", 2.0)], loads)


def test_diagram_closes_beyond_the_right_end():
    rng = np.random.default_rng(1)
    for supports in ([("pin", 0.0), ("roller", 4.0)], [("pin", 1.0), ("roller", 3.0)], [("fixed", 0.0)]):
        loads = random_beam(rng)
        allr = loads + impl.solve_reactions(supports, loads)
        V, M = impl.shear_moment(np.array([4.0 + 1e-12, 5.0]), allr)
        npt.assert_allclose(V, 0.0, atol=1e-8)
        npt.assert_allclose(M, 0.0, atol=1e-7)
        V0, M0 = impl.shear_moment(np.array([-1e-9]), allr)            # nothing to the left of x = 0
        assert V0[0] == 0 and M0[0] == 0


# ---------------------------------------------------------------- independent check: method of sections

def test_singularity_functions_match_method_of_sections():
    rng = np.random.default_rng(2024)
    for supports in ([("pin", 0.0), ("roller", 4.0)], [("pin", 0.5), ("roller", 3.5)], [("fixed", 0.0)]):
        loads = random_beam(rng)
        allr = loads + impl.solve_reactions(supports, loads)
        for x in rng.uniform(0, 4.0, 12):
            V, M = impl.shear_moment(np.array([x]), allr)
            Vc, Mc = impl.sections_check(x, allr)
            assert V[0] == pytest.approx(Vc, abs=1e-6) and M[0] == pytest.approx(Mc, abs=1e-6)


def test_jumps_at_point_force_and_couple():
    loads = [PF(1.0, -2500.0), PM(2.5, 800.0), DL(0, 4.0, -300.0, -300.0)]
    allr = loads + impl.solve_reactions([("pin", 0.0), ("roller", 4.0)], loads)
    e = 1e-9
    V, M = impl.shear_moment(np.array([1.0 - e, 1.0, 2.5 - e, 2.5]), allr)
    assert V[1] - V[0] == pytest.approx(-2500.0, abs=1e-4)
    assert M[1] - M[0] == pytest.approx(0.0, abs=1e-4)              # M continuous at a force
    assert M[3] - M[2] == pytest.approx(-800.0, abs=1e-4)           # ccw couple: M drops by M0
    assert V[3] - V[2] == pytest.approx(0.0, abs=1e-4)              # V continuous at a couple


def test_derivative_relations():
    loads = [DL(0.5, 3.5, -1200.0, -200.0), PF(2.0, -1000.0)]
    allr = loads + impl.solve_reactions([("pin", 0.0), ("roller", 4.0)], loads)
    x = np.linspace(0.7, 1.9, 13)                                  # inside one smooth interval
    h = 1e-5
    Vp, Mp = impl.shear_moment(x + h, allr)
    Vm, Mm = impl.shear_moment(x - h, allr)
    V, _ = impl.shear_moment(x, allr)
    npt.assert_allclose((Mp - Mm) / (2 * h), V, rtol=1e-6, atol=1e-4)                    # dM/dx = V
    w = -1200.0 + 1000.0 * (x - 0.5) / 3.0
    npt.assert_allclose((Vp - Vm) / (2 * h), w, rtol=1e-6)                               # dV/dx = w


# ---------------------------------------------------------------- closed forms and extremes

def test_simply_supported_classics():
    L, w, P = 3.0, -1000.0, -4000.0
    l = [DL(0, L, w, w)]
    (x, M), _ = impl.extreme_moments(l + impl.solve_reactions([("pin", 0.0), ("roller", L)], l), L)
    assert x == pytest.approx(L / 2, abs=1e-9) and M == pytest.approx(-w * L ** 2 / 8)
    a = 1.2
    l = [PF(a, P)]
    (x, M), _ = impl.extreme_moments(l + impl.solve_reactions([("pin", 0.0), ("roller", L)], l), L)
    assert x == pytest.approx(a, abs=1e-6) and M == pytest.approx(-P * a * (L - a) / L, rel=1e-9)
    # triangular load 0 -> w0: M_max = w0 L^2 / (9 sqrt 3) at L / sqrt 3
    l = [DL(0, L, 0.0, w)]
    (x, M), _ = impl.extreme_moments(l + impl.solve_reactions([("pin", 0.0), ("roller", L)], l), L)
    assert x == pytest.approx(L / np.sqrt(3), rel=1e-9) and M == pytest.approx(-w * L ** 2 / (9 * np.sqrt(3)), rel=1e-9)


def test_cantilever_classics():
    L, w, P = 2.0, -800.0, -1500.0
    l = [PF(L, P), DL(0, L, w, w)]
    allr = l + impl.solve_reactions([("fixed", 0.0)], l)
    V0, M0 = impl.shear_moment(np.array([0.0]), allr)
    assert V0[0] == pytest.approx(-(P + w * L))                     # root carries all the load
    assert M0[0] == pytest.approx(P * L + w * L ** 2 / 2)           # hogging for downward loads
    _, (xmin, Mmin) = impl.extreme_moments(allr, L)
    assert xmin == pytest.approx(0.0, abs=1e-6) and Mmin == pytest.approx(M0[0])


def test_symmetry_and_superposition():
    L = 4.0
    sym = [PF(1.0, -1e3), PF(3.0, -1e3), DL(0, L, -500.0, -500.0)]
    allr = sym + impl.solve_reactions([("pin", 0.0), ("roller", L)], sym)
    x = np.linspace(0.05, 1.95, 20)
    V1, M1 = impl.shear_moment(x, allr)
    V2, M2 = impl.shear_moment(L - x, allr)
    npt.assert_allclose(M1, M2, atol=1e-8)                          # symmetric M
    npt.assert_allclose(V1, -V2, atol=1e-8)                         # antisymmetric V
    rng = np.random.default_rng(7)
    a, b = random_beam(rng), random_beam(rng)
    sup = [("pin", 0.0), ("roller", L)]
    xx = rng.uniform(0, L, 15)
    Va, Ma = impl.shear_moment(xx, a + impl.solve_reactions(sup, a))
    Vb, Mb = impl.shear_moment(xx, b + impl.solve_reactions(sup, b))
    Vab, Mab = impl.shear_moment(xx, a + b + impl.solve_reactions(sup, a + b))
    npt.assert_allclose(Va + Vb, Vab, atol=1e-8) and npt.assert_allclose(Ma + Mb, Mab, atol=1e-8)


def test_vanishing_load_gives_zero_diagram():
    R = impl.solve_reactions([("pin", 0.0), ("roller", 3.0)], [PF(1.0, 0.0)])
    V, M = impl.shear_moment(np.linspace(0, 3, 7), [PF(1.0, 0.0)] + R)
    npt.assert_allclose(V, 0.0) and npt.assert_allclose(M, 0.0)


# ---------------------------------------------------------------- the wing spar

def test_elliptical_lift_panels_converge_to_closed_form_root_loads():
    wc = impl.WingCase()
    V_ex, M_ex = wc.exact_lift_root()
    assert V_ex == pytest.approx(0.5 * wc.n * wc.mass * impl.G0)       # half the lift
    assert M_ex / V_ex == pytest.approx(4 * wc.semispan / (3 * np.pi))  # centre of an ellipse quadrant
    errs = []
    for n in (20, 80, 320):
        l = wc.loads(n, relief=False)
        V0, M0 = impl.shear_moment(np.array([0.0]), l + impl.solve_reactions([("fixed", 0.0)], l))
        assert V0[0] < 0 and M0[0] > 0                                  # upward lift: sagging root moment
        errs.append(abs(M0[0] - M_ex) / M_ex)
    orders = np.log(np.array(errs[:-1]) / errs[1:]) / np.log(4)
    assert np.all(orders > 1.3)                                         # ~1.5: sqrt behaviour at the tip
    assert errs[-1] < 3e-4


def test_inertia_relief_reduces_root_loads_and_jumps_at_engine():
    wc = impl.WingCase()
    l = wc.loads(200)
    allr = l + impl.solve_reactions([("fixed", 0.0)], l)
    lo = wc.loads(200, relief=False)
    V0, M0 = impl.shear_moment(np.array([0.0]), allr)
    V1, M1 = impl.shear_moment(np.array([0.0]), lo + impl.solve_reactions([("fixed", 0.0)], lo))
    relief = sum(impl.load_resultant(p)[0] for p in wc.inertia_loads())
    assert V0[0] - V1[0] == pytest.approx(-relief, rel=1e-9)           # less shear by the inertia loads
    assert 0 < M0[0] < M1[0]
    e = 1e-9
    V, _ = impl.shear_moment(np.array([wc.engine_x - e, wc.engine_x]), allr)
    assert V[1] - V[0] == pytest.approx(-wc.n * wc.engine_mass * impl.G0, rel=1e-6)
    # moving the engine outboard relieves more root moment
    l2 = impl.WingCase(engine_x=3.0).loads(200)
    _, M2 = impl.shear_moment(np.array([0.0]), l2 + impl.solve_reactions([("fixed", 0.0)], l2))
    assert M2[0] < M0[0]
    # the exact elliptic lift by sections agrees with the panelled singularity solution
    ex = [impl.DistributedFunction(0.0, wc.semispan, wc.lift)] + wc.inertia_loads()
    allx = ex + impl.solve_reactions([("fixed", 0.0)], ex)
    Vc, Mc = impl.sections_check(1.0, allx)
    V, M = impl.shear_moment(np.array([1.0]), allr)
    assert M[0] == pytest.approx(Mc, rel=1e-3) and V[0] == pytest.approx(Vc, rel=1e-3)


def test_floor_beam_reactions_and_extreme():
    R = impl.solve_reactions([("pin", 0.0), ("roller", impl.FLOOR_SPAN)], impl.FLOOR_LOADS)
    total_down = 4 * 2.0e3 + 500.0 * impl.FLOOR_SPAN
    assert R[0].P + R[1].P == pytest.approx(total_down)
    (x, M), _ = impl.extreme_moments(impl.FLOOR_LOADS + R, impl.FLOOR_SPAN)
    xs = np.linspace(0, impl.FLOOR_SPAN, 20001)
    _, Ms = impl.shear_moment(xs, impl.FLOOR_LOADS + R)
    assert M >= Ms.max() - 1e-6 and M == pytest.approx(Ms.max(), rel=1e-6)


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "vm_diagrams.png").stat().st_size > 10_000
