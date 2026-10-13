"""Public tests for P01.12 Divergence, Stokes & Green's Theorems. They assert identities and convergence, not magic numbers.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-01/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest
import sympy as sp

from conftest import load_impl

impl = load_impl(__file__)

x, y, z = impl.coords()
BOX = ((-0.2, 0.2), (-0.2, 0.2), (1.0, 2.0))
RNG = np.random.default_rng(1211)
PTS = RNG.uniform(-1.0, 1.0, size=(12, 3))

# A smooth test field with non-trivial divergence and curl, and its exact derivatives.
F_SYM = [y * sp.exp(x), z * sp.sin(y), x * sp.cos(z)]
PHI_SYM = x**2 * y + sp.sin(z) * x - y**3


# ---------------------------------------------------------------- symbolic operators and identities

def test_symbolic_divergence_and_curl_of_a_known_field():
    assert sp.simplify(impl.divergence_sym(F_SYM) - (y * sp.exp(x) + z * sp.cos(y) - x * sp.sin(z))) == 0
    c = impl.curl_sym(F_SYM)
    expected = [-sp.sin(y), -sp.cos(z), -sp.exp(x)]
    assert all(sp.simplify(a - b) == 0 for a, b in zip(c, expected))
    with pytest.raises(ValueError):
        impl.divergence_sym([x, y])


def test_curl_of_a_gradient_is_zero_and_divergence_of_a_curl_is_zero():
    assert all(sp.simplify(c) == 0 for c in impl.curl_sym(impl.gradient_sym(PHI_SYM)))
    assert sp.simplify(impl.divergence_sym(impl.curl_sym(F_SYM))) == 0
    # Laplacian via div(grad) of x**2 + y**2 + z**2 is 6.
    assert sp.simplify(impl.divergence_sym(impl.gradient_sym(x**2 + y**2 + z**2)) - 6) == 0


def test_lambdified_field_broadcasts_constants_to_the_input_shape():
    f = impl.lambdify_field([sp.Integer(1), x, sp.Integer(0)])
    u, v, w = f(np.zeros((3, 2)), np.zeros((3, 2)), np.zeros((3, 2)))
    assert u.shape == v.shape == w.shape == (3, 2)
    npt.assert_allclose(u, 1.0)
    npt.assert_allclose(w, 0.0)
    g = impl.lambdify_scalar(sp.Integer(3))
    assert g(np.zeros(5), np.zeros(5), np.zeros(5)).shape == (5,)


# ---------------------------------------------------------------- finite-difference operators

def test_finite_difference_divergence_and_curl_match_sympy_and_converge_at_second_order():
    # F_SYM exercises the divergence; G_SYM = (cos y, sin z, e^x) has a curl whose every term is non-linear
    # in the differentiated coordinate, so the central difference is not accidentally exact.
    G_SYM = [sp.cos(y), sp.sin(z), sp.exp(x)]
    f, g = impl.lambdify_field(F_SYM), impl.lambdify_field(G_SYM)
    div_exact = impl.lambdify_scalar(impl.divergence_sym(F_SYM))(*PTS.T)
    curl_exact = np.stack(impl.lambdify_field(impl.curl_sym(G_SYM))(*PTS.T), axis=-1)
    errs_div, errs_curl = [], []
    for h in (4e-2, 2e-2, 1e-2):
        errs_div.append(np.max(np.abs(impl.divergence_fd(f, *PTS.T, h=h) - div_exact)))
        errs_curl.append(np.max(np.abs(impl.curl_fd(g, *PTS.T, h=h) - curl_exact)))
    for errs in (np.array(errs_div), np.array(errs_curl)):
        rates = np.log2(errs[:-1] / errs[1:])
        assert np.all(rates > 1.8), f"expected second order, got {rates}"
    npt.assert_allclose(impl.divergence_fd(f, *PTS.T, h=1e-4), div_exact, atol=1e-7)
    npt.assert_allclose(impl.curl_fd(g, *PTS.T, h=1e-4), curl_exact, atol=1e-7)
    assert impl.curl_fd(f, *PTS.T).shape == (12, 3)
    # div G = 0 identically; the finite difference must agree to rounding.
    npt.assert_allclose(impl.divergence_fd(g, *PTS.T), 0.0, atol=1e-9)
    with pytest.raises(ValueError):
        impl.divergence_fd(f, 0.0, 0.0, 0.0, h=0.0)


def test_example_fields_have_the_advertised_divergence_and_curl():
    k, omega = 10.0, 0.5
    flow = lambda X, Y, Z: impl.contraction_flow(X, Y, Z, k)  # noqa: E731
    rot = lambda X, Y, Z: impl.rigid_rotation(X, Y, Z, omega)  # noqa: E731
    exp_ = lambda X, Y, Z: impl.uniform_expansion(X, Y, Z, k)  # noqa: E731
    npt.assert_allclose(impl.divergence_fd(flow, *PTS.T), 0.0, atol=1e-8)        # incompressible
    npt.assert_allclose(impl.curl_fd(flow, *PTS.T), 0.0, atol=1e-8)              # irrotational
    npt.assert_allclose(impl.divergence_fd(rot, *PTS.T), 0.0, atol=1e-8)
    npt.assert_allclose(impl.curl_fd(rot, *PTS.T), [[0.0, 0.0, 2 * omega]] * len(PTS), atol=1e-8)
    npt.assert_allclose(impl.divergence_fd(exp_, *PTS.T), 3 * k, rtol=1e-8)
    # Gravity inside a uniform sphere: div g = -4 pi G rho (Gauss's law for gravity).
    rho = 5515.0
    grav = lambda X, Y, Z: impl.uniform_sphere_gravity(X, Y, Z, rho)  # noqa: E731
    npt.assert_allclose(impl.divergence_fd(grav, *PTS.T, h=1.0), -4 * np.pi * impl.G_NEWTON * rho, rtol=1e-6)


# ---------------------------------------------------------------- divergence theorem on a box

def test_face_fluxes_of_the_contraction_flow_balance_exactly():
    k = 10.0
    flow = lambda X, Y, Z: impl.contraction_flow(X, Y, Z, k)  # noqa: E731
    faces = impl.box_face_fluxes(flow, BOX, 8)
    assert set(faces) == {"x0", "x1", "y0", "y1", "z0", "z1"}
    assert all(faces[f] < 0.0 for f in ("x0", "x1", "y0", "y1", "z0"))   # inflow through the sides and bottom
    assert faces["z1"] > 0.0                                             # outflow through the top
    assert faces["z1"] == pytest.approx(k * 2.0 * 0.16)                  # w = k z at z = 2 over 0.16 m**2
    assert abs(impl.box_flux(flow, BOX, 8)) < 1e-12 * abs(faces["z1"])  # net zero: incompressible


def test_divergence_theorem_holds_for_a_uniform_expansion_and_is_linear():
    k = 2.0
    exp_ = lambda X, Y, Z: impl.uniform_expansion(X, Y, Z, k)  # noqa: E731
    three_k = lambda X, Y, Z: np.full(np.shape(X), 3 * k)  # noqa: E731
    flux, vol, diff = impl.divergence_theorem_check(exp_, three_k, BOX, 5)
    volume = 0.4 * 0.4 * 1.0
    assert flux == pytest.approx(3 * k * volume, rel=1e-12)
    assert vol == pytest.approx(3 * k * volume, rel=1e-12)
    assert diff < 1e-12
    twice = lambda X, Y, Z: impl.uniform_expansion(X, Y, Z, 2 * k)  # noqa: E731
    assert impl.box_flux(twice, BOX, 5) == pytest.approx(2 * flux, rel=1e-12)


def test_divergence_theorem_residual_converges_at_second_order_for_a_smooth_field():
    f = impl.lambdify_field(F_SYM)
    d = impl.lambdify_scalar(impl.divergence_sym(F_SYM))
    res = impl.divergence_theorem_convergence(f, d, BOX, [4, 8, 16, 32])
    assert np.all(res > 0.0)
    rates = np.log2(res[:-1] / res[1:])
    assert np.all(rates > 1.8), f"expected second order, got {rates}"
    flux, vol, _ = impl.divergence_theorem_check(f, d, BOX, 32)
    assert flux == pytest.approx(vol, rel=1e-4)


def test_discrete_divergence_theorem_is_exact_when_each_component_is_quadratic_in_its_own_coordinate():
    Q = [x**2 * y, y * z**2, sp.sin(x) * z]
    flux, vol, diff = impl.divergence_theorem_check(impl.lambdify_field(Q), impl.lambdify_scalar(impl.divergence_sym(Q)), BOX, 3)
    assert diff < 1e-12 * max(1.0, abs(flux))


def test_flux_of_a_curl_through_a_closed_box_is_zero():
    curl_f = impl.lambdify_field(impl.curl_sym(F_SYM))
    scale = abs(impl.box_face_fluxes(curl_f, BOX, 16)["z1"])
    assert abs(impl.box_flux(curl_f, BOX, 16)) < 1e-10 * max(scale, 1.0)


def test_box_volume_integral_is_exact_for_linear_integrands_and_rejects_bad_input():
    assert impl.box_volume_integral(lambda X, Y, Z: 2 * X + Z, BOX, 3) == pytest.approx(0.16 * 1.5)   # mean of 2x is 0, of z is 1.5
    with pytest.raises(ValueError):
        impl.box_volume_integral(lambda X, Y, Z: X, BOX, 0)
    with pytest.raises(ValueError):
        impl.box_flux(lambda X, Y, Z: (X, Y, Z), ((0, 1), (1, 0), (0, 1)), 4)


# ---------------------------------------------------------------- Green's theorem

def test_greens_theorem_for_rigid_rotation_gives_twice_the_area_exactly():
    rot2d = lambda X, Y: (-Y, X)  # noqa: E731
    two = lambda X, Y: np.full(np.shape(X), 2.0)  # noqa: E731
    circ, area, diff = impl.greens_theorem_check(rot2d, two, -1.0, 3.0, 0.5, 2.0, 7)
    assert circ == pytest.approx(2 * 4.0 * 1.5, rel=1e-12)
    assert diff < 1e-12


def test_greens_theorem_for_a_general_smooth_field_converges_at_second_order():
    P = lambda X, Y: (np.sin(X) * Y**2, np.exp(X) * np.cos(Y))  # noqa: E731
    curl_z = lambda X, Y: np.exp(X) * np.cos(Y) - 2 * np.sin(X) * Y  # noqa: E731
    diffs = np.array([impl.greens_theorem_check(P, curl_z, 0.0, 1.5, -0.5, 1.0, n)[2] for n in (8, 16, 32, 64)])
    rates = np.log2(diffs[:-1] / diffs[1:])
    assert np.all(rates > 1.8), f"expected second order, got {rates}"


def test_greens_theorem_for_the_rankine_vortex_inside_outside_and_around_the_core():
    gamma, a = 300.0, 1.5
    rv = lambda X, Y: impl.rankine_vortex_2d(X, Y, gamma, a)  # noqa: E731
    wz = lambda X, Y: impl.rankine_vorticity(X, Y, gamma, a)  # noqa: E731
    circ, area, diff = impl.greens_theorem_check(rv, wz, -0.5, 0.5, -0.5, 0.5, 50)       # inside: uniform vorticity
    assert circ == pytest.approx(gamma / (np.pi * a**2) * 1.0, rel=1e-12)
    assert diff < 1e-10
    circ, area, diff = impl.greens_theorem_check(rv, wz, 2.0, 3.0, 2.0, 3.0, 50)         # outside: irrotational
    assert abs(circ) < 1e-3 and abs(area) < 1e-12               # midpoint error only; vorticity exactly zero
    circ, area, diff = impl.greens_theorem_check(rv, wz, -3.0, 3.0, -3.0, 3.0, 400)      # enclosing: all of gamma
    assert circ == pytest.approx(gamma, rel=1e-5)
    assert area == pytest.approx(gamma, rel=5e-3)        # step in vorticity: only first-order accurate


def test_polygon_circulation_flips_sign_with_orientation_and_rejects_bad_input():
    rot2d = lambda X, Y: (-Y, X)  # noqa: E731
    pts = impl.rectangle_boundary(0.0, 1.0, 0.0, 1.0, 3)
    assert impl.polygon_circulation(rot2d, pts[::-1]) == pytest.approx(-impl.polygon_circulation(rot2d, pts))
    with pytest.raises(ValueError):
        impl.polygon_circulation(rot2d, np.zeros((2, 2)))


# ---------------------------------------------------------------- Stokes' theorem

def test_hemisphere_normal_points_outward():
    hemi = impl.hemisphere_surface(1.0)
    pos = lambda X, Y, Z: (X, Y, Z)  # noqa: E731
    flux = impl.surface_flux(pos, hemi, (0.0, np.pi / 2), (0.0, 2 * np.pi), 30, 30)
    assert flux == pytest.approx(2 * np.pi, rel=1e-3)          # r . n = R on the surface; area 2 pi R**2


def test_stokes_theorem_for_rigid_rotation_on_the_hemisphere_and_the_disc():
    omega, R = 0.5, 2.0
    rot = lambda X, Y, Z: impl.rigid_rotation(X, Y, Z, omega)  # noqa: E731
    curl_rot = lambda X, Y, Z: (np.zeros(np.shape(X)), np.zeros(np.shape(X)), np.full(np.shape(X), 2 * omega))  # noqa: E731
    circ, cap, disc = impl.stokes_theorem_check(rot, curl_rot, R, 40)
    exact = 2 * np.pi * omega * R**2
    assert circ == pytest.approx(exact, rel=1e-3)
    assert cap == pytest.approx(exact, rel=1e-3)
    assert disc == pytest.approx(exact, rel=1e-6)
    assert cap == pytest.approx(disc, rel=1e-3)                  # same rim, same flux of curl


def test_stokes_theorem_for_a_general_field_and_its_convergence():
    f = impl.lambdify_field(F_SYM)
    c = impl.lambdify_field(impl.curl_sym(F_SYM))
    R = 1.3
    errs_cap = []
    for n in (10, 20, 40):
        circ, cap, disc = impl.stokes_theorem_check(f, c, R, n)
        errs_cap.append(abs(cap - disc))
    errs_cap = np.array(errs_cap)
    circ, cap, disc = impl.stokes_theorem_check(f, c, R, 80)
    assert circ == pytest.approx(disc, rel=1e-4)
    assert cap == pytest.approx(disc, rel=1e-4)
    rates = np.log2(errs_cap[:-1] / errs_cap[1:])
    assert np.all(rates > 1.5), f"expected roughly second order, got {rates}"


def test_rim_circulation_of_a_gradient_field_is_zero():
    grad = impl.lambdify_field(impl.gradient_sym(PHI_SYM))
    assert abs(impl.rim_circulation(grad, 1.7, 400)) < 1e-6
    with pytest.raises(ValueError):
        impl.rim_circulation(grad, 1.0, 2)
