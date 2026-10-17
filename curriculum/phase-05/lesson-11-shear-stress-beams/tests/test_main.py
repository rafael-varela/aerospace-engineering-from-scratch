"""Public tests for P05.11 Shear Stress in Beams. They assert mechanics, not magic numbers: the
shear stress in a rectangle integrates to V and peaks at 1.5 V/A; the I-section web carries almost all
of V; the thin-walled shear flow is continuous at junctions, vanishes at free edges, has resultant
(Sx, Sy) and is linear in the load; it agrees with an independent route (longitudinal equilibrium of
a slice, integrating the change of bending stress) and with closed forms for the I-beam flanges and
web; the channel's shear centre matches 3 b^2 t_f/(6 b t_f + h t_w) and lies outside the section;
symmetric and antisymmetric sections put the shear centre where symmetry says.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)


# ---------------------------------------------------------------- solid sections

def test_rectangle_shear_stress_integrates_to_V_and_peaks_at_three_halves():
    V, b, h = 1.0e4, 0.02, 0.06
    y = np.linspace(-h / 2, h / 2, 4001)
    tau = impl.rect_shear_stress(V, b, h, y)
    assert np.trapezoid(tau * b, y) == pytest.approx(V, rel=1e-6)
    assert tau.max() == pytest.approx(1.5 * V / (b * h))
    assert tau[0] == pytest.approx(0.0, abs=1e-6) and tau[-1] == pytest.approx(0.0, abs=1e-6)
    # VQ/It at the neutral axis, built by hand
    assert impl.vq_it(V, b * h / 2 * h / 4, b * h ** 3 / 12, b) == pytest.approx(1.5 * V / (b * h))


def test_circle_peak_is_four_thirds_of_average():
    assert impl.circle_max_shear(3.0, 0.5) == pytest.approx(4 / 3 * 3.0 / (np.pi * 0.25))
    # independent: VQ/It with Q = 2 r^3 / 3 and t = 2 r
    r, V = 0.01, 500.0
    assert impl.vq_it(V, 2 * r ** 3 / 3, np.pi * r ** 4 / 4, 2 * r) == pytest.approx(impl.circle_max_shear(V, r))


def test_i_section_web_carries_almost_all_the_shear_with_a_jump_at_the_flange():
    V, h, b, tf, tw = 13974.5, 0.180, 0.060, 0.010, 0.003
    hw = h - 2 * tf
    y = np.linspace(-hw / 2, hw / 2, 4001)
    web = np.trapezoid(impl.i_section_shear_stress(V, h, b, tf, tw, y) * tw, y)
    yf = np.linspace(hw / 2 + 1e-12, h / 2, 801)          # just inside the flange
    fl = 2 * np.trapezoid(impl.i_section_shear_stress(V, h, b, tf, tw, yf) * b, yf)
    assert web + fl == pytest.approx(V, rel=1e-5)                       # total = V
    assert 0.9 < web / V < 1.0                                          # web carries > 90 %
    jump = impl.i_section_shear_stress(V, h, b, tf, tw, hw / 2 - 1e-12) / impl.i_section_shear_stress(V, h, b, tf, tw, hw / 2 + 1e-12)
    assert jump == pytest.approx(b / tw, rel=1e-6)                      # same Q, width ratio
    # web average V / (h_w t_w) is within 1 % of the peak for this deep, thin-webbed spar
    assert impl.i_section_shear_stress(V, h, b, tf, tw, 0.0) == pytest.approx(V / (hw * tw), rel=0.01)


# ---------------------------------------------------------------- thin-walled properties

def test_thin_wall_properties_of_channel_match_hand_formulas():
    h, b, t = 0.1, 0.04, 0.002
    p = impl.thin_wall_properties(impl.channel_segments(h, b, t, t))
    assert p["A"] == pytest.approx(t * (h + 2 * b))
    assert p["xc"] == pytest.approx(b ** 2 / (h + 2 * b))
    assert p["Ixx"] == pytest.approx(t * h ** 3 / 12 + 2 * b * t * (h / 2) ** 2)
    assert abs(p["Ixy"]) < 1e-18
    z = impl.thin_wall_properties(impl.z_segments(0.04, 0.019, 0.002))
    assert abs(z["xc"]) < 1e-15 and abs(z["yc"]) < 1e-15
    # flanges in the first and third quadrants: each gives t (h/2) b^2 / 2 > 0
    assert z["Ixy"] == pytest.approx(2 * 0.002 * 0.02 * 0.019 ** 2 / 2, rel=1e-12)


# ---------------------------------------------------------------- shear flow

@pytest.mark.parametrize("Sx,Sy", [(0.0, 1000.0), (500.0, 0.0), (300.0, -700.0)])
def test_shear_flow_resultant_equals_applied_shear_and_free_edges_are_zero(Sx, Sy):
    for segs in (impl.i_beam_segments(0.17, 0.06, 0.01, 0.003), impl.channel_segments(0.1, 0.04, 0.002, 0.0015),
                 impl.z_segments(0.04, 0.019, 0.002)):
        flows = impl.open_section_shear_flow(segs, Sx, Sy)
        Fx, Fy, _ = impl.shear_flow_resultant(segs, flows)
        scale = abs(Sx) + abs(Sy)
        assert Fx == pytest.approx(Sx, abs=1e-9 * scale)
        assert Fy == pytest.approx(Sy, abs=1e-9 * scale)
        assert abs(flows[0][1][0]) < 1e-12 * scale                      # starts at a free edge
        assert abs(flows[-1][1][-1]) < 1e-6 * scale                     # ends at a free edge: q back to 0


def test_shear_flow_is_continuous_at_junctions_and_linear_in_load():
    segs = impl.i_beam_segments(0.17, 0.06, 0.01, 0.003)
    f1 = impl.open_section_shear_flow(segs, 0.0, 1000.0)
    # top junction: two half-flanges flow in, web flows out
    assert f1[2][1][0] == pytest.approx(f1[0][1][-1] + f1[1][1][-1])
    # bottom junction: web + bottom-left half-flange flow into bottom-right half-flange
    assert f1[4][1][0] == pytest.approx(f1[2][1][-1] + f1[3][1][-1])
    f2 = impl.open_section_shear_flow(segs, 0.0, 2500.0)
    for (s1, q1), (s2, q2) in zip(f1, f2):
        npt.assert_allclose(q2, 2.5 * q1, rtol=1e-12, atol=1e-9)


def test_i_beam_shear_flow_matches_closed_forms():
    hc, b, tf, tw, S = 0.17, 0.06, 0.01, 0.003, 13974.5
    segs = impl.i_beam_segments(hc, b, tf, tw)
    I = tw * hc ** 3 / 12 + 2 * b * tf * (hc / 2) ** 2               # thin-wall I by hand
    flows = impl.open_section_shear_flow(segs, 0.0, S)
    s, q = flows[0]
    npt.assert_allclose(np.abs(q), S * tf * s * (hc / 2) / I, rtol=1e-10, atol=1e-6)   # linear in the flange
    s, q = flows[2]
    y = hc / 2 - s
    q_web = S / I * (b * tf * hc / 2 + tw / 2 * ((hc / 2) ** 2 - y ** 2))             # parabola in the web
    npt.assert_allclose(np.abs(q), q_web, rtol=1e-10)
    assert np.all(q < 0)              # segment runs downward, flow runs upward: the web pushes up with Sy > 0


def test_shear_flow_matches_independent_bending_gradient_route():
    for segs, Sx, Sy in ((impl.i_beam_segments(0.17, 0.06, 0.01, 0.003), 0.0, 13974.5),
                         (impl.channel_segments(0.1, 0.04, 0.002, 0.0015), 400.0, 900.0),
                         (impl.z_segments(0.04, 0.019, 0.002), -250.0, 600.0)):
        f = impl.open_section_shear_flow(segs, Sx, Sy)
        g = impl.shear_flow_from_bending_gradient(segs, Sx, Sy, 0.05, 501)
        scale = max(np.max(np.abs(q)) for _, q in f)
        for (fs, fq), (gs, gq) in zip(f, g):
            npt.assert_allclose(np.interp(fs, gs, gq), fq, atol=1e-8 * scale)


# ---------------------------------------------------------------- shear centre

@pytest.mark.parametrize("h,b,tf,tw", [(0.1, 0.04, 0.002, 0.002), (0.08, 0.03, 0.0025, 0.0016), (0.2, 0.02, 0.003, 0.003)])
def test_channel_shear_centre_matches_closed_form_and_lies_outside(h, b, tf, tw):
    segs = impl.channel_segments(h, b, tf, tw)
    xs, ys = impl.shear_centre(segs)
    e = impl.channel_shear_centre_closed_form(h, b, tf, tw)
    assert xs == pytest.approx(-e, rel=1e-9)
    assert abs(ys) < 1e-9 * h                                           # on the axis of symmetry
    assert xs < 0 < impl.thin_wall_properties(segs)["xc"]               # web between shear centre and centroid
    assert e < b                                                        # within one flange length of the web


def test_shear_centre_on_symmetry_point_for_i_and_z_sections():
    xs, ys = impl.shear_centre(impl.i_beam_segments(0.17, 0.06, 0.01, 0.003))
    assert abs(xs) < 1e-12 and abs(ys) < 1e-12
    xs, ys = impl.shear_centre(impl.z_segments(0.04, 0.019, 0.002))
    assert abs(xs) < 1e-12 and abs(ys) < 1e-12


def test_shear_centre_moment_balance_for_channel():
    h, b, t = 0.1, 0.04, 0.002
    segs = impl.channel_segments(h, b, t, t)
    xs, ys = impl.shear_centre(segs)
    flows = impl.open_section_shear_flow(segs, 0.0, 1000.0)
    _, _, Mz = impl.shear_flow_resultant(segs, flows, about=(xs, ys))
    assert abs(Mz) < 1e-9 * 1000.0 * h                                   # no moment about the shear centre
    _, _, Mw = impl.shear_flow_resultant(segs, flows, about=(0.0, 0.0))
    # by hand: only the flanges have a moment arm about the web centre
    I = t * h ** 3 / 12 + 2 * b * t * (h / 2) ** 2
    F_flange = 1000.0 * t * b ** 2 * h / (4 * I)
    assert abs(Mw) == pytest.approx(F_flange * h, rel=1e-9)


# ---------------------------------------------------------------- fasteners

def test_fastener_pitch_scales_inversely_with_shear_flow():
    assert impl.fastener_pitch(1.0e5, 4000.0) == pytest.approx(0.04)
    assert impl.fastener_pitch(2.0e5, 4000.0) == pytest.approx(0.02)
    assert impl.fastener_pitch(-1.0e5, 4000.0, 1.5) == pytest.approx(0.04 / 1.5)


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "shear_flow.png").stat().st_size > 10_000
