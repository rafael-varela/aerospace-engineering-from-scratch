"""P09.09 Root Locus: reference solution.

An altitude-hold outer loop. An inner pitch-attitude loop (P09.11 designs one) makes pitch follow its
command with a lag; the flight path follows pitch with a further lag; altitude is the integral of
climb rate V gamma. A representative low-order model (typical of a light aircraft or a small UAV at
V = 50 m/s; not one vehicle's data):

    h(s) / theta_c(s) = G(s) = V / (s (T s + 1)(tau s + 1)),   T = 1.0 s (flight-path lag),
                                                               tau = 0.5 s (closed inner pitch loop).

A proportional altitude controller theta_c = K (h_c - h), K in rad per metre, closes the loop. The
closed-loop poles are the roots of den(s) + K num(s) = 0, and as K runs from 0 to infinity they trace
the root locus. This module computes it and its construction rules from scratch:

* the closed-loop roots against K, tracked branch by branch;
* the angle condition arg G(s) = 180 deg (mod 360) and the magnitude condition K |G(s)| = 1;
* the real-axis segments, the asymptotes and their centroid, the breakaway points (dK/ds = 0);
* the Routh array (with the epsilon and auxiliary-polynomial special cases), the number of
  right-half-plane roots from its first column, and the j omega crossings and their gain;
* the gain that puts a complex pair on a damping-ratio line (zeta = 0.5) by searching along the line
  for the angle condition;
* a from-scratch step simulation to see what the dominant-pair prediction misses.

Units: altitude in m, pitch in rad, K in rad/m, time in s, frequencies in rad/s.

Run:  python solution.py [out_dir]     (prints the lesson's tables, writes root_locus.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-09/lesson-09-root-locus
"""
from __future__ import annotations

import itertools
import sys
from pathlib import Path

import numpy as np

V_AIR = 50.0      # m/s
T_PATH = 1.0      # s, flight-path lag (theta -> gamma)
TAU_PITCH = 0.5   # s, closed inner pitch loop


# ---------------------------------------------------------------- given helpers

def poly_add(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Sum of two polynomials (highest power first) of possibly different length."""
    p, q = np.atleast_1d(np.asarray(p, dtype=float)), np.atleast_1d(np.asarray(q, dtype=float))
    n = max(p.size, q.size)
    return np.concatenate([np.zeros(n - p.size), p]) + np.concatenate([np.zeros(n - q.size), q])


def tf_eval(num: np.ndarray, den: np.ndarray, s: np.ndarray | complex) -> np.ndarray:
    """num(s)/den(s) at complex s."""
    return np.polyval(num, s) / np.polyval(den, s)


def altitude_plant(v: float = V_AIR, t_path: float = T_PATH, tau: float = TAU_PITCH) -> tuple[np.ndarray, np.ndarray]:
    """(num, den) of G(s) = V / (s (T s + 1)(tau s + 1))."""
    den = np.convolve(np.convolve([1.0, 0.0], [t_path, 1.0]), [tau, 1.0])
    return np.array([v]), den


def match_roots(prev: np.ndarray, new: np.ndarray) -> np.ndarray:
    """Reorder `new` so that new[i] is the root closest to prev[i] (exhaustive over permutations for
    up to 7 roots, greedy beyond). Keeps each locus branch continuous as K changes."""
    n = len(new)
    if n <= 7:
        best, best_cost = None, np.inf
        for perm in itertools.permutations(range(n)):
            cost = np.sum(np.abs(prev - new[list(perm)]))
            if cost < best_cost:
                best, best_cost = perm, cost
        return new[list(best)]
    out = np.empty_like(new)
    free = list(range(n))
    for i in range(n):
        j = min(free, key=lambda k: abs(prev[i] - new[k]))
        out[i] = new[j]
        free.remove(j)
    return out


def tf_to_companion(b, a):
    """(A, B, C, D) of the controllable companion form of a proper b(s)/a(s)."""
    b, a = np.atleast_1d(np.asarray(b, float)), np.atleast_1d(np.asarray(a, float))
    n = a.size - 1
    b = np.concatenate([np.zeros(a.size - b.size), b]) / a[0]
    a = a / a[0]
    A = np.zeros((n, n))
    A[0, :] = -a[1:]
    A[1:, :-1] = np.eye(n - 1)
    B = np.zeros((n, 1))
    B[0, 0] = 1.0
    return A, B, (b[1:] - b[0] * a[1:]).reshape(1, n), float(b[0])


def step_response(num: np.ndarray, den: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Unit-step response of num/den on the uniform grid t (exact for a step: zero-order-hold steps;
    the matrix exponential by scaling and squaring, as in P02.05)."""
    A, B, C, D = tf_to_companion(num, den)
    n = A.shape[0]
    dt = t[1] - t[0]
    M = np.zeros((n + 1, n + 1))
    M[:n, :n], M[:n, n:] = A, B
    M = M * dt
    s = int(np.ceil(np.log2(max(np.linalg.norm(M, 1), 1e-300) / 0.5))) if np.linalg.norm(M, 1) > 0.5 else 0
    X = M / 2.0 ** s
    term, E = np.eye(n + 1), np.eye(n + 1)
    for k in range(1, 20):
        term = term @ X / k
        E = E + term
    for _ in range(s):
        E = E @ E
    Ad, bd = E[:n, :n], E[:n, n]
    y = np.empty(t.size)
    x = np.zeros(n)
    for k in range(t.size):
        y[k] = C[0] @ x + D
        x = Ad @ x + bd
    return y


# ---------------------------------------------------------------- the locus itself

def char_poly(num: np.ndarray, den: np.ndarray, K: float) -> np.ndarray:
    """Closed-loop characteristic polynomial den(s) + K num(s) (highest power first)."""
    return poly_add(den, K * np.asarray(num, dtype=float))


def locus_roots(num: np.ndarray, den: np.ndarray, gains: np.ndarray) -> np.ndarray:
    """Closed-loop roots for each gain, shape (len(gains), deg den); each column is one continuous branch
    (gains should increase in small steps; the first row is sorted by np.sort_complex)."""
    gains = np.asarray(gains, dtype=float)
    out = np.empty((gains.size, len(den) - 1), dtype=complex)
    prev = None
    for i, K in enumerate(gains):
        r = np.roots(char_poly(num, den, K)).astype(complex)
        r = np.sort_complex(r) if prev is None else match_roots(prev, r)
        out[i] = r
        prev = r
    return out


def angle_condition_error(num: np.ndarray, den: np.ndarray, s: np.ndarray | complex) -> np.ndarray:
    """Deviation of arg G(s) from 180 deg, wrapped to (-180, 180] deg. Zero on the locus (K > 0)."""
    ang = np.degrees(np.angle(tf_eval(num, den, s)))
    return (ang - 180.0 + 180.0) % 360.0 - 180.0


def gain_at_point(num: np.ndarray, den: np.ndarray, s: np.ndarray | complex) -> np.ndarray:
    """The gain that puts a closed-loop root at s, from the magnitude condition K = 1 / |G(s)|."""
    return 1.0 / np.abs(tf_eval(num, den, s))


def real_axis_segments(num: np.ndarray, den: np.ndarray) -> list[tuple[float, float]]:
    """Intervals (lo, hi) of the real axis on the locus for K > 0: points with an odd number of real
    open-loop poles plus zeros to their right. lo may be -inf."""
    crit = np.concatenate([np.roots(den), np.roots(num) if len(num) > 1 else []])
    real = np.sort(np.real(crit[np.abs(np.imag(crit)) < 1e-9]))[::-1]   # descending
    segs = []
    for i in range(0, len(real), 2):
        hi = real[i]
        lo = real[i + 1] if i + 1 < len(real) else -np.inf
        segs.append((float(lo), float(hi)))
    return segs[::-1]


def asymptotes(num: np.ndarray, den: np.ndarray) -> tuple[float, np.ndarray]:
    """Centroid sigma_a = (sum poles - sum zeros)/(n - m) and asymptote angles (2k + 1) 180/(n - m) deg."""
    p = np.roots(den)
    z = np.roots(num) if len(num) > 1 else np.array([])
    q = len(p) - len(z)
    if q <= 0:
        raise ValueError("no asymptotes: as many zeros as poles")
    sigma = float(np.real(np.sum(p) - np.sum(z)) / q)
    return sigma, (2 * np.arange(q) + 1) * 180.0 / q


def breakaway_points(num: np.ndarray, den: np.ndarray) -> np.ndarray:
    """Real breakaway and break-in points: real roots of num(s) den'(s) - num'(s) den(s) = 0 at which
    K = -den/num is positive (so the point lies on the K > 0 locus)."""
    num, den = np.asarray(num, float), np.asarray(den, float)
    d = poly_add(np.convolve(num, np.polyder(den)), -np.convolve(np.polyder(num) if num.size > 1 else [0.0], den))
    cand = np.roots(np.trim_zeros(d, "f"))
    out = []
    for s in cand:
        if abs(s.imag) < 1e-7:
            K = -np.polyval(den, s.real) / np.polyval(num, s.real)
            if K > 0:
                out.append(s.real)
    return np.sort(np.array(out))


# ---------------------------------------------------------------- Routh

def routh_array(coeffs: np.ndarray, eps: float = 1e-9) -> np.ndarray:
    """Routh array of the polynomial (highest power first), rows s^n ... s^0.

    A zero in the first column (rest of the row nonzero) is replaced by eps; a whole row of zeros is
    replaced by the coefficients of the derivative of the auxiliary polynomial formed from the row above.
    """
    c = np.atleast_1d(np.asarray(coeffs, dtype=float))
    n = c.size - 1
    width = n // 2 + 1
    R = np.zeros((n + 1, width))
    R[0, :len(c[0::2])] = c[0::2]
    R[1, :len(c[1::2])] = c[1::2]
    for i in range(2, n + 1):
        if np.all(np.abs(R[i - 1]) < 1e-14):
            # row of zeros: derivative of the auxiliary polynomial from row i-2 (power n - i + 2)
            power = n - i + 2
            aux_pows = power - 2 * np.arange(width)
            R[i - 1] = np.where(aux_pows > 0, R[i - 2] * aux_pows, 0.0)
        if abs(R[i - 1, 0]) < 1e-14:
            R[i - 1, 0] = eps
        for j in range(width - 1):
            R[i, j] = (R[i - 1, 0] * R[i - 2, j + 1] - R[i - 2, 0] * R[i - 1, j + 1]) / R[i - 1, 0]
    return R


def routh_rhp_count(coeffs: np.ndarray) -> int:
    """Number of roots in the open right half plane: sign changes in the first column of the Routh array."""
    col = routh_array(coeffs)[:, 0]
    s = np.sign(col)
    s[s == 0] = 1
    return int(np.sum(s[1:] != s[:-1]))


def jw_crossings(num: np.ndarray, den: np.ndarray) -> list[tuple[float, float]]:
    """Gains K > 0 and frequencies omega >= 0 at which a closed-loop root sits on the imaginary axis.

    On s = j omega, den + K num = 0 needs Im[den(j omega) conj(num(j omega))] = 0 (K real); that is a
    real polynomial in omega, solved here; then K = -den/num must be real and positive.
    """
    num, den = np.asarray(num, float), np.asarray(den, float)

    def at_jw(p: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        # p(j w) = R(w) + j I(w) as real polynomials in w (highest power first)
        deg = p.size - 1
        Rw, Iw = np.zeros(deg + 1), np.zeros(deg + 1)
        for i, ci in enumerate(p):
            k = deg - i
            val = ci * (1j ** k)
            Rw[i] += val.real
            Iw[i] += val.imag
        return Rw, Iw

    Dr, Di = at_jw(den)
    Nr, Ni = at_jw(num)
    poly = poly_add(np.convolve(Di, Nr), -np.convolve(Dr, Ni))      # Im[D conj(N)]
    poly = np.trim_zeros(poly, "f")
    out = []
    for w in np.roots(poly) if poly.size > 1 else []:
        if abs(w.imag) < 1e-7 and w.real >= -1e-12:
            wr = abs(w.real)
            Kc = -np.polyval(den, 1j * wr) / np.polyval(num, 1j * wr)
            if abs(Kc.imag) < 1e-6 * max(1.0, abs(Kc)) and Kc.real > 0:
                out.append((float(Kc.real), float(wr)))
    out = sorted(set((round(k, 12), round(w, 12)) for k, w in out))
    return [(float(k), float(w)) for k, w in out]


def gain_for_damping(num: np.ndarray, den: np.ndarray, zeta: float, wn_max: float = 100.0, n_scan: int = 4000) -> tuple[float, complex]:
    """Smallest K > 0 placing a closed-loop root on the ray s = wn (-zeta + j sqrt(1 - zeta^2)).

    Scans wn for a sign change of the angle-condition error (ignoring the +-180 deg wrap jumps), refines
    by bisection, and returns (K, s) with K from the magnitude condition.
    """
    d = np.array(-zeta + 1j * np.sqrt(1.0 - zeta ** 2))
    wn = np.linspace(wn_max / n_scan, wn_max, n_scan)
    err = angle_condition_error(num, den, wn * d)
    for i in range(wn.size - 1):
        if np.sign(err[i]) != np.sign(err[i + 1]) and abs(err[i] - err[i + 1]) < 90.0:
            lo, hi = wn[i], wn[i + 1]
            elo = err[i]
            for _ in range(100):
                mid = 0.5 * (lo + hi)
                em = angle_condition_error(num, den, mid * d)
                if np.sign(em) == np.sign(elo):
                    lo, elo = mid, em
                else:
                    hi = mid
            s = complex(0.5 * (lo + hi) * d)
            return float(gain_at_point(num, den, s)), s
    raise ValueError("the locus does not cross this damping line for wn < wn_max")


# ---------------------------------------------------------------- the lesson's numbers

def overshoot_and_settling(y: np.ndarray, t: np.ndarray, band: float = 0.02) -> tuple[float, float]:
    """Percent overshoot over the final value y[-1] and the 2 % settling time of a step response."""
    yf = y[-1]
    os_ = 100.0 * (y.max() - yf) / yf
    outside = np.flatnonzero(np.abs(y - yf) > band * abs(yf))
    ts = t[outside[-1] + 1] if outside.size else 0.0
    return float(os_), float(ts)


def pd_design_table(tds=(0.0, 1.0, 2.0, 4.0), zeta: float = 0.5) -> list[dict]:
    """For C = K (Td s + 1) on the altitude plant: the gain for zeta, the predicted dominant pair, all
    closed-loop poles, and the simulated overshoot and settling of a unit altitude step."""
    num_p, den_p = altitude_plant()
    rows = []
    t = np.arange(0.0, 60.0, 0.01)
    for td in tds:
        num = np.convolve(num_p, [td, 1.0]) if td > 0 else num_p
        K, s = gain_for_damping(num, den_p, zeta)
        cl_den = char_poly(num, den_p, K)
        poles = np.roots(cl_den)
        y = step_response(K * num, cl_den, t)
        os_, ts = overshoot_and_settling(y, t)
        rows.append({"td": td, "K": K, "s": s, "poles": poles, "pred_os": 100 * np.exp(-np.pi * zeta / np.sqrt(1 - zeta ** 2)),
                     "pred_ts": 4.0 / (zeta * abs(s)), "sim_os": os_, "sim_ts": ts,
                     "jw": jw_crossings(num, den_p)})
    return rows


def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's tables and save root_locus.png in out_dir."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    num, den = altitude_plant()
    print(f"G(s) = {V_AIR:.0f} / (s ({T_PATH} s + 1)({TAU_PITCH} s + 1)); den =", np.array2string(den, precision=3),
          " poles", np.array2string(np.roots(den).real, precision=3))
    sig, ang = asymptotes(num, den)
    print(f"real-axis segments {real_axis_segments(num, den)}; centroid {sig:.4f}, asymptote angles {ang}")
    nd = poly_add(np.convolve(num, np.polyder(den)), -np.convolve([0.0], den))
    print("breakaway candidates (roots of N D' - N' D):", np.array2string(np.roots(np.trim_zeros(nd, 'f')).real, precision=4),
          "-> on the locus:", np.array2string(breakaway_points(num, den), precision=4),
          f"with K = {float(-np.polyval(den, breakaway_points(num, den)[0]) / np.polyval(num, breakaway_points(num, den)[0])):.6f} rad/m")
    jw = jw_crossings(num, den)
    print(f"j omega crossing: K = {jw[0][0]:.6f} rad/m at omega = {jw[0][1]:.6f} rad/s")
    for K in (0.5 * jw[0][0], 0.999 * jw[0][0], 1.001 * jw[0][0], 2 * jw[0][0]):
        print(f"  Routh at K = {K:.5f}: first column {np.array2string(routh_array(char_poly(num, den, K))[:, 0], precision=4)}, "
              f"RHP roots {routh_rhp_count(char_poly(num, den, K))} (np.roots: {int(np.sum(np.roots(char_poly(num, den, K)).real > 0))})")
    gains = np.concatenate([[0.0], np.logspace(-5, 0, 600)])
    R = locus_roots(num, den, gains)
    err = np.array([np.max(np.abs(angle_condition_error(num, den, R[i]))) for i in range(1, gains.size)])
    mag = np.array([np.max(np.abs(gains[i] * np.abs(tf_eval(num, den, R[i])) - 1)) for i in range(1, gains.size)])
    print(f"{R.size - R.shape[1]} locus points: max angle-condition error {err.max():.2e} deg, max |K|G| - 1| {mag.max():.2e}")
    print(f"sum of closed-loop roots: min {np.sum(R, axis=1).real.min():.6f}, max {np.sum(R, axis=1).real.max():.6f} (= -a2/a3 = {-den[1] / den[0]:.1f})")

    print("\ngain for zeta = 0.5 (P control, then K (Td s + 1)); simulated 1 m altitude step:")
    for r in pd_design_table():
        others = [p for p in r["poles"] if abs(p - r["s"]) > 1e-6 and abs(p - np.conj(r["s"])) > 1e-6]
        print(f"  Td = {r['td']:.1f} s: K = {r['K']:.6f} rad/m ({np.degrees(r['K']):.3f} deg/m), pair {r['s'].real:.4f} +- {abs(r['s'].imag):.4f}j "
              f"(wn {abs(r['s']):.4f}), other poles {np.array2string(np.array(others), precision=4)}; predicted OS {r['pred_os']:.1f} %, "
              f"ts {r['pred_ts']:.2f} s; simulated OS {r['sim_os']:.1f} %, ts {r['sim_ts']:.2f} s; j-omega crossings {r['jw']}")

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8))
    ax = axes[0]
    for j in range(R.shape[1]):
        ax.plot(R[:, j].real, R[:, j].imag, lw=1.2)
    ax.plot(np.roots(den).real, np.roots(den).imag, "kx", ms=8)
    for a_ in ang:
        rr = np.linspace(0, 4, 2)
        ax.plot(sig + rr * np.cos(np.radians(a_)), rr * np.sin(np.radians(a_)), "k:", lw=0.7)
    wl = np.linspace(0, 3, 2)
    ax.plot(-0.5 * wl, np.sqrt(0.75) * wl, "g--", lw=0.8, label="zeta = 0.5")
    r0 = pd_design_table((0.0,))[0]
    ax.plot([r0["s"].real], [r0["s"].imag], "ro", ms=5, label=f"K = {r0['K']:.4f} rad/m")
    ax.axhline(0, color="k", lw=0.4); ax.axvline(0, color="k", lw=0.4)
    ax.set_xlim(-4, 1.5); ax.set_ylim(-3, 3); ax.set_aspect("equal")
    ax.set_xlabel("Re s [rad/s]"); ax.set_ylabel("Im s [rad/s]"); ax.legend(fontsize=7)
    ax.set_title("altitude hold, P control: root locus", fontsize=9)
    t = np.arange(0.0, 40.0, 0.01)
    for r in pd_design_table():
        n_ = np.convolve(num, [r["td"], 1.0]) if r["td"] > 0 else num
        axes[1].plot(t, step_response(r["K"] * n_, char_poly(n_, den, r["K"]), t), label=f"Td = {r['td']:.0f} s")
    axes[1].axhline(1.0, color="k", lw=0.5)
    axes[1].set_xlabel("t [s]"); axes[1].set_ylabel("altitude [m] for a 1 m command"); axes[1].legend(fontsize=7); axes[1].grid(lw=0.3)
    axes[1].set_title("gain chosen for zeta = 0.5", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_dir / "root_locus.png", dpi=130)
    plt.close(fig)
    print(f"saved {out_dir / 'root_locus.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
