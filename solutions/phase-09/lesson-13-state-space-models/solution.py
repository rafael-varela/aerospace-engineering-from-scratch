"""P09.13 State-Space Models: reference solution.

State-space models from ODEs and from transfer functions, and back: the controllable (companion) and
observable canonical forms, the transfer function of any (A, B, C, D) by the Faddeev-LeVerrier algorithm,
similarity transforms and the real modal form, and simulation by the zero-order-hold matrix exponential.
The aerospace model is the longitudinal short period of a small aircraft extended with pitch angle and
altitude, written directly from its physical ODEs, whose altitude transfer function is also derived by
hand; the two must agree, and the simulation must agree with ``scipy.signal.lti``.

Conventions:
    * Single input, single output unless stated; B is (n, 1), C is (1, n), D is (1, 1).
    * Polynomials are NumPy coefficient arrays, highest power first.
    * Pitch-altitude state x = (alpha, q, theta, h): rad, rad/s, rad, m. Input elevator delta [rad],
      positive trailing edge up (nose-up moment).

Run:  python solution.py          (prints the checks, writes state_space.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-09/lesson-13-state-space-models
"""
from __future__ import annotations

from pathlib import Path

import numpy as np


# ---------------------------------------------------------------- given

def pitch_altitude_params() -> dict[str, float]:
    """Return a representative (not aircraft-specific) small-UAV short-period model plus airspeed V = 30 m/s.

    Z_alpha [1/s], M_alpha [1/s^2], M_q [1/s], M_delta [1/s^2] as in P09.11 (no servo here).
    """
    return {"Z_alpha": -1.2, "M_alpha": -6.0, "M_q": -2.5, "M_delta": 12.0, "V": 30.0}


def expm(M: np.ndarray, n_terms: int = 20) -> np.ndarray:
    """Return e^M by scaling and squaring a truncated Taylor series (the method of P02.05)."""
    M = np.asarray(M, dtype=float)
    norm = np.linalg.norm(M, 1)
    s = max(0, int(np.ceil(np.log2(norm / 0.5)))) if norm > 0.5 else 0
    X = M / 2 ** s
    term = np.eye(len(M))
    E = np.eye(len(M))
    for k in range(1, n_terms):
        term = term @ X / k
        E = E + term
    for _ in range(s):
        E = E @ E
    return E


def round_sig(x: np.ndarray, digits: int) -> np.ndarray:
    """Round every entry of x to ``digits`` significant figures (models coefficients typed into flight software)."""
    x = np.asarray(x, dtype=float)
    out = np.zeros_like(x)
    nz = x != 0
    mag = np.floor(np.log10(np.abs(x[nz])))
    out[nz] = np.round(x[nz] / 10 ** mag, digits - 1) * 10 ** mag
    return out


# ---------------------------------------------------------------- transfer function <-> state space

def tf_eval(num: np.ndarray, den: np.ndarray, s: complex | np.ndarray) -> complex | np.ndarray:
    """Return num(s)/den(s)."""
    return np.polyval(num, s) / np.polyval(den, s)


def ss_eval(A: np.ndarray, B: np.ndarray, C: np.ndarray, D: np.ndarray, s: complex | np.ndarray) -> complex | np.ndarray:
    """Return G(s) = C (sI - A)^{-1} B + D (SISO) at a scalar s or at each entry of an array of s."""
    s_arr = np.atleast_1d(np.asarray(s, dtype=complex))
    n = A.shape[0]
    out = np.array([(C @ np.linalg.solve(si * np.eye(n) - A, B) + D)[0, 0] for si in s_arr])
    return out[0] if np.ndim(s) == 0 else out


def tf_to_controllable(num: np.ndarray, den: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return (A, B, C, D) of the controllable canonical (companion) realisation of a proper num/den.

    With den monic, den = s^n + a1 s^{n-1} + ... + an and num padded to length n + 1:
    D = num[0], b = num[1:] - D a, A = [[-a1, -a2, ..., -an], [1, 0, ...], ..., [0, ..., 1, 0]], B = e1, C = b.
    State x1 is the highest derivative of the internal variable; this is the form ``scipy.signal.tf2ss`` returns.
    """
    num = np.atleast_1d(np.asarray(num, dtype=float))
    den = np.atleast_1d(np.asarray(den, dtype=float))
    if len(num) > len(den):
        raise ValueError("improper transfer function: deg num > deg den")
    num, den = num / den[0], den / den[0]
    n = len(den) - 1
    num = np.concatenate([np.zeros(n + 1 - len(num)), num])
    D = num[0]
    b = num[1:] - D * den[1:]
    A = np.zeros((n, n))
    A[0, :] = -den[1:]
    A[1:, :-1] = np.eye(n - 1)
    B = np.zeros((n, 1))
    B[0, 0] = 1.0
    return A, B, b[None, :].copy(), np.array([[D]])


def tf_to_observable(num: np.ndarray, den: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return the observable canonical realisation: the dual (A^T, C^T, B^T, D) of the controllable form."""
    A, B, C, D = tf_to_controllable(num, den)
    return A.T.copy(), C.T.copy(), B.T.copy(), D.copy()


def faddeev_leverrier(A: np.ndarray) -> tuple[np.ndarray, list[np.ndarray]]:
    """Return (c, [N_1, ..., N_n]) with det(sI - A) = sum c_k s^{n-k} (c_0 = 1) and adj(sI - A) = sum N_k s^{n-k}.

    Recursion: N_1 = I, c_1 = -tr(A); N_k = A N_{k-1} + c_{k-1} I, c_k = -tr(A N_k)/k. Uses only matrix products
    and traces; Cayley-Hamilton guarantees A N_n + c_n I = 0, which the tests use as a residual.
    """
    A = np.asarray(A, dtype=float)
    n = A.shape[0]
    c = np.zeros(n + 1)
    c[0] = 1.0
    N = np.eye(n)
    Ns = [N]
    c[1] = -np.trace(A @ N)
    for k in range(2, n + 1):
        N = A @ N + c[k - 1] * np.eye(n)
        Ns.append(N)
        c[k] = -np.trace(A @ N) / k
    return c, Ns


def ss_to_tf(A: np.ndarray, B: np.ndarray, C: np.ndarray, D: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (num, den) of C (sI - A)^{-1} B + D from Faddeev-LeVerrier: den = det(sI - A), num = sum C N_k B s^{n-k} + D den."""
    c, Ns = faddeev_leverrier(A)
    n = A.shape[0]
    num = np.zeros(n + 1)
    for k, N in enumerate(Ns, start=1):
        num[k] = (C @ N @ B)[0, 0]
    num = num + D[0, 0] * c
    return num, c


# ---------------------------------------------------------------- changes of coordinates

def similarity(A: np.ndarray, B: np.ndarray, C: np.ndarray, D: np.ndarray, T: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return the model in new coordinates z with x = T z: (T^{-1} A T, T^{-1} B, C T, D)."""
    Ti = np.linalg.inv(T)
    return Ti @ A @ T, Ti @ B, C @ T, D.copy()


def modal_form(A: np.ndarray, B: np.ndarray, C: np.ndarray, D: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return (Am, Bm, Cm, Dm, T): the real modal form and the transform x = T z (A diagonalisable).

    Each real eigenvalue contributes its eigenvector as a column of T; each complex pair sigma +- j omega with
    eigenvector v = v_r + j v_i contributes the two columns (v_r, v_i), giving the real 2x2 block
    [[sigma, omega], [-omega, sigma]] on the diagonal of Am. Modes are ordered by increasing |Re lambda|.
    """
    lam, V = np.linalg.eig(A)
    order = np.argsort(np.abs(lam.real) + 1e-9 * np.abs(lam.imag))
    lam, V = lam[order], V[:, order]
    cols, used = [], np.zeros(len(lam), dtype=bool)
    for i, li in enumerate(lam):
        if used[i]:
            continue
        if abs(li.imag) < 1e-12 * max(1.0, abs(li)):
            cols.append(V[:, i].real)
            used[i] = True
        else:
            j = next(k for k in range(len(lam)) if not used[k] and k != i and abs(lam[k] - li.conjugate()) < 1e-8 * max(1.0, abs(li)))
            v = V[:, i] if li.imag > 0 else V[:, j]
            cols.extend([v.real, v.imag])
            used[i] = used[j] = True
    T = np.column_stack(cols)
    Am, Bm, Cm, Dm = similarity(A, B, C, D, T)
    return Am, Bm, Cm, Dm, T


# ---------------------------------------------------------------- simulation

def discretize_zoh(A: np.ndarray, B: np.ndarray, dt: float) -> tuple[np.ndarray, np.ndarray]:
    """Return (A_d, B_d) = (e^{A dt}, int_0^dt e^{A s} ds B) from the exponential of [[A, B], [0, 0]] dt."""
    n, m = B.shape
    M = np.zeros((n + m, n + m))
    M[:n, :n] = A
    M[:n, n:] = B
    E = expm(M * dt)
    return E[:n, :n], E[:n, n:]


def simulate(A: np.ndarray, B: np.ndarray, C: np.ndarray, D: np.ndarray, u: np.ndarray, dt: float,
             x0: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Return (y, x): outputs y_k = C x_k + D u_k and states x_k, k = 0..N-1, for u held constant over each dt.

    u has shape (N,) for one input or (N, m). x_{k+1} = A_d x_k + B_d u_k is exact for a zero-order-hold input.
    """
    u = np.asarray(u, dtype=float)
    if u.ndim == 1:
        u = u[:, None]
    n = A.shape[0]
    Ad, Bd = discretize_zoh(A, B, dt)
    x = np.zeros(n) if x0 is None else np.asarray(x0, dtype=float).copy()
    xs = np.empty((len(u), n))
    for k in range(len(u)):
        xs[k] = x
        x = Ad @ x + Bd @ u[k]
    y = xs @ C.T + u @ D.T
    return (y[:, 0] if y.shape[1] == 1 else y), xs


def dc_gain(A: np.ndarray, B: np.ndarray, C: np.ndarray, D: np.ndarray) -> float:
    """Return G(0) = D - C A^{-1} B (A must be invertible: no pole at the origin)."""
    return float((D - C @ np.linalg.solve(A, B))[0, 0])


# ---------------------------------------------------------------- the aircraft model

def pitch_altitude_ss(p: dict[str, float]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return (A, B, C, D) for x = (alpha, q, theta, h), input delta, output h, written from the physical ODEs.

    alpha' = Z_alpha alpha + q;  q' = M_alpha alpha + M_q q + M_delta delta;  theta' = q;
    h' = V (theta - alpha)   (climb rate = airspeed x flight-path angle gamma = theta - alpha, small angles).
    """
    za, ma, mq, md, V = p["Z_alpha"], p["M_alpha"], p["M_q"], p["M_delta"], p["V"]
    A = np.array([[za, 1.0, 0.0, 0.0],
                  [ma, mq, 0.0, 0.0],
                  [0.0, 1.0, 0.0, 0.0],
                  [-V, 0.0, V, 0.0]])
    B = np.array([[0.0], [md], [0.0], [0.0]])
    C = np.array([[0.0, 0.0, 0.0, 1.0]])
    return A, B, C, np.zeros((1, 1))


def pitch_altitude_tf_hand(p: dict[str, float]) -> tuple[np.ndarray, np.ndarray]:
    """Return (num, den) of h/delta derived by hand: -V M_delta Z_alpha / (s^2 (s^2 - (Z_alpha + M_q) s + Z_alpha M_q - M_alpha)).

    From the ODEs: alpha/delta = M_delta/Delta(s), theta/delta = M_delta (s - Z_alpha)/(s Delta), so
    gamma/delta = (theta - alpha)/delta = -M_delta Z_alpha/(s Delta) and h = V gamma / s.
    """
    za, ma, mq, md, V = p["Z_alpha"], p["M_alpha"], p["M_q"], p["M_delta"], p["V"]
    delta = np.array([1.0, -(za + mq), za * mq - ma])
    return np.array([-V * md * za]), np.polymul([1.0, 0.0, 0.0], delta)


# ---------------------------------------------------------------- the example

def main(out_dir: str | Path = ".") -> None:
    """Build the pitch-altitude model three ways, check them against each other and SciPy, save state_space.png."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy import signal

    out_dir = Path(out_dir)
    p = pitch_altitude_params()
    A, B, C, D = pitch_altitude_ss(p)
    num_h, den_h = pitch_altitude_tf_hand(p)
    print("representative small-UAV pitch-altitude model (not a specific aircraft), V = 30 m/s")
    print("eigenvalues of A:", np.array2string(np.linalg.eigvals(A), precision=4, suppress_small=True))
    print(f"hand TF  h/delta = {num_h[0]:.2f} / (s^2 (s^2 + {den_h[1]:.2f} s + {den_h[2]:.2f}))")
    num_f, den_f = ss_to_tf(A, B, C, D)
    print("Faddeev-LeVerrier: den =", np.array2string(den_f, precision=6), " num =", np.array2string(num_f, precision=6, suppress_small=True))
    s_test = np.array([0.5 + 1.0j, 2.0j, -0.3 + 4.0j, 10.0 + 0.0j])
    err = np.max(np.abs(ss_eval(A, B, C, D, s_test) - tf_eval(num_h, den_h, s_test)) / np.abs(tf_eval(num_h, den_h, s_test)))
    print(f"max relative |C(sI-A)^-1 B + D - hand TF| at 4 test points = {err:.2e}")

    Ac, Bc, Cc, Dc = tf_to_controllable(num_h, den_h)
    Ao, Bo, Co, Do = tf_to_observable(num_h, den_h)
    print("companion A first row:", np.array2string(Ac[0], precision=3), " C:", np.array2string(Cc[0], precision=3))

    # modal form of the short-period-plus-theta part (h and theta add poles at the origin: keep the stable part)
    Asp, Bsp, Csp, Dsp = A[:2, :2], B[:2], np.array([[1.0, 0.0]]), np.zeros((1, 1))
    Am, Bm, Cm, Dm, T = modal_form(Asp, Bsp, Csp, Dsp)
    print("short-period modal block:", np.array2string(Am, precision=4, suppress_small=True).replace("\n", ""))

    dt = 0.01
    t = np.arange(0.0, 6.0 + dt / 2, dt)
    u = np.radians(1.0) * np.ones_like(t)          # 1 deg elevator step, trailing edge up
    y_phys, x_phys = simulate(A, B, C, D, u, dt)
    y_comp, _ = simulate(Ac, Bc, Cc, Dc, u, dt)
    y_obs, _ = simulate(Ao, Bo, Co, Do, u, dt)
    _, y_sp = signal.step(signal.lti(A, B, C, D), T=t)
    y_sp = np.radians(1.0) * y_sp
    print(f"1 deg elevator step: h(6 s) = {y_phys[-1]:.3f} m, theta(6 s) = {np.degrees(x_phys[-1, 2]):.3f} deg, alpha(6 s) = {np.degrees(x_phys[-1, 0]):.4f} deg")
    print(f"max |physical - companion| = {np.abs(y_phys - y_comp).max():.2e} m, |physical - observable| = {np.abs(y_phys - y_obs).max():.2e} m,"
          f" |physical - scipy.signal.step| = {np.abs(y_phys - y_sp).max():.2e} m")
    # steady climb: alpha -> alpha_ss, gamma -> gamma_ss, so h grows like V gamma_ss t
    a_ss = dc_gain(Asp, Bsp, Csp, Dsp) * np.radians(1.0)
    q_over = -(p["M_alpha"] * a_ss + p["M_delta"] * np.radians(1.0)) / p["M_q"]
    print(f"steady state of the short period: alpha_ss = {np.degrees(a_ss):.4f} deg, pitch rate q_ss = {np.degrees(q_over):.4f} deg/s"
          f" (from the ODEs with alpha' = q' = 0: q_ss = -Z_alpha alpha_ss)")

    print("\ncoefficient-sensitivity experiment: (s + 1)(s + 2)...(s + n), coefficients rounded to 5 significant figures")
    for n in (4, 6, 8, 10):
        poles = -np.arange(1.0, n + 1.0)
        den = np.real(np.poly(poles))
        lam_c = np.linalg.eigvals(tf_to_controllable(np.array([1.0]), round_sig(den, 5))[0])
        lam_m = np.linalg.eigvals(round_sig(np.diag(poles), 5))
        e_comp = max(np.min(np.abs(lam_c - pk)) for pk in poles)
        e_mod = max(np.min(np.abs(lam_m - pk)) for pk in poles)
        print(f"  n = {n:2d}: largest coefficient {den.max():.3e}; pole error companion {e_comp:.2e} (largest |Im| {np.abs(lam_c.imag).max():.2f}),"
              f" modal {e_mod:.2e}")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    axes[0].plot(t, y_phys, label="physical (alpha, q, theta, h)")
    axes[0].plot(t, y_comp, "--", label="companion form")
    axes[0].plot(t[::25], y_sp[::25], "k.", ms=4, label="scipy.signal.step")
    axes[0].set_xlabel("t [s]"); axes[0].set_ylabel("altitude change h [m]"); axes[0].legend(fontsize=8); axes[0].grid(lw=0.4)
    axes[0].set_title("1 deg elevator step (simulation)", fontsize=9)
    axes[1].plot(t, np.degrees(x_phys[:, 0]), label="alpha")
    axes[1].plot(t, np.degrees(x_phys[:, 2]), label="theta")
    axes[1].plot(t, np.degrees(x_phys[:, 2] - x_phys[:, 0]), label="gamma = theta - alpha")
    axes[1].set_xlabel("t [s]"); axes[1].set_ylabel("angle [deg]"); axes[1].legend(fontsize=8); axes[1].grid(lw=0.4)
    fig.tight_layout()
    fig.savefig(out_dir / "state_space.png", dpi=150)
    plt.close(fig)
    print(f"saved {out_dir / 'state_space.png'}")


if __name__ == "__main__":
    main()
