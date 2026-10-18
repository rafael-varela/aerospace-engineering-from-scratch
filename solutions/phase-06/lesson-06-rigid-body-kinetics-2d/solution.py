"""P06.06 Rigid Body Kinetics in 2D: reference solution.

A landing-gear wheel touches the runway at 70 m/s without turning. For a few hundredths of a second
the tyre slides, Coulomb friction drags the axle back and spins the wheel up, and then the tyre rolls.
This module answers the questions that touchdown raises, with plane rigid-body kinetics
(sum of forces = m a_G, sum of moments about G = I_G alpha) and a rolling constraint v = omega R:

* the warm-up: a wheel on an incline, rolling or sliding, with the friction it needs to roll;
* the closed form of the touchdown spin-up: the friction impulse, the spin-up time, the final speed,
  the slip distance and the energy dissipated by friction, for a constant normal load and for a load
  that ramps up as the strut compresses;
* a hybrid time integration with `scipy.integrate.solve_ivp`: a sliding phase that ends with an event
  when the slip velocity reaches zero, and a rolling phase that ends with an event if braking asks for
  more friction than the tyre can give (the wheel then skids again);
* an energy audit: the energy dissipated by friction computed two ways, as the integral of
  friction force x slip velocity (heat made at the contact) and as the loss of kinetic energy.

Frame: x forward along the runway (inertial, the runway is flat), the axle moves at speed V > 0, the
wheel turns at omega > 0 when rolling forward, so the contact point slips forward at v_s = V - omega R.
The equivalent mass M is the share of the aircraft carried by this wheel; it includes the wheel's own
mass. N is the runway's normal force on the tyre. All values are SI. Numbers are illustrative.

Run:  python solution.py [out_dir]     (prints the lesson's numbers, writes wheel_spinup.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-06/lesson-06-rigid-body-kinetics-2d
"""
from __future__ import annotations

import math
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp

G0 = 9.80665


# ---------------------------------------------------------------- the illustrative main-gear wheel (given)

@dataclass
class Wheel:
    """An illustrative single-aisle airliner main-gear wheel and its share of the aircraft.

    m_wheel   rotating mass (tyre, wheel, brake rotors), kg
    R         rolling radius, m
    k_ratio   radius of gyration / R (0.5 solid disk, 1.0 thin hoop; a tyre on a wheel is about 0.6-0.75)
    M         equivalent translating mass carried by this wheel (aircraft share, includes the wheel), kg
    """
    m_wheel: float = 120.0
    R: float = 0.56
    k_ratio: float = 0.65
    M: float = 15000.0

    @property
    def I(self) -> float:
        """Moment of inertia of the rotating parts about the axle, kg m^2."""
        return self.m_wheel * (self.k_ratio * self.R) ** 2


def ramp_load(N_max: float, t_ramp: float) -> Callable[[float], float]:
    """Normal load that rises linearly from 0 to N_max in t_ramp seconds and then stays (a strut
    compressing at touchdown). t_ramp = 0 gives a constant load."""
    if t_ramp <= 0.0:
        return lambda t: N_max
    return lambda t: N_max * min(max(t, 0.0) / t_ramp, 1.0)


# ---------------------------------------------------------------- warm-up: the incline

def incline_motion(theta: float, m: float, I: float, R: float, mu_s: float, mu_k: float,
                   g: float = G0) -> dict:
    """A round body (mass m, inertia I about its centre, radius R) released on an incline of angle theta.

    Plane kinetics: m a = m g sin(theta) - f, I alpha = f R, N = m g cos(theta).
    If the friction needed to roll, f_roll = m g sin(theta) I / (I + m R^2), is at most mu_s N it rolls:
    a = g sin(theta) / (1 + I/(m R^2)), alpha = a / R. Otherwise it slides with f = mu_k N:
    a = g (sin(theta) - mu_k cos(theta)), alpha = mu_k N R / I.

    Returns {"rolls": bool, "a": m/s^2, "alpha": rad/s^2, "friction": N, "friction_needed": N, "normal": N}.
    """
    N = m * g * math.cos(theta)
    f_need = m * g * math.sin(theta) * I / (I + m * R * R)
    if f_need <= mu_s * N:
        a = g * math.sin(theta) / (1.0 + I / (m * R * R))
        return {"rolls": True, "a": a, "alpha": a / R, "friction": f_need, "friction_needed": f_need, "normal": N}
    f = mu_k * N
    return {"rolls": False, "a": g * math.sin(theta) - f / m, "alpha": f * R / I, "friction": f,
            "friction_needed": f_need, "normal": N}


def min_friction_to_roll(theta: float, m: float, I: float, R: float) -> float:
    """Smallest static friction coefficient that lets the body roll down the incline:
    mu = tan(theta) / (1 + m R^2 / I)."""
    return math.tan(theta) / (1.0 + m * R * R / I)


# ---------------------------------------------------------------- touchdown spin-up in closed form

def effective_mass(M: float, I: float, R: float) -> float:
    """The mass that 'feels' the slip: 1/M_eff = 1/M + R^2/I. A friction force F changes the slip
    velocity v_s = V - omega R at the rate -F/M_eff."""
    return 1.0 / (1.0 / M + R * R / I)


def spin_up_closed_form(V0: float, R: float, I: float, M: float, mu_k: float, N_max: float,
                        t_ramp: float = 0.0, omega0: float = 0.0) -> dict:
    """Closed-form touchdown spin-up with Coulomb friction F = mu_k N(t), N ramping linearly from 0 to
    N_max in t_ramp seconds (t_ramp = 0: constant load). Requires V0 > omega0 R (forward slip).

    Impulse-momentum: M (V0 - Vf) = J, I (omega_f - omega0) = J R, rolling at the end: Vf = omega_f R
      => J = M_eff v_s0 with v_s0 = V0 - omega0 R.
    Spin-up time from J = mu_k * integral of N dt.
    Energy dissipated by friction: E_f = M_eff v_s0^2 / 2 (independent of mu_k and of N(t)).
    Slip distance s = integral of v_s dt, with v_s(t) = v_s0 - (mu_k / M_eff) Q(t), Q = integral of N.

    Returns {"impulse", "t_spin", "V_final", "omega_final", "energy_friction", "slip_distance",
    "peak_drag"} in SI units.
    """
    vs0 = V0 - omega0 * R
    if vs0 <= 0.0:
        raise ValueError("the closed form needs forward slip, V0 > omega0 R")
    Me = effective_mass(M, I, R)
    J = Me * vs0
    if t_ramp <= 0.0:
        ts = J / (mu_k * N_max)
        s = 0.5 * vs0 * ts
        peak = mu_k * N_max
    elif J <= 0.5 * mu_k * N_max * t_ramp:                     # spin-up ends while the load is still rising
        ts = math.sqrt(2.0 * J * t_ramp / (mu_k * N_max))
        # v_s = vs0 - c t^2 with c = mu_k N_max / (2 Me t_ramp); integral = vs0 ts - c ts^3 / 3
        c = mu_k * N_max / (2.0 * Me * t_ramp)
        s = vs0 * ts - c * ts ** 3 / 3.0
        peak = mu_k * N_max * ts / t_ramp
    else:
        ts = J / (mu_k * N_max) + 0.5 * t_ramp
        c = mu_k * N_max / (2.0 * Me * t_ramp)
        s1 = vs0 * t_ramp - c * t_ramp ** 3 / 3.0             # slip distance during the ramp
        v1 = vs0 - c * t_ramp ** 2                            # slip velocity at the end of the ramp
        s = s1 + 0.5 * v1 * (ts - t_ramp)                     # then it falls linearly to zero
        peak = mu_k * N_max
    return {"impulse": J, "t_spin": ts, "V_final": V0 - J / M, "omega_final": omega0 + J * R / I,
            "energy_friction": 0.5 * Me * vs0 ** 2, "slip_distance": s, "peak_drag": peak}


def kinetic_energy_loss(V0: float, R: float, I: float, M: float, omega0: float = 0.0) -> float:
    """The second route to the friction energy: kinetic energy before minus after, with the final state
    from impulse-momentum (no friction coefficient, no load history)."""
    Me = effective_mass(M, I, R)
    J = Me * (V0 - omega0 * R)
    Vf, wf = V0 - J / M, omega0 + J * R / I
    return 0.5 * M * V0 ** 2 + 0.5 * I * omega0 ** 2 - (0.5 * M * Vf ** 2 + 0.5 * I * wf ** 2)


def rolling_friction_for_brake(T_brake: float, I: float, M: float, R: float) -> float:
    """Friction force the runway must supply to keep a braked wheel rolling (no skid):
    M dV/dt = -F, I domega/dt = F R - T_b, dV/dt = R domega/dt  =>  F = T_b M R / (I + M R^2)."""
    return T_brake * M * R / (I + M * R * R)


def max_brake_torque(mu_s: float, N: float, I: float, M: float, R: float) -> float:
    """Largest brake torque that keeps the wheel rolling: F = mu_s N in rolling_friction_for_brake."""
    return mu_s * N * (I + M * R * R) / (M * R)


# ---------------------------------------------------------------- hybrid time integration

@dataclass
class Touchdown:
    """Result of `simulate_touchdown`. Arrays are sampled on t (s).

    V (m/s), omega (rad/s), slip = V - omega R (m/s), friction (N, the runway's force on the tyre, positive
    backward, i.e. a drag on the axle), normal (N), E_friction (J, integral of |friction x slip velocity|),
    E_brake (J, integral of T_brake x omega), slip_distance (m), x (m, axle position), phase (one of
    "sliding", "rolling", "locked" per sample), rolling (bool per sample), t_spin (s, end of the first
    sliding phase, nan if it never ended), transitions (list of (t, phase) at each phase change).
    """
    t: np.ndarray
    x: np.ndarray
    V: np.ndarray
    omega: np.ndarray
    slip: np.ndarray
    friction: np.ndarray
    normal: np.ndarray
    E_friction: np.ndarray
    E_brake: np.ndarray
    slip_distance: np.ndarray
    phase: np.ndarray
    rolling: np.ndarray
    t_spin: float
    transitions: list = field(default_factory=list)


def simulate_touchdown(wheel: Wheel, V0: float, mu_k: float, mu_s: float, normal: Callable[[float], float],
                       t_end: float, omega0: float = 0.0, brake_torque: Callable[[float], float] | None = None,
                       n_samples: int = 2001, rtol: float = 1e-10, atol: float = 1e-10) -> Touchdown:
    """Integrate the touchdown as a hybrid system with three phases, each ended by a solve_ivp event.

    States y = [x, V, omega, E_friction, E_brake, slip_distance]; the brake torque T_b(t) >= 0 opposes
    the wheel's rotation (omega >= 0 throughout).
    sliding (v_s = V - omega R != 0): F = mu_k N sign(v_s); M dV/dt = -F; I domega/dt = F R - T_b.
        Ends when v_s reaches zero (-> rolling) or omega reaches zero under braking (-> locked).
    rolling (v_s = 0): the friction that holds V = omega R is F = T_b M R / (I + M R^2).
        Ends when that F exceeds mu_s N (-> sliding: the tyre skids, v_s > 0 under braking).
    locked (omega = 0): F = mu_k N, the brake holds the wheel while T_b >= mu_k N R.
        Ends when the friction torque beats the brake (-> sliding).
    Friction heats the contact at |F v_s| (zero while rolling); the brake absorbs T_b omega.
    """
    I, R, M = wheel.I, wheel.R, wheel.M
    Tb = brake_torque if brake_torque is not None else (lambda t: 0.0)
    t_eval = np.linspace(0.0, t_end, n_samples)
    y = np.array([0.0, V0, omega0, 0.0, 0.0, 0.0])
    t0 = 0.0
    vs0 = V0 - omega0 * R
    phase = "rolling" if abs(vs0) <= 1e-12 * max(1.0, V0) else "sliding"
    sgn = 1.0 if vs0 >= 0 else -1.0
    out_t, out_y, out_ph, out_F = [np.array([0.0])], [y[:, None].copy()], [np.array([phase])], []
    transitions = [(0.0, phase)]
    t_spin = math.nan

    def friction(ph: str, t: float, s: float) -> float:
        if ph == "rolling":
            return rolling_friction_for_brake(Tb(t), I, M, R)
        if ph == "locked":
            return mu_k * normal(t)
        return mu_k * normal(t) * s

    out_F.append(np.array([friction(phase, 0.0, sgn)]))
    while t0 < t_end:
        s = sgn
        if phase == "sliding":
            def rhs(t, y):
                F = mu_k * normal(t) * s
                vs = y[1] - y[2] * R
                return [y[1], -F / M, (F * R - Tb(t)) / I, abs(F * vs), Tb(t) * y[2], abs(vs)]

            def ev_roll(t, y):
                return y[1] - y[2] * R
            ev_roll.terminal, ev_roll.direction = True, -s

            def ev_lock(t, y):
                return y[2]
            ev_lock.terminal, ev_lock.direction = True, -1
            events = [ev_roll, ev_lock]
        elif phase == "rolling":
            def rhs(t, y):
                F = rolling_friction_for_brake(Tb(t), I, M, R)
                return [y[1], -F / M, -F / (M * R), 0.0, Tb(t) * y[2], 0.0]

            def ev_skid(t, y):
                return mu_s * normal(t) - abs(rolling_friction_for_brake(Tb(t), I, M, R))
            ev_skid.terminal, ev_skid.direction = True, -1
            events = [ev_skid]
        else:                                                       # locked
            def rhs(t, y):
                F = mu_k * normal(t)
                return [y[1], -F / M, 0.0, abs(F * y[1]), 0.0, abs(y[1])]

            def ev_free(t, y):
                return Tb(t) - mu_k * normal(t) * R
            ev_free.terminal, ev_free.direction = True, -1
            events = [ev_free]
        te = t_eval[(t_eval > t0) & (t_eval <= t_end)]
        sol = solve_ivp(rhs, (t0, t_end), y, method="DOP853", t_eval=te, events=events, rtol=rtol, atol=atol,
                        max_step=t_end / 200.0)
        if sol.t.size:
            out_t.append(sol.t)
            out_y.append(sol.y)
            out_ph.append(np.full(sol.t.size, phase))
            out_F.append(np.array([friction(phase, t, s) for t in sol.t]))
        hit = [k for k, te_k in enumerate(sol.t_events) if te_k.size]
        if sol.status != 1 or not hit:
            break
        k = min(hit, key=lambda j: sol.t_events[j][0])
        t0 = float(sol.t_events[k][0])
        y = sol.y_events[k][0].copy()
        if phase == "sliding" and k == 0:
            y[1] = y[2] * R                                         # land exactly on the rolling constraint
            if math.isnan(t_spin):
                t_spin = t0
            phase = "rolling"
        elif phase == "sliding":
            y[2] = 0.0
            phase = "locked"
        elif phase == "rolling":
            sgn = 1.0 if rolling_friction_for_brake(Tb(t0), I, M, R) > 0 else -1.0
            y[2] -= sgn * 1e-12 / R                                 # step off v_s = 0 in the skid direction
            phase = "sliding"
        else:
            sgn = 1.0
            y[2] = 1e-12
            phase = "sliding"
        transitions.append((t0, phase))
    ts = np.concatenate(out_t)
    ys = np.concatenate(out_y, axis=1)
    ph = np.concatenate(out_ph)
    Nn = np.array([normal(t) for t in ts])
    return Touchdown(t=ts, x=ys[0], V=ys[1], omega=ys[2], slip=ys[1] - ys[2] * R, friction=np.concatenate(out_F),
                     normal=Nn, E_friction=ys[3], E_brake=ys[4], slip_distance=ys[5], phase=ph,
                     rolling=ph == "rolling", t_spin=t_spin, transitions=transitions)


def energy_audit(res: Touchdown, wheel: Wheel) -> dict:
    """Kinetic energy, the two dissipations and the residual KE(0) - KE(t) - E_friction - E_brake,
    which must be zero (to the solver tolerance) at every sample."""
    KE = 0.5 * wheel.M * res.V ** 2 + 0.5 * wheel.I * res.omega ** 2
    resid = KE[0] - KE - res.E_friction - res.E_brake
    return {"KE": KE, "residual": resid, "max_rel_residual": float(np.max(np.abs(resid)) / KE[0])}


# ---------------------------------------------------------------- the lesson's output

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save wheel_spinup.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    wh = Wheel()
    I, R, M = wh.I, wh.R, wh.M
    V0, mu_k, mu_s, N_max, t_r = 70.0, 0.55, 0.8, 120e3, 0.10
    print(f"wheel: {wh.m_wheel:.0f} kg rotating, R = {R} m, k = {wh.k_ratio} R, I = {I:.3f} kg m^2; equivalent mass M = {M:.0f} kg; "
          f"M_eff = {effective_mass(M, I, R):.3f} kg")

    print("incline warm-up, 10 deg, mu_s 0.3, mu_k 0.25, R 0.3 m, 10 kg:")
    for name, k in [("hoop", 1.0), ("disk", math.sqrt(0.5)), ("sphere", math.sqrt(0.4))]:
        Ib = 10.0 * (k * 0.3) ** 2
        r = incline_motion(math.radians(10.0), 10.0, Ib, 0.3, 0.3, 0.25)
        print(f"  {name:6s}: rolls {r['rolls']}, a = {r['a']:.4f} m/s^2, friction {r['friction']:.3f} N, "
              f"mu needed {min_friction_to_roll(math.radians(10.0), 10.0, Ib, 0.3):.4f}")

    cf_c = spin_up_closed_form(V0, R, I, M, mu_k, N_max, 0.0)
    cf_r = spin_up_closed_form(V0, R, I, M, mu_k, N_max, t_r)
    dKE = kinetic_energy_loss(V0, R, I, M)
    print(f"touchdown V0 = {V0} m/s, mu_k = {mu_k}, N_max = {N_max / 1e3:.0f} kN")
    for lab, cf in [("constant load", cf_c), (f"load ramp {t_r} s", cf_r)]:
        print(f"  {lab:16s}: impulse {cf['impulse']:.1f} N s, t_spin {cf['t_spin'] * 1e3:.2f} ms, V_f {cf['V_final']:.4f} m/s, "
              f"omega_f {cf['omega_final']:.3f} rad/s ({cf['omega_final'] * 60 / (2 * math.pi):.0f} rpm), slip distance {cf['slip_distance']:.3f} m, "
              f"E_f {cf['energy_friction'] / 1e3:.3f} kJ, peak drag {cf['peak_drag'] / 1e3:.1f} kN")
    print(f"  kinetic-energy loss route: {dKE / 1e3:.3f} kJ; wheel's final KE 0.5 I omega_f^2 = {0.5 * I * cf_c['omega_final'] ** 2 / 1e3:.3f} kJ")

    res = simulate_touchdown(wh, V0, mu_k, mu_s, ramp_load(N_max, t_r), 0.3)
    au = energy_audit(res, wh)
    print(f"  simulated (ramp): t_spin {res.t_spin * 1e3:.4f} ms, E_friction {res.E_friction[-1] / 1e3:.4f} kJ, "
          f"KE loss {(au['KE'][0] - au['KE'][-1]) / 1e3:.4f} kJ, slip distance {res.slip_distance[-1]:.4f} m, "
          f"max relative energy residual {au['max_rel_residual']:.1e}")
    print(f"  differences sim - closed form: t_spin {res.t_spin - cf_r['t_spin']:.1e} s, E_f {res.E_friction[-1] - cf_r['energy_friction']:.1e} J")

    print("energy and spin-up time against mu_k (ramp load):")
    for mk in [0.3, 0.55, 0.8]:
        cf = spin_up_closed_form(V0, R, I, M, mk, N_max, t_r)
        print(f"  mu_k {mk:.2f}: t_spin {cf['t_spin'] * 1e3:.2f} ms, slip {cf['slip_distance']:.3f} m, E_f {cf['energy_friction'] / 1e3:.3f} kJ")
    print("pre-spin before touchdown (fraction of V0/R):")
    for f in [0.0, 0.5, 0.9]:
        cf = spin_up_closed_form(V0, R, I, M, mu_k, N_max, t_r, omega0=f * V0 / R)
        print(f"  {f:.1f}: E_f {cf['energy_friction'] / 1e3:.3f} kJ, t_spin {cf['t_spin'] * 1e3:.2f} ms")
    Tmax = max_brake_torque(mu_s, N_max, I, M, R)
    print(f"max brake torque without skid at N = {N_max / 1e3:.0f} kN, mu_s = {mu_s}: {Tmax / 1e3:.2f} kN m (mu_s N R = {mu_s * N_max * R / 1e3:.2f} kN m)")
    Tb_over = 1.2 * Tmax
    resb = simulate_touchdown(wh, V0, mu_k, mu_s, ramp_load(N_max, t_r), 0.6,
                              brake_torque=lambda t: Tb_over if t >= 0.3 else 0.0)
    aub = energy_audit(resb, wh)
    print(f"brake 1.2 x max applied at 0.3 s: transitions {[(round(t, 4), p) for t, p in resb.transitions]}; "
          f"final V {resb.V[-1]:.3f} m/s, omega {resb.omega[-1]:.2f} rad/s; friction heat {resb.E_friction[-1] / 1e3:.1f} kJ, "
          f"brake {resb.E_brake[-1] / 1e3:.1f} kJ; energy residual {aub['max_rel_residual']:.1e}")
    resg = simulate_touchdown(wh, V0, mu_k, mu_s, ramp_load(N_max, t_r), 0.6,
                              brake_torque=lambda t: 0.9 * Tmax if t >= 0.3 else 0.0)
    aug = energy_audit(resg, wh)
    print(f"brake 0.9 x max applied at 0.3 s: transitions {[(round(t, 4), p) for t, p in resg.transitions]}; "
          f"deceleration {(resg.V[-2] - resg.V[-1]) / (resg.t[-1] - resg.t[-2]):.3f} m/s^2; brake {resg.E_brake[-1] / 1e3:.1f} kJ; "
          f"energy residual {aug['max_rel_residual']:.1e}")

    fig, ax = plt.subplots(1, 3, figsize=(14, 4))
    tm = res.t * 1e3
    ax[0].plot(tm, res.V, label="axle V")
    ax[0].plot(tm, res.omega * R, label="wheel rim speed omega R")
    ax[0].axvline(res.t_spin * 1e3, color="k", ls=":", lw=0.8)
    ax[0].set_xlabel("t (ms)"); ax[0].set_ylabel("m/s"); ax[0].set_xlim(0, 150); ax[0].legend(); ax[0].set_title("spin-up (simulation)")
    ax[1].plot(tm, res.friction / 1e3, label="friction (drag on axle)")
    ax[1].plot(tm, res.normal / 1e3 * mu_k, "--", label="mu_k N")
    ax[1].set_xlabel("t (ms)"); ax[1].set_ylabel("kN"); ax[1].set_xlim(0, 150); ax[1].legend(); ax[1].set_title("spin-up drag load")
    ax[2].plot(tm, res.E_friction / 1e3, label="integral F v_s dt")
    ax[2].plot(tm, (au["KE"][0] - au["KE"]) / 1e3, "--", label="KE(0) - KE(t)")
    ax[2].axhline(cf_r["energy_friction"] / 1e3, color="k", ls=":", lw=0.8, label="closed form")
    ax[2].set_xlabel("t (ms)"); ax[2].set_ylabel("kJ"); ax[2].set_xlim(0, 150); ax[2].legend(); ax[2].set_title("friction energy, two routes")
    fig.tight_layout()
    fig.savefig(out_dir / "wheel_spinup.png", dpi=120)
    plt.close(fig)
    print(f"wrote {out_dir / 'wheel_spinup.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
