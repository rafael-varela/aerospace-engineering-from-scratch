"""P04.08 Electrostatics & Circuits: starter code.

Complete every function marked NotImplementedError. When finished, running this file produces the
lesson's figure and prints the verification tables.

From Coulomb's law to a DC circuit solver, applied to a satellite power bus:

* Coulomb force, the electric field and potential of point charges, superposition, and E = -grad V;
* resistance from resistivity, Ohm's law, power; series and parallel combinations; the voltage divider;
* Kirchhoff's rules written as modified nodal analysis (MNA): one KCL equation per node plus one
  constraint per voltage source, solved as a linear system (P01.13) for every circuit at once,
  including the two-loop Wheatstone bridge that series-parallel reduction cannot handle;
* meter loading: a voltmeter of finite resistance in parallel, an ammeter of finite resistance in series;
* the RC transient: closed form against solve_ivp, the time constant, and the energy accounting
  (charging a capacitor through any resistor wastes exactly half the energy drawn from the source);
* a 28 V satellite power bus with a copper harness: voltage at the loads, wiring loss, wire sizing.

The optional measurement with the p04-electronics-bench kit (P04.15) is described in the lesson; nothing
here depends on it.

Units: SI (V, A, ohm, F, s, m). Node 0 is ground.

Run:  python main.py          (prints the verification tables, writes dc_circuits.png)
Test: pytest curriculum/phase-04/lesson-08-electrostatics-circuits
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

K_COULOMB = 8.9875517923e9           # N m^2 / C^2 (1 / 4 pi eps0)
EPS0 = 8.8541878128e-12              # F/m
E_CHARGE = 1.602176634e-19           # C
RHO_COPPER = 1.68e-8                 # ohm m at 20 C
RHO_ALUMINIUM = 2.65e-8              # ohm m at 20 C


# ---------------------------------------------------------------- electrostatics

def coulomb_force(q1: float, q2: float, r1: np.ndarray, r2: np.ndarray) -> np.ndarray:
    """Force on charge q2 at r2 due to q1 at r1: k q1 q2 (r2 - r1) / |r2 - r1|^3 (repulsive for like charges)."""
    raise NotImplementedError


def electric_field(charges: np.ndarray, positions: np.ndarray, r: np.ndarray) -> np.ndarray:
    """Field at r from point charges (superposition): sum_i k q_i (r - r_i) / |r - r_i|^3.
    ``positions`` has shape (n, d); ``r`` is one point of shape (d,)."""
    raise NotImplementedError


def electric_potential(charges: np.ndarray, positions: np.ndarray, r: np.ndarray) -> float:
    """Potential at r (zero at infinity): sum_i k q_i / |r - r_i|."""
    raise NotImplementedError


def field_from_potential(charges: np.ndarray, positions: np.ndarray, r: np.ndarray, h: float = 1e-6) -> np.ndarray:
    """E = -grad V by central differences of ``electric_potential``, a check on ``electric_field``."""
    raise NotImplementedError


def dipole_moment(charges: np.ndarray, positions: np.ndarray) -> np.ndarray:
    """p = sum_i q_i r_i (independent of origin when the total charge is zero)."""
    raise NotImplementedError


# ---------------------------------------------------------------- resistance, Ohm, power

def resistance(rho: float, length: float, area: float) -> float:
    """R = rho L / A."""
    raise NotImplementedError


def wire_resistance(length: float, diameter: float, rho: float = RHO_COPPER) -> float:
    """Resistance of a round wire of the given length and diameter."""
    raise NotImplementedError


def power_dissipated(V: float, I: float) -> float:
    """P = V I (= I^2 R = V^2 / R for a resistor)."""
    raise NotImplementedError


def series(resistances: np.ndarray | list) -> float:
    """R = sum R_i."""
    raise NotImplementedError


def parallel(resistances: np.ndarray | list) -> float:
    """1 / R = sum 1 / R_i."""
    raise NotImplementedError


def voltage_divider(Vs: float, R1: float, R2: float) -> float:
    """Voltage across R2 when R1 and R2 are in series across Vs: Vs R2 / (R1 + R2)."""
    raise NotImplementedError


# ---------------------------------------------------------------- the circuit solver (modified nodal analysis)

@dataclass
class Resistor:
    n1: int
    n2: int
    R: float
    name: str = ""


@dataclass
class VoltageSource:
    """Ideal source: node n_plus sits V volts above node n_minus."""
    n_plus: int
    n_minus: int
    V: float
    name: str = ""


@dataclass
class CurrentSource:
    """Ideal source: current I flows through the source from n_from to n_to (it is injected into n_to)."""
    n_from: int
    n_to: int
    I: float
    name: str = ""


@dataclass
class Circuit:
    """A list of components between integer nodes; node 0 is the reference (ground)."""
    components: list = field(default_factory=list)

    def add(self, component) -> "Circuit":
        self.components.append(component)
        return self

    @property
    def n_nodes(self) -> int:
        nodes = {0}
        for c in self.components:
            if isinstance(c, Resistor):
                nodes.update((c.n1, c.n2))
            elif isinstance(c, VoltageSource):
                nodes.update((c.n_plus, c.n_minus))
            elif isinstance(c, CurrentSource):
                nodes.update((c.n_from, c.n_to))
        return max(nodes) + 1

    @property
    def resistors(self) -> list:
        return [c for c in self.components if isinstance(c, Resistor)]

    @property
    def voltage_sources(self) -> list:
        return [c for c in self.components if isinstance(c, VoltageSource)]

    @property
    def current_sources(self) -> list:
        return [c for c in self.components if isinstance(c, CurrentSource)]


@dataclass
class DCSolution:
    node_voltages: np.ndarray            # index = node number; node_voltages[0] = 0
    source_currents: np.ndarray          # current through each voltage source, from n_minus to n_plus inside the source


def build_mna(circuit: Circuit) -> tuple[np.ndarray, np.ndarray]:
    """Assemble the MNA system A x = z with unknowns x = [V_1 .. V_{N-1}, I_src_1 .. I_src_M].

    Rows 1..N-1: KCL at each non-ground node, written as sum of conductances times voltages equals
    the injected source current. Rows N..N+M-1: V_plus - V_minus = V for each voltage source, with
    the source current entering the KCL rows through the same incidence pattern (the matrix is
    symmetric for a resistor-and-source network).
    """
    raise NotImplementedError


def solve_dc(circuit: Circuit) -> DCSolution:
    """Solve the MNA system. Raises ValueError for a singular system (a floating node, two ideal sources in parallel)."""
    raise NotImplementedError


def branch_current(r: Resistor, V: np.ndarray) -> float:
    """Current through a resistor from n1 to n2: (V[n1] - V[n2]) / R."""
    raise NotImplementedError


def kcl_residuals(circuit: Circuit, sol: DCSolution) -> np.ndarray:
    """Net current leaving each node (resistors, current sources and voltage sources included); zero at a solution."""
    raise NotImplementedError


def power_balance(circuit: Circuit, sol: DCSolution) -> tuple[float, float]:
    """(power delivered by all sources, power dissipated in all resistors); equal at a solution."""
    raise NotImplementedError


def wheatstone_bridge(Vs: float, R1: float, R2: float, R3: float, Rx: float, Rg: float) -> tuple[Circuit, float]:
    """Bridge: Vs between node 1 and ground; R1 from 1 to 2, R2 from 2 to 0; R3 from 1 to 3, Rx from 3 to 0;
    the galvanometer Rg between 2 and 3. Returns the circuit and the galvanometer current from 2 to 3."""
    raise NotImplementedError


def voltmeter_reading(circuit: Circuit, n1: int, n2: int, R_meter: float) -> float:
    """Voltage V[n1] - V[n2] read by a voltmeter of input resistance R_meter connected across the two nodes."""
    raise NotImplementedError


def ammeter_reading(circuit: Circuit, resistor_index: int, R_meter: float) -> float:
    """Current read by an ammeter of resistance R_meter inserted in series with the resistor at
    ``resistor_index`` (counted among the circuit's resistors)."""
    raise NotImplementedError


# ---------------------------------------------------------------- the RC transient

def rc_time_constant(R: float, C: float) -> float:
    """tau = R C [s]."""
    raise NotImplementedError


def rc_rhs(t: float, y: np.ndarray, Vs: float, R: float, C: float) -> np.ndarray:
    """dV_C/dt = (Vs - V_C) / (R C) for a capacitor charged through R from a source Vs (Vs = 0: discharge)."""
    raise NotImplementedError


def rc_capacitor_voltage(t: np.ndarray | float, Vs: float, R: float, C: float, V0: float = 0.0) -> np.ndarray | float:
    """V_C(t) = Vs + (V0 - Vs) exp(-t / RC)."""
    raise NotImplementedError


def rc_current(t: np.ndarray | float, Vs: float, R: float, C: float, V0: float = 0.0) -> np.ndarray | float:
    """I(t) = (Vs - V_C) / R = (Vs - V0) / R exp(-t / RC)."""
    raise NotImplementedError


def rc_energy_accounting(Vs: float, R: float, C: float, t_end: float, n: int = 20001) -> tuple[float, float, float]:
    """(energy drawn from the source, energy stored in C, energy dissipated in R) over [0, t_end] when charging
    from V0 = 0, by numerical integration of Vs I and I^2 R. As t_end -> inf: C Vs^2, C Vs^2 / 2, C Vs^2 / 2."""
    raise NotImplementedError


# ---------------------------------------------------------------- the satellite bus

def bus_circuit(V_bus: float, R_wire: float, loads: list[float]) -> Circuit:
    """A regulated bus V_bus at node 1, a supply wire R_wire to node 2, the loads (resistances) in parallel
    between node 2 and node 3, and a return wire R_wire from node 3 to ground."""
    raise NotImplementedError


def load_resistance_for_power(P: float, V: float) -> float:
    """Resistance that draws P watts at V volts: V^2 / P."""
    raise NotImplementedError


def wire_diameter_for_drop(I: float, length: float, V_bus: float, max_drop_fraction: float,
                           rho: float = RHO_COPPER) -> float:
    """Smallest round-wire diameter such that the drop over the supply-plus-return run (2 x length) at current I
    is at most max_drop_fraction of V_bus."""
    raise NotImplementedError


# ---------------------------------------------------------------- driver

def main(out_dir: str | Path = ".") -> None:
    """Print the verification tables and write dc_circuits.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)

    # electrostatics
    q = np.array([1e-9, -1e-9])
    pos = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1e-3]])
    r = np.array([0.0, 0.0, 0.05])
    E = electric_field(q, pos, r)
    Efd = field_from_potential(q, pos, r)
    p = dipole_moment(q, pos)
    print(f"dipole p = {np.linalg.norm(p):.2e} C m: E on axis at 50 mm = {E[2]:.4e} V/m, -grad V = {Efd[2]:.4e} V/m, "
          f"far field 2kp/z^3 = {2 * K_COULOMB * np.linalg.norm(p) / 0.05 ** 3:.4e} V/m")
    F = coulomb_force(E_CHARGE, -E_CHARGE, np.zeros(3), np.array([5.29e-11, 0.0, 0.0]))
    print(f"electron-proton force at the Bohr radius: {np.linalg.norm(F):.3e} N (attractive: {F[0] < 0})")

    # resistance and the divider
    d16 = 1.291e-3
    Rw = wire_resistance(10.0, d16)
    print(f"\n10 m of 1.291 mm (AWG 16) copper: R = {Rw:.4f} ohm; same length of 0.5 mm wire: {wire_resistance(10.0, 0.5e-3):.3f} ohm")
    print(f"divider 12 V, 10 k / 4.7 k: V2 = {voltage_divider(12.0, 10e3, 4.7e3):.4f} V; parallel(10k, 4.7k) = {parallel([10e3, 4.7e3]):.1f} ohm")

    # Kirchhoff: the two-loop bridge
    c, ig = wheatstone_bridge(5.0, 1000.0, 2000.0, 1000.0, 2200.0, 500.0)
    sol = solve_dc(c)
    print(f"\nWheatstone bridge 5 V, 1k/2k vs 1k/2.2k, galvanometer 500 ohm: I_g = {ig*1e6:.2f} uA; node voltages {np.round(sol.node_voltages, 5)} V")
    print(f"  KCL residuals max |.| = {np.max(np.abs(kcl_residuals(c, sol))):.1e} A; power delivered {power_balance(c, sol)[0]*1e3:.4f} mW, dissipated {power_balance(c, sol)[1]*1e3:.4f} mW")
    _, ig0 = wheatstone_bridge(5.0, 1000.0, 2000.0, 1000.0, 2000.0, 500.0)
    print(f"  balanced bridge (Rx = 2k): I_g = {ig0:.1e} A")

    # meter loading
    div = Circuit().add(VoltageSource(1, 0, 12.0)).add(Resistor(1, 2, 100e3)).add(Resistor(2, 0, 100e3))
    print(f"\n12 V divider of two 100 k: true V2 = {voltmeter_reading(div, 2, 0, np.inf):.4f} V; a 10 M voltmeter reads "
          f"{voltmeter_reading(div, 2, 0, 10e6):.4f} V, a 1 M meter {voltmeter_reading(div, 2, 0, 1e6):.4f} V, a 100 k meter {voltmeter_reading(div, 2, 0, 100e3):.4f} V")
    led = Circuit().add(VoltageSource(1, 0, 5.0)).add(Resistor(1, 0, 220.0))
    print(f"5 V through 220 ohm: {ammeter_reading(led, 0, 0.0)*1e3:.3f} mA; with a 5 ohm ammeter {ammeter_reading(led, 0, 5.0)*1e3:.3f} mA "
          f"(burden {5.0 - 5.0 * 220 / 225:.3f} V)")

    # RC
    R, C, Vs = 10e3, 100e-6, 5.0
    tau = rc_time_constant(R, C)
    t = np.linspace(0.0, 5 * tau, 501)
    ivp = solve_ivp(rc_rhs, (0.0, 5 * tau), [0.0], args=(Vs, R, C), t_eval=t, method="DOP853", rtol=1e-11, atol=1e-13)
    print(f"\nRC: R = {R/1e3:.0f} k, C = {C*1e6:.0f} uF, tau = {tau:.2f} s; V_C(tau) = {rc_capacitor_voltage(tau, Vs, R, C)/Vs:.4f} Vs, "
          f"V_C(5 tau) = {rc_capacitor_voltage(5*tau, Vs, R, C)/Vs:.4f} Vs; closed form vs DOP853 max |dV| = {np.max(np.abs(rc_capacitor_voltage(t, Vs, R, C) - ivp.y[0])):.1e} V")
    es, est, ed = rc_energy_accounting(Vs, R, C, 20 * tau)
    print(f"  energy to 20 tau: from source {es*1e3:.4f} mJ, stored {est*1e3:.4f} mJ, dissipated {ed*1e3:.4f} mJ (C Vs^2 = {C*Vs**2*1e3:.4f} mJ)")
    es2, est2, ed2 = rc_energy_accounting(Vs, 10 * R, C, 200 * tau)
    print(f"  with 10 R: from source {es2*1e3:.4f} mJ, dissipated {ed2*1e3:.4f} mJ: the lost half does not depend on R")

    # the bus
    V_bus, length = 28.0, 8.0
    loads = [load_resistance_for_power(P, V_bus) for P in (120.0, 60.0, 40.0)]
    print(f"\n28 V bus, loads 120 + 60 + 40 W at nominal voltage over an {length:.0f} m harness (each way):")
    print(f"{'wire d (mm)':>12} {'R_wire (ohm)':>13} {'V_load (V)':>11} {'I (A)':>7} {'P_loads (W)':>12} {'harness loss (W)':>17} {'drop':>7}")
    for d in (0.5e-3, 0.8e-3, 1.291e-3, 2.053e-3):
        Rw = wire_resistance(length, d)
        cb = bus_circuit(V_bus, Rw, loads)
        s = solve_dc(cb)
        Vl = s.node_voltages[2] - s.node_voltages[3]
        I = s.source_currents[0]
        pl = sum(Vl ** 2 / RL for RL in loads)
        print(f"{d*1e3:>12.3f} {Rw:>13.4f} {Vl:>11.3f} {I:>7.3f} {pl:>12.2f} {2 * I**2 * Rw:>17.3f} {(V_bus - Vl)/V_bus:>7.2%}")
    I_nom = 220.0 / V_bus
    d_req = wire_diameter_for_drop(I_nom, length, V_bus, 0.02)
    print(f"  diameter for a 2 % drop at the nominal {I_nom:.2f} A: {d_req*1e3:.2f} mm (copper); aluminium: {wire_diameter_for_drop(I_nom, length, V_bus, 0.02, RHO_ALUMINIUM)*1e3:.2f} mm")

    # figure
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.2))
    ax = axs[0]
    ax.plot(t, rc_capacitor_voltage(t, Vs, R, C), label="closed form")
    ax.plot(t, ivp.y[0], "k--", lw=0.8, label="solve_ivp")
    ax.axvline(tau, color="0.6", ls=":")
    ax.axhline(Vs * (1 - np.exp(-1)), color="0.6", ls=":")
    ax.set_xlabel("t (s)")
    ax.set_ylabel("V_C (V)")
    ax.set_title("RC charging, tau = 1 s")
    ax.legend(fontsize=8)
    ax = axs[1]
    Rm = np.logspace(4, 8, 200)
    ax.semilogx(Rm, [voltmeter_reading(div, 2, 0, rm) for rm in Rm])
    ax.axhline(6.0, color="0.6", ls=":")
    ax.set_xlabel("voltmeter resistance (ohm)")
    ax.set_ylabel("reading (V)")
    ax.set_title("meter loading on a 100 k / 100 k divider")
    ax = axs[2]
    ds = np.linspace(0.4e-3, 2.5e-3, 100)
    drops = []
    for d in ds:
        s = solve_dc(bus_circuit(V_bus, wire_resistance(length, d), loads))
        drops.append(100 * (V_bus - (s.node_voltages[2] - s.node_voltages[3])) / V_bus)
    ax.plot(ds * 1e3, drops)
    ax.axhline(2.0, color="0.6", ls=":")
    ax.set_xlabel("wire diameter (mm)")
    ax.set_ylabel("bus-to-load drop (%)")
    ax.set_title("28 V bus, 220 W, 8 m harness")
    fig.tight_layout()
    fig.savefig(out_dir / "dc_circuits.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main()
