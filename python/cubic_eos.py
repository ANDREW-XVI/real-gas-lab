#!/usr/bin/env python3
"""
cubic_eos.py  -  PVT behaviour, departure (residual) enthalpy and entropy, and
fugacity of PURE real gases from cubic equations of state.

Equations of state (all written in the generic two-parameter form)

    P = R T / (V - b)  -  a(T) / ((V + eps*b)(V + sigma*b))

    EOS   sigma      eps        Omega     Psi       alpha(Tr)      (Omega, Psi shown rounded)
    vdW   0          0          1/8       27/64     1
    RK    1          0          0.08664   0.42748   Tr**-0.5
    SRK   1          0          0.08664   0.42748   [1 + k(1 - sqrt(Tr))]^2 , k = 0.480 + 1.574 w - 0.176 w^2
    PR    1+sqrt2    1-sqrt2    0.07779   0.45724   [1 + k(1 - sqrt(Tr))]^2 , k = 0.37464 + 1.54226 w - 0.26992 w^2

    b = Omega R Tc / Pc          a(T) = Psi alpha(Tr) R^2 Tc^2 / Pc

Compressibility factor (cubic in Z):
    Z^3 + [(e+s)B - (1+B)] Z^2 + [e s B^2 - (1+B)(e+s)B + A] Z - [(1+B) e s B^2 + A B] = 0
    A = a P/(R T)^2 = Psi alpha Pr/Tr^2 ,  B = b P/(R T) = Omega Pr/Tr

With  I = ln[(Z + sigma B)/(Z + eps B)]/(sigma - eps)   (I = B/Z for vdW)  and  q = A/B :

    ln(phi)  = Z - 1 - ln(Z - B) - q I                      (fugacity coefficient, f = phi P)
    H^R/(RT) = Z - 1 + (dln(alpha)/dln(Tr) - 1) q I          (residual enthalpy)
    S^R/R    = ln(Z - B) + dln(alpha)/dln(Tr) q I            (residual entropy)

Total properties = ideal-gas part (Cp^ig/R = A + B T + C T^2 + D/T^2, reference
T0 = 298.15 K, P0 = 1 bar, H^ig = S^ig = 0) + residual part.

Units: T [K], P [bar], V [cm3/mol], H [J/mol], S [J/(mol K)], f [bar].

Usage
-----
    python cubic_eos.py

The program then asks for everything it needs, step by step: the gas (name,
formula or your own Tc, Pc, omega), the equation of state, the kind of
calculation, and the values of T, P (or V).  Temperatures and pressures may carry
units: 25C, 300K, 80F, 50bar, 5MPa, 700psi, 1atm.  The functions in this file
(analyze, saturation, delta_HS, ...) can also be imported from other Python code.
"""
from __future__ import annotations

import csv
import math
import sys
from dataclasses import dataclass
from typing import Optional

import numpy as np
from scipy.optimize import brentq

R_J = 8.314462618          # J/(mol K)
R_V = 83.14462618          # bar cm3/(mol K)
T0, P0 = 298.15, 1.0       # reference state for the ideal-gas H and S


# --------------------------------------------------------------------------- #
#  Component data  (Tc, Pc, omega: Smith, Van Ness & Abbott, Table B.1;
#                   Cp^ig/R = A + B T + C T^2 + D T^-2 : Table C.1, 298-1500 K)
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Species:
    name: str
    formula: str
    Tc: float      # K
    Pc: float      # bar
    omega: float   # acentric factor
    cp: tuple = (3.5, 0.0, 0.0, 0.0)   # (A, B, C, D) for Cp^ig/R


def _S(name, formula, Tc, Pc, w, A, B, C, D):
    return Species(name, formula, Tc, Pc, w, (A, B, C, D))


_DATA = [
    _S("Methane", "CH4", 190.6, 45.99, 0.012, 1.702, 9.081e-3, -2.164e-6, 0.0),
    _S("Ethane", "C2H6", 305.3, 48.72, 0.100, 1.131, 19.225e-3, -5.561e-6, 0.0),
    _S("Propane", "C3H8", 369.8, 42.48, 0.152, 1.213, 28.785e-3, -8.824e-6, 0.0),
    _S("n-Butane", "C4H10", 425.1, 37.96, 0.200, 1.935, 36.915e-3, -11.402e-6, 0.0),
    _S("Isobutane", "iC4H10", 408.1, 36.48, 0.181, 1.677, 37.853e-3, -11.945e-6, 0.0),
    _S("n-Pentane", "C5H12", 469.7, 33.70, 0.252, 2.464, 45.351e-3, -14.111e-6, 0.0),
    _S("n-Hexane", "C6H14", 507.6, 30.25, 0.301, 3.025, 53.722e-3, -16.791e-6, 0.0),
    _S("Ethylene", "C2H4", 282.3, 50.40, 0.087, 1.424, 14.394e-3, -4.392e-6, 0.0),
    _S("Propylene", "C3H6", 365.0, 46.65, 0.140, 1.637, 22.706e-3, -6.915e-6, 0.0),
    _S("Acetylene", "C2H2", 308.3, 61.39, 0.187, 6.132, 1.952e-3, 0.0, -1.299e5),
    _S("Benzene", "C6H6", 562.2, 48.98, 0.210, -0.206, 39.064e-3, -13.301e-6, 0.0),
    _S("Nitrogen", "N2", 126.2, 34.00, 0.038, 3.280, 0.593e-3, 0.0, 0.040e5),
    _S("Oxygen", "O2", 154.6, 50.43, 0.022, 3.639, 0.506e-3, 0.0, -0.227e5),
    _S("Carbon monoxide", "CO", 132.9, 34.99, 0.048, 3.376, 0.557e-3, 0.0, -0.031e5),
    _S("Carbon dioxide", "CO2", 304.2, 73.83, 0.224, 5.457, 1.045e-3, 0.0, -1.157e5),
    _S("Hydrogen", "H2", 33.19, 13.13, -0.216, 3.249, 0.422e-3, 0.0, 0.083e5),
    _S("Water", "H2O", 647.1, 220.55, 0.345, 3.470, 1.450e-3, 0.0, 0.121e5),
    _S("Ammonia", "NH3", 405.7, 112.80, 0.253, 3.578, 3.020e-3, 0.0, -0.186e5),
    _S("Argon", "Ar", 150.9, 48.98, 0.000, 2.500, 0.0, 0.0, 0.0),
    _S("Hydrogen sulfide", "H2S", 373.5, 89.63, 0.094, 3.931, 1.490e-3, 0.0, -0.232e5),
    _S("Sulfur dioxide", "SO2", 430.8, 78.84, 0.245, 5.699, 0.801e-3, 0.0, -1.015e5),
    _S("Methanol", "CH3OH", 512.6, 80.97, 0.564, 2.211, 12.216e-3, -3.450e-6, 0.0),
    _S("Ethanol", "C2H5OH", 513.9, 61.48, 0.645, 3.518, 20.001e-3, -6.002e-6, 0.0),
]

GASES: dict[str, Species] = {}
for _sp in _DATA:
    GASES[_sp.name.lower().replace(" ", "-")] = _sp
    GASES[_sp.formula.lower()] = _sp
GASES["butane"] = GASES["n-butane"]
GASES["pentane"] = GASES["n-pentane"]
GASES["hexane"] = GASES["n-hexane"]


def get_species(name: str) -> Species:
    key = name.strip().lower().replace(" ", "-")
    if key not in GASES:
        raise KeyError(f"Unknown gas '{name}'. Use 'list' to see the database or supply --Tc/--Pc/--omega.")
    return GASES[key]


# --------------------------------------------------------------------------- #
#  Equation-of-state definitions
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class CubicEOS:
    key: str
    name: str
    sigma: float
    eps: float
    Omega: float
    Psi: float

    def alpha(self, Tr: float, w: float) -> tuple[float, float]:
        """Return (alpha, dln(alpha)/dln(Tr))."""
        if self.key == "VDW":
            return 1.0, 0.0
        if self.key == "RK":
            return Tr ** -0.5, -0.5
        if self.key == "SRK":
            k = 0.480 + 1.574 * w - 0.176 * w * w
        else:  # PR
            k = 0.37464 + 1.54226 * w - 0.26992 * w * w
        s = 1.0 + k * (1.0 - math.sqrt(Tr))
        return s * s, -k * math.sqrt(Tr) / s

    def AB(self, sp: Species, T: float, P: float):
        """Dimensionless A, B and dln(alpha)/dln(Tr)."""
        Tr, Pr = T / sp.Tc, P / sp.Pc
        al, dl = self.alpha(Tr, sp.omega)
        return self.Psi * al * Pr / Tr ** 2, self.Omega * Pr / Tr, dl

    def a_b(self, sp: Species, T: float):
        """Dimensional EOS parameters a(T) [bar cm6/mol2] and b [cm3/mol]."""
        al, _ = self.alpha(T / sp.Tc, sp.omega)
        return (self.Psi * al * R_V ** 2 * sp.Tc ** 2 / sp.Pc,
                self.Omega * R_V * sp.Tc / sp.Pc)

    def pressure(self, sp: Species, T: float, V: float) -> float:
        """Explicit P(T, V) [bar]."""
        a, b = self.a_b(sp, T)
        return R_V * T / (V - b) - a / ((V + self.eps * b) * (V + self.sigma * b))


# Exact constants (the textbook values 0.08664/0.42748 and 0.07779/0.45724 are rounded versions)
_OM_RK = (2 ** (1 / 3) - 1) / 3               # 0.0866403...
_PSI_RK = 1 / (9 * (2 ** (1 / 3) - 1))        # 0.4274802...
_OM_PR, _PSI_PR = 0.0777960739038885, 0.457235528921382
EOS = {
    "VDW": CubicEOS("VDW", "van der Waals", 0.0, 0.0, 1 / 8, 27 / 64),
    "RK": CubicEOS("RK", "Redlich-Kwong", 1.0, 0.0, _OM_RK, _PSI_RK),
    "SRK": CubicEOS("SRK", "Soave-Redlich-Kwong", 1.0, 0.0, _OM_RK, _PSI_RK),
    "PR": CubicEOS("PR", "Peng-Robinson", 1 + math.sqrt(2), 1 - math.sqrt(2), _OM_PR, _PSI_PR),
}


def get_eos(key: str) -> CubicEOS:
    k = key.strip().upper().replace("DER", "D").replace("-", "")
    k = {"VANDERWAALS": "VDW", "VAN": "VDW", "REDLICHKWONG": "RK", "SOAVE": "SRK", "PENGROBINSON": "PR"}.get(k, k)
    if k not in EOS:
        raise KeyError(f"Unknown EOS '{key}'. Choose from VDW, RK, SRK, PR.")
    return EOS[k]


# --------------------------------------------------------------------------- #
#  Roots of the cubic and fugacity coefficient
# --------------------------------------------------------------------------- #
def z_roots(eos: CubicEOS, A: float, B: float) -> list[float]:
    """All physical (Z > B) real roots of the cubic, ascending (liquid-like first)."""
    s, e = eos.sigma, eos.eps
    c2 = (e + s) * B - (1 + B)
    c1 = e * s * B * B - (1 + B) * (e + s) * B + A
    c0 = -((1 + B) * e * s * B * B + A * B)
    r = np.roots([1.0, c2, c1, c0])
    real = [z.real for z in r if abs(z.imag) <= 1e-9 * max(1.0, abs(z.real))]
    if not real:                               # near a triple root: take the most real one
        real = [min(r, key=lambda z: abs(z.imag)).real]
    out = []
    for z0 in real:                            # Newton polishing (keeps tiny liquid roots accurate)
        z = z0
        for _ in range(30):
            f = ((z + c2) * z + c1) * z + c0
            fp = (3 * z + 2 * c2) * z + c1
            if fp == 0.0:
                break
            dz = f / fp
            z -= dz
            if abs(dz) <= 1e-15 * abs(z):
                break
        out.append(z if abs(z - z0) <= 1e-6 * max(abs(z0), 1e-300) else z0)
    out = sorted(z for z in out if z > B * (1 + 1e-12))
    merged: list[float] = []
    for z in out:
        if not merged or abs(z - merged[-1]) > 1e-9 * abs(z):
            merged.append(z)
    if not merged:
        raise ValueError("No physical root (Z > B) found for the cubic.")
    return merged


def _I(eos: CubicEOS, Z: float, B: float) -> float:
    if eos.sigma == eos.eps:
        return B / Z
    return math.log((Z + eos.sigma * B) / (Z + eos.eps * B)) / (eos.sigma - eos.eps)


def ln_phi(eos: CubicEOS, Z: float, A: float, B: float) -> float:
    return Z - 1.0 - math.log(Z - B) - (A / B) * _I(eos, Z, B)


# --------------------------------------------------------------------------- #
#  Ideal-gas contribution
# --------------------------------------------------------------------------- #
def dH_ig(sp: Species, Ta: float, Tb: float) -> float:
    """Integral of Cp^ig dT from Ta to Tb  [J/mol]."""
    A, B, C, D = sp.cp
    return R_J * (A * (Tb - Ta) + B / 2 * (Tb ** 2 - Ta ** 2) + C / 3 * (Tb ** 3 - Ta ** 3)
                  + D * (1 / Ta - 1 / Tb))


def dS_ig_T(sp: Species, Ta: float, Tb: float) -> float:
    """Integral of Cp^ig/T dT from Ta to Tb  [J/(mol K)]."""
    A, B, C, D = sp.cp
    return R_J * (A * math.log(Tb / Ta) + B * (Tb - Ta) + C / 2 * (Tb ** 2 - Ta ** 2)
                  + D / 2 * (1 / Ta ** 2 - 1 / Tb ** 2))


def cp_ig(sp: Species, T: float) -> float:
    A, B, C, D = sp.cp
    return R_J * (A + B * T + C * T * T + D / T ** 2)


def H_ig(sp: Species, T: float) -> float:
    return dH_ig(sp, T0, T)


def S_ig(sp: Species, T: float, P: float) -> float:
    return dS_ig_T(sp, T0, T) - R_J * math.log(P / P0)


# --------------------------------------------------------------------------- #
#  State evaluation
# --------------------------------------------------------------------------- #
@dataclass
class Phase:
    label: str
    T: float
    P: float
    Z: float
    V: float          # cm3/mol
    ln_phi: float
    phi: float
    f: float          # bar
    HR: float         # J/mol
    SR: float         # J/(mol K)
    GR: float         # J/mol
    H: float          # J/mol   (ideal-gas reference + residual)
    S: float          # J/(mol K)
    stable: bool = True


def evaluate(sp: Species, eos: CubicEOS, T: float, P: float, Z: float, label: str = "",
             stable: bool = True) -> Phase:
    A, B, dl = eos.AB(sp, T, P)
    I = _I(eos, Z, B)
    q = A / B
    lnp = Z - 1.0 - math.log(Z - B) - q * I
    HR = R_J * T * (Z - 1.0 + (dl - 1.0) * q * I)
    SR = R_J * (math.log(Z - B) + dl * q * I)
    return Phase(label, T, P, Z, Z * R_V * T / P, lnp, math.exp(lnp), P * math.exp(lnp),
                 HR, SR, R_J * T * lnp, H_ig(sp, T) + HR, S_ig(sp, T, P) + SR, stable)


def spinodal_pressures(sp: Species, eos: CubicEOS, T: float) -> tuple[float, float, float, float]:
    """Pressures and volumes at the two extrema of the isotherm: (P_low, P_high, V_low, V_high)."""
    T = float(T)                       # (np.float64 * poly1d would silently return an ndarray)
    a, b = eos.a_b(sp, T)
    s, e = eos.sigma, eos.eps
    poly = (a * np.poly1d([2.0, (e + s) * b]) * np.poly1d([1.0, -b]) ** 2
            - R_V * T * np.poly1d([1.0, e * b]) ** 2 * np.poly1d([1.0, s * b]) ** 2)
    cand = [v for v in poly.roots if v.real > b * (1 + 1e-9)]
    Vs = sorted(v.real for v in cand if abs(v.imag) < 1e-9 * abs(v))
    if len(Vs) < 2 and T < sp.Tc:      # very close to Tc the two extrema nearly merge: accept a near-real pair
        near = sorted(cand, key=lambda v: abs(v.imag))[:2]
        if len(near) == 2 and all(abs(v.imag) < 1e-4 * abs(v) for v in near):
            Vs = sorted(v.real for v in near)
    if len(Vs) < 2:
        raise ValueError("No van der Waals loop at this temperature (T >= Tc for this EOS).")
    Va, Vb = Vs[0], Vs[-1]
    return eos.pressure(sp, T, Va), eos.pressure(sp, T, Vb), Va, Vb


def saturation(sp: Species, eos: CubicEOS, T: float):
    """Vapour pressure from equal fugacity.  Returns (Psat, liquid Phase, vapour Phase)."""
    if T >= sp.Tc:
        raise ValueError(f"T = {T} K is not below Tc = {sp.Tc} K: no vapour-liquid equilibrium.")
    Pl, Ph, Va, Vb = spinodal_pressures(sp, eos, T)
    Vmid = 0.5 * (Va + Vb)
    lo = max(Pl, 1e-12 * sp.Pc)

    def g(lnP):
        P = math.exp(lnP)
        A, B, _ = eos.AB(sp, T, P)
        zs = z_roots(eos, A, B)
        if len(zs) < 2:                         # only one root: decide on which side we are
            return -1.0 if zs[0] * R_V * T / P < Vmid else 1.0
        return ln_phi(eos, zs[0], A, B) - ln_phi(eos, zs[-1], A, B)

    x_lo, x_hi = math.log(lo * (1 + 1e-9)), math.log(Ph * (1 - 1e-9))
    if g(x_lo) * g(x_hi) > 0:
        raise ValueError("Could not bracket the saturation pressure at this temperature.")
    Ps = math.exp(brentq(g, x_lo, x_hi, xtol=1e-14, rtol=1e-14, maxiter=300))
    A, B, _ = eos.AB(sp, T, Ps)
    zs = z_roots(eos, A, B)
    return (Ps, evaluate(sp, eos, T, Ps, zs[0], "saturated liquid"),
            evaluate(sp, eos, T, Ps, zs[-1], "saturated vapour"))


@dataclass
class Result:
    sp: Species
    eos: CubicEOS
    T: float
    P: float
    phase: str                  # description of the stable state
    stable: Phase               # properties of the stable phase
    roots: list                 # all real roots as Phase objects
    Psat: Optional[float]
    a: float
    b: float


def analyze(sp: Species, eos: CubicEOS, T: float, P: float) -> Result:
    """Full (T, P) analysis: all roots, stable root (min fugacity), phase label, properties."""
    T, P = float(T), float(P)
    if T <= 0 or P <= 0:
        raise ValueError("T and P must be positive.")
    A, B, _ = eos.AB(sp, T, P)
    zs = z_roots(eos, A, B)
    phases = [evaluate(sp, eos, T, P, z) for z in zs]
    stable = min(phases, key=lambda p: p.ln_phi)

    Psat = None
    if T < sp.Tc:
        try:
            Psat = saturation(sp, eos, T)[0]
        except ValueError:
            Psat = None
    if T >= sp.Tc:
        desc = "supercritical fluid (T > Tc)" if P >= sp.Pc else "gas (T > Tc)"
    elif Psat is not None and abs(P - Psat) <= 1e-6 * Psat:
        desc = "saturated vapour-liquid equilibrium (P = Psat)"
    elif Psat is not None:
        desc = "superheated vapour / gas" if P < Psat else "compressed (sub-cooled) liquid"
    else:
        desc = "vapour" if stable.Z > 0.3 else "liquid"
    for p in phases:
        p.stable = p is stable
        if len(phases) == 3:
            p.label = "liquid-like root" if p is phases[0] else ("vapour-like root" if p is phases[-1]
                                                                  else "middle root (unstable)")
        else:
            p.label = "only real root"
    a, b = eos.a_b(sp, T)
    return Result(sp, eos, T, P, desc, stable, phases, Psat, a, b)


def analyze_TV(sp: Species, eos: CubicEOS, T: float, V: float) -> Result:
    """(T, V) given: P from the EOS, then the properties of that state."""
    P = eos.pressure(sp, T, V)
    if P <= 0:
        raise ValueError(f"EOS gives P = {P:.4g} bar <= 0 at this (T, V): inside the unstable/tension region.")
    A, B, _ = eos.AB(sp, T, P)
    Z = P * V / (R_V * T)
    res = analyze(sp, eos, T, P)
    st = evaluate(sp, eos, T, P, Z, "given (T,V) state")
    res.stable = st
    # flag mechanically unstable volumes (dP/dV > 0)
    Pl, Ph, Va, Vb = (None,) * 4
    if T < sp.Tc:
        try:
            Pl, Ph, Va, Vb = spinodal_pressures(sp, eos, T)
        except ValueError:
            pass
    if Va is not None and Va < V < Vb:
        res.phase = "inside the spinodal (mechanically unstable) - not a physical single phase"
    return res


def delta_HS(sp: Species, eos: CubicEOS, T1, P1, T2, P2):
    """Property changes between two states: returns (dH, dS, state1, state2)."""
    s1 = analyze(sp, eos, T1, P1).stable
    s2 = analyze(sp, eos, T2, P2).stable
    return s2.H - s1.H, s2.S - s1.S, s1, s2


# --------------------------------------------------------------------------- #
#  Printing
# --------------------------------------------------------------------------- #
def report(res: Result) -> str:
    sp, eos, L = res.sp, res.eos, []
    A, B, _ = eos.AB(sp, res.T, res.P)
    L.append("=" * 72)
    L.append(f" {sp.name} ({sp.formula})   EOS: {eos.name}")
    L.append(f" Tc = {sp.Tc} K   Pc = {sp.Pc} bar   omega = {sp.omega}")
    L.append(f" T = {res.T:.4f} K   P = {res.P:.5g} bar   (Tr = {res.T / sp.Tc:.4f}, Pr = {res.P / sp.Pc:.4f})")
    L.append(f" a(T) = {res.a:.6g} bar cm6/mol2    b = {res.b:.6g} cm3/mol    A = {A:.6g}   B = {B:.6g}")
    L.append("-" * 72)
    L.append(f" Phase (stable): {res.phase}")
    if res.Psat is not None:
        L.append(f" Vapour pressure from the EOS at this T: Psat = {res.Psat:.6g} bar")
    L.append(" Real roots of the cubic:")
    for p in res.roots:
        mark = "  <== stable (lowest fugacity)" if p.stable and len(res.roots) > 1 else ""
        L.append(f"   Z = {p.Z:.6f}  V = {p.V:12.5f} cm3/mol  ln(phi) = {p.ln_phi:+.6f}  [{p.label}]{mark}")
    s = res.stable
    L.append("-" * 72)
    L.append(f" Z            = {s.Z:.6f}")
    L.append(f" V            = {s.V:.5f} cm3/mol      density = {1000.0 / s.V:.5f} mol/L")
    L.append(f" phi          = {s.phi:.6f}     ln(phi) = {s.ln_phi:+.6f}")
    L.append(f" f            = {s.f:.6g} bar")
    L.append(f" H^R          = {s.HR:12.3f} J/mol      S^R = {s.SR:10.5f} J/(mol K)")
    L.append(f" G^R = RT ln(phi) = {s.GR:12.3f} J/mol")
    L.append(f" H (total)    = {s.H:12.3f} J/mol      S   = {s.S:10.5f} J/(mol K)"
             f"   [H = S = 0 for ideal gas at {T0} K, {P0} bar]")
    L.append("=" * 72)
    return "\n".join(L)


def sat_report(sp, eos, T) -> str:
    Ps, liq, vap = saturation(sp, eos, T)
    dH, dS = vap.H - liq.H, vap.S - liq.S
    L = ["=" * 72, f" {sp.name}  {eos.name}  saturation at T = {T:.3f} K (Tr = {T / sp.Tc:.4f})", "-" * 72,
         f" Psat = {Ps:.6g} bar      ln(phi) liquid = {liq.ln_phi:+.8f}   vapour = {vap.ln_phi:+.8f}",
         f" f_sat = {vap.f:.6g} bar     phi_sat = {vap.phi:.6f}",
         f" Z_L = {liq.Z:.6f}   V_L = {liq.V:.4f} cm3/mol        Z_V = {vap.Z:.6f}   V_V = {vap.V:.4f} cm3/mol",
         f" H^R_L = {liq.HR:.2f}  H^R_V = {vap.HR:.2f} J/mol",
         f" S^R_L = {liq.SR:.4f}  S^R_V = {vap.SR:.4f} J/(mol K)",
         f" dH_vap = {dH:.2f} J/mol   dS_vap = {dS:.4f} J/(mol K)   (dH_vap/T = {dH / T:.4f})", "=" * 72]
    return "\n".join(L)


def table(sp, eos, T, Ps, csv_path=None) -> str:
    hdr = ["P/bar", "phase", "Z", "V/cm3mol-1", "phi", "f/bar", "HR/Jmol-1", "SR/JmolK-1", "H/Jmol-1", "S/JmolK-1"]
    rows = []
    for P in Ps:
        r = analyze(sp, eos, T, P)
        s = r.stable
        ph = r.phase.split()[0] if r.phase else ""
        rows.append([P, ph, s.Z, s.V, s.phi, s.f, s.HR, s.SR, s.H, s.S])
    if csv_path:
        with open(csv_path, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(hdr)
            w.writerows(rows)
    out = [f"{sp.name}, {eos.name}, T = {T} K",
           f"{'P/bar':>9} {'phase':>11} {'Z':>9} {'V cm3/mol':>12} {'phi':>9} {'f/bar':>10} {'H^R J/mol':>11} {'S^R J/molK':>11} {'H J/mol':>11} {'S J/molK':>10}"]
    for r in rows:
        out.append(f"{r[0]:9.4g} {r[1][:11]:>11} {r[2]:9.5f} {r[3]:12.4f} {r[4]:9.5f} {r[5]:10.5g} {r[6]:11.2f} {r[7]:11.4f} {r[8]:11.1f} {r[9]:10.3f}")
    return "\n".join(out)


def compare(sp, T, P=None, V=None) -> str:
    out = [f"{sp.name}: T = {T} K, " + (f"P = {P} bar" if P else f"V = {V} cm3/mol"),
           f"{'EOS':>5} {'P/bar':>9} {'Z':>9} {'V cm3/mol':>11} {'phi':>9} {'H^R J/mol':>11} {'S^R J/molK':>11}  phase"]
    for key, eos in EOS.items():
        r = analyze(sp, eos, T, P) if P else analyze_TV(sp, eos, T, V)
        s = r.stable
        out.append(f"{key:>5} {r.P:9.4g} {s.Z:9.5f} {s.V:11.4f} {s.phi:9.5f} {s.HR:11.2f} {s.SR:11.4f}  {r.phase}")
    return "\n".join(out)


# --------------------------------------------------------------------------- #
#  Input helpers
# --------------------------------------------------------------------------- #
def parse_T(x) -> float:
    s = str(x).strip().lower()
    if s.endswith("k"):
        return float(s[:-1])
    if s.endswith("c"):
        return float(s[:-1]) + 273.15
    if s.endswith("f"):
        return (float(s[:-1]) - 32) * 5 / 9 + 273.15
    return float(s)


def parse_P(x) -> float:
    s = str(x).strip().lower()
    for suf, fac in (("mpa", 10.0), ("kpa", 0.01), ("bar", 1.0), ("atm", 1.01325), ("psi", 0.0689475729), ("pa", 1e-5)):
        if s.endswith(suf):
            return float(s[:-len(suf)]) * fac
    return float(s)


# --------------------------------------------------------------------------- #
#  Guided, prompt-driven interface
# --------------------------------------------------------------------------- #
def _ask(prompt, default=None, conv=str):
    """Prompt until the answer can be converted; Enter accepts the default."""
    while True:
        try:
            s = input(f"{prompt}" + (f" [{default}]" if default is not None else "") + ": ").strip()
        except EOFError:
            raise SystemExit("\nInput ended - goodbye.")
        if not s and default is not None:
            s = str(default)
        try:
            return conv(s)
        except Exception:                                  # noqa: BLE001
            print("  Could not understand that value - please try again.")


def _ask_T(label="Temperature T", default=None):
    return _ask(f"{label} [K]  (also accepts e.g. 25C, 80F)", default, parse_T)


def _ask_P(label="Pressure P", default=None):
    return _ask(f"{label} [bar]  (also accepts e.g. 5MPa, 1atm, 700psi)", default, parse_P)


def list_gases():
    print(f"\n  {'name':<18}{'formula':<9}{'Tc/K':>8}{'Pc/bar':>9}{'omega':>8}")
    for sp in _DATA:
        print(f"  {sp.name:<18}{sp.formula:<9}{sp.Tc:8.2f}{sp.Pc:9.2f}{sp.omega:8.3f}")
    print()


def choose_gas(current: Optional[Species] = None) -> Species:
    default = current.name if current is not None else "methane"
    while True:
        s = input("Gas name or formula (e.g. methane, CO2, n-butane; 'list' = show all, "
                  f"'custom' = enter your own data) [{default}]: ").strip()
        if not s:
            return current if current is not None else get_species(default)
        if s.lower() == "list":
            list_gases()
        elif s.lower() == "custom":
            print("Enter the critical properties and (optionally) ideal-gas heat-capacity data of your gas.")
            Tc = _ask("  Critical temperature Tc [K]", conv=float)
            Pc = _ask("  Critical pressure Pc [bar]", conv=float)
            w = _ask("  Acentric factor omega", conv=float)
            cp = _ask("  Cp_ig/R = A + B*T + C*T^2 + D/T^2 -> enter A B C D (needed for total H and S)",
                      "3.5 0 0 0", lambda z: tuple(map(float, z.split()))[:4] if len(z.split()) == 4 else 1 / 0)
            return Species("Custom gas", "-", Tc, Pc, w, cp)
        else:
            try:
                return get_species(s)
            except KeyError as ex:
                print(" ", ex)


def choose_eos(current: Optional[CubicEOS] = None) -> CubicEOS:
    print("Equation of state:  1 = van der Waals   2 = Redlich-Kwong   3 = Soave-Redlich-Kwong   4 = Peng-Robinson")
    keys = {"1": "VDW", "2": "RK", "3": "SRK", "4": "PR"}
    default = current.key if current is not None else "PR"
    while True:
        s = input(f"Choose EOS (number or name) [{default}]: ").strip() or default
        try:
            return get_eos(keys.get(s, s))
        except KeyError as ex:
            print(" ", ex)


MENU = """
What would you like to calculate?
  1) Z, V, fugacity, H and S at a given T and P
  2) Same, but at a given T and molar volume V (P is calculated)
  3) Vapour pressure, saturated fugacity and heat of vaporisation at a given T
  4) Table of properties versus pressure at a fixed T
  5) Change in H and S between two states (T1,P1) -> (T2,P2)
  6) Compare all four equations of state at a given T and P
  7) Change gas      8) Change EOS      9) List gases      0) Quit"""


def interactive():
    print("=" * 72)
    print(" Real-gas properties from cubic equations of state")
    print(" PVT  |  residual and total H and S  |  fugacity and fugacity coefficient")
    print("=" * 72)
    sp = choose_gas()
    eos = choose_eos()
    while True:
        print(f"\n--- Current gas: {sp.name}   EOS: {eos.name} " + "-" * 20)
        print(MENU)
        c = _ask("Enter your choice", "1")
        try:
            if c == "0":
                print("Goodbye.")
                return
            elif c == "1":
                T = _ask_T()
                P = _ask_P()
                print(report(analyze(sp, eos, T, P)))
            elif c == "2":
                T = _ask_T()
                V = _ask("Molar volume V [cm3/mol]", conv=float)
                print(report(analyze_TV(sp, eos, T, V)))
            elif c == "3":
                print(f"(Vapour-liquid equilibrium exists only below Tc = {sp.Tc} K)")
                print(sat_report(sp, eos, _ask_T()))
            elif c == "4":
                T = _ask_T()
                p1 = _ask_P("Lowest pressure")
                p2 = _ask_P("Highest pressure")
                n = _ask("Number of pressure points", 10, int)
                f = input("Save the table to a CSV file? Enter a file name or press Enter to skip: ").strip()
                print(table(sp, eos, T, np.linspace(p1, p2, max(n, 2)), f or None))
                if f:
                    print(f"(saved to {f})")
            elif c == "5":
                print("State 1:")
                T1, P1 = _ask_T("  T1"), _ask_P("  P1")
                print("State 2:")
                T2, P2 = _ask_T("  T2"), _ask_P("  P2")
                dH, dS, s1, s2 = delta_HS(sp, eos, T1, P1, T2, P2)
                print(f"\n  State 1: Z = {s1.Z:.5f}, H^R = {s1.HR:.2f} J/mol, S^R = {s1.SR:.4f} J/(mol K)")
                print(f"  State 2: Z = {s2.Z:.5f}, H^R = {s2.HR:.2f} J/mol, S^R = {s2.SR:.4f} J/(mol K)")
                print(f"  Change:  dH = {dH:.2f} J/mol     dS = {dS:.4f} J/(mol K)")
            elif c == "6":
                T = _ask_T()
                P = _ask_P()
                print(compare(sp, T, P=P))
            elif c == "7":
                sp = choose_gas(sp)
            elif c == "8":
                eos = choose_eos(eos)
            elif c == "9":
                list_gases()
            else:
                print("Please enter a number from 0 to 9.")
        except (ValueError, KeyError) as ex:
            print("  Error:", ex)


def main():
    try:
        interactive()
    except (KeyboardInterrupt, SystemExit) as ex:
        print(ex if isinstance(ex, SystemExit) and ex.code else "\nGoodbye.")


if __name__ == "__main__":
    main()
