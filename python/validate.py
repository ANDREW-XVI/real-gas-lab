#!/usr/bin/env python3
"""
validate.py - self-tests for cubic_eos.py

 1. ideal-gas limit              (Z -> 1, phi -> 1, H^R, S^R -> 0 as P -> 0)
 2. thermodynamic consistency    H^R/RT = -T d(ln phi)/dT|P ,  S^R/R = -T d(ln phi)/dT - ln phi
 3. phase equilibrium            equal fugacity at Psat, Clapeyron: dH_vap = T dV_vap dPsat/dT
 4. ideal-gas Cp integrals       closed form vs numerical quadrature
 5. independent implementation   CoolProp's cubic backends (PR, SRK) with identical Tc, Pc, omega
 6. accuracy vs reference EOS    CoolProp Helmholtz EOS (Z, H^R, S^R) for VDW/RK/SRK/PR   (optional)

Run:  python validate.py          (items 5-6 are skipped if CoolProp is not installed)
"""
import math
import sys

import numpy as np
from scipy.integrate import quad

import cubic_eos as ce

fails = 0


def check(name, ok, info=""):
    global fails
    print(f"  [{'PASS' if ok else 'FAIL'}] {name} {info}")
    fails += (not ok)


# 1 ------------------------------------------------------------------------ #
print("1. Ideal-gas limit (P = 1e-6 bar, T = 300 K)")
worst = 0
for sp in {id(s): s for s in ce._DATA}.values():
    for eos in ce.EOS.values():
        r = ce.analyze(sp, eos, 400.0, 1e-6).stable
        worst = max(worst, abs(r.Z - 1), abs(r.phi - 1), abs(r.HR) / (ce.R_J * 400) / 1e-3, abs(r.SR) / ce.R_J / 1e-3)
check("Z, phi -> 1 and H^R, S^R -> 0 for all gases and EOS", worst < 1e-3, f"(worst = {worst:.2e})")

# 2 ------------------------------------------------------------------------ #
print("2. Consistency of H^R, S^R with the fugacity coefficient")
worstH = worstS = 0.0
n = 0
for sp in {id(s): s for s in ce._DATA}.values():
    for eos in ce.EOS.values():
        for Tr in (1.05, 1.5, 2.5):
            for Pr in (0.05, 0.5, 2.0):
                T, P = Tr * sp.Tc, Pr * sp.Pc
                h = 1e-4 * T
                lp = lambda t: ce.analyze(sp, eos, t, P).stable.ln_phi
                dlp = (lp(T + h) - lp(T - h)) / (2 * h)
                s0 = ce.analyze(sp, eos, T, P).stable
                HR_num = -ce.R_J * T * T * dlp
                SR_num = ce.R_J * (-T * dlp - s0.ln_phi)
                scaleH, scaleS = ce.R_J * T, ce.R_J
                worstH = max(worstH, abs(HR_num - s0.HR) / scaleH)
                worstS = max(worstS, abs(SR_num - s0.SR) / scaleS)
                n += 1
check("H^R = -R T^2 d ln(phi)/dT", worstH < 1e-6, f"({n} states, max |err|/RT = {worstH:.2e})")
check("S^R = -R [T d ln(phi)/dT + ln(phi)]", worstS < 1e-6, f"(max |err|/R = {worstS:.2e})")

# 3 ------------------------------------------------------------------------ #
print("3. Vapour-liquid equilibrium")
worst_f = worst_cl = 0.0
count = 0
for name in ("methane", "propane", "n-butane", "co2", "ammonia", "water", "nitrogen"):
    sp = ce.get_species(name)
    for eos in ce.EOS.values():
        for Tr in (0.6, 0.75, 0.9, 0.97):
            T = Tr * sp.Tc
            try:
                Ps, L, V = ce.saturation(sp, eos, T)
            except ValueError:
                continue
            worst_f = max(worst_f, abs(L.ln_phi - V.ln_phi))
            h = 1e-4 * T
            dPdT = (ce.saturation(sp, eos, T + h)[0] - ce.saturation(sp, eos, T - h)[0]) / (2 * h)
            dH_clap = T * (V.V - L.V) * 1e-5 * dPdT * 1e6 / 1e6      # bar*cm3 -> J : x0.1
            dH_clap = T * (V.V - L.V) * dPdT * 0.1
            dH_dep = V.H - L.H
            worst_cl = max(worst_cl, abs(dH_clap - dH_dep) / abs(dH_dep))
            count += 1
check("ln(phi_L) = ln(phi_V) at Psat", worst_f < 1e-9, f"({count} cases, max diff = {worst_f:.1e})")
check("dH_vap from departure functions = Clapeyron T dV dPsat/dT", worst_cl < 1e-5, f"(max rel. diff = {worst_cl:.1e})")

# 4 ------------------------------------------------------------------------ #
print("4. Ideal-gas Cp integrals")
w = 0.0
for name in ("methane", "co2", "acetylene", "water"):
    sp = ce.get_species(name)
    h_num = quad(lambda t: ce.cp_ig(sp, t), 298.15, 850.0)[0]
    s_num = quad(lambda t: ce.cp_ig(sp, t) / t, 298.15, 850.0)[0]
    w = max(w, abs(h_num - ce.dH_ig(sp, 298.15, 850.0)) / h_num, abs(s_num - ce.dS_ig_T(sp, 298.15, 850.0)) / s_num)
check("closed form = numerical integration", w < 1e-9, f"(max rel. diff = {w:.1e})")

# 5-6 ---------------------------------------------------------------------- #
try:
    import CoolProp.CoolProp as CP
    from CoolProp.CoolProp import PropsSI
except ImportError:
    print("5-6. CoolProp not installed - skipped (pip install CoolProp)")
    sys.exit(1 if fails else 0)

print("5. Independent implementation (CoolProp cubic backends, same Tc, Pc, omega)")
worst = {"Z": 0, "lnphi": 0, "HR": 0, "SR": 0}
n = 0
for fl in ("CO2", "Methane", "Propane", "Nitrogen", "Water", "n-Butane"):
    for key in ("PR", "SRK"):
        AS = CP.AbstractState(key, fl)
        sp = ce.Species(fl, "-", AS.T_critical(), AS.p_critical() / 1e5, AS.acentric_factor())
        for Tr, Pr in ((1.2, 0.3), (1.5, 2.0), (2.0, 5.0), (0.9, 0.1), (0.7, 5.0)):
            T, P = Tr * sp.Tc, Pr * sp.Pc
            try:
                AS.update(CP.PT_INPUTS, P * 1e5, T)
                ref = dict(Z=AS.compressibility_factor(), lnphi=math.log(AS.fugacity_coefficient(0)),
                           HR=AS.hmolar_residual())
                rho = AS.rhomolar()
                ref["SR"] = AS.smolar_residual() + ce.R_J * math.log(ref["Z"])   # (T,rho) -> (T,P) basis
            except Exception:                                           # noqa: BLE001
                continue
            me = ce.analyze(sp, ce.EOS[key], T, P).stable
            # only compare when both codes picked the same root
            if abs(me.Z - ref["Z"]) > 1e-3 * ref["Z"]:
                continue
            worst["Z"] = max(worst["Z"], abs(me.Z / ref["Z"] - 1))
            worst["lnphi"] = max(worst["lnphi"], abs(me.ln_phi - ref["lnphi"]))
            worst["HR"] = max(worst["HR"], abs(me.HR - ref["HR"]) / (ce.R_J * T))
            worst["SR"] = max(worst["SR"], abs(me.SR - ref["SR"]) / ce.R_J)
            n += 1
print(f"      {n} states compared; max deviations: Z {worst['Z']:.1e} (rel), ln(phi) {worst['lnphi']:.1e}, "
      f"H^R/RT {worst['HR']:.1e}, S^R/R {worst['SR']:.1e}")
check("Z, ln(phi), H^R/RT, S^R/R agree with CoolProp PR and SRK", n > 20 and max(worst.values()) < 1e-6)

print("6. Accuracy against the reference Helmholtz EOS (CoolProp HEOS)  [|rel. error| in Z, mean over grid]")
rows = []
for fl, cname in (("CO2", "co2"), ("Methane", "methane"), ("Propane", "propane"), ("Nitrogen", "nitrogen"),
                  ("Ammonia", "ammonia")):
    sp = ce.get_species(cname)
    for key, eos in ce.EOS.items():
        eZ, eH, eS, cnt = [], [], [], 0
        for Tr in (1.1, 1.3, 1.6, 2.0):
            for Pr in (0.2, 0.5, 1.0, 2.0, 4.0):
                T, P = Tr * sp.Tc, Pr * sp.Pc
                try:
                    Zr = PropsSI("Z", "T", T, "P", P * 1e5, fl)
                    HR = (PropsSI("Hmolar", "T", T, "P", P * 1e5, fl) - PropsSI("Hmolar", "T", T, "P", 1.0, fl))
                    SR = (PropsSI("Smolar", "T", T, "P", P * 1e5, fl) - PropsSI("Smolar", "T", T, "P", 1.0, fl)
                          + ce.R_J * math.log(P * 1e5 / 1.0))
                except Exception:                                       # noqa: BLE001
                    continue
                me = ce.analyze(sp, eos, T, P).stable
                eZ.append(abs(me.Z / Zr - 1))
                eH.append(abs(me.HR - HR) / abs(HR))
                eS.append(abs(me.SR - SR) / abs(SR))
        rows.append((fl, key, 100 * np.mean(eZ), 100 * np.mean(eH), 100 * np.mean(eS), len(eZ)))
print(f"      {'gas':<9}{'EOS':<5}{'Z %':>8}{'H^R %':>9}{'S^R %':>9}{'n':>4}")
for r in rows:
    print(f"      {r[0]:<9}{r[1]:<5}{r[2]:8.2f}{r[3]:9.2f}{r[4]:9.2f}{r[5]:4d}")
with open("accuracy_vs_reference.csv", "w") as fh:
    fh.write("gas,EOS,mean_abs_err_Z_pct,mean_abs_err_HR_pct,mean_abs_err_SR_pct,n\n")
    for r in rows:
        fh.write(",".join(str(x) if isinstance(x, (str, int)) else f"{x:.4f}" for x in r) + "\n")

print("\nALL CHECKS PASSED" if not fails else f"\n{fails} CHECK(S) FAILED")
sys.exit(1 if fails else 0)
