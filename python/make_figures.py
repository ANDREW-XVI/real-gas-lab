#!/usr/bin/env python3
"""make_figures.py - figures and worked-example numbers for the report (needs matplotlib; CoolProp optional)."""
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import cubic_eos as ce

try:
    from CoolProp.CoolProp import PropsSI
except ImportError:
    PropsSI = None

os.makedirs("../docs/figures", exist_ok=True)
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": "#52514e", "axes.labelcolor": "#0b0b0b", "xtick.color": "#52514e",
                     "ytick.color": "#52514e", "axes.grid": True, "grid.color": "#e6e5e1", "grid.linewidth": 0.5,
                     "legend.frameon": False, "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb",
                     "savefig.facecolor": "#fcfcfb"})
# categorical slots 1-4 of the validated default palette, fixed order; line style is the secondary encoding
STYLE = {"VDW": ("#2a78d6", "-"), "RK": ("#eb6834", "--"), "SRK": ("#1baf7a", "-."), "PR": ("#eda100", ":")}
INK = "#0b0b0b"
LBL = {"VDW": "vdW", "RK": "RK", "SRK": "SRK", "PR": "PR"}

# ---------------------------------------------------------------- Figure 1
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.5))

sp, eos = ce.get_species("n-butane"), ce.EOS["PR"]
a = ax[0]
cols = ["#2a78d6", "#1baf7a", "#eb6834", "#6250d6"]
for Tr, c in zip((0.85, 0.95, 1.0, 1.15), cols):
    T = Tr * sp.Tc
    _, b = eos.a_b(sp, T)
    V = np.logspace(math.log10(b * 1.02), 3.9, 600)
    P = np.array([eos.pressure(sp, T, v) for v in V])
    a.plot(V, P, color=c, lw=1.3, label=f"$T_r$ = {Tr}")
Ts = np.linspace(0.55, 0.998, 90) * sp.Tc
dome = [ce.saturation(sp, eos, T) for T in Ts]
VL = [d[1].V for d in dome]; VV = [d[2].V for d in dome]; PS = [d[0] for d in dome]
a.plot(VL + VV[::-1], PS + PS[::-1], color=INK, lw=1.0)
Zcc = ce.EOS['PR']
for T in (0.999 * sp.Tc,):
    pass
a.set_xscale("log"); a.set_ylim(0, 60); a.set_xlim(60, 8000)
a.set_xlabel("V (cm$^3$/mol)"); a.set_ylabel("P (bar)"); a.set_title("(a) n-butane, PR: isotherms + dome", fontsize=8)
a.legend(loc="center right", bbox_to_anchor=(1.0, 0.68), fontsize=7)

# (b) Z vs P for CO2 at 350 K
sp, T = ce.get_species("co2"), 350.0
Pg = np.linspace(2, 300, 60)
b = ax[1]
if PropsSI:
    Zr = [PropsSI("Z", "T", T, "P", p * 1e5, "CO2") for p in Pg]
    b.plot(Pg, Zr, "o", color=INK, ms=2.6, label="reference EOS", zorder=5)
for k, (c, ls) in STYLE.items():
    b.plot(Pg, [ce.analyze(sp, ce.EOS[k], T, p).stable.Z for p in Pg], color=c, ls=ls, lw=1.4, label=LBL[k])
b.set_xlabel("P (bar)"); b.set_ylabel("Z"); b.set_title("(b) CO$_2$ at 350 K: Z vs P", fontsize=8)
b.legend(fontsize=6.5, loc="upper right")

# (c) ln(phi) of all roots, n-butane, PR, 350 K  (fugacity 'swallow tail')
sp, eos, T = ce.get_species("n-butane"), ce.EOS["PR"], 350.0
Ps, L, V = ce.saturation(sp, eos, T)
Pl, Ph, _, _ = ce.spinodal_pressures(sp, eos, T)
c = ax[2]
Pp = np.linspace(max(Pl, 0.5), Ph, 400)
curves = {0: [], 1: [], 2: []}
for p in Pp:
    A, B, _ = eos.AB(sp, T, p)
    zs = ce.z_roots(eos, A, B)
    if len(zs) == 3:
        for i, z in enumerate(zs):
            curves[i].append((p, ce.ln_phi(eos, z, A, B)))
for i, (lab, col, ls) in enumerate((("liquid root", "#2a78d6", "-"), ("middle root", "#8a8984", "--"), ("vapour root", "#eb6834", "-"))):
    arr = np.array(curves[i]); c.plot(arr[:, 0], arr[:, 1], color=col, ls=ls, lw=1.4, label=lab)
c.plot([Ps], [L.ln_phi], "o", color=INK, ms=4, zorder=5)
c.annotate(f"$P^{{sat}}$ = {Ps:.2f} bar", (Ps, L.ln_phi), xytext=(Ps + 2.6, 0.45), fontsize=7,
           arrowprops=dict(arrowstyle="-", color="#52514e", lw=0.6))
c.set_xlabel("P (bar)"); c.set_ylabel("ln $\\phi$"); c.set_title("(c) n-butane, PR, 350 K: roots", fontsize=8)
c.set_ylim(-0.6, 1.6); c.legend(fontsize=6.5, loc="upper right")
fig.tight_layout(w_pad=0.8)
fig.savefig("../docs/figures/fig1.png", dpi=300)
plt.close(fig)

# ---------------------------------------------------------------- Figure 2
sp, T = ce.get_species("co2"), 350.0
Pg = np.linspace(2, 300, 60)
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.3))
ref = {}
if PropsSI:
    ref["phi"] = [math.exp(0) for _ in Pg]
    HR, SR, PHI = [], [], []
    for p in Pg:
        H = PropsSI("Hmolar", "T", T, "P", p * 1e5, "CO2") - PropsSI("Hmolar", "T", T, "P", 1.0, "CO2")
        S = (PropsSI("Smolar", "T", T, "P", p * 1e5, "CO2") - PropsSI("Smolar", "T", T, "P", 1.0, "CO2")
             + ce.R_J * math.log(p * 1e5))
        HR.append(H / 1000); SR.append(S)
        # fugacity coefficient from the reference EOS: ln phi = (H^R - T S^R)/(R T)
        PHI.append(math.exp((H - T * S) / (ce.R_J * T)))
    ref = {"HR": HR, "SR": SR, "phi": PHI}
panels = (("phi", "$\\phi$", "(a) fugacity coefficient", lambda s: s.phi),
          ("HR", "$H^R$ (kJ/mol)", "(b) residual enthalpy", lambda s: s.HR / 1000),
          ("SR", "$S^R$ (J/mol K)", "(c) residual entropy", lambda s: s.SR))
for a, (key, yl, ttl, fn) in zip(ax, panels):
    if ref:
        a.plot(Pg, ref[key], "o", color=INK, ms=2.6, label="reference EOS", zorder=5)
    for k, (c, ls) in STYLE.items():
        a.plot(Pg, [fn(ce.analyze(sp, ce.EOS[k], T, p).stable) for p in Pg], color=c, ls=ls, lw=1.4, label=LBL[k])
    a.set_xlabel("P (bar)"); a.set_ylabel(yl); a.set_title(ttl + " (CO$_2$, 350 K)", fontsize=7.5)
ax[0].legend(fontsize=6.5, loc="upper right")
fig.tight_layout(w_pad=0.8)
fig.savefig("../docs/figures/fig2.png", dpi=300)
plt.close(fig)

# ---------------------------------------------------------------- worked numbers for the report
print("\nExample 1: CO2, 350 K, 100 bar  (all EOS + reference)")
for k, e in ce.EOS.items():
    s = ce.analyze(sp, e, 350, 100).stable
    print(f"  {k:4s} Z={s.Z:.4f}  V={s.V:7.2f}  phi={s.phi:.4f}  f={s.f:6.2f}  HR={s.HR:8.1f}  SR={s.SR:7.3f}")
if PropsSI:
    H = PropsSI("Hmolar", "T", 350, "P", 100e5, "CO2") - PropsSI("Hmolar", "T", 350, "P", 1.0, "CO2")
    S = (PropsSI("Smolar", "T", 350, "P", 100e5, "CO2") - PropsSI("Smolar", "T", 350, "P", 1.0, "CO2") + ce.R_J * math.log(100e5))
    print(f"  ref  Z={PropsSI('Z','T',350,'P',100e5,'CO2'):.4f}  HR={H:8.1f}  SR={S:7.3f}  phi={math.exp((H-350*S)/(ce.R_J*350)):.4f}")

print("\nExample 2: n-butane saturation at 350 K")
for k, e in ce.EOS.items():
    Ps, L, V = ce.saturation(ce.get_species("n-butane"), e, 350.0)
    print(f"  {k:4s} Psat={Ps:7.4f}  phi_sat={V.phi:.4f}  f={V.f:.4f}  dHvap={(V.H-L.H)/1000:.3f} kJ/mol  dSvap={V.S-L.S:.2f}  Z_V={V.Z:.4f}")
if PropsSI:
    print(f"  ref  Psat={PropsSI('P','T',350,'Q',1,'n-Butane')/1e5:7.4f}  dHvap={(PropsSI('Hmolar','T',350,'Q',1,'n-Butane')-PropsSI('Hmolar','T',350,'Q',0,'n-Butane'))/1000:.3f} kJ/mol")

print("\nExample 3: N2 compressed 300 K,1 bar -> 400 K,100 bar (PR) and 300 K ->300 K isothermal 1->200 bar methane")
dH, dS, s1, s2 = ce.delta_HS(ce.get_species("nitrogen"), ce.EOS["PR"], 300, 1, 400, 100)
print(f"  N2: dH={dH:.1f} J/mol dS={dS:.3f}   (ideal gas: dH={ce.dH_ig(ce.get_species('nitrogen'),300,400):.1f}, dS={ce.dS_ig_T(ce.get_species('nitrogen'),300,400)-ce.R_J*math.log(100):.3f})")
