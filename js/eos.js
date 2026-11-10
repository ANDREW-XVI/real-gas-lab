/* eos.js - cubic equations of state for pure real gases (vdW, RK, SRK, PR)
 *
 * PVT, fugacity, residual and total enthalpy/entropy, vapour pressure.
 * Same formulation and constants as python/cubic_eos.py (Smith, Van Ness & Abbott, Ch. 3, 6, 10, 13).
 * Units: T [K], P [bar], V [cm3/mol], H [J/mol], S [J/(mol K)], f [bar].
 */
(function (root) {
  'use strict';
  const R_J = 8.314462618;       // J/(mol K)
  const R_V = 83.14462618;       // bar cm3/(mol K)
  const T0 = 298.15, P0 = 1.0;   // reference state of the ideal-gas H and S

  const OM_RK = (Math.cbrt(2) - 1) / 3;
  const PSI_RK = 1 / (9 * (Math.cbrt(2) - 1));
  const EOS = {
    VDW: { key: 'VDW', name: 'van der Waals', short: 'vdW', sigma: 0, eps: 0, Omega: 1 / 8, Psi: 27 / 64 },
    RK:  { key: 'RK',  name: 'Redlich-Kwong', short: 'RK', sigma: 1, eps: 0, Omega: OM_RK, Psi: PSI_RK },
    SRK: { key: 'SRK', name: 'Soave-Redlich-Kwong', short: 'SRK', sigma: 1, eps: 0, Omega: OM_RK, Psi: PSI_RK },
    PR:  { key: 'PR',  name: 'Peng-Robinson', short: 'PR', sigma: 1 + Math.SQRT2, eps: 1 - Math.SQRT2,
           Omega: 0.0777960739038885, Psi: 0.457235528921382 }
  };

  /** [alpha, dln(alpha)/dln(Tr)] */
  function alpha(eos, Tr, w) {
    if (eos.key === 'VDW') return [1, 0];
    if (eos.key === 'RK') return [Math.pow(Tr, -0.5), -0.5];
    const k = eos.key === 'SRK' ? 0.480 + 1.574 * w - 0.176 * w * w
                                : 0.37464 + 1.54226 * w - 0.26992 * w * w;
    const s = 1 + k * (1 - Math.sqrt(Tr));
    return [s * s, -k * Math.sqrt(Tr) / s];
  }

  function AB(eos, sp, T, P) {
    const Tr = T / sp.Tc, Pr = P / sp.Pc, al = alpha(eos, Tr, sp.omega);
    return [eos.Psi * al[0] * Pr / (Tr * Tr), eos.Omega * Pr / Tr, al[1]];
  }
  function ab(eos, sp, T) {
    const al = alpha(eos, T / sp.Tc, sp.omega)[0];
    return [eos.Psi * al * R_V * R_V * sp.Tc * sp.Tc / sp.Pc, eos.Omega * R_V * sp.Tc / sp.Pc];
  }
  function pressure(eos, sp, T, V) {
    const [a, b] = ab(eos, sp, T);
    return R_V * T / (V - b) - a / ((V + eos.eps * b) * (V + eos.sigma * b));
  }

  /** physical (Z > B) real roots of the cubic, ascending */
  function zRoots(eos, A, B) {
    const s = eos.sigma, e = eos.eps;
    const c2 = (e + s) * B - (1 + B);
    const c1 = e * s * B * B - (1 + B) * (e + s) * B + A;
    const c0 = -((1 + B) * e * s * B * B + A * B);
    // depressed cubic t^3 + p t + q = 0, Z = t - c2/3
    const p = c1 - c2 * c2 / 3, q = 2 * c2 * c2 * c2 / 27 - c2 * c1 / 3 + c0;
    const disc = q * q / 4 + p * p * p / 27;
    let ts = [];
    if (disc > 0) {
      const sq = Math.sqrt(disc);
      ts = [Math.cbrt(-q / 2 + sq) + Math.cbrt(-q / 2 - sq)];
    } else {
      const r = Math.sqrt(Math.max(-p * p * p / 27, 0));
      const m = 2 * Math.sqrt(Math.max(-p / 3, 0));
      const ph = Math.acos(Math.max(-1, Math.min(1, r === 0 ? 0 : -q / (2 * r))));
      for (let k = 0; k < 3; k++) ts.push(m * Math.cos((ph + 2 * Math.PI * k) / 3));
    }
    let zs = ts.map(t => t - c2 / 3);
    zs = zs.map(z0 => {                    // Newton polishing
      let z = z0;
      for (let i = 0; i < 40; i++) {
        const f = ((z + c2) * z + c1) * z + c0, fp = (3 * z + 2 * c2) * z + c1;
        if (fp === 0) break;
        const dz = f / fp; z -= dz;
        if (Math.abs(dz) <= 1e-15 * Math.abs(z)) break;
      }
      return Math.abs(z - z0) <= 1e-6 * Math.max(Math.abs(z0), 1e-300) ? z : z0;
    });
    zs = zs.filter(z => z > B * (1 + 1e-12)).sort((x, y) => x - y);
    const out = [];
    for (const z of zs) if (!out.length || Math.abs(z - out[out.length - 1]) > 1e-9 * Math.abs(z)) out.push(z);
    if (!out.length) throw new Error('No physical root (Z > B) for the cubic.');
    return out;
  }

  function Iint(eos, Z, B) {
    return eos.sigma === eos.eps ? B / Z :
      Math.log((Z + eos.sigma * B) / (Z + eos.eps * B)) / (eos.sigma - eos.eps);
  }
  function lnPhi(eos, Z, A, B) { return Z - 1 - Math.log(Z - B) - (A / B) * Iint(eos, Z, B); }

  /* ideal-gas part: Cp/R = A + B T + C T^2 + D T^-2 */
  function dHig(sp, Ta, Tb) {
    const [A, B, C, D] = sp.cp;
    return R_J * (A * (Tb - Ta) + B / 2 * (Tb * Tb - Ta * Ta) + C / 3 * (Tb ** 3 - Ta ** 3) + D * (1 / Ta - 1 / Tb));
  }
  function dSigT(sp, Ta, Tb) {
    const [A, B, C, D] = sp.cp;
    return R_J * (A * Math.log(Tb / Ta) + B * (Tb - Ta) + C / 2 * (Tb * Tb - Ta * Ta) + D / 2 * (1 / (Ta * Ta) - 1 / (Tb * Tb)));
  }
  const Hig = (sp, T) => dHig(sp, T0, T);
  const Sig = (sp, T, P) => dSigT(sp, T0, T) - R_J * Math.log(P / P0);

  function evaluate(sp, eos, T, P, Z, label, stable) {
    const [A, B, dl] = AB(eos, sp, T, P);
    const I = Iint(eos, Z, B), q = A / B;
    const lnp = Z - 1 - Math.log(Z - B) - q * I;
    const HR = R_J * T * (Z - 1 + (dl - 1) * q * I);
    const SR = R_J * (Math.log(Z - B) + dl * q * I);
    return { label: label || '', T, P, Z, V: Z * R_V * T / P, lnPhi: lnp, phi: Math.exp(lnp), f: P * Math.exp(lnp),
             HR, SR, GR: R_J * T * lnp, H: Hig(sp, T) + HR, S: Sig(sp, T, P) + SR, stable: stable !== false };
  }

  /** extrema of the isotherm: {Pl, Ph, Va, Vb}; throws if there is no loop */
  function spinodal(eos, sp, T) {
    const [a, b] = ab(eos, sp, T), s = eos.sigma, e = eos.eps;
    const f = V => -R_V * T / ((V - b) * (V - b)) + a * (2 * V + (e + s) * b) / (Math.pow(V + e * b, 2) * Math.pow(V + s * b, 2));
    const N = 6000, lo = Math.log(b * (1 + 1e-6)), hi = Math.log(b * 1e5), roots = [];
    let v0 = Math.exp(lo), f0 = f(v0);
    for (let i = 1; i <= N; i++) {
      const v1 = Math.exp(lo + (hi - lo) * i / N), f1 = f(v1);
      if (f0 * f1 < 0) {
        let x = v0, y = v1, fx = f0;
        for (let k = 0; k < 200; k++) {
          const m = Math.sqrt(x * y), fm = f(m);
          if (fm * fx > 0) { x = m; fx = fm; } else { y = m; }
          if (Math.abs(y - x) < 1e-14 * x) break;
        }
        roots.push(Math.sqrt(x * y));
      }
      v0 = v1; f0 = f1;
    }
    if (roots.length < 2) throw new Error('No van der Waals loop at this temperature.');
    const Va = roots[0], Vb = roots[roots.length - 1];
    return { Pl: pressure(eos, sp, T, Va), Ph: pressure(eos, sp, T, Vb), Va, Vb };
  }

  /** vapour pressure from equal fugacity: {Psat, liquid, vapour} */
  function saturation(eos, sp, T) {
    if (T >= sp.Tc) throw new Error('T must be below Tc = ' + sp.Tc + ' K for vapour-liquid equilibrium.');
    const sg = spinodal(eos, sp, T), Vmid = 0.5 * (sg.Va + sg.Vb);
    const lo = Math.max(sg.Pl, 1e-12 * sp.Pc);
    const g = lnP => {
      const P = Math.exp(lnP), [A, B] = AB(eos, sp, T, P), zs = zRoots(eos, A, B);
      if (zs.length < 2) return zs[0] * R_V * T / P < Vmid ? -1 : 1;
      return lnPhi(eos, zs[0], A, B) - lnPhi(eos, zs[zs.length - 1], A, B);
    };
    let x = Math.log(lo * (1 + 1e-9)), y = Math.log(sg.Ph * (1 - 1e-9));
    let gx = g(x);
    if (gx * g(y) > 0) throw new Error('Could not bracket the saturation pressure.');
    for (let k = 0; k < 200; k++) {
      const m = 0.5 * (x + y), gm = g(m);
      if (gm * gx > 0) { x = m; gx = gm; } else { y = m; }
      if (Math.abs(y - x) < 1e-14) break;
    }
    const Ps = Math.exp(0.5 * (x + y)), [A, B] = AB(eos, sp, T, Ps), zs = zRoots(eos, A, B);
    return { Psat: Ps, liquid: evaluate(sp, eos, T, Ps, zs[0], 'saturated liquid'),
             vapour: evaluate(sp, eos, T, Ps, zs[zs.length - 1], 'saturated vapour'), spinodal: sg };
  }

  /** full (T, P) analysis */
  function analyze(sp, eos, T, P) {
    if (!(T > 0) || !(P > 0)) throw new Error('T and P must be positive.');
    const [A, B] = AB(eos, sp, T, P), zs = zRoots(eos, A, B);
    const roots = zs.map(z => evaluate(sp, eos, T, P, z));
    const stable = roots.reduce((m, p) => (p.lnPhi < m.lnPhi ? p : m), roots[0]);
    let Psat = null;
    if (T < sp.Tc) { try { Psat = saturation(eos, sp, T).Psat; } catch (e) { Psat = null; } }
    let phase;
    if (T >= sp.Tc) phase = P >= sp.Pc ? 'supercritical fluid' : 'gas (T > Tc)';
    else if (Psat !== null && Math.abs(P - Psat) <= 1e-6 * Psat) phase = 'saturated vapour-liquid equilibrium';
    else if (Psat !== null) phase = P < Psat ? 'superheated vapour' : 'compressed liquid';
    else phase = stable.Z > 0.3 ? 'vapour' : 'liquid';
    roots.forEach(p => {
      p.stable = p === stable;
      p.label = roots.length === 3 ? (p === roots[0] ? 'liquid-like root' : p === roots[2] ? 'vapour-like root' : 'middle root (unstable)') : 'only real root';
    });
    const [a, b] = ab(eos, sp, T);
    return { sp, eos, T, P, phase, stable, roots, Psat, a, b, A, B };
  }

  /** property change between two states */
  function deltaHS(sp, eos, T1, P1, T2, P2) {
    const r1 = analyze(sp, eos, T1, P1), r2 = analyze(sp, eos, T2, P2), s1 = r1.stable, s2 = r2.stable;
    return { dH: s2.H - s1.H, dS: s2.S - s1.S, dH_ig: dHig(sp, T1, T2), dS_ig: dSigT(sp, T1, T2) - R_J * Math.log(P2 / P1),
             s1, s2 };
  }

  const api = { R_J, R_V, T0, P0, EOS, alpha, AB, ab, pressure, zRoots, lnPhi, evaluate, spinodal, saturation, analyze,
                deltaHS, dHig, dSigT, Hig, Sig };
  root.RealGas = api;
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
})(typeof window !== 'undefined' ? window : globalThis);
