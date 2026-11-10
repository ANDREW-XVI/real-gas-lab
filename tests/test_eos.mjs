// Compares the JavaScript engine with the Python program (reference.json is produced by python/cubic_eos.py)
// Run:  node tests/test_eos.mjs
import { createRequire } from 'node:module';
import { readFileSync } from 'node:fs';
const require = createRequire(import.meta.url);
const RG = require('../js/eos.js');
const GASES = require('../js/gases.js');
const ref = JSON.parse(readFileSync(new URL('./reference.json', import.meta.url)));
const byId = Object.fromEntries(GASES.map(g => [g.id, g]));

let n = 0, bad = 0, worst = { Z: 0, lnphi: 0, HR: 0, SR: 0, H: 0, S: 0, Psat: 0, dH: 0 };
const rel = (a, b) => Math.abs(a - b) / Math.max(Math.abs(b), 1e-12);
for (const r of ref) {
  const sp = byId[r.gas], eos = RG.EOS[r.eos];
  try {
    if (r.sat) {
      const s = RG.saturation(eos, sp, r.T);
      worst.Psat = Math.max(worst.Psat, rel(s.Psat, r.Psat));
      worst.dH = Math.max(worst.dH, Math.abs((s.vapour.H - s.liquid.H) - r.dH) / Math.abs(r.dH));
      if (rel(s.Psat, r.Psat) > 1e-7) { bad++; console.log('Psat mismatch', r.gas, r.eos, r.T, s.Psat, r.Psat); }
    } else {
      const a = RG.analyze(sp, eos, r.T, r.P), s = a.stable;
      worst.Z = Math.max(worst.Z, rel(s.Z, r.Z));
      worst.lnphi = Math.max(worst.lnphi, Math.abs(s.lnPhi - r.lnphi));
      worst.HR = Math.max(worst.HR, Math.abs(s.HR - r.HR) / (RG.R_J * r.T));
      worst.SR = Math.max(worst.SR, Math.abs(s.SR - r.SR) / RG.R_J);
      worst.H = Math.max(worst.H, rel(s.H, r.H) * Math.min(1, Math.abs(r.H) / 100));
      worst.S = Math.max(worst.S, Math.abs(s.S - r.S) / RG.R_J);
      if (rel(s.Z, r.Z) > 1e-7 || a.roots.length !== r.nroots) { bad++; console.log('mismatch', r.gas, r.eos, r.T, r.P, s.Z, r.Z, a.roots.length, r.nroots); }
    }
  } catch (e) { bad++; console.log('ERROR', r.gas, r.eos, r.T, r.P, e.message); }
  n++;
}
console.log(`${n} cases, ${bad} problems`);
console.log('worst deviations:', Object.fromEntries(Object.entries(worst).map(([k, v]) => [k, v.toExponential(1)])));
process.exit(bad ? 1 : 0);
