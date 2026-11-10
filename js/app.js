/* app.js - user interface of the Real-gas lab (uses RealGas from eos.js and GASES from gases.js) */
(function () {
  'use strict';
  const RG = window.RealGas, EOS = RG.EOS;
  const $ = (s, r) => (r || document).querySelector(s);
  const $$ = (s, r) => Array.from((r || document).querySelectorAll(s));

  /* ---------- state ---------- */
  const S = { gas: 'carbon-dioxide', eos: 'PR', T: 350, P: 100, tab: 'diagram', metric: 'Z', T2: null, P2: null,
              custom: { Tc: 190.6, Pc: 45.99, w: 0.012, cp: [3.5, 0, 0, 0] } };
  const gasList = window.GASES.slice();
  const species = () => S.gas === 'custom'
    ? { id: 'custom', name: 'Custom gas', formula: 'custom', Tc: S.custom.Tc, Pc: S.custom.Pc, omega: S.custom.w, cp: S.custom.cp }
    : gasList.find(g => g.id === S.gas);

  /* ---------- helpers ---------- */
  const MINUS = '−';
  function nf(x, d) {
    if (!isFinite(x)) return '–';
    const s = Math.abs(x) >= 1e5 ? x.toExponential(3) : x.toLocaleString('en-US', { minimumFractionDigits: d, maximumFractionDigits: d });
    return s.replace('-', MINUS);
  }
  const sig = (x, d) => nf(x, d);
  function parseT(x) {
    const s = String(x).trim().toLowerCase().replace(',', '.');
    if (/k$/.test(s)) return parseFloat(s);
    if (/c$/.test(s)) return parseFloat(s) + 273.15;
    if (/f$/.test(s)) return (parseFloat(s) - 32) * 5 / 9 + 273.15;
    return parseFloat(s);
  }
  function parseP(x) {
    const s = String(x).trim().toLowerCase().replace(',', '.');
    const units = [['mpa', 10], ['kpa', 0.01], ['bar', 1], ['atm', 1.01325], ['psi', 0.0689475729], ['pa', 1e-5]];
    for (const [u, f] of units) if (s.endsWith(u)) return parseFloat(s) * f;
    return parseFloat(s);
  }
  const esc = s => String(s).replace(/[&<>]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;' }[c]));
  const formulaHTML = f => esc(f).replace(/(\d+)/g, '<sub>$1</sub>');

  /* ---------- url hash ---------- */
  function readHash() {
    const h = new URLSearchParams(location.hash.slice(1));
    if (h.get('gas') && (h.get('gas') === 'custom' || gasList.some(g => g.id === h.get('gas')))) S.gas = h.get('gas');
    if (EOS[h.get('eos')]) S.eos = h.get('eos');
    if (isFinite(parseFloat(h.get('T'))) && parseFloat(h.get('T')) > 0) S.T = parseFloat(h.get('T'));
    if (isFinite(parseFloat(h.get('P'))) && parseFloat(h.get('P')) > 0) S.P = parseFloat(h.get('P'));
    if (['diagram', 'properties', 'compare', 'theory'].includes(h.get('tab'))) S.tab = h.get('tab');
    if (S.gas === 'custom') {
      const n = k => parseFloat(h.get(k));
      if (n('Tc') > 0 && n('Pc') > 0 && isFinite(n('w'))) S.custom = { Tc: n('Tc'), Pc: n('Pc'), w: n('w'), cp: (h.get('cp') || '3.5 0 0 0').split('_').map(Number) };
    }
  }
  function writeHash() {
    const h = new URLSearchParams({ gas: S.gas, eos: S.eos, T: +S.T.toPrecision(8), P: +S.P.toPrecision(8), tab: S.tab });
    if (S.gas === 'custom') { h.set('Tc', S.custom.Tc); h.set('Pc', S.custom.Pc); h.set('w', S.custom.w); h.set('cp', S.custom.cp.join('_')); }
    history.replaceState(null, '', '#' + h.toString());
  }

  /* ---------- sliders ---------- */
  const Trange = sp => [0.4 * sp.Tc, 3 * sp.Tc];
  const Prange = sp => [0.1, 5 * sp.Pc];
  function setSliders() {
    const sp = species(), [t0, t1] = Trange(sp), [p0, p1] = Prange(sp);
    $('#Tr').value = Math.max(0, Math.min(1000, (S.T - t0) / (t1 - t0) * 1000));
    $('#Pr').value = Math.max(0, Math.min(1000, Math.log(S.P / p0) / Math.log(p1 / p0) * 1000));
  }

  /* ---------- computation ---------- */
  let res = null;
  function compute() {
    const sp = species();
    res = RG.analyze(sp, EOS[S.eos], S.T, S.P);
  }

  /* ---------- readout ---------- */
  function readout() {
    const s = res.stable, el = $('#readout');
    const warn = /liquid/.test(res.phase) && !/saturated/.test(res.phase) ? '' : '';
    el.innerHTML =
      `<div><span class="pill${warn}">${esc(res.phase)}</span></div>
       <div class="rd"><div class="l">Compressibility Z</div><div class="v">${nf(s.Z, 4)}</div></div>
       <div class="rd"><div class="l">Molar volume</div><div class="v">${nf(s.V, s.V < 1000 ? 1 : 0)}<small>cm³/mol</small></div></div>
       <div class="rd"><div class="l">Fugacity</div><div class="v">${nf(s.f, s.f < 1000 ? 2 : 0)}<small>bar</small></div><div class="s">φ = ${nf(s.phi, 4)} · ln φ = ${nf(s.lnPhi, 4)}</div></div>
       <div class="rd"><div class="l">Residual enthalpy H<sup>R</sup></div><div class="v">${nf(s.HR, 0)}<small>J/mol</small></div></div>
       <div class="rd"><div class="l">Residual entropy S<sup>R</sup></div><div class="v">${nf(s.SR, 2)}<small>J/mol·K</small></div></div>
       <div class="rd sm"><div class="l">Total H · S (ref. ideal gas, 298.15 K, 1 bar)</div><div class="v">${nf(s.H, 0)}<small>J/mol</small></div><div class="v">${nf(s.S, 2)}<small>J/mol·K</small></div></div>`;
    $('#bsub').textContent = `${species().name} · ${EOS[S.eos].name}`;
  }

  /* ---------- svg helpers ---------- */
  function niceStep(span, n) {
    const raw = span / n, p = Math.pow(10, Math.floor(Math.log10(raw))), m = raw / p;
    return (m < 1.5 ? 1 : m < 3.5 ? 2 : m < 7.5 ? 5 : 10) * p;
  }
  function logTicks(a, b) {
    const out = [];
    for (let k = Math.floor(Math.log10(a)); k <= Math.ceil(Math.log10(b)); k++)
      for (const m of [1, 2, 5]) { const v = m * Math.pow(10, k); if (v >= a && v <= b) out.push(v); }
    return out;
  }
  function pathSegments(pts, X, Y, ymax) {
    let d = '', pen = false;
    for (const [v, p] of pts) {
      if (!(p >= 0 && p <= ymax) || !isFinite(p)) { pen = false; continue; }
      d += (pen ? 'L' : 'M') + X(v).toFixed(1) + ' ' + Y(p).toFixed(1); pen = true;
    }
    return d;
  }

  /* ---------- dome cache ---------- */
  const domeCache = new Map();
  function dome(sp, eos) {
    const key = [sp.Tc, sp.Pc, sp.omega, eos.key].join('|');
    if (domeCache.has(key)) return domeCache.get(key);
    const L = [], V = [];
    for (let i = 0; i <= 80; i++) {
      const T = sp.Tc * (0.55 + 0.4495 * Math.pow(i / 80, 0.7));
      try { const s = RG.saturation(eos, sp, T); L.push([s.liquid.V, s.Psat]); V.push([s.vapour.V, s.Psat]); } catch (e) { /* skip */ }
    }
    const out = { L, V };
    domeCache.set(key, out);
    return out;
  }

  /* ---------- P-V diagram ---------- */
  function drawPV() {
    const host = $('#plot'); if (!host) return;
    const W = Math.max(host.clientWidth, 320), H = Math.max(host.clientHeight, 300);
    const sp = species(), eos = EOS[S.eos], st = res.stable;
    const m = { l: 56, r: 18, t: 26, b: 44 };
    const b = res.b, xmin = b * 1.25, xmax = Math.max(b * 120, st.V * 3), ymax = Math.max(1.6 * sp.Pc, S.P * 1.2);
    const X = v => m.l + (Math.log(Math.min(Math.max(v, xmin), xmax)) - Math.log(xmin)) / (Math.log(xmax) - Math.log(xmin)) * (W - m.l - m.r);
    const Y = p => m.t + (1 - p / ymax) * (H - m.t - m.b);
    const inv = (px, py) => ({ V: Math.exp(Math.log(xmin) + (px - m.l) / (W - m.l - m.r) * (Math.log(xmax) - Math.log(xmin))), P: (1 - (py - m.t) / (H - m.t - m.b)) * ymax });
    let s = `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="Pressure versus molar volume with isotherms and saturation dome">`;
    const ys = niceStep(ymax, 5);
    for (let p = 0; p <= ymax; p += ys) {
      s += `<line x1="${m.l}" x2="${W - m.r}" y1="${Y(p)}" y2="${Y(p)}" stroke="#222a2e"/><text x="${m.l - 10}" y="${Y(p) + 4}" text-anchor="end" font-size="12" fill="#7f8c92">${+p.toPrecision(6)}</text>`;
    }
    for (const v of logTicks(xmin, xmax)) {
      s += `<line y1="${m.t}" y2="${H - m.b}" x1="${X(v)}" x2="${X(v)}" stroke="#1a2124"/><text x="${X(v)}" y="${H - m.b + 18}" text-anchor="middle" font-size="12" fill="#7f8c92">${v >= 1000 ? (v / 1000) + 'k' : v}</text>`;
    }
    s += `<text x="${(m.l + W - m.r) / 2}" y="${H - 6}" text-anchor="middle" font-size="12" fill="#7f8c92">Molar volume (cm³/mol, log scale)</text>`;
    s += `<text x="${m.l - 44}" y="12" font-size="12" fill="#7f8c92">P (bar)</text>`;
    // dome
    const d = dome(sp, eos);
    if (d.L.length > 3) {
      const pts = d.L.concat(d.V.slice().reverse()).map(([v, p]) => `${X(v).toFixed(1)},${Y(Math.min(p, ymax)).toFixed(1)}`).join(' ');
      s += `<polygon points="${pts}" fill="#16302f" stroke="#2f6f6b" stroke-width="1.5"/>`;
      const mid = d.L[Math.floor(d.L.length * 0.55)];
      if (mid) s += `<text x="${(X(mid[0]) + X(d.V[Math.floor(d.V.length * 0.55)][0])) / 2}" y="${Y(mid[1] * 0.55)}" text-anchor="middle" font-size="13" fill="#4f9d97">liquid + vapour</text>`;
    }
    // isotherms
    const temps = [[0.85 * sp.Tc, 'ref'], [0.95 * sp.Tc, 'ref'], [sp.Tc, 'tc'], [1.15 * sp.Tc, 'ref']].filter(([t]) => Math.abs(t - S.T) / S.T > 0.012);
    temps.push([S.T, 'sel']);
    for (const [T, kind] of temps) {
      const [, bb] = RG.ab(eos, sp, T), pts = [];
      for (let i = 0; i <= 420; i++) { const v = bb * 1.02 * Math.pow(xmax / (bb * 1.02), i / 420); pts.push([v, RG.pressure(eos, sp, T, v)]); }
      const col = kind === 'sel' ? '#3fd0b4' : kind === 'tc' ? '#8ea1a8' : '#3b4a50';
      s += `<path d="${pathSegments(pts, X, Y, ymax)}" fill="none" stroke="${col}" stroke-width="${kind === 'sel' ? 3 : 1.5}" ${kind === 'tc' ? 'stroke-dasharray="5 4"' : ''} stroke-linejoin="round"/>`;
      const first = pts.find(([, p]) => p >= 0 && p <= ymax);
      if (first && (kind !== 'ref' || true)) {
        const lab = kind === 'tc' ? 'Tc' : Math.round(T) + ' K';
        s += `<text x="${Math.min(X(first[0]) + 6, W - m.r - 40)}" y="${Y(first[1]) + 14}" font-size="12" fill="${kind === 'sel' ? '#3fd0b4' : '#7f8c92'}" font-weight="${kind === 'sel' ? 600 : 400}">${lab}</text>`;
      }
    }
    // tie line at Psat
    if (res.Psat) {
      try {
        const sat = RG.saturation(eos, sp, S.T);
        s += `<line x1="${X(sat.liquid.V)}" x2="${X(sat.vapour.V)}" y1="${Y(sat.Psat)}" y2="${Y(sat.Psat)}" stroke="#e7b04a" stroke-width="1.5" stroke-dasharray="4 3"/>`;
        s += `<circle cx="${X(sat.liquid.V)}" cy="${Y(sat.Psat)}" r="3.5" fill="#e7b04a"/><circle cx="${X(sat.vapour.V)}" cy="${Y(sat.Psat)}" r="3.5" fill="#e7b04a"/>`;
      } catch (e) { /* ignore */ }
    }
    // state
    s += `<line x1="${X(st.V)}" x2="${X(st.V)}" y1="${Y(S.P)}" y2="${H - m.b}" stroke="#3fd0b4" stroke-dasharray="3 4"/>`;
    s += `<line x1="${m.l}" x2="${X(st.V)}" y1="${Y(S.P)}" y2="${Y(S.P)}" stroke="#3fd0b4" stroke-dasharray="3 4"/>`;
    s += `<circle cx="${X(st.V)}" cy="${Y(S.P)}" r="7" fill="#0e1214" stroke="#3fd0b4" stroke-width="3"/>`;
    s += `<rect class="hit" x="${m.l}" y="${m.t}" width="${W - m.l - m.r}" height="${H - m.t - m.b}" fill="transparent" style="cursor:crosshair"/></svg><div class="tip" id="tip"></div>`;
    host.innerHTML = s;
    const svg = $('svg', host), tip = $('#tip');
    const pos = e => { const r = svg.getBoundingClientRect(); return { x: (e.clientX - r.left) * W / r.width, y: (e.clientY - r.top) * H / r.height, cx: e.clientX - r.left, cy: e.clientY - r.top }; };
    svg.addEventListener('pointermove', e => {
      const p = pos(e); if (p.x < m.l || p.x > W - m.r || p.y < m.t || p.y > H - m.b) { tip.style.display = 'none'; return; }
      const q = inv(p.x, p.y);
      tip.textContent = `V ${nf(q.V, 1)} cm³/mol · P ${nf(q.P, 1)} bar (click to set P)`;
      tip.style.display = 'block'; tip.style.left = Math.min(p.cx + 14, host.clientWidth - tip.offsetWidth - 4) + 'px'; tip.style.top = Math.max(p.cy - 34, 4) + 'px';
    });
    svg.addEventListener('pointerleave', () => { tip.style.display = 'none'; });
    svg.addEventListener('click', e => {
      const p = pos(e); if (p.y < m.t || p.y > H - m.b) return;
      const q = inv(p.x, p.y); if (q.P > 0.01) { S.P = +q.P.toPrecision(5); update(); }
    });
  }

  /* ---------- line chart for the compare view ---------- */
  const EOS_STYLE = { VDW: ['#3987e5', '9 4 2 4'], RK: ['#d95926', '2 3'], SRK: ['#c98500', '6 4'], PR: ['#3fd0b4', ''] };
  const METRICS = { Z: ['Z', s => s.Z, 'Compressibility factor Z'], phi: ['φ', s => s.phi, 'Fugacity coefficient φ'],
                    HR: ['HR', s => s.HR / 1000, 'Residual enthalpy H^R (kJ/mol)'], SR: ['SR', s => s.SR, 'Residual entropy S^R (J/mol·K)'] };
  function drawCompare() {
    const host = $('#cmpPlot'); if (!host) return;
    const W = Math.max(host.clientWidth, 320), H = Math.max(host.clientHeight, 260), sp = species();
    const m = { l: 60, r: 44, t: 12, b: 40 }, pmax = 4 * sp.Pc, N = 90;
    const Ps = Array.from({ length: N }, (_, i) => Math.max(0.05, pmax * (i + 0.5) / N));
    const fn = METRICS[S.metric][1], series = {};
    let lo = Infinity, hi = -Infinity;
    for (const k of Object.keys(EOS)) {
      series[k] = Ps.map(p => { try { return fn(RG.analyze(sp, EOS[k], S.T, p).stable); } catch (e) { return NaN; } });
      for (const v of series[k]) if (isFinite(v)) { lo = Math.min(lo, v); hi = Math.max(hi, v); }
    }
    const here = fn(res.stable); lo = Math.min(lo, here); hi = Math.max(hi, here);
    const pad = (hi - lo) * 0.06; lo -= pad; hi += pad;
    const X = p => m.l + p / pmax * (W - m.l - m.r), Y = v => m.t + (1 - (v - lo) / (hi - lo)) * (H - m.t - m.b);
    let s = `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="${METRICS[S.metric][2]} versus pressure">`;
    const st = niceStep(hi - lo, 5);
    for (let v = Math.ceil(lo / st) * st; v <= hi; v += st) s += `<line x1="${m.l}" x2="${W - m.r}" y1="${Y(v)}" y2="${Y(v)}" stroke="#222a2e"/><text x="${m.l - 10}" y="${Y(v) + 4}" text-anchor="end" font-size="12" fill="#7f8c92">${(+v.toPrecision(6)).toString().replace('-', MINUS)}</text>`;
    const xs = niceStep(pmax, 5);
    for (let p = 0; p <= pmax; p += xs) s += `<text x="${X(p)}" y="${H - 20}" text-anchor="middle" font-size="12" fill="#7f8c92">${+p.toPrecision(6)}</text>`;
    s += `<text x="${(m.l + W - m.r) / 2}" y="${H - 4}" text-anchor="middle" font-size="12" fill="#7f8c92">Pressure (bar) at ${nf(S.T, 1)} K</text>`;
    const ends = Object.keys(EOS).map(k => ({ k, y: Y(series[k].filter(isFinite).pop() ?? lo) })).sort((a, b) => a.y - b.y);
    for (let i = 1; i < ends.length; i++) if (ends[i].y - ends[i - 1].y < 14) ends[i].y = ends[i - 1].y + 14;
    const labY = Object.fromEntries(ends.map(e => [e.k, e.y]));
    for (const k of ['VDW', 'RK', 'SRK', 'PR']) {
      const [col, dash] = EOS_STYLE[k], sel = k === S.eos;
      let d = '', pen = false;
      series[k].forEach((v, i) => { if (!isFinite(v)) { pen = false; return; } d += (pen ? 'L' : 'M') + X(Ps[i]).toFixed(1) + ' ' + Y(v).toFixed(1); pen = true; });
      s += `<path d="${d}" fill="none" stroke="${col}" stroke-width="${sel ? 3 : 1.8}" stroke-dasharray="${dash}" stroke-linejoin="round"/>`;
      s += `<text x="${W - m.r + 6}" y="${labY[k] + 4}" font-size="12" fill="${col}" font-weight="${sel ? 650 : 500}">${EOS[k].short}</text>`;
    }
    s += `<line x1="${X(S.P)}" x2="${X(S.P)}" y1="${m.t}" y2="${H - m.b}" stroke="#7f8c92" stroke-dasharray="3 4"/><circle cx="${X(Math.min(S.P, pmax))}" cy="${Y(here)}" r="5.5" fill="#0e1214" stroke="${EOS_STYLE[S.eos][0]}" stroke-width="3"/></svg>`;
    host.innerHTML = s;
  }

  /* ---------- views ---------- */
  function viewDiagram() {
    const sp = species();
    return `<div class="shd"><h1>P–V diagram, ${formulaHTML(sp.formula === 'custom' ? sp.name : sp.formula)}</h1><span>${EOS[S.eos].name}: isotherms and saturation dome. Click the diagram to set the pressure.</span></div>
      <div class="plot" id="plot"></div>
      <div class="legendline"><span><i style="border-color:#3fd0b4"></i>selected T</span><span><i style="border-color:#8ea1a8;border-top-style:dashed"></i>critical isotherm</span><span><i style="border-color:#3b4a50"></i>other isotherms</span><span><i style="border-color:#2f6f6b;background:#16302f;height:8px"></i>two-phase dome</span><span><i style="border-color:#e7b04a;border-top-style:dashed"></i>vapour pressure (tie line)</span></div>`;
  }

  function viewProperties() {
    const s = res.stable, sp = species(), eos = EOS[S.eos];
    const rows = res.roots.map(p => `<tr class="${p.stable && res.roots.length > 1 ? 'sel' : (res.roots.length > 1 && !p.stable ? 'dim' : '')}"><td>${esc(p.label)}${p.stable && res.roots.length > 1 ? ' · stable' : ''}</td><td>${nf(p.Z, 5)}</td><td>${nf(p.V, 2)}</td><td>${nf(p.lnPhi, 5)}</td><td>${nf(p.phi, 5)}</td><td>${nf(p.f, 3)}</td><td>${nf(p.HR, 1)}</td><td>${nf(p.SR, 3)}</td></tr>`).join('');
    let vle = '';
    if (S.T < sp.Tc) {
      try {
        const v = RG.saturation(eos, sp, S.T), dH = v.vapour.H - v.liquid.H, dS = v.vapour.S - v.liquid.S;
        vle = `<div class="card"><h2>Vapour–liquid equilibrium at ${nf(S.T, 2)} K</h2><div class="sub">Pressure at which the liquid and vapour fugacities are equal</div>
          <dl class="kv"><dt>Vapour pressure P<sup>sat</sup></dt><dd>${nf(v.Psat, 4)} bar</dd><dt>Fugacity f<sup>sat</sup> (φ<sup>sat</sup>)</dt><dd>${nf(v.vapour.f, 4)} bar (${nf(v.vapour.phi, 4)})</dd>
          <dt>Z, V of saturated liquid</dt><dd>${nf(v.liquid.Z, 5)}, ${nf(v.liquid.V, 2)} cm³/mol</dd><dt>Z, V of saturated vapour</dt><dd>${nf(v.vapour.Z, 5)}, ${nf(v.vapour.V, 1)} cm³/mol</dd>
          <dt>Heat of vaporisation ΔH<sup>vap</sup></dt><dd>${nf(dH, 0)} J/mol</dd><dt>ΔS<sup>vap</sup> = ΔH<sup>vap</sup>/T</dt><dd>${nf(dS, 3)} J/mol·K</dd></dl></div>`;
      } catch (e) { vle = `<div class="card"><h2>Vapour–liquid equilibrium</h2><div class="sub">Not available at this temperature (${esc(e.message)})</div></div>`; }
    } else {
      vle = `<div class="card"><h2>Vapour–liquid equilibrium</h2><div class="sub">T = ${nf(S.T, 1)} K is above Tc = ${nf(sp.Tc, 1)} K, so only one fluid phase exists. Lower T below Tc to see the vapour pressure.</div></div>`;
    }
    if (S.T2 === null) { S.T2 = +(S.T + 50).toFixed(2); S.P2 = +(S.P / 2).toPrecision(4); }
    return `<div class="shd"><h1>${esc(sp.name)} at ${nf(S.T, 2)} K and ${nf(S.P, 3)} bar</h1><span>${EOS[S.eos].name} · a = ${(res.a).toPrecision(5)} bar cm⁶/mol², b = ${res.b.toFixed(3)} cm³/mol · A = ${res.A.toPrecision(5)}, B = ${res.B.toPrecision(5)}</span></div>
      <div class="card"><h2>Roots of the cubic in Z</h2><div class="sub">When three roots exist, the stable phase is the one with the lowest ln φ (lowest Gibbs energy)</div>
        <div class="tw"><table><tr><th>Root</th><th>Z</th><th>V (cm³/mol)</th><th>ln φ</th><th>φ</th><th>f (bar)</th><th>H<sup>R</sup> (J/mol)</th><th>S<sup>R</sup> (J/mol·K)</th></tr>${rows}</table></div></div>
      <div class="two-col"><div class="card"><h2>Total enthalpy and entropy</h2><div class="sub">Ideal-gas part from Cp<sup>ig</sup>(T) plus the residual part</div>
        <dl class="kv"><dt>H<sup>ig</sup>(T)</dt><dd>${nf(RG.Hig(sp, S.T), 1)} J/mol</dd><dt>H<sup>R</sup></dt><dd>${nf(s.HR, 1)} J/mol</dd><dt><b>H = H<sup>ig</sup> + H<sup>R</sup></b></dt><dd><b>${nf(s.H, 1)} J/mol</b></dd>
        <dt>S<sup>ig</sup>(T, P)</dt><dd>${nf(RG.Sig(sp, S.T, S.P), 3)} J/mol·K</dd><dt>S<sup>R</sup></dt><dd>${nf(s.SR, 3)} J/mol·K</dd><dt><b>S = S<sup>ig</sup> + S<sup>R</sup></b></dt><dd><b>${nf(s.S, 3)} J/mol·K</b></dd>
        <dt>G<sup>R</sup> = RT ln φ</dt><dd>${nf(s.GR, 1)} J/mol</dd></dl>
        <p class="hint">Reference: H = S = 0 for the ideal gas at 298.15 K and 1 bar. Only changes in H and S are physically meaningful.</p></div>
      ${vle}</div>
      <div class="card"><h2>Change in H and S to another state</h2><div class="sub">From the state in the inputs to (T₂, P₂), using the stable phase at each end</div>
        <div class="inl"><div><label for="T2">T₂ (K)</label><input id="T2" inputmode="decimal" value="${S.T2}"></div><div><label for="P2">P₂ (bar)</label><input id="P2" inputmode="decimal" value="${S.P2}"></div><button class="btn" id="go2" type="button">Calculate</button></div>
        <div id="delta"></div></div>`;
  }
  function updateDelta() {
    const host = $('#delta'); if (!host) return;
    const t2 = parseT($('#T2').value), p2 = parseP($('#P2').value);
    if (!(t2 > 0 && p2 > 0)) { host.innerHTML = '<p class="hint">Enter positive values for T₂ and P₂.</p>'; return; }
    S.T2 = t2; S.P2 = p2;
    try {
      const d = RG.deltaHS(species(), EOS[S.eos], S.T, S.P, t2, p2);
      host.innerHTML = `<dl class="kv"><dt>State 2</dt><dd>Z = ${nf(d.s2.Z, 4)}, H<sup>R</sup> = ${nf(d.s2.HR, 0)} J/mol, S<sup>R</sup> = ${nf(d.s2.SR, 2)} J/mol·K</dd>
        <dt><b>ΔH</b> (ideal-gas part ${nf(d.dH_ig, 0)})</dt><dd><b>${nf(d.dH, 1)} J/mol</b></dd><dt><b>ΔS</b> (ideal-gas part ${nf(d.dS_ig, 3)})</dt><dd><b>${nf(d.dS, 3)} J/mol·K</b></dd></dl>`;
    } catch (e) { host.innerHTML = `<p class="err">${esc(e.message)}</p>`; }
  }

  function viewCompare() {
    const sp = species();
    const rows = Object.keys(EOS).map(k => {
      let r; try { r = RG.analyze(sp, EOS[k], S.T, S.P); } catch (e) { return `<tr><td>${EOS[k].name}</td><td colspan="8">${esc(e.message)}</td></tr>`; }
      const s = r.stable;
      return `<tr class="${k === S.eos ? 'sel' : ''}"><td>${EOS[k].name}</td><td>${nf(s.Z, 4)}</td><td>${nf(s.V, 1)}</td><td>${nf(s.phi, 4)}</td><td>${nf(s.f, 2)}</td><td>${nf(s.HR, 0)}</td><td>${nf(s.SR, 2)}</td><td>${nf(s.H, 0)}</td><td>${nf(s.S, 2)}</td><td>${r.Psat ? nf(r.Psat, 3) : '–'}</td></tr>`;
    }).join('');
    const seg = Object.entries(METRICS).map(([k, v]) => `<button type="button" data-metric="${k}" class="${k === S.metric ? 'on' : ''}">${k === 'phi' ? 'φ' : k === 'HR' ? 'H<sup>R</sup>' : k === 'SR' ? 'S<sup>R</sup>' : 'Z'}</button>`).join('');
    return `<div class="shd"><h1>Four equations, one state</h1><span>${esc(sp.name)} at ${nf(S.T, 2)} K and ${nf(S.P, 3)} bar</span></div>
      <div class="card"><div class="tw"><table><tr><th>EOS</th><th>Z</th><th>V (cm³/mol)</th><th>φ</th><th>f (bar)</th><th>H<sup>R</sup> (J/mol)</th><th>S<sup>R</sup> (J/mol·K)</th><th>H (J/mol)</th><th>S (J/mol·K)</th><th>P<sup>sat</sup> (bar)</th></tr>${rows}</table></div></div>
      <div class="shd" style="margin-top:4px"><h1 style="font-size:16px">${METRICS[S.metric][2]} versus pressure</h1><div class="seg" id="seg">${seg}</div></div>
      <div class="plot" id="cmpPlot" style="min-height:300px"></div>
      <div class="legendline">${Object.keys(EOS).map(k => `<span><i style="border-color:${EOS_STYLE[k][0]};${EOS_STYLE[k][1] ? 'border-top-style:dashed' : ''}"></i>${EOS[k].name}</span>`).join('')}<span>the selected equation is drawn thicker; a jump marks condensation at the vapour pressure</span></div>`;
  }

  /* Typeset TeX with the bundled KaTeX; fall back to plain text if it failed to load. */
  function tex(src, display) {
    if (window.katex) {
      try { return window.katex.renderToString(src, { displayMode: !!display, throwOnError: false, output: 'html' }); }
      catch (e) { /* fall through */ }
    }
    return '<code>' + src.replace(/&/g, '&amp;').replace(/</g, '&lt;') + '</code>';
  }
  function eq(src, label) {
    return `<div class="eq">${tex(src, true)}${label ? `<span class="eqn">${label}</span>` : ''}</div>`;
  }

  function viewTheory() {
    const t = (s) => tex(s, false);
    return `<div class="prose"><div class="shd"><h1>Theory and method</h1></div>
    <p>All four equations of state are written in one generic form, so a single set of formulae covers van der Waals, Redlich-Kwong, Soave-Redlich-Kwong and Peng-Robinson. Symbols follow Smith, Van Ness &amp; Abbott, <i>Introduction to Chemical Engineering Thermodynamics</i>.</p>
    <h2>Generic cubic equation</h2>
    ${eq('P = \\frac{RT}{V-b} - \\frac{a(T)}{(V+\\epsilon b)(V+\\sigma b)}', '1')}
    <p>with ${t('b = \\Omega\\,\\dfrac{R T_c}{P_c}')} and ${t('a(T) = \\Psi\\,\\alpha(T_r)\\,\\dfrac{R^2 T_c^{2}}{P_c}')}, where ${t('T_r = T/T_c')}. The constants are:</p>
    <div class="card"><div class="tw"><table class="tx"><tr><th>Equation</th><th>${t('\\sigma')}</th><th>${t('\\epsilon')}</th><th>${t('\\Omega')}</th><th>${t('\\Psi')}</th><th>${t('\\alpha(T_r)')}</th></tr>
    <tr><td>van der Waals</td><td>${t('0')}</td><td>${t('0')}</td><td>${t('1/8')}</td><td>${t('27/64')}</td><td>${t('1')}</td></tr>
    <tr><td>Redlich-Kwong</td><td>${t('1')}</td><td>${t('0')}</td><td>${t('0.08664')}</td><td>${t('0.42748')}</td><td>${t('T_r^{-1/2}')}</td></tr>
    <tr><td>Soave-Redlich-Kwong</td><td>${t('1')}</td><td>${t('0')}</td><td>${t('0.08664')}</td><td>${t('0.42748')}</td><td>${t('\\left[1+\\kappa_{\\mathrm{SRK}}\\left(1-\\sqrt{T_r}\\right)\\right]^2')}</td></tr>
    <tr><td>Peng-Robinson</td><td>${t('1+\\sqrt{2}')}</td><td>${t('1-\\sqrt{2}')}</td><td>${t('0.07779')}</td><td>${t('0.45724')}</td><td>${t('\\left[1+\\kappa_{\\mathrm{PR}}\\left(1-\\sqrt{T_r}\\right)\\right]^2')}</td></tr></table></div></div>
    <p>where ${t('\\omega')} is the acentric factor and</p>
    ${eq('\\kappa_{\\mathrm{SRK}} = 0.480 + 1.574\\,\\omega - 0.176\\,\\omega^2, \\qquad \\kappa_{\\mathrm{PR}} = 0.37464 + 1.54226\\,\\omega - 0.26992\\,\\omega^2', '2')}
    <h2>Compressibility factor</h2>
    <p>With ${t('A = \\dfrac{aP}{(RT)^2}')} and ${t('B = \\dfrac{bP}{RT}')}, the factor ${t('Z = \\dfrac{PV}{RT}')} satisfies a cubic. Below ${t('T_c')} it can have three real roots: the smallest is liquid-like, the largest vapour-like, the middle one is unstable.</p>
    ${eq('Z^3 + \\big[(\\epsilon+\\sigma)B - (1+B)\\big]Z^2 + \\big[\\epsilon\\sigma B^2 - (1+B)(\\epsilon+\\sigma)B + A\\big]Z - \\big[(1+B)\\epsilon\\sigma B^2 + AB\\big] = 0', '3')}
    <h2>Fugacity</h2>
    ${eq('\\ln\\varphi = Z - 1 - \\ln(Z-B) - \\frac{A}{B}\\,I', '4')}
    ${eq('I = \\frac{1}{\\sigma-\\epsilon}\\ln\\frac{Z+\\sigma B}{Z+\\epsilon B}, \\qquad f = \\varphi P', '5')}
    <p>For van der Waals (${t('\\sigma=\\epsilon=0')}) the integral reduces to ${t('I = B/Z')}. The stable phase is the root with the lowest ${t('\\ln\\varphi')}. The vapour pressure ${t('P^{\\mathrm{sat}}')} is found where ${t('\\ln\\varphi')} of the liquid root equals ${t('\\ln\\varphi')} of the vapour root.</p>
    <h2>Residual enthalpy and entropy</h2>
    ${eq('\\frac{H^R}{RT} = Z - 1 + \\left(\\frac{d\\ln\\alpha}{d\\ln T_r} - 1\\right)\\frac{A}{B}\\,I', '6')}
    ${eq('\\frac{S^R}{R} = \\ln(Z-B) + \\frac{d\\ln\\alpha}{d\\ln T_r}\\,\\frac{A}{B}\\,I', '7')}
    <p>Total properties add the ideal-gas contribution, with ${t('\\dfrac{C_p^{ig}}{R} = A_0 + B_0 T + C_0 T^2 + \\dfrac{D_0}{T^2}')}, measured from ${t('H = S = 0')} for the ideal gas at 298.15 K and 1 bar:</p>
    ${eq('H = H^{ig}(T) + H^R, \\qquad S = S^{ig}(T,P) + S^R, \\qquad \\Delta H^{\\mathrm{vap}} = H^V - H^L', '8')}
    <h2>Check of the numbers</h2>
    <ul><li>The JavaScript engine of this page agrees with the Python program in the repository on 2,760 test states, to about ${t('10^{-14}')}.</li>
    <li>The Python program agrees with CoolProp's independent PR and SRK implementations to about ${t('10^{-9}')}. The residual properties obey ${t('H^R = -RT^2\\left(\\dfrac{\\partial \\ln\\varphi}{\\partial T}\\right)_P')} and the Clapeyron equation.</li>
    <li>Gas data (${t('T_c')}, ${t('P_c')}, ${t('\\omega')}, ${t('C_p')}) come from Smith, Van Ness &amp; Abbott, Tables B.1 and C.1. Cubic equations are approximate: they are weakest for liquid volumes, near the critical point and for polar fluids such as water and ammonia.</li></ul>
    <p>Source code, tests and the project report are in the project repository (<code>python/</code>, <code>tests/</code>, <code>docs/report.pdf</code>).</p></div>`;
  }

  /* ---------- render ---------- */
  function renderStage() {
    const st = $('#stage');
    st.innerHTML = S.tab === 'diagram' ? viewDiagram() : S.tab === 'properties' ? viewProperties() : S.tab === 'compare' ? viewCompare() : viewTheory();
    if (S.tab === 'diagram') drawPV();
    if (S.tab === 'compare') drawCompare();
    if (S.tab === 'properties') { updateDelta(); $('#go2').addEventListener('click', updateDelta); ['#T2', '#P2'].forEach(i => $(i).addEventListener('input', updateDelta)); }
    if (S.tab === 'compare') $$('#seg button').forEach(b => b.addEventListener('click', () => { S.metric = b.dataset.metric; renderStage(); }));
    $('.readout').hidden = S.tab === 'theory';
  }

  let raf = 0;
  function update(opts) {
    opts = opts || {};
    cancelAnimationFrame(raf);
    raf = requestAnimationFrame(() => {
      try { compute(); $('#err').hidden = true; }
      catch (e) { const er = $('#err'); er.textContent = e.message; er.hidden = false; return; }
      if (!opts.keepInputs) {
        const ae = document.activeElement;
        if (ae !== $('#T')) $('#T').value = +S.T.toPrecision(7);
        if (ae !== $('#P')) $('#P').value = +S.P.toPrecision(7);
      }
      setSliders(); readout(); renderStage(); writeHash();
    });
  }

  /* ---------- wiring ---------- */
  function init() {
    readHash();
    const sel = $('#gas');
    sel.innerHTML = gasList.map(g => `<option value="${g.id}">${esc(g.name)} (${esc(g.formula)})</option>`).join('') + '<option value="custom">Custom gas…</option>';
    sel.value = S.gas;
    $('#custom').hidden = S.gas !== 'custom';
    $('#cTc').value = S.custom.Tc; $('#cPc').value = S.custom.Pc; $('#cW').value = S.custom.w; $('#cCp').value = S.custom.cp.join(' ');
    $$('#eosBtns button').forEach(b => b.classList.toggle('on', b.dataset.eos === S.eos));
    $$('#tabs button').forEach(b => b.classList.toggle('on', b.dataset.tab === S.tab));

    sel.addEventListener('change', () => {
      S.gas = sel.value; $('#custom').hidden = S.gas !== 'custom';
      const sp = species(); if (S.gas !== 'custom') { S.T = +(1.15 * sp.Tc).toPrecision(5); S.P = +sp.Pc.toPrecision(5); } S.T2 = null; update();
    });
    const readCustom = () => {
      const Tc = parseFloat($('#cTc').value), Pc = parseFloat($('#cPc').value), w = parseFloat($('#cW').value), cp = $('#cCp').value.trim().split(/\s+/).map(Number);
      if (Tc > 0 && Pc > 0 && isFinite(w) && cp.length === 4 && cp.every(isFinite)) { S.custom = { Tc, Pc, w, cp }; update(); }
    };
    ['#cTc', '#cPc', '#cW', '#cCp'].forEach(i => $(i).addEventListener('input', readCustom));
    $$('#eosBtns button').forEach(b => b.addEventListener('click', () => { S.eos = b.dataset.eos; $$('#eosBtns button').forEach(x => x.classList.toggle('on', x === b)); update(); }));
    $$('#tabs button').forEach(b => b.addEventListener('click', () => { S.tab = b.dataset.tab; $$('#tabs button').forEach(x => x.classList.toggle('on', x === b)); update(); }));
    $('#T').addEventListener('input', e => { const v = parseT(e.target.value); if (v > 0) { S.T = v; update({ keepInputs: true }); } });
    $('#P').addEventListener('input', e => { const v = parseP(e.target.value); if (v > 0) { S.P = v; update({ keepInputs: true }); } });
    $('#T').addEventListener('change', () => update());
    $('#P').addEventListener('change', () => update());
    $('#Tr').addEventListener('input', e => { const sp = species(), [a, b] = Trange(sp); S.T = +(a + (b - a) * e.target.value / 1000).toPrecision(5); update(); });
    $('#Pr').addEventListener('input', e => { const sp = species(), [a, b] = Prange(sp); S.P = +(a * Math.pow(b / a, e.target.value / 1000)).toPrecision(4); update(); });
    $('#share').addEventListener('click', () => {
      const b = $('#share'); const done = () => { b.textContent = 'Copied'; setTimeout(() => (b.textContent = 'Copy link'), 1400); };
      (navigator.clipboard ? navigator.clipboard.writeText(location.href).then(done, done) : done());
    });
    if (window.ResizeObserver) new ResizeObserver(() => { if (S.tab === 'diagram') drawPV(); if (S.tab === 'compare') drawCompare(); }).observe($('#stage'));
    window.addEventListener('hashchange', () => { readHash(); update(); });
    update();
  }
  init();
})();
