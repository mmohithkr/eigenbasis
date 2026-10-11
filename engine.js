/* Eigenbasis numerics. Plain functions on typed arrays, no DOM; index.html draws the results.
   Every basis here is orthonormal in L²(dx) on its domain, so cₙ = ∫ f φₙ dx, Sₙ = Σ_{n<N} cₙ φₙ and
   Parseval reads Σ cₙ² = ‖f‖². Integrals use composite Gauss–Legendre on panels that (a) split f where it
   is not smooth and (b) are fine enough to resolve the highest basis function. */
(function (root) {
'use strict';
const PI = Math.PI;

/* ---------- Gauss–Legendre rules on [-1, 1] */
const GLC = new Map();
function legP(m, z) { // returns [P_m(z), P_{m-1}(z)]
  let p1 = 1, p2 = 0;
  for (let j = 1; j <= m; j++) { const p3 = p2; p2 = p1; p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j; }
  return [p1, p2];
}
function gaussLegendre(m) {
  let r = GLC.get(m); if (r) return r;
  const x = new Float64Array(m), w = new Float64Array(m);
  for (let i = 0; i < (m + 1) >> 1; i++) {
    let z = Math.cos(PI * (i + 0.75) / (m + 0.5));
    for (let it = 0; it < 60; it++) {
      const [p1, p2] = legP(m, z), pp = m * (z * p1 - p2) / (z * z - 1), dz = p1 / pp;
      z -= dz; if (Math.abs(dz) < 4e-16) break;
    }
    const [p1, p2] = legP(m, z), pp = m * (z * p1 - p2) / (z * z - 1);
    x[i] = -z; x[m - 1 - i] = z; w[i] = w[m - 1 - i] = 2 / ((1 - z * z) * pp * pp);
  }
  r = { x, w }; GLC.set(m, r); return r;
}

/* ---------- breakpoint helpers */
function uniform(a, b, n) { const o = []; for (let i = 0; i <= n; i++) o.push(a + (b - a) * i / n); o[n] = b; return o; }
function uniq(arr, L) {
  const s = arr.filter(Number.isFinite).sort((p, q) => p - q), o = [];
  const tol = 1e-13 * L;
  for (const v of s) if (!o.length || v - o[o.length - 1] > tol) o.push(v);
  return o;
}

/* ---------- adaptive splitting: cut [a, b] until f looks like a degree-12 polynomial on every panel */
const CM = 16, CT = new Float64Array(CM), CCOS = new Float64Array(CM * CM);
for (let j = 0; j < CM; j++) CT[j] = Math.cos(PI * (j + 0.5) / CM);
for (let k = 0; k < CM; k++) for (let j = 0; j < CM; j++) CCOS[k * CM + j] = Math.cos(PI * k * (j + 0.5) / CM);

class NotFinite extends Error { constructor(x) { super('not finite'); this.x = x; } }

function sampleScale(f, a, b) {
  let G = 0;
  for (let i = 0; i <= 4000; i++) {
    const x = a + (b - a) * (i + 0.5) / 4001, v = f(x);
    if (!Number.isFinite(v)) throw new NotFinite(x);
    if (Math.abs(v) > G) G = Math.abs(v);
  }
  return G;
}

function adaptiveBreaks(f, a, b, seeds, G) {
  const L = b - a, minW = 1e-14 * L, tol = 1e-13 * Math.max(G, 1e-300);
  const start = uniq(uniform(a, b, 32).concat(seeds.filter(s => s > a && s < b)), L);
  const out = [a], sing = [], v = new Float64Array(CM);
  let budget = 9000;
  const rec = (l, r) => {
    const mid = 0.5 * (l + r), half = 0.5 * (r - l);
    for (let j = 0; j < CM; j++) { const x = mid + half * CT[j], y = f(x); if (!Number.isFinite(y)) throw new NotFinite(x); v[j] = y; }
    let tail = 0;
    for (let k = CM - 3; k < CM; k++) { let s = 0; for (let j = 0; j < CM; j++) s += v[j] * CCOS[k * CM + j]; if (Math.abs(s) > tail) tail = Math.abs(s); }
    tail *= 2 / CM;
    if (tail <= tol || budget <= 0) { out.push(r); return; }
    if (r - l <= minW) { out.push(r); sing.push(mid); return; }
    budget--;
    rec(l, mid); rec(mid, r);
  };
  for (let i = 0; i + 1 < start.length; i++) rec(start[i], start[i + 1]);
  // clusters of unresolved panels mark a jump or a kink
  const pts = [];
  for (const s of sing) { const last = pts[pts.length - 1]; if (last && s - last.x1 < 1e-10 * L) { last.x1 = s; last.sum += s; last.n++; } else pts.push({ x1: s, sum: s, n: 1 }); }
  return { breaks: out, sing: pts.map(p => p.sum / p.n) };
}

function quadFromBreaks(br, m, L) {
  const g = gaussLegendre(m), gs = gaussLegendre(6);
  let count = 0;
  for (let p = 0; p + 1 < br.length; p++) count += (br[p + 1] - br[p] < L / 2000) ? 6 : m;
  const X = new Float64Array(count), W = new Float64Array(count);
  let k = 0;
  for (let p = 0; p + 1 < br.length; p++) {
    const l = br[p], r = br[p + 1], mid = 0.5 * (l + r), half = 0.5 * (r - l), q = (r - l < L / 2000) ? gs : g;
    for (let j = 0; j < q.x.length; j++) { X[k] = mid + half * q.x[j]; W[k] = half * q.w[j]; k++; }
  }
  return { X, W };
}

/* ---------- subscripts for labels */
const SUB = '₀₁₂₃₄₅₆₇₈₉';
const sub = n => String(n).split('').map(d => SUB[+d] || d).join('');

/* ---------- the bases */
const ISQ2PI = 1 / Math.sqrt(2 * PI), ISQPI = 1 / Math.sqrt(PI), SQ2 = Math.SQRT2;

const fourier = {
  id: 'fourier', name: 'Fourier', dom: [-PI, PI], quad: [-PI, PI], nMax: 201, m: 20, first: 0, envW: 4,
  basisBreaks: () => uniform(-PI, PI, 64),
  evalAll(x, out, n) {
    if (n > 0) out[0] = ISQ2PI;
    const c1 = Math.cos(x), s1 = Math.sin(x);
    let c = c1, s = s1;
    for (let i = 1, k = 1; i < n; k++) {
      out[i++] = c * ISQPI; if (i < n) out[i++] = s * ISQPI;
      if ((k & 15) === 15) { c = Math.cos((k + 1) * x); s = Math.sin((k + 1) * x); }
      else { const cn = c * c1 - s * s1; s = s * c1 + c * s1; c = cn; }
    }
  },
  label(i) {
    if (i === 0) return '1/√(2π)';
    const k = (i + 1) >> 1, kx = (k === 1 ? '' : k) + 'x';
    return (i & 1 ? 'cos ' : 'sin ') + kx + ' / √π';
  },
  short(i) { if (i === 0) return 'constant'; const k = (i + 1) >> 1; return (i & 1 ? 'cos ' : 'sin ') + (k === 1 ? '' : k) + 'x'; },
  peak(N) { const K = (N - 1) >> 1; return K >= 2 ? PI / (K + 0.5) : 0; },
  hints: [[1, 'a jump in f, or between its two ends'], [1.5, 'a square-root cusp'], [2, 'a kink: f′ jumps, or the end slopes differ'], [3, 'a jump in f″']],
  // a free particle on a ring of circumference 2π: H = −½ d²/dx², E = k²/2 for cos kx and sin kx
  En: i => { const k = (i + 1) >> 1; return 0.5 * k * k; }, period: 4 * PI, periodTxt: '4π', system: 'free particle on a ring',
};

const sine = {
  id: 'sine', name: 'Particle in a box', dom: [0, 1], quad: [0, 1], nMax: 200, m: 20, first: 1, envW: 3, qm: true,
  basisBreaks: () => uniform(0, 1, 64),
  evalAll(x, out, n) {
    const t = PI * x, c1 = Math.cos(t), s1 = Math.sin(t);
    let c = c1, s = s1;
    for (let i = 0; i < n; i++) {
      out[i] = SQ2 * s;
      if ((i & 15) === 15) { c = Math.cos((i + 2) * t); s = Math.sin((i + 2) * t); }
      else { const cn = c * c1 - s * s1; s = s * c1 + c * s1; c = cn; }
    }
  },
  label(i) { const n = i + 1; return '√2 sin ' + (n === 1 ? '' : n) + 'πx'; },
  short(i) { return 'φ' + sub(i + 1); },
  energy(i) { const n = i + 1; return 'E' + sub(n) + ' = ' + n * n + ' E' + sub(1); },
  peak(N) { return N >= 3 ? 1 / N : 0; },
  hints: [[1, 'a jump, or f ≠ 0 at a wall'], [1.5, 'a square-root cusp'], [2, 'a kink, where f′ jumps'], [3, 'a jump in f″']],
  En: i => { const n = i + 1; return 0.5 * PI * PI * n * n; }, E1: 0.5 * PI * PI, period: 4 / PI, periodTxt: '4/π', system: 'particle in a box',
};

const LNRM = new Float64Array(260); for (let n = 0; n < 260; n++) LNRM[n] = Math.sqrt((2 * n + 1) / 2);
const legendre = {
  id: 'legendre', name: 'Legendre', dom: [-1, 1], quad: [-1, 1], nMax: 200, m: 24, first: 0, envW: 3,
  basisBreaks: () => { const o = []; for (let j = 0; j <= 64; j++) o.push(-Math.cos(PI * j / 64)); o[0] = -1; o[64] = 1; return o; },
  evalAll(x, out, n) {
    let p0 = 1, p1 = x;
    if (n > 0) out[0] = LNRM[0];
    if (n > 1) out[1] = x * LNRM[1];
    for (let k = 1; k + 1 < n; k++) { const p2 = ((2 * k + 1) * x * p1 - k * p0) / (k + 1); p0 = p1; p1 = p2; out[k + 1] = p2 * LNRM[k + 1]; }
  },
  label(i) { return '√(' + (2 * i + 1) + '/2) P' + sub(i) + '(x)'; },
  short(i) { return 'P' + sub(i); },
  peak(N, x0) { return N >= 3 ? PI * Math.sqrt(Math.max(1 - x0 * x0, 0)) / N + 2 / (N * N) : 0; },
  hints: [[1, 'a jump inside the interval'], [1.5, 'a square-root cusp'], [2, 'a kink, where f′ jumps'], [3, 'a jump in f″']],
  // a rigid rotor in its m = 0 states, x = cos θ: E_l = l(l + 1)/2, all whole numbers
  En: i => 0.5 * i * (i + 1), period: 2 * PI, periodTxt: '2π', system: 'rigid rotor (m = 0)',
};

const HS1 = new Float64Array(260), HS2 = new Float64Array(260);
for (let n = 0; n < 260; n++) { HS1[n] = Math.sqrt(2 / (n + 1)); HS2[n] = Math.sqrt(n / (n + 1)); }
const PIQ = Math.pow(PI, -0.25);
const hermite = {
  id: 'hermite', name: 'Harmonic oscillator', dom: [-6, 6], quad: [-25, 25], nMax: 200, m: 18, first: 0, envW: 4, envProp: 0.15, qm: true,
  basisBreaks: () => uniform(-25, 25, 250),
  evalAll(x, out, n) {
    let h0 = PIQ * Math.exp(-0.5 * x * x), h1 = SQ2 * x * h0;
    if (n > 0) out[0] = h0;
    if (n > 1) out[1] = h1;
    for (let k = 1; k + 1 < n; k++) { const h2 = HS1[k] * x * h1 - HS2[k] * h0; h0 = h1; h1 = h2; out[k + 1] = h2; }
  },
  label(i) { return 'ψ' + sub(i) + '(x)'; },
  short(i) { return 'ψ' + sub(i); },
  energy(i) { return 'E' + sub(i) + ' = ' + i + '½'; },
  peak(N, x0) { const q = 2 * N + 1 - x0 * x0; return N >= 4 && q > 1 ? PI / Math.sqrt(q) : 0; },
  hints: [[0.75, 'a jump'], [1.25, 'a kink, where f′ jumps']],
  En: i => i + 0.5, period: 2 * PI, periodTxt: '2π', system: 'harmonic oscillator',
};

const haar = {
  id: 'haar', name: 'Haar', dom: [0, 1], quad: [0, 1], nMax: 256, m: 12, first: 0, envW: 0,
  basisBreaks: () => uniform(0, 1, 512),
  evalAll(x, out, n) {
    if (n > 0) out[0] = (x >= 0 && x <= 1) ? 1 : 0;
    for (let i = 1; i < n; i++) out[i] = 0;
    if (!(x >= 0 && x <= 1)) return;
    for (let j = 0, p = 1; p < n; j++, p *= 2) {
      const k = Math.min(Math.floor(x * p), p - 1), i = p + k;
      if (i < n) { const t = x * p - k; out[i] = (t < 0.5 ? 1 : -1) * Math.sqrt(p); }
    }
  },
  label(i) { if (i === 0) return '1'; const j = Math.floor(Math.log2(i)), k = i - (1 << j); return 'h' + sub(j) + ',' + sub(k) + '(x)'; },
  short(i) { if (i === 0) return '1'; const j = Math.floor(Math.log2(i)), k = i - (1 << j); return 'h' + sub(j) + ',' + sub(k); },
  peak(N) { return N >= 3 ? 1 / N : 0; },
  hints: [[0.5, 'a jump'], [1.5, 'a continuous f (Haar only sees slopes)']],
};

const BASES = { fourier, sine, legendre, hermite, haar };

/* ---------- preset functions, written on u ∈ [−1, 1] and stretched onto each basis's window */
const PRESETS = [
  { id: 'square', name: 'Square wave', f: u => u > 0 ? 1 : u < 0 ? -1 : 0, breaks: [0] },
  { id: 'box', name: 'Box', f: u => { const a = Math.abs(u); return a < 0.5 ? 1 : a > 0.5 ? 0 : 0.5; }, breaks: [-0.5, 0.5] },
  { id: 'tent', name: 'Tent', f: u => Math.max(0, 1 - 2 * Math.abs(u)), breaks: [-0.5, 0, 0.5] },
  { id: 'ramp', name: 'Ramp', f: u => u, breaks: [] },
  { id: 'gaussian', name: 'Gaussian', f: u => Math.exp(-18 * u * u), breaks: [], decays: true },
  { id: 'shifted', name: 'Shifted Gaussian', short: 'Shifted', f: u => Math.exp(-18 * (u - 1 / 3) * (u - 1 / 3)), breaks: [], decays: true },
  { id: 'runge', name: 'Runge bump', f: u => 1 / (1 + 25 * u * u), breaks: [] },
  { id: 'exp', name: 'Exponential', f: u => Math.exp(u), breaks: [] },
  { id: 'cusp', name: 'Cusp', f: u => Math.sqrt(Math.abs(u)), breaks: [0] },
];

/* ---------- formula parser: x, numbers, + − × ÷ ^, |…|, comparisons (give 1 or 0), implicit products, the usual functions */
const FN1 = {
  sin: Math.sin, cos: Math.cos, tan: Math.tan, sec: x => 1 / Math.cos(x), csc: x => 1 / Math.sin(x), cot: x => 1 / Math.tan(x),
  asin: Math.asin, acos: Math.acos, atan: Math.atan, arcsin: Math.asin, arccos: Math.acos, arctan: Math.atan,
  sinh: Math.sinh, cosh: Math.cosh, tanh: Math.tanh, sech: x => 1 / Math.cosh(x), csch: x => 1 / Math.sinh(x), coth: x => 1 / Math.tanh(x), exp: Math.exp, log: Math.log, ln: Math.log, log10: Math.log10, log2: Math.log2,
  sqrt: Math.sqrt, cbrt: Math.cbrt, abs: Math.abs, sign: Math.sign, sgn: Math.sign, floor: Math.floor, ceil: Math.ceil, round: Math.round,
  step: x => x > 0 ? 1 : x < 0 ? 0 : 0.5, heaviside: x => x > 0 ? 1 : x < 0 ? 0 : 0.5,
  sinc: x => x === 0 ? 1 : Math.sin(x) / x,
};
const FN2 = { min: Math.min, max: Math.max, atan2: Math.atan2, pow: Math.pow, mod: (a, b) => a - b * Math.floor(a / b) };
const CONST = { pi: PI, e: Math.E, tau: 2 * PI };
const NAMES = ['x', ...Object.keys(FN1), ...Object.keys(FN2), ...Object.keys(CONST)].sort((p, q) => q.length - p.length);

class ParseError extends Error {}
function tokenize(src) {
  const s = src.replace(/[−–—]/g, '-').replace(/[×·⋅∙]/g, '*').replace(/÷/g, '/').replace(/π/g, ' pi ').replace(/τ/g, ' tau ')
    .replace(/√/g, ' sqrt ').replace(/\*\*/g, '^').replace(/²/g, '^2').replace(/³/g, '^3').replace(/≤/g, '<=').replace(/≥/g, '>=');
  const out = [];
  let i = 0;
  while (i < s.length) {
    const ch = s[i];
    if (/\s/.test(ch)) { i++; continue; }
    if (/[0-9.]/.test(ch)) {
      const m = /^(\d+\.?\d*|\.\d+)(e[+-]?\d+)?/i.exec(s.slice(i));
      if (!m) throw new ParseError('I can’t read the number at “' + s.slice(i, i + 6) + '”.');
      out.push({ t: 'num', v: parseFloat(m[0]) }); i += m[0].length; continue;
    }
    if (/[a-zA-Z_]/.test(ch)) {
      const m = /^[a-zA-Z_][a-zA-Z0-9_]*/.exec(s.slice(i))[0]; i += m.length;
      let w = m.toLowerCase();
      if (NAMES.includes(w)) { out.push({ t: 'id', v: w }); continue; }
      // split run-together names such as "xsin" or "pix"
      const parts = [];
      while (w.length) {
        const nm = NAMES.find(n => w.startsWith(n));
        if (!nm) {
          const d = /^\d+/.exec(w);
          if (d && parts.length) { parts.push({ t: 'num', v: +d[0] }); w = w.slice(d[0].length); continue; }
          throw new ParseError('Unknown name “' + m + '”. Use x as the variable.');
        }
        parts.push({ t: 'id', v: nm }); w = w.slice(nm.length);
      }
      out.push(...parts); continue;
    }
    const two = s.slice(i, i + 2);
    if (two === '<=' || two === '>=') { out.push({ t: 'op', v: two }); i += 2; continue; }
    if ('+-*/^(),|<>%[]{}'.includes(ch)) { out.push({ t: 'op', v: ch === '[' || ch === '{' ? '(' : ch === ']' || ch === '}' ? ')' : ch }); i++; continue; }
    throw new ParseError('Unexpected “' + ch + '”.');
  }
  return out;
}

function parseFormula(src) {
  if (!src || !src.trim()) throw new ParseError('Type a formula in x, for example abs(x) or x^2.');
  const T = tokenize(src);
  let p = 0, absDepth = 0;
  const peek = () => T[p], isOp = v => T[p] && T[p].t === 'op' && T[p].v === v;
  const expect = v => { if (!isOp(v)) throw new ParseError(v === ')' ? 'A bracket is not closed.' : v === '|' ? 'An absolute-value bar is not closed.' : 'Expected “' + v + '”.'); p++; };
  const startsPrimary = () => { const t = peek(); if (!t) return false; if (t.t !== 'op') return true; return t.v === '(' || (t.v === '|' && absDepth === 0); };

  function cmp() {
    let a = add();
    while (isOp('<') || isOp('>') || isOp('<=') || isOp('>=')) {
      const op = T[p++].v, b = add(), A = a;
      a = op === '<' ? x => +(A(x) < b(x)) : op === '>' ? x => +(A(x) > b(x)) : op === '<=' ? x => +(A(x) <= b(x)) : x => +(A(x) >= b(x));
    }
    return a;
  }
  function add() {
    let a = mul();
    while (isOp('+') || isOp('-')) { const op = T[p++].v, b = mul(), A = a; a = op === '+' ? x => A(x) + b(x) : x => A(x) - b(x); }
    return a;
  }
  function mul() {
    let a = unary();
    for (;;) {
      if (isOp('*') || isOp('/') || isOp('%')) {
        const op = T[p++].v, b = unary(), A = a;
        a = op === '*' ? x => A(x) * b(x) : op === '/' ? x => A(x) / b(x) : x => { const d = b(x); return A(x) - d * Math.floor(A(x) / d); };
      } else if (startsPrimary()) { const b = unary(), A = a; a = x => A(x) * b(x); }
      else return a;
    }
  }
  function unary() {
    if (isOp('-')) { p++; const a = unary(); return x => -a(x); }
    if (isOp('+')) { p++; return unary(); }
    return pow();
  }
  function pow() {
    const a = primary();
    if (isOp('^')) { p++; const b = unary(); return x => Math.pow(a(x), b(x)); }
    return a;
  }
  function primary() {
    const t = peek();
    if (!t) throw new ParseError('The formula ends too early.');
    if (t.t === 'num') { p++; const v = t.v; return () => v; }
    if (t.t === 'op' && t.v === '(') { p++; const a = cmp(); expect(')'); return a; }
    if (t.t === 'op' && t.v === '|') { p++; absDepth++; const a = cmp(); absDepth--; expect('|'); return x => Math.abs(a(x)); }
    if (t.t === 'id') {
      p++;
      if (t.v === 'x') return x => x;
      if (t.v in CONST) { const v = CONST[t.v]; return () => v; }
      if (t.v in FN2) {
        const fn = FN2[t.v]; expect('(');
        const a = cmp(); if (!isOp(',')) throw new ParseError(t.v + ' needs two arguments, like ' + t.v + '(x, 0).');
        p++; const b = cmp(); expect(')'); return x => fn(a(x), b(x));
      }
      if (t.v in FN1) {
        const fn = FN1[t.v];
        if (isOp('(')) { p++; const a = cmp(); expect(')'); if (isOp('^')) { p++; const e = unary(); return x => Math.pow(fn(a(x)), e(x)); } return x => fn(a(x)); }
        const a = pow(); return x => fn(a(x)); // sin x, sqrt x
      }
    }
    if (t.t === 'op' && t.v === ')') throw new ParseError('There is a “)” without a matching “(”.');
    throw new ParseError('Unexpected “' + t.v + '”.');
  }
  let f;
  try { f = cmp(); }
  catch (e) {
    let depth = 0; for (const t of T) if (t.t === 'op') depth += t.v === '(' ? 1 : t.v === ')' ? -1 : 0;
    if (e instanceof ParseError && depth > 0 && p >= T.length) throw new ParseError('A bracket is not closed.');
    throw e;
  }
  if (p < T.length) { const t = T[p]; throw new ParseError(t.v === ')' ? 'There is a “)” without a matching “(”.' : 'Unexpected “' + t.v + '”.'); }
  return f;
}

/* ---------- a function to expand, placed on a basis */
function makeSource(basis, spec) {
  const [a, b] = basis.dom, toU = x => 2 * (x - a) / (b - a) - 1, toX = u => a + (u + 1) * (b - a) / 2;
  const isH = basis.id === 'hermite';
  let f, breaks = [], windowed = false;
  if (spec.kind === 'preset') {
    const P = PRESETS.find(q => q.id === spec.id) || PRESETS[0];
    breaks = P.breaks.map(toX);
    if (isH && !P.decays) { windowed = true; f = x => (x >= a && x <= b) ? P.f(toU(x)) : 0; }
    else f = x => P.f(toU(x));
  } else if (spec.kind === 'formula') {
    const g = spec.fn || parseFormula(spec.text);
    if (isH) {
      // keep the whole formula if it has died away by |x| = 20; otherwise cut it off outside the plotted window
      let big = 0, tailMax = 0, ok = true;
      for (let i = 0; i <= 2000; i++) {
        const x = -25 + 50 * i / 2000, v = g(x);
        if (!Number.isFinite(v)) { if (Math.abs(x) > 6) { ok = false; break; } else continue; }
        big = Math.max(big, Math.abs(v)); if (Math.abs(x) >= 20) tailMax = Math.max(tailMax, Math.abs(v));
      }
      if (ok && big > 0 && tailMax <= 1e-14 * big) f = g;
      else { windowed = true; f = x => (x >= a && x <= b) ? g(x) : 0; }
    } else f = g;
  } else if (spec.kind === 'draw') {
    const v = spec.vals, K = v.length, h = (b - a) / (K - 1);
    for (let k = 1; k < K - 1; k++) breaks.push(a + k * h);
    f = x => {
      if (x < a || x > b) return 0;
      let t = (x - a) / h, k = Math.floor(t); if (k >= K - 1) k = K - 2; if (k < 0) k = 0;
      t -= k; return v[k] * (1 - t) + v[k + 1] * t;
    };
    if (isH) windowed = true;
  }
  if (windowed) breaks.push(a, b);
  // a starting momentum: the real or imaginary part of f(x)·e^{ik₀x}
  if (spec.k0 && spec.part) { const g0 = f, k0 = spec.k0; f = spec.part === 'sin' ? x => g0(x) * Math.sin(k0 * x) : x => g0(x) * Math.cos(k0 * x); }
  return { f, breaks, windowed, spec };
}

/* ---------- jumps, for the Gibbs readout */
function detectJumps(basis, src, sing, brk, G) {
  const [a, b] = basis.dom, L = b - a, eps = 1e-9 * L, thr = 1e-6 * Math.max(G, 1e-300), f = src.f, out = [];
  const near = (x, y) => Math.abs(x - y) < 1e-12 * L;
  const cands = uniq(src.breaks.concat(sing, brk), L);
  for (const c0 of cands) {
    const edge = near(c0, a) || near(c0, b);
    if (edge && !src.windowed) continue;
    if (!edge && !(c0 > a && c0 < b)) continue;
    if (out.some(j => Math.abs(j.x0 - c0) < 4 * eps)) continue; // a neighbour of a jump already found
    const vL = f(c0 - eps), vR = f(c0 + eps);
    if (!(Math.abs(vR - vL) > thr)) continue;
    // pin the jump down to rounding level by bisection
    let lo = c0 - eps, hi = c0 + eps;
    if (!edge) for (let it = 0; it < 80 && hi - lo > 4e-16 * Math.max(1, Math.abs(c0)); it++) {
      const mid = 0.5 * (lo + hi), vm = f(mid);
      if (Math.abs(vm - vL) <= Math.abs(vm - vR)) lo = mid; else hi = mid;
    }
    const x0 = edge ? c0 : 0.5 * (lo + hi);
    if (out.some(j => Math.abs(j.x0 - x0) < 4 * eps)) continue;
    const fL = f(x0 - eps), fR = f(x0 + eps);
    out.push({ x0, fL, fR, J: fR - fL, sides: [{ x0, dir: -1, high: fL > fR }, { x0, dir: 1, high: fR > fL }] });
  }
  if (basis.id === 'fourier') {
    const fL = f(b - eps), fR = f(a + eps);
    if (Math.abs(fR - fL) > thr) out.push({ x0: b, seam: true, fL, fR, J: fR - fL, sides: [{ x0: b, dir: -1, high: fL > fR }, { x0: a, dir: 1, high: fR > fL }] });
  }
  if (basis.id === 'sine') {
    const f0 = f(a + eps), f1 = f(b - eps);
    if (Math.abs(f0) > thr) out.push({ x0: a, wall: true, fL: -f0, fR: f0, J: 2 * f0, sides: [{ x0: a, dir: 1, high: f0 > 0 }] });
    if (Math.abs(f1) > thr) out.push({ x0: b, wall: true, fL: f1, fR: -f1, J: -2 * f1, sides: [{ x0: b, dir: -1, high: f1 > 0 }] });
  }
  out.sort((p, q) => Math.abs(q.J) - Math.abs(p.J));
  return out;
}

/* ---------- how fast do the coefficients fall? */
function linfit(X, Y) {
  const n = X.length; let sx = 0, sy = 0, sxx = 0, sxy = 0;
  for (let i = 0; i < n; i++) { sx += X[i]; sy += Y[i]; sxx += X[i] * X[i]; sxy += X[i] * Y[i]; }
  const d = n * sxx - sx * sx; if (Math.abs(d) < 1e-300) return null;
  const s = (n * sxy - sx * sy) / d, c = (sy - s * sx) / n;
  let rss = 0; for (let i = 0; i < n; i++) { const r = Y[i] - (c + s * X[i]); rss += r * r; }
  return { s, c, rss };
}
function decayFit(basis, c) {
  const n = c.length, A = Array.from(c, Math.abs);
  const amax = Math.max(...A);
  if (!(amax > 0)) return { kind: 'zero' };
  const floor = 1e-12 * amax;
  let last = -1, sig = 0;
  for (let i = 0; i < n; i++) if (A[i] > floor) { last = i; sig++; }
  if (sig <= 8 && last < n - 1) return { kind: 'exact', terms: last + 1 };
  let X = [], Y = [];
  if (basis.id === 'haar') {
    for (let j = 2, p = 4; p < n; j++, p *= 2) {
      let m = 0; for (let i = p; i < Math.min(2 * p, n); i++) m = Math.max(m, A[i]);
      if (m > 10 * floor) { X.push(Math.log(1.5 * p)); Y.push(Math.log(m)); }
    }
    const fp = X.length >= 3 ? linfit(X, Y) : null;
    return fp ? { kind: 'power', p: -fp.s, c0: fp.c, i0: 4, i1: n - 1 } : { kind: 'unclear' };
  }
  // upper envelope: parity zeros and slow beats would otherwise wreck a straight-line fit
  const env = new Float64Array(n), Wof = i => basis.envW + Math.round((basis.envProp || 0) * i);
  for (let i = 0; i < n; i++) { const W = Wof(i); let m = 0; for (let j = Math.max(0, i - W); j <= Math.min(n - 1, i + W); j++) m = Math.max(m, A[j]); env[i] = m; }
  // where the coefficients have already fallen below 10⁻⁷ of the largest one for good
  let tiny = n; for (let i = n - 1; i >= 0; i--) { if (env[i] > 1e-7 * amax) break; tiny = i; }
  if (last < 0.85 * n) {
    const i0 = Math.max(1, Math.floor(0.3 * last));
    for (let i = i0; i <= last; i++) if (env[i] > 10 * floor) { X.push(i); Y.push(Math.log(env[i])); }
    const fg = X.length >= 4 ? linfit(X, Y) : null;
    return fg && fg.s < 0 ? { kind: 'geom', rate: -fg.s, c0: fg.c, i0: X[0], i1: X[X.length - 1], floorAt: last + 1 } : { kind: 'unclear' };
  }
  const i0 = Math.floor(n / 4), i1 = n - 1 - Wof(n - 1), Xl = [];
  for (let i = i0; i <= i1; i++) if (env[i] > 10 * floor) { X.push(i); Xl.push(Math.log(i + 1)); Y.push(Math.log(env[i])); }
  if (X.length < 6) return { kind: 'unclear' };
  const fp = linfit(Xl, Y), fg = linfit(X, Y);
  const tail = tiny < i0 ? { tinyAt: tiny } : {};
  const span = { i0: X[0], i1: X[X.length - 1] };
  // call it geometric only when a straight line in n is clearly better AND it actually falls a fair way over the window;
  // a nearly flat tail is a slow power law (drawn functions, jumps), not a sign of analyticity
  const geomWins = fg && fg.s < 0 && fp && fg.rss < 0.5 * fp.rss && -fg.s * (span.i1 - span.i0) > 3;
  if (fp && !geomWins) return Object.assign({ kind: 'power', p: -fp.s, c0: fp.c }, span, tail);
  if (fg && fg.s < 0) return Object.assign({ kind: 'geom', rate: -fg.s, c0: fg.c }, span, tail);
  return { kind: 'unclear' };
}
function decayHint(basis, fit) {
  if (fit.kind !== 'power' || !basis.hints) return '';
  for (const [p, t] of basis.hints) if (Math.abs(fit.p - p) < 0.15) return t;
  return '';
}

/* ---------- the whole expansion */
const GRID = 1601;
function expand(basisId, spec) {
  const basis = BASES[basisId] || fourier, src = makeSource(basis, spec), f = src.f;
  const [qa, qb] = basis.quad, QL = qb - qa, nMax = basis.nMax;
  const G = sampleScale(f, qa, qb);
  const ad = adaptiveBreaks(f, qa, qb, src.breaks, G);
  const br = uniq(ad.breaks.concat(basis.basisBreaks()), QL);
  const { X, W } = quadFromBreaks(br, basis.m, QL), M = X.length;
  const F = new Float64Array(M); for (let i = 0; i < M; i++) F[i] = f(X[i]);
  const buf = new Float64Array(nMax), c = new Float64Array(nMax);
  let norm2 = 0;
  for (let i = 0; i < M; i++) {
    basis.evalAll(X[i], buf, nMax);
    const wf = W[i] * F[i]; norm2 += wf * F[i];
    for (let n = 0; n < nMax; n++) c[n] += wf * buf[n];
  }
  // ‖f − S_N‖² by quadrature for every N (no cancellation, so it stays accurate down to rounding level)
  const E = new Float64Array(nMax + 1);
  for (let i = 0; i < M; i++) {
    basis.evalAll(X[i], buf, nMax);
    const w = W[i], fi = F[i]; let s = 0;
    E[0] += w * fi * fi;
    for (let n = 0; n < nMax; n++) { s += c[n] * buf[n]; const d = fi - s; E[n + 1] += w * d * d; }
  }
  const cum = new Float64Array(nMax + 1);
  for (let n = 0; n < nMax; n++) cum[n + 1] = cum[n] + c[n] * c[n];

  // values on the plotting grid: f, every partial sum, every basis function
  const [a, b] = basis.dom, xs = new Float64Array(GRID), fs = new Float64Array(GRID);
  const P = new Float32Array((nMax + 1) * GRID), Phi = new Float32Array(nMax * GRID), sup = new Float64Array(nMax + 1);
  let fmax = 0;
  for (let g = 0; g < GRID; g++) {
    const x = a + (b - a) * g / (GRID - 1); xs[g] = x;
    let fv = f(x); if (!Number.isFinite(fv)) fv = NaN;
    const ok = fv === fv;
    fs[g] = fv; if (ok) fmax = Math.max(fmax, Math.abs(fv));
    basis.evalAll(x, buf, nMax);
    let s = 0; if (ok) sup[0] = Math.max(sup[0], Math.abs(fv));
    for (let n = 0; n < nMax; n++) {
      Phi[n * GRID + g] = buf[n]; s += c[n] * buf[n]; P[(n + 1) * GRID + g] = s;
      if (ok) { const e = Math.abs(fv - s); if (e > sup[n + 1]) sup[n + 1] = e; }
    }
  }
  const jumps = detectJumps(basis, src, ad.sing, ad.breaks, G);
  const fit = decayFit(basis, c);
  return { basis, src, c, norm2, E, cum, xs, fs, P, Phi, sup, fmax, jumps, fit, hint: decayHint(basis, fit), nQuad: M, nPanels: br.length - 1 };
}

/* S_N at any point */
function partialAt(R, x, N, buf) {
  buf = buf || new Float64Array(R.basis.nMax);
  R.basis.evalAll(x, buf, N);
  let s = 0; for (let n = 0; n < N; n++) s += R.c[n] * buf[n];
  return s;
}

/* Gibbs overshoot next to the biggest jump, as a fraction of the jump */
function gibbs(R, N) {
  const { basis, jumps, src } = R;
  if (!jumps.length || N < 3) return null;
  const J0 = jumps[0], d = basis.peak(N, J0.x0);
  if (!d) return null;
  const [a, b] = basis.dom, buf = new Float64Array(basis.nMax), K = 240;
  const others = jumps.slice(1).map(j => j.x0);
  let best = -Infinity, bx = 0, by = 0;
  for (const s of J0.sides) {
    let len = 3 * d;
    // the ringing has to fit before the next jump, or the overshoot is not this jump's alone
    for (const o of others) { const dist = (o - s.x0) * s.dir; if (dist > 0 && dist < len) return null; }
    const outward = (s.dir > 0 && s.x0 >= b - 1e-12) || (s.dir < 0 && s.x0 <= a + 1e-12);
    if (!(src.windowed && outward)) len = Math.min(len, s.dir > 0 ? b - s.x0 : s.x0 - a);
    if (!(len > 0)) continue;
    const sg = s.high ? 1 : -1, ev = x => (partialAt(R, x, N, buf) - src.f(x)) * sg;
    let lb = -Infinity, lx = 0;
    for (let k = 1; k <= K; k++) { const x = s.x0 + s.dir * len * k / K, e = ev(x); if (e > lb) { lb = e; lx = x; } }
    const h = len / K; // refine around the best sample
    for (let k = -20; k <= 20; k++) { const x = lx + h * k / 20; if ((x - s.x0) * s.dir <= 0) continue; const e = ev(x); if (e > lb) { lb = e; lx = x; } }
    if (lb > best) { best = lb; bx = lx; by = partialAt(R, lx, N, buf); }
  }
  if (best === -Infinity) return null;
  return { frac: Math.max(0, best) / Math.abs(J0.J), x: bx, y: by, jump: J0 };
}

/* ---------- completeness: the kernel K_N(x, x′) = Σ_{n<N} φₙ(x) φₙ(x′) */
// one row, K_N(x, x₀) on the plotting grid; returns K_N(x₀, x₀)
function kernelRow(R, x0, N, out) {
  const B = R.basis, G = R.xs.length, w = new Float64Array(B.nMax);
  B.evalAll(x0, w, N);
  out.fill(0);
  let d = 0;
  for (let n = 0; n < N; n++) {
    const a = w[n]; d += a * a; if (!a) continue;
    const base = n * G; for (let g = 0; g < G; g++) out[g] += a * R.Phi[base + g];
  }
  return d;
}
// the whole kernel on an M × M grid over the window
function kernelGrid(R, N, M) {
  const G = R.xs.length, idx = new Int32Array(M), K = new Float64Array(M * M), col = new Float64Array(M);
  for (let j = 0; j < M; j++) idx[j] = Math.round(j * (G - 1) / (M - 1));
  for (let n = 0; n < N; n++) {
    const base = n * G; for (let j = 0; j < M; j++) col[j] = R.Phi[base + idx[j]];
    for (let i = 0; i < M; i++) { const a = col[i]; if (!a) continue; const r = i * M; for (let j = i; j < M; j++) K[r + j] += a * col[j]; }
  }
  for (let i = 0; i < M; i++) for (let j = 0; j < i; j++) K[i * M + j] = K[j * M + i];
  return K;
}

/* ---------- time evolution: ψ(x, t) = Σ cₙ e^{−iEₙt} φₙ(x), with complex cₙ = (re) + i (im), normalised over n < N */
function amplitudes(RE, RI, N) {
  const a = new Float64Array(N), b = new Float64Array(N);
  let nn = 0;
  for (let n = 0; n < N; n++) { a[n] = RE.c[n]; b[n] = RI ? RI.c[n] : 0; nn += a[n] * a[n] + b[n] * b[n]; }
  const s = nn > 0 ? 1 / Math.sqrt(nn) : 0;
  for (let n = 0; n < N; n++) { a[n] *= s; b[n] *= s; }
  return { a, b, norm2: nn };
}
function evolveAt(RE, amp, t, outRe, outIm, idx) {
  // idx: optional list of grid indices to evaluate (for the carpet); otherwise the whole plotting grid
  const B = RE.basis, G = RE.xs.length, N = amp.a.length, M = idx ? idx.length : G;
  outRe.fill(0); outIm.fill(0);
  for (let n = 0; n < N; n++) {
    const a = amp.a[n], b = amp.b[n]; if (!a && !b) continue;
    const th = B.En(n) * t, cs = Math.cos(th), sn = Math.sin(th);
    const ar = a * cs + b * sn, ai = b * cs - a * sn, base = n * G;
    if (idx) for (let j = 0; j < M; j++) { const p = RE.Phi[base + idx[j]]; outRe[j] += ar * p; outIm[j] += ai * p; }
    else for (let g = 0; g < G; g++) { const p = RE.Phi[base + g]; outRe[g] += ar * p; outIm[g] += ai * p; }
  }
}
// |ψ(x, t)|² over one period on a cols × rows grid (t down the rows)
function carpet(RE, amp, cols, rows) {
  const B = RE.basis, G = RE.xs.length, idx = new Int32Array(cols), out = new Float32Array(cols * rows);
  for (let j = 0; j < cols; j++) idx[j] = Math.round(j * (G - 1) / (cols - 1));
  const re = new Float64Array(cols), im = new Float64Array(cols);
  let max = 0;
  for (let r = 0; r < rows; r++) {
    evolveAt(RE, amp, B.period * r / (rows - 1), re, im, idx);
    for (let j = 0; j < cols; j++) { const v = re[j] * re[j] + im[j] * im[j]; out[r * cols + j] = v; if (v > max) max = v; }
  }
  return { data: out, cols, rows, max };
}
// |⟨ψ(0)|ψ(t)⟩|² = |Σ pₙ e^{−iEₙt}|², the chance of finding the starting state again
function returnProb(basis, amp, t) {
  let sr = 0, si = 0;
  for (let n = 0; n < amp.a.length; n++) { const p = amp.a[n] * amp.a[n] + amp.b[n] * amp.b[n]; if (!p) continue; const th = basis.En(n) * t; sr += p * Math.cos(th); si -= p * Math.sin(th); }
  return sr * sr + si * si;
}

const API = { BASES, PRESETS, expand, partialAt, gibbs, parseFormula, ParseError, NotFinite, gaussLegendre, quadFromBreaks, uniform, sub, kernelRow, kernelGrid, amplitudes, evolveAt, carpet, returnProb };
if (typeof module !== 'undefined' && module.exports) module.exports = API;
else root.EB = API;
})(typeof self !== 'undefined' ? self : this);
