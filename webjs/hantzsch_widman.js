"use strict";
// P-22.2.2 Hantzsch-Widman 杂单环命名：3-10 元环的词干、位次与名称组装。
// 未命中保留模板的杂单环由本模块生成词干（饱和/mancude 两形态），
// 位次按 P-22.2.2.1.3 元素优先序排列，省略判定见 P-22.2.2.1.7。

const { C, N, MULT_EN, MULT_ZH, P145_SENIOR, zhNumeral } = require("./constants");
const { kekulized, sssrRings } = require("./layer1_rings");

const HW_ID = "hw_mono";  // 生成式 HW 杂单环骨架 id：不登记进模板表，靠生成器产出词干
const HW_CLASS = "heterocycle";  // 其 naming_class（与 monohetero 分列）
const HW_COMPONENT_PREFIX = "hw:";  // 稠合组分编码前缀：'hw:' + 环序元素符号串

const HW_PREFIX_EN = { 9: "fluora", 17: "chlora", 35: "broma", 53: "ioda",  // Table 2.4 'a' 前缀（按优先性递减）
  8: "oxa", 16: "thia", 34: "selena", 52: "tellura",
  7: "aza", 15: "phospha", 33: "arsa", 51: "stiba", 83: "bisma",
  14: "sila", 32: "germa", 50: "stanna", 82: "plumba",
  5: "bora", 13: "aluma", 31: "galla", 49: "indiga", 81: "thalla" };
const HW_PREFIX_ZH = { 9: "氟杂", 17: "氯杂", 35: "溴杂", 53: "碘杂",  // 表 3-3 中文前缀（长式）
  8: "氧杂", 16: "硫杂", 34: "硒杂", 52: "碲杂",
  7: "氮杂", 15: "磷杂", 33: "砷杂", 51: "锑杂", 83: "铋杂",
  14: "硅杂", 32: "锗杂", 50: "锡杂", 82: "铅杂",
  5: "硼杂", 13: "铝杂", 31: "镓杂", 49: "铟杂", 81: "铊杂" };
const HW_ZH_SHORT = new Map([[508, "噁"], [516, "噻"], [616, "噻"]]);  // 含氮环短式，键为 环大小*100+原子序数
const HW_MAX_VALENCE = { 1: 1, 5: 3, 6: 4, 7: 3, 8: 2, 9: 1, 13: 3, 14: 4, 15: 3, 16: 2, 17: 1,  // Table 2.4 键数，供 mancude 参照
  31: 3, 32: 4, 33: 3, 34: 2, 35: 1, 49: 3, 50: 4, 51: 3, 52: 2, 53: 1, 81: 3, 82: 4, 83: 3 };
const HW_VOWELS = new Set(["a", "e", "i", "o", "u"]);  // 首尾元音省略判据（P-22.2.2.1.1）
const HW_SIX_A = new Set([8, 16, 34, 52, 83]);  // Table 2.5 六元环 A 组
const HW_SIX_B = new Set([7, 14, 32, 50, 82]);  // Table 2.5 六元环 B 组（余者归 C 组）
const HW_UNSAT_SIX = { A: "ine", B: "ine", C: "inine" };  // 六元不饱和词干（P-22.2.2.1.6）
const HW_SAT_SIX = { A: "ane", B: "inane", C: "inane" };  // 六元饱和词干
const HW_UNSAT_TAIL = { 7: "epine", 8: "ocine", 9: "onine", 10: "ecine" };  // Table 2.5 七至十元不饱和
const HW_SAT_TAIL = { 7: "epane", 8: "ocane", 9: "onane", 10: "ecane" };
const HW_ZH_RING = { 3: "环丙", 4: "环丁", 5: "环戊", 6: "环己",  // 中文基干环前缀（表 3-4）
  7: "环庚", 8: "环辛", 9: "环壬", 10: "环癸" };

// 元素符号 ↔ 原子序数（'hw:' 组分键的编解码；词表外元素返回 null）
const Z_BY_SYMBOL = { H: 1, B: 5, C: 6, N: 7, O: 8, F: 9, Al: 13, Si: 14, P: 15, S: 16, Cl: 17,
  Ga: 31, Ge: 32, As: 33, Se: 34, Br: 35, In: 49, Sn: 50, Sb: 51, Te: 52, I: 53, Tl: 81, Pb: 82, Bi: 83 };
const SYMBOL_BY_Z = {};
for (const [sym, z] of Object.entries(Z_BY_SYMBOL)) SYMBOL_BY_Z[z] = sym;

let _scaffold = null;
let _numbering = null;

function scaffoldDeps() {
  // 延迟取 ring_scaffold：本模块被 ring_scaffold 反向引用，直接 require 会成环。
  if (_scaffold === null) {
    try { _scaffold = require("./ring_scaffold"); } catch (e) { _scaffold = false; }
  }
  return _scaffold || null;
}

function numberingDeps() {
  // 延迟取 L4 编号引擎：保证名中位次与 L4 编号同源。
  if (_numbering === null) {
    try { _numbering = require("./numbering_engine"); } catch (e) { _numbering = false; }
  }
  return _numbering || null;
}

// ── 词干与名称组装 ──────────────

function ringHeteros(zs) {
  return zs.filter((z) => z !== C);
}

function sixGroup(zs) {
  // 六元环词干分组：由优先性最低的杂原子（名称紧邻词干者）所属组决定（P-22.2.2.1.6）
  let least = null;
  for (const z of ringHeteros(zs)) {
    if (least === null || P145_SENIOR.indexOf(z) > P145_SENIOR.indexOf(least)) least = z;
  }
  if (HW_SIX_A.has(least)) return "A";
  return HW_SIX_B.has(least) ? "B" : "C";
}

function unsatStem(n, zs) {
  // 不饱和（mancude）词干（Table 2.5、P-22.2.2.1.5.1）
  if (n === 3) return new Set(ringHeteros(zs)).size === 1 && zs.indexOf(N) >= 0 ? "irine" : "irene";
  if (n === 4) return "ete";
  if (n === 5) return "ole";
  if (n === 6) return HW_UNSAT_SIX[sixGroup(zs)];
  return HW_UNSAT_TAIL[n] || null;
}

function satStem(n, zs) {
  // 饱和词干（Table 2.5、P-22.2.2.1.5.2）
  const hasN = zs.indexOf(N) >= 0;
  if (n === 3) return hasN ? "iridine" : "irane";
  if (n === 4) return hasN ? "etidine" : "etane";
  if (n === 5) return hasN ? "olidine" : "olane";
  if (n === 6) return HW_SAT_SIX[sixGroup(zs)];
  return HW_SAT_TAIL[n] || null;
}

function zhStem(n, zs, nDb) {
  // 中文基干：含氮五/六元环取唑/嗪系，其余按环大小与实际环内双键数
  const hasN = zs.indexOf(N) >= 0;
  if (n === 5 && hasN) return nDb ? "唑" : "唑烷";
  if (n === 6 && hasN) return nDb ? "嗪" : "嗪烷";
  const ring = HW_ZH_RING[n];
  if (ring === undefined) return null;
  return ring + zhUnsatTail(nDb);
}

function zhUnsatTail(nDb) {
  // 中文双键词尾：0 烷、1 烯、n 二烯…
  if (nDb <= 0) return "烷";
  if (nDb === 1) return "烯";
  return `${zhNumeral(nDb)}烯`;
}

function elideA(terms) {
  // 按 P-22.2.2.1.1 拼接：前项末尾 'a' 在后项起首元音前省略
  let out = terms[0];
  for (const term of terms.slice(1)) {
    if (out.endsWith("a") && HW_VOWELS.has(term[0])) out = out.slice(0, -1);
    out += term;
  }
  return out;
}

function locantMap(zs) {
  // 按 P-22.2.2.1.3 引用顺序分组的 {元素: 位次升序}（zs 为已编号的环序）
  const seen = new Map();
  zs.forEach((z, i) => {
    if (z === C) return;
    if (!seen.has(z)) seen.set(z, []);
    seen.get(z).push(i + 1);
  });
  const out = new Map();
  for (const z of P145_SENIOR) if (seen.has(z)) out.set(z, seen.get(z));
  return out;
}

function locantString(byZ) {
  // 位次串按引用顺序排列（非数值升序），如 1,6,2-dioxazepane 的 '1,6,2'
  const locs = [];
  for (const v of byZ.values()) locs.push(...v);
  return locs.join(",");
}

function heteroTerms(n, zs, byZ, en) {
  // 元素前缀项（含倍数词）；中文在含氮五/六元环里把 N 折入唑/嗪基干，故中文跳过 N
  const hasN = zs.indexOf(N) >= 0;
  const folded = !en && hasN && (n === 5 || n === 6);
  const short = hasN && (n === 5 || n === 6);  // 噁/噻短式只在唑/嗪系（含氮环）启用
  const terms = [];
  for (const [z, locs] of byZ) {
    if (folded && z === N) continue;
    const mult = (en ? MULT_EN : MULT_ZH)[locs.length];
    let word;
    if (en) word = HW_PREFIX_EN[z];
    else word = (short ? HW_ZH_SHORT.get(n * 100 + z) : undefined) || HW_PREFIX_ZH[z];
    terms.push(mult ? elideA([mult, word]) : word);  // tetra+aza → tetraza
  }
  return terms;
}

function zhFoldedN(n, zs, byZ, zhBase) {
  // 折入唑/嗪的氮数量词：三唑、二嗪烷（仅五/六元含氮环的唑嗪系）
  if (!(zs.indexOf(N) >= 0 && (n === 5 || n === 6))) return zhBase;
  const count = (byZ.get(N) || []).length;
  return count ? MULT_ZH[count] + zhBase : zhBase;
}

function hwNameFromCycle(zs, nDb) {
  // 已编号环序元素表 + 环内双键数 → [英文名, 中文名, 英文位次前缀]；第 i 个原子位次为 i+1
  if (!(zs.length >= 3 && zs.length <= 10)) return null;
  const byZ = locantMap(zs);
  if (!byZ.size) return null;
  const stem = nDb ? unsatStem(zs.length, zs) : satStem(zs.length, zs);
  let zhBase = zhStem(zs.length, zs, nDb);
  if (stem === null || zhBase === null) return null;
  zhBase = zhFoldedN(zs.length, zs, byZ, zhBase);
  const locs = omitLocants(zs.length, zs) ? "" : `${locantString(byZ)}-`;
  const en = `${locs}${elideA(heteroTerms(zs.length, zs, byZ, true).concat([stem]))}`;
  const zh = `${locs}${heteroTerms(zs.length, zs, byZ, false).join("")}${zhBase}`;
  return [en, zh, locs];
}

// ── 位次省略（P-22.2.2.1.7）─────────────

const _arrCache = new Map();

function omitLocants(n, zs) {
  // 单杂原子、或元素组成在环上排布唯一时省略全部位次
  const heteros = ringHeteros(zs);
  if (!heteros.length) return false;
  if (heteros.length === 1) return true;
  return arrangementClasses(n, heteros.slice().sort((a, b) => a - b)) === 1;
}

function combos(n, k) {
  // 从 1..n-1 中取 k 个位置的全部升序组合
  const out = [];
  const walk = (start, cur) => {
    if (cur.length === k) { out.push(cur.slice()); return; }
    for (let v = start; v <= n - (k - cur.length); v++) { cur.push(v); walk(v + 1, cur); cur.pop(); }
  };
  walk(1, []);
  return out;
}

function permUnique(arr) {
  // 多重集全排列（去重）
  const sorted = arr.slice().sort((a, b) => a - b);
  const out = [], used = new Array(sorted.length).fill(false);
  const walk = (cur) => {
    if (cur.length === sorted.length) { out.push(cur.slice()); return; }
    for (let i = 0; i < sorted.length; i++) {
      if (used[i]) continue;
      if (i > 0 && sorted[i] === sorted[i - 1] && !used[i - 1]) continue;
      used[i] = true; cur.push(sorted[i]); walk(cur); cur.pop(); used[i] = false;
    }
  };
  walk([]);
  return out;
}

function rotationKey(seq) {
  // 环上序列的旋转最小表示（元素名定宽，避免 '10' < '2' 的字典序陷阱）
  const pad = seq.map((z) => String(z).padStart(2, "0"));
  let best = null;
  for (let i = 0; i < pad.length; i++) {
    const key = pad.slice(i).concat(pad.slice(0, i)).join("");
    if (best === null || key < best) best = key;
  }
  return best;
}

function arrangementClasses(n, heteros) {
  // 杂原子多重集在 n 元环上的循环排布数（模旋转与翻转）
  const key = `${n}|${heteros.join(",")}`;
  const hit = _arrCache.get(key);
  if (hit !== undefined) return hit;
  const k = heteros.length;
  const classes = new Set();
  for (const spots of combos(n, k - 1)) {
    const pos = [0].concat(spots);
    for (const perm of permUnique(heteros)) {
      const seq = new Array(n).fill(C);
      pos.forEach((p, j) => { seq[p] = perm[j]; });
      const rev = seq.slice().reverse();
      classes.add(rotationKey(seq) < rotationKey(rev) ? rotationKey(seq) : rotationKey(rev));
    }
  }
  const out = classes.size;
  _arrCache.set(key, out);
  return out;
}

// ── 环分析 ──────────────

function isHwScaffold(sid) {
  // sid 是否为生成式 HW 杂单环骨架
  return sid === HW_ID;
}

function inScope(zs) {
  // 环内元素是否都在 HW 词表内（词表外元素不得臆造前缀）
  return zs.every((z) => z === C || HW_PREFIX_EN[z] !== undefined);
}

function isolatedRing(mol, atoms) {
  // 原子集是否为不与他环稠合的单环（稠合/桥环另走 P-25/P-23）
  const ring = new Set(atoms), rings = sssrRings(mol);
  if (!rings.some((r) => r.length === ring.size && r.every((a) => ring.has(a)))) return false;
  for (const a of ring) if (rings.filter((r) => r.indexOf(a) >= 0).length !== 1) return false;
  return true;
}

function isDoubleBond(b) {
  // kekulé 视图下的双键（芳香标记不参与，只看实际键级）
  return b.boDouble === 2 || b.type === "DOUBLE";
}

function ringDoubleBonds(mol, chain) {
  // 环内 Kekulé 双键数（芳香环按 Kekulé 视图计）
  const kek = kekulized(mol) || mol, atoms = new Set(chain);
  let n = 0;
  for (const b of kek.g.bonds) if (isDoubleBond(b) && atoms.has(b.a1) && atoms.has(b.a2)) n += 1;
  return n;
}

function ringNumbering(mol, chain) {
  // P-22.2.2.1.3 环编号：复用 L4 窄化引擎，保证名中位次与 L4 编号同源
  const ne = numberingDeps();
  if (ne === null) return null;
  const cands = ne.narrowHeteroRing(ne.ringCands(chain.slice()), mol, chain.slice(), false);
  return cands.length ? ne.toChain(cands[0]) : null;
}

function ringOrder(mol, atoms) {
  // 按环连接走出原子集的环序；原子集不是单环或输入顺序无关
  const rset = new Set(atoms), g = mol.g, adj = new Map();
  for (const a of rset) {
    const nb = [];
    for (const bi of g.atoms[a].bonds) {
      const b = g.bonds[bi], o = b.a1 === a ? b.a2 : b.a1;
      if (rset.has(o)) nb.push(o);
    }
    if (nb.length !== 2) return null;
    adj.set(a, nb);
  }
  const start = Math.min(...rset);
  const out = [start, adj.get(start)[0]];
  while (out.length < rset.size) {
    const nxt = adj.get(out[out.length - 1]).filter((x) => x !== out[out.length - 2]);
    if (!nxt.length) return null;
    out.push(nxt[0]);
  }
  return out;
}

// ── mancude 参照与加氢位（P-31.2）─────────────

function maxMatchingMask(n, capable) {
  // 环上最大匹配（边 i 连接 i 与 i+1）；n<=10 直接枚举边子集
  let best = -1, bestMask = 0;
  for (let mask = 0; mask < (1 << n); mask++) {
    if (mask & (mask >> 1)) continue;  // 相邻边不可同选
    if ((mask & 1) && ((mask >> (n - 1)) & 1)) continue;  // 首尾边不可同选
    let ok = true;
    for (let i = 0; i < n; i++) if ((mask >> i) & 1) { if (!(capable[i] && capable[(i + 1) % n])) { ok = false; break; } }
    if (!ok) continue;
    let cnt = 0;
    for (let m = mask; m; m >>= 1) cnt += m & 1;
    if (cnt > best) { best = cnt; bestMask = mask; }
  }
  return bestMask;
}

function popcount(mask) {
  let n = 0;
  for (let m = mask; m; m >>= 1) n += m & 1;
  return n;
}

function effectiveValence(atom) {
  // 有效价：Table 2.4 键数 + 形式电荷（N+ 视作 4 价、O- 视作 1 价）
  const base = HW_MAX_VALENCE[atom.z];
  return base === undefined ? null : base + atom.chg;
}

function outerOrderSum(kek, a, rset) {
  // 环外键级和（不含氢）
  let s = 0;
  for (const bi of kek.g.atoms[a].bonds) {
    const b = kek.g.bonds[bi], o = b.a1 === a ? b.a2 : b.a1;
    if (rset.has(o) || kek.g.atoms[o].z === 1) continue;
    s += b.boDouble;
  }
  return s;
}

function mancudeHydrogens(mol, chain) {
  // mancude 参照：[{环原子: 参照氢数}, 参照环内双键数]（P-22.2.2.1.1/P-31.2）
  const order = ringOrder(mol, chain);
  if (order === null) return null;
  const kek = kekulized(mol) || mol;
  const rset = new Set(order), n = order.length;
  const outer = new Map(), caps = new Map(), capable = [];
  for (const a of order) {
    const atom = kek.g.atoms[a], cap = effectiveValence(atom);
    if (cap === null) return null;
    caps.set(a, cap);
    const o = outerOrderSum(kek, a, rset);
    outer.set(a, o);
    capable.push(o + 3 <= cap);  // 环内两键 + 一个 π 键不超价
  }
  const mask = maxMatchingMask(n, capable), ref = new Map();
  order.forEach((a, i) => {
    const matched = ((mask >> i) & 1) === 1 || ((mask >> ((i - 1 + n) % n)) & 1) === 1;
    ref.set(a, caps.get(a) - (matched ? 3 : 2) - outer.get(a));
  });
  return [ref, popcount(mask)];
}

function hydroAtoms(mol, chain) {
  // 加氢位：实际氢数多于 mancude 参照的环原子（P-31.2.2 / P-54.4.1）
  // 环内无重键时取饱和词干，加氢前缀不再使用（否则饱和环会被误标 dihydro）
  if (ringDoubleBonds(mol, chain) === 0) return new Set();
  const ref = mancudeHydrogens(mol, chain);
  if (ref === null) return new Set();
  const out = new Set();
  for (const [a, h] of ref[0]) if (mol.g.atoms[a].nHs > h) out.add(a);
  return out;
}

// ── 稠合组分（P-25.2.2.1.1）─────────────

function parseKey(sid) {
  // 'hw:SiCCO' → 环序元素表（按编号顺序）
  const body = sid.slice(HW_COMPONENT_PREFIX.length), zs = [];
  for (const m of body.matchAll(/[A-Z][a-z]?/g)) {
    const z = Z_BY_SYMBOL[m[0]];
    if (z === undefined || (z !== C && HW_PREFIX_EN[z] === undefined)) return null;
    zs.push(z);
  }
  return zs.length ? zs : null;
}

function mancudeDb(zs) {
  // mancude 环内双键数（纯拓扑最大匹配：价 >= 3 者才可承担 π 键）
  const capable = zs.map((z) => HW_MAX_VALENCE[z] >= 3);
  return popcount(maxMatchingMask(zs.length, capable));
}

function componentKey(mol, atoms) {
  // 稠合组分键 'hw:OCOCC'：仅 3-10 元、含杂、mancude 的孤立环（P-25.2.2.1.1）
  const order = ringOrder(mol, atoms);
  if (order === null || !(order.length >= 3 && order.length <= 10)) return null;
  const zs = order.map((a) => mol.g.atoms[a].z);
  if (!ringHeteros(zs).length || !inScope(zs)) return null;
  const ref = mancudeHydrogens(mol, order);
  if (ref === null || ringDoubleBonds(mol, order) !== ref[1]) return null;  // 非 mancude 环不作组分
  const ordered = ringNumbering(mol, order);
  if (ordered === null) return null;
  return HW_COMPONENT_PREFIX + ordered.map((a) => SYMBOL_BY_Z[mol.g.atoms[a].z]).join("");
}

function componentNames(sid) {
  // 'hw:…' → 稠合母体组分词干 [en, zh]（P-25.3.2.1.2）；附加组分前缀由 fused_namer 通用式给出
  const zs = parseKey(sid);
  if (zs === null || !(zs.length >= 3 && zs.length <= 10)) return null;
  const names = hwNameFromCycle(zs, mancudeDb(zs));
  if (names === null) return null;
  const en = names[0], zh = names[1], locs = names[2];
  if (!locs) return [en, zh];
  const bare = locs.slice(0, -1);  // 稠合组分位次须加方括号（P-25.3.2.1.2）
  return [`[${bare}]${en.slice(locs.length)}`, `[${bare}]${zh.slice(locs.length)}`];
}

// ── 接线入口 ──────────────

function parentNames(sid, mol, chain) {
  // 生成式 HW 母体双语名（含位次前缀）；非 HW 骨架或词表外元素返回 null
  if (!isHwScaffold(sid) || mol == null || !chain || !chain.length) return null;
  const ordered = ringNumbering(mol, chain);
  if (ordered === null) return null;
  const zs = ordered.map((a) => mol.g.atoms[a].z);
  if (!inScope(zs)) return null;
  const names = hwNameFromCycle(zs, ringDoubleBonds(mol, ordered));
  return names ? [names[0], names[1]] : null;
}

function identity(info, skeleton) {
  // 未命中保留模板的 3-10 元孤立含杂单环 → 生成式 HW 骨架身份（P-22.2.2）
  const mol = info.mol, atoms = skeleton.atom_ids.slice();
  if (mol == null || !(atoms.length >= 3 && atoms.length <= 10)) return null;
  const zs = atoms.map((a) => mol.g.atoms[a].z);
  if (!ringHeteros(zs).length || !inScope(zs)) return null;  // 全碳环走通用 carbocycle
  if (!isolatedRing(mol, atoms)) return null;
  return scaffoldDeps().scaffoldIdentity(HW_ID, HW_CLASS, 1, "hetero");
}

module.exports = {
  HW_ID, HW_CLASS, HW_COMPONENT_PREFIX,
  ringHeteros, sixGroup, unsatStem, satStem, zhStem, elideA, locantMap, locantString,
  hwNameFromCycle, omitLocants, isHwScaffold, inScope, isolatedRing, ringDoubleBonds,
  ringNumbering, ringOrder, mancudeHydrogens, hydroAtoms, componentKey, componentNames,
  parentNames, identity, _parse_key: parseKey, _mancude_db: mancudeDb,
};
