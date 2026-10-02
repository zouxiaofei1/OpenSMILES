// 把 webjs/*.js 重新打包回单文件 index.min.html（页面内联载荷 = base85(brotli(JSON))）。
// 用法: node webjs/pack.mjs [--template index.min.html] [--mods webjs] [--out tmp/index.min.html]
import fs from "node:fs";
import path from "node:path";
import zlib from "node:zlib";
import { execFileSync } from "node:child_process";
import { fileURLToPath } from "node:url";

const HERE = path.dirname(fileURLToPath(import.meta.url));
const ROOT = path.resolve(HERE, "..");

function opt(name, def) {
  const i = process.argv.indexOf("--" + name);
  return i >= 0 ? process.argv[i + 1] : def;
}
const TEMPLATE = opt("template", path.join(ROOT, "index.min.html"));
const MODS_DIR = opt("mods", HERE);
const OUT = opt("out", path.join(ROOT, "tmp", "index.min.html"));

// 页面用的自定义 base85（字母表不含引号/反斜杠/尖括号/斜杠，可直接放进 script 文本）
const X85A = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!#$%&()*+,-.:;=?@^_[]{}";
const X85R = new Int16Array(128);
for (let i = 0; i < 85; i++) X85R[X85A.charCodeAt(i)] = i;

// 字节 -> base85 文本：4 字节成 5 字符；末组按页面解码器口径补 'u'(84) 位
function x85Encode(buf) {
  const n = buf.length, full = Math.floor(n / 4), rem = n - full * 4;
  let out = "";
  for (let g = 0; g < full; g++) {
    const v = ((buf[g * 4] << 24) | (buf[g * 4 + 1] << 16) | (buf[g * 4 + 2] << 8) | buf[g * 4 + 3]) >>> 0;
    let s = "", x = v;
    for (let k = 0; k < 5; k++) { s = X85A[x % 85] + s; x = Math.floor(x / 85); }
    out += s;
  }
  if (rem === 0) return out;
  // 末组：k 个实际字节 -> k+1 个字符（解码器按 rem 读字符、余位补 'u'，输出 rem-1 字节）
  const k = rem;
  const remChars = k + 1;
  let top = 0;
  for (let m = 0; m < k; m++) top = (top + buf[full * 4 + m] * 2 ** (24 - 8 * m)) >>> 0;
  const padDigits = 5 - remChars;
  let U = 0;
  for (let j = 0; j < padDigits; j++) U = U * 85 + 84;
  const P = 85 ** padDigits;
  const w = top + (((U - (top % P)) % P) + P) % P;
  let v = Math.floor(w / P), s = "";
  for (let j = 0; j < remChars; j++) { s = X85A[v % 85] + s; v = Math.floor(v / 85); }
  return out + s;
}

// 页面同款解码器（仅用于打包后自检）
function x85Decode(s) {
  s = s.replace(/\s+/g, "");
  const n = s.length, full = Math.floor(n / 5), rem = n - full * 5;
  const out = new Uint8Array(full * 4 + (rem > 0 ? rem - 1 : 0));
  let o = 0, i = 0;
  for (let g = 0; g < full; g++) {
    let v = 0;
    for (let k = 0; k < 5; k++) v = (v * 85 + X85R[s.charCodeAt(i++)]) >>> 0;
    out[o++] = v >>> 24; out[o++] = (v >>> 16) & 255; out[o++] = (v >>> 8) & 255; out[o++] = v & 255;
  }
  if (rem > 0) {
    let w = 0;
    for (let j = 0; j < 5; j++) {
      const d = j < rem ? X85R[s.charCodeAt(i++)] : 84;
      w = (w * 85 + d) >>> 0;
    }
    for (let m = 0; m < rem - 1; m++) out[o++] = (w >>> (24 - 8 * m)) & 255;
  }
  return Buffer.from(out);
}

// 自检：任意长度字节经 encode->decode 必须原样还原
function selfTest() {
  for (let n = 0; n < 40; n++) {
    const b = Buffer.alloc(n);
    for (let i = 0; i < n; i++) b[i] = (i * 37 + n * 11) & 255;
    const rt = x85Decode(x85Encode(b));
    if (!rt.equals(b)) throw new Error(`x85 自检失败 n=${n}`);
  }
}

// esbuild 压缩单文件源码；不可用时保留原源码（载荷随后仍会 brotli）
function minify(src) {
  try {
    return execFileSync("npx", ["--no-install", "esbuild", "--minify", "--loader=js"],
      { input: src, encoding: "utf8", shell: true, maxBuffer: 1 << 26 }).trim();
  } catch (e) {
    console.warn("  minify 失败，保留原源码:", String(e.message).slice(0, 120));
    return src;
  }
}

function loadModules(dir, doMinify) {
  const mods = {};
  for (const f of fs.readdirSync(dir).filter((f) => f.endsWith(".js")).sort()) {
    const src = fs.readFileSync(path.join(dir, f), "utf8");
    mods[f] = doMinify ? minify(src) : src;
  }
  if (!Object.keys(mods).length) throw new Error(`模块目录为空: ${dir}`);
  return mods;
}

selfTest();
const mods = loadModules(MODS_DIR, !process.argv.includes("--no-minify"));
const raw = Buffer.from(JSON.stringify(mods), "utf8");
const packed = zlib.brotliCompressSync(raw, {
  params: { [zlib.constants.BROTLI_PARAM_QUALITY]: 11, [zlib.constants.BROTLI_PARAM_SIZE_HINT]: raw.length },
});
const payload = x85Encode(packed);
// 回读校验：解码 -> brotli -> JSON 必须与写入的模块表一致
const back = JSON.parse(zlib.brotliDecompressSync(x85Decode(payload)).toString("utf8"));
if (JSON.stringify(back) !== JSON.stringify(mods)) throw new Error("回读校验失败");

const html = fs.readFileSync(TEMPLATE, "utf8");
const re = /MODULES_PAYLOAD\s*=\s*"([^"]*)"/;
if (!re.test(html)) throw new Error(`模板缺少 MODULES_PAYLOAD: ${TEMPLATE}`);
// 必须用函数式替换：载荷含 $ & 字符，字符串替换里的 $& 等会被当成替换模式展开
const outHtml = html.replace(re, () => `MODULES_PAYLOAD = "${payload}"`);
fs.writeFileSync(OUT, outHtml, "utf8");
console.log(`modules=${Object.keys(mods).length} raw=${raw.length}B brotli=${packed.length}B payload=${payload.length}B`);
console.log(`template=${path.relative(ROOT, TEMPLATE)} -> out=${path.relative(ROOT, OUT)} (${outHtml.length}B)`);
