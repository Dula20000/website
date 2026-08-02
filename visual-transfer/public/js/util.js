/** Small shared helpers. No dependencies. */

const CRC_TABLE = (() => {
  const t = new Uint32Array(256);
  for (let n = 0; n < 256; n++) {
    let c = n;
    for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1;
    t[n] = c >>> 0;
  }
  return t;
})();

/** CRC-32 (IEEE) over bytes — used to prove a reassembled file is byte-exact. */
export function crc32(bytes) {
  return crc32Final(crc32Update(crc32Init(), bytes));
}

/* Streaming form, so a multi-gigabyte file can be checksummed chunk by chunk
   without ever being held in memory whole. */
export const crc32Init = () => 0xffffffff;

export function crc32Update(state, bytes) {
  let c = state;
  for (let i = 0; i < bytes.length; i++) c = CRC_TABLE[(c ^ bytes[i]) & 0xff] ^ (c >>> 8);
  return c >>> 0;
}

export const crc32Final = (state) => (state ^ 0xffffffff) >>> 0;

export function xorInto(target, source) {
  const n = Math.min(target.length, source.length);
  let i = 0;
  // 4 bytes at a time where alignment allows; XOR is the decoder's hot loop.
  if (n >= 16 && target.byteOffset % 4 === 0 && source.byteOffset % 4 === 0) {
    const words = (n / 4) | 0;
    const a = new Uint32Array(target.buffer, target.byteOffset, words);
    const b = new Uint32Array(source.buffer, source.byteOffset, words);
    for (let w = 0; w < words; w++) a[w] ^= b[w];
    i = words * 4;
  }
  for (; i < n; i++) target[i] ^= source[i];
  return target;
}

export function formatBytes(n) {
  if (!Number.isFinite(n)) return '—';
  if (n < 1024) return `${n} B`;
  const units = ['KB', 'MB', 'GB', 'TB'];
  let v = n / 1024;
  let u = 0;
  while (v >= 1024 && u < units.length - 1) { v /= 1024; u++; }
  return `${v < 10 ? v.toFixed(1) : Math.round(v)} ${units[u]}`;
}

export function formatRate(bytesPerSec) {
  if (!Number.isFinite(bytesPerSec) || bytesPerSec <= 0) return '—';
  return `${formatBytes(bytesPerSec)}/s`;
}

export function formatDuration(seconds) {
  if (!Number.isFinite(seconds) || seconds < 0) return '—';
  if (seconds < 60) return `${Math.ceil(seconds)}s`;
  const m = Math.floor(seconds / 60);
  const s = Math.ceil(seconds % 60);
  return `${m}m ${String(s).padStart(2, '0')}s`;
}

/** Strip path segments and characters that break downloads on some platforms. */
export function safeFilename(name) {
  const base = String(name || 'file').split(/[/\\]/).pop() || 'file';
  const clean = base.replace(/[\x00-\x1f<>:"|?*]/g, '_').trim();
  return clean.slice(0, 180) || 'file';
}

/** Best-effort icon for a file, driven by MIME first then extension. */
export function fileGlyph(name, mime) {
  const m = (mime || '').toLowerCase();
  if (m.startsWith('image/')) return '🖼️';
  if (m.startsWith('video/')) return '🎬';
  if (m.startsWith('audio/')) return '🎵';
  if (m === 'application/pdf') return '📕';
  if (m.startsWith('text/') || m.includes('json') || m.includes('xml')) return '📄';
  const ext = (name || '').split('.').pop().toLowerCase();
  if (['zip', 'gz', 'tar', 'rar', '7z', 'bz2', 'xz'].includes(ext)) return '🗜️';
  if (['doc', 'docx', 'odt', 'rtf', 'pages'].includes(ext)) return '📝';
  if (['xls', 'xlsx', 'csv', 'numbers'].includes(ext)) return '📊';
  if (['ppt', 'pptx', 'key'].includes(ext)) return '📽️';
  if (['js', 'ts', 'py', 'rs', 'go', 'java', 'c', 'cpp', 'rb', 'sh', 'html', 'css'].includes(ext)) return '⌨️';
  if (['stl', 'obj', 'step', 'stp', '3mf', 'gcode'].includes(ext)) return '🧊';
  return '📦';
}

export const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

export function el(sel, root = document) {
  return root.querySelector(sel);
}

export function els(sel, root = document) {
  return [...root.querySelectorAll(sel)];
}
