/**
 * Optical wire format.
 *
 * Every frame is self-describing binary — no state is assumed between frames,
 * because the receiver may start watching at any moment and may miss any frame.
 *
 * META (repeated periodically so a late receiver learns the file quickly)
 *   0      u8   magic 0x56 'V'
 *   1      u8   type 1
 *   2..5   u32  session id
 *   6..9   u32  file size in bytes
 *   10..11 u16  block size
 *   12..15 u32  crc32 of the file
 *   16..17 u16  json length
 *   18..   utf8 json  { n: name, t: mime, i: index, c: count }
 *
 * DATA
 *   0      u8   magic 0x56
 *   1      u8   type 2
 *   2..5   u32  session id
 *   6..9   u32  seed  (defines which blocks are XORed — see fountain.js)
 *   10..   payload, exactly blockSize bytes
 */

export const MAGIC = 0x56;
export const TYPE_META = 1;
export const TYPE_DATA = 2;
export const DATA_HEADER = 10;
export const META_HEADER = 18;

const encoder = new TextEncoder();
const decoder = new TextDecoder();

export function buildMeta({ session, fileSize, blockSize, crc, name, mime, index = 0, count = 1 }) {
  const json = encoder.encode(JSON.stringify({ n: name, t: mime || '', i: index, c: count }));
  const buf = new Uint8Array(META_HEADER + json.length);
  const view = new DataView(buf.buffer);
  buf[0] = MAGIC;
  buf[1] = TYPE_META;
  view.setUint32(2, session, true);
  view.setUint32(6, fileSize, true);
  view.setUint16(10, blockSize, true);
  view.setUint32(12, crc, true);
  view.setUint16(16, json.length, true);
  buf.set(json, META_HEADER);
  return buf;
}

export function buildData({ session, seed, payload }) {
  const buf = new Uint8Array(DATA_HEADER + payload.length);
  const view = new DataView(buf.buffer);
  buf[0] = MAGIC;
  buf[1] = TYPE_DATA;
  view.setUint32(2, session, true);
  view.setUint32(6, seed, true);
  buf.set(payload, DATA_HEADER);
  return buf;
}

/** Parse a scanned frame. Returns null for anything that isn't ours. */
export function parse(bytes) {
  if (!bytes || bytes.length < DATA_HEADER || bytes[0] !== MAGIC) return null;
  const arr = bytes instanceof Uint8Array ? bytes : new Uint8Array(bytes);
  const view = new DataView(arr.buffer, arr.byteOffset, arr.byteLength);
  const type = arr[1];
  const session = view.getUint32(2, true);

  if (type === TYPE_META) {
    if (arr.length < META_HEADER) return null;
    const fileSize = view.getUint32(6, true);
    const blockSize = view.getUint16(10, true);
    const crc = view.getUint32(12, true);
    const jsonLen = view.getUint16(16, true);
    if (blockSize === 0 || arr.length < META_HEADER + jsonLen) return null;
    let info;
    try {
      info = JSON.parse(decoder.decode(arr.subarray(META_HEADER, META_HEADER + jsonLen)));
    } catch {
      return null;
    }
    return {
      type: TYPE_META,
      session,
      fileSize,
      blockSize,
      crc,
      name: info.n || 'file',
      mime: info.t || '',
      index: info.i ?? 0,
      count: info.c ?? 1,
    };
  }

  if (type === TYPE_DATA) {
    return {
      type: TYPE_DATA,
      session,
      seed: view.getUint32(6, true),
      payload: arr.subarray(DATA_HEADER),
    };
  }

  return null;
}
