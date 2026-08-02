/**
 * Codec tests — the parts that must be exactly right, since a bug here means a
 * file that silently arrives corrupt.
 *
 * Run: npm test
 */

import crypto from 'node:crypto';
import { LTEncoder, LTDecoder, blocksForSeed, degreeCDF } from '../public/js/fountain.js';
import { crc32, crc32Init, crc32Update, crc32Final, xorInto } from '../public/js/util.js';
import { buildMeta, buildData, parse, TYPE_META, TYPE_DATA, DATA_HEADER } from '../public/js/packet.js';

let failures = 0;
let checks = 0;

function ok(condition, label, detail = '') {
  checks++;
  if (condition) {
    console.log(`  ok   ${label}${detail ? `  ${detail}` : ''}`);
  } else {
    failures++;
    console.log(`  FAIL ${label}${detail ? `  ${detail}` : ''}`);
  }
}

function group(name) {
  console.log(`\n${name}`);
}

/* ------------------------------------------------------------------ crc32 */

group('crc32');
{
  // Known-answer test: CRC-32 of "123456789" is 0xCBF43926.
  const known = crc32(new TextEncoder().encode('123456789'));
  ok(known === 0xcbf43926, 'known answer vector', `got 0x${known.toString(16).toUpperCase()}`);

  const bytes = new Uint8Array(crypto.randomBytes(10_000));
  let streamed = crc32Init();
  for (let i = 0; i < bytes.length; i += 997) {
    streamed = crc32Update(streamed, bytes.subarray(i, i + 997));
  }
  ok(crc32Final(streamed) === crc32(bytes), 'streaming matches one-shot');
  ok(crc32(new Uint8Array(0)) === 0, 'empty input');
}

/* -------------------------------------------------------------------- xor */

group('xorInto');
{
  const a = new Uint8Array([1, 2, 3, 4, 5]);
  xorInto(a, new Uint8Array([1, 2, 3, 4, 5]));
  ok(a.every((b) => b === 0), 'self-xor zeroes');

  // Exercise the word-aligned fast path plus the byte tail.
  const big = new Uint8Array(crypto.randomBytes(1000));
  const copy = big.slice();
  const other = new Uint8Array(crypto.randomBytes(1000));
  xorInto(big, other);
  xorInto(big, other);
  ok(Buffer.compare(Buffer.from(big), Buffer.from(copy)) === 0, 'xor is its own inverse (1000 B)');

  // Unaligned view: must fall back to the byte loop and still be correct.
  const backing = new Uint8Array(1003);
  backing.set(other, 3);
  const unaligned = backing.subarray(3);
  const before = unaligned.slice();
  xorInto(unaligned, big);
  xorInto(unaligned, big);
  ok(Buffer.compare(Buffer.from(unaligned), Buffer.from(before)) === 0, 'unaligned view');
}

/* ----------------------------------------------------------------- packet */

group('packet framing');
{
  const meta = parse(buildMeta({
    session: 0xdeadbeef,
    fileSize: 123456,
    blockSize: 790,
    crc: 0xcafebabe,
    name: 'holiday photo (1).HEIC',
    mime: 'image/heic',
    index: 2,
    count: 5,
  }));
  ok(meta?.type === TYPE_META, 'meta type');
  ok(meta.session === 0xdeadbeef, 'session survives u32 round-trip');
  ok(meta.fileSize === 123456 && meta.blockSize === 790, 'sizes');
  ok(meta.crc === 0xcafebabe, 'crc');
  ok(meta.name === 'holiday photo (1).HEIC' && meta.mime === 'image/heic', 'name and mime');
  ok(meta.index === 2 && meta.count === 5, 'queue position');

  const payload = new Uint8Array(crypto.randomBytes(790));
  const data = parse(buildData({ session: 7, seed: 0xffffffff, payload }));
  ok(data?.type === TYPE_DATA, 'data type');
  ok(data.seed === 0xffffffff, 'max seed round-trips');
  ok(Buffer.compare(Buffer.from(data.payload), Buffer.from(payload)) === 0, 'payload is byte-exact');
  ok(buildData({ session: 7, seed: 1, payload }).length === payload.length + DATA_HEADER, 'header size');

  ok(parse(new Uint8Array([1, 2, 3])) === null, 'rejects short garbage');
  ok(parse(new Uint8Array(20)) === null, 'rejects wrong magic');
  ok(parse(new TextEncoder().encode('https://example.com/#j=ABC234')) === null, 'rejects a pairing URL');

  // Unicode filenames must survive the length-prefixed JSON block.
  const uni = parse(buildMeta({
    session: 1, fileSize: 10, blockSize: 10, crc: 0,
    name: '日本語のファイル 🎬.mkv', mime: 'video/x-matroska',
  }));
  ok(uni.name === '日本語のファイル 🎬.mkv', 'unicode filename');
}

/* --------------------------------------------------------------- fountain */

group('fountain: encoder/decoder agreement');
{
  // The decoder derives block sets from the seed alone. If the two sides ever
  // disagree, files decode to garbage that still "completes", so pin it down.
  for (const K of [1, 2, 7, 64, 1000]) {
    const cdf = degreeCDF(K);
    let match = true;
    let inRange = true;
    for (let seed = 0; seed < 300; seed++) {
      const a = blocksForSeed(seed, K, cdf);
      const b = blocksForSeed(seed, K, degreeCDF(K));
      if (a.join(',') !== b.join(',')) match = false;
      if (a.some((i) => i < 0 || i >= K) || a.length === 0 || a.length > K) inRange = false;
    }
    ok(match && inRange, `deterministic and in range (K=${K})`);
  }

  const cdf = degreeCDF(50);
  let systematic = true;
  for (let seed = 0; seed < 50; seed++) {
    const blocks = blocksForSeed(seed, 50, cdf);
    if (blocks.length !== 1 || blocks[0] !== seed) systematic = false;
  }
  ok(systematic, 'first K seeds are systematic (fast start for a fresh receiver)');
}

group('fountain: round trips');

function roundTrip(size, blockSize, lossRate, startSeed = 0) {
  const source = new Uint8Array(crypto.randomBytes(size));
  const encoder = new LTEncoder(source, blockSize);
  const decoder = new LTDecoder(encoder.K, blockSize, source.length);

  let seed = startSeed;
  let guard = 0;
  const cap = encoder.K * 100 + 10_000;
  while (!decoder.complete && guard++ < cap) {
    if (Math.random() >= lossRate) decoder.addPacket(seed >>> 0, encoder.encode(seed >>> 0));
    seed++;
  }
  if (!decoder.complete) return { ok: false, detail: `stalled at ${(decoder.progress * 100).toFixed(1)}%` };

  const out = decoder.result();
  const exact = Buffer.compare(Buffer.from(out), Buffer.from(source)) === 0;
  const overhead = ((decoder.packetsUsed / encoder.K - 1) * 100).toFixed(1);
  return {
    ok: exact && crc32(out) === crc32(source),
    detail: `K=${encoder.K} used=${decoder.packetsUsed} overhead=+${overhead}%`,
  };
}

for (const [size, block, loss, label] of [
  [1, 790, 0, '1 byte'],
  [789, 790, 0, 'just under one block'],
  [790, 790, 0, 'exactly one block'],
  [791, 790, 0, 'one block plus one byte'],
  [50 * 1024, 790, 0, '50 KB, clean channel'],
  [50 * 1024, 790, 0.3, '50 KB, 30% frame loss'],
  [50 * 1024, 790, 0.6, '50 KB, 60% frame loss'],
  [50 * 1024, 790, 0.85, '50 KB, 85% frame loss'],
  [250 * 1024, 790, 0.4, '250 KB, 40% frame loss'],
  [1024 * 1024, 1590, 0.35, '1 MB at turbo density, 35% loss'],
]) {
  const r = roundTrip(size, block, loss);
  ok(r.ok, label, r.detail);
}

// Receiver tunes in after the stream has been running: the real case, since the
// user aims the camera whenever they get around to it.
for (const startSeed of [137, 5000, 4294967000]) {
  const r = roundTrip(80 * 1024, 790, 0.3, startSeed);
  ok(r.ok, `late join at seed ${startSeed}`, r.detail);
}

group('fountain: robustness');
{
  // Duplicate packets (a looping stream re-scanned) must not corrupt state.
  const source = new Uint8Array(crypto.randomBytes(20_000));
  const encoder = new LTEncoder(source, 790);
  const decoder = new LTDecoder(encoder.K, 790, source.length);
  let seed = 0;
  while (!decoder.complete && seed < 100_000) {
    const payload = encoder.encode(seed);
    decoder.addPacket(seed, payload);
    decoder.addPacket(seed, payload); // exact duplicate
    seed++;
  }
  ok(
    decoder.complete && Buffer.compare(Buffer.from(decoder.result()), Buffer.from(source)) === 0,
    'duplicate packets are ignored safely',
    `used=${decoder.packetsUsed} of K=${encoder.K}`,
  );

  // A partially-fed decoder must refuse to hand back a half-built file rather
  // than returning zero-padded garbage.
  const partial = new LTDecoder(10, 100, 1000);
  partial.addPacket(0, new Uint8Array(100).fill(1));
  let threw = false;
  try { partial.result(); } catch { threw = true; }
  ok(threw, 'incomplete decode throws instead of returning garbage');

  // An all-zero file is a nice trap for XOR bugs.
  const zeros = new Uint8Array(5000);
  const zEnc = new LTEncoder(zeros, 790);
  const zDec = new LTDecoder(zEnc.K, 790, zeros.length);
  let s = 0;
  while (!zDec.complete && s < 50_000) { zDec.addPacket(s, zEnc.encode(s)); s++; }
  ok(zDec.complete && zDec.result().every((b) => b === 0), 'all-zero file');
}

/* --------------------------------------------------- end-to-end simulation */

group('end-to-end: optical frames through a lossy camera');
{
  // Full pipeline: encode → frame → (drop some) → parse → decode → verify,
  // including metadata discovery the way a real receiver experiences it.
  const source = new Uint8Array(crypto.randomBytes(120_000));
  const blockSize = 790;
  const session = 0x1234abcd;
  const encoder = new LTEncoder(source, blockSize);
  const fileCrc = crc32(source);

  let decoder = null;
  let seed = 0;
  let emitted = 0;
  let frames = 0;
  const cap = 60_000;

  while ((!decoder || !decoder.complete) && frames++ < cap) {
    // Mirror the sender's cadence: metadata leads and repeats every 14 frames.
    const isMeta = emitted < 2 || emitted % 14 === 0;
    const frame = isMeta
      ? buildMeta({ session, fileSize: source.length, blockSize, crc: fileCrc, name: 'clip.mp4', mime: 'video/mp4' })
      : buildData({ session, seed: seed++, payload: encoder.encode(seed - 1) });
    emitted++;

    if (Math.random() < 0.45) continue; // camera missed this frame entirely

    const packet = parse(frame);
    if (packet.type === TYPE_META) {
      if (!decoder) {
        const K = Math.ceil(packet.fileSize / packet.blockSize);
        decoder = new LTDecoder(K, packet.blockSize, packet.fileSize);
      }
    } else if (decoder) {
      decoder.addPacket(packet.seed, packet.payload);
    }
  }

  const done = decoder?.complete;
  const out = done ? decoder.result() : null;
  ok(done, 'completes through a 45%-loss channel', `frames shown=${emitted}`);
  ok(out && crc32(out) === fileCrc, 'crc matches the source');
  ok(out && Buffer.compare(Buffer.from(out), Buffer.from(source)) === 0, 'bytes are identical');
}

console.log(`\n${failures ? 'FAILED' : 'PASSED'} — ${checks - failures}/${checks} checks\n`);
process.exit(failures ? 1 : 0);
