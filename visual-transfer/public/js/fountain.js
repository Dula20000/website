/**
 * Luby-Transform (fountain) coding.
 *
 * Why this and not "frame 1 of 340": a camera watching a screen drops frames
 * unpredictably — motion blur, autofocus hunting, rolling shutter, a hand moving.
 * An ordered scheme has to detect every gap and ask for a resend, which needs a
 * back-channel the optical path does not have.
 *
 * With an LT code the sender emits an endless stream of XOR combinations of the
 * file's blocks. The receiver collects whatever it happens to catch, in any order,
 * and can rebuild the file once it has slightly more than K useful packets.
 * Nothing is ever "missed" — it just keeps watching. That is what makes the
 * receiver feel like a camera rather than a protocol.
 *
 * Both sides derive a packet's block set from its 32-bit seed alone, so a packet
 * is fully self-describing in 4 bytes of header.
 */

import { xorInto } from './util.js';

/** mulberry32 — tiny, fast, and identical across JS engines (important: both sides must agree). */
function mulberry32(seed) {
  let a = seed >>> 0;
  return function next() {
    a = (a + 0x6d2b79f5) | 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/**
 * Robust soliton distribution as a CDF over degrees 1..K.
 * The ideal soliton is optimal in expectation but fragile in practice; the robust
 * variant adds a spike of low-degree packets so decoding actually gets started
 * and rarely stalls.
 */
export function degreeCDF(K, c = 0.03, delta = 0.5) {
  const cdf = new Float64Array(K + 1);
  if (K <= 1) {
    if (K === 1) cdf[1] = 1;
    return cdf;
  }
  const rho = new Float64Array(K + 1);
  rho[1] = 1 / K;
  for (let d = 2; d <= K; d++) rho[d] = 1 / (d * (d - 1));

  const tau = new Float64Array(K + 1);
  const R = c * Math.log(K / delta) * Math.sqrt(K);
  const spike = Math.floor(K / R);
  if (spike >= 1 && spike <= K) {
    for (let i = 1; i <= spike - 1; i++) tau[i] = R / (i * K);
    tau[spike] = (R * Math.log(R / delta)) / K;
  }

  let beta = 0;
  for (let d = 1; d <= K; d++) beta += rho[d] + tau[d];

  let acc = 0;
  for (let d = 1; d <= K; d++) {
    acc += (rho[d] + tau[d]) / beta;
    cdf[d] = acc;
  }
  cdf[K] = 1;
  return cdf;
}

function sampleDegree(cdf, u, K) {
  let lo = 1;
  let hi = K;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (cdf[mid] >= u) hi = mid;
    else lo = mid + 1;
  }
  return Math.min(Math.max(lo, 1), K);
}

/**
 * The block indices a packet combines. Pure function of (seed, K) — this is the
 * shared contract between encoder and decoder.
 *
 * Seeds below K are systematic (one block, plainly). That means a fresh receiver
 * gets usable blocks immediately instead of waiting for the peeling decoder to
 * find a foothold, and it makes the very common small-file case near-optimal.
 */
export function blocksForSeed(seed, K, cdf) {
  if (K <= 1) return [0];
  if (seed < K) return [seed];

  const rand = mulberry32((seed ^ 0x9e3779b9) >>> 0);
  const degree = sampleDegree(cdf, rand(), K);
  if (degree >= K) {
    const all = new Array(K);
    for (let i = 0; i < K; i++) all[i] = i;
    return all;
  }
  const picked = new Set();
  let guard = 0;
  while (picked.size < degree && guard++ < degree * 40) {
    picked.add(Math.floor(rand() * K) % K);
  }
  return [...picked];
}

export class LTEncoder {
  /** @param {Uint8Array} bytes @param {number} blockSize */
  constructor(bytes, blockSize) {
    this.blockSize = blockSize;
    this.byteLength = bytes.length;
    this.K = Math.max(1, Math.ceil(bytes.length / blockSize));
    this.cdf = degreeCDF(this.K);

    // One flat padded buffer; blocks are views into it so there's no per-block copy.
    this.padded = new Uint8Array(this.K * blockSize);
    this.padded.set(bytes);
    this.blocks = new Array(this.K);
    for (let i = 0; i < this.K; i++) {
      this.blocks[i] = this.padded.subarray(i * blockSize, (i + 1) * blockSize);
    }
    this.scratch = new Uint8Array(blockSize);
  }

  /**
   * Encoded payload for a seed. Returns a reused scratch buffer — copy it if you
   * need to hold on to it past the next call.
   */
  encode(seed) {
    const indices = blocksForSeed(seed, this.K, this.cdf);
    const out = this.scratch;
    out.set(this.blocks[indices[0]]);
    for (let i = 1; i < indices.length; i++) xorInto(out, this.blocks[indices[i]]);
    return out;
  }
}

export class LTDecoder {
  /** @param {number} K @param {number} blockSize @param {number} byteLength */
  constructor(K, blockSize, byteLength) {
    this.K = K;
    this.blockSize = blockSize;
    this.byteLength = byteLength;
    this.cdf = degreeCDF(K);

    this.solved = new Array(K).fill(null);
    this.solvedCount = 0;
    this.seen = new Set();
    this.packetsUsed = 0;
    this.complete = false;

    /** Unsolved combinations, plus an index → combinations map so peeling is O(1) per hit. */
    this.pending = [];
    this.byIndex = new Array(K);
    this.readyQueue = [];
  }

  get progress() {
    return this.K === 0 ? 1 : this.solvedCount / this.K;
  }

  /**
   * Feed one packet. Returns true once the whole file is recovered.
   * Duplicate seeds are ignored, so re-scanning a looping stream is free.
   */
  addPacket(seed, payload) {
    if (this.complete || this.seen.has(seed)) return this.complete;
    this.seen.add(seed);

    const indices = blocksForSeed(seed, this.K, this.cdf);
    const data = new Uint8Array(this.blockSize);
    data.set(payload.subarray(0, this.blockSize));

    // Cancel out every block we already know; what's left is the new information.
    const unknown = [];
    for (const i of indices) {
      const known = this.solved[i];
      if (known) xorInto(data, known);
      else unknown.push(i);
    }
    if (unknown.length === 0) return this.complete; // redundant packet

    this.packetsUsed++;
    this.#enqueue(unknown, data);
    this.#peel();
    return this.complete;
  }

  #enqueue(indices, data) {
    if (indices.length === 1) {
      this.readyQueue.push([indices[0], data]);
      return;
    }
    const entry = { indices: new Set(indices), data, dead: false };
    this.pending.push(entry);
    for (const i of indices) {
      (this.byIndex[i] ||= []).push(entry);
    }
  }

  /** Classic peeling: each newly solved block may reduce others to degree 1, cascading. */
  #peel() {
    while (this.readyQueue.length) {
      const [index, data] = this.readyQueue.pop();
      if (this.solved[index]) continue;

      this.solved[index] = data;
      this.solvedCount++;
      if (this.solvedCount === this.K) {
        this.complete = true;
        this.pending.length = 0;
        this.readyQueue.length = 0;
        return;
      }

      const dependents = this.byIndex[index];
      if (!dependents) continue;
      this.byIndex[index] = null;
      for (const entry of dependents) {
        if (entry.dead || !entry.indices.has(index)) continue;
        entry.indices.delete(index);
        xorInto(entry.data, data);
        if (entry.indices.size === 1) {
          entry.dead = true;
          this.readyQueue.push([entry.indices.values().next().value, entry.data]);
        } else if (entry.indices.size === 0) {
          entry.dead = true;
        }
      }
    }
  }

  /** Reassembled bytes, trimmed back to the real file length. */
  result() {
    const out = new Uint8Array(this.byteLength);
    for (let i = 0; i < this.K; i++) {
      const block = this.solved[i];
      if (!block) throw new Error('decode incomplete');
      const offset = i * this.blockSize;
      const take = Math.min(this.blockSize, this.byteLength - offset);
      if (take <= 0) break;
      out.set(block.subarray(0, take), offset);
    }
    return out;
  }
}
