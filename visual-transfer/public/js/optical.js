/**
 * Beam mode — moving bytes with nothing but a screen and a camera.
 *
 * Sender: renders an endless stream of fountain-coded QR frames.
 * Receiver: watches with the camera and rebuilds whatever it sees.
 *
 * There is no back-channel, so the sender never learns what arrived. It handles
 * that by looping forever and rotating through the queue: anything a receiver
 * missed on one pass lands on the next. The receiver decodes every session it
 * sees in parallel and ignores sessions it has already finished, so the user
 * experience is "point the camera and wait" rather than a handshake.
 */

import { LTEncoder, LTDecoder } from './fountain.js';
import { buildMeta, buildData, parse, TYPE_META, TYPE_DATA, DATA_HEADER } from './packet.js';
import { crc32, safeFilename } from './util.js';
import { renderQR } from './qrdraw.js';

/**
 * Frame capacity presets. Bigger payload = fewer frames, but a denser code needs
 * a steadier hand and a better camera. Sizes are chosen to land just under a QR
 * version boundary at error-correction level L.
 */
export const DENSITY = {
  gentle: { label: 'Gentle', bytes: 400, note: 'v13 · older phones, shaky hands' },
  balanced: { label: 'Balanced', bytes: 800, note: 'v20 · recommended' },
  turbo: { label: 'Turbo', bytes: 1600, note: 'v29 · steady mount, good camera' },
};

const QUIET_ZONE = 3; // modules of white margin — scanners need this to lock on

function randomSession() {
  return crypto.getRandomValues(new Uint32Array(1))[0];
}

/* ------------------------------------------------------------------ sender */

export class OpticalSender {
  /**
   * @param {HTMLCanvasElement} canvas
   * @param {{fps?: number, density?: keyof typeof DENSITY, onFrame?: Function}} opts
   */
  constructor(canvas, opts = {}) {
    this.canvas = canvas;
    this.ctx = canvas.getContext('2d', { alpha: false });
    this.fps = opts.fps || 10;
    this.densityKey = opts.density || 'balanced';
    this.onFrame = opts.onFrame || (() => {});

    /** @type {Array<{session:number, name:string, mime:string, size:number, crc:number, encoder:LTEncoder, budget:number}>} */
    this.queue = [];
    this.cursor = 0;
    this.emitted = 0; // frames emitted for the current file this pass
    this.seed = 0;
    this.pass = 1;
    this.totalFrames = 0;
    this.running = false;
    this.timer = null;
  }

  get frameBytes() {
    return DENSITY[this.densityKey].bytes;
  }

  get blockSize() {
    return this.frameBytes - DATA_HEADER;
  }

  get current() {
    return this.queue[this.cursor] || null;
  }

  /** Read files into fountain encoders. Call before start(). */
  async load(files) {
    const blockSize = this.blockSize;
    this.queue = [];
    for (let i = 0; i < files.length; i++) {
      const file = files[i];
      const bytes = new Uint8Array(await file.arrayBuffer());
      const encoder = new LTEncoder(bytes, blockSize);
      this.queue.push({
        session: randomSession(),
        name: safeFilename(file.name),
        mime: file.type || 'application/octet-stream',
        size: bytes.length,
        crc: crc32(bytes),
        index: i,
        count: files.length,
        encoder,
        // Emit ~2x the minimum before rotating away: enough that a receiver with
        // a mediocre view still finishes on the first pass.
        budget: Math.ceil(encoder.K * 2) + 10,
      });
    }
    this.cursor = 0;
    this.emitted = 0;
    this.seed = 0;
    return this.queue;
  }

  /** Rough single-pass estimate, assuming the receiver catches most frames. */
  estimateSeconds() {
    const frames = this.queue.reduce((sum, item) => sum + Math.ceil(item.encoder.K * 1.3) + 4, 0);
    return frames / this.fps;
  }

  start() {
    if (this.running || !this.queue.length) return;
    this.running = true;
    const tick = () => {
      if (!this.running) return;
      this.#renderNext();
      this.timer = setTimeout(tick, 1000 / this.fps);
    };
    tick();
  }

  stop() {
    this.running = false;
    if (this.timer) clearTimeout(this.timer);
    this.timer = null;
  }

  /** Jump to the next queued file immediately. */
  next() {
    this.#advance();
    if (this.running) this.#renderNext();
  }

  #advance() {
    this.emitted = 0;
    this.cursor++;
    if (this.cursor >= this.queue.length) {
      this.cursor = 0;
      this.pass++;
    }
  }

  #renderNext() {
    const item = this.current;
    if (!item) return;

    // Lead each file with metadata, then repeat it periodically so a receiver
    // that tunes in mid-file still learns the name, size and checksum fast.
    const isMeta = this.emitted < 2 || this.emitted % 14 === 0;
    let frame;
    if (isMeta) {
      frame = buildMeta({
        session: item.session,
        fileSize: item.size,
        blockSize: this.blockSize,
        crc: item.crc,
        name: item.name,
        mime: item.mime,
        index: item.index,
        count: item.count,
      });
    } else {
      const payload = item.encoder.encode(this.seed >>> 0);
      frame = buildData({ session: item.session, seed: this.seed >>> 0, payload });
      this.seed = (this.seed + 1) >>> 0;
    }

    this.#draw(frame);
    this.emitted++;
    this.totalFrames++;
    this.onFrame({
      file: item,
      emitted: this.emitted,
      budget: item.budget,
      pass: this.pass,
      totalFrames: this.totalFrames,
      isMeta,
    });

    if (this.emitted >= item.budget) {
      this.#advance();
      this.seed = 0; // fresh file: start with the systematic pass again
    }
  }

  #draw(bytes) {
    // Byte mode keeps the payload binary-clean — no base64 inflation.
    // ECC level L: the fountain code above already handles loss far better than
    // per-frame error correction can, so spend the symbol on payload instead.
    const info = renderQR(this.canvas, bytes, { ecc: 'L', quiet: QUIET_ZONE });
    this.lastVersion = info.version;
    return info.version;
  }
}

/* ---------------------------------------------------------------- receiver */

export class OpticalReceiver {
  /**
   * @param {HTMLVideoElement} video
   * @param {{onFile?:Function, onProgress?:Function, onPairing?:Function, onStatus?:Function}} handlers
   */
  constructor(video, handlers = {}) {
    this.video = video;
    this.handlers = handlers;
    this.worker = null;
    this.stream = null;
    this.canvas = document.createElement('canvas');
    this.ctx = this.canvas.getContext('2d', { alpha: false, willReadFrequently: true });
    this.scanning = false;
    this.inFlight = false;
    this.frameId = 0;
    this.pool = null; // reused pixel buffer, ping-ponged with the worker

    /** @type {Map<number, {decoder:LTDecoder, meta:object, started:number}>} */
    this.sessions = new Map();
    this.finished = new Set();
    this.framesSeen = 0;
    this.decodesSeen = 0;
  }

  async start() {
    if (!navigator.mediaDevices?.getUserMedia) {
      throw new Error(
        window.isSecureContext
          ? 'This browser does not expose a camera to web pages.'
          : 'Camera access needs a secure page (https:// or localhost).',
      );
    }

    this.stream = await navigator.mediaDevices.getUserMedia({
      video: {
        facingMode: { ideal: 'environment' },
        width: { ideal: 1920 },
        height: { ideal: 1080 },
      },
      audio: false,
    });

    // Continuous autofocus if the platform allows it; screens are close subjects
    // and a fixed focus hunts badly on dense codes.
    const [track] = this.stream.getVideoTracks();
    try {
      const caps = track.getCapabilities?.() || {};
      const advanced = [];
      if (caps.focusMode?.includes('continuous')) advanced.push({ focusMode: 'continuous' });
      if (advanced.length) await track.applyConstraints({ advanced });
    } catch {
      /* constraint support is patchy; the defaults are fine */
    }

    this.video.srcObject = this.stream;
    this.video.playsInline = true;
    this.video.muted = true;
    await this.video.play();

    this.worker = new Worker('js/scanworker.js');
    this.worker.onmessage = (e) => this.#onWorkerMessage(e.data);

    this.scanning = true;
    this.#pump();
    return this.stream;
  }

  stop() {
    this.scanning = false;
    if (this.worker) { this.worker.terminate(); this.worker = null; }
    if (this.stream) { for (const t of this.stream.getTracks()) t.stop(); this.stream = null; }
    if (this.video) this.video.srcObject = null;
  }

  /** Grab a frame and hand the pixels to the worker; one decode in flight at a time. */
  #pump() {
    if (!this.scanning) return;
    const again = () => {
      if (this.scanning) requestAnimationFrame(() => this.#pump());
    };
    if (this.inFlight || this.video.readyState < 2) return again();

    const vw = this.video.videoWidth;
    const vh = this.video.videoHeight;
    if (!vw || !vh) return again();

    // Cap the working resolution: above ~1280px wide the decode cost climbs
    // faster than the accuracy does.
    const scale = Math.min(1, 1280 / vw);
    const w = Math.round(vw * scale);
    const h = Math.round(vh * scale);
    if (this.canvas.width !== w || this.canvas.height !== h) {
      this.canvas.width = w;
      this.canvas.height = h;
      this.pool = null;
    }

    this.ctx.drawImage(this.video, 0, 0, w, h);
    const image = this.ctx.getImageData(0, 0, w, h);
    this.framesSeen++;

    let buffer;
    if (this.pool && this.pool.byteLength === image.data.byteLength) {
      new Uint8ClampedArray(this.pool).set(image.data);
      buffer = this.pool;
      this.pool = null;
    } else {
      buffer = image.data.buffer;
    }

    this.inFlight = true;
    this.worker.postMessage({ buffer, width: w, height: h, id: ++this.frameId }, [buffer]);
    again();
  }

  #onWorkerMessage(msg) {
    this.inFlight = false;
    if (msg.buffer) this.pool = msg.buffer;
    if (msg.found) {
      this.decodesSeen++;
      this.#ingest(msg.bytes, msg.text);
    }
    if (this.scanning) this.#pump();
  }

  #ingest(bytes, text) {
    const packet = parse(bytes);
    if (!packet) {
      // Not a Beam frame — could be a pairing link from Link mode, which means
      // the user pointed the scanner at the other kind of code. Hand it up so the
      // app can switch modes instead of showing nothing.
      if (text && /^https?:\/\//i.test(text)) this.handlers.onPairing?.(text);
      return;
    }
    if (this.finished.has(packet.session)) return;

    if (packet.type === TYPE_META) {
      if (!this.sessions.has(packet.session)) {
        if (packet.fileSize === 0 || packet.blockSize === 0) return;
        const K = Math.max(1, Math.ceil(packet.fileSize / packet.blockSize));
        this.sessions.set(packet.session, {
          decoder: new LTDecoder(K, packet.blockSize, packet.fileSize),
          meta: packet,
          started: performance.now(),
        });
        this.handlers.onStatus?.(`Receiving ${packet.name}`);
      }
      this.#report(packet.session);
      return;
    }

    if (packet.type === TYPE_DATA) {
      const session = this.sessions.get(packet.session);
      // Data before metadata: we can't size the file yet, so drop it. Metadata
      // repeats every 14 frames, so this resolves within about a second.
      if (!session) return;
      if (packet.payload.length < session.decoder.blockSize) return;

      const done = session.decoder.addPacket(packet.seed, packet.payload);
      this.#report(packet.session);
      if (done) this.#finish(packet.session);
    }
  }

  #report(sessionId) {
    const session = this.sessions.get(sessionId);
    if (!session) return;
    this.handlers.onProgress?.({
      session: sessionId,
      name: session.meta.name,
      mime: session.meta.mime,
      size: session.meta.fileSize,
      index: session.meta.index,
      count: session.meta.count,
      progress: session.decoder.progress,
      blocks: session.decoder.solvedCount,
      totalBlocks: session.decoder.K,
      elapsed: (performance.now() - session.started) / 1000,
    });
  }

  #finish(sessionId) {
    const session = this.sessions.get(sessionId);
    if (!session) return;
    this.sessions.delete(sessionId);
    this.finished.add(sessionId);

    const bytes = session.decoder.result();
    const actual = crc32(bytes);
    const verified = actual === session.meta.crc;
    const seconds = (performance.now() - session.started) / 1000;

    this.handlers.onFile?.({
      name: session.meta.name,
      mime: session.meta.mime || 'application/octet-stream',
      bytes,
      size: bytes.length,
      verified,
      seconds,
      rate: seconds > 0 ? bytes.length / seconds : 0,
      via: 'beam',
    });
  }
}
