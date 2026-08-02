/**
 * Link mode — the QR code carries a pairing URL, not the file.
 *
 * The scan is still the entire interaction, but the bytes then travel over a
 * direct WebRTC data channel at network speed instead of through the camera.
 * That is what makes real files — a 40 MB video, a folder of RAWs — practical.
 *
 * The relay server only forwards SDP and ICE. File data is peer-to-peer.
 *
 * Roles: the device showing the code hosts; the device that scans joins and makes
 * the offer. After that the channel is symmetric — either side can send.
 */

import { crc32Init, crc32Update, crc32Final, safeFilename } from './util.js';

const ICE_SERVERS = [
  { urls: 'stun:stun.l.google.com:19302' },
  { urls: 'stun:stun.cloudflare.com:3478' },
];

/** Conservative floor that every SCTP implementation accepts. */
const CHUNK_FALLBACK = 16 * 1024;
/** Ceiling regardless of what the transport claims — bigger buys nothing here. */
const CHUNK_MAX = 64 * 1024;
const HIGH_WATER = 4 * 1024 * 1024;
const LOW_WATER = 1024 * 1024;

export class Peer {
  /**
   * @param {{onState?:Function, onCode?:Function, onFile?:Function,
   *          onProgress?:Function, onError?:Function, onSendProgress?:Function}} handlers
   */
  constructor(handlers = {}) {
    this.h = handlers;
    this.ws = null;
    this.pc = null;
    this.channel = null;
    this.code = null;
    this.isHost = false;
    this.state = 'idle';
    this.closed = false;

    /** Incoming file being assembled. */
    this.incoming = null;
    this.sendQueue = [];
    this.sending = false;
  }

  #setState(state, detail) {
    if (this.closed && state !== 'closed') return;
    // Collapse no-op transitions. 'connected' legitimately arrives twice — once
    // from the connection state change and once from the channel opening — and a
    // listener that starts a transfer on it must not be told twice.
    if (this.state === state) return;
    this.state = state;
    this.h.onState?.(state, detail);
  }

  #signalUrl() {
    const proto = location.protocol === 'https:' ? 'wss:' : 'ws:';
    return `${proto}//${location.host}/ws`;
  }

  async #connectSignalling() {
    return new Promise((resolve, reject) => {
      const ws = new WebSocket(this.#signalUrl());
      this.ws = ws;
      // On a static host (GitHub Pages and friends) there is no relay to reach,
      // and there never will be — say so plainly and point at the mode that does
      // work there, rather than implying a transient network problem.
      const failed = () => reject(new Error(
        'No pairing server at this address. If this page is on static hosting, '
        + 'Link mode needs the Node server from the project — Beam mode works here as-is.',
      ));
      ws.onopen = () => resolve(ws);
      ws.onerror = failed;
      ws.onclose = () => {
        if (this.state === 'idle' || this.state === 'waiting') failed();
      };
      ws.onmessage = (event) => this.#onSignal(event);
    });
  }

  async #onSignal(event) {
    let msg;
    try {
      msg = JSON.parse(event.data);
    } catch {
      return;
    }

    if (msg.t === 'hosted') {
      this.code = msg.code;
      this.h.onCode?.(msg.code);
      this.#setState('waiting');
      return;
    }
    if (msg.t === 'joined') {
      this.code = msg.code;
      await this.#makeOffer();
      return;
    }
    if (msg.t === 'peer') {
      this.#setState('pairing');
      return;
    }
    if (msg.t === 'error') {
      this.h.onError?.(new Error(msg.m || 'Pairing failed.'));
      return;
    }
    if (msg.t === 'bye') {
      if (this.state !== 'done') this.#setState('peer-left');
      return;
    }
    if (msg.t === 'sig') {
      await this.#onRemote(msg.d);
    }
  }

  #send(obj) {
    if (this.ws?.readyState === WebSocket.OPEN) this.ws.send(JSON.stringify(obj));
  }

  #newPeerConnection() {
    const pc = new RTCPeerConnection({ iceServers: ICE_SERVERS });
    this.pc = pc;

    pc.onicecandidate = (e) => {
      if (e.candidate) this.#send({ t: 'sig', d: { ice: e.candidate } });
    };
    pc.onconnectionstatechange = () => {
      const s = pc.connectionState;
      if (s === 'connected') this.#setState('connected');
      else if (s === 'failed') this.h.onError?.(new Error('Direct connection failed. Both devices may be on networks that block peer-to-peer traffic.'));
      else if (s === 'disconnected' && this.state !== 'done') this.#setState('peer-left');
    };
    pc.ondatachannel = (e) => this.#bindChannel(e.channel);
    return pc;
  }

  #bindChannel(channel) {
    this.channel = channel;
    channel.binaryType = 'arraybuffer';
    channel.bufferedAmountLowThreshold = LOW_WATER;
    channel.onopen = () => {
      this.#setState('connected');
      // The signalling relay has done its job; drop it so no room lingers.
      if (this.ws?.readyState === WebSocket.OPEN) this.ws.close();
      this.#drain();
    };
    channel.onclose = () => {
      if (this.state !== 'done') this.#setState('peer-left');
    };
    channel.onmessage = (e) => this.#onChannelMessage(e.data);
  }

  /* ------------------------------------------------------------ handshake */

  async host() {
    await this.#connectSignalling();
    this.isHost = true;
    this.#send({ t: 'host' });
  }

  async join(code) {
    await this.#connectSignalling();
    this.isHost = false;
    this.#setState('pairing');
    this.#send({ t: 'join', code: String(code).toUpperCase().trim() });
  }

  async #makeOffer() {
    const pc = this.#newPeerConnection();
    this.#bindChannel(pc.createDataChannel('files', { ordered: true }));
    const offer = await pc.createOffer();
    await pc.setLocalDescription(offer);
    this.#send({ t: 'sig', d: { sdp: pc.localDescription } });
  }

  async #onRemote(payload) {
    if (!payload) return;
    if (payload.sdp) {
      const desc = payload.sdp;
      if (desc.type === 'offer') {
        const pc = this.pc || this.#newPeerConnection();
        await pc.setRemoteDescription(desc);
        const answer = await pc.createAnswer();
        await pc.setLocalDescription(answer);
        this.#send({ t: 'sig', d: { sdp: pc.localDescription } });
      } else if (desc.type === 'answer' && this.pc) {
        await this.pc.setRemoteDescription(desc);
      }
      return;
    }
    if (payload.ice && this.pc) {
      try {
        await this.pc.addIceCandidate(payload.ice);
      } catch {
        /* candidates can arrive before the description is set; harmless */
      }
    }
  }

  /* ---------------------------------------------------------------- send */

  get ready() {
    return this.channel?.readyState === 'open';
  }

  /** Queue files; they go out one at a time, in order. */
  enqueue(files) {
    for (const file of files) this.sendQueue.push(file);
    this.#drain();
  }

  async #drain() {
    if (this.sending || !this.ready) return;
    this.sending = true;
    try {
      while (this.sendQueue.length && this.ready) {
        await this.#sendFile(this.sendQueue.shift());
      }
    } catch (err) {
      this.h.onError?.(err);
    } finally {
      this.sending = false;
      if (this.ready && !this.sendQueue.length) this.h.onSendProgress?.({ idle: true });
    }
  }

  #waitForDrain() {
    if (this.channel.bufferedAmount <= HIGH_WATER) return Promise.resolve();
    return new Promise((resolve) => {
      const onLow = () => {
        this.channel.removeEventListener('bufferedamountlow', onLow);
        resolve();
      };
      this.channel.addEventListener('bufferedamountlow', onLow);
    });
  }

  async #sendFile(file) {
    const name = safeFilename(file.name);
    const mime = file.type || 'application/octet-stream';
    const size = file.size;
    const started = performance.now();

    this.channel.send(JSON.stringify({ t: 'begin', name, mime, size }));

    // Stream from disk rather than buffering the file — this is what lets a
    // multi-gigabyte file go through on a phone without blowing up memory.
    // Ask the transport how big a message it will actually carry instead of
    // assuming the legacy 16 KB floor; on current browsers this is far larger,
    // and bigger chunks mean proportionally less per-message overhead.
    const negotiated = this.pc?.sctp?.maxMessageSize;
    const chunkSize = Math.max(
      CHUNK_FALLBACK,
      Math.min(CHUNK_MAX, Number.isFinite(negotiated) && negotiated > 0 ? negotiated : CHUNK_FALLBACK),
    );

    const reader = file.stream().getReader();
    let sent = 0;
    let crc = crc32Init();
    let lastReport = 0;

    while (true) {
      const { done, value } = await reader.read();
      if (done) break;
      let offset = 0;
      while (offset < value.length) {
        if (!this.ready) throw new Error('Connection closed mid-transfer.');
        const slice = value.subarray(offset, offset + chunkSize);
        crc = crc32Update(crc, slice);
        // Copy: the channel takes ownership asynchronously and the view is
        // backed by a buffer the reader may reuse.
        this.channel.send(slice.slice().buffer);
        sent += slice.length;
        offset += slice.length;
        await this.#waitForDrain();

        // Report about ten times a second, not once per 16 KB chunk — at LAN
        // speeds that is hundreds of DOM updates a second and it becomes the
        // bottleneck rather than the network.
        const now = performance.now();
        if (now - lastReport > 100 || sent === size) {
          lastReport = now;
          const seconds = (now - started) / 1000;
          this.h.onSendProgress?.({
            name, size, sent,
            progress: size ? sent / size : 1,
            rate: seconds > 0 ? sent / seconds : 0,
            remaining: this.sendQueue.length,
          });
        }
      }
    }

    this.channel.send(JSON.stringify({ t: 'end', name, crc: crc32Final(crc) }));
  }

  /* ------------------------------------------------------------- receive */

  #onChannelMessage(data) {
    if (typeof data === 'string') {
      let msg;
      try {
        msg = JSON.parse(data);
      } catch {
        return;
      }
      if (msg.t === 'begin') {
        this.incoming = {
          name: safeFilename(msg.name),
          mime: msg.mime || 'application/octet-stream',
          size: Number(msg.size) || 0,
          parts: [],
          received: 0,
          crc: crc32Init(),
          started: performance.now(),
          lastReport: 0,
        };
        this.h.onProgress?.({ ...this.#snapshot(), progress: 0 });
      } else if (msg.t === 'end') {
        this.#completeIncoming(msg.crc);
      }
      return;
    }

    const file = this.incoming;
    if (!file) return;
    const bytes = new Uint8Array(data);
    file.parts.push(bytes);
    file.received += bytes.length;
    file.crc = crc32Update(file.crc, bytes);

    // Same reasoning as the send path: cap UI updates rather than emitting one
    // per arriving chunk.
    const now = performance.now();
    if (now - file.lastReport > 100 || file.received >= file.size) {
      file.lastReport = now;
      this.h.onProgress?.(this.#snapshot());
    }
  }

  #snapshot() {
    const f = this.incoming;
    const seconds = (performance.now() - f.started) / 1000;
    return {
      name: f.name,
      mime: f.mime,
      size: f.size,
      received: f.received,
      progress: f.size ? Math.min(1, f.received / f.size) : 0,
      rate: seconds > 0 ? f.received / seconds : 0,
    };
  }

  #completeIncoming(expectedCrc) {
    const f = this.incoming;
    if (!f) return;
    this.incoming = null;
    const bytes = new Uint8Array(f.received);
    let offset = 0;
    for (const part of f.parts) {
      bytes.set(part, offset);
      offset += part.length;
    }
    const seconds = (performance.now() - f.started) / 1000;
    const actual = crc32Final(f.crc);
    this.h.onFile?.({
      name: f.name,
      mime: f.mime,
      bytes,
      size: bytes.length,
      verified: expectedCrc === undefined ? true : actual === expectedCrc,
      seconds,
      rate: seconds > 0 ? bytes.length / seconds : 0,
      via: 'link',
    });
  }

  close() {
    this.closed = true;
    try { this.channel?.close(); } catch { /* already gone */ }
    try { this.pc?.close(); } catch { /* already gone */ }
    try { this.ws?.close(); } catch { /* already gone */ }
    this.#setState('closed');
  }
}
