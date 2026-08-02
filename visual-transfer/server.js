#!/usr/bin/env node
/**
 * Visual Transfer — static host + WebSocket pairing relay.
 *
 * Runs two ways, same code:
 *   local  : `npm start`  -> HTTP on :8080 and (if openssl exists) HTTPS on :8443
 *            HTTPS matters locally because getUserMedia() needs a secure context,
 *            so a phone hitting http://192.168.x.x cannot open its camera.
 *   hosted : `PORT=$PORT node server.js` -> HTTP only; the platform terminates TLS.
 *
 * The relay only ever passes WebRTC signalling blobs between two peers in a room.
 * File bytes never touch this process.
 */
'use strict';

const http = require('node:http');
const https = require('node:https');
const fs = require('node:fs');
const path = require('node:path');
const os = require('node:os');
const crypto = require('node:crypto');
const { execFileSync } = require('node:child_process');
const { WebSocketServer } = require('ws');

const ROOT = path.join(__dirname, 'public');
const PORT = Number(process.env.PORT || 8080);
const HTTPS_PORT = Number(process.env.HTTPS_PORT || 8443);
const CERT_DIR = path.join(__dirname, '.certs');
const ENABLE_HTTPS = process.env.NO_HTTPS !== '1' && !process.env.RENDER && !process.env.VERCEL;

/* ------------------------------------------------------------------ static */

const MIME = {
  '.html': 'text/html; charset=utf-8',
  '.js': 'text/javascript; charset=utf-8',
  '.css': 'text/css; charset=utf-8',
  '.json': 'application/json; charset=utf-8',
  '.svg': 'image/svg+xml',
  '.ico': 'image/x-icon',
  '.png': 'image/png',
  '.webmanifest': 'application/manifest+json',
};

function serve(req, res) {
  let urlPath;
  try {
    urlPath = decodeURIComponent(new URL(req.url, 'http://x').pathname);
  } catch {
    res.writeHead(400).end('bad request');
    return;
  }
  if (urlPath === '/' || urlPath === '') urlPath = '/index.html';
  if (urlPath === '/healthz') {
    res.writeHead(200, { 'content-type': 'text/plain' }).end('ok');
    return;
  }

  // Resolve inside ROOT only — blocks ../ traversal.
  const filePath = path.join(ROOT, path.normalize(urlPath));
  if (!filePath.startsWith(ROOT + path.sep)) {
    res.writeHead(403).end('forbidden');
    return;
  }

  fs.stat(filePath, (err, st) => {
    if (err || !st.isFile()) {
      // SPA-ish fallback so /#j=CODE deep links always land on the app.
      if (!path.extname(filePath)) return sendFile(res, path.join(ROOT, 'index.html'));
      res.writeHead(404, { 'content-type': 'text/plain' }).end('not found');
      return;
    }
    sendFile(res, filePath, st);
  });
}

function sendFile(res, filePath, st) {
  const ext = path.extname(filePath).toLowerCase();
  const headers = {
    'content-type': MIME[ext] || 'application/octet-stream',
    // The whole point is scanning a camera at a live peer — never serve stale JS.
    'cache-control': 'no-cache',
    // getUserMedia + WebRTC only; no framing.
    'x-frame-options': 'DENY',
    'referrer-policy': 'no-referrer',
  };
  if (st) headers['content-length'] = st.size;
  res.writeHead(200, headers);
  fs.createReadStream(filePath).pipe(res).on('error', () => res.destroy());
}

/* ----------------------------------------------------------------- pairing */

const CODE_ALPHABET = '23456789ABCDEFGHJKLMNPQRSTUVWXYZ'; // no 0/O/1/I/L
const CODE_LEN = 6;
const ROOM_TTL_MS = 15 * 60 * 1000;

/** @type {Map<string, {peers: Set<import('ws').WebSocket>, created: number}>} */
const rooms = new Map();

function newCode() {
  for (let attempt = 0; attempt < 50; attempt++) {
    const bytes = crypto.randomBytes(CODE_LEN);
    let code = '';
    for (let i = 0; i < CODE_LEN; i++) code += CODE_ALPHABET[bytes[i] % CODE_ALPHABET.length];
    if (!rooms.has(code)) return code;
  }
  return crypto.randomBytes(8).toString('hex').toUpperCase();
}

function send(ws, obj) {
  if (ws.readyState === ws.OPEN) ws.send(JSON.stringify(obj));
}

function leave(ws) {
  const code = ws._room;
  if (!code) return;
  ws._room = null;
  const room = rooms.get(code);
  if (!room) return;
  room.peers.delete(ws);
  for (const peer of room.peers) send(peer, { t: 'bye' });
  if (room.peers.size === 0) rooms.delete(code);
}

function attachSignalling(server) {
  const wss = new WebSocketServer({ server, path: '/ws', maxPayload: 256 * 1024 });

  wss.on('connection', (ws) => {
    ws._room = null;
    ws.isAlive = true;
    ws.on('pong', () => { ws.isAlive = true; });

    ws.on('message', (raw) => {
      let msg;
      try {
        msg = JSON.parse(raw.toString());
      } catch {
        return;
      }
      if (!msg || typeof msg.t !== 'string') return;

      if (msg.t === 'host') {
        leave(ws);
        const code = newCode();
        rooms.set(code, { peers: new Set([ws]), created: Date.now() });
        ws._room = code;
        send(ws, { t: 'hosted', code });
        return;
      }

      if (msg.t === 'join') {
        const code = String(msg.code || '').toUpperCase().trim();
        const room = rooms.get(code);
        if (!room) return send(ws, { t: 'error', m: 'That code has expired. Ask the sender for a fresh one.' });
        if (room.peers.size >= 2) return send(ws, { t: 'error', m: 'That transfer already has two devices connected.' });
        leave(ws);
        room.peers.add(ws);
        ws._room = code;
        send(ws, { t: 'joined', code });
        for (const peer of room.peers) if (peer !== ws) send(peer, { t: 'peer' });
        return;
      }

      if (msg.t === 'sig') {
        const room = rooms.get(ws._room);
        if (!room) return;
        for (const peer of room.peers) if (peer !== ws) send(peer, { t: 'sig', d: msg.d });
        return;
      }
    });

    ws.on('close', () => leave(ws));
    ws.on('error', () => leave(ws));
  });

  // Drop half-open sockets so rooms don't wedge with a ghost peer.
  const heartbeat = setInterval(() => {
    for (const ws of wss.clients) {
      if (!ws.isAlive) { ws.terminate(); continue; }
      ws.isAlive = false;
      ws.ping();
    }
    const cutoff = Date.now() - ROOM_TTL_MS;
    for (const [code, room] of rooms) {
      if (room.created < cutoff && room.peers.size < 2) {
        for (const ws of room.peers) send(ws, { t: 'error', m: 'Pairing code expired.' });
        rooms.delete(code);
      }
    }
  }, 30_000);
  wss.on('close', () => clearInterval(heartbeat));
  return wss;
}

/* -------------------------------------------------------------------- tls */

function lanAddresses() {
  const out = [];
  for (const list of Object.values(os.networkInterfaces())) {
    for (const ni of list || []) {
      if (ni.family === 'IPv4' && !ni.internal) out.push(ni.address);
    }
  }
  return out;
}

/** Self-signed cert covering localhost + this machine's LAN IPs, so phones can use the camera. */
function ensureCert() {
  const keyPath = path.join(CERT_DIR, 'key.pem');
  const certPath = path.join(CERT_DIR, 'cert.pem');
  const ips = lanAddresses();
  const wantSan = ['DNS:localhost', 'IP:127.0.0.1', ...ips.map((ip) => `IP:${ip}`)].join(',');
  const stampPath = path.join(CERT_DIR, 'san.txt');

  const fresh =
    fs.existsSync(keyPath) &&
    fs.existsSync(certPath) &&
    fs.existsSync(stampPath) &&
    fs.readFileSync(stampPath, 'utf8') === wantSan;

  if (!fresh) {
    fs.mkdirSync(CERT_DIR, { recursive: true });
    execFileSync('openssl', [
      'req', '-x509', '-newkey', 'rsa:2048', '-nodes',
      '-keyout', keyPath, '-out', certPath,
      '-days', '825', '-subj', '/CN=visual-transfer.local',
      '-addext', `subjectAltName=${wantSan}`,
    ], { stdio: 'ignore' });
    fs.writeFileSync(stampPath, wantSan);
  }
  return { key: fs.readFileSync(keyPath), cert: fs.readFileSync(certPath) };
}

/* ------------------------------------------------------------------- boot */

const httpServer = http.createServer(serve);
attachSignalling(httpServer);
httpServer.listen(PORT, () => {
  console.log(`\n  Visual Transfer`);
  console.log(`  http  →  http://localhost:${PORT}`);
});

if (ENABLE_HTTPS) {
  try {
    const creds = ensureCert();
    const httpsServer = https.createServer(creds, serve);
    attachSignalling(httpsServer);
    httpsServer.listen(HTTPS_PORT, () => {
      console.log(`  https →  https://localhost:${HTTPS_PORT}`);
      for (const ip of lanAddresses()) {
        console.log(`           https://${ip}:${HTTPS_PORT}   ← open this on your phone`);
      }
      console.log(`\n  The phone will warn about the self-signed certificate.`);
      console.log(`  Tap "Advanced" → "Proceed" once; cameras need HTTPS to turn on.\n`);
    });
  } catch (err) {
    console.log(`  https →  unavailable (${err.code === 'ENOENT' ? 'openssl not found' : err.message})`);
    console.log(`           Camera access from other devices needs HTTPS.\n`);
  }
} else {
  console.log('');
}

for (const sig of ['SIGINT', 'SIGTERM']) {
  process.on(sig, () => process.exit(0));
}
