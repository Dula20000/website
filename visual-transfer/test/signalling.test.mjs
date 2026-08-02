/**
 * Pairing relay tests. Boots the real server on a scratch port and drives it with
 * WebSocket clients.
 *
 * Run: node test/signalling.test.mjs
 */

import { spawn } from 'node:child_process';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import WebSocket from 'ws';

const here = path.dirname(fileURLToPath(import.meta.url));
const PORT = 8791;
const BASE = `ws://127.0.0.1:${PORT}/ws`;

let failures = 0;
let checks = 0;

function ok(condition, label, detail = '') {
  checks++;
  console.log(`  ${condition ? 'ok  ' : 'FAIL'} ${label}${detail ? `  ${detail}` : ''}`);
  if (!condition) failures++;
}

const open = () => new Promise((resolve, reject) => {
  const ws = new WebSocket(BASE);
  ws.queue = [];
  ws.waiters = [];
  ws.on('message', (raw) => {
    const msg = JSON.parse(raw.toString());
    const waiter = ws.waiters.shift();
    if (waiter) waiter(msg);
    else ws.queue.push(msg);
  });
  ws.next = () => new Promise((res, rej) => {
    if (ws.queue.length) return res(ws.queue.shift());
    const timer = setTimeout(() => rej(new Error('timed out waiting for a message')), 4000);
    ws.waiters.push((m) => { clearTimeout(timer); res(m); });
  });
  ws.on('open', () => resolve(ws));
  ws.on('error', reject);
});

const server = spawn(process.execPath, [path.join(here, '..', 'server.js')], {
  env: { ...process.env, PORT: String(PORT), NO_HTTPS: '1' },
  stdio: ['ignore', 'ignore', 'inherit'],
});

async function waitForServer() {
  for (let i = 0; i < 60; i++) {
    try {
      const res = await fetch(`http://127.0.0.1:${PORT}/healthz`);
      if (res.ok) return;
    } catch { /* not up yet */ }
    await new Promise((r) => setTimeout(r, 100));
  }
  throw new Error('server never came up');
}

try {
  await waitForServer();
  console.log('\npairing relay');

  // Static hosting and traversal guard.
  const index = await fetch(`http://127.0.0.1:${PORT}/`);
  ok(index.ok && (await index.text()).includes('Visual Transfer'), 'serves the app');
  const escape = await fetch(`http://127.0.0.1:${PORT}/../server.js`);
  ok(escape.status === 404 || escape.status === 403, 'refuses to serve outside public/', `status ${escape.status}`);

  // Host assigns a code.
  const host = await open();
  host.send(JSON.stringify({ t: 'host' }));
  const hosted = await host.next();
  ok(hosted.t === 'hosted' && /^[2-9A-HJ-NP-Z]{6}$/.test(hosted.code), 'host gets a 6-char code', hosted.code);

  // Join notifies the host and confirms to the joiner.
  const guest = await open();
  guest.send(JSON.stringify({ t: 'join', code: hosted.code }));
  const [joined, notified] = await Promise.all([guest.next(), host.next()]);
  ok(joined.t === 'joined', 'joiner is admitted');
  ok(notified.t === 'peer', 'host is told a peer arrived');

  // Codes are case-insensitive: users retype them from a screen.
  const host2 = await open();
  host2.send(JSON.stringify({ t: 'host' }));
  const hosted2 = await host2.next();
  const guest2 = await open();
  guest2.send(JSON.stringify({ t: 'join', code: hosted2.code.toLowerCase() }));
  ok((await guest2.next()).t === 'joined', 'lowercase code works');
  await host2.next();

  // Signalling is relayed to the other peer only.
  guest.send(JSON.stringify({ t: 'sig', d: { sdp: { type: 'offer', sdp: 'v=0 test' } } }));
  const relayed = await host.next();
  ok(relayed.t === 'sig' && relayed.d.sdp.sdp === 'v=0 test', 'relays signalling to the peer');
  host.send(JSON.stringify({ t: 'sig', d: { ice: 'candidate:1' } }));
  ok((await guest.next()).d.ice === 'candidate:1', 'relays in both directions');

  // A third device cannot muscle into an occupied room.
  const third = await open();
  third.send(JSON.stringify({ t: 'join', code: hosted.code }));
  const rejected = await third.next();
  ok(rejected.t === 'error' && /two devices/i.test(rejected.m), 'rejects a third device');

  // Unknown codes fail cleanly rather than opening an empty room.
  third.send(JSON.stringify({ t: 'join', code: 'ZZZZZZ' }));
  const unknown = await third.next();
  ok(unknown.t === 'error' && /expired/i.test(unknown.m), 'unknown code is rejected');

  // Disconnect tells the other side, so the UI can stop pretending.
  guest.close();
  ok((await host.next()).t === 'bye', 'peer departure is announced');

  // Malformed input must not take the relay down.
  const noisy = await open();
  noisy.send('not json at all');
  noisy.send(JSON.stringify({ nope: true }));
  noisy.send(JSON.stringify({ t: 'sig', d: { x: 1 } })); // no room
  noisy.send(JSON.stringify({ t: 'host' }));
  ok((await noisy.next()).t === 'hosted', 'survives malformed messages');

  for (const ws of [host, host2, guest2, third, noisy]) ws.close();
} catch (err) {
  failures++;
  console.log(`  FAIL harness: ${err.message}`);
} finally {
  server.kill('SIGKILL');
}

console.log(`\n${failures ? 'FAILED' : 'PASSED'} — ${checks - failures}/${checks} checks\n`);
process.exit(failures ? 1 : 0);
