/**
 * App shell: view routing, file staging, and the two transports.
 *
 * The seamless bits worth knowing about:
 *  - `#j=CODE` in the URL auto-joins, so scanning the Link code with the phone's
 *    own camera app is the entire receiving flow — nothing to install or press.
 *  - The in-app scanner accepts either kind of code. Point it at a pairing code
 *    and it switches to Link; point it at a data frame and it starts decoding.
 *  - Received files download on arrival and stay listed with a preview, because a
 *    silent auto-download that a browser blocks is worse than a visible card.
 */

import { OpticalSender, OpticalReceiver, DENSITY } from './optical.js';
import { Peer } from './rtc.js';
import { renderQR } from './qrdraw.js';
import {
  formatBytes, formatRate, formatDuration, fileGlyph, el, els,
} from './util.js';

const app = {
  view: 'home',
  staged: [],
  mode: 'link',
  sender: null,
  receiver: null,
  peer: null,
  objectUrls: [],
  receivedCount: 0,
};

/* ────────────────────────────────────────────────────────── plumbing */

function showView(name) {
  // Tear down anything holding a camera, a socket, or a render loop.
  if (app.view !== name) {
    if (app.view === 'beam') stopBeam();
    if (app.view === 'receive') stopScanner();
    if ((app.view === 'link' || app.view === 'join') && name !== 'join') closePeer();
  }
  app.view = name;
  for (const section of els('.view')) section.hidden = section.dataset.view !== name;
  window.scrollTo({ top: 0, behavior: 'instant' in window ? 'instant' : 'auto' });
}

function toast(message, kind = '') {
  const node = document.createElement('div');
  node.className = `toast${kind ? ` is-${kind}` : ''}`;
  node.textContent = message;
  el('#toasts').append(node);
  setTimeout(() => {
    node.style.transition = 'opacity .3s, transform .3s';
    node.style.opacity = '0';
    node.style.transform = 'translateY(6px)';
    setTimeout(() => node.remove(), 320);
  }, kind === 'bad' ? 6000 : 3600);
}

function setProgress(prefix, { name, progress, meta, show = true }) {
  const block = el(`#${prefix}Progress`);
  if (!block) return;
  block.hidden = !show;
  if (!show) return;
  if (name !== undefined) el(`#${prefix}ProgName`).textContent = name;
  if (progress !== undefined) {
    const pct = Math.round(progress * 100);
    el(`#${prefix}ProgBar`).style.width = `${pct}%`;
    el(`#${prefix}ProgPct`).textContent = `${pct}%`;
  }
  if (meta !== undefined) el(`#${prefix}ProgMeta`).textContent = meta;
}

/* ─────────────────────────────────────────────────── staging files */

function stageFiles(list) {
  const incoming = [...list].filter((f) => f && (f.size > 0 || f.type));
  if (!incoming.length) return;
  app.staged.push(...incoming);
  renderStaged();
}

function renderStaged() {
  const list = el('#fileList');
  list.innerHTML = '';
  for (const [i, file] of app.staged.entries()) {
    const row = document.createElement('li');
    row.className = 'file-row';
    row.innerHTML = `
      <span class="file-glyph">${fileGlyph(file.name, file.type)}</span>
      <span class="file-name"></span>
      <span class="file-size">${formatBytes(file.size)}</span>
      <button class="file-drop" type="button" aria-label="Remove">×</button>`;
    row.querySelector('.file-name').textContent = file.name;
    row.querySelector('.file-drop').addEventListener('click', () => {
      app.staged.splice(i, 1);
      renderStaged();
    });
    list.append(row);
  }

  const any = app.staged.length > 0;
  el('#sendConfig').hidden = !any;
  updateBeamEstimate();
}

function totalStaged() {
  return app.staged.reduce((sum, f) => sum + f.size, 0);
}

function updateBeamEstimate() {
  const density = DENSITY[el('#densitySelect').value] || DENSITY.balanced;
  const fps = Number(el('#fpsRange').value);
  const payload = density.bytes - 10;
  const bytes = totalStaged();
  if (!bytes) { el('#beamEstimate').textContent = ''; return; }

  const blocks = app.staged.reduce((sum, f) => sum + Math.ceil(f.size / payload), 0);
  const seconds = (blocks * 1.3 + app.staged.length * 4) / fps;
  const throughput = payload * fps;
  el('#beamEstimate').textContent =
    `${formatBytes(bytes)} · about ${formatDuration(seconds)} at ${formatRate(throughput)} · ${blocks} blocks`;

  if (bytes > 3 * 1024 * 1024) {
    el('#beamEstimate').textContent += ' — consider Link for something this size.';
  }
}

/* ───────────────────────────────────────────────── link (sender) ── */

async function startLink() {
  showView('link');
  el('#linkStatus').textContent = 'Opening a pairing slot…';
  el('#codeRow').hidden = true;
  setProgress('link', { show: false });

  const peer = new Peer({
    onCode: (code) => {
      const url = `${location.origin}/#j=${code}`;
      renderQR(el('#linkCanvas'), url, { ecc: 'M', quiet: 3 });
      el('#pairCode').textContent = code;
      el('#hostHint').textContent = location.host;
      el('#codeRow').hidden = false;
      el('#linkStatus').textContent = 'Scan this with the other device’s camera.';
    },
    onState: (state) => {
      if (state === 'pairing') el('#linkStatus').textContent = 'Device found — connecting…';
      if (state === 'connected') {
        el('#linkStatus').textContent = `Connected. Sending ${app.staged.length} file${app.staged.length === 1 ? '' : 's'}…`;
        el('#codeRow').hidden = true;
        peer.enqueue(app.staged);
      }
      if (state === 'peer-left') {
        el('#linkStatus').textContent = 'The other device disconnected.';
        setProgress('link', { show: false });
      }
    },
    onSendProgress: (p) => {
      if (p.idle) {
        el('#linkStatus').textContent = 'All files sent.';
        setProgress('link', { name: 'Done', progress: 1, meta: 'Every file delivered and verified.' });
        toast('All files sent.', 'good');
        return;
      }
      setProgress('link', {
        name: p.name,
        progress: p.progress,
        meta: `${formatBytes(p.sent)} of ${formatBytes(p.size)} · ${formatRate(p.rate)}${p.remaining ? ` · ${p.remaining} queued` : ''}`,
      });
    },
    // The sender can also receive, so the channel is useful in both directions.
    onProgress: (p) => setProgress('link', {
      name: `↓ ${p.name}`,
      progress: p.progress,
      meta: `${formatBytes(p.received)} of ${formatBytes(p.size)} · ${formatRate(p.rate)}`,
    }),
    onFile: (file) => acceptFile(file, 'receivedList'),
    onError: (err) => {
      el('#linkStatus').textContent = err.message;
      toast(err.message, 'bad');
    },
  });

  app.peer = peer;
  try {
    await peer.host();
  } catch (err) {
    el('#linkStatus').textContent = err.message;
    toast(err.message, 'bad');
  }
}

function closePeer() {
  app.peer?.close();
  app.peer = null;
}

/* ───────────────────────────────────────────────── link (receiver) */

async function joinCode(code) {
  showView('join');
  el('#joinStatus').textContent = 'Connecting…';
  el('#joinPulse').classList.remove('is-live');
  setProgress('join', { show: false });
  el('#reverseSend').hidden = true;

  const peer = new Peer({
    onState: (state) => {
      if (state === 'connected') {
        el('#joinStatus').textContent = 'Connected — waiting for files.';
        el('#joinPulse').classList.add('is-live');
        el('#reverseSend').hidden = false;
      }
      if (state === 'peer-left') {
        el('#joinStatus').textContent = 'The other device disconnected.';
        el('#joinPulse').classList.remove('is-live');
      }
    },
    onProgress: (p) => {
      el('#joinStatus').textContent = 'Receiving…';
      setProgress('join', {
        name: p.name,
        progress: p.progress,
        meta: `${formatBytes(p.received)} of ${formatBytes(p.size)} · ${formatRate(p.rate)}`,
      });
    },
    onSendProgress: (p) => {
      if (p.idle) { toast('Sent.', 'good'); return; }
      setProgress('join', {
        name: `↑ ${p.name}`,
        progress: p.progress,
        meta: `${formatBytes(p.sent)} of ${formatBytes(p.size)} · ${formatRate(p.rate)}`,
      });
    },
    onFile: (file) => {
      acceptFile(file, 'joinReceivedList');
      el('#joinStatus').textContent = 'Connected — waiting for files.';
      setProgress('join', { show: false });
    },
    onError: (err) => {
      el('#joinStatus').textContent = err.message;
      toast(err.message, 'bad');
    },
  });

  app.peer = peer;
  try {
    await peer.join(code);
  } catch (err) {
    el('#joinStatus').textContent = err.message;
    toast(err.message, 'bad');
  }
}

/* ──────────────────────────────────────────────────────────── beam */

async function startBeam() {
  showView('beam');
  const canvas = el('#beamCanvas');
  const density = el('#densitySelect').value;
  const fps = Number(el('#fpsRange').value);

  const sender = new OpticalSender(canvas, {
    fps,
    density,
    onFrame: ({ file, emitted, budget, pass, totalFrames }) => {
      el('#beamFile').textContent = `${file.name} (${file.index + 1}/${file.count})`;
      el('#beamFrames').textContent = String(totalFrames);
      el('#beamPass').textContent = String(pass);
      el('#beamBar').style.width = `${Math.min(100, (emitted / budget) * 100)}%`;
    },
  });
  app.sender = sender;

  el('#beamRate').textContent = `${DENSITY[density].bytes} B · v${density === 'gentle' ? 13 : density === 'balanced' ? 20 : 29}`;

  try {
    await sender.load(app.staged);
  } catch (err) {
    toast(`Could not read the files: ${err.message}`, 'bad');
    showView('send');
    return;
  }

  el('#beamHint').textContent =
    `One full pass takes about ${formatDuration(sender.estimateSeconds())}. `
    + 'The stream repeats forever, so nothing is lost if the camera looks away.';
  sender.start();
}

function stopBeam() {
  app.sender?.stop();
  app.sender = null;
}

/* ───────────────────────────────────────────────────────── scanner */

async function startScanner() {
  showView('receive');
  el('#cameraNote').textContent = '';
  el('#scanHint').textContent = 'Starting camera…';
  el('#rxProgress').hidden = true;

  const receiver = new OpticalReceiver(el('#video'), {
    onStatus: (text) => { el('#scanHint').textContent = text; },
    onProgress: (p) => {
      el('#viewport').classList.add('is-locked');
      el('#rxProgress').hidden = false;
      el('#rxName').textContent = p.count > 1 ? `${p.name} (${p.index + 1}/${p.count})` : p.name;
      const pct = Math.round(p.progress * 100);
      el('#rxBar').style.width = `${pct}%`;
      el('#rxPct').textContent = `${pct}%`;
      const rate = p.elapsed > 0 ? (p.blocks * p.size) / p.totalBlocks / p.elapsed : 0;
      el('#rxMeta').textContent =
        `${p.blocks} of ${p.totalBlocks} blocks · ${formatBytes(p.size)} · ${formatRate(rate)}`;
      el('#scanHint').textContent = `Reading ${p.name} — keep it in frame`;
    },
    onFile: (file) => {
      acceptFile(file, 'receivedList');
      el('#rxProgress').hidden = true;
      el('#viewport').classList.remove('is-locked');
      el('#scanHint').textContent = 'Looking for a code…';
    },
    onPairing: (url) => {
      // They aimed the in-app scanner at a Link code. Follow it rather than
      // making them figure out that this is the other mode.
      const code = parsePairingUrl(url);
      if (!code) return;
      toast('Pairing code found — connecting over the network.');
      stopScanner();
      joinCode(code);
    },
  });

  app.receiver = receiver;
  try {
    await receiver.start();
    el('#scanHint').textContent = 'Looking for a code…';
    await setupTorch(receiver);
  } catch (err) {
    const denied = err.name === 'NotAllowedError';
    el('#scanHint').textContent = denied ? 'Camera permission denied' : 'Camera unavailable';
    el('#cameraNote').textContent = denied
      ? 'Allow camera access in your browser’s site settings, then reload. You can still connect with a code above.'
      : `${err.message} You can still connect with a code above.`;
    toast(denied ? 'Camera permission denied.' : err.message, 'bad');
  }
}

function stopScanner() {
  app.receiver?.stop();
  app.receiver = null;
  el('#viewport').classList.remove('is-locked');
}

/** Torch helps a lot in a bright room where the sending screen is washed out. */
async function setupTorch(receiver) {
  const [track] = receiver.stream?.getVideoTracks() || [];
  const caps = track?.getCapabilities?.() || {};
  if (!caps.torch) return;
  const btn = el('#torchBtn');
  btn.hidden = false;
  let on = false;
  btn.onclick = async () => {
    on = !on;
    try {
      await track.applyConstraints({ advanced: [{ torch: on }] });
      btn.textContent = on ? 'Light off' : 'Light';
    } catch {
      btn.hidden = true;
    }
  };
}

/* ─────────────────────────────────────────────── received files ── */

function parsePairingUrl(text) {
  try {
    const url = new URL(text);
    const match = /[#&?]j=([A-Z0-9]{4,12})/i.exec(url.hash + url.search);
    if (!match) return null;
    // Only follow codes for this deployment — a scanned code is untrusted input,
    // and a foreign origin here would mean sending someone else's traffic.
    if (url.host !== location.host) return null;
    return match[1].toUpperCase();
  } catch {
    return null;
  }
}

function acceptFile(file, listId) {
  const blob = new Blob([file.bytes], { type: file.mime });
  const url = URL.createObjectURL(blob);
  app.objectUrls.push(url);
  app.receivedCount++;

  const card = document.createElement('li');
  card.className = 'rx-card';
  card.innerHTML = `
    <span class="rx-thumb">${fileGlyph(file.name, file.mime)}</span>
    <span class="rx-body">
      <span class="rx-name"></span>
      <span class="rx-meta"></span>
    </span>
    <span class="rx-actions"><a download>Save</a></span>`;

  card.querySelector('.rx-name').textContent = file.name;
  const badge = file.verified
    ? '<span class="rx-badge">✓ verified</span>'
    : '<span class="rx-badge is-bad">⚠ checksum mismatch</span>';
  card.querySelector('.rx-meta').innerHTML =
    `${formatBytes(file.size)} · ${file.via} · ${formatRate(file.rate)} · ${badge}`;

  const link = card.querySelector('a');
  link.href = url;
  link.download = file.name;

  if (file.mime.startsWith('image/')) {
    const img = document.createElement('img');
    img.src = url;
    img.alt = '';
    card.querySelector('.rx-thumb').replaceChildren(img);
  }

  const list = el(`#${listId}`);
  list.prepend(card);

  const isJoin = listId.startsWith('join');
  el(isJoin ? '#joinReceivedWrap' : '#receivedWrap').hidden = false;
  el(isJoin ? '#joinReceivedCount' : '#receivedCount').textContent = String(list.children.length);

  // Try to save straight away; browsers may block this without a gesture, which
  // is exactly why the card keeps a Save button.
  try {
    link.click();
  } catch {
    /* the visible Save button covers it */
  }

  toast(
    file.verified ? `${file.name} received and verified.` : `${file.name} received, but the checksum did not match.`,
    file.verified ? 'good' : 'bad',
  );
}

/* ──────────────────────────────────────────────────────── wiring ── */

function initSendUi() {
  const select = el('#densitySelect');
  for (const [key, preset] of Object.entries(DENSITY)) {
    const option = document.createElement('option');
    option.value = key;
    option.textContent = `${preset.label} — ${preset.note}`;
    if (key === 'balanced') option.selected = true;
    select.append(option);
  }
  select.addEventListener('change', updateBeamEstimate);

  const fps = el('#fpsRange');
  fps.addEventListener('input', () => {
    el('#fpsOut').textContent = fps.value;
    updateBeamEstimate();
  });

  for (const btn of els('.seg-btn')) {
    btn.addEventListener('click', () => {
      app.mode = btn.dataset.mode;
      for (const other of els('.seg-btn')) {
        const active = other === btn;
        other.classList.toggle('is-active', active);
        other.setAttribute('aria-checked', String(active));
      }
      el('#beamOptions').hidden = app.mode !== 'beam';
      el('#startSend').textContent = app.mode === 'beam' ? 'Start beaming' : 'Show the code';
      updateBeamEstimate();
    });
  }

  const dropzone = el('#dropzone');
  const input = el('#fileInput');
  dropzone.addEventListener('click', () => input.click());
  dropzone.addEventListener('keydown', (e) => {
    if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); input.click(); }
  });
  input.addEventListener('change', () => {
    stageFiles(input.files);
    input.value = '';
  });

  for (const type of ['dragenter', 'dragover']) {
    dropzone.addEventListener(type, (e) => {
      e.preventDefault();
      dropzone.classList.add('is-over');
    });
  }
  for (const type of ['dragleave', 'drop']) {
    dropzone.addEventListener(type, () => dropzone.classList.remove('is-over'));
  }
  dropzone.addEventListener('drop', (e) => {
    e.preventDefault();
    // Stop the bubble: the window-level fallback below would otherwise stage the
    // same drop a second time.
    e.stopPropagation();
    if (e.dataTransfer?.files?.length) stageFiles(e.dataTransfer.files);
  });

  // Dropping anywhere on the send view should work, not just on the target.
  window.addEventListener('dragover', (e) => e.preventDefault());
  window.addEventListener('drop', (e) => {
    if (app.view !== 'send') return;
    e.preventDefault();
    if (e.dataTransfer?.files?.length) stageFiles(e.dataTransfer.files);
  });

  window.addEventListener('paste', (e) => {
    if (app.view !== 'send' && app.view !== 'home') return;
    const files = [...(e.clipboardData?.files || [])];
    if (!files.length) return;
    if (app.view === 'home') showView('send');
    stageFiles(files);
    toast(`${files.length} file${files.length === 1 ? '' : 's'} pasted.`);
  });

  el('#startSend').addEventListener('click', () => {
    if (!app.staged.length) return toast('Add a file first.', 'bad');
    if (app.mode === 'beam') startBeam();
    else startLink();
  });
}

function initBeamControls() {
  el('#beamNext').addEventListener('click', () => app.sender?.next());

  const pause = el('#beamPause');
  pause.addEventListener('click', () => {
    const sender = app.sender;
    if (!sender) return;
    if (sender.running) { sender.stop(); pause.textContent = 'Resume'; }
    else { sender.start(); pause.textContent = 'Pause'; }
  });

  el('#beamFullscreen').addEventListener('click', async () => {
    const stage = el('#beamStage');
    try {
      if (document.fullscreenElement) await document.exitFullscreen();
      else await stage.requestFullscreen();
    } catch {
      toast('Fullscreen was refused by the browser.', 'bad');
    }
  });
}

function initReceiveUi() {
  const join = () => {
    const code = el('#manualCode').value.trim().toUpperCase();
    if (code.length < 4) return toast('Enter the 6-character code.', 'bad');
    stopScanner();
    joinCode(code);
  };
  el('#manualJoin').addEventListener('click', join);
  el('#manualCode').addEventListener('keydown', (e) => { if (e.key === 'Enter') join(); });

  const reverseInput = el('#reverseInput');
  el('#reverseBtn').addEventListener('click', () => reverseInput.click());
  reverseInput.addEventListener('change', () => {
    if (reverseInput.files.length && app.peer?.ready) {
      app.peer.enqueue([...reverseInput.files]);
      toast(`Sending ${reverseInput.files.length} file${reverseInput.files.length === 1 ? '' : 's'} back.`);
    }
    reverseInput.value = '';
  });
}

function initShell() {
  el('#brand').addEventListener('click', () => { location.hash = ''; showView('home'); });
  el('#goSend').addEventListener('click', () => showView('send'));
  el('#goReceive').addEventListener('click', () => startScanner());
  for (const btn of els('[data-back]')) {
    btn.addEventListener('click', () => {
      if (location.hash) location.hash = '';
      showView(app.staged.length && app.view !== 'receive' ? 'send' : 'home');
    });
  }

  const sheet = el('#helpSheet');
  el('#helpBtn').addEventListener('click', () => sheet.showModal());
  el('#helpClose').addEventListener('click', () => sheet.close());

  // Warn early if this page can't use a camera at all — far better than a
  // confusing permission failure after the user has picked files.
  if (!window.isSecureContext) {
    const badge = el('#netBadge');
    badge.hidden = false;
    badge.textContent = 'insecure — no camera';
    badge.title = 'Browsers only allow camera access on https:// or localhost. Link mode still works.';
  }
}

function route() {
  const match = /[#&]j=([A-Z0-9]{4,12})/i.exec(location.hash);
  if (match) {
    joinCode(match[1].toUpperCase());
    return true;
  }
  return false;
}

window.addEventListener('hashchange', () => {
  if (!route() && (app.view === 'join')) showView('home');
});

window.addEventListener('pagehide', () => {
  stopBeam();
  stopScanner();
  closePeer();
  for (const url of app.objectUrls) URL.revokeObjectURL(url);
});

initShell();
initSendUi();
initBeamControls();
initReceiveUi();
route();
