/* Clawd Run — a pixelated offline-dino homage starring Clawd.
   All art is generated from pixel maps below; no external assets. */
(() => {
'use strict';

const W = 320, H = 180, GROUND = 150;
const BG = '#f0e9dc';

const cv = document.getElementById('game');
const ctx = cv.getContext('2d');
ctx.imageSmoothingEnabled = false;

/* ---------- palette ---------- */
const PAL = {
  K: '#2a2019', // ink
  W: '#ffffff',
  C: '#faf6ec', // cloud
  y: '#f2c14e', // coin main
  h: '#c9922e', // coin shade
};

function sprite(rows, over) {
  const p = over ? Object.assign({}, PAL, over) : PAL;
  const w = Math.max(...rows.map(r => r.length));
  const c = document.createElement('canvas');
  c.width = w; c.height = rows.length;
  const g = c.getContext('2d');
  rows.forEach((row, yy) => {
    for (let xx = 0; xx < row.length; xx++) {
      const k = row[xx];
      if (k === '.' || k === ' ') continue;
      const col = p[k];
      if (!col) continue;
      g.fillStyle = col;
      g.fillRect(xx, yy, 1, 1);
    }
  });
  return { c, w, h: rows.length };
}

/* ---------- Clawd (26 x 20), built to the mosaic reference:
   head slab with rectangular eyes over a wider body slab with side
   wings, four stubby legs. Faint seams keep the tiled-mosaic look. */
const CLAWD_MAIN = '#cc7f5e', CLAWD_SEAM = '#c07152', CLAWD_EYE = '#161311';
function buildClawd(legMode, dead) {
  const c = document.createElement('canvas');
  c.width = 26; c.height = 20;
  const g = c.getContext('2d');
  const rect = (x, y, w, h, col) => { g.fillStyle = col; g.fillRect(x, y, w, h); };
  const slab = (x, y, w, h) => {
    rect(x, y, w, h, CLAWD_MAIN);
    for (let sx = x + 3; sx < x + w - 1; sx += 4) rect(sx, y, 1, h, CLAWD_SEAM);
  };
  slab(4, 0, 18, 7);                    // head slab
  slab(2, 7, 22, 7);                    // body slab
  slab(0, 7, 2, 4); slab(24, 7, 2, 4);  // side wings
  if (dead) {
    for (const ex of [7, 15]) for (let i = 0; i < 4; i++) {
      rect(ex + i, 3 + i, 1, 1, CLAWD_EYE);
      rect(ex + 3 - i, 3 + i, 1, 1, CLAWD_EYE);
    }
  } else {
    rect(8, 3, 2, 4, CLAWD_EYE); rect(16, 3, 2, 4, CLAWD_EYE);
  }
  [6, 10, 16, 20].forEach((lx, i) => {
    const long = legMode === 'A' ? i % 2 === 0 : legMode === 'B' ? i % 2 === 1 : false;
    rect(lx, 14, 2, long ? 6 : 3, CLAWD_MAIN);
  });
  return { c, w: 26, h: 20 };
}
const sprRunA = buildClawd('A');
const sprRunB = buildClawd('B');
const sprJump = buildClawd('tuck');
const sprDead = buildClawd('tuck', true);

/* ---------- CAPTCHA obstacle: the reCAPTCHA logo, pixelated.
   Three-colour broken ring rendered pixel-by-pixel (no anti-aliasing)
   with the two arrowheads; the large variant gets the wordmark. */
const RC_LIGHT = '#5b9bf5', RC_DARK = '#2148a8', RC_GRAY = '#a9a9a9';
function buildLogo(S, label) {
  const lw = label ? label.length * 4 - 1 : 0;
  const w = Math.max(S + 4, lw + 2);
  const h = S + (label ? 8 : 0);
  const c = document.createElement('canvas');
  c.width = w; c.height = h;
  const g = c.getContext('2d');
  const ox = Math.floor((w - S) / 2);
  const cx = (S - 1) / 2, cy = (S - 1) / 2;
  const rO = S * 0.5, rI = S * 0.3;
  for (let y = 0; y < S; y++) for (let x = 0; x < S; x++) {
    const dx = x - cx, dy = y - cy, r = Math.hypot(dx, dy);
    if (r < rI || r > rO) continue;
    const a = Math.atan2(dy, dx) * 180 / Math.PI;
    let col = dy < 0 ? (dx < 0 ? RC_LIGHT : RC_DARK) : RC_GRAY;
    if ((a > 20 && a < 62) || a < -162) col = null; // breaks before the arrow tails
    if (col) { g.fillStyle = col; g.fillRect(ox + x, y, 1, 1); }
  }
  const rM = Math.round((rO + rI) / 2);
  const tri = (tx, ty, dir, col) => { // 4-row arrowhead, dir 1 = down, -1 = up
    g.fillStyle = col;
    for (let i = 0; i < 4; i++) g.fillRect(ox + tx - (3 - i), ty + i * dir, 2 * (3 - i) + 1, 1);
  };
  tri(Math.round(cx + rM), Math.round(cy) + 1, 1, RC_DARK);  // right head, clockwise down
  tri(Math.round(cx - rM), Math.round(cy), -1, RC_GRAY);     // left head, clockwise up
  if (label) {
    g.fillStyle = RC_GRAY;
    let tx = Math.floor((w - lw) / 2);
    for (const ch of label) {
      const gl = FONT[ch];
      if (gl) for (let r = 0; r < 5; r++) for (let cc = 0; cc < 3; cc++)
        if (gl[r][cc] === '1') g.fillRect(tx + cc, S + 2 + r, 1, 1);
      tx += 4;
    }
  }
  return { c, w, h };
}

/* ---------- token coins (diamond, tier-coloured) ---------- */
const COIN0 = [
  '....y....',
  '...yWy...',
  '..yyyyy..',
  '.yyKKKyy.',
  'yyyyKyyyh',
  '.yyyKyyh.',
  '..yyyyh..',
  '...yyh...',
  '....h....',
];
const COIN1 = [
  '....y....',
  '....y....',
  '...yWy...',
  '...yKy...',
  '...yKy...',
  '...yKy...',
  '...yyy...',
  '....h....',
  '....h....',
];
const TIER_COLS = [
  { y: '#f2c14e', h: '#c9922e' }, // gold
  { y: '#58c470', h: '#2e9e4f' }, // emerald
  { y: '#5ba0f2', h: '#2d6fc9' }, // sapphire
  { y: '#b07df2', h: '#7e4fc9' }, // amethyst
];
const coinSpr = TIER_COLS.map(tc => [sprite(COIN0, tc), sprite(COIN1, tc)]);

const sprCloud = sprite([
  '....CCCC........',
  '..CCCCCCCC.CCC..',
  '.CCCCCCCCCCCCCC.',
  'CCCCCCCCCCCCCCCC',
  '.CCCCCCCCCCCCC..',
]);

/* ---------- 3x5 pixel font ---------- */
const FONT = {
  A:['010','101','111','101','101'], B:['110','101','110','101','110'],
  C:['011','100','100','100','011'], D:['110','101','101','101','110'],
  E:['111','100','110','100','111'], F:['111','100','110','100','100'],
  G:['111','100','101','101','111'], H:['101','101','111','101','101'],
  I:['111','010','010','010','111'], J:['001','001','001','101','010'],
  K:['101','101','110','101','101'], L:['100','100','100','100','111'],
  M:['101','111','101','101','101'], N:['110','101','101','101','101'],
  O:['111','101','101','101','111'], P:['111','101','111','100','100'],
  Q:['111','101','101','111','001'], R:['111','101','110','101','101'],
  S:['011','100','010','001','110'], T:['111','010','010','010','010'],
  U:['101','101','101','101','111'], V:['101','101','101','101','010'],
  W:['101','101','101','111','101'], X:['101','101','010','101','101'],
  Y:['101','101','010','010','010'], Z:['111','001','010','100','111'],
  0:['111','101','101','101','111'], 1:['010','110','010','010','111'],
  2:['111','001','111','100','111'], 3:['111','001','111','001','111'],
  4:['101','101','111','001','001'], 5:['111','100','111','001','111'],
  6:['111','100','111','101','111'], 7:['111','001','001','010','010'],
  8:['111','101','111','101','111'], 9:['111','101','111','001','111'],
  '+':['000','010','111','010','000'], '-':['000','000','111','000','000'],
  '!':['010','010','010','000','010'], '.':['000','000','000','000','010'],
  ':':['000','010','000','010','000'], '/':['001','001','010','100','100'],
};
const sprLogoS = buildLogo(20);
const sprLogoL = buildLogo(28, 'RECAPTCHA');
function text(str, x, y, col, sc = 1) {
  ctx.fillStyle = col;
  str = String(str).toUpperCase();
  for (let i = 0; i < str.length; i++) {
    const g = FONT[str[i]];
    if (!g) continue;
    for (let r = 0; r < 5; r++)
      for (let c = 0; c < 3; c++)
        if (g[r][c] === '1') ctx.fillRect(x + i * 4 * sc + c * sc, y + r * sc, sc, sc);
  }
}
function textC(str, cx, y, col, sc = 1) {
  const w = String(str).length * 4 * sc - sc;
  text(str, Math.round(cx - w / 2), y, col, sc);
}

/* ---------- audio ---------- */
let AC = null;
function ensureAudio() {
  if (!AC) { try { AC = new (window.AudioContext || window.webkitAudioContext)(); } catch (e) {} }
  if (AC && AC.state === 'suspended') AC.resume();
}
function tone(f0, f1, dur, delay = 0, type = 'square', vol = 0.045) {
  if (!AC) return;
  const t = AC.currentTime + delay;
  const o = AC.createOscillator(), g = AC.createGain();
  o.type = type;
  o.frequency.setValueAtTime(f0, t);
  if (f1) o.frequency.exponentialRampToValueAtTime(f1, t + dur);
  g.gain.setValueAtTime(vol, t);
  g.gain.exponentialRampToValueAtTime(0.0001, t + dur);
  o.connect(g); g.connect(AC.destination);
  o.start(t); o.stop(t + dur + 0.03);
}
const sJump = () => tone(340, 720, 0.13);
const sCoin = () => { tone(988, 0, 0.055, 0, 'square', 0.035); tone(1319, 0, 0.07, 0.05, 'square', 0.035); };
const sTier = () => { tone(659, 0, 0.07); tone(784, 0, 0.07, 0.07); tone(988, 0, 0.1, 0.14); };
const sDie = () => tone(300, 70, 0.35, 0, 'sawtooth', 0.05);
const sMile = () => tone(880, 0, 0.06, 0, 'square', 0.025);

/* ---------- state ---------- */
const GRV = 0.26, JUMPV = -5.0;
const TIER_MS = 15000; // coin value goes up every 15s survived
let state = 'title', frameN = 0, paused = false, jumpHeld = false;
let hi = +(localStorage.getItem('clawdrun_hi') || 0);
let bank = +(localStorage.getItem('clawdrun_bank') || 0);
let speed, dist, score, elapsed, tokens, coinValue, newHi;
let obstacles, coins, popups, nextObs, nextCoin, milestone, scoreFlash;
let flashMsg = '', flashT = 0, deadAt = 0;
const clawd = { x: 24, y: GROUND - 20, vy: 0, onGround: true, run: 0 };

const clouds = [
  { x: 46, y: 24, v: 0.24 },
  { x: 158, y: 54, v: 0.18 },
  { x: 262, y: 14, v: 0.3 },
];
const PEB_COLS = ['#cfc3ab', '#b9ab90'];
const pebbles = [];
for (let i = 0; i < 40; i++) pebbles.push({
  x: Math.random() * W,
  y: GROUND + 5 + Math.random() * 22,
  s: Math.random() < 0.3 ? 2 : 1,
  c: PEB_COLS[i % 2],
});

function reset() {
  speed = 2; dist = 0; score = 0; elapsed = 0; tokens = 0; coinValue = 1;
  obstacles = []; coins = []; popups = [];
  nextObs = 420; nextCoin = 240; milestone = 0; scoreFlash = 0;
  flashMsg = ''; flashT = 0; newHi = false; paused = false;
  clawd.y = GROUND - 20; clawd.vy = 0; clawd.onGround = true;
}
reset();

function start() { reset(); state = 'run'; }

function jump() {
  if (clawd.onGround) {
    clawd.vy = JUMPV; clawd.onGround = false; jumpHeld = true; sJump();
  }
}
function action() {
  ensureAudio();
  if (state === 'title') start();
  else if (state === 'over') { if (performance.now() - deadAt > 450) start(); }
  else if (state === 'run') { if (paused) { paused = false; return; } jump(); }
}
function die() {
  state = 'over'; deadAt = performance.now();
  if (score > hi) { hi = score; newHi = true; localStorage.setItem('clawdrun_hi', hi); }
  bank += tokens; localStorage.setItem('clawdrun_bank', bank);
  sDie();
}

const overlap = (a, b) => a.x < b.x + b.w && a.x + a.w > b.x && a.y < b.y + b.h && a.y + a.h > b.y;

function spawnObstacle() {
  const r = Math.random();
  let parts;
  if (r < 0.45) parts = [{ s: sprLogoS, dx: 0 }];
  else if (r < 0.65 && speed > 2.6) parts = [{ s: sprLogoS, dx: 0 }, { s: sprLogoS, dx: 26 }];
  else parts = [{ s: sprLogoL, dx: 0 }];
  const w = Math.max(...parts.map(p => p.dx + p.s.w));
  const h = Math.max(...parts.map(p => p.s.h));
  let x = W + 8;
  for (const c of coins) if (c.x + 25 > x) x = c.x + 25; // clear coins near the spawn edge
  obstacles.push({ x, y: GROUND - h, w, h, parts });
}

function spawnCoins() {
  let bx = W + 12;
  for (const o of obstacles) if (o.x + o.w + 34 > bx) bx = o.x + o.w + 34;
  const r = Math.random();
  if (r < 0.4) {           // ground row
    const n = 3 + ((Math.random() * 2) | 0);
    for (let i = 0; i < n; i++) coins.push({ x: bx + i * 13, y: GROUND - 16 });
  } else if (r < 0.75) {   // jump arc
    [-14, -32, -42, -32, -14].forEach((dy, i) =>
      coins.push({ x: bx + i * 12, y: GROUND - 16 + dy }));
  } else {                 // high row (needs a full jump)
    for (let i = 0; i < 3; i++) coins.push({ x: bx + i * 13, y: GROUND - 52 });
  }
}

function update(dt) {
  speed = Math.min(speed + 0.0011 * dt, 5.2);
  dist += speed * dt;
  elapsed += dt * (1000 / 60);
  score = Math.floor(dist / 6);

  if (Math.floor(score / 100) > milestone) {
    milestone = Math.floor(score / 100); scoreFlash = 24; sMile();
  }
  if (scoreFlash > 0) scoreFlash -= dt;

  const cvNew = 1 + Math.floor(elapsed / TIER_MS);
  if (cvNew !== coinValue) {
    coinValue = cvNew;
    flashMsg = 'TOKEN VALUE UP! +' + coinValue + ' PER COIN';
    flashT = 130; sTier();
  }
  if (flashT > 0) flashT -= dt;

  if (!clawd.onGround) {
    clawd.vy += GRV * (clawd.vy < 0 && jumpHeld ? 0.52 : 1) * dt;
    clawd.y += clawd.vy * dt;
    if (clawd.y >= GROUND - 20) { clawd.y = GROUND - 20; clawd.vy = 0; clawd.onGround = true; }
  } else clawd.run += dt;

  nextObs -= speed * dt;
  if (nextObs <= 0) { spawnObstacle(); nextObs = 42 * speed + 68 + Math.random() * 150; }
  nextCoin -= speed * dt;
  if (nextCoin <= 0) { spawnCoins(); nextCoin = 300 + Math.random() * 450; }

  obstacles.forEach(o => o.x -= speed * dt);
  coins.forEach(c => c.x -= speed * dt);
  obstacles = obstacles.filter(o => o.x + o.w > -4);
  coins = coins.filter(c => c.x > -12 && !c.got);

  clouds.forEach(cl => {
    cl.x -= cl.v * dt * (speed / 2);
    if (cl.x < -18) { cl.x = W + 10; cl.y = 12 + Math.random() * 56; }
  });
  pebbles.forEach(p => {
    p.x -= speed * dt;
    if (p.x < -3) { p.x = W + Math.random() * 30; p.y = GROUND + 5 + Math.random() * 22; }
  });
  popups.forEach(p => { p.y -= 0.5 * dt; p.t += dt; });
  popups = popups.filter(p => p.t < 50);

  const hb = { x: clawd.x + 4, y: clawd.y + 2, w: 18, h: 17 };
  for (const c of coins) {
    if (!c.got && overlap(hb, { x: c.x, y: c.y, w: 9, h: 9 })) {
      c.got = true; tokens += coinValue;
      popups.push({ x: c.x - 2, y: c.y - 6, txt: '+' + coinValue, t: 0 });
      sCoin();
    }
  }
  for (const o of obstacles) {
    if (overlap(hb, { x: o.x + 1, y: o.y + 1, w: o.w - 2, h: o.h - 2 })) { die(); break; }
  }
}

/* ---------- render ---------- */
function drawScene() {
  ctx.fillStyle = BG;
  ctx.fillRect(0, 0, W, H);
  clouds.forEach(cl => ctx.drawImage(sprCloud.c, Math.round(cl.x), cl.y));
  ctx.fillStyle = '#8f8272';
  ctx.fillRect(0, GROUND, W, 2);
  pebbles.forEach(p => { ctx.fillStyle = p.c; ctx.fillRect(Math.round(p.x), Math.round(p.y), p.s, 1); });
}
function drawClawd() {
  let s;
  if (state === 'over') s = sprDead;
  else if (!clawd.onGround) s = sprJump;
  else s = Math.floor(clawd.run / 6) % 2 ? sprRunB : sprRunA;
  ctx.drawImage(s.c, Math.round(clawd.x), Math.round(clawd.y));
}
function drawHUD() {
  const tier = Math.min(coinValue - 1, 3);
  ctx.drawImage(coinSpr[tier][0].c, 5, 4);
  text('TOKENS ' + tokens, 18, 6, '#2a2019');
  const right = 'HI ' + String(hi).padStart(5, '0') + '  ' + String(score).padStart(5, '0');
  const flash = scoreFlash > 0 && Math.floor(frameN / 4) % 2 === 0;
  text(right, W - right.length * 4 - 5, 6, flash ? '#d97757' : '#6f675c');
  if (flashT > 0 && Math.floor(frameN / 8) % 2 === 0) textC(flashMsg, W / 2, 54, '#d97757');
}

function render() {
  drawScene();

  const tier = Math.min(coinValue - 1, 3);
  const cf = Math.floor(frameN / 14) % 2;
  if (state !== 'title') {
    coins.forEach(c => ctx.drawImage(coinSpr[tier][cf].c, Math.round(c.x), Math.round(c.y)));
    obstacles.forEach(o => o.parts.forEach(p =>
      ctx.drawImage(p.s.c, Math.round(o.x + p.dx), GROUND - p.s.h)));
  }
  drawClawd();
  popups.forEach(p => text(p.txt, Math.round(p.x), Math.round(p.y), '#a0622d'));

  if (state === 'title') {
    ctx.drawImage(sprLogoL.c, 246, GROUND - sprLogoL.h);
    for (let i = 0; i < 3; i++) ctx.drawImage(coinSpr[i][cf].c, 190 + i * 14, GROUND - 16);
    textC('CLAWD RUN', W / 2, 34, '#2a2019', 3);
    textC('DODGE THE CAPTCHAS. GRAB TOKEN COINS.', W / 2, 66, '#6f675c');
    textC('COINS ARE WORTH MORE THE LONGER YOU SURVIVE', W / 2, 76, '#6f675c');
    if (Math.floor(frameN / 24) % 2 === 0) textC('PRESS SPACE OR TAP TO RUN', W / 2, 96, '#d97757');
    textC('TOKEN BANK ' + bank + '   HI ' + String(hi).padStart(5, '0'), W / 2, 116, '#a0988a');
  } else {
    drawHUD();
  }

  if (state === 'over') {
    ctx.fillStyle = '#2a2019';
    ctx.fillRect(56, 36, 208, 82);
    ctx.fillStyle = '#faf6ec';
    ctx.fillRect(58, 38, 204, 78);
    textC('CAPTCHA FAILED!', W / 2, 44, '#b3402e', 2);
    if (newHi && Math.floor(frameN / 10) % 2 === 0) textC('NEW HI!', W / 2, 58, '#d97757');
    textC('SCORE ' + String(score).padStart(5, '0') + '  HI ' + String(hi).padStart(5, '0'), W / 2, 70, '#2a2019');
    textC('TOKENS EARNED +' + tokens, W / 2, 82, '#a0622d');
    textC('TOKEN BANK ' + bank, W / 2, 92, '#a0622d');
    if (Math.floor(frameN / 24) % 2 === 0) textC('PRESS SPACE TO RUN AGAIN', W / 2, 106, '#6f675c');
  }
  if (paused && state === 'run') {
    ctx.fillStyle = 'rgba(240,233,220,0.7)';
    ctx.fillRect(0, 0, W, H);
    textC('PAUSED', W / 2, 80, '#2a2019', 2);
    textC('PRESS P OR SPACE', W / 2, 98, '#6f675c');
  }
}

/* ---------- main loop ---------- */
let last = performance.now();
function loop(t) {
  requestAnimationFrame(loop);
  let dt = (t - last) / 16.667;
  last = t;
  if (dt <= 0) return;
  if (dt > 3) dt = 3;
  frameN++;
  if (state === 'run' && !paused) update(dt);
  else if (state === 'title') {
    clawd.run += dt;
    clouds.forEach(cl => { cl.x -= cl.v * dt; if (cl.x < -18) cl.x = W + 10; });
  }
  render();
}
requestAnimationFrame(loop);

/* ---------- input ---------- */
addEventListener('keydown', e => {
  if (e.code === 'Space' || e.code === 'ArrowUp' || e.code === 'KeyW') {
    e.preventDefault();
    action();
  } else if (e.code === 'KeyP' && state === 'run') {
    paused = !paused;
  } else if (e.code === 'ArrowDown' && state === 'run' && !clawd.onGround) {
    e.preventDefault();
    clawd.vy = Math.max(clawd.vy, 4); // fast fall
    jumpHeld = false;
  }
});
addEventListener('keyup', e => {
  if (e.code === 'Space' || e.code === 'ArrowUp' || e.code === 'KeyW') jumpHeld = false;
});
cv.addEventListener('pointerdown', e => { e.preventDefault(); action(); });
addEventListener('pointerup', () => { jumpHeld = false; });
addEventListener('blur', () => { if (state === 'run') paused = true; });

/* ---------- favicon from the Clawd sprite ---------- */
(() => {
  const f = document.createElement('canvas');
  f.width = f.height = 32;
  const g = f.getContext('2d');
  g.imageSmoothingEnabled = false;
  g.drawImage(sprRunA.c, 3, 6);
  const l = document.createElement('link');
  l.rel = 'icon';
  l.href = f.toDataURL();
  document.head.appendChild(l);
})();

/* hooks for automated testing */
window.__clawd = {
  start: () => { if (state !== 'run') start(); },
  step: (n = 1) => {
    for (let i = 0; i < n && state === 'run'; i++) { frameN++; update(1); }
    render();
  },
  die: () => { if (state === 'run') die(); },
  jump,
  addTime: ms => { elapsed += ms; },
  get state() { return state; },
  get obstacles() { return obstacles.map(o => ({ x: o.x, w: o.w, h: o.h })); },
  get coinsPos() { return coins.map(c => ({ x: c.x, y: c.y })); },
  get onGround() { return clawd.onGround; },
  get score() { return score; },
  get tokens() { return tokens; },
  get coinValue() { return coinValue; },
};
})();
