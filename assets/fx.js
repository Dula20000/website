/* ════════════════════════════════════════════════════════════
   DULRAN — shared FX engine
   Provides on every page:
     • custom lag cursor (dot + ring, grows on links/[data-hover])
     • pointer-parallax depth layers      [data-depth="0.05"]
     • 3D tilt cards with glare           [data-tilt]
     • magnetic pull elements             [data-magnetic]
     • scroll reveals                     .reveal → .visible
     • scroll progress bar                .dk-progress
     • DK.stars(canvas, opts) — starfield with pointer parallax
       and scroll-velocity warp (aerospace signature)
═════════════════════════════════════════════════════════════ */
(function () {
  'use strict';

  const DK = (window.DK = window.DK = {});
  const fine = window.matchMedia('(hover: hover) and (pointer: fine)').matches;
  const reduced = window.matchMedia('(prefers-reduced-motion: reduce)').matches;

  /* shared pointer state, normalized -1..1 from screen center */
  const pointer = (DK.pointer = { x: 0, y: 0, px: 0, py: 0 });
  window.addEventListener('pointermove', (e) => {
    pointer.px = e.clientX;
    pointer.py = e.clientY;
    pointer.x = (e.clientX / window.innerWidth) * 2 - 1;
    pointer.y = (e.clientY / window.innerHeight) * 2 - 1;
  }, { passive: true });

  /* ── custom cursor ── */
  function initCursor() {
    if (!fine || reduced) return;
    const dot = document.createElement('div');
    const ring = document.createElement('div');
    dot.className = 'dk-cursor-dot';
    ring.className = 'dk-cursor-ring';
    document.body.append(dot, ring);

    let rx = innerWidth / 2, ry = innerHeight / 2;
    (function loop() {
      rx += (pointer.px - rx) * 0.16;
      ry += (pointer.py - ry) * 0.16;
      dot.style.transform = `translate(${pointer.px - 3}px, ${pointer.py - 3}px)`;
      ring.style.transform = `translate(${rx - 17}px, ${ry - 17}px)`;
      requestAnimationFrame(loop);
    })();

    document.addEventListener('pointerover', (e) => {
      if (e.target.closest('a, button, [data-hover]'))
        document.body.classList.add('dk-cursor-hover');
    });
    document.addEventListener('pointerout', (e) => {
      if (e.target.closest('a, button, [data-hover]'))
        document.body.classList.remove('dk-cursor-hover');
    });
  }

  /* ── pointer-parallax depth layers ── */
  function initDepth() {
    const els = document.querySelectorAll('[data-depth]');
    if (!els.length || !fine || reduced) return;
    (function loop() {
      els.forEach((el) => {
        const d = parseFloat(el.dataset.depth) || 0.05;
        const x = -pointer.x * d * 100;
        const y = -pointer.y * d * 100;
        el.style.transform = `translate3d(${x}px, ${y}px, 0)`;
      });
      requestAnimationFrame(loop);
    })();
  }

  /* ── 3D tilt with glare ── */
  function initTilt() {
    if (!fine || reduced) return;
    document.querySelectorAll('[data-tilt]').forEach((el) => {
      const max = parseFloat(el.dataset.tilt) || 10;
      if (!el.querySelector('.tilt-glare')) {
        const glare = document.createElement('div');
        glare.className = 'tilt-glare';
        el.appendChild(glare);
      }
      el.addEventListener('pointermove', (e) => {
        const r = el.getBoundingClientRect();
        const nx = (e.clientX - r.left) / r.width - 0.5;
        const ny = (e.clientY - r.top) / r.height - 0.5;
        el.style.transform =
          `perspective(800px) rotateY(${nx * max}deg) rotateX(${-ny * max}deg) translateZ(6px)`;
        el.style.setProperty('--gx', `${(nx + 0.5) * 100}%`);
        el.style.setProperty('--gy', `${(ny + 0.5) * 100}%`);
      });
      el.addEventListener('pointerleave', () => {
        el.style.transition = 'transform 0.6s cubic-bezier(0.16,1,0.3,1)';
        el.style.transform = 'perspective(800px) rotateX(0) rotateY(0)';
        setTimeout(() => (el.style.transition = ''), 600);
      });
    });
  }

  /* ── magnetic elements ── */
  function initMagnetic() {
    if (!fine || reduced) return;
    document.querySelectorAll('[data-magnetic]').forEach((el) => {
      const strength = parseFloat(el.dataset.magnetic) || 0.3;
      el.addEventListener('pointermove', (e) => {
        const r = el.getBoundingClientRect();
        const x = (e.clientX - r.left - r.width / 2) * strength;
        const y = (e.clientY - r.top - r.height / 2) * strength;
        el.style.transform = `translate(${x}px, ${y}px)`;
      });
      el.addEventListener('pointerleave', () => {
        el.style.transition = 'transform 0.5s cubic-bezier(0.16,1,0.3,1)';
        el.style.transform = 'translate(0,0)';
        setTimeout(() => (el.style.transition = ''), 500);
      });
    });
  }

  /* ── scroll reveals ── */
  function initReveal() {
    const els = document.querySelectorAll('.reveal');
    if (!els.length) return;
    const io = new IntersectionObserver((entries) => {
      entries.forEach((en) => {
        if (en.isIntersecting) {
          en.target.classList.add('visible');
          io.unobserve(en.target);
        }
      });
    }, { threshold: 0.15 });
    els.forEach((el) => io.observe(el));
  }

  /* ── scroll progress bar ── */
  function initProgress() {
    const bar = document.querySelector('.dk-progress');
    if (!bar) return;
    const update = () => {
      const max = document.documentElement.scrollHeight - innerHeight;
      bar.style.transform = `scaleX(${max > 0 ? scrollY / max : 0})`;
    };
    window.addEventListener('scroll', update, { passive: true });
    update();
  }

  /* ── starfield: pointer parallax + scroll-velocity warp ── */
  DK.stars = function (canvas, opts = {}) {
    const ctx = canvas.getContext('2d');
    const COUNT = opts.count || 220;
    const COLOR = opts.color || '255, 255, 255';
    const ACCENTS = opts.accents || ['255, 122, 41', '232, 176, 75'];
    let w, h, stars = [];
    let lastScroll = scrollY, warp = 0;

    function resize() {
      w = canvas.width = canvas.offsetWidth * devicePixelRatio;
      h = canvas.height = canvas.offsetHeight * devicePixelRatio;
    }
    window.addEventListener('resize', resize);
    resize();

    for (let i = 0; i < COUNT; i++) {
      stars.push({
        x: Math.random() * w,
        y: Math.random() * h,
        z: Math.random() * 0.9 + 0.1,           // depth 0.1..1
        r: Math.random() * 1.4 + 0.3,
        tw: Math.random() * Math.PI * 2,
        c: Math.random() < 0.12
          ? ACCENTS[Math.floor(Math.random() * ACCENTS.length)]
          : COLOR,
      });
    }

    window.addEventListener('scroll', () => {
      warp += Math.abs(scrollY - lastScroll) * 0.06;
      lastScroll = scrollY;
    }, { passive: true });

    let t = 0;
    (function frame() {
      t += 0.016;
      warp *= 0.92;
      ctx.clearRect(0, 0, w, h);
      for (const s of stars) {
        const px = pointer.x * s.z * 28 * devicePixelRatio;
        const py = pointer.y * s.z * 28 * devicePixelRatio;
        s.y += (0.08 + warp * 2.4) * s.z * devicePixelRatio;
        if (s.y > h + 4) { s.y = -4; s.x = Math.random() * w; }
        const a = (0.25 + 0.55 * s.z) * (0.7 + 0.3 * Math.sin(t * 2 + s.tw));
        ctx.beginPath();
        const stretch = 1 + warp * 6 * s.z;
        ctx.ellipse(s.x - px, s.y - py, s.r * devicePixelRatio,
                    s.r * stretch * devicePixelRatio, 0, 0, Math.PI * 2);
        ctx.fillStyle = `rgba(${s.c}, ${a})`;
        ctx.fill();
      }
      requestAnimationFrame(frame);
    })();
  };

  /* ── floating ember/gold dust particles ── */
  DK.dust = function (canvas, opts = {}) {
    const ctx = canvas.getContext('2d');
    const COUNT = opts.count || 60;
    const COLORS = opts.colors || ['232, 176, 75', '255, 122, 41', '255, 255, 255'];
    let w, h, ps = [];
    function resize() {
      w = canvas.width = canvas.offsetWidth * devicePixelRatio;
      h = canvas.height = canvas.offsetHeight * devicePixelRatio;
    }
    window.addEventListener('resize', resize);
    resize();
    for (let i = 0; i < COUNT; i++) {
      ps.push({
        x: Math.random() * w, y: Math.random() * h,
        vx: (Math.random() - 0.5) * 0.18, vy: -(Math.random() * 0.3 + 0.08),
        r: Math.random() * 1.8 + 0.5,
        a: Math.random() * 0.5 + 0.15,
        tw: Math.random() * Math.PI * 2,
        c: COLORS[Math.floor(Math.random() * COLORS.length)],
      });
    }
    let t = 0;
    (function frame() {
      t += 0.016;
      ctx.clearRect(0, 0, w, h);
      for (const p of ps) {
        p.x += (p.vx + pointer.x * 0.25) * devicePixelRatio;
        p.y += p.vy * devicePixelRatio;
        if (p.y < -6) { p.y = h + 6; p.x = Math.random() * w; }
        if (p.x < -6) p.x = w + 6;
        if (p.x > w + 6) p.x = -6;
        const a = p.a * (0.6 + 0.4 * Math.sin(t * 1.6 + p.tw));
        ctx.beginPath();
        ctx.arc(p.x, p.y, p.r * devicePixelRatio, 0, Math.PI * 2);
        ctx.fillStyle = `rgba(${p.c}, ${a})`;
        ctx.fill();
      }
      requestAnimationFrame(frame);
    })();
  };

  /* boot */
  function boot() {
    initCursor();
    initDepth();
    initTilt();
    initMagnetic();
    initReveal();
    initProgress();
  }
  if (document.readyState === 'loading')
    document.addEventListener('DOMContentLoaded', boot);
  else boot();
})();
