/* Propulsion Research Series — shared Plotly theming helper.
   Load Plotly first:
     <script src="https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js"></script>
     <script src="../assets/plot.js"></script>

   Usage:
     PR.plot('divId', traces, { xaxis: { title: 'x/L' }, yaxis: { title: 'Mach' } });
     PR.series(i)        -> i-th categorical colour from paper.css (--s1..--s6)
     PR.loadJSON(url)    -> fetch + parse, returns a Promise
   Every plot re-themes itself automatically when the OS light/dark mode flips. */
(function () {
  const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
  const registry = new Map();

  function series(i) { return css('--s' + ((i % 6) + 1)); }

  function baseLayout() {
    const ink = css('--ink'), ink2 = css('--ink-2'), line = css('--line');
    const axis = {
      gridcolor: line, zerolinecolor: line, linecolor: css('--line-2'),
      tickfont: { family: 'JetBrains Mono, monospace', size: 11, color: ink2 },
      title: { font: { family: 'DM Sans, sans-serif', size: 13, color: ink2 } },
      ticks: 'outside', ticklen: 4, tickcolor: css('--line-2'),
      automargin: true,
    };
    return {
      paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
      font: { family: 'DM Sans, sans-serif', color: ink, size: 12 },
      margin: { l: 64, r: 22, t: 28, b: 52 },
      xaxis: { ...axis }, yaxis: { ...axis },
      legend: { orientation: 'h', x: 0, y: 1.02, yanchor: 'bottom', bgcolor: 'rgba(0,0,0,0)',
                font: { size: 12, color: ink2 } },
      hoverlabel: { font: { family: 'JetBrains Mono, monospace', size: 12 } },
      colorway: [1, 2, 3, 4, 5, 6].map((k) => css('--s' + k)),
    };
  }

  function deepMerge(a, b) {
    const out = { ...a };
    for (const k in b) {
      const v = b[k];
      if (v && typeof v === 'object' && !Array.isArray(v) && a[k] && typeof a[k] === 'object') out[k] = deepMerge(a[k], v);
      else out[k] = v;
    }
    return out;
  }

  // Axis titles may be given as plain strings.
  function normalise(layout) {
    const l = { ...layout };
    for (const k of Object.keys(l)) {
      if (/^[xy]axis\d*$/.test(k) && l[k] && typeof l[k].title === 'string') l[k] = { ...l[k], title: { text: l[k].title } };
    }
    return l;
  }

  function plot(id, traces, layout = {}, config = {}) {
    const el = typeof id === 'string' ? document.getElementById(id) : id;
    const lay = deepMerge(baseLayout(), normalise(layout));
    for (const k of Object.keys(lay)) {
      if (/^[xy]axis\d+$/.test(k)) lay[k] = deepMerge(baseLayout().xaxis, lay[k]);
    }
    const cfg = { responsive: true, displaylogo: false, modeBarButtonsToRemove: ['lasso2d', 'select2d'], ...config };
    registry.set(el, { layout, config });
    return Plotly.react(el, traces, lay, cfg);
  }

  function retheme() {
    for (const [el, { layout, config }] of registry) {
      if (!el.isConnected) { registry.delete(el); continue; }
      plot(el, el.data, layout, config);
    }
  }
  window.matchMedia('(prefers-color-scheme: light)').addEventListener('change', retheme);

  async function loadJSON(url) {
    const r = await fetch(url);
    if (!r.ok) throw new Error(url + ': ' + r.status);
    return r.json();
  }

  window.PR = { plot, series, css, loadJSON, retheme };
})();
