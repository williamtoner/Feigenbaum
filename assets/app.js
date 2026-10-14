/* Feigenbaum constants: interactive page. Data: window.FEIG (assets/data.js). */
(() => {
  "use strict";
  const F = window.FEIG;
  const DEG = [2, 4, 6, 8];
  const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
  const seriesColor = (z) => css(`--s${z}`);
  const num = (key) => parseFloat(F.digits[key].slice(0, 18));
  const DELTA = Object.fromEntries(DEG.map((z) => [z, num(`delta${z}`)]));
  const ALPHA = Object.fromEntries(DEG.map((z) => [z, num(`alpha${z}`)]));
  const MUINF = Object.fromEntries(DEG.map((z) => [z, parseFloat(F.superstable[z].mu_inf)]));
  const SUP = { 2: "²", 4: "⁴", 6: "⁶", 8: "⁸" };
  const fmt = (n) => n.toLocaleString("en-US");

  // ------------------------------------------------------------ Feigenbaum function g (low precision, for plots)
  function gFun(z) {
    const c = F.cheb[z];
    return (x) => {
      const y = 2 * Math.pow(Math.abs(x), z) - 1;
      let b1 = 0, b2 = 0;
      for (let k = c.length - 1; k >= 1; k--) { const t = c[k] + 2 * y * b1 - b2; b2 = b1; b1 = t; }
      return c[0] / 2 + y * b1 - b2;
    };
  }

  // ------------------------------------------------------------ plotly theming
  function baseLayout(extra) {
    const text = css("--text-2"), grid = css("--grid");
    return Object.assign({
      paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)",
      font: { family: "Inter, system-ui, sans-serif", size: 13, color: text },
      margin: { l: 58, r: 24, t: 16, b: 48 },
      hoverlabel: { bgcolor: css("--surface-2"), bordercolor: css("--border"), font: { color: css("--text"), family: "Inter" } },
      legend: { orientation: "h", y: 1.08, x: 0, font: { size: 12 } },
      xaxis: { gridcolor: grid, zerolinecolor: grid, linecolor: css("--border"), ticks: "outside", tickcolor: css("--border") },
      yaxis: { gridcolor: grid, zerolinecolor: grid, linecolor: css("--border"), ticks: "outside", tickcolor: css("--border") },
    }, extra);
  }
  const PCONF = { displayModeBar: false, responsive: true };

  // ------------------------------------------------------------ hero: digit streams + faint diagram
  function heroStreams() {
    const items = [["stream-delta", F.digits.delta2], ["stream-alpha", F.digits.alpha2]];
    const els = items.map(([id, s]) => ({ el: document.getElementById(id), s, pos: 0 }));
    const reduce = window.matchMedia("(prefers-reduced-motion: reduce)").matches;
    const show = () => els.forEach((o) => {
      const start = Math.floor(o.pos);
      o.el.textContent = o.s.slice(start, start + 140);
    });
    if (reduce) { show(); return; }
    let last = performance.now();
    const t0 = last + 3000;
    function tick(t) {
      const dt = (t - last) / 1000; last = t;
      if (t < t0) { requestAnimationFrame(tick); return; }
      els.forEach((o, i) => { o.pos = (o.pos + dt * (i ? 5.5 : 7)) % (o.s.length - 160); });
      show();
      requestAnimationFrame(tick);
    }
    show();
    requestAnimationFrame(tick);
  }

  function densityRender(ctx, W, H, win, z, flip, quality, inkVar) {
    const cols = Math.max(80, Math.floor(W * quality));
    const rows = Math.max(60, Math.floor(H * quality));
    const dens = new Uint32Array(cols * rows);
    const span = win.m1 - win.m0, xs = win.x1 - win.x0;
    const depth = Math.max(0, Math.log(2 / span) / Math.log(4.0));
    const transient = Math.min(30000, Math.floor(300 + 160 * Math.pow(2, depth)));
    const samples = Math.min(20000, Math.floor((240 + 90 * Math.pow(2, depth)) * Math.max(quality, 0.5)));
    let max = 1;
    for (let i = 0; i < cols; i++) {
      const mu = win.m0 + ((i + 0.5) / cols) * span;
      let x = 0;
      for (let t = 0; t < transient; t++) x = 1 - mu * Math.pow(Math.abs(x), z);
      for (let t = 0; t < samples; t++) {
        x = 1 - mu * Math.pow(Math.abs(x), z);
        let r = (win.x1 - x) / xs;
        if (flip) r = 1 - r;
        if (r >= 0 && r < 1) {
          const k = Math.floor(r * rows) * cols + i;
          const v = ++dens[k];
          if (v > max) max = v;
        }
      }
    }
    const img = ctx.createImageData(cols, rows);
    const ink = css(inkVar).split(",").map((s) => parseInt(s, 10));
    const lmax = Math.log(1 + max * 0.25);
    for (let k = 0; k < dens.length; k++) {
      if (!dens[k]) continue;
      const a = Math.min(1, Math.log(1 + dens[k]) / lmax);
      img.data[4 * k] = ink[0]; img.data[4 * k + 1] = ink[1]; img.data[4 * k + 2] = ink[2];
      img.data[4 * k + 3] = Math.floor(40 + 215 * a);
    }
    const off = document.createElement("canvas");
    off.width = cols; off.height = rows;
    off.getContext("2d").putImageData(img, 0, 0);
    ctx.imageSmoothingEnabled = true;
    ctx.drawImage(off, 0, 0, W, H);
  }

  function heroCanvas() {
    const cv = document.getElementById("hero-canvas");
    const draw = () => {
      const r = cv.getBoundingClientRect(), dpr = Math.min(2, window.devicePixelRatio || 1);
      cv.width = r.width * dpr; cv.height = r.height * dpr;
      const ctx = cv.getContext("2d");
      ctx.clearRect(0, 0, cv.width, cv.height);
      densityRender(ctx, cv.width, cv.height, { m0: 0.55, m1: 2.0, x0: -1.05, x1: 1.05 }, 2, false, 0.6, "--hero-ink");
      // fade toward the text side
      const grd = ctx.createLinearGradient(0, 0, cv.width, 0);
      grd.addColorStop(0, css("--bg")); grd.addColorStop(0.55, "rgba(0,0,0,0)");
      ctx.globalCompositeOperation = "destination-out";
      ctx.fillStyle = grd; ctx.fillRect(0, 0, cv.width, cv.height);
      ctx.globalCompositeOperation = "source-over";
    };
    draw();
    return draw;
  }

  // ------------------------------------------------------------ interactive bifurcation diagram
  function bifurcation() {
    const cv = document.getElementById("bif"), tip = document.getElementById("bif-tip");
    const levelEl = document.getElementById("bif-level");
    let z = 2, level = 0, win = null, flip = false, anim = null, drag = null;

    const levelWindow = (zz, k) => {
      const mi = MUINF[zz], d = Math.pow(DELTA[zz], -k), a = Math.pow(ALPHA[zz], -k);
      const L = mi - 0.2, R = 2.0 - mi;
      return { m0: mi - L * d, m1: mi + R * d, x0: -1.05 * a, x1: 1.05 * a };
    };
    const size = () => {
      const r = cv.getBoundingClientRect(), dpr = Math.min(2, window.devicePixelRatio || 1);
      cv.width = Math.floor(r.width * dpr); cv.height = Math.floor(r.height * dpr);
    };
    function draw(quality = 1) {
      const ctx = cv.getContext("2d");
      ctx.fillStyle = css("--surface"); ctx.fillRect(0, 0, cv.width, cv.height);
      densityRender(ctx, cv.width, cv.height, win, z, flip, quality, "--bif-ink");
      // superstable parameters as ticks along the bottom (only those inside the window)
      const mus = F.superstable[z].mu, dpr = cv.width / cv.getBoundingClientRect().width;
      ctx.font = `${11 * dpr}px Inter, sans-serif`;
      ctx.fillStyle = css("--muted"); ctx.strokeStyle = css("--muted");
      let shown = 0;
      for (let n = 1; n < mus.length; n++) {
        const px = ((mus[n] - win.m0) / (win.m1 - win.m0)) * cv.width;
        if (px < 0 || px > cv.width) continue;
        if (shown > 0 && n > 1) {
          const prev = ((mus[n - 1] - win.m0) / (win.m1 - win.m0)) * cv.width;
          if (px - prev < 14 * dpr) continue;
        }
        ctx.globalAlpha = 0.5;
        ctx.beginPath(); ctx.moveTo(px, cv.height); ctx.lineTo(px, cv.height - 10 * dpr); ctx.stroke();
        ctx.globalAlpha = 0.9;
        ctx.fillText(`2${n > 1 ? superscript(n) : ""}`, px + 3 * dpr, cv.height - 4 * dpr);
        shown++;
      }
      ctx.globalAlpha = 1;
      const mi = MUINF[z];
      const pxi = ((mi - win.m0) / (win.m1 - win.m0)) * cv.width;
      if (pxi > 0 && pxi < cv.width) {
        ctx.setLineDash([4 * dpr, 4 * dpr]); ctx.strokeStyle = css("--accent"); ctx.globalAlpha = 0.6;
        ctx.beginPath(); ctx.moveTo(pxi, 0); ctx.lineTo(pxi, cv.height); ctx.stroke();
        ctx.setLineDash([]); ctx.globalAlpha = 1;
        ctx.fillStyle = css("--accent");
        ctx.fillText("μ∞", pxi + 4 * dpr, 14 * dpr);
      }
      levelEl.textContent = level === null ? "custom window" :
        `level ${level}: window ÷ ${fmt(+Math.pow(DELTA[z], level).toPrecision(4))} horizontally, ÷ ${fmt(+Math.pow(ALPHA[z], level).toPrecision(4))} vertically${flip ? ", flipped" : ""}`;
    }
    function superscript(n) {
      const map = "⁰¹²³⁴⁵⁶⁷⁸⁹";
      return String(n).split("").map((d) => map[+d]).join("");
    }
    function setLevel(k, animate) {
      const target = levelWindow(z, k);
      const targetFlip = k % 2 === 1;
      if (!animate || !win) { win = target; flip = targetFlip; level = k; draw(1); return; }
      const from = Object.assign({}, win), t0 = performance.now(), dur = 1100;
      cancelAnimationFrame(anim);
      const mi = MUINF[z];
      const step = (t) => {
        const u = Math.min(1, (t - t0) / dur), e = u < 0.5 ? 2 * u * u : 1 - Math.pow(-2 * u + 2, 2) / 2;
        // interpolate log of distances from (mu_inf, 0)
        const lerpLog = (a, b) => Math.sign(b) * Math.exp(Math.log(Math.abs(a)) * (1 - e) + Math.log(Math.abs(b)) * e);
        win = {
          m0: mi + lerpLog(from.m0 - mi, target.m0 - mi), m1: mi + lerpLog(from.m1 - mi, target.m1 - mi),
          x0: lerpLog(from.x0, target.x0), x1: lerpLog(from.x1, target.x1),
        };
        if (u >= 0.5) flip = targetFlip;
        level = k;
        draw(u < 1 ? 0.35 : 1);
        if (u < 1) anim = requestAnimationFrame(step);
      };
      anim = requestAnimationFrame(step);
    }
    document.getElementById("bif-deeper").addEventListener("click", () => {
      if (level === null) { setLevel(0, false); return; }
      setLevel(Math.min(level + 1, 6), true);
    });
    document.getElementById("bif-reset").addEventListener("click", () => setLevel(0, false));
    document.querySelectorAll("#bif-degree button").forEach((b) => b.addEventListener("click", () => {
      document.querySelectorAll("#bif-degree button").forEach((o) => { o.classList.toggle("on", o === b); o.setAttribute("aria-checked", o === b); });
      z = +b.dataset.z; setLevel(0, false);
    }));

    // drag to zoom, hover readout
    const toData = (ev) => {
      const r = cv.getBoundingClientRect();
      const u = (ev.clientX - r.left) / r.width, v = (ev.clientY - r.top) / r.height;
      const mu = win.m0 + u * (win.m1 - win.m0);
      const vv = flip ? 1 - v : v;
      return { u, v, mu, x: win.x1 - vv * (win.x1 - win.x0), px: ev.clientX - r.left, py: ev.clientY - r.top };
    };
    cv.addEventListener("pointerdown", (ev) => { drag = toData(ev); cv.setPointerCapture(ev.pointerId); });
    cv.addEventListener("pointermove", (ev) => {
      const p = toData(ev);
      tip.hidden = false;
      tip.style.left = `${Math.min(p.px + 18, cv.clientWidth - 150)}px`;
      tip.style.top = `${p.py + 14}px`;
      tip.textContent = `μ = ${p.mu.toPrecision(9)}   x = ${p.x.toPrecision(5)}`;
      if (drag) {
        draw(1);
        const ctx = cv.getContext("2d"), s = cv.width / cv.clientWidth;
        ctx.strokeStyle = css("--accent"); ctx.lineWidth = 1.5 * s;
        ctx.strokeRect(drag.px * s, drag.py * s, (p.px - drag.px) * s, (p.py - drag.py) * s);
      }
    });
    cv.addEventListener("pointerleave", () => { tip.hidden = true; });
    cv.addEventListener("pointerup", (ev) => {
      if (!drag) return;
      const p = toData(ev), a = drag; drag = null;
      if (Math.abs(p.px - a.px) < 6 || Math.abs(p.py - a.py) < 6) { draw(1); return; }
      win = { m0: Math.min(a.mu, p.mu), m1: Math.max(a.mu, p.mu), x0: Math.min(a.x, p.x), x1: Math.max(a.x, p.x) };
      flip = false; level = null; draw(1);
    });
    cv.addEventListener("dblclick", () => setLevel(0, false));

    size(); setLevel(0, false);
    return () => { size(); draw(1); };
  }

  // ------------------------------------------------------------ delta: ratio convergence
  function ratioPlot() {
    const traces = [];
    DEG.forEach((z) => {
      const m = F.superstable[z].mu, xs = [], ys = [];
      for (let n = 1; n + 1 < m.length; n++) { xs.push(n); ys.push((m[n] - m[n - 1]) / (m[n + 1] - m[n])); }
      const col = seriesColor(z);
      traces.push({ x: xs, y: ys, name: `x${SUP[z]}`, mode: "lines+markers", line: { color: col, width: 2 },
        marker: { size: 8, color: col, line: { color: css("--surface"), width: 2 } },
        hovertemplate: `x${SUP[z]}, n=%{x}<br>ratio %{y:.10f}<extra></extra>` });
      traces.push({ x: [0.5, xs.length + 0.5], y: [DELTA[z], DELTA[z]], mode: "lines", showlegend: false, hoverinfo: "skip",
        line: { color: col, width: 1.5, dash: "dash" } });
    });
    const ann = DEG.map((z) => ({ x: 1, xref: "paper", xanchor: "left", y: DELTA[z], text: `δ = ${DELTA[z].toFixed(6)}…`, showarrow: false,
      font: { size: 12, color: css("--text-2") } }));
    Plotly.react("ratio-plot", traces, baseLayout({
      margin: { l: 58, r: 120, t: 16, b: 48 }, annotations: ann,
      xaxis: Object.assign(baseLayout().xaxis, { title: { text: "doubling n" }, dtick: 1 }),
      yaxis: Object.assign(baseLayout().yaxis, { title: { text: "gap ratio" }, range: [0.5, 12.2] }),
    }), PCONF);
  }

  // ------------------------------------------------------------ renormalisation demo
  const EPS = [0, 1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2];
  const epsLabel = (e) => (e === 0 ? "0" : `10${String(Math.round(Math.log10(e))).replace("-", "⁻").replace(/\d/g, (d) => "⁰¹²³⁴⁵⁶⁷⁸⁹"[+d])}`);
  function renorm() {
    let z = 2;
    const kIn = document.getElementById("ren-k"), eIn = document.getElementById("ren-eps");
    const render = () => {
      const K = +kIn.value, eps = EPS[+eIn.value];
      document.getElementById("ren-k-val").textContent = K;
      document.getElementById("ren-eps-val").textContent = epsLabel(eps);
      const mu = MUINF[z] + eps;
      // F_0 = f, F_k(x) = F_{k-1}(F_{k-1}(lam_{k-1} x)) / lam_{k-1}
      const Fk = [(x) => 1 - mu * Math.pow(Math.abs(x), z)], lam = [Fk[0](1)];
      for (let k = 1; k <= K; k++) {
        const prev = Fk[k - 1], l = lam[k - 1];
        Fk.push((x) => prev(prev(l * x)) / l);
        lam.push(Fk[k](1));
      }
      const xs = []; for (let i = 0; i <= 240; i++) xs.push(-1 + i / 120);
      const g = gFun(z), col = seriesColor(z), traces = [];
      for (let k = 0; k <= K; k++) {
        const ys = xs.map((x) => { const v = Fk[k](x); return Number.isFinite(v) && Math.abs(v) < 3 ? v : null; });
        const cur = k === K;
        traces.push({ x: xs, y: ys, mode: "lines", name: cur ? `𝒯${k > 0 ? superscriptK(k) : "⁰"}f` : `𝒯^${k} f`, showlegend: cur,
          line: { color: col, width: cur ? 3 : 1.2 }, opacity: cur ? 1 : 0.18 + 0.5 * (k / Math.max(K, 1)),
          hovertemplate: `𝒯^${k} f(%{x:.3f}) = %{y:.6f}<extra></extra>` });
      }
      traces.push({ x: xs, y: xs.map(g), mode: "lines", name: "fixed point g", line: { color: css("--text"), width: 1.5, dash: "dash" },
        hovertemplate: "g(%{x:.3f}) = %{y:.6f}<extra></extra>" });
      Plotly.react("ren-plot", traces, baseLayout({
        xaxis: Object.assign(baseLayout().xaxis, { title: { text: "x" }, range: [-1, 1] }),
        yaxis: Object.assign(baseLayout().yaxis, { range: [-1.6, 1.4] }),
      }), PCONF);
    };
    function superscriptK(n) { return String(n).split("").map((d) => "⁰¹²³⁴⁵⁶⁷⁸⁹"[+d]).join(""); }
    kIn.addEventListener("input", render);
    eIn.addEventListener("input", render);
    document.querySelectorAll("#ren-degree button").forEach((b) => b.addEventListener("click", () => {
      document.querySelectorAll("#ren-degree button").forEach((o) => { o.classList.toggle("on", o === b); o.setAttribute("aria-checked", o === b); });
      z = +b.dataset.z; render();
    }));
    render();
    return render;
  }

  function gPlot() {
    const xs = []; for (let i = 0; i <= 240; i++) xs.push(-1 + i / 120);
    const traces = DEG.map((z) => ({
      x: xs, y: xs.map(gFun(z)), mode: "lines", name: `x${SUP[z]}: α = ${ALPHA[z].toFixed(6)}…`,
      line: { color: seriesColor(z), width: 2 }, hovertemplate: `degree ${z}: g(%{x:.3f}) = %{y:.6f}<extra></extra>`,
    }));
    Plotly.react("g-plot", traces, baseLayout({
      xaxis: Object.assign(baseLayout().xaxis, { title: { text: "x" } }),
      yaxis: Object.assign(baseLayout().yaxis, { title: { text: "g(x)" } }),
    }), PCONF);
  }

  // ------------------------------------------------------------ records
  const RECORDS = [
    { key: "delta2", z: 2, sym: "δ", prev: 1019, src: "Broadhurst 1999 (OEIS b-file)", reported: 10022, repSrc: "Molteni 2016, reported, digits unpublished", oeis: "A006890", now: 10000 },
    { key: "alpha2", z: 2, sym: "α", prev: 1018, src: "Broadhurst 1999 (OEIS b-file)", reported: 10026, repSrc: "Molteni 2016, reported, digits unpublished", oeis: "A006891", now: 10000 },
    { key: "delta4", z: 4, sym: "δ", prev: 325, src: "Burbanks, Osbaldestin & Thurlby 2021 (proven)", oeis: "new", now: 5000 },
    { key: "alpha4", z: 4, sym: "α", prev: 331, src: "Burbanks, Osbaldestin & Thurlby 2021 (proven)", oeis: "A195102", now: 5000 },
    { key: "delta6", z: 6, sym: "δ", prev: 46, src: "Briggs 1991", oeis: "new", now: 5000 },
    { key: "alpha6", z: 6, sym: "α", prev: 59, src: "Mathar 2010", oeis: "A396218", now: 5000 },
    { key: "delta8", z: 8, sym: "δ", prev: 38, src: "Briggs 1991", oeis: "new", now: 5000 },
    { key: "alpha8", z: 8, sym: "α", prev: 42, src: "Briggs 1991", oeis: "A396219", now: 5000 },
  ];
  const label = (r) => `${r.sym}, degree ${r.z}`;
  function recordsPlot() {
    const ys = RECORDS.map(label), traces = [];
    RECORDS.forEach((r) => traces.push({ x: [r.prev, r.now], y: [label(r), label(r)], mode: "lines", showlegend: false, hoverinfo: "skip",
      line: { color: css("--border"), width: 3 } }));
    traces.push({ x: RECORDS.map((r) => r.prev), y: ys, mode: "markers", name: "previous best published",
      marker: { size: 11, color: css("--surface"), line: { color: css("--muted"), width: 2 } },
      customdata: RECORDS.map((r) => r.src), hovertemplate: "%{y}<br>%{x:,} digits, %{customdata}<extra></extra>" });
    const rep = RECORDS.filter((r) => r.reported);
    traces.push({ x: rep.map((r) => r.reported), y: rep.map(label), mode: "markers", name: "reported, unpublished",
      marker: { size: 20, symbol: "diamond-open", color: css("--muted"), line: { width: 2 } },
      customdata: rep.map((r) => r.repSrc), hovertemplate: "%{y}<br>%{x:,} digits, %{customdata}<extra></extra>" });
    traces.push({ x: RECORDS.map((r) => r.now), y: ys, mode: "markers+text", name: "this work",
      marker: { size: 14, color: RECORDS.map((r) => seriesColor(r.z)), line: { color: css("--surface"), width: 2 } },
      text: RECORDS.map((r) => `${r.reported ? "public " : ""}×${Math.round(r.now / r.prev)}`), textposition: "middle right", textfont: { color: css("--text-2"), size: 12 },
      hovertemplate: "%{y}<br>%{x:,} digits, this work<extra></extra>" });
    Plotly.react("records-plot", traces, baseLayout({
      margin: { l: 120, r: 80, t: 30, b: 48 },
      xaxis: Object.assign(baseLayout().xaxis, { type: "log", title: { text: "decimal digits (log scale)" }, range: [1.3, 4.45],
        tickvals: [30, 100, 300, 1000, 3000, 10000], ticktext: ["30", "100", "300", "1,000", "3,000", "10,000"] }),
      yaxis: Object.assign(baseLayout().yaxis, { autorange: "reversed", gridcolor: "rgba(0,0,0,0)" }),
    }), PCONF);
    const tb = document.getElementById("records-body");
    tb.innerHTML = RECORDS.map((r) => {
      const v = F.digits[r.key].slice(0, 14) + "…";
      const oeis = r.oeis === "new" ? `<span class="badge new">new sequence</span>` :
        `<a href="https://oeis.org/${r.oeis}">${r.oeis}</a>`;
      return `<tr><td><span class="swatch" style="background:${seriesColor(r.z)}"></span>${label(r)}</td><td class="mono">${v}</td>
        <td class="num">${fmt(r.prev)}${r.reported ? `<br><span class="muted small">(${fmt(r.reported)} reported)</span>` : ""}</td>
        <td>${r.src}</td><td class="num"><b>${fmt(r.now)}</b></td><td>${oeis}</td></tr>`;
    }).join("");
  }

  // ------------------------------------------------------------ method: convergence + verification
  function convPlot() {
    const traces = DEG.map((z) => ({
      x: F.convergence[z].map((p) => p[0]), y: F.convergence[z].map((p) => p[1]), name: `x${SUP[z]}`, mode: "lines+markers",
      line: { color: seriesColor(z), width: 2 }, marker: { size: 8, color: seriesColor(z), line: { color: css("--surface"), width: 2 } },
      hovertemplate: `degree ${z}: %{x} nodes → %{y:.1f} correct digits<extra></extra>`,
    }));
    Plotly.react("conv-plot", traces, baseLayout({
      xaxis: Object.assign(baseLayout().xaxis, { title: { text: "collocation nodes n" } }),
      yaxis: Object.assign(baseLayout().yaxis, { title: { text: "correct digits of α" }, rangemode: "tozero" }),
    }), PCONF);
    const tb = document.getElementById("verify-body");
    const prevPub = { delta2: "all 1,019 (OEIS)", alpha2: "all 1,018 (OEIS)", delta4: "all 325 proven + all 60 (Briggs)", alpha4: "all 331 proven + all 64 (Briggs)",
      delta6: "all 46 (Briggs)", alpha6: "all 59 (Mathar) + 50 (Briggs)", delta8: "all 38 (Briggs)", alpha8: "all 42 (Briggs) + 28 (OEIS)" };
    tb.innerHTML = RECORDS.map((r) => {
      const run = F.runs[r.z], agree = r.sym === "δ" ? run.agree_delta : run.agree_alpha;
      const indep = r.sym === "δ" ? run.indep_delta : run.indep_alpha;
      return `<tr><td><span class="swatch" style="background:${seriesColor(r.z)}"></span>${label(r)}</td><td class="num">${fmt(r.now)}</td>
        <td class="num">${fmt(agree)}</td><td class="num">${fmt(indep)}</td><td>${prevPub[r.key]}</td></tr>`;
    }).join("");
  }

  // ------------------------------------------------------------ digit explorer
  function digits() {
    const sel = document.getElementById("dig-select"), box = document.getElementById("dig-box");
    const search = document.getElementById("dig-search"), found = document.getElementById("dig-found");
    RECORDS.forEach((r) => { const o = document.createElement("option"); o.value = r.key; o.textContent = `${label(r)}  (${fmt(r.now)} digits)`; sel.appendChild(o); });
    const render = () => {
      const s = F.digits[sel.value], dot = s.indexOf("."), ip = s.slice(0, dot), fr = s.slice(dot + 1);
      const q = search.value.replace(/\D/g, "");
      const hits = [];
      if (q.length >= 1) { let i = fr.indexOf(q); while (i !== -1 && hits.length < 500) { hits.push(i); i = fr.indexOf(q, i + 1); } }
      const marked = new Uint8Array(fr.length);
      hits.forEach((h) => { for (let j = h; j < h + q.length; j++) marked[j] = 1; });
      const out = [`<div class="row"><span class="idx"></span><span>${ip}.</span></div>`];
      for (let row = 0; row < fr.length; row += 50) {
        let line = "";
        for (let g = row; g < Math.min(row + 50, fr.length); g += 10) {
          let grp = "";
          for (let j = g; j < Math.min(g + 10, fr.length); j++) grp += marked[j] ? `<mark>${fr[j]}</mark>` : fr[j];
          line += `<span class="grp">${grp.replace(/<\/mark><mark>/g, "")}</span>`;
        }
        out.push(`<div class="row"><span class="idx">${row + 1}</span><span>${line}</span></div>`);
      }
      box.innerHTML = out.join("");
      found.textContent = q ? (hits.length ? `${hits.length}${hits.length === 500 ? "+" : ""} match${hits.length === 1 ? "" : "es"}, first at decimal place ${fmt(hits[0] + 1)}` : "not found") : "";
      if (hits.length) { const m = box.querySelector("mark"); if (m) box.scrollTop = m.offsetTop - box.clientHeight / 3; }
    };
    sel.addEventListener("change", render);
    search.addEventListener("input", render);
    document.getElementById("dig-copy").addEventListener("click", async (ev) => {
      try { await navigator.clipboard.writeText(F.digits[sel.value]); ev.target.textContent = "Copied"; } catch { ev.target.textContent = "Copy failed"; }
      setTimeout(() => { ev.target.textContent = "Copy"; }, 1400);
    });
    document.getElementById("dig-download").addEventListener("click", () => {
      const r = RECORDS.find((x) => x.key === sel.value), s = F.digits[r.key], dot = s.indexOf(".");
      const ds = s.replace(".", ""), lines = [];
      for (let i = 0; i < ds.length; i++) lines.push(`${i + dot} ${ds[i]}`);
      const name = r.oeis === "new" ? `${r.key.replace(/\d/, "")}_degree${r.z}.txt` : `b${r.oeis.slice(1)}.txt`;
      const blob = new Blob([lines.join("\n") + "\n"], { type: "text/plain" });
      const a = document.createElement("a"); a.href = URL.createObjectURL(blob); a.download = name; a.click();
      setTimeout(() => URL.revokeObjectURL(a.href), 2000);
    });
    render();
  }

  // ------------------------------------------------------------ boot
  heroStreams();
  const redrawHero = heroCanvas();
  const redrawBif = bifurcation();
  const redrawRen = renorm();
  ratioPlot(); gPlot(); recordsPlot(); convPlot(); digits();
  const rethemed = () => { redrawHero(); redrawBif(); redrawRen(); ratioPlot(); gPlot(); recordsPlot(); convPlot(); };
  window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", rethemed);
  let rt = null;
  window.addEventListener("resize", () => { clearTimeout(rt); rt = setTimeout(() => { redrawHero(); redrawBif(); }, 200); });
})();
