/* Visualizations: protein map, sequence view, 3D structure, alignment. */
/* References
 *
 * Cheng, J., Novati, G., Pan, J., Bycroft, C., Žemgulytė, A., Applebaum, T., Pritzel, A., Wong,
 *     L. H., Zielinski, M., Sargeant, T., Schneider, R. G., Senior, A. W., Jumper, J., Hassabis, D.,
 *     Kohli, P., & Avsec, Ž. (2023). Accurate proteome-wide missense variant effect prediction with
 *     AlphaMissense. Science, 381(6664), Article eadg7492. https://doi.org/10.1126/science.adg7492
 *
 * Jumper, J., Evans, R., Pritzel, A., Green, T., Figurnov, M., Ronneberger, O., Tunyasuvunakool,
 *     K., Bates, R., Žídek, A., Potapenko, A., Bridgland, A., Meyer, C., Kohl, S. A. A., Ballard, A.
 *     J., Cowie, A., Romera-Paredes, B., Nikolov, S., Jain, R., Adler, J., . . . Hassabis, D. (2021).
 *     Highly accurate protein structure prediction with AlphaFold. Nature, 596(7873), 583–589.
 *     https://doi.org/10.1038/s41586-021-03819-2
 *
 * Rego, N., & Koes, D. (2015). 3Dmol.js: Molecular visualization with WebGL. Bioinformatics,
 *     31(8), 1322–1324. https://doi.org/10.1093/bioinformatics/btu829
 */

function useWidth(ref) {
  const [w, setW] = React.useState(0);
  React.useLayoutEffect(() => {
    setW(ref.current.getBoundingClientRect().width);
    const ro = new ResizeObserver(([e]) => setW(e.contentRect.width));
    ro.observe(ref.current);
    return () => ro.disconnect();
  }, []);
  return w;
}
const cssVar = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();

/** Every candidate peptide span [start, end] of a site, one per repeat. */
function siteSpans(s) {
  if (!s.peptide_start || !s.positions || !s.positions.length) return [];
  const offset = s.positions[0] - s.peptide_start;
  const len = s.peptide_end - s.peptide_start + 1;
  return s.positions.map((p) => [p - offset, p - offset + len - 1]);
}

/** Greedy interval packing -> row index per item (items sorted by start). */
function packRows(items, gap) {
  const ends = [];
  return items.map((it) => {
    let r = ends.findIndex((e) => it.x0 > e + gap);
    if (r < 0) { r = ends.length; ends.push(0); }
    ends[r] = it.x1;
    return r;
  });
}

function niceStep(span, px) {
  const raw = span / Math.max(2, px / 90);
  const mag = 10 ** Math.floor(Math.log10(raw));
  return [1, 2, 5, 10].map((m) => m * mag).find((s) => s >= raw);
}

/* ---------------------------- Protein map ---------------------------- */
const fitLabel = (text, px) => { const n = Math.floor(px / 5.6); return text.length <= n ? text : text.slice(0, Math.max(1, n - 1)) + "…"; };
const GUTTER = 124, RIGHT = 16;

function ProteinMap({ length, sequence, sites, focus, annotations, onSelect }) {
  const wrap = useRef();
  const width = useWidth(wrap);
  const [view, setView] = useState([0, length]);
  const [hover, setHover] = useState(null);
  const [brush, setBrush] = useState(null);
  useEffect(() => setView([0, length]), [length]);

  const plotW = Math.max(10, width - GUTTER - RIGHT);
  const [v0, v1] = view;
  const k = plotW / (v1 - v0);
  const x = (p) => GUTTER + (p - v0) * k;          // p in residue coordinates (residue i spans i-1..i)
  const visible = (a, b) => b >= v0 && a - 1 <= v1;
  const focusSite = sites.find((s) => s.id === focus);
  const focusPos = focusSite && focusSite.positions[0];

  // One lollipop per position.
  const groups = useMemo(() => {
    const m = new Map();
    sites.forEach((s) => s.positions.forEach((p) => {
      const g = m.get(p) || { pos: p, saaps: [], obs: 0 };
      g.saaps.push(s); g.obs += s.n_observations || 0; m.set(p, g);
    }));
    return [...m.values()].sort((a, b) => a.pos - b.pos);
  }, [sites]);
  const maxObs = Math.max(1, ...groups.map((g) => g.obs));

  const feats = (annotations && annotations.uniprot.features) || [];
  const af = annotations && annotations.alphafold;
  const plddt = af && af.plddt;
  const amMean = af && af.am_mean;

  const zoomTo = (a, b) => {
    const span = Math.max(30, b - a);
    const c = (a + b) / 2;
    const lo = Math.max(0, Math.min(length - span, c - span / 2));
    setView([lo, Math.min(length, lo + span)]);
  };
  const zoom = (f) => zoomTo(v0 + (v1 - v0) * (1 - f) / 2, v1 - (v1 - v0) * (1 - f) / 2);
  const pan = (d) => { const s = (v1 - v0) * d; const lo = Math.max(0, Math.min(length - (v1 - v0), v0 + s)); setView([lo, lo + (v1 - v0)]); };
  const posAt = (mx) => Math.min(length, Math.max(1, Math.floor(v0 + (mx - GUTTER) / k) + 1));
  const localX = (e) => e.clientX - wrap.current.getBoundingClientRect().left;

  /* ---- tracks ---- */
  const tracks = [];
  let y = 0;
  const add = (label, h, render) => { tracks.push({ label, y, h }); const node = render(y); y += h; return node; };
  const nodes = [];

  // axis
  nodes.push(add("", 24, (ty) => {
    const step = niceStep(v1 - v0, plotW);
    const ticks = [];
    for (let t = Math.ceil(v0 / step) * step || step; t <= v1; t += step) ticks.push(t);
    return <g key="axis">{ticks.map((t) => (
      <g key={t}><line className="tick" x1={x(t - 0.5)} x2={x(t - 0.5)} y1={ty + 16} y2={ty + 22} />
        <text className="axis" x={x(t - 0.5)} y={ty + 12} textAnchor="middle">{t}</text></g>))}</g>;
  }));

  // SAAP sites
  nodes.push(add("SAAP sites", 70, (ty) => {
    const base = ty + 66;
    const items = groups.filter((g) => visible(g.pos, g.pos));
    const sorted = [...items].sort((a, b) => (a.saaps.some((s) => s.id === focus) ? 1 : 0) - (b.saaps.some((s) => s.id === focus) ? 1 : 0));
    return <g key="sites">
      <line className="grid" x1={GUTTER} x2={GUTTER + plotW} y1={base} y2={base} />
      {sorted.map((g) => {
        const isFocus = g.saaps.some((s) => s.id === focus);
        const cls = isFocus ? " focus" : focus ? " dim" : "";
        const cx = x(g.pos - 0.5);
        const h = 12 + 44 * Math.log1p(g.obs) / Math.log1p(maxObs);
        const r = isFocus ? 5.5 : 3 + Math.min(3, Math.log2(g.saaps.length));
        return <g key={g.pos}>
          <line className={"stem" + cls} x1={cx} x2={cx} y1={base} y2={base - h} />
          <circle className={"head" + cls} cx={cx} cy={base - h} r={r}
                  onMouseDown={(e) => e.stopPropagation()}
                  onClick={() => onSelect && onSelect(g)} />
        </g>;
      })}
    </g>;
  }));

  // backbone + peptide span + residues when zoomed in
  nodes.push(add("Sequence", 22, (ty) => {
    const showRes = k >= 9 && sequence;
    const res = [];
    if (showRes) for (let i = Math.max(1, Math.floor(v0)); i <= Math.min(length, Math.ceil(v1)); i++) {
      res.push(<text key={i} className={"res" + (i === focusPos ? " focus" : "")} x={x(i - 0.5)} y={ty + 15}>{sequence[i - 1]}</text>);
    }
    return <g key="bb">
      <rect className="backbone" x={x(Math.max(0, v0))} width={x(Math.min(length, v1)) - x(Math.max(0, v0))} y={ty + 3} height={16} rx={2} opacity={showRes ? 0.5 : 1} />
      {focusSite && siteSpans(focusSite).map(([a, b], i) => (
        <rect key={i} className="span" x={x(a - 1)} width={Math.max(2, (b - a + 1) * k)} y={ty + (showRes ? 19 : 8)} height={showRes ? 2 : 6} rx={1} />))}
      {res}
    </g>;
  }));

  // UniProt span features, one packed track per category
  [["domain", "Domains"], ["region", "Regions"], ["topology", "Topology"]].forEach(([track, label]) => {
    const items = feats.filter((f) => f.track === track && f.end > f.start && visible(f.start, f.end))
      .map((f) => ({ ...f, x0: x(f.start - 1), x1: x(f.end) }));
    if (!items.length) return;
    const rows = packRows(items, 2);
    const n = Math.max(...rows) + 1;
    nodes.push(add(label, n * 18 + 4, (ty) => <g key={track}>{items.map((f, i) => {
      const w = Math.max(1, f.x1 - f.x0);
      const fy = ty + 2 + rows[i] * 18;
      const name = f.description || f.type;
      return <g key={i}>
        <rect className={"feat " + track} x={f.x0} y={fy} width={w} height={14} rx={2} />
        {w > 24 && <text className="feat-label" x={f.x0 + 4} y={fy + 10.5}>{fitLabel(name, w - 8)}</text>}
      </g>;
    })}</g>));
  });

  // Point features
  const points = feats.filter((f) => (f.track === "ptm" || f.track === "site") && visible(f.start, f.start));
  if (points.length) {
    // Label a point only when its text cannot collide with the previous label.
    let lastEnd = -Infinity;
    const labels = points.map((f) => {
      const text = fitLabel((f.description || f.type).replace(/^\(.*?\)\s*/, ""), 90);
      const cx = x(f.start - 0.5), w = text.length * 5.6;
      if (k < 3 || cx - w / 2 < lastEnd + 6) return null;
      lastEnd = cx + w / 2;
      return text;
    });
    const any = labels.some(Boolean);
    nodes.push(add("PTMs & sites", any ? 30 : 16, (ty) => <g key="ptm">{points.map((f, i) => (
      <g key={i}><path className={"ptm " + f.track} d={`M${x(f.start - 0.5)},${ty + 3} l3.5,8 h-7 z`} />
        {labels[i] && <text className="feat-label" x={x(f.start - 0.5)} y={ty + 24} textAnchor="middle">{labels[i]}</text>}</g>))}</g>));
  }
  const vars = feats.filter((f) => f.track === "variant" && visible(f.start, f.end));
  if (vars.length) nodes.push(add("Known variants", 16, (ty) => <g key="var">{vars.map((f, i) => (
    <line key={i} className="var" x1={x(f.start - 0.5)} x2={x(f.start - 0.5)} y1={ty + 3} y2={ty + 13} />))}</g>));

  // Per-residue heat tracks, binned to ~2px when zoomed out
  const heat = (key, label, values, color) => {
    if (!values) return;
    nodes.push(add(label, 16, (ty) => {
      const lo = Math.max(0, Math.floor(v0)), hi = Math.min(length, Math.ceil(v1));
      const per = Math.max(1, Math.ceil((hi - lo) / (plotW / 2)));
      const rects = [];
      for (let i = lo; i < hi; i += per) {
        const vals = values.slice(i, i + per).filter((v) => v != null);
        if (!vals.length) continue;
        const mean = vals.reduce((a, b) => a + b, 0) / vals.length;
        rects.push(<rect key={i} x={x(i)} width={Math.max(1, per * k) + 0.3} y={ty + 3} height={10} fill={color(mean)} />);
      }
      return <g key={key}>{rects}</g>;
    }));
  };
  heat("plddt", "AlphaFold pLDDT", plddt, plddtHex);
  heat("am", "AlphaMissense", amMean, amColor);  // average pathogenicity per residue

  const H = y + 4;

  /* ---- interaction ---- */
  const onMove = (e) => {
    const mx = localX(e);
    if (brush) setBrush({ ...brush, x1: mx });
    setHover(mx >= GUTTER && mx <= GUTTER + plotW ? { mx, my: e.clientY - wrap.current.getBoundingClientRect().top, pos: posAt(mx) } : null);
  };
  const onDown = (e) => { const mx = localX(e); if (mx >= GUTTER) setBrush({ x0: mx, x1: mx }); };
  const onUp = () => {
    if (brush && Math.abs(brush.x1 - brush.x0) > 4) {
      const a = posAt(Math.min(brush.x0, brush.x1)) - 1, b = posAt(Math.max(brush.x0, brush.x1));
      zoomTo(a, b);
    }
    setBrush(null);
  };

  const tip = hover && !brush && (() => {
    const p = hover.pos;
    const here = groups.find((g) => g.pos === p);
    const fs = feats.filter((f) => f.start <= p && p <= f.end && f.track !== "variant");
    const vs = feats.filter((f) => f.track === "variant" && f.start === p);
    return (
      <div className="tip" style={{ left: Math.min(hover.mx + 12, width - 300), top: hover.my + 14 }}>
        <b>{sequence ? sequence[p - 1] : ""}{p}</b>
        {plddt && plddt[p - 1] != null && <span className="t"> · pLDDT {plddt[p - 1].toFixed(1)}</span>}
        {amMean && amMean[p - 1] != null && <span className="t"> · AlphaMissense {amMean[p - 1].toFixed(2)} (mean of 19 substitutions)</span>}
        {here && here.saaps.slice(0, 6).map((s) => <div key={s.id}>{s.aa_sub} · {s.mtp_seq} · {plural(s.n_observations, "observation")}</div>)}
        {here && here.saaps.length > 6 && <div className="t">+{here.saaps.length - 6} more</div>}
        {fs.slice(0, 4).map((f, i) => <div key={i} className="t">{f.type}: {f.description || "—"}</div>)}
        {vs.slice(0, 3).map((f, i) => <div key={"v" + i} className="t">Variant {f.ref}→{f.alt} {f.description.split(";")[0]}</div>)}
      </div>
    );
  })();

  const zoomed = v0 > 0 || v1 < length;
  return (
    <div>
      <div className="card-head">
        <h2>Protein map</h2>
        <span className="muted">{length.toLocaleString()} aa · {plural(groups.length, "site")}{zoomed ? ` · ${Math.floor(v0) + 1}–${Math.ceil(v1)}` : ""}</span>
        <div className="map-tools">
          {focusPos && <button className="sm ghost" onClick={() => zoomTo(focusPos - 30, focusPos + 30)}>Zoom to site</button>}
          {zoomed && <button className="ghost icon sm" onClick={() => pan(-0.5)} title="Pan left"><Icon name="left" /></button>}
          {zoomed && <button className="ghost icon sm" onClick={() => pan(0.5)} title="Pan right"><Icon name="right" /></button>}
          <button className="ghost icon sm" onClick={() => zoom(0.5)} title="Zoom in"><Icon name="zoomIn" /></button>
          <button className="ghost icon sm" onClick={() => zoom(2)} disabled={!zoomed} title="Zoom out"><Icon name="zoomOut" /></button>
          <button className="ghost icon sm" onClick={() => setView([0, length])} disabled={!zoomed} title="Reset"><Icon name="reset" /></button>
        </div>
      </div>
      <div className="pmap" ref={wrap} onMouseMove={onMove} onMouseLeave={() => { setHover(null); setBrush(null); }}
           onMouseDown={onDown} onMouseUp={onUp} onDoubleClick={() => setView([0, length])} style={{ padding: "8px 0" }}>
        {width > 0 && (
          <svg height={H} style={{ cursor: "crosshair" }}>
            {tracks.filter((t) => t.label).map((t) => (
              <text key={t.label} className="tlabel" x={16} y={t.y + t.h / 2 + 4}>{t.label}</text>))}
            {nodes}
            {hover && !brush && <line className="cross" x1={x(hover.pos - 0.5)} x2={x(hover.pos - 0.5)} y1={20} y2={H} />}
            {brush && <rect className="brush" x={Math.min(brush.x0, brush.x1)} width={Math.abs(brush.x1 - brush.x0)} y={0} height={H} />}
          </svg>
        )}
        {tip}
      </div>
      <MapLegend annotations={annotations} />
    </div>
  );
}

function MapLegend({ annotations }) {
  if (!annotations) return <div className="legend"><span><Spinner />&nbsp;Loading UniProt and AlphaFold tracks</span></div>;
  const { uniprot, alphafold } = annotations;
  return (
    <div className="legend">
      <span>Drag to zoom · double-click to reset</span>
      {alphafold.plddt && <span className="lg"><b>Model confidence</b>{PLDDT_BANDS.map(([t, name], i) => <span key={t}><i className="dot" style={{ background: PLDDT_HEX[i] }} />{name} ({PLDDT_RANGES[i]})</span>)}</span>}
      {alphafold.am_mean && <span className="lg"><b>AlphaMissense</b>{[[0.1, "Likely benign"], [0.45, "Ambiguous"], [0.9, "Likely pathogenic"]].map(([v, name]) => (
        <span key={name}><i className="dot" style={{ background: amColor(v) }} />{name}</span>))}</span>}
      {uniprot.issue && <span>UniProt features: {uniprot.issue}</span>}
      {alphafold.issue && <span>AlphaFold: {alphafold.issue}</span>}
    </div>
  );
}

/* ---------------------------- Sequence view ---------------------------- */
function SequenceView({ sequence, sites, focus, plddt }) {
  const [mode, setMode] = useState("sites");
  const marks = useMemo(() => {
    const m = new Array(sequence.length).fill("");
    const set = (i, c) => { if (i >= 1 && i <= m.length) m[i - 1] = c; };
    sites.forEach((s) => siteSpans(s).forEach(([a, b]) => { for (let i = a; i <= b; i++) if (!m[i - 1]) set(i, "pep"); }));
    sites.forEach((s) => s.positions.forEach((p) => set(p, "site")));
    const f = sites.find((s) => s.id === focus);
    if (f) {
      siteSpans(f).forEach(([a, b]) => { for (let i = a; i <= b; i++) if (m[i - 1] !== "site") set(i, "fpep"); });
      f.positions.forEach((p) => set(p, "focus"));
    }
    return m;
  }, [sequence, sites, focus]);

  const lines = [];
  for (let s = 0; s < sequence.length; s += 60) {
    const blocks = [];
    for (let b = s; b < Math.min(s + 60, sequence.length); b += 10) {
      const res = [];
      for (let i = b; i < Math.min(b + 10, sequence.length); i++) {
        const style = mode === "plddt" && plddt && plddt[i] != null ? { background: plddtHex(plddt[i]) + "66" } : null;
        res.push(<span key={i} className={"r " + (mode === "sites" ? marks[i] : "")} style={style}
                       title={`${sequence[i]}${i + 1}${plddt && plddt[i] != null ? ` · pLDDT ${plddt[i].toFixed(1)}` : ""}`}>{sequence[i]}</span>);
      }
      blocks.push(<Fragment key={b}>{res}{" "}</Fragment>);
    }
    lines.push(<div className="ln" key={s}><span className="n">{s + 1}</span>{blocks}</div>);
  }
  return (
    <div>
      <div className="card-head">
        <h2>Sequence</h2><span className="muted">{sequence.length.toLocaleString()} aa</span>
        <span className="spacer" />
        <Segmented value={mode} onChange={setMode} options={[["sites", "SAAP sites"], ...(plddt ? [["plddt", "pLDDT"]] : []), ["plain", "Plain"]]} />
      </div>
      <div className="seqview">{lines}</div>
    </div>
  );
}

/* ---------------------------- 3D structure ---------------------------- */
function StructureView({ accession, sites, focus }) {
  const el = useRef();
  const viewer = useRef();
  const [status, setStatus] = useState("loading");
  const [color, setColor] = useState("plddt");
  const focusSite = sites.find((s) => s.id === focus);
  const focusPos = focusSite && focusSite.positions[0];

  const paint = () => {
    const v = viewer.current;
    if (!v) return;
    v.removeAllLabels();
    v.setStyle({}, { cartoon: color === "plddt" ? { colorfunc: (a) => plddtHex(a.b) } : { color: "#c8c8cc" } });
    const others = [...new Set(sites.flatMap((s) => s.positions))].filter((p) => p !== focusPos);
    if (others.length) v.addStyle({ resi: others, atom: "CA" }, { sphere: { radius: 1.0, color: cssVar("--site") } });
    if (focusPos) {
      v.addStyle({ resi: focusPos }, { stick: { radius: 0.3, color: cssVar("--mark") } });
      v.addStyle({ resi: focusPos, atom: "CA" }, { sphere: { radius: 1.6, color: cssVar("--mark") } });
      v.addLabel(`${focusSite.ref}${focusPos}${focusSite.alt}`, {
        position: { resi: focusPos, atom: "CA" }, fontSize: 12, fontColor: "white",
        backgroundColor: cssVar("--mark"), backgroundOpacity: 0.95, borderRadius: 3,
      });
    }
    v.render();
  };
  const frame = () => {
    const v = viewer.current;
    if (!v) return;
    if (focusPos) { v.zoomTo({ resi: focusPos }); v.zoom(0.25); } else v.zoomTo();
    v.render();
  };

  useEffect(() => {
    let live = true;
    setStatus("loading");
    api.structure(accession).then((pdb) => {
      if (!live) return;
      const v = $3Dmol.createViewer(el.current, { backgroundColor: cssVar("--surface"), antialias: true });
      v.addModel(pdb, "pdb");
      viewer.current = v;
      paint(); frame();
      setStatus("ready");
    }).catch((e) => live && setStatus(/getParameter|WebGL/i.test(e.message) ? "3D view requires WebGL" : e.message));
    return () => { live = false; if (viewer.current) viewer.current.clear(); viewer.current = null; };
  }, [accession]);
  useEffect(paint, [color, focus, sites]);

  return (
    <div className="structure">
      <div className="viewer" ref={el} />
      {status === "loading" && <Loading label="Loading AlphaFold model" />}
      {status !== "loading" && status !== "ready" && <div className="empty">{status}</div>}
      {status === "ready" && (
        <Fragment>
          <div className="tools">
            <Segmented value={color} onChange={setColor} options={[["plddt", "pLDDT"], ["plain", "Plain"]]} />
            <button className="sm" onClick={frame}>Reset view</button>
          </div>
          <div className="overlay">
            {color === "plddt" && PLDDT_BANDS.map(([t, name], i) => <span key={t}><i className="dot" style={{ background: PLDDT_HEX[i] }} />{name} ({PLDDT_RANGES[i]})</span>)}
            <span><i className="dot" style={{ background: cssVar("--mark"), borderRadius: 4 }} />This SAAP</span>
            <span><i className="dot" style={{ background: cssVar("--site"), borderRadius: 4 }} />Other sites</span>
          </div>
        </Fragment>
      )}
    </div>
  );
}

/* ------------------------------ Alignment ------------------------------ */
/** BLAST-style alignment of BP and SAAP against the protein context. */
function Alignment({ saap, sequence }) {
  const bp = saap.bp_seq || "", mtp = saap.mtp_seq || "";
  const start = sequence ? saap.peptide_start : null;
  const FLANK = 12;
  const lo = start ? Math.max(1, start - FLANK) : 1;
  const hi = start ? Math.min(sequence.length, start + bp.length - 1 + FLANK) : 0;
  const pad = " ".repeat(start ? start - lo : 0);
  const head = (label, n) => label.padEnd(9) + String(n ?? "").padStart(5) + "  ";
  const end = (n) => (n ? <span className="end">{n}</span> : null);
  const last = (seq) => (start ? start + seq.length - 1 : null);
  return (
    <div className="aln">
      {start && <div>{head("Protein", lo)}
        <span className="ctx">{sequence.slice(lo - 1, start - 1)}</span>
        <span className="pep">{sequence.slice(start - 1, start - 1 + bp.length)}</span>
        <span className="ctx">{sequence.slice(start - 1 + bp.length, hi)}</span>{end(hi)}</div>}
      <div>{head("BP", start)}{pad}{bp}{end(last(bp))}</div>
      <div className="idents">{head("", null)}{pad}{[...bp].map((c, i) => (c === mtp[i] ? "|" : " ")).join("")}</div>
      <div>{head("SAAP", start)}{pad}{[...mtp].map((c, i) => (c !== bp[i] ? <span key={i} className="mm">{c}</span> : c))}{end(last(mtp))}</div>
    </div>
  );
}

/* ------------------------- Substitution matrix ------------------------- */
// Residues grouped by side-chain chemistry, so related swaps sit together.
const MATRIX_ORDER = "GAVLIMPFWYSTCNQDEKRH";
const MATRIX_GROUPS = [[0, 7, "Hydrophobic"], [7, 10, "Aromatic"], [10, 15, "Polar"], [15, 17, "Acidic"], [17, 20, "Basic"]];

function SubstitutionMatrix({ data, onSelect }) {
  const [hover, setHover] = useState(null);
  const wrap = useRef();
  if (!data) return <Loading />;
  const aa = MATRIX_ORDER.split("");
  const n = (r, a) => (data.counts[r] && data.counts[r][a]) || 0;
  const max = Math.max(1, ...aa.flatMap((r) => aa.map((a) => n(r, a))));
  const t = (v) => Math.log1p(v) / Math.log1p(max);  // log scale: H->D dwarfs everything else
  const mix = (f) => `color-mix(in oklab, var(--accent) ${Math.round(12 + 88 * f)}%, var(--surface))`;

  // Geometry: 26px cells, an 8px gutter between chemical classes.
  const C = 22, G = 8, L = 178, T = 92;
  const group = (i) => MATRIX_GROUPS.findIndex(([s, e]) => i >= s && i < e);
  const pos = (i) => i * C + group(i) * G;
  const span = pos(19) + C;
  const legendX = L + span + 30, W = legendX + 100, H = T + span + 12;
  const ticks = [1, 10, 100, 1000, 10000].filter((v) => v < max * 0.8).concat([max]);
  const ty = (v) => T + span - t(v) * span;
  const hoverAt = (e, r, a, v) => {
    const box = wrap.current.getBoundingClientRect();
    setHover({ r, a, v, x: e.clientX - box.left, y: e.clientY - box.top });
  };

  return (
    <div className="matrix" ref={wrap} onMouseLeave={() => setHover(null)}>
      <svg viewBox={`0 0 ${W} ${H}`} style={{ maxWidth: W }}>
        <defs>
          <linearGradient id="mx-ramp" x1="0" y1="1" x2="0" y2="0">
            {[0, 0.25, 0.5, 0.75, 1].map((f) => <stop key={f} offset={f} style={{ stopColor: mix(f) }} />)}
          </linearGradient>
        </defs>
        {/* axis titles */}
        <text className="mx-title" x={L + span / 2} y={16}>Substituted residue</text>
        <text className="mx-title" transform={`translate(14 ${T + span / 2}) rotate(-90)`}>Reference residue</text>
        {/* class names and residue letters */}
        {MATRIX_GROUPS.map(([s, e, name]) => (
          <Fragment key={name}>
            <text className="mx-group" x={L + (pos(s) + pos(e - 1) + C) / 2} y={44}>{name}</text>
            <line className="mx-bracket" x1={L + pos(s) + 2} x2={L + pos(e - 1) + C - 2} y1={52} y2={52} />
            <text className="mx-group left" x={L - 42} y={T + (pos(s) + pos(e - 1) + C) / 2 + 5}>{name}</text>
            <line className="mx-bracket" x1={L - 32} x2={L - 32} y1={T + pos(s) + 2} y2={T + pos(e - 1) + C - 2} />
          </Fragment>
        ))}
        {aa.map((a, j) => <text key={"c" + a} className="mx-axis" x={L + pos(j) + C / 2} y={T - 12}>{a}</text>)}
        {aa.map((r, i) => <text key={"r" + r} className="mx-axis" x={L - 15} y={T + pos(i) + C / 2 + 5}>{r}</text>)}
        {/* cells */}
        {aa.map((r, i) => aa.map((a, j) => {
          const v = n(r, a), diag = r === a;
          return <rect key={r + a} x={L + pos(j) + 1} y={T + pos(i) + 1} width={C - 2} height={C - 2} rx={3}
                       className={"mx-cell" + (diag ? " diag" : v ? " has" : "")} style={diag ? null : { fill: v ? mix(t(v)) : "var(--surface-2)" }}
                       onMouseMove={(e) => !diag && hoverAt(e, r, a, v)}
                       onClick={() => !diag && v && onSelect && onSelect(r, a)} />;
        }))}
        {/* legend: vertical log-scale bar beside the matrix */}
        <rect x={legendX} y={T} width={16} height={span} rx={4} fill="url(#mx-ramp)" />
        {ticks.map((v) => (
          <Fragment key={v}>
            <line className="mx-tick" x1={legendX + 16} x2={legendX + 22} y1={ty(v)} y2={ty(v)} />
            <text className="mx-legend-label" x={legendX + 27} y={ty(v) + 4.5}>{v.toLocaleString()}</text>
          </Fragment>
        ))}
        <text className="mx-title" transform={`translate(${legendX + 90} ${T + span / 2}) rotate(90)`}># SAAPs</text>
      </svg>
      {hover && <div className="tip" style={{ left: hover.x + 14, top: hover.y + 14 }}><b>{hover.r} → {hover.a}</b> · {plural(hover.v, "SAAP")}</div>}
    </div>
  );
}

/* ------------------------------- Small ------------------------------- */
function MiniMap({ length, positions, width = 160 }) {
  if (!length) return DASH;
  const sx = (p) => 1 + ((p - 0.5) / length) * (width - 2);
  return (
    <svg className="minimap" width={width} height={14}>
      <line className="bb" x1={1} x2={width - 1} y1={10} y2={10} />
      {positions.map((p) => <line key={p} className="s" x1={sx(p)} x2={sx(p)} y1={2} y2={10} />)}
    </svg>
  );
}

function BarList({ items }) {
  if (!items || !items.length) return <div className="empty">No data</div>;
  const max = Math.max(...items.map((i) => i.n), 1);
  return (
    <div className="bars">
      {items.map((i) => (
        <div className="bar-row" key={i.label}>
          <div className="bar-label" title={i.label}>{i.label}</div>
          <div className="bar-track"><div className="bar-fill" style={{ width: `${(i.n / max) * 100}%` }} /></div>
          <div className="bar-value">{i.n.toLocaleString()}</div>
        </div>
      ))}
    </div>
  );
}
