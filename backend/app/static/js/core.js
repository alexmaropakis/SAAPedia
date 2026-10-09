/* Shared primitives: API client, formatting, routing, small UI components. */
/* References
 *
 * Cheng, J., Novati, G., Pan, J., Bycroft, C., Žemgulytė, A., Applebaum, T., Pritzel, A., Wong,
 *     L. H., Zielinski, M., Sargeant, T., Schneider, R. G., Senior, A. W., Jumper, J., Hassabis, D.,
 *     Kohli, P., & Avsec, Ž. (2023). Accurate proteome-wide missense variant effect prediction with
 *     AlphaMissense. Science, 381(6664), Article eadg7492. https://doi.org/10.1126/science.adg7492
 *
 * Grantham, R. (1974). Amino acid difference formula to help explain protein evolution. Science,
 *     185(4154), 862–864. https://doi.org/10.1126/science.185.4154.862
 *
 * Henikoff, S., & Henikoff, J. G. (1992). Amino acid substitution matrices from protein blocks.
 *     Proceedings of the National Academy of Sciences, 89(22), 10915–10919.
 *     https://doi.org/10.1073/pnas.89.22.10915
 *
 * Jumper, J., Evans, R., Pritzel, A., Green, T., Figurnov, M., Ronneberger, O., Tunyasuvunakool,
 *     K., Bates, R., Žídek, A., Potapenko, A., Bridgland, A., Meyer, C., Kohl, S. A. A., Ballard, A.
 *     J., Cowie, A., Romera-Paredes, B., Nikolov, S., Jain, R., Adler, J., . . . Hassabis, D. (2021).
 *     Highly accurate protein structure prediction with AlphaFold. Nature, 596(7873), 583–589.
 *     https://doi.org/10.1038/s41586-021-03819-2
 *
 * Li, W.-H., Wu, C.-I., & Luo, C.-C. (1984). Nonrandomness of point mutation as reflected in
 *     nucleotide substitutions in pseudogenes and its evolutionary implications. Journal of Molecular
 *     Evolution, 21(1), 58–71. https://doi.org/10.1007/BF02100628
 *
 * Varadi, M., Bertoni, D., Magana, P., Paramval, U., Pidruchna, I., Radhakrishnan, M., Tsenkov,
 *     M., Nair, S., Mirdita, M., Yeo, J., Kovalevskiy, O., Tunyasuvunakool, K., Laydon, A., Žídek,
 *     A., Tomlinson, H., Hariharan, D., Abrahamson, J., Green, T., Jumper, J., . . . Velankar, S.
 *     (2024). AlphaFold Protein Structure Database in 2024: Providing structure coverage for over 214
 *     million protein sequences. Nucleic Acids Research, 52(D1), D368–D375.
 *     https://doi.org/10.1093/nar/gkad1011
 */
const { useState, useEffect, useCallback, useMemo, useRef, Fragment } = React;

/* ------------------------------- API ------------------------------- */
async function request(url, opts = {}) {
  const r = await fetch(url, opts);
  if (!r.ok) {
    const body = await r.json().catch(() => ({}));
    throw new Error(body.detail || `Request failed (${r.status})`);
  }
  return r;
}
const getJSON = async (url) => (await request(url)).json();
const postJSON = async (url, body) => request(url, {
  method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body || {}),
});
const qs = (params) => {
  const q = new URLSearchParams();
  Object.entries(params).forEach(([k, v]) => { if (v !== "" && v != null) q.set(k, v); });
  return q.toString();
};

const api = {
  stats: () => getJSON("/api/stats"),
  facets: () => getJSON("/api/facets"),
  datasets: () => getJSON("/api/datasets"),
  curation: () => getJSON("/api/curation"),
  list: (params) => getJSON("/api/saap?" + qs(params)),
  saap: (id) => getJSON(`/api/saap/${id}`),
  proteins: (params) => getJSON("/api/proteins?" + qs(params)),
  protein: (acc) => getJSON(`/api/proteins/${acc}`),
  proteinAnnotations: (acc) => getJSON(`/api/proteins/${acc}/annotations`),
  structure: async (acc) => (await request(`/api/proteins/${acc}/structure.pdb`)).text(),
  annotateStatus: () => getJSON("/api/annotate/status"),
  annotate: async (payload) => (await postJSON("/api/annotate", payload)).json(),
  remove: async (payload) => (await postJSON("/api/saap/delete", payload)).json(),
  upload: async (file) => {
    const fd = new FormData();
    fd.append("file", file);
    return (await request("/api/upload", { method: "POST", body: fd })).json();
  },
  exportFile: async (kind, payload, reference) => {
    let r;
    if (kind === "fasta") {
      const fd = new FormData();
      fd.append("payload", JSON.stringify(payload));
      if (reference) fd.append("reference", reference);
      r = await request("/api/export/fasta", { method: "POST", body: fd });
    } else {
      r = await postJSON(`/api/export/${kind}`, payload);
    }
    const name = (r.headers.get("Content-Disposition") || "").match(/filename="(.+)"/);
    return { blob: await r.blob(), name: name ? name[1] : `saap.${kind}`, skipped: Number(r.headers.get("X-SAAP-Skipped") || 0) };
  },
};

function download(blob, name) {
  const url = URL.createObjectURL(blob);
  const a = Object.assign(document.createElement("a"), { href: url, download: name });
  document.body.appendChild(a); a.click(); a.remove();
  URL.revokeObjectURL(url);
}

/* ---------------------------- formatting --------------------------- */
const DASH = <span className="dash">—</span>;
const fmt = {
  int: (v) => (v == null ? DASH : Number(v).toLocaleString()),
  num: (v, d = 2) => (v == null ? DASH : Number(v).toFixed(d)),
  sci: (v) => (v == null ? DASH : v === 0 ? "0" : Number(v).toExponential(1)),
  bool: (v) => (v === true ? <span className="yes">Yes</span> : v === false ? <span className="no">No</span> : DASH),
  text: (v) => (v == null || v === "" ? DASH : v),
};
const plural = (n, word) => `${Number(n).toLocaleString()} ${word}${n === 1 ? "" : "s"}`;
const firstOf = (v) => (v || "").split(";").map((s) => s.trim()).find(Boolean) || "";

const AA3 = { A: "Ala", R: "Arg", N: "Asn", D: "Asp", C: "Cys", Q: "Gln", E: "Glu", G: "Gly", H: "His", I: "Ile",
  L: "Leu", K: "Lys", M: "Met", F: "Phe", P: "Pro", S: "Ser", T: "Thr", W: "Trp", Y: "Tyr", V: "Val" };

/** Substituted peptide with the changed residue marked. */
function SubSeq({ saap, bp }) {
  if (!saap) return DASH;
  if (!bp || bp.length !== saap.length) return saap;
  const diff = [...saap].findIndex((c, i) => c !== bp[i]);
  if (diff < 0) return saap;
  return <Fragment>{saap.slice(0, diff)}<span className="sub">{saap[diff]}</span>{saap.slice(diff + 1)}</Fragment>;
}

/* -------------------------- score semantics ------------------------ */
const PLDDT_BANDS = [
  [90, "Very high", "var(--plddt-vh)"], [70, "Confident", "var(--plddt-h)"],
  [50, "Low", "var(--plddt-l)"], [0, "Very low", "var(--plddt-vl)"],
];
const PLDDT_RANGES = ["pLDDT > 90", "90 > pLDDT > 70", "70 > pLDDT > 50", "pLDDT < 50"];
const plddtBand = (v) => PLDDT_BANDS.find(([t]) => v >= t);
const PLDDT_HEX = ["#0053d6", "#65cbf3", "#ffdb13", "#ff7d45"];
const plddtHex = (v) => PLDDT_HEX[PLDDT_BANDS.findIndex(([t]) => v >= t)];

const AM_CLASS = { LBen: "Likely benign", Amb: "Ambiguous", LPath: "Likely pathogenic" };
/** AlphaMissense pathogenicity -> colour by its published classes:
 *  likely benign < 0.34 (blue), ambiguous 0.34–0.564 (grey), likely pathogenic > 0.564 (red). */
const AM_THRESHOLDS = [0.34, 0.564];
function amColor(v) {
  if (v == null) return "transparent";
  const mix = (a, b, t) => `rgb(${a.map((x, i) => Math.round(x + (b[i] - x) * t)).join(",")})`;
  const [blue, grey, red] = [[59, 130, 246], [190, 190, 196], [220, 38, 38]];
  if (v < AM_THRESHOLDS[0]) return mix(blue, grey, (v / AM_THRESHOLDS[0]) * 0.6);
  if (v <= AM_THRESHOLDS[1]) return `rgb(${grey.join(",")})`;
  return mix(grey, red, 0.4 + 0.6 * (v - AM_THRESHOLDS[1]) / (1 - AM_THRESHOLDS[1]));
}
/** Grantham classes per Li et al. 1984. */
const granthamClass = (g) => g == null ? "" : g <= 50 ? "Conservative" : g <= 100 ? "Moderately conservative"
  : g <= 150 ? "Moderately radical" : "Radical";

/* ------------------------------ routing ---------------------------- */
function parseHash() {
  const [path, query] = (location.hash.replace(/^#\/?/, "") || "home").split("?");
  const [page, param] = path.split("/");
  return { page, param: param ? decodeURIComponent(param) : null, query: new URLSearchParams(query || "") };
}
function useRoute() {
  const [route, setRoute] = useState(parseHash);
  useEffect(() => {
    const on = () => { setRoute(parseHash()); window.scrollTo(0, 0); };
    window.addEventListener("hashchange", on);
    return () => window.removeEventListener("hashchange", on);
  }, []);
  return route;
}
const go = (path) => { location.hash = "#/" + path; };

/** Fetch on mount / when deps change; returns [data, error, reload].
 *  `keep` holds the previous result while reloading (tables), instead of
 *  clearing it (detail pages, where stale data would be misleading). */
function useAsync(fn, deps, keep = false) {
  const [state, setState] = useState({ data: null, error: null });
  const [tick, setTick] = useState(0);
  const prevDeps = useRef(deps);
  useEffect(() => {
    let live = true;
    const sameTarget = deps.every((d, i) => d === prevDeps.current[i]);
    prevDeps.current = deps;
    setState((s) => ({ data: keep || sameTarget ? s.data : null, error: null }));
    fn().then((data) => live && setState({ data, error: null }))
      .catch((error) => live && setState({ data: null, error }));
    return () => { live = false; };
  }, [...deps, tick]);
  return [state.data, state.error, () => setTick((t) => t + 1)];
}

function useStored(key, initial) {
  const [v, setV] = useState(() => {
    try { const s = localStorage.getItem(key); return s ? JSON.parse(s) : initial; } catch { return initial; }
  });
  useEffect(() => { try { localStorage.setItem(key, JSON.stringify(v)); } catch { /* storage unavailable */ } }, [key, v]);
  return [v, setV];
}

/* ----------------------------- components -------------------------- */
const ICONS = {
  search: "M11 19a8 8 0 1 0 0-16 8 8 0 0 0 0 16zM21 21l-4.3-4.3",
  filter: "M3 5h18M6 12h12M10 19h4",
  x: "M18 6 6 18M6 6l12 12",
  down: "M6 9l6 6 6-6",
  left: "M15 18l-6-6 6-6",
  right: "M9 18l6-6-6-6",
  columns: "M4 4h16v16H4zM10 4v16M16 4v16",
  download: "M12 3v12M7 10l5 5 5-5M5 21h14",
  external: "M14 4h6v6M20 4l-9 9M18 14v5a1 1 0 0 1-1 1H5a1 1 0 0 1-1-1V7a1 1 0 0 1 1-1h5",
  zoomIn: "M11 19a8 8 0 1 0 0-16 8 8 0 0 0 0 16zM21 21l-4.3-4.3M11 8v6M8 11h6",
  zoomOut: "M11 19a8 8 0 1 0 0-16 8 8 0 0 0 0 16zM21 21l-4.3-4.3M8 11h6",
  reset: "M3 12a9 9 0 1 0 3-6.7L3 8M3 3v5h5",
  upload: "M12 21V9M7 14l5-5 5 5M5 3h14",
  spark: "M12 3v4M12 17v4M3 12h4M17 12h4M5.6 5.6l2.8 2.8M15.6 15.6l2.8 2.8M5.6 18.4l2.8-2.8M15.6 8.4l2.8-2.8",
};
/** "SAAPedia" with SAAP in the accent; one letter, picked once per page load, is the "substitution". */
const SUB_LETTER = Math.floor(Math.random() * 8);
const Wordmark = () => (
  <span className="wordmark">{[..."SAAPedia"].map((c, i) => (
    <span key={i} className={i === SUB_LETTER ? "wm-sub" : i < 4 ? "wm-saap" : ""}>{c}</span>
  ))}</span>
);

const Icon = ({ name }) => <svg className="i" viewBox="0 0 24 24"><path d={ICONS[name]} /></svg>;

const Spinner = () => <span className="spin" />;
const Loading = ({ label = "Loading" }) => <div className="loading"><Spinner />{label}</div>;
const ErrorNote = ({ error }) => <div className="empty">{error.message}</div>;

function Menu({ label, children, disabled, className = "" }) {
  const [open, setOpen] = useState(false);
  const ref = useRef();
  useEffect(() => {
    if (!open) return;
    const close = (e) => { if (!ref.current.contains(e.target)) setOpen(false); };
    document.addEventListener("mousedown", close);
    return () => document.removeEventListener("mousedown", close);
  }, [open]);
  return (
    <div className="menu-wrap" ref={ref}>
      <button className={className} disabled={disabled} onClick={() => setOpen((o) => !o)}>{label}<Icon name="down" /></button>
      {open && <div className="menu" onClick={(e) => { if (e.target.closest("button")) setOpen(false); }}>{children}</div>}
    </div>
  );
}

function Segmented({ value, options, onChange }) {
  return (
    <div className="segmented">
      {options.map(([v, label]) => (
        <button key={v} className={value === v ? "on" : ""} onClick={() => onChange(v)}>{label}</button>
      ))}
    </div>
  );
}

function Modal({ title, onClose, children, footer }) {
  useEffect(() => {
    const esc = (e) => e.key === "Escape" && onClose();
    document.addEventListener("keydown", esc);
    return () => document.removeEventListener("keydown", esc);
  }, [onClose]);
  return (
    <div className="backdrop" onMouseDown={onClose}>
      <div className="modal" onMouseDown={(e) => e.stopPropagation()}>
        <div className="card-head"><h2>{title}</h2><span className="spacer" />
          <button className="ghost icon" onClick={onClose}><Icon name="x" /></button></div>
        <div className="card-body">{children}</div>
        {footer && <div className="modal-foot">{footer}</div>}
      </div>
    </div>
  );
}

function Tabs({ value, onChange, tabs }) {
  return (
    <div className="tabs">
      {tabs.map(([v, label, count]) => (
        <button key={v} className={value === v ? "on" : ""} onClick={() => onChange(v)}>
          {label}{count != null && <span className="count">&nbsp;{count}</span>}
        </button>
      ))}
    </div>
  );
}

function Pager({ page, pageSize, total, onPage, onPageSize }) {
  const pages = Math.max(1, Math.ceil(total / pageSize));
  const from = total ? (page - 1) * pageSize + 1 : 0;
  return (
    <div className="pager">
      <span>{from.toLocaleString()}–{Math.min(page * pageSize, total).toLocaleString()} of {total.toLocaleString()}</span>
      <span className="spacer" />
      <span>Rows</span>
      <select value={pageSize} onChange={(e) => onPageSize(Number(e.target.value))}>
        {[25, 50, 100, 200, 500].map((n) => <option key={n} value={n}>{n}</option>)}
      </select>
      <button className="ghost icon sm" disabled={page <= 1} onClick={() => onPage(page - 1)}><Icon name="left" /></button>
      <span>{page} / {pages}</span>
      <button className="ghost icon sm" disabled={page >= pages} onClick={() => onPage(page + 1)}><Icon name="right" /></button>
    </div>
  );
}

/** Sortable table header cell. */
function Th({ col, sort, onSort }) {
  const active = sort.key === col.key;
  return (
    <th className={(col.sortable ? "sort " : "") + (col.num ? "num" : "")} title={col.title}
        onClick={() => col.sortable && onSort(col.key)}>
      {col.label}{active && <span className="arrow">{sort.order === "asc" ? "↑" : "↓"}</span>}
    </th>
  );
}
const toggleSort = (key) => (s) => (s.key === key ? { key, order: s.order === "asc" ? "desc" : "asc" } : { key, order: "desc" });

const ToastContext = React.createContext(() => {});
const useToast = () => React.useContext(ToastContext);
