const { useState, useEffect, useCallback, useRef } = React;

/* ------------------------------ write auth ------------------------------ */
let writePassword = sessionStorage.getItem("saap_write_pw") || "";
function setWritePassword(p) {
  writePassword = p || "";
  sessionStorage.setItem("saap_write_pw", writePassword);
}
function hasWritePassword() { return !!writePassword; }

// fetch for mutating calls: attaches the password, and on 401 prompts once and retries.
async function writeFetch(url, opts = {}) {
  const build = () => {
    const headers = Object.assign({}, opts.headers, { "X-Write-Password": writePassword });
    return Object.assign({}, opts, { headers });
  };
  let r = await fetch(url, build());
  if (r.status === 401) {
    const p = window.prompt("Enter the write password to make changes:");
    if (p === null) return r;
    setWritePassword(p);
    r = await fetch(url, build());
  }
  return r;
}

/* ----------------------------- API helpers ----------------------------- */
const api = {
  async config() { return (await fetch("/api/config")).json(); },
  async stats() { return (await fetch("/api/stats")).json(); },
  async facets() { return (await fetch("/api/facets")).json(); },
  async datasets() { return (await fetch("/api/datasets")).json(); },
  async list(params) {
    const qs = new URLSearchParams();
    Object.entries(params).forEach(([k, v]) => {
      if (v !== "" && v !== null && v !== undefined) qs.set(k, v);
    });
    const r = await fetch("/api/saap?" + qs.toString());
    if (!r.ok) throw new Error("Failed to load SAAP list");
    return r.json();
  },
  async detail(id) { return (await fetch("/api/saap/" + id)).json(); },
  async upload(file, doiMap) {
    const fd = new FormData();
    fd.append("file", file);
    if (doiMap && Object.keys(doiMap).length) fd.append("dataset_doi_map", JSON.stringify(doiMap));
    const r = await writeFetch("/api/upload", { method: "POST", body: fd });
    const body = await r.json();
    if (!r.ok) throw new Error(body.detail || "Upload failed");
    return body;
  },
  async saveDatasetDois(map) {
    const r = await writeFetch("/api/datasets", {
      method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ map }),
    });
    if (!r.ok) throw new Error((await r.json().catch(() => ({}))).detail || "Failed to save DOIs");
    return r.json();
  },
  async deleteSaap(payload) {
    const r = await writeFetch("/api/saap/delete", {
      method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload),
    });
    const body = await r.json();
    if (!r.ok) throw new Error(body.detail || "Delete failed");
    return body;
  },
  async exportFasta(payload) {
    const r = await fetch("/api/export/fasta", {
      method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload),
    });
    if (!r.ok) {
      const b = await r.json().catch(() => ({}));
      throw new Error(b.detail || "Export failed");
    }
    return r.blob();
  },
};

/* ------------------------------ formatting ------------------------------ */
const fmt = {
  num(v, d = 2) {
    if (v === null || v === undefined) return "—";
    if (Math.abs(v) !== 0 && (Math.abs(v) < 1e-3 || Math.abs(v) >= 1e6)) return v.toExponential(2);
    return Number(v).toFixed(d);
  },
  pep(v) { return v === null || v === undefined ? "—" : v.toExponential(1); },
  bool(v) {
    if (v === true) return <span className="flag-yes">yes</span>;
    if (v === false) return <span className="flag-no">no</span>;
    return <span className="flag-no">—</span>;
  },
};

/* -------------------------------- Toast --------------------------------- */
function useToast() {
  const [toast, setToast] = useState(null);
  const show = useCallback((msg, err = false) => {
    setToast({ msg, err });
    setTimeout(() => setToast(null), 3500);
  }, []);
  const node = toast ? <div className={"toast" + (toast.err ? " err" : "")}>{toast.msg}</div> : null;
  return [node, show];
}

/* ------------------------------- Columns -------------------------------- */
const COLUMNS = [
  { key: "mtp_seq", label: "SAAP", sortable: true, cls: "seq" },
  { key: "bp_seq", label: "Base peptide", sortable: true, cls: "seq" },
  { key: "aa_sub", label: "AAS", sortable: true },
  { key: "source_gene", label: "Genes", sortable: true },
  { key: "ref_proteins", label: "RefProteins", sortable: true },
  { key: "source_accession", label: "UniProt", sortable: true },
  { key: "species", label: "Species" },
  { key: "n_datasets", label: "Datasets", sortable: true },
  { key: "digests", label: "Digest" },
  { key: "acquisition_types", label: "Acquisition" },
  { key: "n_observations", label: "Obs", sortable: true, cls: "num" },
  { key: "best_saap_pep", label: "Best PEP", sortable: true, cls: "num" },
  { key: "max_positional_probability", label: "Max PosProb", sortable: true, cls: "num" },
  { key: "max_evidence_fragments", label: "Max Frags", sortable: true, cls: "num" },
  { key: "immunoglobulin", label: "Ig" },
  { key: "trypsin", label: "Tryp" },
  { key: "missed_cleavage", label: "MissClv" },
  { key: "aas_at_peptide_terminus", label: "TermAAS" },
  { key: "greater_than_shared", label: ">Shared" },
];

function chips(arr, cls) {
  return arr && arr.length ? arr.map((v) => <span key={v} className={cls}>{v}</span>) : "—";
}

function renderCell(col, row, ctx) {
  switch (col.key) {
    case "mtp_seq":
      return <span className="link" onClick={() => ctx.openDetail(row.id)}>{row.mtp_seq}</span>;
    case "aa_sub":
      return <span className="sub-chip">{row.aa_sub || "—"}</span>;
    case "n_datasets":
      if (!row.datasets.length) return "—";
      return row.datasets.map((d) => {
        const url = ctx.datasetUrl[d];
        return url
          ? <a key={d} className="ds-chip linked" href={url} target="_blank" rel="noreferrer">{d}</a>
          : <span key={d} className="ds-chip">{d}</span>;
      });
    case "digests": return chips(row.digests, "sub-chip");
    case "species": return chips(row.species, "sub-chip");
    case "acquisition_types": return chips(row.acquisition_types, "sub-chip");
    case "best_saap_pep": return fmt.pep(row.best_saap_pep);
    case "max_positional_probability": return fmt.num(row.max_positional_probability, 3);
    case "max_evidence_fragments": return row.max_evidence_fragments ?? "—";
    case "immunoglobulin": return fmt.bool(row.immunoglobulin);
    case "trypsin": return fmt.bool(row.trypsin);
    case "missed_cleavage": return fmt.bool(row.missed_cleavage);
    case "aas_at_peptide_terminus": return fmt.bool(row.aas_at_peptide_terminus);
    case "greater_than_shared": return fmt.bool(row.greater_than_shared);
    default: return row[col.key] ?? "—";
  }
}

/* ------------------------------ filter bits ----------------------------- */
function Facet({ label, value, opts, onChange }) {
  return (
    <div className="field">
      <label>{label}</label>
      <select value={value} onChange={(e) => onChange(e.target.value)}>
        <option value="">All</option>
        {(opts || []).map((o) => <option key={o} value={o}>{o}</option>)}
      </select>
    </div>
  );
}

function TriState({ label, value, onChange }) {
  return (
    <div className="field">
      <label>{label}</label>
      <select value={value} onChange={(e) => onChange(e.target.value)}>
        <option value="">Any</option><option value="true">Yes</option><option value="false">No</option>
      </select>
    </div>
  );
}

/* --------------------------------- App ---------------------------------- */
function App() {
  const [tab, setTab] = useState("browse");
  const [stats, setStats] = useState(null);
  const [datasets, setDatasets] = useState([]);
  const [facets, setFacets] = useState({ datasets: [], digests: [], species: [], acquisition_types: [], aa_subs: [] });
  const [cfg, setCfg] = useState({ write_protected: false });
  const [unlocked, setUnlocked] = useState(hasWritePassword());
  const [toastNode, showToast] = useToast();

  const refreshMeta = useCallback(() => {
    api.stats().then(setStats).catch(() => {});
    api.facets().then(setFacets).catch(() => {});
    api.datasets().then(setDatasets).catch(() => {});
  }, []);
  useEffect(() => { refreshMeta(); api.config().then(setCfg).catch(() => {}); }, [refreshMeta]);

  const unlock = () => {
    const p = window.prompt("Enter the write password to enable editing:");
    if (p === null) return;
    setWritePassword(p);
    setUnlocked(!!p);
    if (p) showToast("Editing unlocked");
  };
  const lock = () => { setWritePassword(""); setUnlocked(false); showToast("Editing locked"); };

  const datasetUrl = {};
  datasets.forEach((d) => { if (d.url) datasetUrl[d.name] = d.url; });

  return (
    <div className="app">
      <header className="top">
        <div className="brand">
          <h1>SAAP Database</h1>
          {cfg.write_protected && (
            unlocked
              ? <span className="lockpill unlocked" onClick={lock} title="Click to lock">Editing unlocked · lock</span>
              : <span className="lockpill" onClick={unlock} title="Click to enter write password">Read-only · unlock to edit</span>
          )}
        </div>
        <div className="stats">
          <div className="stat"><div className="num">{stats ? stats.n_saap : "—"}</div><div className="lbl">SAAP</div></div>
          <div className="stat"><div className="num">{stats ? stats.n_observations : "—"}</div><div className="lbl">Observations</div></div>
          <div className="stat"><div className="num">{stats ? stats.n_datasets : "—"}</div><div className="lbl">Datasets</div></div>
          <div className="stat"><div className="num">{stats ? stats.n_genes : "—"}</div><div className="lbl">Genes</div></div>
        </div>
      </header>

      <div className="tabs">
        <button className={"tab" + (tab === "browse" ? " active" : "")} onClick={() => setTab("browse")}>Browse</button>
        <button className={"tab" + (tab === "import" ? " active" : "")} onClick={() => setTab("import")}>Import</button>
        <button className={"tab" + (tab === "datasets" ? " active" : "")} onClick={() => setTab("datasets")}>Datasets</button>
      </div>

      {tab === "browse" && (
        <BrowseTab facets={facets} datasetUrl={datasetUrl} onDataChanged={refreshMeta} showToast={showToast} />
      )}
      {tab === "import" && (
        <ImportTab knownDatasets={facets.datasets} onIngested={() => { refreshMeta(); }} showToast={showToast} />
      )}
      {tab === "datasets" && (
        <DatasetsTab datasets={datasets} onChanged={refreshMeta} showToast={showToast} />
      )}

      {toastNode}
    </div>
  );
}

/* ------------------------------ Browse tab ------------------------------ */
function BrowseTab({ facets, datasetUrl, onDataChanged, showToast }) {
  const [data, setData] = useState({ items: [], total: 0, page: 1, page_size: 50 });
  const [loading, setLoading] = useState(false);
  const EMPTY_FILTERS = {
    q: "", dataset: "", digest: "", species: "", acquisition_type: "", aa_sub: "",
    immunoglobulin: "", trypsin: "", min_pos_prob: "", max_pep: "",
  };
  const [filters, setFilters] = useState(EMPTY_FILTERS);
  const [sort, setSort] = useState({ key: "n_observations", order: "desc" });
  const [page, setPage] = useState(1);
  const pageSize = 50;
  const [selected, setSelected] = useState(() => new Set());
  const [detailId, setDetailId] = useState(null);
  const [showExport, setShowExport] = useState(false);

  const load = useCallback(() => {
    setLoading(true);
    api.list({ ...filters, sort: sort.key, order: sort.order, page, page_size: pageSize })
      .then(setData).catch((e) => showToast(e.message, true)).finally(() => setLoading(false));
  }, [filters, sort, page, showToast]);
  useEffect(() => { load(); }, [load]);

  const setFilter = (k, v) => { setFilters((f) => ({ ...f, [k]: v })); setPage(1); };
  const clearFilters = () => { setFilters(EMPTY_FILTERS); setPage(1); };
  const toggleSort = (col) => {
    if (!col.sortable) return;
    setSort((s) => s.key === col.key ? { key: col.key, order: s.order === "asc" ? "desc" : "asc" } : { key: col.key, order: "desc" });
  };
  const toggleRow = (id) => setSelected((prev) => {
    const next = new Set(prev); next.has(id) ? next.delete(id) : next.add(id); return next;
  });
  const allOnPage = data.items.length > 0 && data.items.every((r) => selected.has(r.id));
  const togglePage = () => setSelected((prev) => {
    const next = new Set(prev);
    if (allOnPage) data.items.forEach((r) => next.delete(r.id));
    else data.items.forEach((r) => next.add(r.id));
    return next;
  });

  const deleteSelected = async () => {
    if (!selected.size) return;
    if (!window.confirm(`Delete ${selected.size} SAAP and all their observations? This cannot be undone.`)) return;
    try {
      const res = await api.deleteSaap({ ids: Array.from(selected) });
      showToast(`Deleted ${res.deleted} SAAP`);
      setSelected(new Set());
      onDataChanged();
      load();
    } catch (e) { showToast(e.message, true); }
  };

  const totalPages = Math.max(1, Math.ceil(data.total / pageSize));
  const activeFilters = Object.values(filters).filter((v) => v !== "").length;
  const ctx = { openDetail: setDetailId, datasetUrl };

  return (
    <React.Fragment>
      <div className="card">
        <div className="filters">
          <div className="field">
            <label>Search</label>
            <input type="search" placeholder="peptide, gene, protein, UniProt"
                   value={filters.q} onChange={(e) => setFilter("q", e.target.value)} />
          </div>
          <Facet label="Dataset" value={filters.dataset} opts={facets.datasets} onChange={(v) => setFilter("dataset", v)} />
          <Facet label="Digest" value={filters.digest} opts={facets.digests} onChange={(v) => setFilter("digest", v)} />
          <Facet label="Species" value={filters.species} opts={facets.species} onChange={(v) => setFilter("species", v)} />
          <Facet label="Acquisition" value={filters.acquisition_type} opts={facets.acquisition_types} onChange={(v) => setFilter("acquisition_type", v)} />
          <Facet label="AAS" value={filters.aa_sub} opts={facets.aa_subs} onChange={(v) => setFilter("aa_sub", v)} />
          <TriState label="Ig" value={filters.immunoglobulin} onChange={(v) => setFilter("immunoglobulin", v)} />
          <TriState label="Trypsin" value={filters.trypsin} onChange={(v) => setFilter("trypsin", v)} />
          <div className="field">
            <label>Min PosProb</label>
            <input type="number" step="0.01" min="0" max="1" style={{ width: 90 }}
                   value={filters.min_pos_prob} onChange={(e) => setFilter("min_pos_prob", e.target.value)} />
          </div>
          <div className="field">
            <label>Max PEP</label>
            <input type="number" step="any" min="0" style={{ width: 90 }}
                   value={filters.max_pep} onChange={(e) => setFilter("max_pep", e.target.value)} />
          </div>
          <button className="ghost" onClick={clearFilters} disabled={activeFilters === 0}>Clear</button>
        </div>
      </div>

      <div className="actionbar">
        <div className="selcount"><b>{selected.size}</b> selected</div>
        <div className="btn-group">
          <button className="danger" disabled={selected.size === 0} onClick={deleteSelected}>Delete selected</button>
          <button className="ghost" disabled={selected.size === 0} onClick={() => setShowExport({ mode: "selected" })}>Export selected</button>
          <button disabled={data.total === 0} onClick={() => setShowExport({ mode: "filtered" })}>
            Export {activeFilters ? "filtered " : "all "}({data.total})
          </button>
        </div>
      </div>

      <div className="table-wrap">
        <table>
          <thead>
            <tr>
              <th className="checkcol"><input type="checkbox" checked={allOnPage} onChange={togglePage} /></th>
              {COLUMNS.map((c) => (
                <th key={c.key} className={c.sortable ? "sortable" : ""} onClick={() => toggleSort(c)}>
                  {c.label}{sort.key === c.key && <span className="arrow"> {sort.order === "asc" ? "▲" : "▼"}</span>}
                </th>
              ))}
            </tr>
          </thead>
          <tbody>
            {loading ? (
              <tr><td colSpan={COLUMNS.length + 1} className="spinner">Loading…</td></tr>
            ) : data.items.length === 0 ? (
              <tr><td colSpan={COLUMNS.length + 1} className="empty">No SAAP found. Import a CSV, or adjust filters.</td></tr>
            ) : data.items.map((r) => (
              <tr key={r.id} className={selected.has(r.id) ? "selected" : ""}>
                <td className="checkcol"><input type="checkbox" checked={selected.has(r.id)} onChange={() => toggleRow(r.id)} /></td>
                {COLUMNS.map((c) => <td key={c.key} className={c.cls || ""}>{renderCell(c, r, ctx)}</td>)}
              </tr>
            ))}
          </tbody>
        </table>
      </div>

      <div className="pager">
        <button className="ghost small" disabled={page <= 1} onClick={() => setPage((p) => p - 1)}>Prev</button>
        <span>Page {data.page} of {totalPages} · {data.total} SAAP</span>
        <button className="ghost small" disabled={page >= totalPages} onClick={() => setPage((p) => p + 1)}>Next</button>
      </div>

      {detailId && <DetailDrawer id={detailId} datasetUrl={datasetUrl} onClose={() => setDetailId(null)} />}
      {showExport && (
        <ExportModal mode={showExport.mode} selected={selected} filters={filters}
          total={showExport.mode === "filtered" ? data.total : selected.size}
          onClose={() => setShowExport(false)} showToast={showToast} />
      )}
    </React.Fragment>
  );
}

/* ------------------------------ Import tab ------------------------------ */
function ImportTab({ knownDatasets, onIngested, showToast }) {
  const [drag, setDrag] = useState(false);
  const [busy, setBusy] = useState(false);
  const [result, setResult] = useState(null);
  const [rows, setRows] = useState([{ dataset: "", doi: "" }]);
  const inputRef = useRef();

  const setRow = (i, field, val) => setRows((rs) => rs.map((r, idx) => idx === i ? { ...r, [field]: val } : r));
  const addRow = () => setRows((rs) => [...rs, { dataset: "", doi: "" }]);
  const removeRow = (i) => setRows((rs) => rs.filter((_, idx) => idx !== i));

  const buildMap = () => {
    const map = {};
    rows.forEach((r) => { if (r.dataset.trim()) map[r.dataset.trim()] = r.doi.trim(); });
    return map;
  };

  const handleFile = async (file) => {
    if (!file) return;
    setBusy(true);
    try {
      const res = await api.upload(file, buildMap());
      setResult(res);
      showToast(`Imported ${res.rows_read} rows → ${res.saap_created} new SAAP`);
      onIngested();
    } catch (e) { showToast(e.message, true); }
    finally { setBusy(false); }
  };

  return (
    <React.Fragment>
      <div className="card">
        <h2>Import CSV</h2>
        <div className="desc">Columns are auto-mapped and de-duplicated. Peptides without a UniProt ID are dropped on import. "N datasets" is computed from the datasets each SAAP appears in.</div>
        <div className={"drop" + (drag ? " drag" : "")}
             onClick={() => inputRef.current.click()}
             onDragOver={(e) => { e.preventDefault(); setDrag(true); }}
             onDragLeave={() => setDrag(false)}
             onDrop={(e) => { e.preventDefault(); setDrag(false); handleFile(e.dataTransfer.files[0]); }}>
          {busy ? "Importing…" : <span><strong>Drop a .csv or .xlsx here</strong> or click to browse</span>}
          <input ref={inputRef} type="file" accept=".csv,.tsv,.txt,.xlsx" hidden onChange={(e) => handleFile(e.target.files[0])} />
        </div>

        {result && (
          <div className="result">
            <span className="good">Imported {result.filename}</span>
            <div className="badge-row">
              <span className="mini">rows read <b>{result.rows_read}</b></span>
              <span className="mini">new SAAP <b>{result.saap_created}</b></span>
              <span className="mini">observations <b>{result.observations_created}</b></span>
              <span className="mini">duplicates skipped <b>{result.duplicate_observations_skipped}</b></span>
              <span className="mini">DOIs saved <b>{result.dataset_dois_saved}</b></span>
              {result.saap_removed_no_uniprot > 0 && <span className="mini warn">removed (no UniProt ID) <b>{result.saap_removed_no_uniprot}</b></span>}
            </div>
            {result.columns_unmapped && result.columns_unmapped.length > 0 &&
              <div className="badge-row"><span className="mini warn">unmapped columns: {result.columns_unmapped.join(", ")}</span></div>}
          </div>
        )}
      </div>

      <div className="card">
        <h2>Dataset → paper DOI (optional)</h2>
        <div className="desc">Attach a paper to each dataset. Applied when you import above; datasets can also be edited later under the Datasets tab.</div>
        <div className="doi-rows">
          {rows.map((r, i) => (
            <div className="doi-row" key={i}>
              <input className="ds-name" list="known-datasets" placeholder="dataset name"
                     value={r.dataset} onChange={(e) => setRow(i, "dataset", e.target.value)} />
              <input placeholder="DOI or URL (e.g. 10.1038/s41586-020-0000-0)"
                     value={r.doi} onChange={(e) => setRow(i, "doi", e.target.value)} />
              <button className="ghost small" onClick={() => removeRow(i)} disabled={rows.length === 1}>Remove</button>
            </div>
          ))}
          <datalist id="known-datasets">{knownDatasets.map((d) => <option key={d} value={d} />)}</datalist>
          <div><button className="ghost small" onClick={addRow}>Add row</button></div>
        </div>
      </div>
    </React.Fragment>
  );
}

/* ----------------------------- Datasets tab ----------------------------- */
function DatasetsTab({ datasets, onChanged, showToast }) {
  const [edits, setEdits] = useState({});
  const [saving, setSaving] = useState(false);

  useEffect(() => {
    const init = {};
    datasets.forEach((d) => { init[d.name] = d.doi || ""; });
    setEdits(init);
  }, [datasets]);

  const save = async () => {
    setSaving(true);
    try {
      await api.saveDatasetDois(edits);
      showToast("Saved dataset DOIs");
      onChanged();
    } catch (e) { showToast(e.message, true); }
    finally { setSaving(false); }
  };

  const wipeAll = async () => {
    if (!window.confirm("Delete ALL SAAP and observations from the database? Dataset DOIs are kept. This cannot be undone.")) return;
    try {
      const res = await api.deleteSaap({ all: true });
      showToast(`Cleared ${res.deleted} SAAP`);
      onChanged();
    } catch (e) { showToast(e.message, true); }
  };

  return (
    <React.Fragment>
      <div className="card">
        <h2>Datasets &amp; papers</h2>
        <div className="desc">Each dataset's DOI links its rows to the source paper.</div>
        {datasets.length === 0 ? (
          <div className="empty">No datasets yet. Import a CSV to populate them.</div>
        ) : (
          <React.Fragment>
            <div className="table-wrap">
              <table>
                <thead><tr><th>Dataset</th><th className="num">SAAP</th><th className="num">Obs</th><th>DOI / URL</th><th>Link</th></tr></thead>
                <tbody>
                  {datasets.map((d) => (
                    <tr key={d.name}>
                      <td>{d.name}</td>
                      <td className="num">{d.n_saap}</td>
                      <td className="num">{d.n_observations}</td>
                      <td style={{ minWidth: 320 }}>
                        <input style={{ width: "100%" }} placeholder="DOI or URL"
                               value={edits[d.name] ?? ""} onChange={(e) => setEdits((s) => ({ ...s, [d.name]: e.target.value }))} />
                      </td>
                      <td>{d.url ? <a className="link" href={d.url} target="_blank" rel="noreferrer">open</a> : "—"}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
            <div style={{ marginTop: 14 }}><button onClick={save} disabled={saving}>{saving ? "Saving…" : "Save DOIs"}</button></div>
          </React.Fragment>
        )}
      </div>

      <div className="danger-zone">
        <h2>Danger zone</h2>
        <div className="desc">Remove all imported SAAP and observations. Dataset DOIs are preserved.</div>
        <button className="danger" onClick={wipeAll}>Clear all data</button>
      </div>
    </React.Fragment>
  );
}

/* ---------------------------- Detail drawer ----------------------------- */
function DetailDrawer({ id, datasetUrl, onClose }) {
  const [data, setData] = useState(null);
  useEffect(() => { api.detail(id).then(setData); }, [id]);

  const OBS_COLS = [
    ["dataset", "Dataset"], ["tmt_tissue", "TMT/Tissue"], ["digest", "Digest"],
    ["species", "Species"], ["acquisition_type", "Acquisition"],
    ["saap_pep", "PEP"], ["positional_probability", "PosProb"], ["n_evidence_fragments", "Frags"],
  ];

  return (
    <React.Fragment>
      <div className="drawer-backdrop" onClick={onClose} />
      <div className="drawer">
        <button className="ghost close" onClick={onClose}>Close</button>
        {!data ? <div className="spinner">Loading…</div> : (
          <React.Fragment>
            <h3>{data.saap.mtp_seq}</h3>
            <div className="mono-small">base {data.saap.bp_seq} · {data.saap.aa_sub}</div>
            <div className="kv">
              <div className="k">Genes</div><div className="v">{data.saap.source_gene || "—"}</div>
              <div className="k">RefProteins</div><div className="v">{data.saap.ref_proteins || "—"}</div>
              <div className="k">UniProt</div><div className="v">{data.saap.source_accession || "—"}</div>
              <div className="k">Immunoglobulin</div><div className="v">{String(data.saap.immunoglobulin)}</div>
              <div className="k">Trypsin</div><div className="v">{String(data.saap.trypsin)}</div>
              <div className="k">Missed cleavage</div><div className="v">{String(data.saap.missed_cleavage)}</div>
              <div className="k">AAS at terminus</div><div className="v">{String(data.saap.aas_at_peptide_terminus)}</div>
              <div className="k">Greater than shared</div><div className="v">{String(data.saap.greater_than_shared)}</div>
            </div>
            <h2 style={{ margin: "18px 0 10px" }}>{data.observations.length} observations</h2>
            <div className="table-wrap">
              <table>
                <thead><tr>{OBS_COLS.map(([k, l]) => <th key={k}>{l}</th>)}</tr></thead>
                <tbody>
                  {data.observations.map((o) => (
                    <tr key={o.id}>
                      {OBS_COLS.map(([k]) => (
                        <td key={k} className={typeof o[k] === "number" ? "num" : ""}>
                          {k === "dataset" && datasetUrl[o[k]]
                            ? <a className="link" href={datasetUrl[o[k]]} target="_blank" rel="noreferrer">{o[k]}</a>
                            : o[k] === null || o[k] === undefined ? "—"
                            : typeof o[k] === "number" ? fmt.num(o[k], 3) : o[k]}
                        </td>
                      ))}
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </React.Fragment>
        )}
      </div>
    </React.Fragment>
  );
}

/* ----------------------------- Export modal ----------------------------- */
const DEFAULT_TEMPLATE =
  ">saap|{accession}|{entry_name} {protein} SAAP variant ({aa_sub}) OS={species} GN={gene} BP={bp_seq}";

function ExportModal({ mode, selected, filters, total, onClose, showToast }) {
  const [template, setTemplate] = useState(DEFAULT_TEMPLATE);
  const [busy, setBusy] = useState(false);

  const doExport = async () => {
    setBusy(true);
    try {
      const payload = { header_template: template };
      if (mode === "selected") payload.ids = Array.from(selected);
      else payload.filters = filters;
      const blob = await api.exportFasta(payload);
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url; a.download = `saap_variants_${total}.fasta`;
      document.body.appendChild(a); a.click(); a.remove();
      URL.revokeObjectURL(url);
      showToast(`Exported ${total} SAAP to FASTA`);
      onClose();
    } catch (e) { showToast(e.message, true); }
    finally { setBusy(false); }
  };

  return (
    <div className="modal-backdrop" onClick={onClose}>
      <div className="modal" onClick={(e) => e.stopPropagation()}>
        <h3>Export {total} SAAP to FASTA</h3>
        <div className="desc">Each entry's sequence is the variant (substituted) peptide. Species (OS=) is taken from each SAAP's Species column.</div>
        <div className="field">
          <label>Header template</label>
          <textarea rows={3} value={template} onChange={(e) => setTemplate(e.target.value)} />
          <div className="hint">Fields: {"{accession} {entry_name} {gene} {protein} {aa_sub} {sub_compact} {bp_seq} {mtp_seq} {species} {id}"}</div>
        </div>
        <div className="actions">
          <button className="ghost" onClick={onClose} disabled={busy}>Cancel</button>
          <button onClick={doExport} disabled={busy || total === 0}>{busy ? "Generating…" : "Download .fasta"}</button>
        </div>
      </div>
    </div>
  );
}

ReactDOM.createRoot(document.getElementById("root")).render(<App />);
