const { useState, useEffect, useCallback, useRef } = React;

/* ----------------------------- API helpers ----------------------------- */
const api = {
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
    const r = await fetch("/api/upload", { method: "POST", body: fd });
    const body = await r.json();
    if (!r.ok) throw new Error(body.detail || "Upload failed");
    return body;
  },
  async saveDatasetDois(map) {
    const r = await fetch("/api/datasets", {
      method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ map }),
    });
    if (!r.ok) throw new Error((await r.json().catch(() => ({}))).detail || "Failed to save DOIs");
    return r.json();
  },
  async deleteSaap(payload) {
    const r = await fetch("/api/saap/delete", {
      method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload),
    });
    const body = await r.json();
    if (!r.ok) throw new Error(body.detail || "Delete failed");
    return body;
  },
  async exportFasta(payload, refFile) {
    const fd = new FormData();
    fd.append("payload", JSON.stringify(payload));
    if (refFile) fd.append("reference", refFile);
    const r = await fetch("/api/export/fasta", { method: "POST", body: fd });
    if (!r.ok) {
      const b = await r.json().catch(() => ({}));
      throw new Error(b.detail || "Export failed");
    }
    return r.blob();
  },
  async exportCsv(payload) {
    const r = await fetch("/api/export/csv", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
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
  { key: "mtp_seq", label: "SAAP", sortable: true, cls: "seq", w: 170 },
  { key: "bp_seq", label: "Base peptide", sortable: true, cls: "seq", w: 170 },
  { key: "aa_sub", label: "AAS", sortable: true, w: 90 },
  { key: "source_gene", label: "Genes", sortable: true, w: 120 },
  { key: "ref_proteins", label: "RefProteins", sortable: true, w: 200 },
  { key: "source_accession", label: "UniProt", sortable: true, w: 110 },
  { key: "species", label: "Species", w: 130 },
  { key: "n_datasets", label: "Datasets", sortable: true, w: 170 },
  { key: "digests", label: "Digest", w: 140 },
  { key: "acquisition_types", label: "Acquisition", w: 120 },
  { key: "n_observations", label: "Obs", sortable: true, cls: "num", w: 70 },
  { key: "best_saap_pep", label: "Best PEP", sortable: true, cls: "num", w: 100 },
  { key: "max_positional_probability", label: "Max PosProb", sortable: true, cls: "num", w: 110 },
  { key: "max_evidence_fragments", label: "Max Frags", sortable: true, cls: "num", w: 90 },
  { key: "immunoglobulin", label: "Ig", w: 64 },
  { key: "trypsin", label: "Tryp", w: 70 },
  { key: "missed_cleavage", label: "MissClv", w: 84 },
  { key: "aas_at_peptide_terminus", label: "TermAAS", w: 90 },
  { key: "greater_than_shared", label: ">Shared", w: 90 },
];
const DEFAULT_COL_W = 120;

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
  const [toastNode, showToast] = useToast();

  const refreshMeta = useCallback(() => {
    api.stats().then(setStats).catch(() => {});
    api.facets().then(setFacets).catch(() => {});
    api.datasets().then(setDatasets).catch(() => {});
  }, []);
  useEffect(() => { refreshMeta(); }, [refreshMeta]);

  const datasetUrl = {};
  datasets.forEach((d) => { if (d.url) datasetUrl[d.name] = d.url; });

  return (
    <div className="app">
      <header className="top">
        <div className="brand">
          <a className="logo-link" href="https://github.com/alexmaropakis/SAAPedia" target="_blank" rel="noopener noreferrer" aria-label="SAAPedia on GitHub">
          <svg className="logo" viewBox="0 0 64 64" width="62" height="62" aria-hidden="true">
            {/* mRNA strand (threads through the ribosome; only the ends show) */}
            <path className="logo-mrna-strand" d="M5 39 H59" stroke="var(--logo-mrna)" strokeWidth="3" strokeLinecap="round" fill="none" />
            <circle cx="8" cy="39" r="1.8" fill="var(--logo-mrna)" />
            <circle cx="52" cy="39" r="1.8" fill="var(--logo-mrna)" />
            {/* large + small subunits */}
            <ellipse className="logo-large" cx="28" cy="29" rx="17" ry="14" fill="var(--logo-large)" />
            <ellipse className="logo-small" cx="28" cy="44" rx="13" ry="8.5" fill="var(--logo-small)" />
            {/* growing polypeptide: backbone + amino-acid beads */}
            <path d="M38 19 Q48 8 61 13" stroke="var(--logo-mrna)" strokeWidth="2" strokeLinecap="round" fill="none" />
            <circle className="logo-bead b1" cx="38" cy="19" r="3.5" fill="var(--logo-b1)" />
            <circle className="logo-bead b2" cx="46" cy="12" r="3.5" fill="var(--logo-b2)" />
            <circle className="logo-bead b3" cx="54" cy="11" r="3.5" fill="var(--logo-b3)" />
            <circle className="logo-bead b4" cx="61" cy="13" r="3.5" fill="var(--logo-b4)" />
          </svg>
          </a>
          <h1>SAAPedia</h1>
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
        <ImportTab onIngested={() => { refreshMeta(); }} showToast={showToast} />
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
  const [pageSize, setPageSize] = useState(50);
  const PAGE_SIZES = [25, 50, 100, 200, 500];
  const [selected, setSelected] = useState(() => new Set());
  const [detailId, setDetailId] = useState(null);
  const [showExport, setShowExport] = useState(false);
  const [colWidths, setColWidths] = useState({});

  const colW = (c) => colWidths[c.key] ?? c.w ?? DEFAULT_COL_W;
  const tableWidth = 40 + COLUMNS.reduce((s, c) => s + colW(c), 0);
  const startResize = (e, key) => {
    e.preventDefault();
    e.stopPropagation();
    const startX = e.clientX;
    const col = COLUMNS.find((c) => c.key === key);
    const startW = colWidths[key] ?? (col && col.w) ?? DEFAULT_COL_W;
    const onMove = (ev) => {
      const w = Math.max(50, startW + (ev.clientX - startX));
      setColWidths((prev) => ({ ...prev, [key]: w }));
    };
    const onUp = () => {
      document.removeEventListener("mousemove", onMove);
      document.removeEventListener("mouseup", onUp);
      document.body.style.cursor = "";
      document.body.style.userSelect = "";
    };
    document.addEventListener("mousemove", onMove);
    document.addEventListener("mouseup", onUp);
    document.body.style.cursor = "col-resize";
    document.body.style.userSelect = "none";
  };

  const load = useCallback(() => {
    setLoading(true);
    api.list({ ...filters, sort: sort.key, order: sort.order, page, page_size: pageSize })
      .then(setData).catch((e) => showToast(e.message, true)).finally(() => setLoading(false));
  }, [filters, sort, page, pageSize, showToast]);
  useEffect(() => { load(); }, [load]);

  const changePageSize = (n) => { setPageSize(n); setPage(1); };
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

  const pager = (
    <div className="pager">
      <label htmlFor="pagesize" className="pager-lbl">Rows</label>
      <select id="pagesize" value={pageSize} onChange={(e) => changePageSize(Number(e.target.value))}>
        {PAGE_SIZES.map((n) => <option key={n} value={n}>{n}</option>)}
      </select>
      <button className="ghost small" disabled={page <= 1} onClick={() => setPage((p) => p - 1)}>Prev</button>
      <span>Page {data.page} of {totalPages} · {data.total} SAAP</span>
      <button className="ghost small" disabled={page >= totalPages} onClick={() => setPage((p) => p + 1)}>Next</button>
    </div>
  );

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

      {pager}

      <div className="table-wrap">
        <table className="table-resizable" style={{ width: tableWidth }}>
          <colgroup>
            <col style={{ width: 40 }} />
            {COLUMNS.map((c) => <col key={c.key} style={{ width: colW(c) }} />)}
          </colgroup>
          <thead>
            <tr>
              <th className="checkcol"><input type="checkbox" checked={allOnPage} onChange={togglePage} /></th>
              {COLUMNS.map((c) => (
                <th key={c.key} className={c.sortable ? "sortable" : ""} onClick={() => toggleSort(c)}>
                  {c.label}{sort.key === c.key && <span className="arrow"> {sort.order === "asc" ? "▲" : "▼"}</span>}
                  <span className="col-resizer" onMouseDown={(e) => startResize(e, c.key)} onClick={(e) => e.stopPropagation()} />
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

      {pager}

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
function ImportTab({ onIngested, showToast }) {
  const [drag, setDrag] = useState(false);
  const [busy, setBusy] = useState(false);
  const [result, setResult] = useState(null);
  const inputRef = useRef();

  const handleFile = async (file) => {
    if (!file) return;
    setBusy(true);
    try {
      const res = await api.upload(file);
      setResult(res);
      showToast(`Imported ${res.rows_read} rows → ${res.saap_created} new SAAP`);
      onIngested();
    } catch (e) { showToast(e.message, true); }
    finally { setBusy(false); }
  };

  return (
    <div className="card">
      <h2>Import CSV</h2>
      <div className="desc">Columns are auto-mapped and de-duplicated. Peptides without a UniProt ID are dropped on import. "N datasets" is computed from the datasets each SAAP appears in. Attach source-paper DOIs to datasets later under the Datasets tab.</div>
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
            {result.saap_removed_no_uniprot > 0 && <span className="mini warn">removed (no UniProt ID) <b>{result.saap_removed_no_uniprot}</b></span>}
          </div>
          {result.columns_unmapped && result.columns_unmapped.length > 0 &&
            <div className="badge-row"><span className="mini warn">unmapped columns: {result.columns_unmapped.join(", ")}</span></div>}
        </div>
      )}
    </div>
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
        <div className="desc" style={{ color: "var(--danger)", fontWeight: 700 }}>
          Remove all imported SAAP and observations. Dataset DOIs are preserved.
        </div>
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
  ">sp|{accession}-{mid}-{tok}|{gene}-mut {gene} mistranslated {mid} OS={species} OX={taxid} GN={gene} PE=1 SV=1";

function ExportModal({ mode, selected, filters, total, onClose, showToast }) {
  const [format, setFormat] = useState("fasta");
  const [decoys, setDecoys] = useState(false);
  const [refFile, setRefFile] = useState(null);
  const [template, setTemplate] = useState(DEFAULT_TEMPLATE);
  const [busy, setBusy] = useState(false);

  const download = (blob, name) => {
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url; a.download = name;
    document.body.appendChild(a); a.click(); a.remove();
    URL.revokeObjectURL(url);
  };

  const doExport = async () => {
    setBusy(true);
    try {
      const payload = {};
      if (mode === "selected") payload.ids = Array.from(selected);
      else payload.filters = filters;
      if (format === "csv") {
        download(await api.exportCsv(payload), `saap_export_${total}.csv`);
        showToast(`Exported ${total} SAAP to CSV`);
      } else {
        const blob = await api.exportFasta({ ...payload, decoys, header_template: template }, refFile);
        download(blob, `saap_variants_${total}.fasta`);
        showToast(`Exported ${total} SAAP to FASTA`);
      }
      onClose();
    } catch (e) { showToast(e.message, true); }
    finally { setBusy(false); }
  };

  return (
    <div className="modal-backdrop" onClick={onClose}>
      <div className="modal" onClick={(e) => e.stopPropagation()}>
        <h3>Export {total} SAAP as {format === "csv" ? ".CSV" : "FASTA"}</h3>
        <div className="field">
          <label>Format</label>
          <div className="segmented">
            <button className={format === "fasta" ? "on" : ""} onClick={() => setFormat("fasta")}>FASTA</button>
            <button className={format === "csv" ? "on" : ""} onClick={() => setFormat("csv")}>CSV</button>
          </div>
        </div>
        {format === "fasta" ? (
          <React.Fragment>
            <div className="field">
              <label>Reference proteome FASTA (optional)</label>
              <input type="file" accept=".fasta,.fa,.faa,.txt"
                     onChange={(e) => setRefFile(e.target.files[0] || null)} />
            </div>
            <div className="field">
              <label style={{ display: "flex", alignItems: "center", gap: 8, textTransform: "none", cursor: "pointer" }}>
                <input type="checkbox" checked={decoys} onChange={(e) => setDecoys(e.target.checked)} style={{ width: "auto" }} />
                Append <code>rev_</code> decoys.
              </label>
            </div>
            <div className="field">
              <label>Header template</label>
              <textarea rows={3} value={template} onChange={(e) => setTemplate(e.target.value)} />
            </div>
          </React.Fragment>
        ) : null}
        <div className="actions">
          <button className="ghost" onClick={onClose} disabled={busy}>Cancel</button>
          <button onClick={doExport} disabled={busy || total === 0}>
            {busy ? "Generating…" : format === "csv" ? "Download .csv" : "Download .fasta"}
          </button>
        </div>
      </div>
    </div>
  );
}

ReactDOM.createRoot(document.getElementById("root")).render(<App />);
