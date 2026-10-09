/* Pages. */

const dataChanged = () => window.dispatchEvent(new Event("saap:changed"));
function useDataChanged(fn) {
  useEffect(() => {
    window.addEventListener("saap:changed", fn);
    return () => window.removeEventListener("saap:changed", fn);
  }, [fn]);
}
function useSession(key, initial) {
  const [v, setV] = useState(() => {
    try { const s = sessionStorage.getItem(key); return s ? JSON.parse(s) : initial; } catch { return initial; }
  });
  useEffect(() => { try { sessionStorage.setItem(key, JSON.stringify(v)); } catch { /* unavailable */ } }, [key, v]);
  return [v, setV];
}
const shortSpecies = (s) => s.replace(/^(\w)\w*\s+(\w+)/, "$1. $2");
const chips = (arr, map = (x) => x, max = Infinity) => {
  if (!arr || !arr.length) return DASH;
  const shown = arr.slice(0, max).map((v) => <span key={v} className="chip">{map(v)}</span>);
  return arr.length > max ? <span title={arr.join("\n")}>{shown}<span className="muted">+{arr.length - max}</span></span> : shown;
};
const fmtAF = (v) => (v == null ? DASH : v === 0 ? "0" : v.toExponential(1));
const variantLabel = (s) => (s.ref && s.alt && s.positions && s.positions.length ? `${s.ref}${s.positions[0]}${s.alt}` : null);

function annotateSummary(r) {
  const bits = [`${r.positioned.toLocaleString()} of ${r.requested.toLocaleString()} positioned`];
  [["resolved_by_gene", "via gene"], ["resolved_by_peptide", "via known peptide"], ["resolved_by_reference", "via reference proteome"], ["resolved_by_sequence", "via peptide search"],
   ["species_corrected", "species-corrected"], ["aas_filled", "AAS filled"], ["merged_duplicates", "duplicates merged"],
   ["not_found", "not found"], ["unmatched_peptide", "peptide not in protein"], ["failed", "failed — retry"]]
    .forEach(([k, label]) => r[k] && bits.push(`${r[k].toLocaleString()} ${label}`));
  return bits.join(" · ");
}

/* ================================ Home ================================ */
const ACCESSION = /^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-\d+)?$/;

/** Split-flap board: tiles flip through random residues, settle on a real base
 *  peptide, then the substituted tile flips again to the SAAP residue. */
const RESIDUES = "ACDEFGHIKLMNPQRSTVWY";
const FALLBACK_PAIRS = [{ bp: "SAVTALWGK", saap: "SAVTALDGK" }];

function FlapBoard() {
  const [pairs, setPairs] = useState(FALLBACK_PAIRS);
  const [frame, setFrame] = useState({ text: FALLBACK_PAIRS[0].saap, diff: 6, marked: true, settled: [] });
  useEffect(() => {
    api.list({ ...DEFAULT_FILTERS, sort: "n_observations", order: "desc", page_size: 60 }).then((d) => {
      const ok = d.items.filter((r) => r.bp_seq && r.mtp_seq && r.bp_seq.length === r.mtp_seq.length && r.bp_seq.length <= 14
        && [...r.bp_seq].filter((c, i) => c !== r.mtp_seq[i]).length === 1);
      if (ok.length) setPairs(ok.map((r) => ({ bp: r.bp_seq, saap: r.mtp_seq })));
    }).catch(() => {});
  }, []);
  useEffect(() => {
    if (window.matchMedia("(prefers-reduced-motion: reduce)").matches) {
      const { bp, saap } = pairs[0];
      setFrame({ text: saap, diff: [...bp].findIndex((c, i) => c !== saap[i]), marked: true, settled: [] });
      return;
    }
    let n = 0, start = performance.now();
    const rand = () => RESIDUES[Math.floor(Math.random() * RESIDUES.length)];
    const tick = setInterval(() => {
      const { bp, saap } = pairs[n % pairs.length];
      const diff = [...bp].findIndex((c, i) => c !== saap[i]);
      const t = performance.now() - start;
      const settleAt = (i) => 250 + i * 110;          // phase 1: tiles land left to right
      const flipStart = settleAt(bp.length) + 1100;   // phase 2: the substituted tile flips again
      const flipEnd = flipStart + 650;
      const done = flipEnd + 2400;
      if (t > done) { n += 1; start = performance.now(); return; }
      const text = [...bp].map((c, i) => {
        if (t < settleAt(i)) return rand();
        if (i === diff && t >= flipStart) return t < flipEnd ? rand() : saap[i];
        return c;
      }).join("");
      setFrame({ text, diff, marked: t >= flipEnd, settled: [...bp].map((_, i) => t >= settleAt(i) && !(i === diff && t >= flipStart && t < flipEnd)) });
    }, 55);
    return () => clearInterval(tick);
  }, [pairs]);
  return (
    <div className="flap" aria-label="Substituted peptide">
      {[...frame.text].map((c, i) => (
        <span key={i} className={"tile" + (i === frame.diff && frame.marked ? " sub" : "") + (frame.settled[i] === false ? " spinning" : "")}>{c}</span>
      ))}
    </div>
  );
}

function HomePage({ stats }) {
  const [q, setQ] = useState("");
  const [proteins] = useAsync(() => api.proteins({ sort: "n_saap", order: "desc", page_size: 8 }), []);
  const [overview] = useAsync(api.datasets, []);
  const search = (e) => {
    e.preventDefault();
    const term = q.trim();
    if (!term) return go("browse");
    if (ACCESSION.test(term.toUpperCase())) return go(`protein/${term.toUpperCase()}`);
    try {
      sessionStorage.setItem("browse.filters", JSON.stringify({ ...DEFAULT_FILTERS, q: term }));
      sessionStorage.setItem("browse.page", "1");
    } catch { /* storage unavailable */ }
    go("browse");
  };
  return (
    <div className="page stack home">
      <div className="hero">
        <FlapBoard />
        <h1><Wordmark /></h1>
        <p className="muted">Functional database of amino acid substitutions arising from alternate RNA decoding identified via mass spectrometry-based proteomics</p>
        <form className="search hero-search" onSubmit={search}>
          <Icon name="search" />
          <input autoFocus type="search" placeholder="Peptide, gene, protein or UniProt accession" value={q} onChange={(e) => setQ(e.target.value)} />
          <button className="primary" type="submit">Search</button>
        </form>
      </div>
      {stats && (
        <div className="metrics">
          <Metric k="SAAP" v={stats.n_saap.toLocaleString()} />
          <Metric k="Proteins" v={stats.n_proteins.toLocaleString()} />
          <Metric k="Observations" v={stats.n_observations.toLocaleString()} />
          <Metric k="Tissues & cell types" v={stats.n_tissues.toLocaleString()} />
        </div>
      )}
      <div className="grid-2">
        <div className="card">
          <div className="card-head"><h2>Most substituted proteins</h2><span className="spacer" /><a href="#/proteins">All proteins</a></div>
          {!proteins ? <Loading /> : (
            <table><tbody>{proteins.items.map((p) => (
              <tr key={p.protein_accession} className="link" onClick={() => go(`protein/${p.protein_accession}`)}>
                <td><b>{fmt.text(p.gene)}</b></td>
                <td className="trunc muted" title={p.description}>{fmt.text(p.description)}</td>
                <td className="num">{plural(p.n_saap, "SAAP")}</td>
                <td><MiniMap length={p.length} positions={p.positions} width={120} /></td>
              </tr>
            ))}</tbody></table>
          )}
        </div>
        <div className="card">
          <div className="card-head"><h2>Top substitutions</h2><span className="spacer" /><a href="#/tissues">Tissues</a></div>
          <div className="card-body">{overview ? <BarList items={overview.top_substitutions.slice(0, 8)} /> : <Loading />}</div>
        </div>
      </div>
    </div>
  );
}

/* =============================== Browse =============================== */
const COLUMNS = [
  { key: "mtp_seq", label: "SAAP", sortable: true, cls: "seq", fixed: true, render: (r) => <SubSeq saap={r.mtp_seq} bp={r.bp_seq} /> },
  { key: "bp_seq", label: "Base peptide", sortable: true, cls: "seq" },
  { key: "aa_sub", label: "AAS", sortable: true },
  { key: "source_gene", label: "Gene", sortable: true, render: (r) => fmt.text(firstOf(r.source_gene)) },
  { key: "ref_proteins", label: "Protein", sortable: true, cls: "trunc",
    render: (r) => { const p = firstOf(r.protein_description || r.ref_proteins); return p ? <span title={p}>{p}</span> : DASH; } },
  { key: "source_accession", label: "UniProt", sortable: true, hidden: true },
  { key: "ensembl_gene", label: "Ensembl gene", sortable: true, hidden: true },
  { key: "position_in_protein", label: "Position", sortable: true, num: true,
    render: (r) => (r.positions_all ? (r.n_positions > 1 ? <span title={r.positions_all}>{r.positions_all.split(",")[0]} +{r.n_positions - 1}</span> : r.positions_all) : fmt.text(r.position_in_protein)) },
  { key: "species", label: "Species", render: (r) => chips(r.species, shortSpecies) },
  { key: "n_tissues", label: "Tissues", sortable: true, title: "Tissues / cell types (per species)", render: (r) => chips(r.tissues, undefined, 2) },
  { key: "datasets", label: "Datasets", hidden: true, private: true, render: (r) => chips(r.datasets, undefined, 2) },
  { key: "digests", label: "Digest", hidden: true, render: (r) => chips(r.digests) },
  { key: "acquisition_types", label: "Acquisition", hidden: true, render: (r) => chips(r.acquisition_types) },
  { key: "n_observations", label: "Obs", sortable: true, num: true, title: "Observations" },
  { key: "best_saap_pep", label: "Best PEP", sortable: true, num: true, render: (r) => fmt.sci(r.best_saap_pep) },
  { key: "max_positional_probability", label: "PosProb", sortable: true, num: true, title: "Best positional probability (observations and source data)",
    render: (r) => (r.max_positional_probability == null ? DASH : <span className={r.max_positional_probability < 0.9 ? "no" : ""}>{r.max_positional_probability.toFixed(3)}</span>) },
  { key: "max_evidence_fragments", label: "Fragments", sortable: true, num: true, hidden: true, title: "Max evidence fragments" },
  { key: "trypsin", label: "Trypsin", hidden: true, render: (r) => fmt.bool(r.trypsin) },
  { key: "missed_cleavage", label: "Missed cleavage", hidden: true, render: (r) => fmt.bool(r.missed_cleavage) },
  { key: "aas_at_peptide_terminus", label: "Terminal AAS", hidden: true, render: (r) => fmt.bool(r.aas_at_peptide_terminus) },
  { key: "greater_than_shared", label: "> Shared", hidden: true, render: (r) => fmt.bool(r.greater_than_shared) },
  { key: "at_cleavage_site", label: "Cleavage site", title: "Computed: the reference or substituted residue is a cut site for the observed digest",
    render: (r) => fmt.bool(r.at_cleavage_site) },
  { key: "gnomad_af", label: "gnomAD AF", sortable: true, num: true,
    title: "gnomAD v4 allele frequency of germline variants producing this exact substitution (human)",
    render: (r) => (r.gnomad_status === "present" ? fmtAF(r.gnomad_af) : r.gnomad_status === "absent" ? <span className="no">Absent</span> : DASH) },
];

const FILTERS = [
  { key: "tissue", label: "Tissue / cell type", facet: "tissues" },
  { key: "dataset", label: "Dataset", facet: "datasets", private: true },
  { key: "digest", label: "Digest", facet: "digests" },
  { key: "species", label: "Species", facet: "species" },
  { key: "acquisition_type", label: "Acquisition", facet: "acquisition_types" },
  { key: "aa_sub", label: "AAS", facet: "aa_subs" },
  { key: "at_cleavage_site", label: "Cleavage site", bool: true },
  { key: "in_gnomad", label: "In gnomAD", bool: true },
  { key: "min_pos_prob", label: "Min PosProb", number: true, step: 0.01 },
  { key: "max_pep", label: "Max PEP", number: true, step: "any" },
];
const EMPTY_FILTERS = Object.fromEntries([["q", ""], ...FILTERS.map((f) => [f.key, ""])]);
// Noisiest artifact classes hidden by default; shown as removable chips.
const DEFAULT_FILTERS = { ...EMPTY_FILTERS, at_cleavage_site: "false" };
const selection = { ids: new Set() };  // survives navigating into a SAAP and back

function BrowsePage({ facets }) {
  const toast = useToast();
  const [stored, setFilters] = useSession("browse.filters", DEFAULT_FILTERS);
  // Drop keys from older sessions that are no longer filters.
  const filters = useMemo(() => Object.fromEntries(Object.keys(EMPTY_FILTERS).map((k) => [k, stored[k] ?? ""])), [stored]);
  const [sort, setSort] = useSession("browse.sort", { key: "n_observations", order: "desc" });
  const [page, setPage] = useSession("browse.page", 1);
  const [pageSize, setPageSize] = useStored("browse.pageSize", 50);
  const [hidden, setHidden] = useStored("browse.hiddenColumns", COLUMNS.filter((c) => c.hidden).map((c) => c.key));
  const [showFilters, setShowFilters] = useState(false);
  const [selected, setSelected] = useState(selection.ids);
  const [exporting, setExporting] = useState(null);
  const [busy, setBusy] = useState(false);
  const [q, setQ] = useState(filters.q);

  useEffect(() => { selection.ids = selected; }, [selected]);
  useEffect(() => { const t = setTimeout(() => q !== filters.q && setFilter("q", q), 250); return () => clearTimeout(t); }, [q]);

  const [data, error, reload] = useAsync(
    () => api.list({ ...filters, sort: sort.key, order: sort.order, page, page_size: pageSize }),
    [JSON.stringify(filters), sort.key, sort.order, page, pageSize], true);
  useDataChanged(reload);

  const setFilter = (k, v) => { setFilters((f) => ({ ...f, [k]: v })); setPage(1); };
  const replaceFilters = (f) => { setFilters(f); setQ(f.q); setPage(1); };
  const priv = Boolean(facets.datasets);  // dataset names are only served in private mode
  const visibleColumns = COLUMNS.filter((c) => priv || !c.private);
  const visibleFilters = FILTERS.filter((f) => priv || !f.private);
  const cols = visibleColumns.filter((c) => c.fixed || !hidden.includes(c.key));
  const items = data ? data.items : [];
  const allOnPage = items.length > 0 && items.every((r) => selected.has(r.id));
  const toggle = (id) => setSelected((s) => { const n = new Set(s); n.has(id) ? n.delete(id) : n.add(id); return n; });
  const togglePage = () => setSelected((s) => { const n = new Set(s); items.forEach((r) => (allOnPage ? n.delete(r.id) : n.add(r.id))); return n; });
  const active = visibleFilters.filter((f) => filters[f.key] !== "");
  const isDefault = JSON.stringify({ ...filters, q: "" }) === JSON.stringify({ ...DEFAULT_FILTERS, q: "" });

  const annotateSelected = async () => {
    setBusy(true);
    try { toast(annotateSummary(await api.annotate({ ids: [...selected] }))); dataChanged(); }
    catch (e) { toast(e.message, true); } finally { setBusy(false); }
  };
  const deleteSelected = async () => {
    if (!confirm(`Delete ${plural(selected.size, "SAAP")} and their observations? This cannot be undone.`)) return;
    try { const r = await api.remove({ ids: [...selected] }); toast(`Deleted ${plural(r.deleted, "SAAP")}`); setSelected(new Set()); dataChanged(); }
    catch (e) { toast(e.message, true); }
  };

  const filterValue = (f) => (f.bool ? (filters[f.key] === "true" ? "Yes" : "No") : filters[f.key]);

  return (
    <div className="stack">
      <div className="toolbar">
        <div className="search"><Icon name="search" />
          <input type="search" placeholder="Search peptide, gene, protein or ID" value={q} onChange={(e) => setQ(e.target.value)} /></div>
        <button className={showFilters ? "on" : ""} onClick={() => setShowFilters((s) => !s)}>
          <Icon name="filter" />Filters{active.length > 0 && <span className="badge">{active.length}</span>}</button>
        <span className="spacer" />
        {selected.size > 0 ? (
          <div className="selbar">
            <b>{selected.size.toLocaleString()}</b>&nbsp;selected
            <button onClick={() => setExporting({ ids: [...selected], count: selected.size })}><Icon name="download" />Export</button>
            {priv && <button onClick={annotateSelected} disabled={busy}>{busy ? <Spinner /> : <Icon name="spark" />}Annotate</button>}
            {priv && <button onClick={deleteSelected}>Delete</button>}
            <button className="icon" onClick={() => setSelected(new Set())} title="Clear selection"><Icon name="x" /></button>
          </div>
        ) : (
          <Fragment>
            <Menu label={<Fragment><Icon name="columns" />Columns</Fragment>}>
              {visibleColumns.filter((c) => !c.fixed).map((c) => (
                <label key={c.key}><input type="checkbox" checked={!hidden.includes(c.key)}
                  onChange={() => setHidden((h) => (h.includes(c.key) ? h.filter((k) => k !== c.key) : [...h, c.key]))} />{c.label}</label>
              ))}
            </Menu>
            <button className="primary" disabled={!data || !data.total}
                    onClick={() => setExporting({ filters, count: data.total })}><Icon name="download" />Export {data ? data.total.toLocaleString() : ""}</button>
          </Fragment>
        )}
      </div>

      <div className="card">
        {showFilters && (
          <div className="filter-panel">
            {visibleFilters.map((f) => (
              <div className="field" key={f.key}>
                <label>{f.label}</label>
                {f.number ? (
                  <input type="number" step={f.step} min="0" value={filters[f.key]} onChange={(e) => setFilter(f.key, e.target.value)} />
                ) : (
                  <select value={filters[f.key]} onChange={(e) => setFilter(f.key, e.target.value)}>
                    <option value="">Any</option>
                    {f.bool ? <Fragment><option value="true">Yes</option><option value="false">No</option></Fragment>
                      : (facets[f.facet] || []).map((o) => <option key={o} value={o}>{o}</option>)}
                  </select>
                )}
              </div>
            ))}
          </div>
        )}
        {(active.length > 0 || !isDefault) && (
          <div className="chips-row">
            {active.map((f) => (
              <span key={f.key} className="chip filter"><b>{f.label}</b>{filterValue(f)}
                <button onClick={() => setFilter(f.key, "")}><Icon name="x" /></button></span>
            ))}
            {active.length > 0 && <button className="ghost sm" onClick={() => replaceFilters({ ...EMPTY_FILTERS, q: filters.q })}>Clear all</button>}
            {!isDefault && <button className="ghost sm" onClick={() => replaceFilters({ ...DEFAULT_FILTERS, q: filters.q })}>Reset to defaults</button>}
          </div>
        )}
        <div className="table-wrap">
          <table>
            <thead><tr>
              <th className="check"><input type="checkbox" checked={allOnPage} onChange={togglePage} /></th>
              {cols.map((c) => <Th key={c.key} col={c} sort={sort} onSort={(k) => { setSort(toggleSort(k)); setPage(1); }} />)}
            </tr></thead>
            <tbody>
              {error ? <tr><td colSpan={cols.length + 1}><ErrorNote error={error} /></td></tr>
                : !data ? <tr><td colSpan={cols.length + 1}><Loading /></td></tr>
                : !items.length ? <tr><td colSpan={cols.length + 1} className="empty">No SAAP match these filters</td></tr>
                : items.map((r) => (
                  <tr key={r.id} className={"link" + (selected.has(r.id) ? " selected" : "")} onClick={() => go(`saap/${r.id}`)}>
                    <td className="check" onClick={(e) => e.stopPropagation()}><input type="checkbox" checked={selected.has(r.id)} onChange={() => toggle(r.id)} /></td>
                    {cols.map((c) => <td key={c.key} className={(c.cls || "") + (c.num ? " num" : "")}>{c.render ? c.render(r) : fmt.text(r[c.key])}</td>)}
                  </tr>
                ))}
            </tbody>
          </table>
        </div>
        {data && <Pager page={page} pageSize={pageSize} total={data.total} onPage={setPage} onPageSize={(n) => { setPageSize(n); setPage(1); }} />}
      </div>
      {exporting && <ExportModal scope={exporting} onClose={() => setExporting(null)} />}
    </div>
  );
}

/* ============================== Proteins ============================== */
const PROTEIN_COLUMNS = [
  { key: "gene", label: "Gene", sortable: true },
  { key: "protein_accession", label: "Accession", sortable: true },
  { key: "description", label: "Protein" },
  { key: "length", label: "Length", sortable: true, num: true },
  { key: "n_saap", label: "SAAP", sortable: true, num: true },
  { key: "n_sites", label: "Sites", sortable: true, num: true },
  { key: "n_observations", label: "Obs", sortable: true, num: true, title: "Observations" },
  { key: "map", label: "Site map" },
];

function ProteinsPage() {
  const [q, setQ] = useState("");
  const [query, setQuery] = useState("");
  const [sort, setSort] = useSession("proteins.sort", { key: "n_saap", order: "desc" });
  const [page, setPage] = useSession("proteins.page", 1);
  const [pageSize, setPageSize] = useStored("proteins.pageSize", 50);
  useEffect(() => { const t = setTimeout(() => { setQuery(q); setPage(1); }, 250); return () => clearTimeout(t); }, [q]);
  const [data, error] = useAsync(() => api.proteins({ q: query, sort: sort.key, order: sort.order, page, page_size: pageSize }),
    [query, sort.key, sort.order, page, pageSize], true);

  return (
    <div className="stack">
      <div className="toolbar">
        <div className="search"><Icon name="search" />
          <input type="search" placeholder="Search gene, accession or protein" value={q} onChange={(e) => setQ(e.target.value)} /></div>
      </div>
      <div className="card">
        <div className="table-wrap">
          <table>
            <thead><tr>{PROTEIN_COLUMNS.map((c) => <Th key={c.key} col={c} sort={sort} onSort={(k) => { setSort(toggleSort(k)); setPage(1); }} />)}</tr></thead>
            <tbody>
              {error ? <tr><td colSpan={8}><ErrorNote error={error} /></td></tr>
                : !data ? <tr><td colSpan={8}><Loading /></td></tr>
                : !data.items.length ? <tr><td colSpan={8} className="empty">No annotated proteins</td></tr>
                : data.items.map((p) => (
                  <tr key={p.protein_accession} className="link" onClick={() => go(`protein/${p.protein_accession}`)}>
                    <td><b>{fmt.text(p.gene)}</b></td>
                    <td className="mono">{p.protein_accession}</td>
                    <td className="trunc" title={p.description}>{fmt.text(p.description)}</td>
                    <td className="num">{fmt.int(p.length)}</td>
                    <td className="num">{fmt.int(p.n_saap)}</td>
                    <td className="num">{fmt.int(p.n_sites)}</td>
                    <td className="num">{fmt.int(p.n_observations)}</td>
                    <td><MiniMap length={p.length} positions={p.positions} /></td>
                  </tr>
                ))}
            </tbody>
          </table>
        </div>
        {data && <Pager page={page} pageSize={pageSize} total={data.total} onPage={setPage} onPageSize={(n) => { setPageSize(n); setPage(1); }} />}
      </div>
    </div>
  );
}

/* ====================== Protein + SAAP shared bits ===================== */
const ExtLink = ({ href, children }) => <a href={href} target="_blank" rel="noopener noreferrer">{children}<Icon name="external" /></a>;

function Metric({ k, v, d, dClass, mono, title, children }) {
  return (
    <div className="metric" title={title}>
      <div className="k">{k}</div>
      <div className={"v" + (mono ? " mono" : "")}>{v ?? DASH}</div>
      {d && <div className={"d " + (dClass || "")}>{d}</div>}
      {children}
    </div>
  );
}

function SiteScores({ s, site }) {
  return (
    <Fragment>
      <td className="num">{fmt.text(s.blosum62)}</td>
      <td className="num">{fmt.text(s.grantham)}</td>
      <td className="num">{site && site.plddt != null ? <span><i className="dot" style={{ background: plddtHex(site.plddt) }} />{site.plddt.toFixed(1)}</span> : DASH}</td>
      <td className="num">{site && site.am_pathogenicity != null ? <span title={AM_CLASS[site.am_class]}><i className="dot" style={{ background: amColor(site.am_pathogenicity) }} />{site.am_pathogenicity.toFixed(2)}</span> : DASH}</td>
      <td>{site && site.known_variants.length ? site.known_variants.map((v, i) => (
        <Fragment key={i}>{i > 0 && ", "}{v.alt.split(",").includes(s.alt) ? <b className="sub">{v.ref}→{v.alt}</b> : `${v.ref}→${v.alt}`}</Fragment>)) : DASH}</td>
    </Fragment>
  );
}

function ObservationsTable({ observations }) {
  const priv = observations.length > 0 && "dataset" in observations[0];  // private mode only
  const cols = [["tissue", "Tissue / cell type"], ["species", "Species"], ["digest", "Digest"],
    ["acquisition_type", "Acquisition"], ["saap_pep", "PEP"], ["positional_probability", "PosProb"], ["n_evidence_fragments", "Fragments"],
    ...(priv ? [["dataset", "Dataset"], ["tmt_tissue", "TMT / tissue"], ["source_file", "File"]] : [])];
  const numeric = { saap_pep: fmt.sci, positional_probability: (v) => fmt.num(v, 3), n_evidence_fragments: fmt.int };
  return (
    <div className="table-wrap">
      <table>
        <thead><tr>{cols.map(([k, l]) => <th key={k} className={numeric[k] ? "num" : ""}>{l}</th>)}</tr></thead>
        <tbody>{observations.map((o) => (
          <tr key={o.id}>{cols.map(([k]) => <td key={k} className={numeric[k] ? "num" : ""}>{numeric[k] ? numeric[k](o[k]) : fmt.text(o[k])}</td>)}</tr>
        ))}</tbody>
      </table>
    </div>
  );
}

/* =============================== Protein ============================== */
function FeatureTable({ features, sites }) {
  if (!features) return <Loading />;
  const rows = features.filter((f) => f.track !== "variant");
  if (!rows.length) return <div className="empty">No UniProt features</div>;
  const positions = sites.flatMap((s) => s.positions);
  return (
    <div className="table-wrap">
      <table>
        <thead><tr><th>Type</th><th>Description</th><th className="num">Start</th><th className="num">End</th><th className="num" title="SAAP sites inside this feature">SAAP sites</th></tr></thead>
        <tbody>{rows.map((f, i) => {
          const n = new Set(positions.filter((p) => f.start <= p && p <= f.end)).size;
          return <tr key={i}><td>{f.type}</td><td className="trunc" title={f.description}>{fmt.text(f.description)}</td>
            <td className="num">{f.start}</td><td className="num">{f.end}</td><td className="num">{n || DASH}</td></tr>;
        })}</tbody>
      </table>
    </div>
  );
}

function ProteinPage({ accession, initialTab }) {
  const [view, error] = useAsync(() => api.protein(accession), [accession]);
  const [ann] = useAsync(() => api.proteinAnnotations(accession), [accession]);
  const [tab, setTab] = useState(initialTab || "sites");
  const [pos, setPos] = useState(null);

  if (error) return <div className="page"><ErrorNote error={error} /></div>;
  if (!view) return <Loading />;
  const { protein: p, sites } = view;
  const af = ann && ann.alphafold;
  const shown = pos ? sites.filter((s) => s.positions.includes(pos)) : sites;
  const nObs = sites.reduce((a, s) => a + (s.n_observations || 0), 0);
  const nSites = new Set(sites.flatMap((s) => s.positions)).size;

  return (
    <div className="page stack">
      <div className="crumbs"><a href="#/proteins">Proteins</a><span>/</span><span>{p.gene || accession}</span></div>
      <div>
        <div className="title-row"><h1>{p.gene || accession}</h1><span className="muted" style={{ fontSize: 15 }}>{p.description}</span></div>
        <div className="subtitle" style={{ marginTop: 6 }}>
          <ExtLink href={`https://www.uniprot.org/uniprotkb/${accession}`}>UniProt {accession}</ExtLink>
          {af && af.model_id && <ExtLink href={`https://alphafold.ebi.ac.uk/entry/${accession}`}>AlphaFold</ExtLink>}
          {p.ensembl_gene && <ExtLink href={`https://www.ensembl.org/id/${p.ensembl_gene}`}>{p.ensembl_gene}</ExtLink>}
          {ann && ann.uniprot.organism && <i>{ann.uniprot.organism}</i>}
        </div>
      </div>
      <div className="metrics">
        <Metric k="Length" v={`${p.length.toLocaleString()} aa`} />
        <Metric k="SAAP" v={sites.length.toLocaleString()} />
        <Metric k="Sites" v={nSites.toLocaleString()} d={`${(100 * nSites / p.length).toFixed(1)}% of residues`} />
        <Metric k="Observations" v={nObs.toLocaleString()} />
        <Metric k="AlphaFold pLDDT" v={af && af.global_plddt != null ? af.global_plddt.toFixed(1) : ann ? null : <Spinner />}
                d={af && af.global_plddt != null ? plddtBand(af.global_plddt)[1] : af && af.issue} />
      </div>
      <div className="card">
        <ProteinMap length={p.length} sequence={p.sequence} sites={sites} annotations={ann}
                    onSelect={(g) => { setPos(g.pos); setTab("sites"); }} />
      </div>
      <div className="card"><SequenceView sequence={p.sequence} sites={sites} plddt={af && af.plddt} /></div>
      <div className="card">
        <Tabs value={tab} onChange={setTab} tabs={[["sites", "SAAP", shown.length], ["features", "Features", ann ? ann.uniprot.features.filter((f) => f.track !== "variant").length : null], ...(af && af.has_structure ? [["structure", "Structure"]] : [])]} />
        {tab === "sites" && (
          <Fragment>
            {pos && <div className="chips-row"><span className="chip filter"><b>Position</b>{pos}<button onClick={() => setPos(null)}><Icon name="x" /></button></span></div>}
            <div className="table-wrap">
              <table>
                <thead><tr>
                  <th>Variant</th><th>SAAP</th><th>AAS</th>
                  <th className="num" title="BLOSUM62 substitution score">BLOSUM62</th>
                  <th className="num" title="Grantham chemical distance">Grantham</th>
                  <th className="num" title={PLDDT_TITLE}>pLDDT</th>
                  <th className="num" title="AlphaMissense pathogenicity">AlphaMissense pathogenicity</th>
                  <th title="UniProt natural variants at this position; a match to this substitution suggests a genetic variant">UniProt variant</th>
                  <th className="num">Obs</th><th className="num">Best PEP</th><th>Tissues</th>
                </tr></thead>
                <tbody>{shown.map((s) => (
                  <tr key={s.id} className="link" onClick={() => go(`saap/${s.id}`)}>
                    <td className="mono">{fmt.text(variantLabel(s))}</td>
                    <td className="seq"><SubSeq saap={s.mtp_seq} bp={s.bp_seq} /></td>
                    <td>{s.aa_sub}</td>
                    <SiteScores s={s} site={ann && ann.sites[s.id]} />
                    <td className="num">{fmt.int(s.n_observations)}</td>
                    <td className="num">{fmt.sci(s.best_saap_pep)}</td>
                    <td>{chips(s.tissues, undefined, 2)}</td>
                  </tr>
                ))}</tbody>
              </table>
            </div>
          </Fragment>
        )}
        {tab === "features" && <FeatureTable features={ann ? ann.uniprot.features : null} sites={sites} />}
        {tab === "structure" && af && af.has_structure && <StructureView accession={accession} sites={sites} />}
      </div>
    </div>
  );
}

/* ================================ SAAP ================================ */
const PLDDT_TITLE = "pLDDT: AlphaFold's per-residue confidence in its predicted structure (0–100)";

function KnownVariantMetric({ s, site }) {
  if (!site) return <Metric k="UniProt variant" v={null} />;
  const vs = site.known_variants;
  const same = vs.filter((v) => v.alt.split(",").includes(s.alt));
  const title = vs.map((v) => `${v.id} ${v.ref}→${v.alt} ${v.description}`).join("\n");
  if (same.length) return <Metric k="UniProt variant" v={same.map((v) => `${v.ref}→${v.alt}`).join(", ")} d="Known genetic variant" dClass="bad" title={title} />;
  return <Metric k="UniProt variant" v="None" d={vs.length ? `Other substitutions known here: ${vs.map((v) => v.alt).join(", ")}` : "Not a known variant"} dClass="ok" title={title} />;
}

function GnomadMetric({ s, human }) {
  const title = "gnomAD v4 (exomes + genomes): germline variants producing this exact substitution";
  if (!s.gnomad_status) return <Metric k="gnomAD" v={null} d={human ? "Not checked yet" : "Human SAAPs only"} title={title} />;
  if (s.gnomad_status === "unmapped") return <Metric k="gnomAD" v="Not mapped" d="No Ensembl transcript identical to this protein" title={title} />;
  if (s.gnomad_status === "absent") return <Metric k="gnomAD" v="Absent" d="Not a known polymorphism" dClass="ok" title={title} />;
  const ids = (s.gnomad_variants || "").split(",").filter(Boolean);
  return (
    <Metric k="gnomAD" v={`AF ${fmtAF(s.gnomad_af)}`} dClass="bad" title={title + "\n" + ids.join("\n")}
            d={<Fragment>Rare polymorphism · {ids.slice(0, 2).map((id, i) => (
              <Fragment key={id}>{i > 0 && ", "}<a href={`https://gnomad.broadinstitute.org/variant/${id}?dataset=gnomad_r4`} target="_blank" rel="noopener noreferrer">{id}</a></Fragment>))}</Fragment>} />
  );
}

/** UniProt features overlapping the substitution site. */
function SiteContext({ pos, features }) {
  const here = features.filter((f) => f.track !== "variant" && f.start <= pos && pos <= f.end)
    .sort((a, b) => (a.end - a.start) - (b.end - b.start));
  return (
    <div className="card">
      <div className="card-head"><h2>Site context</h2><span className="muted">UniProt features at position {pos}</span></div>
      {here.length ? (
        <table><tbody>{here.map((f, i) => (
          <tr key={i}><td className="muted" style={{ width: 180 }}>{f.type}</td><td>{f.description || DASH}</td>
            <td className="num mono">{f.start === f.end ? f.start : `${f.start}–${f.end}`}</td></tr>
        ))}</tbody></table>
      ) : <div className="notice">No annotated feature covers this position</div>}
    </div>
  );
}

function SaapPage({ id, initialTab }) {
  const [detail, error] = useAsync(() => api.saap(id), [id]);
  const acc = detail && detail.saap.protein_accession;
  const [view] = useAsync(() => (acc ? api.protein(acc) : Promise.resolve(null)), [acc]);
  const [ann] = useAsync(() => (acc ? api.proteinAnnotations(acc) : Promise.resolve(null)), [acc]);
  const [tab, setTab] = useState(initialTab || "observations");

  if (error) return <div className="page"><ErrorNote error={error} /></div>;
  if (!detail || detail.saap.id !== Number(id)) return <Loading />;
  const s = detail.saap;
  const site = ann && ann.sites[s.id];
  const af = ann && ann.alphafold;
  const gene = firstOf(s.source_gene);
  const variant = variantLabel(s);
  const best = detail.observations.reduce((m, o) => (o.saap_pep != null && (m == null || o.saap_pep < m) ? o.saap_pep : m), null);
  const species = [...new Set(detail.observations.map((o) => o.species).filter(Boolean))];
  const protRes = view && s.positions.length ? view.protein.sequence[s.positions[0] - 1] : null;
  const refMismatch = protRes && s.ref && protRes !== s.ref && !(protRes + s.ref).split("").every((c) => "IL".includes(c));

  return (
    <div className="page stack">
      <div className="crumbs">
        <a href="#/browse">Browse</a><span>/</span>
        {acc ? <a href={`#/protein/${acc}`}>{gene || acc}</a> : <span>{gene || "Unannotated"}</span>}
        <span>/</span><span>{variant || s.aa_sub}</span>
      </div>
      <div>
        <div className="title-row"><h1 className="seq"><SubSeq saap={s.mtp_seq} bp={s.bp_seq} /></h1></div>
        <div className="subtitle" style={{ marginTop: 6 }}>
          <span className="chip">{s.aa_sub || "—"}</span>
          {gene && <b>{gene}</b>}
          {s.protein_description && <span>{s.protein_description}</span>}
          {acc && <ExtLink href={`https://www.uniprot.org/uniprotkb/${acc}`}>{acc}</ExtLink>}
          {species.map((x) => <i key={x}>{x}</i>)}
          <span>{plural(detail.observations.length, "observation")}</span>
          {best != null && <span>Best PEP {best === 0 ? 0 : best.toExponential(1)}</span>}
          {s.max_positional_probability != null && <span className={s.max_positional_probability < 0.9 ? "sub" : ""}>PosProb {s.max_positional_probability.toFixed(2)}</span>}
        </div>
      </div>

      <div className="metrics">
        <Metric k="Variant" mono v={variant} dClass={refMismatch ? "bad" : ""}
                d={!variant ? "Not positioned" : refMismatch ? `Protein has ${protRes}${s.positions[0]}, not ${s.ref}`
                  : `p.${AA3[s.ref]}${s.positions[0]}${AA3[s.alt]}${s.positions.length > 1 ? ` · ${s.positions.length} repeat sites` : ""}`} />
        <Metric k="BLOSUM62" v={s.blosum62} d={s.blosum62 == null ? null : s.blosum62 > 0 ? "Frequent substitution" : s.blosum62 === 0 ? "Neutral" : "Rare substitution"} title="NCBI BLOSUM62 log-odds score" />
        <Metric k="Grantham" v={s.grantham} d={granthamClass(s.grantham)} title="Grantham (1974) chemical distance" />
        <Metric k="pLDDT at site" v={site && site.plddt != null ? site.plddt.toFixed(1) : acc && !ann ? <Spinner /> : null}
                d={site && site.plddt != null ? plddtBand(site.plddt)[1] : af && af.issue} title={PLDDT_TITLE} />
        <Metric k="AlphaMissense" v={site && site.am_pathogenicity != null ? site.am_pathogenicity.toFixed(3) : acc && !ann ? <Spinner /> : null}
                d={site && site.am_class ? AM_CLASS[site.am_class] : protRes && protRes === s.alt ? `Not scored: protein already has ${s.alt}` : af && !af.am_mean && !af.issue ? "Not available for this proteome" : af && af.issue}
                dClass={site && site.am_class === "LPath" ? "bad" : ""} title="AlphaMissense pathogenicity · likely benign < 0.34 · ambiguous 0.34–0.564 · likely pathogenic > 0.564" />
        <KnownVariantMetric s={s} site={site} />
        <GnomadMetric s={s} human={species.includes("Homo sapiens")} />
      </div>

      <div className="card">
        <div className="card-head"><h2>Alignment</h2><span className="muted">{s.at_cleavage_site ? "Substitution sits at a protease cleavage site" : ""}</span></div>
        <Alignment saap={s} sequence={view && view.protein.sequence} />
      </div>

      {view && (
        <div className="card">
          <ProteinMap length={view.protein.length} sequence={view.protein.sequence} sites={view.sites} focus={s.id} annotations={ann}
                      onSelect={(g) => { const top = [...g.saaps].sort((a, b) => b.n_observations - a.n_observations)[0]; if (top.id !== s.id) go(`saap/${top.id}`); }} />
        </div>
      )}

      {view && s.positions.length > 0 && ann && <SiteContext pos={s.positions[0]} features={ann.uniprot.features} />}
      {view && <div className="card"><SequenceView sequence={view.protein.sequence} sites={view.sites} focus={s.id} plddt={af && af.plddt} /></div>}

      <div className="card">
        <Tabs value={tab} onChange={setTab} tabs={[
          ["observations", "Observations", detail.observations.length],
          ...(af && af.has_structure ? [["structure", "Structure"]] : []),
          ["details", "Details"],
        ]} />
        {tab === "observations" && <ObservationsTable observations={detail.observations} />}
        {tab === "structure" && view && af && af.has_structure && <StructureView accession={acc} sites={view.sites} focus={s.id} />}
        {tab === "details" && (
          <div className="kv">
            {[["Base peptide", s.bp_seq], ["UniProt (source)", s.source_accession], ["Protein accession", s.protein_accession],
              ["Gene", s.source_gene], ["RefProteins", s.ref_proteins], ["Ensembl gene", s.ensembl_gene],
              ["Ensembl transcript", s.ensembl_transcript], ["Ensembl protein", s.ensembl_protein],
              ["Positions", s.positions_all || s.position_in_protein], ["Peptide start", s.peptide_start],
              ["Protein length", s.protein_length], ["Annotation source", s.annotation_source]]
              .map(([k, v]) => <Fragment key={k}><div className="k">{k}</div><div className="v">{fmt.text(v)}</div></Fragment>)}
            {[["Trypsin", s.trypsin], ["Missed cleavage", s.missed_cleavage],
              ["AAS at peptide terminus", s.aas_at_peptide_terminus], ["Greater than shared", s.greater_than_shared],
              ["At cleavage site (computed)", s.at_cleavage_site]]
              .map(([k, v]) => <Fragment key={k}><div className="k">{k}</div><div className="v">{fmt.bool(v)}</div></Fragment>)}
          </div>
        )}
      </div>
    </div>
  );
}

/* ============================== Datasets ============================== */
/** Browse filtered to one tissue / cell type (an organ includes its sub-sites). */
const browseTissue = (name, species) => {
  try {
    sessionStorage.setItem("browse.filters", JSON.stringify({ ...DEFAULT_FILTERS, tissue: name, species: species || "" }));
    sessionStorage.setItem("browse.page", "1");
  } catch { /* unavailable */ }
  go("browse");
};

function GroupTable({ rows, label, clickable }) {
  // Rows sorted by species, then organ; an organ heading precedes its sub-sites.
  const sorted = [...rows].sort((a, b) => a.species.localeCompare(b.species) || a.name.localeCompare(b.name));
  return (
    <div className="table-wrap">
      <table>
        <thead><tr><th>{label}</th><th>Species</th><th className="num">SAAP</th><th className="num">Observations</th><th>Positioned</th><th>Digest</th><th>Acquisition</th></tr></thead>
        <tbody>{sorted.map((d) => {
          const pct = d.n_saap ? Math.round((100 * d.n_annotated) / d.n_saap) : 0;
          const [org, sub] = d.name.split(" (");
          return (
            <tr key={d.name + d.species} className={clickable ? "link" : ""} onClick={() => clickable && browseTissue(d.name, d.species)}>
              <td><b>{org}</b>{sub && <span className="muted"> ({sub}</span>}{d.sample_type === "cell type" && <span className="chip" style={{ marginLeft: 8 }}>cell type</span>}</td>
              <td><i>{shortSpecies(d.species || "")}</i></td>
              <td className="num">{d.n_saap.toLocaleString()}</td>
              <td className="num">{d.n_observations.toLocaleString()}</td>
              <td><span className="inline-bar"><span className="bar-track"><span className="bar-fill" style={{ width: `${pct}%`, display: "block" }} /></span>{pct}%</span></td>
              <td>{chips(d.digests)}</td><td>{chips(d.acquisition_types)}</td>
            </tr>
          );
        })}</tbody>
      </table>
    </div>
  );
}

function TissuesPage({ stats }) {
  const toast = useToast();
  const [ov, error, reload] = useAsync(api.datasets, []);
  useDataChanged(reload);
  const wipe = async () => {
    if (!confirm("Delete ALL SAAP and observations? This cannot be undone.")) return;
    try { const r = await api.remove({ all: true }); toast(`Cleared ${plural(r.deleted, "SAAP")}`); dataChanged(); }
    catch (e) { toast(e.message, true); }
  };
  if (error) return <ErrorNote error={error} />;
  if (!ov) return <Loading />;
  if (!ov.tissues.length) return <div className="card empty">No data yet. <a href="#/import">Import a file</a>.</div>;
  // Organ-level totals for the trend chart (sub-sites grouped under their organ).
  const organs = Object.values(ov.tissues.reduce((m, t) => {
    const o = t.name.split(" (")[0];
    m[o] = m[o] || { label: o, n: 0 }; m[o].n += t.n_saap; return m;
  }, {})).sort((a, b) => b.n - a.n).slice(0, 12);

  return (
    <div className="page stack">
      {stats && <div className="metrics">
        <Metric k="SAAP" v={stats.n_saap.toLocaleString()} />
        <Metric k="Observations" v={stats.n_observations.toLocaleString()} />
        <Metric k="Tissues & cell types" v={stats.n_tissues.toLocaleString()} />
        <Metric k="Proteins" v={stats.n_proteins.toLocaleString()} />
        <Metric k="Genes" v={stats.n_genes.toLocaleString()} />
      </div>}
      <div className="grid-2">
        {[["Top substitutions", ov.top_substitutions], ["Organs", organs],
          ["Species", ov.species_distribution], ["Acquisition", ov.acquisition_distribution]].map(([t, items]) => (
          <div className="card" key={t}><div className="card-head"><h2>{t}</h2><span className="muted">distinct SAAP</span></div>
            <div className="card-body"><BarList items={items} /></div></div>
        ))}
      </div>
      <div className="card">
        <div className="card-head"><h2>Tissues &amp; cell types</h2><span className="muted">{ov.tissues.length} by species</span></div>
        <GroupTable rows={ov.tissues} label="Tissue / cell type" clickable />
      </div>
      {ov.datasets && (
        <div className="card">
          <div className="card-head"><h2>Datasets</h2><span className="muted">private · {ov.datasets.length}</span></div>
          <GroupTable rows={ov.datasets} label="Dataset" />
        </div>
      )}
      {stats && stats.private && <div className="row"><span className="spacer" /><button className="ghost danger sm" onClick={wipe}>Clear all data</button></div>}
    </div>
  );
}

/* =============================== Import =============================== */
function ImportPage() {
  const toast = useToast();
  const [drag, setDrag] = useState(false);
  const [busy, setBusy] = useState(null);
  const [result, setResult] = useState(null);
  const [lastRun, setLastRun] = useState(null);
  const [status, , reloadStatus] = useAsync(api.annotateStatus, []);
  const [curation, , reloadCuration] = useAsync(api.curation, []);
  useDataChanged(useCallback(() => { reloadStatus(); reloadCuration(); }, []));
  const input = useRef();

  const upload = async (file) => {
    if (!file) return;
    setBusy("upload");
    try {
      const r = await api.upload(file); setResult(r);
      toast(`Imported ${plural(r.rows_read, "row")}${r.rows_skipped_immunoglobulin ? ` · ${r.rows_skipped_immunoglobulin.toLocaleString()} immunoglobulin rows skipped` : ""}`);
      dataChanged();
    }
    catch (e) { toast(e.message, true); } finally { setBusy(null); input.current.value = ""; }
  };
  const annotate = async (overwrite) => {
    if (overwrite && !confirm("Re-fetch every SAAP from UniProt and overwrite existing annotation? This can take a while.")) return;
    setBusy(overwrite ? "all" : "new");
    try { const r = await api.annotate(overwrite ? { overwrite: true } : {}); setLastRun(r); toast(annotateSummary(r), r.failed > 0 && !r.positioned); dataChanged(); }
    catch (e) { toast(e.message, true); } finally { setBusy(null); }
  };

  return (
    <div className="page stack" style={{ maxWidth: 880 }}>
      <div className="card">
        <div className="card-head"><h2>Import</h2><span className="muted">CSV, TSV or XLSX</span></div>
        <div className="card-body stack">
          <div className={"drop" + (drag ? " drag" : "")} onClick={() => input.current.click()}
               onDragOver={(e) => { e.preventDefault(); setDrag(true); }} onDragLeave={() => setDrag(false)}
               onDrop={(e) => { e.preventDefault(); setDrag(false); upload(e.dataTransfer.files[0]); }}>
            {busy === "upload" ? <Loading label="Importing" /> : <span><strong>Drop a file</strong> or click to browse</span>}
            <input ref={input} type="file" accept=".csv,.tsv,.txt,.xlsx" hidden onChange={(e) => upload(e.target.files[0])} />
          </div>
          {result && (
            <div className="stack" style={{ gap: 6 }}>
              <div className="stats-line">
                <span><b>{result.filename}</b></span>
                <span>Rows <b>{result.rows_read.toLocaleString()}</b></span>
                <span>New SAAP <b>{result.saap_created.toLocaleString()}</b></span>
                <span>Observations <b>{result.observations_created.toLocaleString()}</b></span>
                <span>Duplicates skipped <b>{result.duplicate_observations_skipped.toLocaleString()}</b></span>
              </div>
              {result.columns_unmapped.length > 0 && <div className="warn-text">Unmapped columns: {result.columns_unmapped.join(", ")}</div>}
            </div>
          )}
        </div>
      </div>

      <div className="card">
        <div className="card-head"><h2>Annotation</h2><span className="muted">UniProt · Ensembl</span></div>
        <div className="card-body stack">
          {status ? (
            <div className="stats-line">
              <span>Positioned <b>{status.n_with_position.toLocaleString()}</b> of {status.n_saap.toLocaleString()}</span>
              <span>Awaiting annotation <b>{status.n_needs_annotation.toLocaleString()}</b></span>
            </div>
          ) : <Loading />}
          <div className="row">
            <button className="primary" disabled={!!busy || !status || !status.n_needs_annotation} onClick={() => annotate(false)}>
              {busy === "new" ? <Spinner /> : <Icon name="spark" />}Annotate {status ? status.n_needs_annotation.toLocaleString() : ""} pending</button>
            <button className="ghost" disabled={!!busy || !status || !status.n_saap} onClick={() => annotate(true)}>
              {busy === "all" && <Spinner />}Re-annotate all</button>
          </div>
          {lastRun && <div className="muted">{annotateSummary(lastRun)}</div>}
        </div>
      </div>

      <div className="card">
        <div className="card-head"><h2>Curation log</h2></div>
        {curation ? (
          <div className="card-body stack" style={{ gap: 8 }}>
            <div className="stats-line">
              <span>PosProb &lt; {curation.min_positional_probability} <b>{curation.below_min_positional_probability.toLocaleString()}</b></span>
              <span>No PosProb <b>{curation.no_positional_probability.toLocaleString()}</b></span>
              <span>gnomAD absent <b>{(curation.gnomad.absent || 0).toLocaleString()}</b></span>
              <span>rare <b>{(curation.gnomad.present || 0).toLocaleString()}</b></span>
            </div>
            <div className="muted">Removed: immunoglobulins, contaminants, genome-encoded peptides, gnomAD AF ≥ 0.0001, unconfirmed known variants</div>
          </div>
        ) : <Loading />}
      </div>
    </div>
  );
}

/* =============================== Export =============================== */
const HEADER_TEMPLATES = {
  peptide: ">sp|{accession}-{mid}-{tok}|{gene}-mut {gene} substituted {mid} OS={species} OX={taxid} GN={gene} PE=1 SV=1",
  protein: ">sp|{accession}-{mid}-{tok}|{gene}-mut {gene} substituted {mid} {sub_compact}@{position} OS={species} OX={taxid} GN={gene} PE=1 SV=1",
  base: ">sp|{accession}-{bid}-{tok}|{gene}-base {gene} base peptide {bid} OS={species} OX={taxid} GN={gene} PE=1 SV=1",
};

function ExportModal({ scope, onClose }) {
  const toast = useToast();
  const [format, setFormat] = useState("fasta");
  const [entry, setEntry] = useState("peptide");
  const [decoys, setDecoys] = useState(false);
  const [base, setBase] = useState(false);
  const [reference, setReference] = useState(null);
  const [template, setTemplate] = useState(HEADER_TEMPLATES.peptide);
  const [baseTemplate, setBaseTemplate] = useState(HEADER_TEMPLATES.base);
  const [busy, setBusy] = useState(false);
  const refInput = useRef();

  const changeEntry = (m) => {
    setEntry(m);
    setTemplate((t) => (t === HEADER_TEMPLATES.peptide || t === HEADER_TEMPLATES.protein ? HEADER_TEMPLATES[m] : t));
  };
  const run = async () => {
    setBusy(true);
    try {
      const target = scope.ids ? { ids: scope.ids } : { filters: scope.filters };
      const payload = format === "fasta"
        ? { ...target, entry_mode: entry, decoys, base_peptides: entry === "peptide" && base, header_template: template, base_header_template: baseTemplate }
        : target;
      const { blob, name, skipped } = await api.exportFile(format, payload, format === "fasta" ? reference : null);
      download(blob, name);
      toast(`Exported ${plural(scope.count, "SAAP")}${skipped ? ` · ${skipped.toLocaleString()} skipped (no gene, species mismatch, or not positioned)` : ""}`, skipped > 0 && skipped >= scope.count);
      onClose();
    } catch (e) { toast(e.message, true); } finally { setBusy(false); }
  };

  return (
    <Modal title={`Export ${plural(scope.count, "SAAP")}`} onClose={onClose}
           footer={<Fragment><button className="ghost" onClick={onClose}>Cancel</button>
             <button className="primary" onClick={run} disabled={busy}>{busy ? <Spinner /> : <Icon name="download" />}Download</button></Fragment>}>
      <div className="field"><label>Format</label>
        <Segmented value={format} onChange={setFormat} options={[["fasta", "FASTA"], ["csv", "CSV"], ["pairs", "SAAP–BP pairs"]]} /></div>
      {format === "fasta" && (
        <Fragment>
          <div className="field"><label>Entries</label>
            <Segmented value={entry} onChange={changeEntry} options={[["peptide", "Substituted peptides"], ["protein", "Full-length proteins"]]} /></div>
          <div className="stack" style={{ gap: 8 }}>
            {entry === "peptide" && <label className="check"><input type="checkbox" checked={base} onChange={(e) => setBase(e.target.checked)} />Include base peptides</label>}
            <label className="check"><input type="checkbox" checked={decoys} onChange={(e) => setDecoys(e.target.checked)} />Append reversed decoys (<code>rev_</code>)</label>
          </div>
          <div className="field"><label>Reference proteome (optional)</label>
            <div className="drop compact" onClick={() => refInput.current.click()}
                 onDragOver={(e) => e.preventDefault()} onDrop={(e) => { e.preventDefault(); setReference(e.dataTransfer.files[0] || null); }}>
              {reference ? <span><strong>{reference.name}</strong> · <a onClick={(e) => { e.stopPropagation(); setReference(null); refInput.current.value = ""; }}>remove</a></span>
                : <span><strong>Drop a FASTA</strong> or click to browse</span>}
              <input ref={refInput} type="file" accept=".fasta,.fa,.faa,.txt" hidden onChange={(e) => setReference(e.target.files[0] || null)} />
            </div>
          </div>
          <details><summary>Header templates</summary>
            <div className="stack" style={{ gap: 8 }}>
              <textarea rows={3} value={template} onChange={(e) => setTemplate(e.target.value)} />
              {entry === "peptide" && base && <textarea rows={3} value={baseTemplate} onChange={(e) => setBaseTemplate(e.target.value)} />}
            </div>
          </details>
        </Fragment>
      )}
    </Modal>
  );
}
