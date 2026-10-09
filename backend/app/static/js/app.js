/* App shell: navigation, global stats, toast, routing. */
const NAV = [["browse", "Browse"], ["proteins", "Proteins"], ["datasets", "Datasets"], ["import", "Import"]];
const NAV_PARENT = { saap: "browse", protein: "proteins" };

function App() {
  const route = useRoute();
  const [stats, , reloadStats] = useAsync(api.stats, [], true);
  const [facets, , reloadFacets] = useAsync(api.facets, [], true);
  const [toast, setToast] = useState(null);
  const refresh = useCallback(() => { reloadStats(); reloadFacets(); }, []);
  useDataChanged(refresh);

  const timer = useRef();
  const showToast = useCallback((msg, err = false) => {
    clearTimeout(timer.current);
    setToast({ msg, err });
    timer.current = setTimeout(() => setToast(null), 5000);
  }, []);

  const section = NAV_PARENT[route.page] || route.page;
  useEffect(() => {
    const label = { saap: "SAAP", protein: route.param }[route.page] || (NAV.find(([k]) => k === route.page) || [])[1];
    document.title = label ? `${label} · SAAPedia` : "SAAPedia";
  }, [route.page, route.param]);

  let page;
  switch (route.page) {
    case "proteins": page = <ProteinsPage />; break;
    case "protein": page = <ProteinPage accession={route.param} initialTab={route.query.get("tab")} />; break;
    case "saap": page = <SaapPage id={route.param} initialTab={route.query.get("tab")} />; break;
    case "datasets": page = <DatasetsPage stats={stats} />; break;
    case "import": page = <ImportPage />; break;
    default: page = <BrowsePage facets={facets || {}} />;
  }

  return (
    <ToastContext.Provider value={showToast}>
      <header className="topbar">
        <a className="brand" href="#/browse">
          <svg viewBox="0 0 64 64" aria-hidden="true">
            <path d="M3 41 Q16 38 26 41 T50 41 Q56 41 61 39" />
            <ellipse cx="27" cy="42" rx="14" ry="8.5" />
            <path d="M14 32 Q13 17 29 14 Q44 13 45 27 Q46 36 34 38 Q20 40 14 32 Z" />
            <path d="M37 18 Q46 9 58 12" />
            <circle cx="45" cy="11" r="3" /><circle cx="58" cy="12" r="3" />
          </svg>
          SAAPedia
        </a>
        <nav className="nav">
          {NAV.map(([k, label]) => <a key={k} href={`#/${k}`} className={section === k ? "on" : ""}>{label}</a>)}
        </nav>
        {stats && <div className="meta">{stats.n_saap.toLocaleString()} SAAP · {stats.n_proteins.toLocaleString()} proteins · {stats.n_observations.toLocaleString()} observations</div>}
        <a className="meta-link" href="https://github.com/alexmaropakis/SAAPedia" target="_blank" rel="noopener noreferrer" title="GitHub"><Icon name="external" /></a>
      </header>
      <main key={route.page + (route.param || "")}>{page}</main>
      {toast && <div className={"toast" + (toast.err ? " err" : "")} onClick={() => setToast(null)}>{toast.msg}</div>}
    </ToastContext.Provider>
  );
}

ReactDOM.createRoot(document.getElementById("root")).render(<App />);
