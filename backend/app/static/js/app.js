/* App shell: navigation, global stats, toast, routing. */
const NAV = [["home", "Home"], ["browse", "Browse"], ["proteins", "Proteins"], ["tissues", "Tissues"], ["import", "Import"]];
const NAV_PARENT = { saap: "browse", protein: "proteins", datasets: "tissues" };

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
    case "home": page = <HomePage stats={stats} />; break;
    case "browse": page = <BrowsePage facets={facets || {}} />; break;
    case "proteins": page = <ProteinsPage />; break;
    case "protein": page = <ProteinPage accession={route.param} initialTab={route.query.get("tab")} />; break;
    case "saap": page = <SaapPage id={route.param} initialTab={route.query.get("tab")} />; break;
    case "tissues": case "datasets": page = <TissuesPage stats={stats} />; break;
    case "import": page = <ImportPage />; break;
    default: page = <HomePage stats={stats} />;
  }

  return (
    <ToastContext.Provider value={showToast}>
      <header className="topbar">
        <a className="brand" href="#/">
          <svg className="logo" viewBox="0 0 32 32" aria-hidden="true">
            <rect className="logo-tile" x="1" y="1" width="30" height="30" rx="7" />
            <text className="logo-letter" x="16" y="22.5" textAnchor="middle">S</text>
            <line className="logo-split" x1="1" y1="16" x2="31" y2="16" />
            <path className="logo-mark" d="M22 1 H24 A7 7 0 0 1 31 8 V10 Z" />
          </svg>
          <Wordmark />
        </a>
        <nav className="nav">
          {NAV.filter(([k]) => k !== "import" || (stats && stats.private))
            .map(([k, label]) => <a key={k} href={`#/${k}`} className={section === k ? "on" : ""}>{label}</a>)}
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
