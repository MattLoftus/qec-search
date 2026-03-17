import { useEffect } from "react";
import { useStore } from "./store";
import { api } from "./api";
import { StatsPanel } from "./components/StatsPanel";
import { BestCodesTable } from "./components/BestCodesTable";
import { CampaignList } from "./components/CampaignList";
import { BoundsTable } from "./components/BoundsTable";
import { Guide } from "./components/Guide";
import { DistancePlot } from "./components/DistancePlot";
import type { TabId, StatusSummary, Campaign, CodeEntry, BenchmarkResult, BoundsEntry } from "./types";

const TABS: { id: TabId; label: string }[] = [
  { id: "dashboard", label: "Dashboard" },
  { id: "campaigns", label: "Campaigns" },
  { id: "bounds", label: "Bounds" },
  { id: "guide", label: "Guide" },
];

function App() {
  const activeTab = useStore((s) => s.activeTab);
  const setActiveTab = useStore((s) => s.setActiveTab);
  const status = useStore((s) => s.status);
  const campaigns = useStore((s) => s.campaigns);
  const bestCodes = useStore((s) => s.bestCodes);
  const benchmark = useStore((s) => s.benchmark);
  const bounds = useStore((s) => s.bounds);
  const activeCampaignIndex = useStore((s) => s.activeCampaignIndex);
  const selectedCodeHash = useStore((s) => s.selectedCodeHash);
  const loading = useStore((s) => s.loading);
  const error = useStore((s) => s.error);

  const setStatus = useStore((s) => s.setStatus);
  const setCampaigns = useStore((s) => s.setCampaigns);
  const setBestCodes = useStore((s) => s.setBestCodes);
  const setBenchmark = useStore((s) => s.setBenchmark);
  const setBounds = useStore((s) => s.setBounds);
  const setActiveCampaign = useStore((s) => s.setActiveCampaign);
  const setSelectedCode = useStore((s) => s.setSelectedCode);
  const setLoading = useStore((s) => s.setLoading);
  const setError = useStore((s) => s.setError);

  useEffect(() => {
    async function loadAll() {
      setLoading(true);
      setError(null);
      try {
        const [statusData, campaignsData, bestData, benchmarkData, boundsData] =
          await Promise.all([
            api.status() as Promise<StatusSummary>,
            api.campaigns() as Promise<Campaign[]>,
            api.best() as Promise<CodeEntry[]>,
            api.benchmark() as Promise<BenchmarkResult>,
            api.bounds() as Promise<BoundsEntry[]>,
          ]);
        setStatus(statusData);
        setCampaigns(campaignsData);
        setBestCodes(bestData);
        setBenchmark(benchmarkData);
        setBounds(boundsData);
      } catch (e) {
        setError(e instanceof Error ? e.message : "Failed to connect to API");
      } finally {
        setLoading(false);
      }
    }
    loadAll();
  }, []);

  return (
    <div className="min-h-screen bg-[var(--color-bg-primary)] text-[var(--color-text-primary)]">
      <header className="border-b border-[var(--color-border-dim)] bg-[var(--color-bg-secondary)]">
        <div className="max-w-7xl mx-auto px-4 py-3 flex items-center justify-between">
          <div className="flex items-center gap-3">
            <h1 className="text-lg font-bold tracking-tight">QEC Search</h1>
            <span className="text-xs text-[var(--color-text-muted)]">
              Quantum Error Correction Code Search
            </span>
          </div>
          <nav className="flex gap-1">
            {TABS.map((tab) => (
              <button
                key={tab.id}
                onClick={() => setActiveTab(tab.id)}
                className={`px-3 py-1.5 rounded text-sm transition-colors duration-100 ${
                  activeTab === tab.id
                    ? "bg-[var(--color-accent-blue)] text-white"
                    : "text-[var(--color-text-secondary)] hover:text-[var(--color-text-primary)] hover:bg-white/[0.05]"
                }`}
              >
                {tab.label}
              </button>
            ))}
          </nav>
        </div>
      </header>

      {error && (
        <div className="bg-[var(--color-accent-red)]/10 border-b border-[var(--color-accent-red)]/30 px-4 py-2">
          <p className="text-sm text-[var(--color-accent-red)] max-w-7xl mx-auto">
            {error}
            <span className="text-[var(--color-text-muted)] ml-2">
              — Start the API server: python3 engine/serve.py
            </span>
          </p>
        </div>
      )}

      {loading && (
        <div className="text-center py-8 text-[var(--color-text-muted)]">
          Loading...
        </div>
      )}

      <main className="max-w-7xl mx-auto">
        {activeTab === "dashboard" && (
          <>
            <StatsPanel
              status={status}
              benchmark={benchmark}
              bestCodesCount={bestCodes.length}
              onNavigate={setActiveTab}
            />
            <DistancePlot codes={bestCodes} bounds={bounds} />
            <BestCodesTable
              codes={bestCodes}
              bounds={bounds}
              selectedHash={selectedCodeHash}
              onSelect={setSelectedCode}
            />
          </>
        )}
        {activeTab === "campaigns" && (
          <CampaignList
            campaigns={campaigns}
            activeIndex={activeCampaignIndex}
            onSelect={setActiveCampaign}
          />
        )}
        {activeTab === "bounds" && (
          <BoundsTable bounds={bounds} bestCodes={bestCodes} />
        )}
        {activeTab === "guide" && <Guide />}
      </main>
    </div>
  );
}

export default App;
