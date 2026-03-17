import { create } from "zustand";
import type {
  Campaign,
  CodeEntry,
  BoundsEntry,
  BenchmarkResult,
  StatusSummary,
  TabId,
} from "./types";

interface AppState {
  // Data
  status: StatusSummary | null;
  campaigns: Campaign[];
  bestCodes: CodeEntry[];
  benchmark: BenchmarkResult | null;
  bounds: BoundsEntry[];

  // UI state
  activeTab: TabId;
  activeCampaignIndex: number;
  selectedCodeHash: string | null;
  loading: boolean;
  error: string | null;

  // Actions
  setStatus: (s: StatusSummary) => void;
  setCampaigns: (c: Campaign[]) => void;
  setBestCodes: (c: CodeEntry[]) => void;
  setBenchmark: (b: BenchmarkResult) => void;
  setBounds: (b: BoundsEntry[]) => void;
  setActiveTab: (t: TabId) => void;
  setActiveCampaign: (i: number) => void;
  setSelectedCode: (hash: string | null) => void;
  setLoading: (l: boolean) => void;
  setError: (e: string | null) => void;
}

export const useStore = create<AppState>((set) => ({
  status: null,
  campaigns: [],
  bestCodes: [],
  benchmark: null,
  bounds: [],

  activeTab: "dashboard",
  activeCampaignIndex: 0,
  selectedCodeHash: null,
  loading: false,
  error: null,

  setStatus: (s) => set({ status: s }),
  setCampaigns: (c) => set({ campaigns: c }),
  setBestCodes: (c) => set({ bestCodes: c }),
  setBenchmark: (b) => set({ benchmark: b }),
  setBounds: (b) => set({ bounds: b }),
  setActiveTab: (t) => set({ activeTab: t }),
  setActiveCampaign: (i) => set({ activeCampaignIndex: i }),
  setSelectedCode: (hash) => set({ selectedCodeHash: hash }),
  setLoading: (l) => set({ loading: l }),
  setError: (e) => set({ error: e }),
}));

// Derived selectors
export const useActiveCampaign = () => {
  const campaigns = useStore((s) => s.campaigns);
  const idx = useStore((s) => s.activeCampaignIndex);
  return campaigns[idx] ?? null;
};

export const useSelectedCode = () => {
  const bestCodes = useStore((s) => s.bestCodes);
  const hash = useStore((s) => s.selectedCodeHash);
  return bestCodes.find((c) => c.hash === hash) ?? null;
};
