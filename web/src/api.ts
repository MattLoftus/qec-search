const API_BASE = "http://localhost:8789";
const STATIC_BASE = import.meta.env.BASE_URL + "data";

async function fetchJson<T>(path: string, staticFile: string): Promise<T> {
  // Try local API first, fall back to static data
  try {
    const res = await fetch(`${API_BASE}${path}`, { signal: AbortSignal.timeout(2000) });
    if (res.ok) return res.json();
  } catch {
    // API unavailable — use static data
  }
  const res = await fetch(`${STATIC_BASE}/${staticFile}`);
  if (!res.ok) throw new Error(`Static ${staticFile}: ${res.status}`);
  return res.json();
}

export const api = {
  status: () => fetchJson("/api/status", "status.json"),
  campaigns: () => fetchJson("/api/campaigns", "campaigns.json"),
  campaign: (name: string) => fetchJson(`/api/campaign/${name}`, `campaign_${name}.json`),
  best: () => fetchJson("/api/best", "best.json"),
  benchmark: () => fetchJson("/api/benchmark", "benchmark.json"),
  bounds: () => fetchJson("/api/bounds", "bounds.json"),
};
