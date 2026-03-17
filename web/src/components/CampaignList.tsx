import type { Campaign } from "../types";

interface Props {
  campaigns: Campaign[];
  activeIndex: number;
  onSelect: (index: number) => void;
}

export function CampaignList({ campaigns, activeIndex, onSelect }: Props) {
  if (campaigns.length === 0) {
    return (
      <div className="p-6 text-center text-[var(--color-text-muted)]">
        <p className="mb-2">No campaigns yet.</p>
        <p className="text-xs">
          Run: <code className="bg-white/[0.05] px-1 rounded">python3 engine/evolver.py --mode random --n 10 --k 1</code>
        </p>
      </div>
    );
  }

  return (
    <div className="p-4">
      <h2 className="text-lg font-semibold text-[var(--color-text-primary)] mb-3">
        Search Campaigns
      </h2>
      <div className="space-y-3">
        {campaigns.map((campaign, idx) => {
          const isActive = idx === activeIndex;
          const bestCodes = campaign.best_codes ?? [];
          const bestD = bestCodes.length > 0
            ? Math.max(...bestCodes.map((c) => c.d ?? 0))
            : 0;

          return (
            <div
              key={campaign.name}
              onClick={() => onSelect(idx)}
              className={`bg-[var(--color-bg-card)] border rounded-lg p-4 cursor-pointer transition-colors duration-100 ${
                isActive
                  ? "border-[var(--color-accent-blue)]"
                  : "border-[var(--color-border-dim)] hover:border-[var(--color-border-default)]"
              }`}
            >
              <div className="flex items-center justify-between mb-2">
                <h3 className="font-medium text-[var(--color-text-primary)]">
                  {campaign.name}
                </h3>
                <span className="text-xs text-[var(--color-text-muted)] bg-white/[0.05] px-2 py-0.5 rounded">
                  {campaign.strategy}
                </span>
              </div>
              <div className="grid grid-cols-2 md:grid-cols-4 gap-2 text-xs">
                <div>
                  <span className="text-[var(--color-text-muted)]">Tested: </span>
                  <span className="text-[var(--color-text-secondary)]">
                    {campaign.total_tested.toLocaleString()}
                  </span>
                </div>
                <div>
                  <span className="text-[var(--color-text-muted)]">Valid: </span>
                  <span className="text-[var(--color-text-secondary)]">
                    {campaign.codes_with_k}
                  </span>
                </div>
                <div>
                  <span className="text-[var(--color-text-muted)]">Best d: </span>
                  <span className="text-[var(--color-text-primary)] font-bold">
                    {bestD > 0 ? bestD : "—"}
                  </span>
                </div>
                <div>
                  <span className="text-[var(--color-text-muted)]">Time: </span>
                  <span className="text-[var(--color-text-secondary)]">
                    {campaign.elapsed_sec < 60
                      ? `${campaign.elapsed_sec.toFixed(1)}s`
                      : `${(campaign.elapsed_sec / 60).toFixed(1)}m`}
                  </span>
                </div>
              </div>
              {bestCodes.length > 0 && (
                <div className="mt-2 flex flex-wrap gap-1">
                  {bestCodes.slice(0, 8).map((code, ci) => (
                    <span
                      key={ci}
                      className={`text-xs px-1.5 py-0.5 rounded ${
                        code.beats_known
                          ? "bg-[var(--color-accent-green)]/20 text-[var(--color-accent-green)]"
                          : "bg-white/[0.05] text-[var(--color-text-muted)]"
                      }`}
                    >
                      [[{code.n},{code.k},{code.d ?? "?"}]]
                    </span>
                  ))}
                  {bestCodes.length > 8 && (
                    <span className="text-xs text-[var(--color-text-muted)]">
                      +{bestCodes.length - 8} more
                    </span>
                  )}
                </div>
              )}
            </div>
          );
        })}
      </div>
    </div>
  );
}
