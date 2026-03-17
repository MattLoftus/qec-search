import type { StatusSummary, BenchmarkResult, TabId } from "../types";

interface StatCardProps {
  label: string;
  value: string | number;
  color?: string;
  onClick?: () => void;
}

function StatCard({ label, value, color, onClick }: StatCardProps) {
  return (
    <div
      onClick={onClick}
      className={`bg-[var(--color-bg-card)] border border-[var(--color-border-dim)] rounded-lg p-4 ${
        onClick
          ? "cursor-pointer hover:border-[var(--color-border-default)] transition-colors duration-100"
          : ""
      }`}
    >
      <div className="text-xs text-[var(--color-text-muted)] uppercase tracking-wide mb-1">
        {label}
      </div>
      <div className={`text-2xl font-bold ${color ?? "text-[var(--color-text-primary)]"}`}>
        {value}
      </div>
    </div>
  );
}

interface StatsPanelProps {
  status: StatusSummary | null;
  benchmark: BenchmarkResult | null;
  bestCodesCount: number;
  onNavigate: (tab: TabId) => void;
}

export function StatsPanel({ status, benchmark, bestCodesCount, onNavigate }: StatsPanelProps) {
  if (!status) {
    return (
      <div className="text-[var(--color-text-muted)] p-4">
        Loading status...
      </div>
    );
  }

  return (
    <div className="grid grid-cols-2 md:grid-cols-4 lg:grid-cols-6 gap-3 p-4">
      <StatCard
        label="Codes Tested"
        value={status.total_codes_tested.toLocaleString()}
      />
      <StatCard
        label="Parameters Explored"
        value={bestCodesCount}
      />
      <StatCard
        label="Matches"
        value={benchmark?.summary.matches_best ?? 0}
        color={
          (benchmark?.summary.matches_best ?? 0) > 0
            ? "text-[var(--color-accent-blue)]"
            : "text-[var(--color-text-primary)]"
        }
      />
      <StatCard
        label="Improvements"
        value={benchmark?.summary.improvements ?? 0}
        color={
          (benchmark?.summary.improvements ?? 0) > 0
            ? "text-[var(--color-accent-green)]"
            : "text-[var(--color-text-primary)]"
        }
      />
      <StatCard
        label="Open Gaps"
        value={status.open_gaps}
        color="text-[var(--color-accent-amber)]"
        onClick={() => onNavigate("bounds")}
      />
      <StatCard
        label="Campaigns"
        value={status.total_campaigns}
        onClick={() => onNavigate("campaigns")}
      />
    </div>
  );
}
