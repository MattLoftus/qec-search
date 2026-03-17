import type { CodeEntry, BoundsEntry } from "../types";

interface Props {
  codes: CodeEntry[];
  bounds: BoundsEntry[];
  selectedHash: string | null;
  onSelect: (hash: string | null) => void;
}

export function BestCodesTable({ codes, bounds, selectedHash, onSelect }: Props) {
  // Index bounds for lookup
  const boundsMap = new Map<string, BoundsEntry>();
  for (const b of bounds) {
    boundsMap.set(`${b.n},${b.k}`, b);
  }

  if (codes.length === 0) {
    return (
      <div className="p-6 text-center text-[var(--color-text-muted)]">
        No codes found yet. Run a search campaign to populate results.
      </div>
    );
  }

  return (
    <div className="p-4">
      <h2 className="text-lg font-semibold text-[var(--color-text-primary)] mb-3">
        Best Codes Found
      </h2>
      <div className="overflow-x-auto">
        <table className="w-full text-sm">
          <thead>
            <tr className="border-b border-[var(--color-border-default)]">
              <th className="text-left p-2 text-[var(--color-text-secondary)]">Code</th>
              <th className="text-left p-2 text-[var(--color-text-secondary)]">n</th>
              <th className="text-left p-2 text-[var(--color-text-secondary)]">k</th>
              <th className="text-left p-2 text-[var(--color-text-secondary)]">d</th>
              <th className="text-left p-2 text-[var(--color-text-secondary)]">Known</th>
              <th className="text-left p-2 text-[var(--color-text-secondary)]">Status</th>
              <th className="text-left p-2 text-[var(--color-text-secondary)]">Rate</th>
              <th className="text-left p-2 text-[var(--color-text-secondary)]">Max Wt</th>
              <th className="text-left p-2 text-[var(--color-text-secondary)]">Strategy</th>
            </tr>
          </thead>
          <tbody>
            {codes.map((code) => {
              const bound = boundsMap.get(`${code.n},${code.k}`);
              const isSelected = code.hash === selectedHash;

              let status = "—";
              let statusColor = "text-[var(--color-text-muted)]";
              if (bound && code.d !== null) {
                if (code.d > bound.d_lower) {
                  status = "IMPROVEMENT";
                  statusColor = "text-[var(--color-accent-green)]";
                } else if (code.d === bound.d_lower) {
                  status = "MATCHES";
                  statusColor = "text-[var(--color-accent-blue)]";
                } else {
                  status = "BELOW";
                  statusColor = "text-[var(--color-accent-amber)]";
                }
              }

              return (
                <tr
                  key={code.hash || `${code.n}-${code.k}-${code.d}`}
                  onClick={() => onSelect(isSelected ? null : code.hash)}
                  className={`border-b border-[var(--color-border-dim)] cursor-pointer transition-colors duration-100 ${
                    isSelected
                      ? "bg-[var(--color-accent-blue)]/10"
                      : "hover:bg-white/[0.02]"
                  }`}
                >
                  <td className="p-2 font-mono text-[var(--color-text-primary)]">
                    [[{code.n}, {code.k}, {code.d ?? "?"}]]
                  </td>
                  <td className="p-2 text-[var(--color-text-secondary)]">{code.n}</td>
                  <td className="p-2 text-[var(--color-text-secondary)]">{code.k}</td>
                  <td className="p-2 text-[var(--color-text-primary)] font-bold">{code.d ?? "?"}</td>
                  <td className="p-2 text-[var(--color-text-muted)]">
                    {bound ? `[${bound.d_lower}, ${bound.d_upper}]` : "—"}
                  </td>
                  <td className={`p-2 text-xs font-medium ${statusColor}`}>{status}</td>
                  <td className="p-2 text-[var(--color-text-secondary)]">
                    {code.rate.toFixed(3)}
                  </td>
                  <td className="p-2 text-[var(--color-text-secondary)]">{code.max_weight}</td>
                  <td className="p-2 text-[var(--color-text-muted)] text-xs">{code.strategy ?? "—"}</td>
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>

      {/* Selected code detail */}
      {selectedHash && (() => {
        const code = codes.find((c) => c.hash === selectedHash);
        if (!code) return null;
        return (
          <div className="mt-4 bg-[var(--color-bg-card)] border border-[var(--color-border-dim)] rounded-lg p-4">
            <h3 className="font-medium text-[var(--color-text-primary)] mb-2">
              Code Detail: [[{code.n}, {code.k}, {code.d ?? "?"}]]
            </h3>
            <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
              <div>
                <h4 className="text-xs text-[var(--color-text-muted)] uppercase mb-1">Hx ({code.Hx.length} rows)</h4>
                <pre className="text-xs text-[var(--color-text-secondary)] bg-white/[0.02] p-2 rounded overflow-x-auto">
                  {code.Hx.map((row) => row.join(" ")).join("\n")}
                </pre>
              </div>
              <div>
                <h4 className="text-xs text-[var(--color-text-muted)] uppercase mb-1">Hz ({code.Hz.length} rows)</h4>
                <pre className="text-xs text-[var(--color-text-secondary)] bg-white/[0.02] p-2 rounded overflow-x-auto">
                  {code.Hz.map((row) => row.join(" ")).join("\n")}
                </pre>
              </div>
            </div>
            <div className="mt-2 text-xs text-[var(--color-text-muted)]">
              Hash: {code.hash} | Avg weight: {code.avg_weight.toFixed(1)} | Max weight: {code.max_weight}
            </div>
          </div>
        );
      })()}
    </div>
  );
}
