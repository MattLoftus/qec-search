import type { BoundsEntry, CodeEntry } from "../types";

interface Props {
  bounds: BoundsEntry[];
  bestCodes: CodeEntry[];
}

export function BoundsTable({ bounds, bestCodes }: Props) {
  // Index our best codes
  const ourBest = new Map<string, number>();
  for (const c of bestCodes) {
    if (c.d !== null) {
      const key = `${c.n},${c.k}`;
      const existing = ourBest.get(key);
      if (existing === undefined || c.d > existing) {
        ourBest.set(key, c.d);
      }
    }
  }

  if (bounds.length === 0) {
    return (
      <div className="p-6 text-center text-[var(--color-text-muted)]">
        Loading bounds data...
      </div>
    );
  }

  // Group by n
  const byN = new Map<number, BoundsEntry[]>();
  for (const b of bounds) {
    const arr = byN.get(b.n) ?? [];
    arr.push(b);
    byN.set(b.n, arr);
  }

  return (
    <div className="p-4">
      <h2 className="text-lg font-semibold text-[var(--color-text-primary)] mb-1">
        Known Bounds (codetables.de)
      </h2>
      <p className="text-sm text-[var(--color-text-muted)] mb-4">
        Amber rows have open gaps (d_lower &lt; d_upper) — finding a code with d &gt; d_lower
        in these regimes is a publishable result.
      </p>
      <div className="overflow-x-auto">
        <table className="w-full text-sm">
          <thead>
            <tr className="border-b border-[var(--color-border-default)]">
              <th className="text-left p-2 text-[var(--color-text-secondary)]">n</th>
              <th className="text-left p-2 text-[var(--color-text-secondary)]">k</th>
              <th className="text-left p-2 text-[var(--color-text-secondary)]">d lower</th>
              <th className="text-left p-2 text-[var(--color-text-secondary)]">d upper</th>
              <th className="text-left p-2 text-[var(--color-text-secondary)]">Gap</th>
              <th className="text-left p-2 text-[var(--color-text-secondary)]">Our Best</th>
              <th className="text-left p-2 text-[var(--color-text-secondary)]">Status</th>
            </tr>
          </thead>
          <tbody>
            {Array.from(byN.entries())
              .sort(([a], [b]) => a - b)
              .map(([_n, entries]) =>
                entries
                  .sort((a, b) => a.k - b.k)
                  .map((bound) => {
                    const key = `${bound.n},${bound.k}`;
                    const ourD = ourBest.get(key);
                    const hasGap = bound.d_lower < bound.d_upper;

                    let status = "";
                    let statusColor = "";
                    if (ourD !== undefined) {
                      if (ourD > bound.d_lower) {
                        status = "BEATS";
                        statusColor = "text-[var(--color-accent-green)]";
                      } else if (ourD === bound.d_lower) {
                        status = "MATCHES";
                        statusColor = "text-[var(--color-accent-blue)]";
                      } else {
                        status = `−${bound.d_lower - ourD}`;
                        statusColor = "text-[var(--color-accent-amber)]";
                      }
                    }

                    return (
                      <tr
                        key={key}
                        className={`border-b border-[var(--color-border-dim)] ${
                          hasGap ? "bg-[var(--color-accent-amber)]/5" : ""
                        }`}
                      >
                        <td className="p-2 text-[var(--color-text-primary)]">{bound.n}</td>
                        <td className="p-2 text-[var(--color-text-secondary)]">{bound.k}</td>
                        <td className="p-2 text-[var(--color-text-primary)] font-bold">
                          {bound.d_lower}
                        </td>
                        <td className="p-2 text-[var(--color-text-secondary)]">
                          {bound.d_upper}
                        </td>
                        <td className="p-2">
                          {hasGap ? (
                            <span className="text-[var(--color-accent-amber)]">
                              {bound.d_upper - bound.d_lower}
                            </span>
                          ) : (
                            <span className="text-[var(--color-accent-green)]">0</span>
                          )}
                        </td>
                        <td className="p-2 text-[var(--color-text-secondary)]">
                          {ourD ?? "—"}
                        </td>
                        <td className={`p-2 text-xs font-medium ${statusColor}`}>
                          {status || "—"}
                        </td>
                      </tr>
                    );
                  })
              )}
          </tbody>
        </table>
      </div>
    </div>
  );
}
