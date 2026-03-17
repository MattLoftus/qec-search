import Plot from "./PlotlyWrapper";
import type { CodeEntry, BoundsEntry } from "../types";

interface DistancePlotProps {
  codes: CodeEntry[];
  bounds: BoundsEntry[];
}

export function DistancePlot({ codes, bounds }: DistancePlotProps) {
  // Filter codes that have a computed distance
  const validCodes = codes.filter((c) => c.d !== null);

  // Build the scatter trace for our found codes
  const codesTrace: Plotly.Data = {
    x: validCodes.map((c) => c.n),
    y: validCodes.map((c) => c.d),
    mode: "markers" as const,
    type: "scatter" as const,
    name: "Found codes",
    text: validCodes.map(
      (c) =>
        `[[${c.n}, ${c.k}, ${c.d}]]<br>rate: ${c.rate.toFixed(3)}<br>max weight: ${c.max_weight}<br>strategy: ${c.strategy ?? "unknown"}`
    ),
    hoverinfo: "text" as const,
    marker: {
      color: "#3b82f6",
      size: 9,
      symbol: "circle",
      line: { color: "#60a5fa", width: 1 },
    },
  };

  // Group bounds by k value to draw connected lines
  const kValues = [...new Set(bounds.map((b) => b.k))].sort((a, b) => a - b);

  const lowerTraces: Plotly.Data[] = [];
  const upperTraces: Plotly.Data[] = [];

  for (const k of kValues) {
    const kBounds = bounds
      .filter((b) => b.k === k)
      .sort((a, b) => a.n - b.n);

    if (kBounds.length === 0) continue;

    lowerTraces.push({
      x: kBounds.map((b) => b.n),
      y: kBounds.map((b) => b.d_lower),
      mode: "lines+markers" as const,
      type: "scatter" as const,
      name: `Lower bound (k=${k})`,
      text: kBounds.map(
        (b) =>
          `Lower bound<br>n=${b.n}, k=${b.k}<br>d_lower=${b.d_lower}<br>optimal: ${b.optimal ? "yes" : "no"}`
      ),
      hoverinfo: "text" as const,
      marker: {
        color: "#10b981",
        size: 5,
        symbol: "triangle-up",
      },
      line: {
        color: "#10b981",
        width: 1.5,
        dash: "dot",
      },
      legendgroup: `lower-k${k}`,
    });

    upperTraces.push({
      x: kBounds.map((b) => b.n),
      y: kBounds.map((b) => b.d_upper),
      mode: "lines+markers" as const,
      type: "scatter" as const,
      name: `Upper bound (k=${k})`,
      text: kBounds.map(
        (b) =>
          `Upper bound<br>n=${b.n}, k=${b.k}<br>d_upper=${b.d_upper}<br>optimal: ${b.optimal ? "yes" : "no"}`
      ),
      hoverinfo: "text" as const,
      marker: {
        color: "#ef4444",
        size: 5,
        symbol: "triangle-down",
      },
      line: {
        color: "#ef4444",
        width: 1.5,
        dash: "dot",
      },
      legendgroup: `upper-k${k}`,
    });
  }

  const data = [...lowerTraces, ...upperTraces, codesTrace];

  const layout: Partial<Plotly.Layout> = {
    title: {
      text: "Code Distance vs Block Length",
      font: { color: "#e8edf5", size: 16 },
    },
    xaxis: {
      title: { text: "n (block length)", font: { color: "#8899b0" } },
      color: "#8899b0",
      gridcolor: "#1e2d3d",
      zerolinecolor: "#2a3a4f",
    },
    yaxis: {
      title: { text: "d (distance)", font: { color: "#8899b0" } },
      color: "#8899b0",
      gridcolor: "#1e2d3d",
      zerolinecolor: "#2a3a4f",
    },
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "#0a0f1a",
    font: { color: "#e8edf5" },
    legend: {
      bgcolor: "rgba(26,35,50,0.8)",
      bordercolor: "#1e2d3d",
      borderwidth: 1,
      font: { color: "#8899b0", size: 11 },
    },
    margin: { t: 48, r: 24, b: 56, l: 56 },
    hovermode: "closest" as const,
  };

  const config: Partial<Plotly.Config> = {
    displayModeBar: true,
    displaylogo: false,
    responsive: true,
    modeBarButtonsToRemove: [
      "lasso2d",
      "select2d",
      "sendDataToCloud",
      "toggleSpikelines",
    ] as Plotly.ModeBarDefaultButtons[],
  };

  if (validCodes.length === 0 && bounds.length === 0) {
    return (
      <div className="text-[var(--color-text-muted)] p-4 text-center">
        No data available for plotting.
      </div>
    );
  }

  return (
    <div className="px-4 pb-2">
      <div className="bg-[var(--color-bg-card)] border border-[var(--color-border-dim)] rounded-lg p-3">
        <Plot
          data={data}
          layout={layout}
          config={config}
          useResizeHandler
          style={{ width: "100%", height: "420px" }}
        />
      </div>
    </div>
  );
}
