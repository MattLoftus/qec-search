// Core data types for the QEC search dashboard

export interface CodeEntry {
  type: "css";
  n: number;
  k: number;
  d: number | null;
  Hx: number[][];
  Hz: number[][];
  hash: string;
  max_weight: number;
  avg_weight: number;
  rate: number;
  generation?: number;
  strategy?: string;
  timestamp?: number;
  beats_known?: boolean;
}

export interface Campaign {
  name: string;
  strategy: string;
  total_tested: number;
  valid_codes: number;
  codes_with_k: number;
  generations: number;
  elapsed_sec: number;
  best_distance: Record<string, number>; // "n,k" -> d
  best_codes: CodeEntry[];
  timestamp: number;
}

export interface BoundsEntry {
  n: number;
  k: number;
  d_lower: number;
  d_upper: number;
  optimal: boolean;
}

export interface BenchmarkComparison {
  n: number;
  k: number;
  d_found: number;
  d_known_lower: number | null;
  d_known_upper: number | null;
  status: "IMPROVEMENT" | "MATCHES_BEST" | "BELOW_BEST" | "NO_REFERENCE";
  gap: number | null;
}

export interface BenchmarkResult {
  summary: {
    total: number;
    matches_best: number;
    improvements: number;
    below_best: number;
    no_reference: number;
  };
  comparisons: BenchmarkComparison[];
}

export interface StatusSummary {
  total_campaigns: number;
  total_best_codes: number;
  total_codes_tested: number;
  known_bounds_loaded: number;
  open_gaps: number;
}

export interface SimulationPoint {
  physical_error_rate: number;
  logical_error_rate: number;
  logical_errors: number;
  shots: number;
  time_sec: number;
}

export type TabId = "dashboard" | "campaigns" | "bounds" | "guide";
