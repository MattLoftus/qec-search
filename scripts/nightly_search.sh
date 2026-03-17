#!/bin/bash
# QEC Code Search — Nightly automation
# Runs weight-4 evolutionary search, BB sweep, and benchmarks
# Install: launchctl load ~/Library/LaunchAgents/com.qec-search.nightly.plist
#
# Resource constraint: ≤25% CPU/RAM (nice -n 15, reduced workloads)
# Total runtime: ~3.5 hours max

set -euo pipefail

PROJECT_DIR="$HOME/workspace/qec-search"
ENGINE_DIR="$PROJECT_DIR/engine"
VENV_PYTHON="$PROJECT_DIR/venv/bin/python3"
SYS_PYTHON="/usr/bin/python3"
LOG_DIR="$PROJECT_DIR/logs"
NTFY_TOPIC="qec-search-7k2m"

mkdir -p "$LOG_DIR"
LOGFILE="$LOG_DIR/nightly-$(date +%Y-%m-%d).log"

exec > >(tee -a "$LOGFILE") 2>&1
echo "=== QEC Nightly Search: $(date) ==="

# All compute runs at low priority (nice -n 15)
NICE="nice -n 15"
cd "$ENGINE_DIR"

# Phase 1: Weight-4 evolutionary search (most productive) — 1.5 hour cap
# This found [[20,2,6]] MATCHING known bounds
echo "Phase 1: Weight-4 evolutionary polynomial search..."
timeout 5400 $NICE "$VENV_PYTHON" evolve_polys.py --pop 40 --gens 200 --weight 4 --min-order 9 --max-order 30 2>&1 || true

# Phase 2: Weight-3 evolutionary search (different region) — 45 min cap
echo "Phase 2: Weight-3 evolutionary polynomial search..."
timeout 2700 $NICE "$VENV_PYTHON" evolve_polys.py --pop 30 --gens 100 --weight 3 --min-order 15 --max-order 30 2>&1 || true

# Phase 3: BB code sweep (qLDPC) — 45 min cap
echo "Phase 3: BB code polynomial sweep..."
timeout 2700 $NICE "$VENV_PYTHON" qldpc_search.py --mode bb --max-n 50 2>&1 || true

# Phase 4: Benchmark
echo "Phase 4: Benchmarking..."
"$SYS_PYTHON" evolver.py --benchmark 2>&1

# Phase 5: Status + notification
echo "Phase 5: Status..."
STATUS=$("$SYS_PYTHON" evolver.py --status 2>&1)
echo "$STATUS"

# Count improvements
IMPROVEMENTS=$(echo "$STATUS" | grep -c "IMPROVEMENT" || echo "0")
MATCHES=$(echo "$STATUS" | grep -c "MATCHES" || echo "0")
BEST=$(echo "$STATUS" | grep "All-time best codes:" | grep -o '[0-9]*' || echo "?")

# Send notification
SUMMARY="QEC nightly: ${BEST} best codes, ${MATCHES} matches, ${IMPROVEMENTS} improvements"
PRIORITY=3
if [ "$IMPROVEMENTS" -gt 0 ] 2>/dev/null; then
    PRIORITY=5
    SUMMARY="$SUMMARY"
fi

curl -s \
    -H "Title: QEC Search Nightly" \
    -H "Priority: $PRIORITY" \
    -d "$SUMMARY" \
    "https://ntfy.sh/$NTFY_TOPIC" > /dev/null 2>&1 || true

echo "=== Done: $(date) ==="
