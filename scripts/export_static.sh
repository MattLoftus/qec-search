#!/bin/bash
# Export current search data as static JSON for web deployment.
# Run before deploying to Vercel.

set -euo pipefail

PROJECT_DIR="$HOME/workspace/qec-search"
ENGINE_DIR="$PROJECT_DIR/engine"
OUT_DIR="$PROJECT_DIR/web/public/data"
SYS_PYTHON="/usr/bin/python3"

mkdir -p "$OUT_DIR"

echo "Exporting static data for web deployment..."

cd "$ENGINE_DIR"

# Export all API endpoints as static JSON
"$SYS_PYTHON" -c "
import sys, json, os
sys.path.insert(0, '.')
from serve import APIHandler
import http.server

# Create a mock handler to reuse the data loading logic
class MockHandler:
    pass

h = APIHandler.__new__(APIHandler)

# Status
data = h._get_status()
with open('$OUT_DIR/status.json', 'w') as f:
    json.dump(data, f, default=str)
print(f'  status.json: {json.dumps(data, default=str)[:80]}...')

# Best codes
data = h._get_best()
with open('$OUT_DIR/best.json', 'w') as f:
    json.dump(data, f, default=str)
print(f'  best.json: {len(data)} codes')

# Benchmark
data = h._get_benchmark()
with open('$OUT_DIR/benchmark.json', 'w') as f:
    json.dump(data, f, default=str)
print(f'  benchmark.json: {json.dumps(data.get(\"summary\", {}), default=str)[:80]}')

# Bounds
data = h._get_bounds()
with open('$OUT_DIR/bounds.json', 'w') as f:
    json.dump(data, f, default=str)
print(f'  bounds.json: {len(data)} entries')

# Campaigns
data = h._get_campaigns()
with open('$OUT_DIR/campaigns.json', 'w') as f:
    json.dump(data, f, default=str)
print(f'  campaigns.json: {len(data)} campaigns')
"

echo "Static data exported to $OUT_DIR"
ls -la "$OUT_DIR"
