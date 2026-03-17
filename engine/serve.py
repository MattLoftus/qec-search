"""
REST API server for the QEC search dashboard.

Endpoints:
    GET /api/status     — search status summary
    GET /api/campaigns  — list all campaigns
    GET /api/campaign/<name> — full campaign data
    GET /api/best       — all-time best codes
    GET /api/benchmark  — benchmark comparison
    GET /api/bounds     — known bounds table
"""

import http.server
import json
import os
import sys
import argparse

RESULTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
CAMPAIGNS_DIR = os.path.join(RESULTS_DIR, "campaigns")
BEST_CODES_PATH = os.path.join(RESULTS_DIR, "best_codes.json")
BENCHMARK_PATH = os.path.join(RESULTS_DIR, "benchmark.json")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from known_codes import KNOWN_BOUNDS


class APIHandler(http.server.BaseHTTPRequestHandler):
    def do_GET(self):
        # CORS
        self.send_response(200)
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Content-Type", "application/json")
        self.end_headers()

        path = self.path.rstrip("/")
        try:
            if path == "/api/status":
                self._send_json(self._get_status())
            elif path == "/api/campaigns":
                self._send_json(self._get_campaigns())
            elif path.startswith("/api/campaign/"):
                name = path.split("/api/campaign/")[1]
                self._send_json(self._get_campaign(name))
            elif path == "/api/best":
                self._send_json(self._get_best())
            elif path == "/api/benchmark":
                self._send_json(self._get_benchmark())
            elif path == "/api/bounds":
                self._send_json(self._get_bounds())
            else:
                self._send_json({"error": "not found", "endpoints": [
                    "/api/status", "/api/campaigns", "/api/campaign/<name>",
                    "/api/best", "/api/benchmark", "/api/bounds",
                ]})
        except Exception as e:
            self._send_json({"error": str(e)})

    def do_OPTIONS(self):
        self.send_response(200)
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "GET")
        self.send_header("Access-Control-Allow-Headers", "*")
        self.end_headers()

    def _send_json(self, data):
        self.wfile.write(json.dumps(data, default=str).encode())

    def _get_status(self):
        campaigns = self._get_campaigns()
        best = self._get_best()
        return {
            "total_campaigns": len(campaigns),
            "total_best_codes": len(best),
            "total_codes_tested": sum(c.get("total_tested", 0) for c in campaigns),
            "known_bounds_loaded": len(KNOWN_BOUNDS),
            "open_gaps": len([(n, k) for (n, k), (dl, du) in KNOWN_BOUNDS.items() if dl < du]),
        }

    def _get_campaigns(self):
        if not os.path.exists(CAMPAIGNS_DIR):
            return []
        campaigns = []
        for f in sorted(os.listdir(CAMPAIGNS_DIR)):
            if f.endswith('.json'):
                with open(os.path.join(CAMPAIGNS_DIR, f)) as fh:
                    campaigns.append(json.load(fh))
        return campaigns

    def _get_campaign(self, name):
        path = os.path.join(CAMPAIGNS_DIR, f"{name}.json")
        if not os.path.exists(path):
            return {"error": f"Campaign '{name}' not found"}
        with open(path) as f:
            return json.load(f)

    def _get_best(self):
        if not os.path.exists(BEST_CODES_PATH):
            return []
        with open(BEST_CODES_PATH) as f:
            return json.load(f)

    def _get_benchmark(self):
        if not os.path.exists(BENCHMARK_PATH):
            return {"summary": {}, "comparisons": []}
        with open(BENCHMARK_PATH) as f:
            return json.load(f)

    def _get_bounds(self):
        return [
            {"n": n, "k": k, "d_lower": dl, "d_upper": du, "optimal": dl == du}
            for (n, k), (dl, du) in sorted(KNOWN_BOUNDS.items())
        ]

    def log_message(self, format, *args):
        # Quieter logging
        pass


def main():
    parser = argparse.ArgumentParser(description="QEC Search API Server")
    parser.add_argument("--port", type=int, default=8789)
    args = parser.parse_args()

    server = http.server.HTTPServer(("", args.port), APIHandler)
    print(f"QEC Search API server: http://localhost:{args.port}")
    print(f"Endpoints: /api/status, /api/campaigns, /api/best, /api/benchmark, /api/bounds")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nShutting down.")


if __name__ == "__main__":
    main()
