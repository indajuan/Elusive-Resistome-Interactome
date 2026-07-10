#!/usr/bin/env bash
# Serves the ARG Pipeline Explorer front-end at http://localhost:8010
set -euo pipefail
cd "$(dirname "$0")/webapp"
python3 -m http.server 8010
