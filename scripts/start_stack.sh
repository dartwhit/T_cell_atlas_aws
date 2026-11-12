#!/usr/bin/env bash
# Start (or restart) the production stack and make sure nginx is actually running.
# Usage (on EC2, from the repo root): bash scripts/start_stack.sh
set -uo pipefail
cd "$(dirname "$0")/.."
docker-compose up -d --remove-orphans || echo "⚠️ compose up returned non-zero (likely a slow app healthcheck); continuing"
docker-compose up -d nginx
docker-compose ps
curl -sS -o /dev/null -w "HTTP %{http_code} from http://localhost/atlas/\n" http://localhost/atlas/ || true
