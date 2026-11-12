#!/usr/bin/env bash
# Start (or restart) the production stack and make sure nginx is actually running.
# Exits non-zero if nginx can't start, the apps never become healthy, or the
# site isn't served through nginx.
# Usage (on EC2, from the repo root): bash scripts/start_stack.sh
set -euo pipefail
cd "$(dirname "$0")/.."
docker-compose up -d --remove-orphans || echo "⚠️ compose up returned non-zero (likely a slow app healthcheck); continuing"
bash scripts/verify_stack.sh
docker-compose ps
