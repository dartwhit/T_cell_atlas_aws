#!/usr/bin/env bash
# Make sure nginx is running, then wait for all three Shiny replicas to be healthy
# and probe the proxied HTTPS endpoint. Exits non-zero (with diagnostics) if the
# backends never become ready or the site isn't served through nginx.
# Used by .github/workflows/deploy.yml and scripts/start_stack.sh. Kept as a
# standalone script so appleboy/ssh-action's line-splitting can't corrupt it.
# Usage (from the repo root): bash scripts/verify_stack.sh
set -euo pipefail
cd "$(dirname "$0")/.."

# Healthcheck window in docker-compose.yml is start_period 300s + 10 retries x 15s = 450s.
HEALTH_TIMEOUT="${HEALTH_TIMEOUT:-480}"
SITE_HOST="${SITE_HOST:-ssccellatlas.dartmouth.edu}"
APPS="tcell_atlas_app_1 tcell_atlas_app_2 tcell_atlas_app_3"

diagnostics() {
  echo "Container status:"
  docker-compose ps || true
  echo "Nginx logs:"
  docker-compose logs --tail=30 nginx || true
  echo "App 1 logs:"
  docker-compose logs --tail=30 shiny-app-1 || true
}

# Start nginx first so the site is reachable (502 per backend) while apps warm up.
echo "🌐 Ensuring nginx is running..."
docker-compose up -d nginx

echo "⏳ Waiting for app containers to be healthy (max ${HEALTH_TIMEOUT}s)..."
START=$(date +%s)
while true; do
  HEALTHY=0
  for c in $APPS; do
    status=$(docker inspect -f '{{if .State.Health}}{{.State.Health.Status}}{{else}}none{{end}}' "$c" 2>/dev/null || echo missing)
    if [ "$status" = "healthy" ]; then HEALTHY=$((HEALTHY + 1)); fi
  done
  ELAPSED=$(( $(date +%s) - START ))
  if [ "$HEALTHY" -eq 3 ]; then
    echo "✓ All 3 app containers healthy after ${ELAPSED}s"
    break
  fi
  if [ "$ELAPSED" -ge "$HEALTH_TIMEOUT" ]; then
    echo "❌ Only $HEALTHY/3 app containers healthy after ${ELAPSED}s"
    diagnostics
    exit 1
  fi
  echo "  [${ELAPSED}s] Healthy: $HEALTHY/3"
  sleep 5
done

# Probe through nginx over HTTPS: port 80 only redirects (301), which would
# "pass" curl -f without ever reaching Shiny. --resolve keeps SNI/Host correct;
# -k because this checks backend reachability, not certificate validity.
echo "🩺 Probing https://${SITE_HOST}/atlas/ via nginx on localhost..."
for i in 1 2 3 4 5; do
  if curl -sSfk -o /dev/null --max-time 15 --resolve "${SITE_HOST}:443:127.0.0.1" "https://${SITE_HOST}/atlas/"; then
    echo "✓ HTTPS check through nginx passed"
    exit 0
  fi
  echo "  Attempt $i/5 failed, retrying..."
  sleep 3
done

echo "❌ HTTPS check through nginx failed after 5 attempts"
diagnostics
exit 1
