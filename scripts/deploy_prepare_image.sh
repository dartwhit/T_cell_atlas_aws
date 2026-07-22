#!/usr/bin/env bash
# Prepare the production image tcell_atlas:latest for docker-compose.
# Reuses an already-built image (e.g. the PR image on this host) whose
# build-input fingerprint matches, else builds once. Runs as a standalone
# script so appleboy/ssh-action's line-splitting can't corrupt the logic.
set -euo pipefail

FP="${1:?fingerprint required}"
FORCE_REBUILD="${2:-false}"
IMAGE="tcell_atlas:latest"

echo "  Build fingerprint: $FP"

if [ "$FORCE_REBUILD" = "true" ]; then
  echo "  force_rebuild — full rebuild with a fresh base image"
  docker build --pull -t "$IMAGE" --label "atlas.build_fp=$FP" .
  exit 0
fi

# First image (incl. current latest) whose label matches this fingerprint.
REUSE_ID="$(docker images --filter "label=atlas.build_fp=$FP" --format '{{.ID}}' | head -n1 || true)"

if [ -n "$REUSE_ID" ]; then
  echo "  Reusing image $REUSE_ID (fingerprint match)"
  docker tag "$REUSE_ID" "$IMAGE"
else
  echo "  No matching image — building once (layer cache preserved)"
  docker build -t "$IMAGE" --label "atlas.build_fp=$FP" .
fi
