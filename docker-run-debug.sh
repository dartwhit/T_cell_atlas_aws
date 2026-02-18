#!/bin/bash

# Simple single-container debug setup
# Usage: ./docker-run-debug.sh

set -e

CONTAINER_NAME="atlas-debug"
DATA_PATH="/home/ubuntu/atlas/data"
PORT=3838

echo "🐳 Starting single Shiny container for debugging..."
echo ""

# Stop and remove any existing debug container
if docker ps -a --format '{{.Names}}' | grep -q "^${CONTAINER_NAME}$"; then
    echo "Stopping existing container: $CONTAINER_NAME"
    docker stop $CONTAINER_NAME || true
    docker rm $CONTAINER_NAME || true
fi

# Build fresh image
echo "📦 Building Docker image..."
docker build -t atlas-debug:latest .

# Run single container with interactive logging
echo "🚀 Starting container..."
docker run \
    --name $CONTAINER_NAME \
    --rm \
    -it \
    -p $PORT:3838 \
    -v $DATA_PATH:/srv/shiny-server/atlas/data:ro \
    -v /tmp/atlas-debug-logs:/var/log/shiny-server \
    -e SHINY_LOG_LEVEL=DEBUG \
    atlas-debug:latest

echo ""
echo "Container stopped."
