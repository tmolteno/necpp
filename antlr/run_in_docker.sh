#!/bin/bash
# Run a command inside the necpp-antlr4 Docker container.
# The project root is volume-mounted at /workspace.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
docker run --rm -i --user "$(id -u):$(id -g)" \
  -v "$PROJECT_DIR:/workspace" -w /workspace/antlr \
  necpp-antlr4 "$@"
