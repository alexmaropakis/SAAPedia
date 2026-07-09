#!/usr/bin/env bash
# Launch the SAAP Database app (backend serves the frontend too).
set -euo pipefail

cd "$(dirname "$0")/backend"

if [ ! -d ".venv" ]; then
  echo "Creating virtual environment..."
  python3 -m venv .venv
fi
# shellcheck disable=SC1091
source .venv/bin/activate

echo "Installing dependencies..."
pip install --quiet --upgrade pip
pip install --quiet -r requirements.txt

HOST="${HOST:-127.0.0.1}"
PORT="${PORT:-8000}"
echo ""
echo "SAAP Database running at: http://${HOST}:${PORT}"
echo "Press Ctrl+C to stop."
echo ""
exec uvicorn app.main:app --host "$HOST" --port "$PORT" --reload
