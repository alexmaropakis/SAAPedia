#!/usr/bin/env bash

# Launch SAAPedia

cd "$(dirname "$0")/backend"

if [ ! -d ".venv" ]; then
  echo "Creating virtual environment..."
  python3 -m venv .venv
fi
# shellcheck disable=SC1091
source .venv/bin/activate

echo "venv initialized, installing dependencies..."
pip install --quiet --upgrade pip
pip install --quiet -r requirements.txt

# Local copy: show dataset names and allow import/annotate/delete.
# A public deployment leaves this unset (tissue/species only, read-only).
export SAAP_PRIVATE="${SAAP_PRIVATE:-1}"

HOST="${HOST:-127.0.0.1}"
PORT="${PORT:-8000}"
echo ""
echo "All ready!"
echo "SAAP Database running at: http://${HOST}:${PORT}"
echo "Press Ctrl+C to stop."
echo ""
exec uvicorn app.main:app --host "$HOST" --port "$PORT" --reload
