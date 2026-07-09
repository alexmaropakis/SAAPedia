#!/usr/bin/env bash
# Production/shared launcher (e.g. on Northeastern Explorer). Unlike run.sh, this
# binds to the network, uses a shared database file, and enables the write
# password. Configure via environment variables:
#
#   SAAP_DB_PATH         path to the shared SQLite file (REQUIRED here)
#   SAAP_WRITE_PASSWORD  password required for import/delete/DOI edits
#   HOST                 bind address (default 0.0.0.0 = all interfaces)
#   PORT                 port (default 8000)
#
# Example:
#   SAAP_DB_PATH=/work/mylab/saap/saap.db \
#   SAAP_WRITE_PASSWORD='choose-a-strong-password' \
#   PORT=8000 ./serve.sh
set -euo pipefail

cd "$(dirname "$0")/backend"

if [ ! -d ".venv" ]; then
  python3 -m venv .venv
fi
# shellcheck disable=SC1091
source .venv/bin/activate
pip install --quiet --upgrade pip
pip install --quiet -r requirements.txt

HOST="${HOST:-0.0.0.0}"
PORT="${PORT:-8000}"

: "${SAAP_DB_PATH:?Set SAAP_DB_PATH to the shared database file path (e.g. /work/yourlab/saap/saap.db)}"
mkdir -p "$(dirname "$SAAP_DB_PATH")"

if [ -z "${SAAP_WRITE_PASSWORD:-}" ]; then
  echo "WARNING: SAAP_WRITE_PASSWORD is not set — imports/deletes will be OPEN to" >&2
  echo "         anyone who can reach this server. Set it to protect writes." >&2
fi

echo "Serving SAAP Database on http://${HOST}:${PORT}"
echo "Shared database: ${SAAP_DB_PATH}"
echo "Write-protected: $([ -n "${SAAP_WRITE_PASSWORD:-}" ] && echo yes || echo NO)"
# One worker: SQLite serializes writes; WAL handles concurrent reads.
exec uvicorn app.main:app --host "$HOST" --port "$PORT" --workers 1
