#!/usr/bin/env bash
set -euo pipefail

# --- Require QBHOME and source qbenv.sh ---
if [ -z "${QBHOME:-}" ]; then
  echo "ERROR: QBHOME is not set. Please export QBHOME=/path/to/QuantumBio/root and re-run." >&2
  exit 1
fi
# shellcheck source=/dev/null
source "$QBHOME/etc/qbenv.sh"

# --- Locate this script, then Tutorials directory ---
if command -v readlink >/dev/null 2>&1 && readlink -f / >/dev/null 2>&1; then
  SCRIPT_DIR="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")"
else
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
export TUTORIALS="$SCRIPT_DIR/../Tutorials"

# --- Parallelism: half the cores, min 1 ---
NPROC=$(( $(nproc) / 2 ))
if [ "$NPROC" -lt 1 ]; then
  NPROC=1
fi

# --- Count total .sh scripts in Tutorials ---
TOTAL=$(find "$TUTORIALS" -maxdepth 1 -type f -name '*.sh' \
        | wc -l | awk '{print $1}')
COUNTER_FILE=$(mktemp)
echo 0 > "$COUNTER_FILE"

export COUNTER_FILE TOTAL

# --- Run all tutorial scripts in parallel with progress output ---
find "$TUTORIALS" -maxdepth 1 -type f -name '*.sh'  \
  | xargs -n1 -P"$NPROC" -I{} bash -c '
      script="{}"
      base=$(basename "$script" .sh)

      echo "[START]  $base ($(date +%H:%M:%S))"
      mkdir -p "$base"
      (
        cd "$base" && bash "$script" > OUT.SCREEN 2>&1
      )
      status=$?

      # Update counter atomically
      {
        flock 200
        count=$(<"$COUNTER_FILE")
        count=$((count + 1))
        echo "$count" > "$COUNTER_FILE"
      } 200>"$COUNTER_FILE.lock"

      if [ $status -eq 0 ]; then
        echo "[DONE]   [$count/$TOTAL] $base ($(date +%H:%M:%S))"
      else
        echo "[ERROR]  [$count/$TOTAL] $base ($(date +%H:%M:%S)) (exit code $status)"
      fi
  ' _
