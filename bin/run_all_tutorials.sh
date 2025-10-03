#!/usr/bin/env bash
#
# Run all tutorial shell scripts in the Tutorials directory in parallel,
# creating a separate working directory per tutorial and capturing output.
#
# Produces a final summary including which tutorials FAILED and which SUCCEEDED.
#
# Environment:
#   QBHOME      (required) Root of QuantumBio installation (must contain etc/qbenv.sh)
#   QB_NPROC    (optional) Explicit number of parallel jobs to use.
#
# Exit codes:
#   0 on completion (even if some scripts failed; adjust near end if desired)
#   1 on configuration / environment errors
#
set -euo pipefail

# -------- Colors (disabled if not a terminal) ---------------------------------
if [ -t 1 ]; then
  COLOR_GREEN=$'\033[32m'
  COLOR_RED=$'\033[31m'
  COLOR_YELLOW=$'\033[33m'
  COLOR_CYAN=$'\033[36m'
  COLOR_DIM=$'\033[2m'
  COLOR_RESET=$'\033[0m'
else
  COLOR_GREEN=""; COLOR_RED=""; COLOR_YELLOW=""; COLOR_CYAN=""; COLOR_DIM=""; COLOR_RESET=""
fi

log()  { printf '%s%s%s\n' "${COLOR_CYAN}" "$*" "${COLOR_RESET}"; }
warn() { printf '%s%s%s\n' "${COLOR_YELLOW}" "$*" "${COLOR_RESET}" >&2; }
err()  { printf '%s%s%s\n' "${COLOR_RED}" "$*" "${COLOR_RESET}" >&2; }

# -------- Require QBHOME ------------------------------------------------------
if [ -z "${QBHOME:-}" ]; then
  err "ERROR: QBHOME is not set. Please export QBHOME=/path/to/QuantumBio/root and re-run."
  exit 1
fi
if [ ! -f "$QBHOME/etc/qbenv.sh" ]; then
  err "ERROR: Cannot find $QBHOME/etc/qbenv.sh"
  exit 1
fi
# shellcheck source=/dev/null
source "$QBHOME/etc/qbenv.sh"

# -------- Resolve script directory (portable) ---------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TUTORIALS="$SCRIPT_DIR/../Tutorials"
export TUTORIALS
if [ ! -d "$TUTORIALS" ]; then
  err "ERROR: Tutorials directory not found at: $TUTORIALS"
  exit 1
fi

# -------- Determine parallelism -----------------------------------------------
detect_cpus() {
  if [ -n "${QB_NPROC:-}" ]; then
    echo "$QB_NPROC"; return
  fi
  if command -v nproc >/dev/null 2>&1; then
    nproc; return
  fi
  if command -v sysctl >/dev/null 2>&1; then
    sysctl -n hw.ncpu 2>/dev/null && return
  fi
  if command -v getconf >/dev/null 2>&1; then
    getconf _NPROCESSORS_ONLN 2>/dev/null && return
  fi
  echo 1
}
CPUS="$(detect_cpus || echo 1)"
HALF=$(( (CPUS + 1) / 2 ))
[ "$HALF" -lt 1 ] && HALF=1
NPROC="$HALF"
export NPROC
log "Detected CPUs: $CPUS  -> Using parallel jobs: $NPROC"

# -------- Gather tutorial scripts ---------------------------------------------
IFS=$'\n' read -r -d '' -a SCRIPTS < <(find "$TUTORIALS" -maxdepth 1 -type f -name '*.sh' -print && printf '\0' || true)
TOTAL=${#SCRIPTS[@]}
if [ "$TOTAL" -eq 0 ]; then
  warn "No tutorial scripts (*.sh) found in $TUTORIALS"
  exit 0
fi
export TOTAL
log "Found $TOTAL tutorial script(s)."

# -------- Counter + result tracking -------------------------------------------
COUNTER_FILE=$(mktemp -t qb_tutorial_counter.XXXXXX)
LOCK_FILE="${COUNTER_FILE}.lock"
RESULTS_FILE=$(mktemp -t qb_tutorial_results.XXXXXX)
RESULTS_LOCK="${RESULTS_FILE}.lock"

echo 0 > "$COUNTER_FILE"

cleanup() {
  rm -f "$COUNTER_FILE" "$LOCK_FILE" "$RESULTS_FILE" "$RESULTS_LOCK" 2>/dev/null || true
}
trap cleanup EXIT INT TERM

increment_counter() {
  while :; do
    if ( set -o noclobber; : > "$LOCK_FILE" ) 2>/dev/null; then
      local count
      count=$(<"$COUNTER_FILE")
      count=$((count + 1))
      echo "$count" > "$COUNTER_FILE"
      rm -f "$LOCK_FILE"
      printf '%s' "$count"
      return 0
    fi
    sleep 0.05
  done
}

append_result() {
  # Args: status (OK/FAIL), name, elapsed_secs, exit_code
  local status name elapsed code
  status="$1"; name="$2"; elapsed="$3"; code="$4"
  while :; do
    if ( set -o noclobber; : > "$RESULTS_LOCK" ) 2>/dev/null; then
      printf '%s\t%s\t%s\t%s\n' "$status" "$name" "$elapsed" "$code" >> "$RESULTS_FILE"
      rm -f "$RESULTS_LOCK"
      return 0
    fi
    sleep 0.05
  done
}

export COUNTER_FILE LOCK_FILE RESULTS_FILE RESULTS_LOCK
export -f increment_counter append_result

# -------- Worker function -----------------------------------------------------
worker() {
  script="$1"
  base=$(basename "$script" .sh)
  start_ts=$(date +%s)

  printf '%s[START]%s %s (%s)\n' "${COLOR_DIM}" "${COLOR_RESET}" "$base" "$(date +%H:%M:%S)"

  mkdir -p "$base"
  (
    cd "$base"
    bash "$script" > OUT.SCREEN 2>&1
  )
  status=$?
  end_ts=$(date +%s)
  elapsed=$(( end_ts - start_ts ))

  count=$(increment_counter)

  if [ $status -eq 0 ]; then
    printf '%s[DONE ]%s  [%s/%s] %s (%s) +%ss\n' \
      "${COLOR_GREEN}" "${COLOR_RESET}" "$count" "$TOTAL" "$base" "$(date +%H:%M:%S)" "$elapsed"
    append_result "OK" "$base" "$elapsed" "0"
  else
    printf '%s[ERROR]%s [%s/%s] %s (%s) (exit %d, %ss)\n' \
      "${COLOR_RED}" "${COLOR_RESET}" "$count" "$TOTAL" "$base" "$(date +%H:%M:%S)" "$status" "$elapsed"
    append_result "FAIL" "$base" "$elapsed" "$status"
  fi
  return 0
}
export -f worker

# -------- Parallel execution --------------------------------------------------
printf '%s\0' "${SCRIPTS[@]}" \
| xargs -0 -n1 -P"$NPROC" bash -c 'worker "$@"' _

# -------- Summary -------------------------------------------------------------
COMPLETED=$(<"$COUNTER_FILE")
log "All queued tutorials finished. Completed: $COMPLETED / $TOTAL"

if [ ! -s "$RESULTS_FILE" ]; then
  warn "No results recorded (unexpected)."
  exit 1
fi

# Parse results
OK_LIST=()
FAIL_LIST=()
while IFS=$'\t' read -r status name elapsed code; do
  case "$status" in
    OK)   OK_LIST+=("$name $elapsed"s) ;;
    FAIL) FAIL_LIST+=("$name (exit $code, ${elapsed}s)") ;;
  esac
done < "$RESULTS_FILE"

OK_COUNT=${#OK_LIST[@]}
FAIL_COUNT=${#FAIL_LIST[@]}

echo
echo "==================== SUMMARY ===================="
echo "Succeeded: $OK_COUNT"
echo "Failed:    $FAIL_COUNT"
echo "Total:     $TOTAL"
echo

if [ $FAIL_COUNT -gt 0 ]; then
  echo "Failed Tutorials:"
  # Sort for stable output
  printf '  %s\n' "${FAIL_LIST[@]}" | sort
  echo
fi

echo "Succeeded Tutorials:"
printf '  %s\n' "${OK_LIST[@]}" | sort
echo "================================================="
echo

# OPTIONAL: If you want a non-zero exit code when failures occurred, uncomment:
# if [ $FAIL_COUNT -gt 0 ]; then
#   exit 2
# fi

exit 0