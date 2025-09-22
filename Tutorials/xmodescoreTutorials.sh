#!/usr/bin/env bash
# Refactored: Menu-driven version (default: run ALL tutorials non-interactively)
# XModeScore / qbphenix / qmechanic / qbbuster demonstrations
#
# CHANGELOG (fix for unbound variable):
#   * Ensured script exits early if not executed by Bash.

#  // BEGIN COPYRIGHT
#  /***********************************************************************
#     Copyright (c) 2021 QuantumBio Inc.
#  ***********************************************************************/
#  // END COPYRIGHT

# Require Bash explicitly (defensive)
[ -n "${BASH_VERSION:-}" ] || { echo "ERROR: This script must be run with bash."; exit 1; }

set -u -o pipefail

: "${QBHOME:?ERROR: QBHOME is not set! Source /path/to/DivConSuite/etc/qbenv.sh first.}"

DIVCON_BIN="${QBHOME}/bin/qmechanic"
QBPHENIX_BIN="${QBHOME}/bin/qbphenix"
QBBUSTER_BIN="${QBHOME}/bin/qbbuster"

if [[ ! -x "${DIVCON_BIN}" ]]; then
  echo "ERROR: ${DIVCON_BIN} not executable."
  exit 1
fi
[[ -x "${QBPHENIX_BIN}" ]] || echo "WARNING: qbphenix not found (some tutorials will skip)."
[[ -x "${QBBUSTER_BIN}" ]] || echo "WARNING: qbbuster not found (Buster tutorials will skip)."

WORKDIR="${PWD}"
DATE_FMT="+%Y-%m-%d %H:%M:%S"

# Always declare the array to avoid unbound errors with set -u
ENGINE_DIVCON_ARGS=""
if [[ -n "${ENFORCE_DIVCON:-}" ]]; then
  ENGINE_DIVCON_ARGS="--protonation DivCon --engine DivCon --protonateTautomers DivCon"
fi

log() { printf "[%s] %s\n" "$(date "${DATE_FMT}")" "$*"; }

section() {
  echo
  echo "======================================================================"
  echo "$*"
  echo "======================================================================"
}

safe_cd_root() {
  cd "${WORKDIR}" || { echo "ERROR: Cannot cd back to ${WORKDIR}"; exit 1; }
}

clean_make_cd() {
  local d="$1"
  rm -rf "${WORKDIR:?}/${d}"
  mkdir -p "${WORKDIR}/${d}"
  cd "${WORKDIR}/${d}" || { echo "ERROR: Cannot enter ${WORKDIR}/${d}"; exit 1; }
}

fetch() {
  local url="$1"
  local f
  f="$(basename "${url}")"
  if [[ -f "${f}" ]]; then
    log "Already present: ${f}"
    return 0
  fi
  if command -v wget >/dev/null 2>&1; then
    wget -q "${url}" || { echo "ERROR: wget failed: ${url}"; return 1; }
  elif command -v curl >/dev/null 2>&1; then
    curl -fsSL -O "${url}" || { echo "ERROR: curl failed: ${url}"; return 1; }
  else
    echo "ERROR: Need wget or curl"
    return 1
  fi
}

# ---------------------------
# Tutorial functions
# ---------------------------

tutorial_1() {
  section "Tutorial #1: XModeScore (commented example, no execution)"
  safe_cd_root
  clean_make_cd "xmodeScore_4YJR"
  log "Original qbphenix example was commented out in source; nothing to run."
}

tutorial_2a() {
  section "Tutorial 2a: XModeScore 3HS4 user-provided"
  safe_cd_root
  clean_make_cd "xmodeScore_3HS4"
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/3HS4-tautomers/3HS4+H.pdb
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/3HS4-tautomers/3HS4.mtz
  if [[ -x "${QBPHENIX_BIN}" ]]; then
    "${QBPHENIX_BIN}" --pdbFile 3HS4+H.pdb --dataFile 3HS4.mtz \
      --XModeScore --protomers "-1..1" --mmMethod amberff14sb --qmMethod pm6 \
      --protonation skip --nproc 20 --dir xmodeData \
      --selection "chain A resname AZM resid 701" \
      ${ENGINE_DIVCON_ARGS} \
      || { echo "Tutorial 2a failed"; return 1; }
  else
    echo "SKIP: qbphenix not available."
  fi
}

tutorial_2b() {
  section "Tutorial 2b: XModeScore 3HS4 pre-generated tautomers"
  safe_cd_root
  clean_make_cd "xmodeScore_pregen3HS4"
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/3HS4-tautomers/3HS4+H.pdb
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/3HS4-tautomers/3HS4.mtz
  mkdir 3HS4-tautomers
  cd 3HS4-tautomers
  for f in 0_0_0_0_-1.pdb 1_1_0_0_-1.pdb 3_1_0_2_-1.pdb 4_0_0_0_0.pdb 6_1_0_1_0.pdb; do
    fetch "https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/3HS4-tautomers/$f"
  done
  cd ..
  if [[ -x "${QBPHENIX_BIN}" ]]; then
    "${QBPHENIX_BIN}" --pdbFile 3HS4+H.pdb --dataFile 3HS4.mtz \
      --XModeScore 3HS4-tautomers \
      --selection "chain A resname AZM resid 701" \
      --mmMethod amberff14sb --qmMethod pm6 \
      --nproc 20 --dir xmodeData \
      ${ENGINE_DIVCON_ARGS} \
      || { echo "Tutorial 2b failed"; return 1; }
  else
    echo "SKIP: qbphenix not available."
  fi
}

tutorial_3() {
  section "Tutorial #3: XModeScore using Buster (2BSM)"
  safe_cd_root
  clean_make_cd "xmodeScore_external_2BSM"
  if [[ -x "${QBBUSTER_BIN}" ]]; then
    # shellcheck disable=SC1091
    source /share/apps/GlobalPhasing/linux-x86_64/BUSTER_snapshot_20221121/setup.sh 2>/dev/null || \
      echo "WARNING: BUSTER setup script not found."
    export CSDHOME=/share/apps/CCDC/CSD_2022
    fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/2BSM+H.pdb
    fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/2BSM.mtz
    qbbuster --pdbFile 2BSM+H.pdb --sfFile 2BSM.mtz \
      --XModeScore --protomers "-1..1" --protonation skip --makeCIF divcon \
      --mmMethod amberff14sb --qmMethod pm6 --qmWeight 5.0 --ncycles 1 \
      --resname BSM --engine buster --np 8 -v 2 \
      || { echo "Tutorial 3 failed"; return 1; }
  else
    echo "SKIP: qbbuster not available."
  fi
}

tutorial_4() {
  section "Tutorial #4: XModeScore flip/chiral via Buster (4wq6)"
  safe_cd_root
  clean_make_cd "xmodeScore_4wq6"
  if [[ -x "${QBBUSTER_BIN}" ]]; then
    # shellcheck disable=SC1091
    source /share/apps/GlobalPhasing/linux-x86_64/BUSTER_snapshot_20221121/setup.sh 2>/dev/null || \
      echo "WARNING: BUSTER setup script not found."
    fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/4wq6+H.pdb
    fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/4wq6.mtz
    qbbuster --pdbFile 4wq6+H.pdb --sfFile 4wq6.mtz \
      --XModeScore --protomers "-1..1" --exploreFlip --exploreChiral \
      --protonation skip --makeCIF grade --mmMethod amberff14sb \
      --qmMethod pm6 --engine buster --np 8 -v 2 \
      || { echo "Tutorial 4 failed"; return 1; }
  else
    echo "SKIP: qbbuster not available."
  fi
}

tutorial_5() {
  section "Tutorial #5: Automated Docking + XModeScore (5C3K)"
  safe_cd_root
  clean_make_cd "xmodeScore_autoDock"
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/5C3K.mtz
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/5C3K-H_refine_001.pdb
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/5C3K-H_refine_001.mtz
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/4XF_A_402.pdb
  if [[ -x "${QBPHENIX_BIN}" ]]; then
    "${QBPHENIX_BIN}" \
      --pdbFile 5C3K-H_refine_001.pdb \
      --dataFile 5C3K.mtz \
      --densityFile 5C3K-H_refine_001.mtz \
      --ligandFile 4XF_A_402.pdb \
      --dock dockFolderResults \
      --protonation skip \
      --Xmodescore \
      --protomers "-1..1" \
      --mmMethod amberff14sb --qmMethod pm6 \
      --nproc 16 --dir runDock \
      ${ENGINE_DIVCON_ARGS} \
      || { echo "Tutorial 5 failed"; return 1; }
  else
    echo "SKIP: qbphenix not available."
  fi
}

tutorial_6() {
  section "Tutorial #6: Real-Space Ligand Statistics"
  safe_cd_root
  clean_make_cd "RealSpaceStatictics"
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/RealSpaceStats/4wq6+H_2_1_0_0_0_F2_C1_refined.pdb
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/RealSpaceStats/4wq6+H_2_1_0_0_0_F2_C1_refined.mtz
  "${DIVCON_BIN}" 4wq6+H_2_1_0_0_0_F2_C1_refined.pdb 4wq6+H_2_1_0_0_0_F2_C1_refined.mtz \
      --xstats "/A/3TQ/601//" \
      -v2 --np 2 \
      || { echo "Tutorial 6 failed"; return 1; }
}

tutorial_7() {
  section "Tutorial #7: XModeScore macrocycle states (8alx)"
  safe_cd_root
  clean_make_cd "xmodeScore_8alx"
  if [[ -x "${QBBUSTER_BIN}" ]]; then
    # shellcheck disable=SC1091
    source /share/apps/GlobalPhasing/linux-x86_64/BUSTER_snapshot_20221121/setup.sh 2>/dev/null || \
      echo "WARNING: BUSTER setup script not found."
    fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/8alx+H.pdb
    fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/8alx+H.mtz
    mkdir 8alx-tautomers_depot
    cd 8alx-tautomers_depot
    fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/8alx+H_state1.pdb
    fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/8alx+H_state2.pdb
    cd ..
    "${QBBUSTER_BIN}" \
      --pdbFile 8alx+H.pdb \
      --sfFile 8alx+H.mtz \
      --XmodeScore 8alx-tautomers_depot \
      --makecif divcon \
      --qmWeight 7.0 \
      --protonation skip \
      --region-radius 0.0 \
      --buffer-radius 0.0 \
      --mmMethod amberff14sb \
      --qmMethod pm6 \
      --engine buster \
      --np 8 -v 2 \
      || { echo "Tutorial 7 failed"; return 1; }
  else
    echo "SKIP: qbbuster not available."
  fi
}

# ---------------------------
# Menu & Dispatch
# ---------------------------

print_menu() {
  cat <<'EOF'
Interactive XModeScore Tutorial Menu (-i to enable):
  1  Tutorial #1  : XModeScore (no-op/example)
  2  Tutorial 2a  : XModeScore 3HS4 user-provided
  3  Tutorial 2b  : XModeScore 3HS4 pre-generated tautomers
  4  Tutorial #3  : XModeScore via Buster (2BSM)
  5  Tutorial #4  : XModeScore flip/chiral via Buster (4wq6)
  6  Tutorial #5  : Automatic docking + XModeScore (5C3K)
  7  Tutorial #6  : Real-Space Ligand Statistics
  8  Tutorial #7  : XModeScore macrocycle states (8alx)
  0 / A           : Run All
  Q               : Quit
EOF
}

run_all() {
  tutorial_1
  tutorial_2a
  tutorial_2b
  tutorial_3
  tutorial_4
  tutorial_5
  tutorial_6
  tutorial_7
}

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Default (no options): Run ALL tutorials non-interactively.

Options:
  -i    Interactive menu mode (choose specific tutorials)
  -l    List tutorials only
  -h    Help
EOF
}

dispatch() {
  case "$1" in
    1) tutorial_1 ;;
    2) tutorial_2a ;;
    3) tutorial_2b ;;
    4) tutorial_3 ;;
    5) tutorial_4 ;;
    6) tutorial_5 ;;
    7) tutorial_6 ;;
    8) tutorial_7 ;;
    0|A|a) run_all ;;
    Q|q) log "Quit requested"; exit 0 ;;
    *) echo "WARNING: Unknown selection: $1" ;;
  esac
}

main() {
  log "BEGIN XModeScore Tutorial Batch (DIVCON=${DIVCON_BIN})"

  local mode="all"
  while getopts ":ilh" opt; do
    case "${opt}" in
      i) mode="interactive" ;;
      l) print_menu; exit 0 ;;
      h) usage; exit 0 ;;
      *) usage; exit 1 ;;
    esac
  done
  shift $((OPTIND-1))

  if [[ "${mode}" == "all" ]]; then
    run_all
  else
    print_menu
    echo
    read -r -p "Select tutorials (ENTER=All): " sels
    if [[ -z "${sels}" ]]; then
      run_all
    else
      for s in ${sels}; do
        dispatch "${s}"
      done
    fi
  fi

  log "END XModeScore Tutorial Batch"
}

main "$@"