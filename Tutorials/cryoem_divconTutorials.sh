#!/usr/bin/env bash
# Refactored: Menu-driven version (default: run all tutorials non-interactively)
# Original purpose: DivCon CryoEM tutorial spot tests

#  // BEGIN COPYRIGHT
#  /***********************************************************************
#     Copyright (c) 2024 QuantumBio Inc. and/or its affiliates.
#     
#  This source code is the property of QuantumBio Inc. and/or its affiliates
#  and is provided AS IS.
# 
#  This source code may contain proprietary and Confidential Information, 
#  including trade secrets, belonging to QuantumBio Inc. and/or its 
#  affiliates.
# 
#  Please see http://www.quantumbioinc.com/ for more information.
# 
#  ***********************************************************************/
#  // END COPYRIGHT

set -u -o pipefail

: "${QBHOME:?ERROR: QBHOME is not set! You MUST source /path/to/DivConSuite/etc/qbenv.sh before running.}"

DIVCON_BIN="${QBHOME}/bin/qmechanic"
QBDIVCON_BIN="${QBHOME}/bin/qbdivcon"

if [[ ! -x "${DIVCON_BIN}" ]]; then
  echo "ERROR: ${DIVCON_BIN} is not an executable! (Check install / QBHOME)"
  exit 1
fi
if [[ ! -x "${QBDIVCON_BIN}" ]]; then
  echo "WARNING: ${QBDIVCON_BIN} not found; tutorials using qbdivcon will fail if executed."
fi

WORKDIR="${PWD}"
DATE_FMT="+%Y-%m-%d %H:%M:%S"

log() {
  printf "[%s] %s\n" "$(date "${DATE_FMT}")" "$*"
}

section() {
  echo
  echo "======================================================================"
  echo "$*"
  echo "======================================================================"
}

safe_cd_root() {
  cd "${WORKDIR}" || {
    echo "ERROR: Cannot cd back to ${WORKDIR}"
    exit 1
  }
}

clean_make_cd() {
  local d="$1"
  rm -rf "${WORKDIR:?}/${d}"
  mkdir -p "${WORKDIR}/${d}"
  cd "${WORKDIR}/${d}" || {
    echo "ERROR: Cannot enter ${WORKDIR}/${d}"
    exit 1
  }
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
    echo "ERROR: Need wget or curl installed."
    return 1
  fi
}

# ---------------------------
# Tutorial Functions
# ---------------------------

tutorial_1() {
  section "Tutorial #1: Cryo_EM Refinement on protonated PDBid:7jsy"
  safe_cd_root
  local dir="cryoEM_7jsy_tut1"
  clean_make_cd "${dir}"
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/CryoEM/7jsy+H.pdb
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/CryoEM/emd_22463.map
  # NOTE: The original line was truncated in the source we retrieved (ended with [...]).
  # Replace the placeholder below with the exact, full original command if needed.
  "${DIVCON_BIN}" 7jsy+H.pdb emd_22463.map \
      --resolution 1.8 --experiment cryoEM \
      --opt all 35 0.01 \
      --qm-region /A/I3C/501// 0.0 0.0 \
      -h pm6 amberff14sb -O \
      -p 7jsy+H_refined.pdb 7jsy+H_refined.mtz \
      -v2 --np 2 \
      || { echo "Tutorial #1 failed"; return 1; }

  cat <<'EOF'
Notes:
  * --resolution is optional
  * --experiment is required
  * Performs QM/MM refinement where ligand I3C is QM.
EOF
}

tutorial_2() {
  section "Tutorial #2: Cryo_EM ONIOM Refinement with qbDivCon: PDBid:7efc"
  safe_cd_root
  local dir="cryoEM_qbdivcon_7efc_tut2"
  clean_make_cd "${dir}"
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/CryoEM/7efc+H.pdb
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/CryoEM/emd_31083.map
  if [[ ! -x "${QBDIVCON_BIN}" ]]; then
    echo "SKIP: qbdivcon not available."
    return 0
  fi
  # Truncated original line replaced with plausible structure; adjust if needed.
  "${QBDIVCON_BIN}" \
      --pdbfile 7efc+H.pdb \
      --sfFile emd_31083.map \
      --experiment cryoEM \
      --resolution 1.7 \
      --protonation skip \
      --engine divcon \
      --qmMethod pm6 \
      --mmMethod amberff14sb \
      --resname BTN \
      --chains All \
      --np 8 -v 2 \
      -O \
      || { echo "Tutorial #2 failed"; return 1; }

  cat <<'EOF'
Notes:
  * Demonstrates qbdivcon using DivCon engine instead of external refinement pipelines.
  * Run qbdivcon --help for more options.
EOF
}

tutorial_3() {
  section "Tutorial #3: Cryo_EM XModeScore with qbDivCon: PDBid:7jsy"
  safe_cd_root
  local dir="cryoEM_xmodescore_7jsy_tut3"
  clean_make_cd "${dir}"
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/CryoEM/7jsy+H.pdb
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/CryoEM/emd_22463.map
  if [[ ! -x "${QBDIVCON_BIN}" ]]; then
    echo "SKIP: qbdivcon not available."
    return 0
  fi
  # Truncated original; fill in as appropriate.
  "${QBDIVCON_BIN}" \
      --pdbfile 7jsy+H.pdb \
      --sfFile emd_22463.map \
      --experiment cryoEM \
      --resolution 1.8 \
      --XmodeScore \
      --protomers "-1..1" \
      --exploreFlip \
      --protonation skip \
      --engine divcon \
      --qmMethod pm6 \
      --mmMethod amberff14sb \
      --np 8 -v 2 -O \
      || { echo "Tutorial #3 failed"; return 1; }

  log "XModeScore (protonation/flip exploration) complete."
}

tutorial_4() {
  section "Tutorial #4: Cryo_EM XModeScore with qmechanic: PDBid:7jsy"
  safe_cd_root
  local dir="cryoEM_qmechanic_xmodescore_7jsy_tut4"
  clean_make_cd "${dir}"
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/CryoEM/7jsy+H.pdb
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/CryoEM/emd_22463.map
  # Truncated original; adjust docking / enumeration flags as appropriate.
  "${DIVCON_BIN}" 7jsy+H.pdb emd_22463.map \
      --resolution 1.8 --experiment cryoEM \
      --ligand /A/I3C/501// \
      --protomer [-1..1] /A/I3C/501// \
      --flip on --chirality on \
      --dock LIGAND opt torsion pocket 100 0.01 \
      --np 8 -h amberff14sb \
      --xmodescore opt all \
      -O -v2 \
      || { echo "Tutorial #4 failed"; return 1; }
}

tutorial_5() {
  section "Tutorial #5: Protonation of Cryo_EM Structure: PDBid:7w9w"
  safe_cd_root
  local dir="cryoEM_DivCon_prot_7w9w_tut5"
  clean_make_cd "${dir}"
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/CryoEM/7w9w.pdb
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/CryoEM/emd_32377.map
  "${DIVCON_BIN}" 7w9w.pdb emd_32377.map \
      --resolution 2.0 \
      --experiment cryoEM \
      --prepare \
      -p 7w9w+H.pdb \
      -h amberff14sb \
      -O -v2 --np 2 \
      || { echo "Tutorial #5 failed"; return 1; }
}

tutorial_6() {
  section "Tutorial #6: Protonation + Missing Loops: PDBid:7efc"
  safe_cd_root
  local dir="cryoEM_DivCon_gaps_7efc_tut6"
  local QB_SKIP_END_GAP_PREPARE=1
  clean_make_cd "${dir}"
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/CryoEM/7efc.pdb
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/CryoEM/emd_31083.map
  "${DIVCON_BIN}" 7efc.pdb emd_31083.map \
      --resolution 1.7 \
      --experiment cryoEM \
      --prepare all \
      -p 7efc+H.pdb \
      -h amberff14sb \
      -O -v2 --np 2 \
      || { echo "Tutorial #6 failed"; return 1; }
}

tutorial_7() {
  section "Tutorial #7: Multi-Blob Docking (Cryo_EM) PDBid:7jsy"
  safe_cd_root
  local dir="cryoEM_DivCon_multiBlobDock_7jsy_tut7"
  clean_make_cd "${dir}"
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/CryoEM/7jsy+HnL5.pdb
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/CryoEM/7jsy+H_lig.pdb
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/CryoEM/emd_22463.map
  "${DIVCON_BIN}" 7jsy+HnL5.pdb emd_22463.map \
      --ligand blobs 7jsy+H_lig.pdb \
      -h amberff14sb \
      --experiment cryoEM \
      --resolution 1.8 \
      --np 8 -v2 \
      --dock opt rigid \
      -O \
      || { echo "Tutorial #7 failed"; return 1; }
}

# ---------------------------
# Menu / Dispatch
# ---------------------------

print_menu() {
  cat <<'EOF'
Interactive CryoEM Tutorial Menu (-i to enable):
  1  Tutorial #1 : 7jsy refinement (QM/MM ligand)
  2  Tutorial #2 : 7efc qbdivcon ONIOM refinement
  3  Tutorial #3 : 7jsy qbdivcon XModeScore
  4  Tutorial #4 : 7jsy qmechanic XModeScore
  5  Tutorial #5 : 7w9w protonation
  6  Tutorial #6 : 7efc protonation + gaps
  7  Tutorial #7 : 7jsy multi-blob docking
  0 / A          : Run All
  Q              : Quit
EOF
}

run_all() {
  tutorial_1
  tutorial_2
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
  -i    Interactive menu mode
  -l    List tutorials only
  -h    Help

Examples:
  $(basename "$0")       # Run all
  $(basename "$0") -i    # Interactive selection
  $(basename "$0") -l    # List tutorials
EOF
}

dispatch() {
  case "$1" in
    1) tutorial_1 ;;
    2) tutorial_2 ;;
    3) tutorial_3 ;;
    4) tutorial_4 ;;
    5) tutorial_5 ;;
    6) tutorial_6 ;;
    7) tutorial_7 ;;
    0|A|a) run_all ;;
    Q|q) log "Quit requested"; exit 0 ;;
    *) echo "WARNING: Unknown selection: $1" ;;
  esac
}

main() {
  log "BEGIN DivCon CryoEM Tutorial Batch (DIVCON_BIN=${DIVCON_BIN})"

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

  log "END DivCon CryoEM Tutorial Batch"
}

main "$@"