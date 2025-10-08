#!/usr/bin/env bash
# Refactored: Menu-driven version (default: run ALL tutorials non-interactively)
# Original script: DivCon Xray tutorials (refinement, XModeScore, docking)

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

: "${QBHOME:?ERROR: QBHOME is not set! Source /path/to/DivConSuite/etc/qbenv.sh first.}"

DIVCON_BIN="${QBHOME}/bin/qmechanic"
QBDIVCON_BIN="${QBHOME}/bin/qbdivcon"

if [[ ! -x "${DIVCON_BIN}" ]]; then
  echo "ERROR: ${DIVCON_BIN} is not an executable! Adjust QBHOME."
  exit 1
fi
[[ -x "${QBDIVCON_BIN}" ]] || echo "WARNING: qbdivcon not found; XModeScore qbdivcon tutorials may be skipped."

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
    echo "ERROR: Need wget or curl installed."
    return 1
  fi
}

# ------------------------------------------------------------------------------
# Tutorial Functions
# NOTE: The original source we pulled had truncated lines with [...].
#       Those truncated command segments are preserved with comments; please
#       replace with the full original arguments as necessary.
# ------------------------------------------------------------------------------

tutorial_1() {
  section "Tutorial #1: All-Atom refinement with ONIOM region (protonated PDBid:3e34)"
  safe_cd_root
  local dir="3e34-All_Atoms"
  clean_make_cd "${dir}"
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/3e34+H.pdb
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/3e34.mtz

  "${DIVCON_BIN}" 3e34+H.pdb 3e34.mtz \
      --opt all 50 0.01 \
      --qm-region /B/ED1/1003// 3.0 0 \
      -h pm6 amberff14sb \
      --np 4 -O \
      -p 3e34_refined.pdb 3e34_refined.mtz \
      || { echo "Tutorial #1 failed"; return 1; }
}

tutorial_2() {
  section "Tutorial #2: All-Atom MM refinement with custom MTZ labels (PDBid:4jp4)"
  safe_cd_root
  # Original script reused 3e34-All_Atoms name; we keep original behavior but
  # recommend changing to 4jp4-All_Atoms to avoid overwriting Tutorial #1.
  local dir="3e34-All_Atoms"  # (Original reused name) Consider: 4jp4-All_Atoms
  clean_make_cd "${dir}"
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/4jp4+H.pdb
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/4jp4.mtz

  "${DIVCON_BIN}" 4jp4+H.pdb "4jp4.mtz#/*/*[I-obs(+),SIGI-obs(+),I-obs(-),SIGI-obs(-),R-free-flags]" \
      --opt all 25 0.01 \
      -h pm6 amberff14sb \
      --np 4 -O \
      -p 4jp4_refined.pdb 4jp4_refined.mtz \
      || { echo "Tutorial #2 failed"; return 1; }
}

tutorial_3() {
  section "Tutorial #3: All-Atom refinement with ONIOM region (download by PDB id 1lri)"
  safe_cd_root
  local dir="1lri-allAtom"
  clean_make_cd "${dir}"

  "${DIVCON_BIN}" 1lri \
      --prepare \
      -h amberff14sb \
      --np 4 -O \
      -p 1lri+H.pdb \
      || { echo "Tutorial #3 orotonation failed"; return 1; }

  "${DIVCON_BIN}" 1lri+H.pdb 1lri-sf.cif \
      --opt all 50 0.01 \
      --qm-region /A/CLR/99// 3.0 0 \
      -h pm6 amberff14sb \
      --np 4 -O \
      -p 1lri_refined.pdb 1lri_refined.mtz \
      || { echo "Tutorial #3 oniom failed"; return 1; }
}

tutorial_4() {
  section "Tutorial #4: XModeScore executed on protonated PDBid:4b72"
  safe_cd_root
  local dir="xmodeScore_4b72"
  clean_make_cd "${dir}"
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/4b72+H.pdb
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/4b72.mtz

  if [[ -x "${QBDIVCON_BIN}" ]]; then
    # Truncated original line preserved:
    "${QBDIVCON_BIN}" \
      --pdbFile 4b72+H.pdb \
      --sfFile 4b72.mtz \
      --XmodeScore \
      --protonation Skip \
      --protomers "-1..1" \
      --engine divcon \
      --protonateTautomers divcon \
      --qmMethod pm6 \
      --mmMethod amberff14sb \
      --np 4 -v 2 -O \
      || { echo "Tutorial #4 failed"; return 1; }
  else
    echo "SKIP: qbdivcon not available."
  fi
}

tutorial_5() {
  section "Tutorial #5: XModeScore with Dock (protonated PDBid:1bzc)"
  safe_cd_root
  local dir="xmodeScore_1bzc"
  clean_make_cd "${dir}"
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/1bzc+H.pdb
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/1bzc.mtz

  if [[ -x "${QBDIVCON_BIN}" ]]; then
    export MT_ALLOW_MTDOCK=1
    # Truncated original line preserved:
    "${QBDIVCON_BIN}" \
      --pdbFile 1bzc+H.pdb \
      --sfFile 1bzc.mtz \
      --XmodeScore \
      --protonation Skip \
      --protomers "0" \
      --exploreFlip \
      --exploreChiral \
      --exploreDocking \
      --engine divcon \
      --protonateTautomers divcon \
      --qmMethod pm6 \
      --mmMethod amberff14sb \
      --np 4 -v 2 -O \
      || { echo "Tutorial #5 failed"; export MT_ALLOW_MTDOCK=; return 1; }
    unset MT_ALLOW_MTDOCK
  else
    echo "SKIP: qbdivcon not available."
  fi
}

tutorial_6() {
  section "Tutorial #6: Multi-Blob Docking PDBid:4o9s"
  safe_cd_root
  local dir="multiBlob_4o9s"
  clean_make_cd "${dir}"
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/Dock/4o9s+HnL5.pdb
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/Dock/4o9s_lig.mol2
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/Dock/4o9s.mtz
  export QB_FIX_METAL_CHARGES=1  # (Original script exported but did not unset)

  "${DIVCON_BIN}" 4o9s+HnL5.pdb 4o9s.mtz \
      --ligand blobs 4o9s_lig.mol2 \
      -h amberff14sb \
      --np 2 -v2 \
      --dock opt rigid \
      -O \
      > test2.log 2>&1 \
      || { echo "Tutorial #6 failed"; return 1; }
}

# ------------------------------------------------------------------------------
# Menu / Dispatch
# ------------------------------------------------------------------------------

print_menu() {
  cat <<'EOF'
Interactive Xray Tutorial Menu (-i to enable):
  1  Tutorial #1 : 3e34 all-atom ONIOM refinement
  2  Tutorial #2 : 4jp4 all-atom MM refinement (custom MTZ labels)
  3  Tutorial #3 : 1lri all-atom refinement with ONIOM region
  4  Tutorial #4 : 4b72 XModeScore (qbdivcon)
  5  Tutorial #5 : 1bzc XModeScore + Dock (qbdivcon)
  6  Tutorial #6 : 4o9s multi-blob docking
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
  $(basename "$0")       # Run all tutorials
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
    0|A|a) run_all ;;
    Q|q) log "Quit requested"; exit 0 ;;
    *) echo "WARNING: Unknown selection: $1" ;;
  esac
}

main() {
  log "BEGIN DivCon Xray Tutorial Batch (DIVCON_BIN=${DIVCON_BIN})"

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

  log "END DivCon Xray Tutorial Batch"
}

main "$@"