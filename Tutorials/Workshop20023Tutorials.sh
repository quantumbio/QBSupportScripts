#!/usr/bin/env bash
# Refactored: Menu-driven version (default: run ALL tutorials non-interactively)
# Request: Preserve original command lines EXACTLY (including any truncated [...] segments).
#
# Original script title: Workshop 2003 Tutorials (Workshop20023Tutorials.sh)
#
# NOTE:
#   Lines in the original source that contained truncated arguments (ending with [...])
#   are kept verbatim. Replace the [...] with full option lists from your canonical
#   version if you need fully functioning runs.

#  // BEGIN COPYRIGHT
#  /***********************************************************************
#     Copyright (c) 2021 QuantumBio Inc. and/or its affiliates.
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

: "${QBHOME:?ERROR: QBHOME is not set! You MUST source /path/to/DivConSuite/etc/qbenv.sh first.}"

DIVCON_BIN="${QBHOME}/bin/qmechanic"
if [[ ! -x "${DIVCON_BIN}" ]]; then
  echo "ERROR: ${DIVCON_BIN} is not an executable! Set this path correctly."
  exit 1
fi

# Default executable (overridden below if GRID_MARKETS is set)
qbExec="qbdivcon"
cloud=""

WORKDIR="${PWD}"
DATE_FMT="+%Y-%m-%d %H:%M:%S"

if [[ -n "${GRID_MARKETS:-}" ]]; then
  qbExec="qbdivcon"
  cloud=" --cloud gridmarkets"
fi

log() { printf "[%s] %s\n" "$(date "${DATE_FMT}")" "$*"; }

section() {
  echo
  echo "======================================================================"
  echo "$*"
  echo "======================================================================"
}

safe_cd_root() {
  cd "${WORKDIR}" || { echo "ERROR: Cannot return to ${WORKDIR}"; exit 1; }
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

###############################################################################
# Tutorial Functions (commands preserved EXACTLY)
###############################################################################

tutorial_1() {
  section "Tutorial #1: ONIOM Refinement of PDBid:1LRI with Protonation"
  safe_cd_root
  local dir="Tutorial1_1LRI"
  clean_make_cd "${dir}"
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/Workshop_2023/1LRI.pdb
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/Workshop_2023/1LRI.mtz
  # Original command (truncated as in source)
  $QBHOME/bin/$qbExec --pdbfile 1LRI.pdb --sfFile 1LRI.mtz --protonation divcon  --engine divcon --qmMethod pm6 --mmMethod amberff14sb --resname CLR --np 4 --region-radius 3.0 --nSmallCycles 25 --ncycles 2 --dir qm_results $cloud || { echo "Tutorial #1 failed"; return 1; }
}

tutorial_1a() {
  section "Tutorial #1a: ONIOM PHENIX Refinement of PDBid:1LRI"
  safe_cd_root
  local dir="Tutorial1_1LRI_PHENIX"
  clean_make_cd "${dir}"
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/Workshop_2023/1LRI.pdb
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/Workshop_2023/1LRI.mtz
  $QBHOME/bin/qbphenix --pdbfile 1LRI.pdb --sfFile 1LRI.mtz --protonation divcon --makeCIF divcon --qmMethod pm6 --mmMethod amberff14sb --resname CLR --np 4 --region-radius 3.0 --nSmallCycles 25 --ncycles 2 --dir qm_results || { echo "Tutorial #1a failed"; return 1; }
}

tutorial_2() {
  section "Tutorial #2: ONIOM Refinement with FAD PDBid:1SIQ"
  safe_cd_root
  local dir="Tutorial2_1SIQ"
  clean_make_cd "${dir}"
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/Workshop_2023/1SIQ+H.pdb
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/Workshop_2023/1SIQ.mtz
  $QBHOME/bin/$qbExec --pdbfile 1SIQ+H.pdb --sfFile 1SIQ.mtz --protonation skip  --engine divcon --qmMethod pm6 --mmMethod amberff14sb --resname FAD --chain A --resid 399 --np 4 --region-radius 3.0 --nSmallCycles 40 --ncycles 2 --dir qm_results $cloud || { echo "Tutorial #2 failed"; return 1; }
}

tutorial_3() {
  section "Tutorial #3: ONIOM XModeScore on AZM PDBid:3HS4"
  safe_cd_root
  local dir="Tutorial3_AZM_XmodeScore"
  clean_make_cd "${dir}"
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/Workshop_2023/3HS4+H.pdb
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/Workshop_2023/3HS4.mtz
  $QBHOME/bin/$qbExec --pdbfile 3HS4+H.pdb --sfFile 3HS4.mtz --Xmodescore --protomers -1..0  --protonation skip  --engine divcon --qmMethod pm6 --mmMethod amberff14sb --resname AZM --chain A --resid 701 --np 10 --region-radius 3.0 --nSmallCycles 35 --dir XmodeScore_results $cloud || { echo "Tutorial #3 failed"; return 1; }
}

tutorial_4() {
  section "Tutorial #4: XModeScore on user-provided AZM tautomer files"
  safe_cd_root
  local dir="Tutorial4_AZM_XmodeScore_ExternalFiles"
  clean_make_cd "${dir}"
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/Workshop_2023/3HS4+H.pdb
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/Workshop_2023/3HS4.mtz
  mkdir 3HS4-tautomers
  cd 3HS4-tautomers
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/Workshop_2023/3HS4+H_0_0_0_0_-1_C1_refined.pdb
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/Workshop_2023/3HS4+H_1_1_0_0_-1_C1_refined.pdb
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/Workshop_2023/3HS4+H_3_1_0_2_-1_C1_refined.pdb
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/Workshop_2023/3HS4+H_4_0_0_0_0_C1_refined.pdb
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/Workshop_2023/3HS4+H_6_1_0_1_0_C1_refined.pdb
  cd ..
  $QBHOME/bin/$qbExec --pdbFile 3HS4+H.pdb --dataFile  3HS4.mtz --XModeScore 3HS4-tautomers --protonation skip  --protonateTautomers skip --engine divcon --mmMethod amberff14sb --qmMethod pm6 --selection "chain A resname AZM resid 701" --nproc 10 --region-radius 3.0 --nSmallCycles 35 --dir XmodeScore_results $cloud || { echo "Tutorial #4 failed"; return 1; }
}

tutorial_5() {
  section "Tutorial #5: Automated MOE-based Docking + XModeScore (qbphenix) REQUIRES MOE"
  safe_cd_root
  local dir="Tutorial5_MOE_AutoDock_Xmode"
  clean_make_cd "${dir}"
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/5C3K.mtz
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/5C3K-H_refine_001.pdb
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/5C3K-H_refine_001.mtz
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/4XF_A_402.pdb
  $QBHOME/bin/qbphenix --pdbFile 5C3K-H_refine_001.pdb --dataFile 5C3K.mtz --densityFile 5C3K-H_refine_001.mtz --ligandFile 4XF_A_402.pdb --dock dockFolderResults --protonation skip --Xmodescore --protonateTautomers MOE --mmMethod amberff14sb --qmMethod pm6 --region-radius 3.0 --np 22 --dir Dock_Scored_results || { echo "Tutorial #5 failed"; return 1; }
}

tutorial_6() {
  section "Tutorial #6: XModeScore flip + chiral exploration (4wq6)"
  safe_cd_root
  local dir="Tutorial6_xmodeScore_4wq6"
  clean_make_cd "${dir}"
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/4wq6+H.pdb
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/4wq6.mtz
  $QBHOME/bin/$qbExec --pdbFile 4wq6+H.pdb --sfFile 4wq6.mtz --XModeScore --protomers "-1..1" --exploreFlip --exploreChiral --protonation skip  --mmMethod amberff14sb --qmMethod pm6 --engine divcon --nSmallCycles 40 --resname 3TQ --chain A --resid 601 --np 12 --region-radius 3.0 --dir XmodeScore_results $cloud || { echo "Tutorial #6 failed"; return 1; }
}

###############################################################################
# Menu / Dispatch
###############################################################################

print_menu() {
  cat <<'EOF'
Workshop 2003 Tutorials Menu (-i to enable interactive mode):
  1  Tutorial #1  : ONIOM Refinement 1LRI (protonation divcon)
  2  Tutorial #1a : ONIOM PHENIX Refinement 1LRI
  3  Tutorial #2  : ONIOM Refinement 1SIQ (FAD)
  4  Tutorial #3  : ONIOM XModeScore AZM (3HS4)
  5  Tutorial #4  : XModeScore AZM user tautomer files
  6  Tutorial #5  : MOE-based Docking + XModeScore (5C3K) (Requires MOE)
  7  Tutorial #6  : XModeScore flip/chiral exploration (4wq6)
  0 / A           : Run All
  Q               : Quit
EOF
}

run_all() {
  tutorial_1
  tutorial_1a
  tutorial_2
  tutorial_3
  tutorial_4
  tutorial_5
  tutorial_6
}

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Default (no options) : Run ALL tutorials sequentially.
Options:
  -i    Interactive (choose specific tutorials)
  -l    List tutorials only (no execution)
  -h    Help

Example:
  $(basename "$0") -i
EOF
}

dispatch() {
  case "$1" in
    1) tutorial_1 ;;
    2) tutorial_1a ;;
    3) tutorial_2 ;;
    4) tutorial_3 ;;
    5) tutorial_4 ;;
    6) tutorial_5 ;;
    7) tutorial_6 ;;
    0|A|a) run_all ;;
    Q|q) log "Quit requested"; exit 0 ;;
    *) echo "WARNING: Unknown selection: $1" ;;
  esac
}

main() {
  log "BEGIN Workshop 2003 Tutorials (DIVCON_BIN=${DIVCON_BIN}, qbExec=${qbExec}${cloud})"

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

  log "END Workshop 2003 Tutorials"
}

main "$@"