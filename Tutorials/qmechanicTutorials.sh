#!/usr/bin/env bash
# Menu-driven refactor (patterned after mtscoreTutorials.sh)
# Original content: qmechanicTutorials.sh (single point, decomposition, optimization, ONIOM, protonation, macrocycle, covalent examples)
#
# Changes:
#  - Added robust argument parsing (-i interactive, -l list, -h help)
#  - Added safety helpers (fetch, clean directories, logging with timestamps)
#  - Converted legacy select/menu into numbered dispatcher similar to mtscoreTutorials.sh
#  - Kept original qmechanic command lines functionally identical
#  - Removed stray 'grep' placeholder lines that had no arguments (no functional effect)
#  - Replaced duplicate nested directory creation (macrocycle / bound_ligand) with a cleaner layout
#
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

# (Original script had set -e commented out; we keep similar leniency)
set -u -o pipefail

: "${QBHOME:?ERROR: QBHOME is not set! Must export QBHOME before running.}"

DIVCON_BIN="${QBHOME}/bin/qmechanic"
if [[ ! -x "${DIVCON_BIN}" ]]; then
  echo "ERROR: ${DIVCON_BIN} is not an executable! Adjust QBHOME or install qmechanic."
  exit 1
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
    echo "ERROR: Cannot return to root work directory: ${WORKDIR}"
    exit 1
  }
}

clean_make_cd() {
  local dir="$1"
  rm -rf "${WORKDIR:?}/${dir}"
  mkdir -p "${WORKDIR}/${dir}"
  cd "${WORKDIR}/${dir}" || {
    echo "ERROR: Cannot cd into ${WORKDIR}/${dir}"
    exit 1
  }
}

fetch() {
  local url="$1"
  local file
  file="$(basename "${url}")"
  if [[ -f "${file}" ]]; then
    log "Already downloaded: ${file}"
    return 0
  fi
  if command -v wget >/dev/null 2>&1; then
    wget -q "${url}" && return 0
  elif command -v curl >/dev/null 2>&1; then
    curl -fsSL -O "${url}" && return 0
  else
    echo "ERROR: Neither wget nor curl is available."
    return 1
  fi
  echo "ERROR: Failed to download ${url}"
  return 1
}

###############################################################################
# Tutorial Functions
###############################################################################

tutorial_1() {
  section "Tutorial 1: Single Point Calculation"
  safe_cd_root
  clean_make_cd "SinglePointCalculation"
  fetch http://downloads.quantumbioinc.com/media/tutorials/cli/1TOW-H.pdb
  mv 1TOW-H.pdb 1TOW.pdb
  "${DIVCON_BIN}" 1TOW.pdb --np 8 -h pm6 -v 2
}

tutorial_2() {
  section "Tutorial 2: Interaction Energy Decomposition"
  safe_cd_root
  clean_make_cd "InteractionEnergyDecomposition"
  fetch http://downloads.quantumbioinc.com/media/tutorials/cli/1LRI-addH.pdb
  "${DIVCON_BIN}" 1LRI-addH.pdb -h pm6 --np 8 -v 2 --decompose "/A/CLR/99//" --dc
}

tutorial_3() {
  section "Tutorial 3: Structure Optimization"
  safe_cd_root
  clean_make_cd "StructureOptimization"
  fetch http://downloads.quantumbioinc.com/media/tutorials/cli/1LRI-addH.pdb
  "${DIVCON_BIN}" 1LRI-addH.pdb -h pm6 --np 8 -v 2 --opt all 3 0.01 --symmetry -p pdb
}

tutorial_4() {
  section "Tutorial 4: Active Site Structure Optimization"
  safe_cd_root
  clean_make_cd "ActiveSiteStructureOptimization"
  fetch http://downloads.quantumbioinc.com/media/tutorials/cli/1TOW-H.pdb
  "${DIVCON_BIN}" 1TOW-H.pdb --opt /*/CRZ/*// 3.0 0.0 --np 8 -h pm6 -v 2 -p pdb
}

tutorial_5() {
  section "Tutorial 5: ONIOM (mixed QM/MM) Simulations"
  safe_cd_root
  clean_make_cd "ONIOM_Simulations"
  fetch http://downloads.quantumbioinc.com/media/tutorials/cli/1LRI-addH.pdb
  "${DIVCON_BIN}" 1LRI-addH.pdb -h pm6 amberff14sb --opt 25 0.01 --qm-region /A/CLR/99// 3.0 0.0 --np 8 -v 1 -p pdb
}

tutorial_6() {
  section "Tutorial 6: Protonation (Prepare)"
  safe_cd_root
  clean_make_cd "Protonation"
  # Original used bare PDB ID 4EK4 relying on internal fetch; replicate directly
  "${DIVCON_BIN}" 4EK4 --prepare --np 8 -v 2 -p pdb -h amberff14sb
}

tutorial_7() {
  section "Tutorial 7: Macrocycle Protonation (JSON bond definition)"
  safe_cd_root
  clean_make_cd "macrocycle"
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/cli/3p72_macrocycle_bond.json
  # Original command (kept same flags/order; writing outputs 3p72+H.pdb 3p72.mtz)
  "${DIVCON_BIN}" 3p72 --standards 3p72_macrocycle_bond.json -h amberff14sb --prepare all --np 4 -v2 -p 3p72+H.pdb 3p72.mtz -O 2>&1
}

tutorial_8() {
  section "Tutorial 8: Single Point Covalently Bound Ligand (JSON bond definition)"
  safe_cd_root
  clean_make_cd "bound_ligand"
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/cli/5y41+H.pdb
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/cli/5y41_covalent_bond.json
  "${DIVCON_BIN}" 5y41+H.pdb --standards 5y41_covalent_bond.json -h amberff14sb --np 4 -v2 -O 2>&1
}

tutorial_9() {
  section "Tutorial 9: Template Driven Prepare"
  safe_cd_root
  clean_make_cd "templated_prepare"
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/cli/2VMF_1_model_1_relaxed-dssp.pdb
  "${DIVCON_BIN}" 2vmf --template 2VMF_1_model_1_relaxed-dssp.pdb  --prepare all -h amberff14sb --np 2 -v 2 -O 2>&1
}

###############################################################################
# Menu / Dispatch
###############################################################################

print_menu() {
  cat <<'EOF'
Interactive qmechanic Tutorial Menu (-i to enable):
  1  Single Point Calculation
  2  Interaction Energy Decomposition
  3  Structure Optimization
  4  Active Site Structure Optimization
  5  ONIOM (mixed QM/MM) Simulations
  6  Protonation (Prepare)
  7  Macrocycle Protonation (JSON bond def)
  8  Covalently Bound Ligand (JSON bond def)
  9  Template Driven Prepare
  0 / A  Run All
  Q      Quit
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
  tutorial_8
  tutorial_9
}

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Default (no options): run ALL tutorials sequentially.

Options:
  -i    Interactive menu mode
  -l    List tutorials only (no execution)
  -h    Show this help

Examples:
  $(basename "$0")
  $(basename "$0") -i
  $(basename "$0") -l
EOF
}

dispatch_number() {
  case "$1" in
    1) tutorial_1 ;;
    2) tutorial_2 ;;
    3) tutorial_3 ;;
    4) tutorial_4 ;;
    5) tutorial_5 ;;
    6) tutorial_6 ;;
    7) tutorial_7 ;;
    8) tutorial_8 ;;
    9) tutorial_9 ;;
    0|A|a) run_all ;;
    Q|q) log "User requested quit."; exit 0 ;;
    *) echo "WARNING: Unknown selection: $1" ;;
  esac
}

main() {
  log "BEGIN Tutorial Test using ${DIVCON_BIN}"

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
    read -r -p "Select tutorials (ENTER=All): " selections
    if [[ -z "${selections}" ]]; then
      run_all
    else
      for sel in ${selections}; do
        dispatch_number "${sel}"
      done
    fi
  fi

  log "END Tutorial Test using ${DIVCON_BIN}"
}

main "$@"