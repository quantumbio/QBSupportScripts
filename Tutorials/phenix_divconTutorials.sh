#!/usr/bin/env bash
# Refactored: Menu-driven version (default: run ALL tutorials non-interactively)
# Requirement: Preserve ORIGINAL command lines EXACTLY (including truncated [...])
# Source: Original phenix_divconTutorials.sh
#
# Tutorials:
#  1  Naked PHENIX refine (2WOR)
#  2a qbphenix MOE single ligand (3IX1)
#  2b qbphenix ReadySet single ligand (3IX1)
#  3  qbphenix MOE all ligands (3IX1)
#  4  qbphenix MOE single ligand (1LRI)
#  5  ONIOM QM/MM refinement + Clash (1NAV)
#  6  qbphenix MOE covalent ligand (3NCK)
#
# NOTE: Lines truncated in the original repository with [...] are kept verbatim.
#       Replace the [...] portions with the full original arguments if you have them.
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

set -u -o pipefail

: "${QBHOME:?ERROR: QBHOME is not set! You MUST source /path/to/DivConSuite/etc/qbenv.sh THEN call this script}"
: "${PHENIX:?ERROR: PHENIX is not set! You MUST source /path/to/phenix-1.XX.XXX/phenix_env.sh THEN call this script}"

DIVCON_BIN="${QBHOME}/bin/qmechanic"
if [[ ! -x "${DIVCON_BIN}" ]]; then
  echo "ERROR: ${DIVCON_BIN} is not an executable! Set this path to the qmechanic executable"
  exit 1
fi

WORKDIR="$PWD"

# Match original logic for optional DivCon enforcement
if [[ -v ENFORCE_DIVCON ]]; then
  ENGINE_DIVCON=" --protonation DivCon --engine DivCon"
else
  ENGINE_DIVCON=""
fi

DATE_FMT="+%Y-%m-%d %H:%M:%S"
log() { printf "[%s] %s\n" "$(date "${DATE_FMT}")" "$*"; }

section() {
  echo
  echo "======================================================================"
  echo "$*"
  echo "======================================================================"
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
# Tutorial functions (command lines preserved)
###############################################################################

tutorial_1() {
  section "Tutorial #1 (Running Phenix/DivCon without the qbphenix wrapper): 2WOR"
  local PDBID=2WOR
  clean_make_cd "NakedPHENIX"
  fetch "https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data//XModeScore/${PDBID}.pdb"
  fetch "https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data//XModeScore/${PDBID}-sf.cif"
  ${PHENIX}/build/bin/phenix.ready_set ${PDBID}.pdb add_h_to_water=true
  # Original truncated refine line:
  ${PHENIX}/build/bin/phenix.refine ${PDBID}.updated.pdb ${PDBID}-sf.cif output.write_geo_file=False output.write_eff_file=False output.write_def_file=False refinement.refine.strategy=individual_sites+individual_adp main.number_of_macro_cycles=2 qblib=True qblib_method=pm6 qblib_region_selection="/A/2AN/1098/" qblib_region_radius=3.0 qblib_buffer_radius=2.5 qblib_np=4
}

tutorial_2a() {
  section "Tutorial #2a (qbphenix + MOE single ligand): 3IX1"
  clean_make_cd "qbphenix_moe"
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data//XModeScore/3ix1.pdb
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data//XModeScore/3ix1-sf.cif
  $QBHOME/bin/qbphenix --dataFile 3ix1-sf.cif --pdbFile 3ix1.pdb --selection "chain A resname NFM resid 401"  --phenixOptions "main.number_of_macro_cycles=2" --protonation MOE --qmMethod pm6 --region-radius 3.0 --buffer-radius 2.5 --Nproc 4 --v 1 $ENGINE_DIVCON
}

tutorial_2b() {
  section "Tutorial #2b (qbphenix + phenix.ready_set single ligand): 3IX1"
  clean_make_cd "qbphenix_readyset"
  $QBHOME/bin/qbphenix --pdbID 3ix1 --selection "chain A resname NFM resid 401" --phenixOptions "main.number_of_macro_cycles=2" --protonation ReadySet  --qmMethod pm6 --region-radius 3.0 --buffer-radius 2.5  --Nproc 4
}

tutorial_3() {
  section "Tutorial #3 (qbphenix + MOE all ligands): 3IX1"
  clean_make_cd "qbphenix_moe_all"
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data//XModeScore/3ix1.pdb
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data//XModeScore/3ix1-sf.cif
  $QBHOME/bin/qbphenix --dataFile 3ix1-sf.cif --pdbFile 3ix1.pdb --selection "resname NFM" --phenixOptions "main.number_of_macro_cycles=1" --protonation MOE  --qmMethod pm6 --region-radius 3.0 --buffer-radius 2.5  --Nproc 4
}

tutorial_4() {
  section "Tutorial #4 (qbphenix + MOE single ligand): 1LRI"
  clean_make_cd "qbphenix_moe_1lri"
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data//XModeScore/1lri.pdb
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data//XModeScore/1lri-sf.cif
  $QBHOME/bin/qbphenix --dataFile 1lri-sf.cif --pdbFile 1lri.pdb --selection "chain A resname CLR resid 99"  --phenixOptions "main.number_of_macro_cycles=2" --protonation MOE --qmMethod pm6 --mmMethod amberff14sb --region-radius 3.0 --buffer-radius 2.5 --Nproc 4 $ENGINE_DIVCON
}

tutorial_5() {
  section "Tutorial #5: ONIOM (QM/MM) Based X-ray Refinement & Clash Score (1NAV)"
  clean_make_cd "Protonation_1NAV"
  # Original command with $ENGINE_DIVCON and redirect:
  $QBHOME/bin/qbphenix --pdbID 1NAV --mmMethod amberff14sb --qmMethod pm6 --selection "resname IH5" --np 4 --region-radius 3.0 --buffer-radius 0.0 --protonation MOE $ENGINE_DIVCON >& OUT.screen
}

tutorial_6() {
  section "Tutorial #6 (qbphenix MOE covalent ligand): 3NCK"
  clean_make_cd "Protonation_3NCK"
  $QBHOME/bin/qbphenix --pdbID 3NCK --selection "resname NFF" --phenixOptions "main.number_of_macro_cycles=2" --protonation divcon  --qmMethod pm6 --mmMethod amberff14sb --region-radius 3.0 --Nproc 4  --scriptName run $ENGINE_DIVCON
  # Original script invoked a generated run script:
  if [[ -x ./run ]]; then
    ./run || echo "WARNING: ./run returned non-zero status"
  else
    # Still mimic original behavior (attempt) – keep silent if absent
    :
  fi
}

###############################################################################
# Menu / Dispatch
###############################################################################

print_menu() {
  cat <<'EOF'
PHENIX/DivCon Tutorial Menu (-i for interactive):
  1  Tutorial #1  : Naked PHENIX refine (2WOR)
  2  Tutorial #2a : qbphenix MOE single ligand (3IX1)
  3  Tutorial #2b : qbphenix ReadySet single ligand (3IX1)
  4  Tutorial #3  : qbphenix MOE all ligands (3IX1)
  5  Tutorial #4  : qbphenix MOE single ligand (1LRI)
  6  Tutorial #5  : ONIOM QM/MM refinement + Clash (1NAV)
  7  Tutorial #6  : qbphenix MOE covalent ligand (3NCK)
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
}

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Default (no options): Run ALL tutorials sequentially.

Options:
  -i    Interactive menu mode
  -l    List tutorials only (no execution)
  -h    Help

Examples:
  $(basename "$0")
  $(basename "$0") -i
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
    0|A|a) run_all ;;
    Q|q) log "Quit requested"; exit 0 ;;
    *) echo "WARNING: Unknown selection: $1" ;;
  esac
}

main() {
  local startDate
  startDate="$(date)"
  echo "BEGIN PHENIX/DivCon Tutorial Test at ${startDate} using ${DIVCON_BIN}"

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

  local endDate
  endDate="$(date)"
  echo "END Tutorial Test at ${endDate} using ${DIVCON_BIN}"
}

main "$@"