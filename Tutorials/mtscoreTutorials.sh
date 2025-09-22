#!/usr/bin/env bash
# set -e  # (Left intentionally commented—individual commands handle errors)
set -u -o pipefail

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

#############################################
# Configuration / Globals
#############################################
: "${QBHOME:?ERROR: QBHOME is not set! Must export QBHOME before running.}"

DIVCON_BIN="${QBHOME}/bin/qmechanic"
if [[ ! -x "${DIVCON_BIN}" ]]; then
  echo "ERROR: ${DIVCON_BIN} is not an executable! Adjust QBHOME or install qmechanic."
  exit 1
fi

MOEBATCH_BIN="${MOEBATCH_BIN:-$(command -v moebatch || true)}"
HAVE_MOEBATCH="yes"
if [[ -z "${MOEBATCH_BIN}" ]]; then
  HAVE_MOEBATCH="no"
fi

WORKDIR="${PWD}"
DATE_FMT="+%Y-%m-%d %H:%M:%S"

#############################################
# Utility Functions
#############################################
log() {
  printf "[%s] %s\n" "$(date "${DATE_FMT}")" "$*"
}

section() {
  echo
  echo "======================================================================"
  echo "$*"
  echo "======================================================================"
}

need_moebatch() {
  if [[ "${HAVE_MOEBATCH}" != "yes" ]]; then
    log "SKIP (moebatch not found): $*"
    return 1
  fi
  return 0
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
  local attempt
  for attempt in 1 2; do
    if command -v wget >/dev/null 2>&1; then
      wget -q "${url}" && return 0
    elif command -v curl >/dev/null 2>&1; then
      curl -fsSL -O "${url}" && return 0
    else
      echo "ERROR: Neither wget nor curl is available."
      return 1
    fi
    log "Download failed (attempt ${attempt}) for ${url}; retrying..."
    sleep 2
  done
  echo "ERROR: Failed to download ${url}"
  return 1
}

#############################################
# Tutorial Functions
#############################################
tutorial_1a() {
  section "First A Tutorial: Command Line: MTConfSearch (SMILES)"
  safe_cd_root
  clean_make_cd "MTCS-SMILES"
  "${DIVCON_BIN}" "c1ccccc1CCCC" -p sdf -h amberff14sb -v2
}

tutorial_1b() {
  section "First B Tutorial: Command Line: MTConfSearch (MOL2)"
  safe_cd_root
  local tutorFolder="MTCS-MOL2"
  clean_make_cd "${tutorFolder}"
  fetch "http://downloads.quantumbioinc.com/media/tutorials/MT/4k18_ligand.mol2"
  "${DIVCON_BIN}" --ligand 4k18_ligand.mol2 --mtcs -O -h garf
}

tutorial_1c() {
  section "First C Tutorial: Command Line: MTConfSearch Energy Cutoff"
  safe_cd_root
  local tutorFolder="MTCS-MOL2_CUTOFF"
  clean_make_cd "${tutorFolder}"
  fetch "http://downloads.quantumbioinc.com/media/tutorials/MT/4k18_ligand.mol2"
  "${DIVCON_BIN}" --ligand 4k18_ligand.mol2 --mtcs 10kcal -O -h garf
}

tutorial_1d() {
  section "First D Tutorial: Command Line: MTConfSearch Energy Cutoff With Opt"
  safe_cd_root
  local tutorFolder="MTCS-MOL2_CUTOFF_OPT"
  clean_make_cd "${tutorFolder}"
  fetch "http://downloads.quantumbioinc.com/media/tutorials/MT/4k18_ligand.mol2"
  "${DIVCON_BIN}" --ligand 4k18_ligand.mol2 --mtcs 10kcal input opt torsion 100 0.1 -O -h amberff14sb
}

tutorial_2() {
  section "Second Tutorial: Command Line: MTScore (Endstate)"
  safe_cd_root
  clean_make_cd "MTScoreES"
  fetch "http://downloads.quantumbioinc.com/media/tutorials/MT/4w7t_protein.pdb"
  fetch "http://downloads.quantumbioinc.com/media/tutorials/MT/4w7t_ligand.mol2"
  "${DIVCON_BIN}" 4w7t_protein.pdb --ligand 4w7t_ligand.mol2 -h amberff14sb --mtscore endstate
}

tutorial_3() {
  section "Third Tutorial: 3Step MTCS+MTScoreE (Ensemble) with External Docker"
  safe_cd_root
  clean_make_cd "3step"
  fetch "http://downloads.quantumbioinc.com/media/tutorials/MT/Bace_030215_CAT_4p.tar.gz"
  tar xvf Bace_030215_CAT_4p.tar.gz >/dev/null
  cp Bace_030215_CAT_4p/Bace_030215.pdb .
  cp Bace_030215_CAT_4p/Bace_030215_CAT_4p.mol2 .
  "${DIVCON_BIN}" Bace_030215.pdb --ligand Bace_030215_CAT_4p.mol2 --mtcs input -h amberff14sb --np 2 -v 2
  if need_moebatch "Third Tutorial docking step"; then
    "${MOEBATCH_BIN}" -licwait -run "${QBHOME}/svl/run/qbDockPair.svl" \
      -rec Bace_030215.pdb -lig Bace_030215_CAT_4p.mol2 \
      -conf Bace_030215_CAT_4p_conf.sdf -protonate -delwat
  fi
  "${DIVCON_BIN}" Bace_030215.pdb --ligand Bace_030215_CAT_4p.mol2 -h amberff14sb \
    --mtdock Bace_030215_CAT_4p-dock.sdf --mtscore ensemble --np 4 -v 2 -O
}

tutorial_4() {
  section "Fourth Tutorial: 2Step MTScoreE with External Docker (no-MTCS)"
  safe_cd_root
  clean_make_cd "2Step"
  fetch "http://downloads.quantumbioinc.com/media/tutorials/MT/Bace_030215_CAT_4p.tar.gz"
  tar xvf Bace_030215_CAT_4p.tar.gz >/dev/null
  cp Bace_030215_CAT_4p/Bace_030215.pdb .
  cp Bace_030215_CAT_4p/Bace_030215_CAT_4p.mol2 .
  if need_moebatch "Fourth Tutorial docking step"; then
    "${MOEBATCH_BIN}" -licwait -run "${QBHOME}/svl/run/qbDockPair.svl" \
      -rec Bace_030215.pdb -lig Bace_030215_CAT_4p.mol2 -delwat -protonate
  fi
  "${DIVCON_BIN}" pro_Bace_030215_CAT_4p_predock.pdb \
    --ligand lig_Bace_030215_CAT_4p_predock.mol2 -h amberff14sb \
    --mtdock CAT_4p-dock.sdf --mtscore ensemble --np 4 -v 2
}

tutorial_5() {
  section "Fifth Tutorial: MTScoreE with induced fit docking"
  safe_cd_root
  clean_make_cd "InducedFit"
  fetch "http://downloads.quantumbioinc.com/media/tutorials/MT/Bace_030215_CAT_4p.tar.gz"
  tar xvf Bace_030215_CAT_4p.tar.gz >/dev/null
  cp Bace_030215_CAT_4p/Bace_030215.pdb .
  cp Bace_030215_CAT_4p/Bace_030215_CAT_4p.mol2 .
  if need_moebatch "Fifth Tutorial induced-fit docking"; then
    "${MOEBATCH_BIN}" -licwait -run "${QBHOME}/svl/run/qbDockPair.svl" \
      -rec Bace_030215.pdb -lig Bace_030215_CAT_4p.mol2 -inducedfit -delwat -protonate
  fi
  "${DIVCON_BIN}" pro_Bace_030215_CAT_4p_predock.pdb \
    --ligand lig_Bace_030215_CAT_4p_predock.mol2 -h amberff14sb \
    --mtdock *-dock*.pdb --mtscore ensemble --np 4 -v 2
}

tutorial_6() {
  section "Sixth Tutorial: Manipulating cutoffs to maximize predictions"
  safe_cd_root
  clean_make_cd "Cutoffs"
  fetch "http://downloads.quantumbioinc.com/media/tutorials/MT/Bace_030215_CAT_4p.tar.gz"
  tar xvf Bace_030215_CAT_4p.tar.gz >/dev/null
  cp Bace_030215_CAT_4p/Bace_030215.pdb .
  cp Bace_030215_CAT_4p/Bace_030215_CAT_4p.mol2 .
  if need_moebatch "Sixth Tutorial docking"; then
    "${MOEBATCH_BIN}" -licwait -run "${QBHOME}/svl/run/qbDockPair.svl" \
      -rec Bace_030215.pdb -lig Bace_030215_CAT_4p.mol2 -delwat -protonate
  fi
  "${DIVCON_BIN}" pro_Bace_030215_CAT_4p_predock.pdb \
    --ligand lig_Bace_030215_CAT_4p_predock.mol2 -h amberff14sb \
    --mtdock Bace_030215_CAT_4p-dock.sdf --mtscore ensemble --nb-cutoff 8.0 --np 4 -v 2
}

tutorial_7() {
  section "Seventh Tutorial: MTScoreE external third party docked poses (no Opt)"
  safe_cd_root
  local tutorFolder="externalSDF_noOpt"
  clean_make_cd "${tutorFolder}"
  fetch "http://downloads.quantumbioinc.com/media/tutorials/MT/pro_4wiv_ligand_predock.pdb"
  fetch "http://downloads.quantumbioinc.com/media/tutorials/MT/lig_4wiv_ligand_predock.mol2"
  fetch "http://downloads.quantumbioinc.com/media/tutorials/MT/4wiv_ligand_dock.sdf"
  "${DIVCON_BIN}" pro_4wiv_ligand_predock.pdb \
    --ligand lig_4wiv_ligand_predock.mol2 \
    --mtdock 4wiv_ligand_dock.sdf opt off --mtscore -h garf -O --np 2 -v 2
}

tutorial_8() {
  section "Eighth Tutorial: MTScoreE external docked poses with Rigid Pose Opt"
  safe_cd_root
  local tutorFolder="externalSDF_RigidOpt"
  clean_make_cd "${tutorFolder}"
  fetch "http://downloads.quantumbioinc.com/media/tutorials/MT/pro_4wiv_ligand_predock.pdb"
  fetch "http://downloads.quantumbioinc.com/media/tutorials/MT/lig_4wiv_ligand_predock.mol2"
  fetch "http://downloads.quantumbioinc.com/media/tutorials/MT/4wiv_ligand_dock.sdf"
  "${DIVCON_BIN}" pro_4wiv_ligand_predock.pdb \
    --ligand lig_4wiv_ligand_predock.mol2 \
    --mtdock 4wiv_ligand_dock.sdf opt rigid 100 0.01 --mtscore -h garf -O --np 2 -v 2
}

tutorial_9() {
  section "Ninth Tutorial: MTScoreE external docked poses with Torsion Pose Opt"
  safe_cd_root
  local tutorFolder="externalSDF_TorsionOpt"
  clean_make_cd "${tutorFolder}"
  fetch "http://downloads.quantumbioinc.com/media/tutorials/MT/pro_4wiv_ligand_predock.pdb"
  fetch "http://downloads.quantumbioinc.com/media/tutorials/MT/lig_4wiv_ligand_predock.mol2"
  fetch "http://downloads.quantumbioinc.com/media/tutorials/MT/4wiv_ligand_dock.sdf"
  "${DIVCON_BIN}" pro_4wiv_ligand_predock.pdb \
    --ligand lig_4wiv_ligand_predock.mol2 \
    --mtdock 4wiv_ligand_dock.sdf opt torsion 100 0.01 --mtscore -h garf -O --np 2 -v 2
}

tutorial_10() {
  section "Tenth Tutorial: MTScoreE external docked poses with All Atom Pose Opt"
  safe_cd_root
  local tutorFolder="externalSDF_AllOpt"
  clean_make_cd "${tutorFolder}"
  fetch "http://downloads.quantumbioinc.com/media/tutorials/MT/pro_4wiv_ligand_predock.pdb"
  fetch "http://downloads.quantumbioinc.com/media/tutorials/MT/lig_4wiv_ligand_predock.mol2"
  fetch "http://downloads.quantumbioinc.com/media/tutorials/MT/4wiv_ligand_dock.sdf"
  "${DIVCON_BIN}" pro_4wiv_ligand_predock.pdb \
    --ligand lig_4wiv_ligand_predock.mol2 \
    --mtdock 4wiv_ligand_dock.sdf opt all 100 0.01 --mtscore -h garf -O --np 2 -v 2
}

#############################################
# Menu + Dispatch
#############################################
print_menu() {
  cat <<'EOF'
Interactive Tutorial Menu (only shown with -i):
  1  First A   : MTConfSearch (SMILES)
  2  First B   : MTConfSearch (MOL2)
  3  First C   : MTConfSearch Energy Cutoff
  4  First D   : MTConfSearch Energy Cutoff With Opt
  5  Second    : MTScore (Endstate)
  6  Third     : 3Step MTCS+MTScoreE (Ensemble) External Docker
  7  Fourth    : 2Step MTScoreE External Docker (no-MTCS)
  8  Fifth     : MTScoreE induced fit docking
  9  Sixth     : Manipulating cutoffs
 10  Seventh   : External SDF no Opt
 11  Eighth    : External SDF Rigid Pose Opt
 12  Ninth     : External SDF Torsion Pose Opt
 13  Tenth     : External SDF All Atom Pose Opt
  0 / A        : Run All
  Q            : Quit
EOF
}

run_all() {
  tutorial_1a
  tutorial_1b
  tutorial_1c
  tutorial_1d
  tutorial_2
  tutorial_3
  tutorial_4
  tutorial_5
  tutorial_6
  tutorial_7
  tutorial_8
  tutorial_9
  tutorial_10
}

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Default behavior (no options): run ALL tutorials non-interactively.

Options:
  -i        Interactive menu mode (lets you pick specific tutorials)
  -l        List tutorials (no execution)
  -h        Show this help

Examples:
  $(basename "$0")        # Run all tutorials
  $(basename "$0") -i     # Show interactive menu
  $(basename "$0") -l     # List tutorials then exit
EOF
}

dispatch_number() {
  case "$1" in
    1) tutorial_1a ;;
    2) tutorial_1b ;;
    3) tutorial_1c ;;
    4) tutorial_1d ;;
    5) tutorial_2 ;;
    6) tutorial_3 ;;
    7) tutorial_4 ;;
    8) tutorial_5 ;;
    9) tutorial_6 ;;
    10) tutorial_7 ;;
    11) tutorial_8 ;;
    12) tutorial_9 ;;
    13) tutorial_10 ;;
    0|A|a) run_all ;;
    Q|q) log "User requested quit."; exit 0 ;;
    *) echo "WARNING: Unknown selection: $1" ;;
  esac
}

#############################################
# Main
#############################################
main() {
  log "BEGIN Tutorial Test using ${DIVCON_BIN}"
  log "MOEBATCH available: ${HAVE_MOEBATCH}"

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