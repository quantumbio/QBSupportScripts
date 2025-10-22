#!/usr/bin/env bash
# Refactored: Menu-driven version (default: run ALL tutorials non-interactively)
# Original script: BUSTER/DivCon tutorials (qm-region setup, qbbuster XModeScore scenarios)
#
# NOTE:
#   - Original command lines (including truncated segments ending with [...]) are
#     preserved EXACTLY to avoid altering intended flags. Replace [...] with full
#     arguments from your canonical internal script if needed.
#   - Added interactive menu (-i), list (-l), help (-h).
#   - Added safety helpers (fetch, directory creation, logging).
#   - The original script printed multiple END markers mid-way; here we log a single
#     final END after all selected tutorials complete.
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

: "${QBHOME:?ERROR: QBHOME is not set! Source /path/to/DivConSuite/etc/qbenv.sh THEN run this script.}"
: "${BDG_home:?ERROR: BDG_home is not set! Source your BUSTER setup.sh THEN run this script.}"

DIVCON_BIN="${QBHOME}/bin/qmechanic"
if [[ ! -x "${DIVCON_BIN}" ]]; then
  echo "ERROR: ${DIVCON_BIN} is not an executable! Adjust QBHOME."
  exit 1
fi

# qbbuster executable (normally on PATH after BUSTER setup)
if ! command -v qbbuster >/dev/null 2>&1; then
  echo "WARNING: 'qbbuster' not found in PATH. Tutorials requiring it will likely fail."
fi

WORKDIR="${PWD}"
DATE_FMT="+%Y-%m-%d %H:%M:%S"

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
    echo "ERROR: Need wget or curl installed."
    return 1
  fi
}

###############################################################################
# Tutorial Functions (commands preserved exactly; truncated args kept)
###############################################################################

tutorial_1() {
  section "Tutorial #1 (Running BUSTER/DivCon without the qbbuster wrapper): 2BSM"
  safe_cd_root
  clean_make_cd "NakedBUSTER"

  # Create configuration files as in original script
  tee divcon.ini >/dev/null <<'EOF'
qm-region = /A/BSM/1224// 3 0
hamiltonian = pm6 amberff14sb
np = 4
EOF

  tee 2BSM.qm >/dev/null <<'EOF'
NOTE BUSTER_SET QM01 =  ALL
NOTE BUSTER_QM_CHARGE_01 0
NOTE BUSTER_QM_MULTIP_01 1
NOTE BUSTER_QM_WRESTR_01 0.0
NOTE BUSTER_QM_WRITE_PDB_FOR_HELPER ON
NOTE BUSTER_QM_MAXDISP 1.0E05
NOTE BUSTER_QM_HELPER $QBHOME/scripts/BusterQMHelper.sh
NOTE BUSTER_QM_METHOD PM6
NOTE BUSTER_QM_WEIGHT 3.0
EOF

  fetch http://downloads.quantumbioinc.com/media/tutorials/XModeScore/2BSM+H.pdb
  fetch http://downloads.quantumbioinc.com/media/tutorials/XModeScore/2BSM.mtz
  fetch http://downloads.quantumbioinc.com/media/tutorials/XModeScore/BSM.cif

  # Original refine command was commented out; preserved here:
  # refine -Gelly 2BSM.qm -p 2BSM+H.pdb -m 2BSM.mtz -l BSM.cif -nbig 2  RunGellySanityCheck=no RunGellyScreen=no -qm_weight 3.0
}

tutorial_2() {
  section "Tutorial #2 (qbbuster + Grade & DivCon protonator): 2BSM"
  safe_cd_root
  clean_make_cd "qbbuster_divcon"

  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/2bsm.pdb
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/2BSM.mtz

  # Preserved truncated command line exactly:
  $QBHOME/bin/qbbuster --pdbFile 2bsm.pdb --sfFile 2BSM.mtz --protonation divcon --makeCIF divcon --mmMethod amberff14sb --qmMethod pm6 --qmWeight 5.0 --ncycles 1  --nSmallCycles 15 --selection "resname BSM" --region-radius 3.0 --np 4 --dir qmRun || { echo "Tutorial #2 failed"; return 1; }
}

tutorial_3() {
  section "Tutorial #3 (XModeScore; structure with multiple ligand copies): 4ntk"
  safe_cd_root
  clean_make_cd "qbbuster_4ntk"

  fetch http://downloads.quantumbioinc.com/media/tutorials/XModeScore/4ntk.pdb
  fetch http://downloads.quantumbioinc.com/media/tutorials/XModeScore/4ntk.mtz

  qbbuster --pdbFile 4ntk.pdb --sfFile 4ntk.mtz --XModeScore --protomers "-1" --protonation divcon --protonateTautomers divcon --makeCIF divcon --mmMethod amberff14sb --qmMethod pm6 --qmWeight 3.0 --engine buster --ncycles 1 --nSmallCycles 20 --selection "resname ZSP and resid 202 and chain A" --np 32 --buffer-radius 0.0 --region-radius 3.0 --dir xmodeScore_results || { echo "Tutorial #3 failed"; return 1; }
}

tutorial_4() {
  section "Tutorial #4 (XModeScore; multiple chiral centers): 4gr0"
  safe_cd_root
  clean_make_cd "qbbuster_4gr0"

  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data//XModeScore/4gr0+H.pdb
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data//XModeScore/4gr0.mtz

  # C31 chiral atom noted in original comment
  qbbuster --pdbFile 4gr0+H.pdb --sfFile 4gr0.mtz --XModeScore --protomers "0" --exploreChiral "C31" --protonation skip --protonateTautomers MOE --makeCIF moe --mmMethod amberff14sb --qmMethod pm6 --qmWeight 3.0 --engine buster --ncycles 1 --nSmallCycles 20 --selection "resname R4B and resid 306 and chain A" --np 32 --region-radius 3.0 --dir xmodeScore_results || { echo "Tutorial #4 failed"; return 1; }
}

tutorial_5() {
  section "Tutorial #5 (XModeScore; multiple flips): 3r4p"
  safe_cd_root
  clean_make_cd "qbbuster_3r4p"

  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data//XModeScore/3r4p+H.pdb
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data//XModeScore/3r4p.mtz

  # N6 flip exploration noted in original comment
  qbbuster --pdbFile 3r4p+H.pdb --sfFile 3r4p.mtz --XModeScore --protomers "0" --exploreFLip "N6" --protonation skip --protonateTautomers MOE --makeCIF moe --mmMethod amberff14sb --qmMethod pm6 --qmWeight 3.0 --engine buster --ncycles 1 --nSmallCycles 20 --selection "resname FU7 and resid 901 and chain A" --np 32 --region-radius 3.0 --dir xmodeScore_results || { echo "Tutorial #5 failed"; return 1; }
}

tutorial_6() {
  section "Tutorial #6 (XModeScore; PBS cluster submission example): 5kcv"
  safe_cd_root
  clean_make_cd "qbbuster_5kcv"

  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data//XModeScore/5kcv+H.pdb
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data//XModeScore/5kcv.mtz
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data//XModeScore/phenix_pbs.tmpl

  # Original notes about cluster template and submission command retained:
  # --cluster phenix_pbs.tmpl  (Set PBS template)
  # --clusterSub qsub          (Set PBS submission command)

  qbbuster --pdbFile 5kcv+H.pdb --sfFile 5kcv.mtz --XModeScore --protomers "1" --exploreChiral --exploreFlip --protonation skip --makeCIF divcon --mmMethod amberff14sb --qmMethod pm6 --selection "resname 6S1 and resid 501 and chain A" --ncycles 1 --nSmallCycles 10 --qmWeight 3.0 --engine divcon --np 4 --region-radius 3.0 --dir xmodeScore_results --cluster phenix_pbs.tmpl --clusterSub qsub || { echo "Tutorial #6 failed"; return 1; }
}

###############################################################################
# Menu / Dispatch
###############################################################################

print_menu() {
  cat <<'EOF'
BUSTER/DivCon Tutorial Menu (-i for interactive):
  1  Tutorial #1 : Naked BUSTER / DivCon (2BSM)
  2  Tutorial #2 : qbbuster Grade & DivCon protonation (2BSM)
  3  Tutorial #3 : XModeScore multiple ligand copies (4ntk)
  4  Tutorial #4 : XModeScore chiral centers (4gr0)
  5  Tutorial #5 : XModeScore flips (3r4p)
  6  Tutorial #6 : XModeScore PBS cluster example (5kcv)
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

Default (no options): Run ALL tutorials.

Options:
  -i    Interactive menu mode
  -l    List tutorials (no execution)
  -h    Help

Examples:
  $(basename "$0")
  $(basename "$0") -i
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
  log "BEGIN BUSTER/DivCon Tutorial Batch (DIVCON_BIN=${DIVCON_BIN})"

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

  log "END BUSTER/DivCon Tutorial Batch"
}

main "$@"