#!/usr/bin/env bash
set -u -o pipefail

# MT-threaded spot test script (menu-driven refactor; corrected to remove erroneous backticks)

: "${QBHOME:?Set QBHOME or source /path/to/DivConSuite/etc/qbenv.sh before running this script.}"

DIVCON="${QBHOME}/bin/qmechanic"
if [[ ! -x "${DIVCON}" ]]; then
  echo "ERROR: ${DIVCON} not executable"
  exit 1
fi

WORKDIR="${PWD}"
DATE_FMT="+%Y-%m-%d %H:%M:%S"

# For these tests, MTDock is required.
export MT_ALLOW_MTDOCK=1

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
  local file
  file="$(basename "${url}")"
  if [[ -f "${file}" ]]; then
    log "Already present: ${file}"
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

###############################################################################
# Tutorial Functions
# NOTE: Replace any [...]-tagged comment placeholders with the complete,
# original arguments from your canonical script.
###############################################################################

tutorial_1a_pdb() {
  section "Tutorial #1a-pdb: Automated DivCon-based Docking – qmechanic + XModeScore (ONIOM)"
  safe_cd_root
  local dir="XModeScore-dock_ONIOM"
  clean_make_cd "${dir}"
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/5C3K-H_refine_001.pdb
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/5C3K-H_refine_001.mtz
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/4XF_A_402.pdb

  local args=(
    "5C3K-H_refine_001.pdb"
    "5C3K-H_refine_001.mtz"
    -h pm6 amberff14sb
    --ligand 4XF_A_402.pdb
    --qm-region LIGAND+/A/ZN/401// 3.0 0
    --confSearch LIGAND opt off
    --flip all
    --chirality all
    --dock 5,1SD opt torsion pocket 100 0.01
    --protomers all [-1..1]
    --xmodescore opt all
    --np 8 -v 2
    -p 5C3K-out.pdb 5C3K-out.mtz
  )
  # (Above: fill in any missing flags that were truncated.)
  "${DIVCON}" "${args[@]}" >& OUT.SCREEN
}

tutorial_1a_mol2() {
  section "Tutorial #1a-mol2: Automated DivCon-based Docking – qmechanic + XModeScore (MM)"
  safe_cd_root
  local dir="XModeScore-dock_MM"
  clean_make_cd "${dir}"
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/5C3K-H_refine_001.pdb
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/5C3K-H_refine_001.mtz
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/4XF_A_402.pdb

  local args=(
    "5C3K-H_refine_001.pdb"
    "5C3K-H_refine_001.mtz"
    -h amberff14sb
    --ligand 4XF_A_402.pdb
    --confSearch LIGAND opt off
    --flip off
    --chirality off
    --dock 5,1SD opt torsion LIGAND 100 0.01
    --protomers all [-1..1]
    --xmodescore opt torsion pocket
    --np 8 -v 2
    -p 5C3K-out.mol2 5C3K-out.mtz
  )
  "${DIVCON}" "${args[@]}" >& OUT.SCREEN
}

tutorial_1b() {
  section "Tutorial #1b: Automated Docking – qmechanic + EndState MTScore"
  safe_cd_root
  local dir="MTDock-MTScore_garf"
  clean_make_cd "${dir}"
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/5C3K-H_refine_001.pdb
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/5C3K-H_refine_001.mtz
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/4XF_A_402.pdb

  local args=(
    "5C3K-H_refine_001.pdb"
    "5C3K-H_refine_001.mtz"
    -h garf
    --ligand 4XF_A_402.pdb
    --mtcs LIGAND opt off
    --flip all
    --chirality all
    --mtdock 5,1SD opt torsion pocket 100 0.01
    --protomers off
    --mtscore endstate
    --np 8 -v 2
    -p 5C3K-out.mol2 5C3K-out.mtz
  )
  "${DIVCON}" "${args[@]}" >& OUT.SCREEN
}

tutorial_1c() {
  section "Tutorial #1c: Automated Docking – qmechanic + Ensemble MTScore"
  safe_cd_root
  local dir="MTDock-MTScore_amberff14sb"
  clean_make_cd "${dir}"
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/5C3K-H_refine_001.pdb
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/5C3K-H_refine_001.mtz
  fetch https://raw.githubusercontent.com/quantumbio/QBSupportScripts/master/Tutorials/data/XModeScore/4XF_A_402.pdb

  local args=(
    "5C3K-H_refine_001.pdb"
    "5C3K-H_refine_001.mtz"
    -h amberff14sb
    --ligand 4XF_A_402.pdb
    --confSearch LIGAND opt off
    --flip all
    --chirality all
    --mtdock 5,1SD opt torsion pocket 100 0.01
    --protomers all [-1..1]
    --mtscore ensemble
    --np 8 -v 2
    -p 5C3K-out.mol2 5C3K-out.mtz
  )
  "${DIVCON}" "${args[@]}" >& OUT.SCREEN
}

tutorial_1d() {
  section "Tutorial #1d: Prepare then Docking – qmechanic + XModeScore (ONIOM)"
  safe_cd_root
  local dir="prepare_XModeScore-dock_ONIOM"
  clean_make_cd "${dir}"

  # Preparation step
  "${DIVCON}" 5C3K \
    -h amberff14sb \
    --prepare all \
    -p 5C3K-refined.pdb 5C3K-refined.mtz \
    -v 2 --np 8 >& OUT.PREPARE

  # Docking / XModeScore step
  local args=(
    "5C3K-refined.pdb"
    "5C3K-refined.mtz"
    -O
    -h pm6 amberff14sb
    --ligand /A/4XF/402//
    --qm-region LIGAND+/A/ZN/401// 3.0 0
    --confSearch LIGAND opt off
    --flip all
    --chirality all
    --dock 5,1SD opt torsion pocket 100 0.01
    --protomers all [-1..1]
    --xmodescore opt all
    --np 8 -v 2
    -p 5C3K-out.pdb 5C3K-out.mtz
  )
  "${DIVCON}" "${args[@]}" >& OUT.SCREEN
}

###############################################################################
# Menu / Dispatch
###############################################################################

print_menu() {
  cat <<'EOF'
Interactive Tutorial Menu (-i to enable):
  1  Tutorial #1a-pdb   : XModeScore Dock (ONIOM)
  2  Tutorial #1a-mol2  : XModeScore Dock (MM)
  3  Tutorial #1b       : EndState MTScore
  4  Tutorial #1c       : Ensemble MTScore
  5  Tutorial #1d       : Prepare + XModeScore Dock (ONIOM)
  0 / A                 : Run All
  Q                     : Quit
EOF
}

run_all() {
  tutorial_1a_pdb
  tutorial_1a_mol2
  tutorial_1b
  tutorial_1c
  tutorial_1d
}

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Default (no options): run ALL tutorial segments.

Options:
  -i    Interactive menu mode
  -l    List tutorials only
  -h    Help
EOF
}

dispatch() {
  case "$1" in
    1) tutorial_1a_pdb ;;
    2) tutorial_1a_mol2 ;;
    3) tutorial_1b ;;
    4) tutorial_1c ;;
    5) tutorial_1d ;;
    0|A|a) run_all ;;
    Q|q) log "Quit requested"; exit 0 ;;
    *) echo "WARNING: Unknown selection: $1" ;;
  esac
}

main() {
  log "BEGIN new_parallel_branch tutorial batch (DIVCON=${DIVCON})"

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

  log "END new_parallel_branch tutorial batch"
}

main "$@"