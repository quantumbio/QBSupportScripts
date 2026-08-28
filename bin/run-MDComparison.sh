#!/bin/bash

set -euo pipefail

# =============================================================================
# DivCon/OpenMM C++ vs Python validation
# =============================================================================

usage() {
    echo "Usage: $0 PDBID" >&2
    echo "       $0 PDBID_LIST [percent]" >&2
}

if (( $# < 1 || $# > 2 )); then
    usage
    exit 2
fi

INPUT=$1
PCT=${2:-100}

: "${DIVCON:?DIVCON must be set before running this script}"
: "${QBSUPPORT:?QBSUPPORT must be set before running this script}"

if ! command -v python >/dev/null 2>&1; then
    echo "ERROR: python is not available in PATH." >&2
    exit 1
fi

if ! python -c 'import openmm' >/dev/null 2>&1; then
    echo "ERROR: OpenMM is not available in the current Python environment." >&2
    echo "       Python: $(command -v python)" >&2
    exit 1
fi

PYTHON_MD=${QBSUPPORT}/python/runMD-fromDivCon.py
ANALYSIS_MD=${QBSUPPORT}/python/analysisMD.py

PYTHON=python

# These values are supplied identically to both MD implementations.
#
# EQUILIBRATION_STEPS is used for BOTH NVT and NPT equilibration.
# PRODUCTION_STEPS is the production MD length.
#
# Use defaults only when the variables have not already been supplied
# by the calling environment.
: "${MINIMIZATION_STEPS:=5000}"
: "${EQUILIBRATION_STEPS:=500}"
: "${PRODUCTION_STEPS:=2000}"
: "${REPORT_INTERVAL:=500}"

export MINIMIZATION_STEPS
export EQUILIBRATION_STEPS
export PRODUCTION_STEPS
export REPORT_INTERVAL

# =============================================================================
# Input mode
# =============================================================================

# If the first argument names a regular file, treat it as a PDBID list.  The
# optional percentage selects a random subset, rounded up as in run-prepareALL.sh.
# Each selected PDBID is then processed by a separate invocation of this script,
# preserving the existing single-PDB workflow unchanged.
if [[ -f "${INPUT}" && -z "${MD_COMPARISON_SINGLE:-}" ]]; then
    [[ ${PCT} =~ ^[0-9]+$ ]] || {
        echo "ERROR: percent must be an integer from 0 through 100." >&2
        exit 2
    }
    (( PCT >= 0 && PCT <= 100 )) || {
        echo "ERROR: percent must be from 0 through 100." >&2
        exit 2
    }

    TOTAL=$(awk 'NF && $0 !~ /^[[:space:]]*#/ { n++ } END { print n + 0 }' "${INPUT}")
    TAKE=$(( (TOTAL * PCT + 99) / 100 ))

    # Match run-prepareALL.sh by honoring PBS_NUM_PPN when supplied.  The
    # MD-specific override is useful when each comparison itself uses many CPU
    # threads and a lower outer concurrency is desired.
    N=${MD_PARALLEL_JOBS:-${PBS_NUM_PPN:-4}}
    [[ ${N} =~ ^[1-9][0-9]*$ ]] || {
        echo "ERROR: MD_PARALLEL_JOBS/PBS_NUM_PPN must be a positive integer." >&2
        exit 2
    }

    SCRIPT_PATH="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)/$(basename -- "${BASH_SOURCE[0]}")"

    echo "======================================================================"
    echo "DivCon/OpenMM validation list mode"
    echo "Input list          : ${INPUT}"
    echo "Structures available: ${TOTAL}"
    echo "Percent requested   : ${PCT}%"
    echo "Structures selected : ${TAKE}"
    echo "Parallel jobs       : ${N}"
    echo "======================================================================"

    if (( TAKE == 0 )); then
        exit 0
    fi

    export MD_COMPARISON_SINGLE=1

    if command -v shuf >/dev/null 2>&1; then
        awk 'NF && $0 !~ /^[[:space:]]*#/ { print }' "${INPUT}" \
            | shuf -n "${TAKE}" \
            | xargs -P "${N}" -n1 bash "${SCRIPT_PATH}"
    else
        # Portable fallback when GNU shuf is unavailable.
        awk 'BEGIN { srand() }
             NF && $0 !~ /^[[:space:]]*#/ { printf "%.12f\t%s\n", rand(), $0 }' "${INPUT}" \
            | sort -k1,1 \
            | cut -f2- \
            | sed -n "1,${TAKE}p" \
            | xargs -P "${N}" -n1 bash "${SCRIPT_PATH}"
    fi

    exit $?
fi

if (( $# == 2 )); then
    echo "ERROR: an optional percentage is only valid when the first argument is a PDBID list file." >&2
    usage
    exit 2
fi

PDBID=${INPUT}

echo "======================================================================"
echo "DivCon/OpenMM validation: ${PDBID}"
echo "Started: $(date)"
echo
echo "Simulation protocol:"
echo "  NVT equilibration : ${EQUILIBRATION_STEPS} steps"
echo "  NPT equilibration : ${EQUILIBRATION_STEPS} steps"
echo "  Production        : ${PRODUCTION_STEPS} steps"
echo "======================================================================"

# -----------------------------------------------------------------------------
# DivCon environment
# -----------------------------------------------------------------------------

# Do not source qbenv.sh here.  Keep the validation environment independent
# of any library/PATH changes made by the DivCon environment setup.
export DIVCON

# Development switches still required by the current MD implementation.
export FORCE_INPUT_EXE=1
export MT_ALLOW_MTDOCK=1
export DEV_SHOW_PERCEPTION=1
export QB_SKIP_END_GAP_PREPARE=1
export DEV_MD_MINIMIZE=1
export DEV_MIN_STEPS=5000

# -----------------------------------------------------------------------------
# Primary working directory
# -----------------------------------------------------------------------------

mkdir -p "${PDBID}"
cd "${PDBID}"

ROOT_WORKDIR=$(pwd)

mkdir -p cpp python

# Remove root-level files generated by previous runs.
rm -f \
    *.h5 *.h5.* \
    *.pdb *.cif *.mol2 *.mtz *.map *.json \
    *.qb.log *.config \
    *.xml \
    output.dcd output.prmtop output.inpcrd \
    prepare.screenout analysis.screenout

# Remove previous analysis products for this system.
rm -f \
    "${PDBID}"_*.pdf \
    "${PDBID}"_*.csv \
    "${PDBID}"_summary.json

# Clear the two isolated MD working directories.
rm -rf cpp/*
rm -rf python/*

# =============================================================================
# 1. Prepare common input system in ROOT directory
# =============================================================================

echo
echo "======================================================================"
echo "Preparing common waterboxed system"
echo "Working directory: ${ROOT_WORKDIR}"
echo "Screen output    : ${ROOT_WORKDIR}/prepare.screenout"
echo "======================================================================"

"${DIVCON}/bin/qmechanic" "${PDBID}" \
    --prepare sidechains \
    -h amberff14sb \
    --waterbox \
    --np 16 \
    -O \
    --symmetry off \
    --nb-cutoff 10 \
    -p \
        "${PDBID}+H+wb.pdb" \
        "${PDBID}+wb.parm7" \
        "${PDBID}+wb.inpcrd" \
    > prepare.screenout 2>&1
    
# -----------------------------------------------------------------------------
# Give BOTH implementations private copies of exactly the same prepared system.
# -----------------------------------------------------------------------------

for RUN_DIR in cpp python; do
    cp "${PDBID}+H+wb.pdb" "${RUN_DIR}/"
    cp "${PDBID}+wb.parm7" "${RUN_DIR}/"
    cp "${PDBID}+wb.inpcrd" "${RUN_DIR}/"
done

# =============================================================================
# 2. Launch DivCon/OpenMM C++ MD
# =============================================================================

echo
echo "======================================================================"
echo "Launching DivCon/OpenMM C++ MD"
echo "Working directory: ${ROOT_WORKDIR}/cpp"
echo "Screen output    : ${ROOT_WORKDIR}/cpp/md.screenout"
echo "======================================================================"

(
    cd "${ROOT_WORKDIR}/cpp"

    export DEV_MIN_STEPS="${MINIMIZATION_STEPS}"
    export DEV_MD_REPORT_INTERVAL="${REPORT_INTERVAL}"
    
    "${DIVCON}/bin/qmechanic" "${PDBID}+H+wb.pdb" \
        -h amberff14sb \
        --np 20 -v 3 \
        -O \
        --symmetry off \
        --nb-cutoff 10 \
        --opt all 1 \
        --moleculardynamics \
            "${EQUILIBRATION_STEPS},${PRODUCTION_STEPS}" \
        -p "${PDBID}_md.pdb"
    
) > "${ROOT_WORKDIR}/cpp/md.screenout" 2>&1 &

CPP_PID=$!

# =============================================================================
# 3. Launch Python/OpenMM MD
# =============================================================================

echo
echo "======================================================================"
echo "Launching Python/OpenMM MD"
echo "Working directory: ${ROOT_WORKDIR}/python"
echo "Screen output    : ${ROOT_WORKDIR}/python/md.screenout"
echo "======================================================================"

(
    cd "${ROOT_WORKDIR}/python"

    # Keep the Conda/OpenMM environment changes local to the Python process.

    "${PYTHON}" "${PYTHON_MD}" \
        --skip-waterbox \
        --minimization-steps "${MINIMIZATION_STEPS}" \
        --equilibration-steps "${EQUILIBRATION_STEPS}" \
        --production-steps "${PRODUCTION_STEPS}" \
        --report-interval "${REPORT_INTERVAL}" \
        "${PDBID}+wb.parm7" \
        "${PDBID}+wb.inpcrd"
    
) > "${ROOT_WORKDIR}/python/md.screenout" 2>&1 &

PYTHON_PID=$!

echo
echo "Both MD calculations are running in parallel."
echo "  C++ PID    : ${CPP_PID}"
echo "  Python PID : ${PYTHON_PID}"

# =============================================================================
# 4. Wait for BOTH calculations
# =============================================================================

CPP_STATUS=0
PYTHON_STATUS=0

wait "${CPP_PID}" || CPP_STATUS=$?
wait "${PYTHON_PID}" || PYTHON_STATUS=$?

echo
echo "======================================================================"
echo "MD calculations completed"
echo "  C++ status    : ${CPP_STATUS}"
echo "  Python status : ${PYTHON_STATUS}"
echo "======================================================================"

if (( CPP_STATUS != 0 )); then
    echo "ERROR: DivCon/OpenMM C++ MD failed."
    echo "See ${ROOT_WORKDIR}/cpp/md.screenout"
    exit "${CPP_STATUS}"
fi

if (( PYTHON_STATUS != 0 )); then
    echo "ERROR: Python/OpenMM MD failed."
    echo "See ${ROOT_WORKDIR}/python/md.screenout"
    exit "${PYTHON_STATUS}"
fi

# =============================================================================
# 5. Analyze trajectories in ROOT working directory
# =============================================================================

cd "${ROOT_WORKDIR}"

echo
echo "======================================================================"
echo "Analyzing C++ vs Python trajectories"
echo "Working directory: ${ROOT_WORKDIR}"
echo "Screen output    : ${ROOT_WORKDIR}/analysis.screenout"
echo "======================================================================"

"${PYTHON}" "${ANALYSIS_MD}" \
    cpp/output.dcd \
    "cpp/${PDBID}+wb.parm7" \
    python/output.dcd \
    python/output.prmtop \
    --label1 "C++" \
    --label2 "Python" \
    --out1 cpp/md.screenout \
    --out2 python/md.screenout \
    -o "${PDBID}" \
    > analysis.screenout 2>&1

echo
echo "======================================================================"
echo "Completed: $(date)"
echo
echo "Screen outputs:"
echo "  Preparation : ${ROOT_WORKDIR}/prepare.screenout"
echo "  C++ MD      : ${ROOT_WORKDIR}/cpp/md.screenout"
echo "  Python MD   : ${ROOT_WORKDIR}/python/md.screenout"
echo "  Analysis    : ${ROOT_WORKDIR}/analysis.screenout"
echo
echo "Trajectory outputs:"
echo "  C++         : ${ROOT_WORKDIR}/cpp/output.dcd"
echo "  Python      : ${ROOT_WORKDIR}/python/output.dcd"
echo "======================================================================"
