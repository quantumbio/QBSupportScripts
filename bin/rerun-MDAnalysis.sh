#!/bin/bash

set -euo pipefail

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

: "${QBSUPPORT:?QBSUPPORT must be set before running this script}"

PYTHON_BIN=${PYTHON_BIN:-python}
ANALYSIS_MD=${ANALYSIS_MD:-${QBSUPPORT}/python/analysisMD.py}

if ! command -v "${PYTHON_BIN}" >/dev/null 2>&1; then
    echo "ERROR: ${PYTHON_BIN} is not available in PATH." >&2
    exit 1
fi

if ! "${PYTHON_BIN}" -c 'import openmm, mdtraj' >/dev/null 2>&1; then
    echo "ERROR: OpenMM/MDTraj are not available in the selected Python environment." >&2
    echo "       Python: $(command -v "${PYTHON_BIN}")" >&2
    exit 1
fi

if [[ ! -f "${ANALYSIS_MD}" ]]; then
    echo "ERROR: analysisMD.py not found: ${ANALYSIS_MD}" >&2
    exit 1
fi

# Match run-MDComparison.sh list mode, but run only the existing trajectory/XML analysis.
if [[ -f "${INPUT}" && -z "${MD_ANALYSIS_SINGLE:-}" ]]; then
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
    N=${MD_ANALYSIS_JOBS:-${MD_PARALLEL_JOBS:-${PBS_NUM_PPN:-4}}}
    [[ ${N} =~ ^[1-9][0-9]*$ ]] || {
        echo "ERROR: MD_ANALYSIS_JOBS/MD_PARALLEL_JOBS/PBS_NUM_PPN must be a positive integer." >&2
        exit 2
    }

    SCRIPT_PATH="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)/$(basename -- "${BASH_SOURCE[0]}")"

    echo "======================================================================"
    echo "OpenMM validation analysis-only rerun"
    echo "Input list          : ${INPUT}"
    echo "Structures available: ${TOTAL}"
    echo "Percent requested   : ${PCT}%"
    echo "Structures selected : ${TAKE}"
    echo "Parallel jobs       : ${N}"
    echo "======================================================================"

    if (( TAKE == 0 )); then
        exit 0
    fi

    export MD_ANALYSIS_SINGLE=1

    if command -v shuf >/dev/null 2>&1; then
        awk 'NF && $0 !~ /^[[:space:]]*#/ { print }' "${INPUT}" \
            | shuf -n "${TAKE}" \
            | xargs -P "${N}" -n1 bash "${SCRIPT_PATH}"
    else
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

if [[ ! -d "${PDBID}" ]]; then
    echo "ERROR: comparison directory not found: ${PDBID}" >&2
    exit 1
fi

cd "${PDBID}"
ROOT_WORKDIR=$(pwd)

required_files=(
    "cpp/output.dcd"
    "cpp/${PDBID}+wb.parm7"
    "cpp/md.screenout"
    "python/output.dcd"
    "python/output.prmtop"
    "python/md.screenout"
)

missing=0
for path in "${required_files[@]}"; do
    if [[ ! -f "${path}" ]]; then
        echo "ERROR: missing required analysis input: ${ROOT_WORKDIR}/${path}" >&2
        missing=1
    fi
done
if (( missing != 0 )); then
    exit 1
fi

# Remove only products created by analysisMD.py.  Preserve all MD/preparation inputs,
# trajectories, screenouts, and serialized OpenMM XML files.
rm -f \
    analysis.screenout \
    "${PDBID}"_*.pdf \
    "${PDBID}"_*.csv \
    "${PDBID}"_summary.json

echo "======================================================================"
echo "Reanalyzing existing C++ vs Python MD outputs: ${PDBID}"
echo "Working directory: ${ROOT_WORKDIR}"
echo "Analysis script  : ${ANALYSIS_MD}"
echo "======================================================================"

if "${PYTHON_BIN}" "${ANALYSIS_MD}" \
    cpp/output.dcd \
    "cpp/${PDBID}+wb.parm7" \
    python/output.dcd \
    python/output.prmtop \
    --label1 "C++" \
    --label2 "Python" \
    --out1 cpp/md.screenout \
    --out2 python/md.screenout \
    -o "${PDBID}" \
    > analysis.screenout 2>&1; then
    :
else
    status=$?
    echo "ERROR: analysis failed for ${PDBID}; see ${ROOT_WORKDIR}/analysis.screenout" >&2
    exit "${status}"
fi

echo "Completed analysis: ${PDBID}"
echo "  ${ROOT_WORKDIR}/analysis.screenout"
echo "  ${ROOT_WORKDIR}/${PDBID}_summary.json"
