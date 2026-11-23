#!/bin/bash

if [[ -z "${SPAECIES_ROOT:-}" ]]; then
    SPAECIES_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi

: "${SPAECIES_BUILD_DIR:=${SPAECIES_ROOT}/build}"
: "${RAINSHAFT_DIR:=${SPAECIES_BUILD_DIR}/rainshaft}"
: "${RAINSHAFT_EXE:=${RAINSHAFT_DIR}/${RAINSHAFT_EXE_NAME:-rainshaft}}"
: "${SPAECIES_INPUT_FILE:=${SPAECIES_ROOT}/input/random_rainshaft_samples_12mo.nc}"
: "${IC_FILE:=${SPAECIES_INPUT_FILE}}"
: "${SPAECIES_RESULTS_DIR:=${SPAECIES_ROOT}/results}"

if [[ -z "${SAVE_DIR:-}" ]]; then
    if [[ -n "${SAVE_SUBDIR:-}" ]]; then
        SAVE_DIR="${SPAECIES_RESULTS_DIR}/${SAVE_SUBDIR}"
    else
        SAVE_DIR="${SPAECIES_RESULTS_DIR}"
    fi
fi

if [[ ! -d "${RAINSHAFT_DIR}" ]]; then
    echo "RAINSHAFT_DIR does not exist: ${RAINSHAFT_DIR}" >&2
    echo "Set SPAECIES_BUILD_DIR or RAINSHAFT_DIR to the rainshaft build directory." >&2
    exit 1
fi

if [[ ! -x "${RAINSHAFT_EXE}" ]]; then
    echo "RAINSHAFT_EXE is not executable: ${RAINSHAFT_EXE}" >&2
    echo "Set RAINSHAFT_EXE or RAINSHAFT_EXE_NAME to the rainshaft executable." >&2
    exit 1
fi

if [[ ! -f "${IC_FILE}" ]]; then
    echo "IC_FILE does not exist: ${IC_FILE}" >&2
    echo "Set SPAECIES_INPUT_FILE or IC_FILE to the input NetCDF file." >&2
    exit 1
fi

mkdir -p "${SAVE_DIR}"
