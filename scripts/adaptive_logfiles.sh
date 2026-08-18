#!/bin/bash

RAINSHAFT_EXE_NAME="rainshaft_sedonly_logging"
SAVE_SUBDIR="logfiles/adaptive_sedonly"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/common_paths.sh"

# simulation length
FINAL_TIME=1000

# tolerance for adaptive stepper/nls
RELTOL=(1.e-1 1.e-3 1.e-5 1.e-7)

# order of method
ORDERS=2

# process set to run
PROCESSES="all"

# number of runs to do for averaging. processing more than 1 run is currently not supported.
NUMRUNS=1

# how many time steps to save data for
STEPS_PER_OUTPUT=-1

# number of cases from IC_FILE to run (-1 for all cases)
NUM_CASES=1

# E3SM case number to run
CASE_IDX=

# toggle limiting
POSTPROCESS="false"

# toggle lookup tables
LOOKUP_FLAG="false"

# toggle legacy rain shape parameter
USE_ZERO_MUR="false"

# toggle q_sat_dry regularization
REGULARIZE_QSAT="true"
REGULARIZE_LAMBDAR="true"
EPSILON_QSAT_FAC=1e-7
QSMALL=1e-14

# type of integration (options: explicit, imex, mri)
INTEGRATION_TYPE="explicit"

SETTINGS_NAME="explicit_adaptive_withlog"

printf -- "---------------------------------- RUN %u ----------------------------------\n" $k
# loop over columns
for k in {0..999}
do
    for i in $(seq 0 $((${#RELTOL[@]} - 1)))
    do
        OUTPUT_FILE="${RAINSHAFT_DIR}/rainshaft.nc"

        LOG_FILE="${SAVE_DIR}/icase${k}_reltol${RELTOL[i]}.log"
        export SUNLOGGER_ERROR_FILENAME="${LOG_FILE}"
        export SUNLOGGER_WARNING_FILENAME="${LOG_FILE}"
        export SUNLOGGER_INFO_FILENAME="${LOG_FILE}"
        export SUNLOGGER_DEBUG_FILENAME="${LOG_FILE}"

        printf "# [Integrator settings]\n" > "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
        printf "order       = ${ORDERS}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
        printf "dt          = 0.0\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
        printf "dt_partition_1 = 0.0\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
        printf "dt_partition_2 = 0.0\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
        printf "type        = ${INTEGRATION_TYPE}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
        printf "final_time  = ${FINAL_TIME}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
        printf "rel_tol     = ${RELTOL[i]}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
        printf "postprocess = ${POSTPROCESS}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
        printf "use_lookup  = ${LOOKUP_FLAG}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
        printf "use_zero_mur = ${USE_ZERO_MUR}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"

        printf "\n# [Save settings]\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
        printf "steps    = ${STEPS_PER_OUTPUT}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
        printf "filename = ${OUTPUT_FILE}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"

        printf "\n# [Initial conditions]\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
        printf "ic_file   = ${IC_FILE}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
        printf "num_cases = ${NUM_CASES}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
        printf "case_idx  = ${k}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"

        printf "\n# [Process settings]\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
        printf "processes = ${PROCESSES}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
        printf "regularize_qsat = ${REGULARIZE_QSAT}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
        printf "regularize_lambdar = ${REGULARIZE_LAMBDAR}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
        printf "qsmall = ${QSMALL}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
        printf "epsilon_qsat_fac = ${EPSILON_QSAT_FAC}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"

        "${RAINSHAFT_EXE}" --i "settings_${SETTINGS_NAME}.ini"
        # ./rainshaft --order ${ORDERS[k]} --dt ${TIMESTEPS[i]} --type $INTEGRATION_TYPE --simname $SIMULATION_NAME
    done
done
