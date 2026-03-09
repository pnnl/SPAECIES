#!/bin/bash

IC_FILE="/home/dong9/SPAECIES-input-data/random_rainshaft_samples_12mo.nc"
RAINSHAFT_EXE="/home/dong9/SPAECIES-settings-for-paper/build/rainshaft/rainshaft_sedonly_logging"
RAINSHAFT_DIR="/home/dong9/SPAECIES-settings-for-paper/build/rainshaft"
SAVE_DIR="/home/dong9/postprocessing-SPAECIES/results/final-runs-nov7/logfiles/adaptive_sedonly"

# simulation length
FINAL_TIME=1000

# tolerance for adaptive stepper/nls
RELTOL=(1.e-1 1.e-3 1.e-5 1.e-7)

# order of method
ORDERS=2

# toggle nudging
NUDGING_FLAG="false"

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

        export SUNLOGGER_ERROR_FILENAME=/home/dong9/postprocessing-SPAECIES/results/final-runs-nov7/logfiles/adaptive_sedonly/icase${k}_reltol${RELTOL[i]}.log
        export SUNLOGGER_WARNING_FILENAME=/home/dong9/postprocessing-SPAECIES/results/final-runs-nov7/logfiles/adaptive_sedonly/icase${k}_reltol${RELTOL[i]}.log
        export SUNLOGGER_INFO_FILENAME=/home/dong9/postprocessing-SPAECIES/results/final-runs-nov7/logfiles/adaptive_sedonly/icase${k}_reltol${RELTOL[i]}.log
        export SUNLOGGER_DEBUG_FILENAME=/home/dong9/postprocessing-SPAECIES/results/final-runs-nov7/logfiles/adaptive_sedonly/icase${k}_reltol${RELTOL[i]}.log

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

        printf "\n# [Save settings]\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
        printf "steps    = ${STEPS_PER_OUTPUT}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
        printf "filename = ${OUTPUT_FILE}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"

        printf "\n# [Initial conditions]\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
        printf "ic_file   = ${IC_FILE}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
        printf "num_cases = ${NUM_CASES}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
        printf "case_idx  = ${k}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"

        printf "\n# [Process settings]\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
        printf "nudging = ${NUDGING_FLAG}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
        printf "regularize_qsat = ${REGULARIZE_QSAT}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
        printf "regularize_lambdar = ${REGULARIZE_LAMBDAR}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
        printf "qsmall = ${QSMALL}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
        printf "epsilon_qsat_fac = ${EPSILON_QSAT_FAC}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"

        ${RAINSHAFT_EXE} --i "settings_${SETTINGS_NAME}.ini"
        # ./rainshaft --order ${ORDERS[k]} --dt ${TIMESTEPS[i]} --type $INTEGRATION_TYPE --simname $SIMULATION_NAME
    done
done