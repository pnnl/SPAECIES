#!/bin/bash

# RAINSHAFT_EXE="/home/dong9/SPAECIES-final-version-for-paper/build_production/rainshaft/rainshaft_allprocesses_fixedstep_implicit"
# RAINSHAFT_DIR="/home/dong9/SPAECIES-final-version-for-paper/build_production/rainshaft"
# SAVE_DIR="/home/dong9/postprocessing-SPAECIES/results/final-runs-mar21/fixedstep_lookupon_limiteroff"
IC_FILE="/home/dong9/SPAECIES-final-version-for-paper/data/random_rainshaft_samples_12mo.nc"

RAINSHAFT_EXE="/home/dong9/SPAECIES-regularization/build/rainshaft/rainshaft_regularized"
RAINSHAFT_DIR="/home/dong9/SPAECIES-regularization/build/rainshaft"
SAVE_DIR="/home/dong9/postprocessing-SPAECIES/results/final-runs-sept23/reference_solutions"

# simulation length
FINAL_TIME=300

# fixed steps to try
TIMESTEP=0.0
RELTOLS=(1e-12)

# order of method
ORDERS=(3)

# qsmalls 
QSMALLS=(1e-5 1e-6 1e-7 1e-8 1e-9 1e-10 1e-11 1e-12 1e-13 1e-14)

# epsilon_qsat
EPSILON_QSAT_FACS=(1e-5 1e-4 1e-3 1e-2 1e-1 0.25 1.0)

# toggle nudging
NUDGING_FLAG="false"

# number of runs to do for averaging. processing more than 1 run is currently not supported.
NUMRUNS=1

# how many time steps to save data for
STEPS_PER_OUTPUT=-1

# number of cases from IC_FILE to run (-1 for all cases)
NUM_CASES=-1

# E3SM case number to run
CASE_IDX=

# toggle limiting
POSTPROCESS="false"

# toggle lookup tables
LOOKUP_FLAG="false"

# toggle q_sat_dry regularization
REGULARIZE_QSAT="true"

# type of integration (options: explicit, imex, mri)
INTEGRATION_TYPE="explicit"

# name for this collection of simulations. to be used in plot_workprecision.py to gather the data
SIMULATION_NAME="allprocesses_allcolumns"

# SETTINGS_NAME="imp_skipcrashes_abstol0_reltol1e-12_realjacobian_kvaerno423_fixedstep_lookupon_limiteroff"
SETTINGS_NAME="imp_skipcrashes_lookupon_limiteroff"
# SETTINGS_NAME="imp_skipcrashes_abstol0_reltol1e-12_FDjacobian_fixedstep_lookupon_limiteroff"

for kk in $(seq 1 $NUMRUNS) 
do
    printf -- "---------------------------------- RUN %u ----------------------------------\n" $k
    # loop over requested orders
    for k in $(seq 0 $((${#ORDERS[@]} - 1)))
    do
        for i in $(seq 0 $((${#RELTOLS[@]} - 1)))
        do
            for j in $(seq 0 $((${#QSMALLS[@]} - 1)))
            do
                for m in $(seq 0 $((${#EPSILON_QSAT_FACS[@]} - 1)))
                do
                    OUTPUT_FILE="${SAVE_DIR}/rainshaft_${SIMULATION_NAME}_tf${FINAL_TIME}_qsmall${QSMALLS[j]}_EPSILONQSATFAC${EPSILON_QSAT_FACS[m]}_order${ORDERS[k]}_reltol${RELTOLS[i]}_every1sec_reference.nc"

                    printf "# [Integrator settings]\n" > "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
                    printf "order       = ${ORDERS[k]}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
                    printf "dt          = ${TIMESTEP}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
                    printf "type        = ${INTEGRATION_TYPE}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
                    printf "final_time  = ${FINAL_TIME}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
                    printf "rel_tol     = ${RELTOLS[i]}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
                    printf "postprocess = ${POSTPROCESS}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
                    printf "use_lookup  = ${LOOKUP_FLAG}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"

                    printf "\n# [Save settings]\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
                    printf "steps    = ${STEPS_PER_OUTPUT}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
                    printf "filename = ${OUTPUT_FILE}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"

                    printf "\n# [Initial conditions]\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
                    printf "ic_file   = ${IC_FILE}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
                    printf "num_cases = ${NUM_CASES}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
                    printf "case_idx  = ${CASE_IDX}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"

                    printf "\n# [Process settings]\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
                    printf "nudging = ${NUDGING_FLAG}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
                    printf "regularize_qsat = ${REGULARIZE_QSAT}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
                    printf "qsmall = ${QSMALLS[j]}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
                    printf "epsilon_qsat_fac = ${EPSILON_QSAT_FACS[m]}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"

                    ${RAINSHAFT_EXE} --i "settings_${SETTINGS_NAME}.ini"
                done
            done
        done
    done
done