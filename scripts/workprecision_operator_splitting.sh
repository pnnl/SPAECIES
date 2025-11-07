#!/bin/bash

IC_FILE="/home/dong9/SPAECIES-input-data/random_rainshaft_samples_12mo.nc"
RAINSHAFT_EXE="/home/dong9/SPAECIES-settings-for-paper/build/rainshaft/rainshaft"
RAINSHAFT_DIR="/home/dong9/SPAECIES-settings-for-paper/build/rainshaft"
SAVE_DIR="/home/dong9/postprocessing-SPAECIES/results/final-runs-sept23/operator_splitting"

# simulation length
FINAL_TIME=300

# fixed steps to try
TIMESTEPS=(0.0128 0.0256 0.0512 0.1024 0.2048 0.4096 0.8192 1.6384 3.2768 6.5536 13.1072 26.2144 52.4288 104.8576 209.7152 300.0)

RELTOL=1.e99

# order of method
ORDERS=(1 2 3)

# toggle nudging
NUDGING_FLAG="false"

# how many time steps to save data for
STEPS_PER_OUTPUT=-1

# number of cases from IC_FILE to run (-1 for all cases)
NUM_CASES=-1

# E3SM case number to run
CASE_IDX=

# toggle limiting
POSTPROCESS="true"

# toggle lookup tables
LOOKUP_FLAG="true"

# toggle q_sat_dry regularization
REGULARIZE_QSAT="true"
REGULARIZE_LAMBDAR="true"
EPSILON_QSAT_FACS=(1e-5 1e-2)
QSMALLS=(1e-14 1e-10)

# use CFL substepping for sedimentation
CFL_SUBSTEP="true"

# type of integration (options: explicit, imex, mri)
INTEGRATION_TYPES=("splitting")

cd ${RAINSHAFT_DIR}

for kk in $(seq 0 $((${#INTEGRATION_TYPES[@]} - 1)))
do
    # loop over requested orders
    for k in $(seq 0 $((${#ORDERS[@]} - 1)))
    do
        for i in $(seq 0 $((${#TIMESTEPS[@]} - 1)))
        do
            for ii in $(seq 0 $((${#QSMALLS[@]} - 1)))
            do
                for j in $(seq 0 $((${#EPSILON_QSAT_FACS[@]} - 1)))
                do 
                    if [[ $REGULARIZE_LAMBDAR == "true" && $REGULARIZE_QSAT == "true" ]]
                    then
                        # name for this collection of simulations. to be used in plot_workprecision.py to gather the data
                        SIMULATION_NAME="${INTEGRATION_TYPES[kk]}_regularized"
                        OUTPUT_FILE="${SAVE_DIR}/rainshaft_${SIMULATION_NAME}_finaltime${FINAL_TIME}_qsmall${QSMALLS[ii]}_epsilonqsat${EPSILON_QSAT_FACS[j]}_order${ORDERS[k]}_dt${TIMESTEPS[i]}.nc"
                    else 
                        SIMULATION_NAME="${INTEGRATION_TYPES[kk]}_unregularized"
                        OUTPUT_FILE="${SAVE_DIR}/rainshaft_${SIMULATION_NAME}_finaltime${FINAL_TIME}_order${ORDERS[k]}_dt${TIMESTEPS[i]}.nc"
                    fi

                    printf "# [Integrator settings]\n" > "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
                    printf "order       = ${ORDERS[k]}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
                    printf "dt          = ${TIMESTEPS[i]}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
                    printf "dt_partition_1 = 0.0\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
                    printf "dt_partition_2 = 1.e99\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
                    printf "type        = ${INTEGRATION_TYPES[kk]}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
                    printf "final_time  = ${FINAL_TIME}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
                    printf "rel_tol     = ${RELTOL}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
                    printf "postprocess = ${POSTPROCESS}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
                    printf "use_lookup  = ${LOOKUP_FLAG}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
                    printf "cfl_substep = ${CFL_SUBSTEP}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"

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
                    printf "regularize_lambdar = ${REGULARIZE_LAMBDAR}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
                    printf "qsmall = ${QSMALLS[ii]}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"
                    printf "epsilon_qsat_fac = ${EPSILON_QSAT_FACS[j]}\n" >> "${RAINSHAFT_DIR}/settings_${SETTINGS_NAME}.ini"

                    ${RAINSHAFT_EXE} --i "settings_${SETTINGS_NAME}.ini"
                done
            done
        done
    done
done