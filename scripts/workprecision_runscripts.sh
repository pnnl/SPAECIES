#!/bin/bash

RAINSHAFT_EXE="/home/dong9/SPAECIES-evaluation-against-e3sm/build/rainshaft/rainshaft"
RAINSHAFT_DIR="/home/dong9/SPAECIES-evaluation-against-e3sm/build/rainshaft"
SAVE_DIR="/home/dong9/SPAECIES-evaluation-against-e3sm/build/rainshaft/results"
IC_FILE="/home/dong9/SPAECIES-evaluation-against-e3sm/build/rainshaft/random_rainshaft_samples.nc"

# simulation length
FINAL_TIME=1800.0

# (inner) time steps to try
TIMESTEPS=(0.1024 0.2048 0.4096 0.8192 1.6384 3.2768 6.5536 13.1072 26.2144 52.4288 104.8576 209.7152 419.4304 838.8608 1677.7216)

# tolerance for adaptive stepper/nls
REL_TOL=1.e-4

# slow time steps to try for MRI methods
SLOW_TIMESTEPS=(-2 -4 -8 -16)

# order of method
ORDERS=(2 3 4 5)

# toggle nudging
NUDGING_FLAG="true"

# number of runs to do for averaging. processing more than 1 run is currently not supported.
NUMRUNS=1

# how many time steps to save data for
STEPS_PER_OUTPUT=-1

# number of cases from IC_FILE to run (-1 for all cases)
NUM_CASES=-1

# E3SM case number to run
CASE_IDX=31

# type of integration (options: explicit, imex, mri)
INTEGRATION_TYPE="mri"

# name for this collection of simulations. to be used in plot_workprecision.py to gather the data
SIMULATION_NAME="mri"

for kk in $(seq 1 $NUMRUNS) 
do
    printf -- "---------------------------------- RUN %u ----------------------------------\n" $k
    # loop over requested orders
    for k in $(seq 0 $((${#ORDERS[@]} - 1)))
    do
        # if using MRI, specify slow time steps
        if [[ "$INTEGRATION_TYPE" == "mri" ]]; then
            # loop over slow time steps
            for j in $(seq 0 $((${#SLOW_TIMESTEPS[@]} - 1)))
            do
                # loop over fast time steps
                for i in $(seq 0 $((${#TIMESTEPS[@]} - 1)))
                do
                    OUTPUT_FILE="${SAVE_DIR}/rainshaft_${SIMULATION_NAME}_order${ORDERS[k]}_dtslow${SLOW_TIMESTEPS[j]}_dt${TIMESTEPS[i]}.nc"

                    printf "# [Integrator settings]\n" > "${RAINSHAFT_DIR}/settings_i.ini"
                    printf "order      = ${ORDERS[k]}\n" >> "${RAINSHAFT_DIR}/settings_i.ini"
                    printf "dt         = ${TIMESTEPS[i]}\n" >> "${RAINSHAFT_DIR}/settings_i.ini"
                    printf "dt_slow    = ${SLOW_TIMESTEPS[j]}\n" >> "${RAINSHAFT_DIR}/settings_i.ini"
                    printf "type       = ${INTEGRATION_TYPE}\n" >> "${RAINSHAFT_DIR}/settings_i.ini"
                    printf "final_time = ${FINAL_TIME}\n" >> "${RAINSHAFT_DIR}/settings_i.ini"
                    printf "rel_tol    = ${REL_TOL}\n" >> "${RAINSHAFT_DIR}/settings_i.ini"

                    printf "\n# [Save settings]\n" >> "${RAINSHAFT_DIR}/settings_i.ini"
                    printf "steps    = ${STEPS_PER_OUTPUT}\n" >> "${RAINSHAFT_DIR}/settings_i.ini"
                    printf "filename = ${OUTPUT_FILE}\n" >> "${RAINSHAFT_DIR}/settings_i.ini"

                    printf "\n# [Initial conditions]\n" >> "${RAINSHAFT_DIR}/settings_i.ini"
                    printf "ic_file   = ${IC_FILE}\n" >> "${RAINSHAFT_DIR}/settings_i.ini"
                    printf "num_cases = ${NUM_CASES}\n" >> "${RAINSHAFT_DIR}/settings_i.ini"
                    printf "case_idx  = ${CASE_IDX}\n" >> "${RAINSHAFT_DIR}/settings_i.ini"

                    printf "\n# [Process settings]\n" >> "${RAINSHAFT_DIR}/settings_i.ini"
                    printf "nudging = ${NUDGING_FLAG}\n" >> "${RAINSHAFT_DIR}/settings_i.ini"

                    ${RAINSHAFT_EXE} --i "settings_i.ini"
                    # ./rainshaft --order ${ORDERS[k]} --dt_fast ${TIMESTEPS[i]} --dt ${SLOW_TIMESTEPS[j]} --type $INTEGRATION_TYPE --simname $SIMULATION_NAME
                done
            done
        else # otherwise, only loop over the single time steps
            for i in $(seq 0 $((${#TIMESTEPS[@]} - 1)))
            do
                OUTPUT_FILE="${SAVE_DIR}/rainshaft_${SIMULATION_NAME}_order${ORDERS[k]}_dt${TIMESTEPS[i]}.nc"

                printf "# [Integrator settings]\n" > "${RAINSHAFT_DIR}/settings_i.ini"
                printf "order      = ${ORDERS[k]}\n" >> "${RAINSHAFT_DIR}/settings_i.ini"
                printf "dt         = ${TIMESTEPS[i]}\n" >> "${RAINSHAFT_DIR}/settings_i.ini"
                printf "type       = ${INTEGRATION_TYPE}\n" >> "${RAINSHAFT_DIR}/settings_i.ini"
                printf "final_time = ${FINAL_TIME}\n" >> "${RAINSHAFT_DIR}/settings_i.ini"
                printf "rel_tol    = ${REL_TOL}\n" >> "${RAINSHAFT_DIR}/settings_i.ini"

                printf "\n# [Save settings]\n" >> "${RAINSHAFT_DIR}/settings_i.ini"
                printf "steps    = ${STEPS_PER_OUTPUT}\n" >> "${RAINSHAFT_DIR}/settings_i.ini"
                printf "filename = ${OUTPUT_FILE}\n" >> "${RAINSHAFT_DIR}/settings_i.ini"

                printf "\n# [Initial conditions]\n" >> "${RAINSHAFT_DIR}/settings_i.ini"
                printf "ic_file   = ${IC_FILE}\n" >> "${RAINSHAFT_DIR}/settings_i.ini"
                printf "num_cases = ${NUM_CASES}\n" >> "${RAINSHAFT_DIR}/settings_i.ini"
                printf "case_idx  = ${CASE_IDX}\n" >> "${RAINSHAFT_DIR}/settings_i.ini"

                printf "\n# [Process settings]\n" >> "${RAINSHAFT_DIR}/settings_i.ini"
                printf "nudging = ${NUDGING_FLAG}\n" >> "${RAINSHAFT_DIR}/settings_i.ini"

                ${RAINSHAFT_EXE} --i "settings_i.ini"
                # ./rainshaft --order ${ORDERS[k]} --dt ${TIMESTEPS[i]} --type $INTEGRATION_TYPE --simname $SIMULATION_NAME
            done
        fi
    done
done