#!/bin/bash

make -j

# (inner) time steps to try
TIMESTEPS=(0.1024 0.2048 0.4096 0.8192 1.6384 3.2768 6.5536 13.1072 26.2144 52.4288 104.8576 209.7152 419.4304 838.8608 1677.7216)

# slow time steps to try for MRI methods
SLOW_TIMESTEPS=(-2 -4 -8 -16)
# SLOW_TIMESTEPS=-2

# order of method
ORDERS=(2 3 4 5)

# number of runs to do for averaging. processing more than 1 run is currently not supported.
NUMRUNS=1

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
                    ./rainshaft --order ${ORDERS[k]} --dt ${TIMESTEPS[i]} --dt_slow ${SLOW_TIMESTEPS[j]} --type $INTEGRATION_TYPE --simname $SIMULATION_NAME
                done
            done
        else # otherwise, only loop over fast time steps
            for i in $(seq 0 $((${#TIMESTEPS[@]} - 1)))
            do
                ./rainshaft --order ${ORDERS[k]} --dt ${TIMESTEPS[i]} --type $INTEGRATION_TYPE --simname $SIMULATION_NAME
            done
        fi
    done
done