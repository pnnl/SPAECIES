#!/usr/bin/env bash

RAINSHAFT_EXE="./rainshaft/rainshaft"

DATA_LOC="/home/sant720/data/e3sm-microphysics"

INITIAL_CONDITION="${DATA_LOC}/random_rainshaft_samples.nc"

if [ "$INITIAL_CONDITION" != "adiabatic" ]; then
    NUM_CASES="100"
else
    NUM_CASES=""
fi

SIM_LEN=300.

NUDGE_OPT="FALSE"

METHOD_TYPE="original"

if [ "$METHOD_TYPE" != "original" ]; then
    ORDER=2
    STEPS_PER_OUTPUT=0
else
    ORDER=""
    STEPS_PER_OUTPUT=""
fi

for dt in 1 2 4 8 16 32 64 128 256 300; do
    if [ "$METHOD_TYPE" == "mri" ]; then
        FAST_FAC=5
        dt_fast=`bc -l <<< "${dt}/${FAST_FAC}"`
        fast_suffix="_${FAST_FAC}xfast"
    else
        dt_fast=""
        fast_suffix=""
    fi
    OUTPUT_FILE="./rainshaft_${METHOD_TYPE}${ORDER}_${dt}s${fast_suffix}.nc"
    cmd="$RAINSHAFT_EXE $INITIAL_CONDITION $NUM_CASES $SIM_LEN $NUDGE_OPT $METHOD_TYPE $dt_fast $dt $ORDER $STEPS_PER_OUTPUT $OUTPUT_FILE"
    echo "$cmd"
    eval "$cmd"
done
