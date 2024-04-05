#!/bin/bash

make -j

# list of slow time step factors: dt_slow = SLOWFAC * dt_fast
SLOWFAC=(4 10 100 300)
NUMRUNS=0
for k in $(seq 0 $NUMRUNS) # ensemble runs
do
	printf -- "---------------------------------- RUN %u ----------------------------------\n" $k
	for j in $(seq 0 ${#SLOWFAC[@]}) # number of different slow time step factors
	do
		for i in {10..10} # # number of fast time steps in convergence plot: time step is (1.28e-5) * pow(2,i)
		do
		    ./rainshaft $i ${SLOWFAC[j]} $k
		done
	done
done