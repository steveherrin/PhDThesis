#!/bin/bash

for RUNNUM in `cat golden_runs_20130213.txt`; do

    # Find how many files there are for the run, spread over all the places
    # they can be (yes sometimes runs are split over directories).
    # The goldenrun_muons script handles finding the actual files if it's
    # told which number to grab.
    NFILES=0
    for ((I=2;I<=5;I+=1)); do
	DIR="/nfs/slac/g/exo_data${I}/exo_data/data/WIPP/root/${RUNNUM}"
	if [[ -d "${DIR}" ]]; then
	    N=`find ${DIR} -name "*.root" -print | wc -l`
	    NFILES=`expr ${NFILES} + ${N}`
	fi
    done
    
    # If we couldn't find that run, exit
    if [[ "${NFILES}" -eq 0 ]]; then
	echo "Couldn't find run ${RUNNUM}"
	continue
    fi

    #echo "Found ${NFILES} files for run ${RUNNUM}"

    for ((FILENUM=0;FILENUM<${NFILES};FILENUM+=1)); do
	bsub -q long ./goldenrun_muons.sh ${RUNNUM} ${FILENUM}
    done

done
