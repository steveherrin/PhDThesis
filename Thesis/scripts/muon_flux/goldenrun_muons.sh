#!/bin/bash

MAXEVENTS=-1
#PRINTMODULO=$((${MAXEVENTS}/5))
PRINTMODULO=1000

ulimit -c 1024
finalRC=1
set -e

EXODIR=/nfs/slac/g/exo
#EXOBASE=${EXODIR}/sherrin/software_muonsim/exoout
#EXOBASE=${EXODIR}/software/hudson/builds-rhel5/trunk
#EXOBASE=${EXODIR}/software/builds/trunk
EXOBASE=${EXODIR}/sherrin/software_muonsim/exoout
#EXOBASE=/tmp/sherrin/exoout

RUNNUM=$1
FILENUM=$2

# Find where the file actually lives
EXODATANUM=0
for ((I=2;I<=5;I+=1)); do
    INPUT_DIR=/nfs/slac/g/exo_data${I}/exo_data/data/WIPP/root
    FILE_NAME=`printf "run%08i-%03i.root" ${RUNNUM} ${FILENUM}`
    INPUT_FILE=${INPUT_DIR}/${RUNNUM}/${FILE_NAME}
    if [[ -f "${INPUT_FILE}" ]]; then
	EXODATANUM=${I}
	break
    fi
done

# If we couldn't find that run, exit
if [[ "${EXODATANUM}" -eq 0 ]]; then
    echo "Couldn't find ${INPUT_FILE}."
    exit 1
fi

INPUT_DIR=/nfs/slac/g/exo_data${EXODATANUM}/exo_data/data/WIPP/root
FILE_NAME=`printf "run%08i-%03i.root" ${RUNNUM} ${FILENUM}`
INPUT_FILE=${INPUT_DIR}/${RUNNUM}/${FILE_NAME}

OUTPUT_DIR=/nfs/slac/g/exo_data3/sherrin/muon_flux/golden_runs
OUTPUT_FILE=${OUTPUT_DIR}/muon_${RUNNUM}_${FILENUM}.root

export SCRATCH_DIR=/scratch/sherrin/muon_${RUNNUM}
mkdir -p ${SCRATCH_DIR}
gotEXIT()
{
    rm -rf ${SCRATCH_DIR}
    exit $finalRC
}
trap gotEXIT EXIT

source ${EXOBASE}/setup.sh

cat > ${SCRATCH_DIR}/myexo.exo << EOF
use input alphanoisetag muontrack evsel toutput

/input/file ${INPUT_FILE}

/alphanoisetag/verbose false

#/muontrack/nwirecounts 80.0
#/muontrack/napdcounts 5000.0
#/muontrack/minlightdt 0.0
#/muontrack/driftspeed 1.7
/muontrack/indvhists false
/muontrack/verbose true

/evsel/muon_include 1
/evsel/filter 1

/toutput/writeSignals false
/toutput/file ${SCRATCH_DIR}/output.root

printmodulo ${PRINTMODULO}
maxevents ${MAXEVENTS}
begin
exit
EOF

#cp -pv ${SCRATCH_DIR}/myexo.exo ${OUTPUT_DIR}/myexo.exo

EXOAnalysis ${SCRATCH_DIR}/myexo.exo

mkdir -p ${OUTPUT_DIR}
cp -pv ${SCRATCH_DIR}/output.root ${OUTPUT_FILE}

rm -rf ${SCRATCH_DIR}

finalRC=0
