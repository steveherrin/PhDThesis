#!/bin/bash

MAXEVENTS=2000
PRINTMODULO=$((${MAXEVENTS}/5))
#PRINTMODULO=10

ulimit -c 1024
finalRC=1
set -e

EXODIR=/nfs/slac/g/exo
#EXOBASE=${EXODIR}/sherrin/software_muonsim/exoout
#EXOBASE=${EXODIR}/software/hudson/builds-rhel5/trunk
#EXOBASE=${EXODIR}/software/builds/trunk
EXOBASE=${EXODIR}/sherrin/software_muonsim/exoout
#EXOBASE=/tmp/sherrin/exoout

RUNNUM=$1;

OUTPUT_DIR=/nfs/slac/g/exo_data3/sherrin/muon_flux/simulations
OUTPUT_FILE=${OUTPUT_DIR}/muonsim_${RUNNUM}.root

export SCRATCH_DIR=/scratch/sherrin/muonsim_${RUNNUM}
mkdir -p ${SCRATCH_DIR}
gotEXIT()
{
    rm -rf ${SCRATCH_DIR}
    exit $finalRC
}
trap gotEXIT EXIT

source ${EXOBASE}/setup.sh

cat > ${SCRATCH_DIR}/myexo.exo << EOF
load ${EXOBASE}/plugins/EXOGeant4Module.*
use exosim digitizer alphanoisetag muontrack toutput
#use exosim toutput

/exosim/macro ${SCRATCH_DIR}/mymac.mac
/exosim/initial_seed ${RUNNUM}
/exosim/run_number ${RUNNUM}
/exosim/physics_MuNuclear false
/exosim/physics_Neutron false
/exosim/physics_Optics false
/exosim/SkipEmptyEvents false

#/digitizer/useRealNoise true
#/digitizer/LXeEnergyRes 0.12
/digitizer/driftVelocity 0.171
/digitizer/electronLifetime 2500 microsecond
/digitizer/setTriggerTime 512 microsecond
/digitizer/setDigitizationTime 1024

/alphanoisetag/verbose false

#/muontrack/nwirecounts 80.0
#/muontrack/napdcounts 5000.0
#/muontrack/minlightdt 0.0
#/muontrack/driftspeed 1.7
/muontrack/indvhists false
/muontrack/verbose false

#/rec/pattern_recognition 4
#/rec/filter true
#/rec/drift_velocity_mm_per_ns 0.00180

/toutput/writeSignals false
/toutput/file ${SCRATCH_DIR}/output.root

printmodulo ${PRINTMODULO}
maxevents ${MAXEVENTS}
begin
exit
EOF

cat > ${SCRATCH_DIR}/mymac.mac <<EOF
/event/LXeEventsOnly false
/event/digitizeWires true
/event/digitizeAPDs true
/event/printParticleInfo false
/event/totalEventWindowTime 1024.0 microsecond

/grdm/analogueMC 1

/generator/setGenerator muon

EOF

#cp -pv ${SCRATCH_DIR}/myexo.exo ${OUTPUT_DIR}/myexo.exo
#cp -pv ${SCRATCH_DIR}/mymac.mac ${OUTPUT_DIR}/mymac.mac

EXOAnalysis ${SCRATCH_DIR}/myexo.exo

mkdir -p ${OUTPUT_DIR}
cp -pv ${SCRATCH_DIR}/output.root ${OUTPUT_FILE}

rm -rf ${SCRATCH_DIR}

finalRC=0
