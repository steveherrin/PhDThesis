#!/bin/bash

N=999

#echo "DID YOU ./movestuf???"
#cat

for I in `seq 0 ${N}`; do
  NUM=`printf %03i ${I}`
  #bjobs | grep "ns.sh ${NUM}" | awk '{print $1}' |xargs bkill
  bsub -q xlong ./sim_muons.sh ${NUM}
done
