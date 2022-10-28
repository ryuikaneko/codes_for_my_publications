#!/bin/bash

L=24
#L=240
Niter=500
mixing=0.5

prog=../../prog/HF_2sub.out
init_state=3
name_state=TMI

for U in \
1.00
#`seq -f%.2f 0.00 0.25 6.001`
do
  for V in \
  `seq -f%.2f 0.00 0.25 3.001`
  do
    dat=dat_L${L}U${U}V${V}_${name_state}
    ./${prog} -U ${U} -V ${V} -k ${L} -i ${Niter} -m ${mixing} -s ${init_state} > ${dat}
#    ./${prog} -U ${U} -V ${V} -k ${L} -i ${Niter} -m ${mixing} -s ${init_state} -b > ${dat}
  done
done
