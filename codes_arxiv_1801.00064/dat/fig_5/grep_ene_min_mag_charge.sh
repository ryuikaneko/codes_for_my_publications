#!/bin/bash

L=24
#L=120

for V in \
4.00
#$(seq 0 0.25 3.01)
do
  output=dat_V${V}
  echo -ne "" > ${output}
  for U in \
  $(seq 0 0.25 10.01)
#  $(seq 0 0.25 6.01)
  do
    echo -ne "${U} ${V} " >> ${output}
    awk \
      '$1=='${U}'{print $4,\
      0.5*($8-$14),0.5*($9-$15),0.5*($10-$16),0.5*($11-$17),0.5*($12-$18),0.5*($13-$19),\
      ($8+$14),($9+$15),($10+$16),($11+$17),($12+$18),($13+$19)\
      }' \
      ./catdat_L${L}V${V}* | sort -g | head -n 1 >> ${output}
  done
done
