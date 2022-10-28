#!/bin/bash

for th in \
0.05 0.10
#`seq -f %.2f 0.00 0.05 1.00`
do
  echo -ne "${th} "
  if [ -e ../dat_t${th} ]; then
    tail -n 434 ../dat_t${th} | head -n 288 > dat
    ./BBB_calc_quad.py > dat_sq_t${th}
  fi
  echo -ne "\n"
done

rm dat
