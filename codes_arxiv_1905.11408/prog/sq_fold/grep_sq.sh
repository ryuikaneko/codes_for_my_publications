#!/bin/bash

if [ -e ../dat ]; then
  tail -n 434 ../dat | head -n 288 > dat
#  tail -n 29 ../dat | head -n 18 > dat
  ./BBB_calc_quad.py > dat_sq
fi

rm dat
