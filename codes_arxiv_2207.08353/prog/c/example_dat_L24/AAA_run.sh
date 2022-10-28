#!/bin/bash

nmax=64

for L in \
24
do

for n in \
`seq 1 $((nmax))`
do
  t=`echo ${L} ${n} ${nmax} | awk '{printf("%.10f",2.0*$1*$2/$3)}'`
  echo ${L} ${n} ${t}
  ../main -L ${L} -t ${t}
done

done
