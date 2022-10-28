#!/bin/bash
#PBS -q batch
#PBS -l nodes=gr5
#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
#PBS -N test

source scl_source enable python27

export OMP_NUM_THREADS=$PBS_NP
cd $PBS_O_WORKDIR
date

th=0.02

Gp=0.00
./main.out -t ${th} -p ${Gp} > dat

date

cd sq_fold

if [ -e ../dat ]; then
  tail -n 434 ../dat | head -n 288 > dat
#  tail -n 29 ../dat | head -n 18 > dat
#  ./BBB_calc_quad.py > dat_sq
  python2.7 BBB_calc_quad.py > dat_sq
fi

rm dat

cd -
