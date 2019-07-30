#!/bin/bash

echo "####"
echo "compiling mean field codes"
cd ./prog
./Compile.sh
cd -
echo "####"

#----

echo "####"
echo "making data files"
cd ./dat
echo "####"

echo "####"
echo "for fig 5"
cd ./fig_5
for i in \
make_dat_*
do
  echo "running a mean field code: ${i}"
  ./${i}
done
for i in \
make_catdat_*
do
  echo "gathering data: ${i}"
  ./${i}
done
echo "finding lowest-energy states"
./grep_ene_min_mag_charge.sh
cd -
echo "####"

echo "####"
echo "for fig 6"
echo "making DOS data"
cd ./fig_6
for i in \
make_dat_*
do
  echo "running a mean field code: ${i}"
  ./${i}
done
cd -
echo "####"

echo "####"
echo "for fig 7"
cd ./fig_7
for i in \
make_dat_*
do
  echo "running a mean field code: ${i}"
  ./${i}
done
for i in \
make_catdat_*
do
  echo "gathering data: ${i}"
  ./${i}
done
cd -
echo "####"

cd ..

#----

echo "####"
echo "making figures"
cd ./fig
for i in \
gnuplot_fig_*.sh
do
  echo "running a gnuplot script: ${i}"
  ./${i}
done
cd -
echo "####"
