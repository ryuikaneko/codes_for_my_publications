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
echo "for fig 4"
cd ./fig_4
for i in \
make_dat_*
do
  echo "running a mean field code: ${i}"
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
