#!/bin/bash

for name in \
HF_2sub \
HF_6sub
do
  gcc -O3 ${name}.c -lm -llapack -o ${name}.out -Wall
done

rm -f  *~
