#!/bin/bash

#L=12
L=24
#L=120
UpV=0
intkmax=${L}
Niter=2000
mixing=0.9

declare -a dir_array=("normal_metal" "AF_insulator" "CO_AF_insulator" "FM_metal")
declare -a state_array=("0" "2" "4" "3")

for num in \
0 1 2 3
do
  dir=dat_${dir_array[${num}]}
  state=${state_array[${num}]}
  echo ${dir}
  mkdir -p ${dir}

  prog=../../prog/HF_4sub.out
  catdat=./${dir}/catdat_L${L}UpV${UpV}
#  catdatgap=./${dir}/catdatgap_L${L}UpV${UpV}
  echo -ne "# U Up=V Niter " > ${catdat}
  echo -ne " ene ene_0 ene_U HalfN DeltaNc[0] DeltaNc[1] DeltaNf[0] DeltaNf[1] Chi[0] Chi[1] DeltaChi[0] DeltaChi[1]\n" >> ${catdat}
#  echo -ne "" > ${catdatgap}

  for U in \
  `seq -f%.1f 0.0 0.5 3.51` \
  `seq -f%.1f 4.0 0.1 6.01` \
  `seq -f%.1f 6.5 0.5 10.01`
  do
    dat=./${dir}/dat_L${L}U${U}UpV${UpV}
    ./${prog} -U ${U} -P ${UpV} -V ${UpV} -k ${intkmax} -i ${Niter} -m ${mixing} -s ${state} > ${dat}
#    ./${prog} -U ${U} -P ${UpV} -V ${UpV} -k ${intkmax} -i ${Niter} -m ${mixing} -s ${state} -b > ${dat}
  done

  for U in \
  `seq -f%.1f 0.0 0.5 3.51` \
  `seq -f%.1f 4.0 0.1 6.01` \
  `seq -f%.1f 6.5 0.5 10.01`
  do
    dat=./${dir}/dat_L${L}U${U}UpV${UpV}
    echo -ne "${U} ${UpV} " >> ${catdat}
    grep "converged" -A1 ${dat} | tail -n 1 >> ${catdat}
  done

#  for U in \
#  `seq -f%.1f 0.0 0.5 3.51` \
#  `seq -f%.1f 4.0 0.1 6.01` \
#  `seq -f%.1f 6.5 0.5 10.01`
#  do
#    dat=./${dir}/dat_L${L}U${U}UpV${UpV}
#    echo -ne "${U} ${UpV} " >> ${catdatgap}
#    grep Egap ${dat} | sed 's/.*=//g' >> ${catdatgap}
#  done
done
