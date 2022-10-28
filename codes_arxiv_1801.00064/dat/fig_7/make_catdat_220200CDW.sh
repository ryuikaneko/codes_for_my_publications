#!/bin/bash

L=24
#L=240
name_state=220200CDW

for U in \
1.00
#`seq -f%.2f 0.00 0.25 6.001`
do
  catdat=catdat_L${L}U${U}_${name_state}
#  catdatgap=catdatgap_L${L}U${U}_${name_state}
  echo -ne "# U V Niter " > ${catdat}
  echo -ne " ene ene_0 ene_U \n" >> ${catdat}
#  echo -ne "" > ${catdatgap}

  for V in \
  `seq -f%.2f 0.00 0.25 3.001`
  do
    dat=dat_L${L}U${U}V${V}_${name_state}

    echo -ne "${U} ${V} " >> ${catdat}
    grep "converged" -A1 ${dat} | tail -n 1 >> ${catdat}

#    echo -ne "${U} ${V} " >> ${catdatgap}
#    grep Egap ${dat} | sed 's/.*=//g' >> ${catdatgap}
  done
done
