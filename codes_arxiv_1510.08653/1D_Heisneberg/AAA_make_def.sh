#!/bin/bash

L=6
#L=10
#L=14

#p=0.0
p=1.0
#p=2.0

./make_coulombinter.py -L ${L}
./make_exchange.py -L ${L}
./make_greenone.py -L ${L}
./make_greentwo.py -L ${L}
./make_gutzwilleridx.py -L ${L}
./make_hund.py -L ${L}
./make_jastrowidx.py -L ${L}
./make_locspn.py -L ${L}
./make_modpara.py -L ${L}
./make_modpara_aft.py -L ${L}
./make_namelist.py
./make_namelist_aft.py
./make_orbitalidx.py -L ${L}
./make_orbitalidxpara.py -L ${L}
./make_qptransidx.py -L ${L}
./make_trans.py -L ${L}
./make_zqp_init_fij.py -L ${L} -p ${p}
#./make_zqp_init_fij_uRVB.py -L ${L}
