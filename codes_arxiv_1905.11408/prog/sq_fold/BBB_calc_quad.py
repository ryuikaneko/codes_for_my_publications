#!/usr/bin/env python
# -*- coding: utf-8 -*-

## tail -864 dat_t-0.50 > dat

from __future__ import print_function

import numpy as np
import scipy as scipy

sizeL = 12
#sizeL = 3
sizeK = sizeL*2
sizeOrb = 2
sizeNs = sizeL*sizeL*sizeOrb

dat_ind = np.loadtxt('dat', usecols=[0,1,2], dtype='int')
dat_val = np.loadtxt('dat', usecols=[3,4,5,6,7], dtype='float')

dist = np.zeros((sizeL,sizeL,sizeOrb,2),dtype='float')
spin = np.zeros((sizeL,sizeL,sizeOrb,3),dtype='float')
quad = np.zeros((sizeL,sizeL,sizeOrb,5),dtype='float')

sizeNs = sizeL*sizeL*sizeOrb
for i in range(sizeNs):
    ind = dat_ind[i]
    val = dat_val[i]
    dist[ind[0],ind[1],ind[2],0] = val[0]
    dist[ind[0],ind[1],ind[2],1] = val[1]
    spin[ind[0],ind[1],ind[2],0] = val[2]
    spin[ind[0],ind[1],ind[2],1] = val[3]
    spin[ind[0],ind[1],ind[2],2] = val[4]
    quad[ind[0],ind[1],ind[2],0] = (3.0*val[4]*val[4] - 1.0)/np.sqrt(3.0)
    quad[ind[0],ind[1],ind[2],1] = val[2]*val[2] - val[3]*val[3]
    quad[ind[0],ind[1],ind[2],2] = 2.0*val[2]*val[3]
    quad[ind[0],ind[1],ind[2],3] = 2.0*val[3]*val[4]
    quad[ind[0],ind[1],ind[2],4] = 2.0*val[4]*val[2]

quad2=np.sum(np.sum(np.sum(quad,axis=0),axis=0),axis=0)
quad2=quad2*quad2/sizeNs/sizeNs
quad2=np.sum(quad2,axis=0)
print ("# quad2: ",quad2)

for ix in range(sizeL):
    for iy in range(sizeL):
        for iorb in range(sizeOrb):
            print(ix,iy,iorb,
                dist[ix,iy,iorb,0],dist[ix,iy,iorb,1],
                spin[ix,iy,iorb,0],spin[ix,iy,iorb,1],spin[ix,iy,iorb,2],
                quad[ix,iy,iorb,0],quad[ix,iy,iorb,1],quad[ix,iy,iorb,2],
                quad[ix,iy,iorb,3],quad[ix,iy,iorb,4])

expiqr = [[[[[\
    np.exp(1.0j*2.0*np.pi/sizeL*(\
    + (qx-sizeK)/np.sqrt(3.0) * dist[ix,iy,iorb,0]\
    + (-(qx-sizeK)-(qy-sizeK)*2.0)/3.0 * dist[ix,iy,iorb,1])\
    )/sizeNs \
    for iorb in range(sizeOrb)] \
    for iy in range(sizeL)] \
    for ix in range(sizeL)] \
    for qy in range(sizeK*2+1)] \
    for qx in range(sizeK*2+1)]

sq_qx_qy = np.tensordot(spin,expiqr,axes=((0,1,2),(2,3,4)))
sqsmq = np.einsum('ijk,ijk->ijk',sq_qx_qy,np.conjugate(sq_qx_qy))
sqsmq = np.sum(sqsmq,axis=0)

qq_qx_qy = np.tensordot(quad,expiqr,axes=((0,1,2),(2,3,4)))
qqqmq = np.einsum('ijk,ijk->ijk',qq_qx_qy,np.conjugate(qq_qx_qy))
qqqmq = np.sum(qqqmq,axis=0)

print ()
print ()
sqargmax=np.argmax(np.real(sqsmq))
sqarg1=sqargmax/(sizeK*2+1)
sqarg2=sqargmax%(sizeK*2+1)
print ("# s(q) argmax max: ",sqarg1,sqarg2,np.max(np.real(sqsmq)))
qqargmax=np.argmax(np.real(qqqmq))
qqarg1=qqargmax/(sizeK*2+1)
qqarg2=qqargmax%(sizeK*2+1)
print ("# q(q) argmax max: ",qqarg1,qqarg2,np.max(np.real(qqqmq)))
for qx in range(sizeK*2+1):
    for qy in range(sizeK*2+1):
        print (qx,qy,
            + (qx-sizeK)/np.sqrt(3.0) * 2.0*np.pi/sizeL,
            + (-(qx-sizeK)-(qy-sizeK)*2.0)/3.0 * 2.0*np.pi/sizeL,
            np.real(sqsmq[qx,qy]),
            np.real(qqqmq[qx,qy])
            )
    print ()
