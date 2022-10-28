#!/usr/bin/env python

# coding:utf-8
import numpy as np
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='prepare zqp_init.dat')
    parser.add_argument('-L',metavar='L',dest='L',type=int,default=16,help='set L')
    parser.add_argument('-Nsub',metavar='Nsub',dest='Nsub',type=int,default=1,help='set Nsub')
    parser.add_argument('-output',metavar='output_file',dest='output_file',type=str,default='zqp_init.dat',help='set output file name')
    return parser.parse_args()

## output header
def print_header(output_file,L,Nsub):
    f = open(output_file,'w')
    f.write('0.0 0.0 0.0\n') # ene
    f.write('0.0 0.0 0.0\n') # var
    f.write('0.0 0.0 0.0\n') # Pg
    f.write('0.0 0.0 0.0\n') # Pj
    f.close()

## output orbital
def print_orb(output_file,L,Nsub):
    f = open(output_file,'a')
    Ns = L
    Norb = Nsub*(Ns/2+1)
    Nup = int(Ns/2)
    m0 = int(int(Nup)/2)
    fij = np.ones(Ns,dtype=float) # fij = 1
    for x in range(Norb):
        for m in range(1,m0+1):
            k = 2.0*np.pi*m/Ns
            fij[x] += 2.0*np.cos(k*x) # fij += 2cos(kx)
    for x in range(Norb):
        fij[x] *= 4.0/Ns
#    fij[0] = 0.0
    for x in range(Norb):
        f.write('%.16e 0.0 0.0\n' % fij[x])
    f.write('0.0 0.0 0.0\n') # orbitalidxpara up
    f.write('0.0 0.0 0.0\n') # orbitalidxpara down
    f.close()

## main
def main():
    args = parse_args()
    L = args.L
    Nsub = args.Nsub
    output_file = args.output_file

    print_header(output_file,L,Nsub)
    print_orb(output_file,L,Nsub)

if __name__ == "__main__":
    main()
