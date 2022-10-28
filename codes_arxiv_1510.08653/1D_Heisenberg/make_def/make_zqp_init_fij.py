#!/usr/bin/env python

# coding:utf-8
import numpy as np
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='prepare zqp_init.dat')
    parser.add_argument('-L',metavar='L',dest='L',type=int,default=16,help='set L')
    parser.add_argument('-p',metavar='p',dest='p',type=float,default=1.0,help='set p')
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
def print_orb(output_file,L,p,Nsub):
    f = open(output_file,'a')
    Ns = L
    Norb = Nsub*(Ns/2+1)
    kx = np.array([2.0*np.pi*x/L for x in range(L)])
    enek = np.array([-2.0*np.cos(x) for x in kx])
    delk = np.array([+pow(np.abs(x),p) if x<=0 else -pow(np.abs(x),p) for x in enek])
    denom = enek + np.sqrt(enek*enek + delk*delk)
    fk = delk/denom
    fij = np.real(np.fft.ifft(fk))
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
    p = args.p
    Nsub = args.Nsub
    output_file = args.output_file

    print_header(output_file,L,Nsub)
    print_orb(output_file,L,p,Nsub)

if __name__ == "__main__":
    main()
