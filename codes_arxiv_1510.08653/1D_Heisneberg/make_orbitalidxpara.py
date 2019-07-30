#!/usr/bin/env python

# coding:utf-8
import numpy as np
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='prepare orbitalidxpara.def')
    parser.add_argument('-L',metavar='L',dest='L',type=int,default=16,help='set L')
    parser.add_argument('-Nsub',metavar='Nsub',dest='Nsub',type=int,default=1,help='set Nsub')
    parser.add_argument('-output',metavar='output_file',dest='output_file',type=str,default='orbitalidxpara.def',help='set output file name')
    return parser.parse_args()

## output header
def print_header(output_file,L,Nsub):
    f = open(output_file,'w')
    Ns = L
    Norb = 1
    f.write('========================\n')
    f.write('NOrbitalIdxPara %4d\n' % Norb)
#    f.write('ComplexType %4d\n' % 0)
    f.write('ComplexType %4d\n' % 1)
    f.write('========================\n')
    f.write('========================\n')
    f.close()

## output int
def print_int(output_file,L,Nsub):
    f = open(output_file,'a')
    Ns = L
    Norb = Nsub*(Ns/2+1)
    for i in range(Ns):
        for j in range(i+1,Ns):
            f.write('%3d %3d %3d\n' % (i,j,0))
    f.write('%3d %3d\n' % (0,0))
    f.close()

## main
def main():
    args = parse_args()
    L = args.L
    Nsub = args.Nsub
    output_file = args.output_file

    print_header(output_file,L,Nsub)
    print_int(output_file,L,Nsub)

if __name__ == "__main__":
    main()
