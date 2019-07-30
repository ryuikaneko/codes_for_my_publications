#!/usr/bin/env python

# coding:utf-8
import numpy as np
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='prepare greenone.def')
    parser.add_argument('-L',metavar='L',dest='L',type=int,default=16,help='set L')
    parser.add_argument('-Nsub',metavar='Nsub',dest='Nsub',type=int,default=1,help='set Nsub')
    parser.add_argument('-output',metavar='output_file',dest='output_file',type=str,default='greenone.def',help='set output file name')
    return parser.parse_args()

## output header
def print_header(output_file,L,Nsub):
    f = open(output_file,'w')
    Ns = L
    num = Nsub*4
    f.write('========================\n')
    f.write('NCisAjs               %d\n' % num)
    f.write('========================\n')
    f.write('= green one ============\n')
    f.write('========================\n')
    f.close()

## output int
def print_int(output_file,L,Nsub):
    f = open(output_file,'a')
    for i in range(Nsub):
        j=i
        f.write('%3d %3d %3d %3d\n' % (i,0,j,0))
        f.write('%3d %3d %3d %3d\n' % (i,0,j,1))
        f.write('%3d %3d %3d %3d\n' % (i,1,j,0))
        f.write('%3d %3d %3d %3d\n' % (i,1,j,1))
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
