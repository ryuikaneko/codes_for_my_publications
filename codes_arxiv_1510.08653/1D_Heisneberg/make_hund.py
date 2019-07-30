#!/usr/bin/env python

# coding:utf-8
import numpy as np
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='prepare hund.def')
    parser.add_argument('-L',metavar='L',dest='L',type=int,default=16,help='set L')
    parser.add_argument('-J1',metavar='J1',dest='J1',type=int,default=1.0,help='set J1')
    parser.add_argument('-J2',metavar='J2',dest='J2',type=int,default=0.5,help='set J2')
    parser.add_argument('-interaction_name',metavar='interaction_name',dest='interaction_name',type=str,default='hund',help='set interaction_name')
    parser.add_argument('-output',metavar='output_file',dest='output_file',type=str,default='hund.def',help='set output file name')
    return parser.parse_args()

## output header
def print_header(output_file,L,interaction_name,num):
    f = open(output_file,'w')
    f.write('========================\n')
    f.write('N%s %d\n' % (interaction_name,num))
    f.write('========================\n')
    f.write('%s\n' % interaction_name)
    f.write('========================\n')
    f.close()

## output int
def print_int(output_file,L,J1,J2):
    f = open(output_file,'a')
    parJ1 = -0.5*J1
    parJ2 = -0.5*J2
    for i in range(L):
        f.write('%3d %3d %.16e\n' % (i,(i+1)%L,parJ1))
#        f.write('%3d %3d %.16e\n' % (i,(i+2)%L,parJ2))
    f.close()

## main
def main():
    args = parse_args()
    L = args.L
    J1 = args.J1
    J2 = args.J2
    interaction_name = args.interaction_name
    output_file = args.output_file

    num = L
#    num = 2*L

    print_header(output_file,L,interaction_name,num)
    print_int(output_file,L,J1,J2)

if __name__ == "__main__":
    main()
