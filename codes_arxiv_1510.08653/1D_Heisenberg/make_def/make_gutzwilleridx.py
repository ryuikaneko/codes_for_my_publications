#!/usr/bin/env python

# coding:utf-8
import numpy as np
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='prepare gutzwilleridx.def')
    parser.add_argument('-L',metavar='L',dest='L',type=int,default=16,help='set L')
    parser.add_argument('-output',metavar='output_file',dest='output_file',type=str,default='gutzwilleridx.def',help='set output file name')
    return parser.parse_args()

## output header
def print_header(output_file):
    f = open(output_file,'w')
    f.write('========================\n')
    f.write('NGutzwillerIdx %d\n' % 1)
    f.write('ComplexType %d\n' % 0)
    f.write('========================\n')
    f.write('========================\n')
    f.close()

## output int
def print_int(output_file,L):
    f = open(output_file,'a')
    for i in range(L):
        f.write('%3d %3d\n' % (i,0))
    f.write('%3d %3d\n' % (0,0))
    f.close()

## main
def main():
    args = parse_args()
    L = args.L
    output_file = args.output_file

    print_header(output_file)
    print_int(output_file,L)

if __name__ == "__main__":
    main()
