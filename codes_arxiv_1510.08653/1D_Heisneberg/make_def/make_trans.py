#!/usr/bin/env python

# coding:utf-8
import numpy as np
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='prepare trans.def')
    parser.add_argument('-L',metavar='L',dest='L',type=int,default=16,help='set L')
    parser.add_argument('-output',metavar='output_file',dest='output_file',type=str,default='trans.def',help='set output file name')
    return parser.parse_args()

## output header
def print_header(output_file):
    f = open(output_file,'w')
    f.write('========================\n')
    f.write('NTransfer %d\n' % 0)
    f.write('========================\n')
    f.write('i_j_s_tijs\n')
    f.write('========================\n')
    f.close()

## main
def main():
    args = parse_args()
    L = args.L
    output_file = args.output_file

    print_header(output_file)

if __name__ == "__main__":
    main()
