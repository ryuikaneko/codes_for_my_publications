#!/usr/bin/env python

# coding:utf-8
import numpy as np
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='prepare modpara.def')
    parser.add_argument('-L',metavar='L',dest='L',type=int,default=16,help='set L')
    parser.add_argument('-output',metavar='output_file',dest='output_file',type=str,default='modpara.def',help='set output file name')
    return parser.parse_args()

## output name of files
def print_files(output_file,L):
    f = open(output_file,'w')
    f.write('--------------------\n')
    f.write('Model_Parameters   0\n')
    f.write('--------------------\n')
    f.write('VMC_Cal_Parameters\n')
    f.write('--------------------\n')
    f.write('CDataFileHead  zvo\n')
    f.write('CParaFileHead  zqp\n')
    f.write('--------------------\n')
    f.write('NVMCCalMode    0\n')
    f.write('--------------------\n')
    f.write('NDataIdxStart  1\n')
    f.write('NDataQtySmp    1\n')
    f.write('--------------------\n')
    f.write('Nsite          %d\n' % L)
    f.write('Ncond          0\n')
#    f.write('2Sz            0\n')
#    f.write('NSPGaussLeg    1\n')
#    f.write('NSPStot        0\n')
    f.write('NMPTrans       1\n')
    f.write('NSROptItrStep  2000\n')
    f.write('NSROptItrSmp   200\n')
    f.write('DSROptRedCut   0.0010000000\n')
    f.write('DSROptStaDel   0.0200000000\n')
    f.write('DSROptStepDt   0.0200000000\n')
    f.write('NVMCWarmUp     10\n')
    f.write('NVMCInterval   1\n')
    f.write('NVMCSample     2000\n')
    f.write('NExUpdatePath  2\n')
    f.write('RndSeed        123456789\n')
    f.write('NSplitSize     1\n')
    f.write('NStore         1\n')
    f.write('NSRCG          0\n')
    f.close()

## main
def main():
    args = parse_args()
    L = args.L
    output_file = args.output_file

    print_files(output_file,L)

if __name__ == "__main__":
    main()
