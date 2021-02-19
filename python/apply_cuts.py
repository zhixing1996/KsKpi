#!/usr/bin/env python
"""
Apply cuts on root files
"""

__author__ = "Maoqiang JING <jingmq@ihep.ac.cn>"
__copyright__ = "Copyright (c) Maoqiang JING"
__created__ = "[2020-01-07 Tue 23:30]"

import math
from array import array
import ROOT
from ROOT import TCanvas, gStyle, TLorentzVector, TTree
from ROOT import TFile, TH1F, TLegend, TArrow, TChain, TVector3
import sys, os
import logging
from math import *
from tools import *
logging.basicConfig(level=logging.DEBUG, format=' %(asctime)s - %(levelname)s- %(message)s')

def usage():
    sys.stdout.write('''
NAME
    apply_cuts.py

SYNOPSIS
    ./apply_cuts.py [file_in] [file_out]

AUTHOR
    Maoqiang JING <jingmq@ihep.ac.cn>

DATE
    January 2020
\n''')

def apply(file_in, file_out):
    try:
        chain = TChain('vfit')
        chain.Add(file_in)
    except:
        logging.error(file_in + ' is invalid!')
        sys.exit()

    cut = ''
    cut_dl = '(vfit2_dl/vfit2_dle<2.3) && '
    cut_ks = '(fabs(vfit2_mks-0.497614)>0.008) && '
    cut_kst =  '(var_kstar>=0.12)'
    cut = cut_dl + cut_ks + cut_kst

    t = chain.CopyTree(cut)
    t.SaveAs(file_out)

def main():
    args = sys.argv[1:]
    if len(args)<2:
        return usage()
    file_in = args[0]
    file_out = args[1]

    print '--> Begin to process file: ' + file_in
    apply(file_in, file_out)
    print '--> End of processing file: ' + file_in

if __name__ == '__main__':
    main()
