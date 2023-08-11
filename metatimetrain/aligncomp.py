#!/usr/bin/env python
# coding: utf-8

import sys
import subprocess
import os
import numpy as np
import pandas as pd
import argparse
import scipy.stats as ss

from src import icfun
import config


#scicDIR = '../analysis/20210109_scICA/'
#OUTDIR = '../analysis/20210109_scICA/'
#indir = config.SCICDIR
#outdir = config.SCICDIR

def rotateComps( indir, outdir,filenames=[], n_components=100):
    """
    Read in all component files and output skew aligned ones
    """
    if( len(filenames)<1):
        filenames = [t for t in os.listdir(indir) if '.c'+str(n_components)+'.txt' in t]
        filenames = [t for t in filenames if 'rotate' not in filenames]
    
    names = [t.replace('.c'+str(n_components)+'.txt', '') for t in filenames]

    components_dict = {}
    for name, filename in zip( names, filenames):
        components = pd.read_table(os.path.join(indir , filename ), index_col = 0)
        components_dict.update( { name: components })
        
    for name, filename in zip( names, filenames):
        print('[Log] aligning components ', filename)
        components = components_dict[ name ]
        components, skews = icfun.ic_rotate(  components, n_components )
        #skewfile = os.path.join( OUTDIR, name+'.sc.c'+str(n_components)+'.icskew.txt')
        rotate_components_file = os.path.join( outdir, name +'.sc.rotate.c'+str(n_components)+'.txt' )
        #icskews = pd.DataFrame({ 'ic_skew': skews})
        #icskews.to_csv( skewfile, sep='\t')
        components.to_csv( rotate_components_file, sep='\t' )


__doc__ = """"""
if __name__ == '__main__':
    # python aligncomp.py -c ../test/analysis/decompose -o ../test/analysis/decompose_align/ -k 100
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-c','--compdir', help="input dir of component files")
    parser.add_argument('-l','--listfile', default='', help="list with filenames")
    parser.add_argument('-t','--threads',default=4,type=int,help='Number of threads')
    parser.add_argument('-o', '--outdir', default='../',type=str, help="output directory")
    parser.add_argument('-k','--kcomps', default=100, type=int, help='number of components, needs to be the same as decompose.py input')

    args = parser.parse_args()

    if not os.path.isdir(args.outdir): #generate output folder if not there
        os.makedirs(args.outdir)
    if not os.path.isfile(args.listfile): 
        filenames = [t for t in os.listdir(args.compdir) if '.c'+str(args.kcomps)+'.txt' in t]
        filenames = [t for t in filenames if 'rotate' not in filenames]
    else:
        filenames = pd.read_table(args.listfile, names=['file'])['file'].values.tolist()
     
    rotateComps( args.compdir, args.outdir, filenames=filenames, n_components=args.kcomps)