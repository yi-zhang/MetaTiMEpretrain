#!/usr/bin/env python
# coding: utf-8

# Pre-processing all dataset from 10x triplet format to h5ad. Also accept: 10x triplets, h5file, in folder or in subfolder (1 sample). If in subfolder, the result file will be named by the subfolder name. Read into h5ad, compatible with .gz compressed files


import sys
import subprocess
import os
import argparse
import numpy as np
import pandas as pd
import importlib as imp
import gc
import scanpy as sc
from multiprocessing import Pool
import gzip
import shutil
from src import tischfun #import tischfun
import config

def decompressgz(datadir, outdir, filename, mode ='nodecompress'):
    file = os.path.join(datadir, filename)
    if (mode =='decompress' and filename[-3:] =='.gz'):
        
        newfilename = filename[:-3]
        newfile = os.path.join(outdir, newfilename)
        with gzip.open(file, 'rb') as f_in:
            with open(newfile, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    else:
        
        shutil.copy(file, outdir)


        
def decompress10x(datadir, outdir='../TMP/'):
    
    # triplet files
    prefix10s = [ t.replace('matrix.mtx.gz','') for t in os.listdir(datadir) if ('matrix.mtx.gz' in t and os.path.isfile(t)) ] 
    for prefix10 in prefix10s:
        # two types iof feature suffix
        fea_fn = prefix10 + 'features.tsv.gz'
        if(not os.path.isfile(os.path.join(datadir, fea_fn))):
            fea_fn = prefix10 + 'genes.tsv.gz'
            if(not os.path.isfile(os.path.join(datadir, fea_fn))):
                print('[Warning] Skip ', prefix10)
                continue
        print('[Log] Collecting ', prefix10)
        mtx_fn = prefix10 + 'matrix.mtx.gz'
        bar_fn = prefix10 + 'barcodes.tsv.gz'
        #[todo] check compression need
        decompressgz(datadir, outdir, mtx_fn)
        decompressgz(datadir, outdir, fea_fn)
        decompressgz(datadir, outdir, bar_fn)
    return(prefix10s)

if __name__ == '__main__':
    """
    Example
    ----------
    ls /liulab/yiz/DATA/scRNAmm/202306/geo_breast| while read folder; do echo $folder; python scformat.py -d /liulab/yiz/DATA/scRNAmm/202306/geo_breast/$folder -o ../TMP/ ; done
    """
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d','--datadir', help="input dir of 10x triplet files")
    parser.add_argument('-l','--listfile', default='', help="list with filenames")
    parser.add_argument('-t','--threads',default=4,type=int,help='Number of threads')
    parser.add_argument('-o', '--outdir', default='../TMP/',type=str, help="output directory")

    args = parser.parse_args()
    if not os.path.isdir(args.outdir): #generate output folder if not there
        os.makedirs(args.outdir)
        
    # for 10x triplets gz in folder directly
    prefix10s = decompress10x(args.datadir, outdir = args.outdir) #[todo] a few format not comptible. see  GSM5645904_Unknown_Carc_
    ## getting 10x prefix
    for prefix10 in prefix10s:
        print('[Log] converting ', prefix10)
        outfile = os.path.join(args.outdir, prefix10+'.h5ad')
        if(os.path.isfile(outfile)):
            print('[Log File Exist]', outfile)
            continue
        adata = sc.read_10x_mtx(path = args.datadir, prefix = prefix10)
        adata.write(outfile)
       
    # if this is a args.datadir where prefix shall be subfolder names
    prefix10dirs = [t for t in os.listdir(args.datadir) if os.path.isdir(t)] 
    for prefix10dir in prefix10dirs:
        if 'matrix.mtx.gz' in os.listdir( os.path.join(args.datadir, prefix10dir)): #the mtx usually does not have a name
            print('[Log] converting 10x triplet in subfolder', prefix10dir)
            outfile = os.path.join(args.outdir, prefix10dir+'.h5ad')
            if(os.path.isfile(outfile)):
                print('[Log File Exist]', outfile)
                continue
            adata = sc.read_10x_mtx(path = os.path.join(args.datadir, prefix10dir) )
            adata.write(outfile)
            
    # h5 file in the folder directly
    h5files = [t for t in os.listdir(args.datadir) if ( os.path.isfile(t) and t.split('.')[-1]=='h5') ]
    prefixh5s = [t[:-3] for t in h5files]
    for prefixh5 in prefixh5s:
        print('[Log] converting h5 in folder', prefixh5)
        outfile = os.path.join(args.outdir, prefixh5+'.h5ad')
        if(os.path.isfile(outfile)):
            print('[Log File Exist]', outfile)
            continue
        adata = sc.read_10x_h5(os.path.join(args.datadir, prefixh5+'.h5'))
        adata.write(outfile)
        
    
                                                         
    