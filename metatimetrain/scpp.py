#!/usr/bin/env python
# coding: utf-8

# Pre-processing all dataset provided in list.


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
    
from src import tischfun #import tischfun


def worker( filename, indir, outdir ): 
    """
    Parallel processing single datasets, read in {indir}/{filename} and output {outdir}/{filename}.pp.h5ad. autodetect count layer and integer in adata.X. simple filtering.

    """
    #filename, indir, outdir = args
    file = os.path.join(indir, filename)
    #file = os.path.join( config.SCPPDIR, filename )
    outfilename = filename+'.pp.h5ad'
    outfile = os.path.join(outdir, outfilename)

    if(os.path.isfile(outfile)):
        print(f'[Log]File exists; skipping preprocessing {outfile}')
        return(True)
    print(f'[Log]preprocessing: {file}')
    adata = tischfun.readTischAnn( file , suffix='',  preprocessing=True)
    adata.write( outfile )
    print(f'[Log]pp done: {outfile}')

    import gc
    gc.collect()

'''def ppDatasets( job_list_args, n_process = 4):
    """ 
    paralell pre-processing scrna dataset
    ##### UNIT TEST #####
    ## parse
    job_list = config.testnames 
    input_dir = '/liulab/yiz/proj/MetaTiMEpretrain/test/data/icb/'
    #job_list = [ 'BCC_GSE123813_aPD1' ]
    ##### RUN ALL #####
    #job_list = config.cohort_inIC 

    """'''
    
    

__doc__ = """
Example
----------
# Unit test
python scpp.py --datadir ../test/data/icb/ --listfile ../test/data/datalst1.tsv --outdir ../test/analysis/pp
# Example, mmmetatime.
python scpp.py --datadir ../mmmetatime202307cohort/gsedata/ --outdir ../mmmetatime202307cohort/pp/
"""
if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d','--datadir', help="input dir of h5ad  files")
    parser.add_argument('-l','--listfile', default='', help="list with filenames")
    parser.add_argument('-t','--threads',default=4,type=int,help='Number of threads')
    parser.add_argument('-o', '--outdir', default='../TMP/',type=str, help="output directory")

    args = parser.parse_args()

    if not os.path.isdir(args.outdir): #generate output folder if not there
        os.makedirs(args.outdir)
    if not os.path.isfile(args.listfile): 
        filenames = [ t for t in os.listdir(args.datadir) if t.split('.')[-1] == 'h5ad']
    else:
        filenames = pd.read_table(args.listfile, names=['file'])['file'].values.tolist()
        # []Error
        #print('Error, no input filelist found')
        
    
    job_list_args = [
        (filename, args.datadir, args.outdir) for filename in filenames
    ]
    
    with Pool(processes = min( len(job_list_args), args.threads) ) as pool:
        pool.starmap( worker, job_list_args)



