#!/usr/bin/env python
# coding: utf-8

import sys
import subprocess
import os
import numpy as np
import pandas as pd
import argparse
import joblib
import sklearn
import sklearn.decomposition as decomp
import scanpy as sc
from multiprocessing import Pool
from src import tischfun, icfun


def decomposeData(filename, indir, outdir, n_components=100):
    """
    processing decomposition of single dataset.
    """
    #filename = 'AML_GSE116256_res.pp.h5ad'
    file = os.path.join(indir, filename)
    outfilename =  filename+'.sc.c'+ str(n_components)+'.txt'
    outfile = os.path.join(outdir, outfilename)
    adata = tischfun.readTischAnn(file , 
                            preprocessing = False, 
                            DATADIR = '',
                            suffix = '')
    if os.path.isfile( outfile ) :
        print(f'Exist: {outfile}')
        sys.exit(0)


    # [todo]simple malignant cell filtering
    #adata = adata[~adata.obs['assign_level3_anno'].str.contains('alignant')] # simple non-malignant

    # decomposition
    x = adata.to_df() #[todo] layer = 'norm_data'
    X = x.copy()
    print('[Log] Decomposing ', filename)
    model, c = icfun.runIC(X, n_components = n_components )
    c.to_csv( outfile , sep = '\t')
    #joblib.dump( model, os.path.join(OUTDIR, name+'.sc.c'+ str(n_components)+'.joblib')) 


def decomposeDatasets( job_list_args, n_process = 4):
    """DNU.
    Parallele processing decomposition of single dataset.
    job_list_args -> (filename, args.datadir, args.outdir)
    """
    n_parallel = min( len(job_list_args), n_process)
    print(n_parallel, ' nodes')
    with Pool(processes = n_parallel ) as pool:
        pool.starmap( decomposeData, job_list_args)


__doc__ = """"""

if __name__ == '__main__':
    """
    Example
    ----------
    python decompose.py -d ../test1/ -t 4 -o ../decompose/ -k 100
    """
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d','--datappdir', help="input dir of preprocessed scRNA h5ad files")
    parser.add_argument('-l','--listfile', default='', help="list with filenames")
    parser.add_argument('-t','--threads',default=4,type=int,help='Number of threads')
    parser.add_argument('-o', '--outdir', default='../decompose/',type=str, help="output directory of component files")
    parser.add_argument('-k','--kcomps', default=100, type=int, help='number of components')

    args = parser.parse_args()

    if not os.path.isdir(args.outdir): #generate output folder if not there
        os.makedirs(args.outdir)
    if not os.path.isfile(args.listfile): 
        filenames = [ t for t in os.listdir(args.datappdir) if t.split('.')[-1] == 'h5ad']
    else:
        # read in a list
        filenames = pd.read_table(args.listfile, names=['file'])['file'].values.tolist()
        # []Error
        #print('Error, no input filelist found')
        
    
    job_list_args = [
        (filename, args.datappdir, args.outdir) for filename in filenames
    ]

    n_parallel = min( len(job_list_args), args.threads)
    print(n_parallel, ' nodes')
    with Pool(processes = n_parallel ) as pool:
        pool.starmap( decomposeData, job_list_args)

    #decomposeDatasets( job_list_args, n_process = )

