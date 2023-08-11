## DNU
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
import config



def worker( gseid, datadir, outdir, prefixs): 
    """
    Parallel processing single datasets, read in {indir}/{filename} and output {outdir}/{filename}.pp.h5ad. autodetect count layer and integer in adata.X. simple filtering.

    """

    # output
    outfilename = gseid+'.h5ad'
    adatalst = []
    for prefix in prefixs:

        print('[Log] loading ', gseid, prefix)
        adata = sc.read(os.path.join(datadir, prefix +'.h5ad'))
        adatalst.append(adata)

    print('[Log] merging ', gseid)
    adata = adatalst[0].concatenate(adatalst[1:], join='inner' ,
            batch_categories = prefixs)

    del adatalst
    import gc
    gc.collect()

    adata.write(os.path.join(outdir, outfilename))






__doc__ = """
## A step to take per GSM 10x triplet samples and merge GSMs from the same GSE.

Example
----------
python mergecohort.py --datadir ../mmmetatime202306/data --metafile ../mmmetatime202307cohort/breast_final_downloaded_metatable_gsm.csv --outdir ../mmmetatime202307cohort/gsedata -t 8

"""
if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d','--datadir', help="input dir of h5ad  files")
    parser.add_argument('-m','--metafile', default='', help="metatable file with GSM and corresponding GSE.")
#    parser.add_argument('--samplecol',default='id',type=str,help='column name of GSM')
    parser.add_argument('--cohortcol',default='GSE',type=str,help='column name of GSE')
    parser.add_argument('-t','--threads',default=4,type=int,help='Number of threads')
    parser.add_argument('-o', '--outdir', default='../TMP/',type=str, help="output directory")
#    parser.add_argument('-l','--listfile', default='', help="list with filenames")

    args = parser.parse_args()

    if not os.path.isdir(args.outdir): #generate output folder if not there
        os.makedirs(args.outdir)


    #datadir = '/liulab/yiz/DATA/scRNAmm/202306/geo_breast'
#    datadir = '/liulab/yiz/proj/MetaTiMEpretrain/mmmetatime202306/data'
#    metafile = '/liulab/yiz/proj/MetaTiMEpretrain/mmmetatime202307cohort/breast_final_downloaded_metatable_gsm.csv'
    metatable = pd.read_table(args.metafile, sep=',')
    metatable['prefix10'] = metatable['mtx_fn'].str.replace('matrix.mtx.gz','') # matrix.mtx.gz

    prefixes = [t.replace('.h5ad','') for t in os.listdir(args.datadir) if '.h5ad' in t]
    print('metatable size', args.metafile, metatable.shape[0])
    metatable = metatable[metatable['prefix10'].isin(prefixes)]
    print('metatable size after  matching directory ',args.datadir, metatable.shape[0])

    ### group gses for the filenames.
    
    print('[Log]sample count')
    print(metatable.groupby(args.cohortcol).count()['id'])    

    allgseids = metatable[args.cohortcol].drop_duplicates().values.tolist()
    print('[Log] All GSEs')
    print(allgseids)

    gsegsmdict = {}
    for gseid in allgseids:
        prefix = metatable[metatable[args.cohortcol]==gseid]['prefix10'].values.tolist()
        gsegsmdict.update( {gseid: prefix})
    
    job_list_args = [
        (gseid, args.datadir, args.outdir, gsegsmdict[gseid]) for gseid in allgseids
    ]
    
    with Pool(processes = min( len(job_list_args), args.threads) ) as pool:
        pool.starmap( worker, job_list_args)

