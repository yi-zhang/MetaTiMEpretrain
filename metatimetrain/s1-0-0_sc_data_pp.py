#!/usr/bin/env python
# coding: utf-8

# Pre-processing all TISCH data in scanpy.


import sys
import subprocess
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#%matplotlib inline
import seaborn as sns
import scipy.stats as ss
import importlib as imp
mypyPath = '/homes6/yiz/yipy/'
if(mypyPath not in sys.path):
    sys.path.append(mypyPath)
    
mypyPath = '../scICA/'
if(mypyPath not in sys.path):
    sys.path.append(mypyPath)

import gc
import scanpy as sc
import tischfun
from multiprocessing import Pool

import importlib as imp
import config
imp.reload(config)


# global
n_process = 12 


##### UNIT TEST #####
job_list = config.cohort_inIC[:2] 
#job_list = [ 'BCC_GSE123813_aPD1' ]
##### RUN ALL #####
#job_list = config.cohort_inIC 
job_list = config.allnames 
job_list = ['LSCC_GSE150321']


def worker( cohortname ): 
    filename = cohortname + '_res.pp.h5ad'
    file = os.path.join( config.SCPPDIR, filename )
    if(os.path.isfile(file)):
        print(f'[Log]File exists: {config.SCPPDIR}/{cohortname}')
        return(True)
    adata = tischfun.readTischAnn( cohortname , preprocessing=False)
    #if('norm_data' not in adata.layers or 'umap_cell_embeddings' not in adata.obsm.keys() or 'X_umap' not in adata.obsm.keys() ):
    print(f'[Log]pp: {cohortname}')
    adata = tischfun.readTischAnn( cohortname , preprocessing=True)
    adata.write( file )
    print(f'[Log]pp done: {config.SCPPDIR}/{cohortname}')
        
    # save adatap.obs
    import gc
    gc.collect()


with Pool(processes = min( len(job_list), n_process) ) as pool:
    pool.map( worker, job_list)







