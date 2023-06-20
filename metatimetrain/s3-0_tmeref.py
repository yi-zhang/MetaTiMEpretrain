#!/usr/bin/env python
# coding: utf-8

# ### calculate per-cell projection for all tisch data.

# In[12]:


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
import tischfun
from multiprocessing import Pool
import importlib as imp
import config; imp.reload(config)
import anndata
import metacomp; imp.reload(metacomp)
import mecmapper;imp.reload(mecmapper)

# In[35]:


# global
n_process = 4 
OUTDIR = config.SCMECTISCHDIR 


mecnamedict = metacomp.load_mecname()
mecres = metacomp.load_mec_precomputed()
mec_topg = mecres['mec_topg']
mec_rank = mecres['mec_rank']
mec_score = mecres['mec_score']
mec_score_topg = mecres['mec_score_topg']
mec_cluster = mecres['mec_cluster']

# In[36]:



if not os.path.exists( OUTDIR ):
    os.makedirs( OUTDIR )
#name = 'BCC_GSE123813_aPD1'
# job dictionary
#job_list = [ 'BCC_GSE123813_aPD1' ]
##### UNIT TEST #####
job_list = config.cohort_inIC[:2] 
##### RUN ALL #####
#job_list = config.cohort_inIC 
job_list = config.allnames 

print(f'[Log] Processing {len(job_list)} datasets')

def worker( cohortname ): 
    print('[Log] calculating ', cohortname)
    filename = cohortname + '.mec.scoretopg.proj.noscale.txt'
    zfilename = cohortname + '.mec.scoretopg.proj.zscale.txt'
    obsfilename = cohortname + '.mec.obs.txt'
    file = os.path.join( OUTDIR, filename )
    zfile = os.path.join( OUTDIR, zfilename )
    obsfile = os.path.join( OUTDIR, obsfilename )
    if(os.path.isfile(file)):
        print(f'[Log]File exists: {cohortname}')
        return(True)
    adata = tischfun.readTischAnn( cohortname , preprocessing=False, DATADIR = config.SCPPDIR, suffix = '_res.pp.h5ad')
    if('norm_data' not in adata.layers):
        print(f'[Log]Calculating norm_data: {cohortname}')
        adata = tischfun.readTischAnn( cohortname , preprocessing=True)
    #adatap = mecmapper.projectMecAnn( adata, mecscore_topg, scaling = False, )
    pdataz = mecmapper.projectMecAnn(adata, mec_score_topg, sigscaling=True, scaling=False, addon=False)
    pdata = mecmapper.projectMecAnn(adata, mec_score_topg, sigscaling=False, scaling=False, addon=False)
    #adatap = ic100.projectModuleAnn( adata, mod1k )
    print(f'[Log]Projected: {cohortname}')
        
    # save adatap.obs
    pdataz.to_df().to_csv( zfile, sep='\t')
    pdata.to_df().to_csv( file, sep='\t')
    pdataz.obs.to_csv( obsfile , sep='\t')
    # memory release
    import gc
    gc.collect()
    return(True)


with Pool(processes = min( len(job_list), n_process) ) as pool:
    pool.map( worker, job_list)

