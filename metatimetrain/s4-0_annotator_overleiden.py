## high reso leiden clustering . records top1 prediction. 


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
import scanpy as sc
# In[35]:


# global
n_process = 12
OUTDIR = config.SCANNOLEIDENDIR3


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
##### UNIT TEST #####
#job_list = config.cohort_inIC[:2] +['PBMC_8K']
#job_list = ['BRCA_Alex']
##### RUN ALL #####
#job_list = config.cohort_inIC
#job_list = config.allnames
#job_list = config.cohort_ICB
job_list = ['BCC_GSE123813_aPD1']
print(f'[Log] Processing {len(job_list)} datasets')

def annotator( projmat , gcol, mecnamedict):
    """ Annotator based on scores.
    ----
    projmat:
    dataframe with scores and one column gcol (grouping column) e.g. leiden class
    gcol:
    grouping column
    mecnamedict: for renaming projmat columns into functional names
    ----
    Return:
    projmat: matrix with predicted label column
    gpreddict: dictionary about how is each gcol mapped to a label
    
    """
    # median score for each group
    gpred = projmat.groupby( gcol ).median() # gpred: leiden cluster by mec.
    # top score for each group, with filtering
    gpreddict = projmat.groupby( gcol ).median().idxmax(axis=1).apply(lambda t: mecnamedict[t])# cluster max quantile 75% larger than zscore-1
    pass_maxmedian = (projmat.groupby( gcol ).quantile(q=0.75) >1 ).max(axis=1)
    gpreddict[ ~ pass_maxmedian] = 'Others'
    # generate top1 prediction
    projmat['assign_pred_'+gcol] = projmat[ gcol ].apply(lambda t: gpreddict[t])
    
    return( projmat, gpred, gpreddict )

def worker( cohortname ):
    """ leiden label prediction """
    print('[Log] calculating ', cohortname)
    h5file = os.path.join( OUTDIR, cohortname + '.proj.h5ad')
    leidenscorefile = os.path.join( OUTDIR, cohortname + '.mecscore.leiden.txt')

    if(os.path.isfile(h5file)):
        print(f'[Log]File exists: {cohortname}')
        return(True)

    # read data
    adata = tischfun.readTischAnn( cohortname , preprocessing=False, DATADIR = config.SCPPDIR, suffix = '_res.pp.h5ad')
    if('norm_data' not in adata.layers):
        print(f'[Log]Calculating norm_data: {cohortname}')
        adata = tischfun.readTischAnn( cohortname , preprocessing=True)

    # batch correction
    batchcols = list(set(config.batchcols).intersection(set(adata.obs.columns.tolist())))
    if(len(batchcols))>=1:
        print(f'[Log ] {cohortname} with batch ', batchcols)
        try:
            sc.external.pp.harmony_integrate(adata, key = batchcols )
            sc.pp.neighbors(adata, use_rep = 'X_pca_harmony')
        except:
            print('[Log ] {cohortname} harmony problem, skipping batch correction')
            sc.pp.neighbors(adata)
    else:
        print(f'[Log ] {cohortname} no batch label detected')
        sc.pp.neighbors(adata)

    # clustering and umap
    print(f'[Log ] {cohortname} Clustering ' )
    #sc.tl.leiden(adata, resolution=5)
    sc.tl.leiden(adata, resolution=8)
    print(f'[Log ] {cohortname} umap ' )
    sc.tl.umap(adata)

    print(f'[Log] {cohortname} Projecting')

    # projection
    pdata = mecmapper.projectMecAnn(adata, mec_score_topg, sigscaling=True, scaling=False, addon=False)

    # cell state prediction
    
    # per group prediction, maximum mean score, renaming
    gcol = 'leiden'
    projmat = adata.obs.join( pdata.to_df() )
    projmat = projmat[config.scorecols + [ gcol] ]
    
    projmat, gpred, gpreddict = annotator( projmat , gcol, mecnamedict)
    
    # write in pdata
    pdata.obsm['X_umap'] = adata.obsm['X_umap']
    if('X_pca_harmony' in adata.obsm):
        pdata.obsm['X_pca_harmony'] = adata.obsm['X_pca_harmony']
    pdata.obs['leiden'] = adata.obs['leiden']
    pdata.obs['assign_pred_'+gcol] = projmat['assign_pred_'+gcol]

    # save gpred and pdata.
    print(f'[Log] {cohortname} write to {h5file}')
    gpred.to_csv(leidenscorefile, sep='\t')
    pdata.write(h5file )

    # memory release
    del(adata)
    import gc
    gc.collect()
    return(True)


with Pool(processes = min( len(job_list), n_process) ) as pool:
    pool.map( worker, job_list)
    


