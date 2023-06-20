import sys
import subprocess
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#%matplotlib inline
import seaborn as sns
import scipy.stats as ss
mypyPath = '/homes6/yiz/yipy/'
if(mypyPath not in sys.path):
    sys.path.append(mypyPath)
    
mypyPath = '../scICA/'
if(mypyPath not in sys.path):
    sys.path.append(mypyPath)

import joblib
import sklearn
import sklearn.decomposition as decomp
import scanpy as sc
import config
import icfun
import tischfun



name = sys.argv[1]
OUTDIR = sys.argv[2]
#name = 'AML_GSE116256'
#OUTDIR = '../TMP/'
n_components =  100
OUTfile = os.path.join(OUTDIR, name+'.sc.c'+ str(n_components)+'.txt')

if os.path.isfile( OUTfile ) :
    print(f'Exist: {OUTfile}')
    sys.exit(0)

adata = tischfun.readTischAnn(name , preprocessing = False, 
                           DATADIR = config.SCPPDIR,
                            suffix = '_res.pp.h5ad')

adata = adata[~adata.obs['assign_level3_anno'].str.contains('alignant')] # simple non-malignant
x = adata.to_df(layer = 'norm_data')

TEST = False 
if(TEST):

    X = x[:100].copy()
    n_components = 20; model, c = icfun.runIC(X, n_components = n_components )
    c.to_csv(os.path.join(OUTDIR, name+'.sc.c'+ str(n_components)+'.txt'), sep = '\t')
    joblib.dump( model, os.path.join( OUTDIR, name+'.sc.c'+ str(n_components)+'.joblib'))
    sys.exit()
else:
    X = x.copy()
    n_components = 100; model, c = icfun.runIC(X, n_components = n_components )
    c.to_csv( OUTfile , sep = '\t')
    #joblib.dump( model, os.path.join(OUTDIR, name+'.sc.c'+ str(n_components)+'.joblib')) 


