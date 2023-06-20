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
import scanpy as sc
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
import config
import icfun

comp_concat = icfun.concat_ic( compDIR = config.SCICDIR , OUTDIR = config.SCICDIR )
comp_concat = icfun.concat_ic( compDIR = config.SCICDIR , OUTDIR = config.SCICDIR )
comp_concat_z = icfun.qc_comp_ztrans(comp_concat)
comp_concat_z_topg = icfun.qc_comp_topg(comp_concat_z )
ic_sims = icfun.comps_similarity(comp_concat_z_topg)

comp_concat_z.to_csv( os.path.join(config.SCICDIR, 'components_all_concat.txt' ), sep='\t')
comp_concat_z_topg.to_csv( os.path.join(config.SCICDIR, 'components_topg_concat.txt' ), sep='\t')
ic_sims.to_csv(os.path.join(config.SCICDIR, 'components_similarity.txt' ), sep='\t')

print('comp_concat.shape', comp_concat.shape)
print('comp_concat_z.shape', comp_concat_z.shape)
print('comp_concat_z_topg.shape', comp_concat_z_topg.shape)

