# annotate metacomponents
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
import interpretmec 
from multiprocessing import Pool

import importlib as imp
import config
imp.reload(config)


# global
n_process = 12 
metacomp_weightaverage = pd.read_table( os.path.join( config.SCMECDIR, 'MeC_allgene_average-weights.txt'), sep='\t',index_col = 0)
##### UNIT TEST #####
#job_list = ['7']
##### RUN ALL #####
job_list = metacomp_weightaverage.columns.tolist( )

mecmodule = interpretmec.mec_gene_top( metacomp_weightaverage  )

mecanno = mecmodule.copy()
def worker( i ):
    print(i,'\n')
    glist = mecanno.loc[i, 'TopGene100' ].split(',')
    mec_enrich_file = os.path.join(config.SCMEC_enrich_DIR, 'enrich.mec.'+str(i)+'.pickle')
    enrspd = interpretmec.enrich( glist , printing=False, savefile = mec_enrich_file )
    enrspd.to_csv( os.path.join( config.SCMEC_enrich_DIR, 'enrich.mec.'+str(i)+'.txt') ,sep='\t')
    #mecanno.loc[i,  'Pathway'  ] = ';  '.join( enrspd.head(3)['Term'].tolist())
    return({ i : enrspd })

with Pool( processes = 12 ) as pool:
    pool.map( worker, job_list )

