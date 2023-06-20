import sys
import subprocess
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# %matplotlib inline
import seaborn as sns
import scipy.stats as ss
mypyPath = '/homes6/yiz/yipy/'
if(mypyPath not in sys.path):
    sys.path.append(mypyPath)
    
mypyPath = '../scICA/'
if(mypyPath not in sys.path):
    sys.path.append(mypyPath)
    
import icfun
import config

n_components = 100 
#scicDIR = '../analysis/20210109_scICA/'
#OUTDIR = '../analysis/20210109_scICA/'
scicDIR = config.SCICDIR
OUTDIR = config.SCICDIR

files = [t for t in os.listdir(scicDIR) if 'sc.c'+str(n_components)+'.txt' in t]
names = [file.split('.')[0] for file in files]

components_dict = {}
for file in files:
    name = file.split('.')[0]
    components = pd.read_table(os.path.join(scicDIR , file ), index_col = 0)
    components_dict.update( { name: components })
    
for name in names:
    components = components_dict[ name ]

    components, skews = icfun.ic_rotate(  components, n_components )
    #skewfile = os.path.join( OUTDIR, name+'.sc.c'+str(n_components)+'.icskew.txt')
    rotate_components_file = os.path.join( OUTDIR, name +'.sc.rotate.c'+str(n_components)+'.txt' )
    #icskews = pd.DataFrame({ 'ic_skew': skews})
    #icskews.to_csv( skewfile, sep='\t')
    components.to_csv( rotate_components_file, sep='\t' )

