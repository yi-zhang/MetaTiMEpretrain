#!/usr/bin/env python
# coding: utf-8


import os
import pandas as pd

THIS_DIR = os.path.dirname( os.path.abspath(__file__))

def getalltfs( mode = 'all'):
    alltf = pd.read_csv(
        os.path.join(THIS_DIR, 'AllTF_Human.txt'),
        sep='\t')
    alltfs = alltf['Symbol'].drop_duplicates().values.tolist()
    alltfs = [t for t in alltfs if isinstance(t,str) ]# remove nan
    return(alltfs)
