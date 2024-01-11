#!/usr/bin/env python
# coding: utf-8

import sys
import subprocess
import os
import argparse
import numpy as np
import pandas as pd
import scipy.stats as ss
import importlib as imp
import scanpy as sc
from multiprocessing import Pool    
import gc
from src import  icfun

def concatcomp(compdir, outdir):
    """pull together all components"""

    print('[Log] Concatenating all components')
    comp_concat = icfun.concat_ic( compDIR = compdir , OUTDIR = outdir )
    print('[Log] ztransform all components')
    comp_concat_z = icfun.qc_comp_ztrans(comp_concat)
    print('[Log] Getting top genes on all components')
    comp_concat_z_topg = icfun.qc_comp_topg(comp_concat_z )
    print('[Log] Eval similarity of all components')
    ic_sims = icfun.comps_similarity(comp_concat_z_topg, file = 'sim_pairwise_components.txt', OUTDIR = outdir)
    print('[Log] Only keep components with correlation >=0.3 with at least one compoennt from another cohort.  ')
    rcp_z = icfun.qc_comp_maxcor(comp_concat_z_topg, ic_sims, CORRCUT = 0.3)
    print('[Log] Getting top genes on selected components')
    rcp_topg = icfun.qc_comp_topg( rcp_z)
    
    comp_concat_z.to_csv( os.path.join( outdir, 'components_all_concat.txt' ), sep='\t')
    comp_concat_z_topg.to_csv( os.path.join( outdir, 'components_topg_concat.txt' ), sep='\t')
    ic_sims.to_csv(os.path.join( outdir, 'components_similarity.txt' ), sep='\t')
    rcp_z.to_csv( os.path.join( outdir, 'rcp_z.txt' ), sep='\t')
    rcp_topg.to_csv( os.path.join( outdir, 'rcp_topg.txt' ), sep='\t')
    
    print(f'[Log] output in {outdir}')
    print('comp_concat.shape', comp_concat.shape)
    print('comp_concat_z.shape', comp_concat_z.shape)
    print('comp_concat_z_topg.shape', comp_concat_z_topg.shape)
    print('rcp_z.shape', rcp_z.shape)
    print('rcp_topg.shape', rcp_topg.shape)
    

__doc__ = """"""
if __name__ == '__main__':
    # python pullcomp.py -c ../test/analysis/decompose_align/ -o ../test/analysis/decompose_pull/ 
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-c','--compdir', help="input dir of component weight files")
    parser.add_argument('-l','--listfile', default='', help="list with filenames")
    parser.add_argument('-t','--threads',default=4,type=int,help='Number of threads')
    parser.add_argument('-o', '--outdir', default='../',type=str, help="output directory")
    parser.add_argument('-s','--filesuffix', default='txt')
    args = parser.parse_args()

    if not os.path.isdir(args.outdir): #generate output folder if not there
        os.makedirs(args.outdir)
    if not os.path.isfile(args.listfile): 
        filenames = [ t for t in os.listdir(args.compdir) if t.split('.')[-len(args.filesuffix):] == args.filesuffix ]
    else:
        filenames = pd.read_table(args.listfile, names=['file'])['file'].values.tolist()
        print('Inputted list file: ', filenames)

        # if only performing on the subset, create a temp folder

        import shutil
        import os
        subsetcompdir = args.compdir + '/' + 'subset_to_call_mec/'
        if not os.path.isdir(subsetcompdir):
            os.makedirs(subsetcompdir)
        
        for filename in filenames:
            print(filename)
            realfilename = [t for t in os.listdir(args.compdir) if filename in t][0]
            shutil.copyfile(os.path.join(args.compdir, realfilename),
                            os.path.join(subsetcompdir, realfilename)
                            )
        args.compdir = subsetcompdir
        # []Error
        #print('Error, no input filelist found')
    concatcomp(args.compdir, args.outdir)   