#!/usr/bin/env python
# coding: utf-8

# # meta-components calling

import sys
import argparse
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

from src import tischfun, icfun, metacomp, complexClustermap
import sklearn
import gc
from multiprocessing import Pool
#import config
from sknetwork.clustering import Louvain # need scikitnetwork >=0.28.2

plt.rcParams['figure.dpi'] = 200
## for test: min_cluster_size = 2
#OUTDIR = '../analysis/20211001_scICA/'



# ### load concatenated components
def callmec(compconcatDIR, OUTDIR, min_cluster_size = 2):
    ic100res = icfun.load_components_precomputed( compconcatDIR)

    comp_concat_z = ic100res['comp_concat_z']
    comp_concat_z_topg = ic100res['comp_concat_z_topg']
    ic_sims = ic100res['ic_sims']
    rcp_concat_z = ic100res['rcp_concat_z']
    rcp_concat_z_topg = ic100res['rcp_concat_z_topg']


    ### statis
    print( f"comp_concat_z_topg: {len(comp_concat_z_topg.columns.str.split('.').str.get(0).drop_duplicates())} cohorts, {len(comp_concat_z_topg.index)} genes" )
    print( f"rcp_concat_z: {len(rcp_concat_z.columns.str.split('.').str.get(0).drop_duplicates())} cohorts, {len(rcp_concat_z.index)} genes" )
    print( f"rcp_concat_z_topg: {len(rcp_concat_z_topg.columns.str.split('.').str.get(0).drop_duplicates())} cohorts, {len(rcp_concat_z_topg.index)} genes" )
    '''

    # ## Model selection

    # In[10]:


    params = np.arange(1,1.35,0.05)


    # In[11]:


    score = []
    for reso in params:
        comp_cluster_louvain_pd, n = metacomp.cluster_louvain( rcp_concat_z_topg.T , resolution = reso, min_cluster_size = 5 )
        s = metacomp.eval_comp_cluster( rcp_concat_z_topg.T.loc[comp_cluster_louvain_pd.index], 
                                    comp_cluster_louvain_pd, metric='cosine')
        score.append([n,s, reso])


    # In[23]:


    score = np.array(score)

    fig, ax = plt.subplots(2, figsize = (6,4), sharex=True)
    ax[1].plot( np.array(score)[:,2], np.array(score)[:,1],'o--' , label = 'silhouette_clustering_cosine_', color = 'orange')
    ax[1].set_ylabel('Silhouette score')

    ax[0].plot( np.array(score)[:,2], np.array(score)[:,0],'o--' , label = 'n_mod')
    ax[0].set_ylabel('Number of MeC')

    ax[1].set_xlabel(' Louvain resolution')
    plt.savefig('../figure/model_selection_louvain_1.25.pdf',dpi = 300)


    # ## set best model
    best_reso = params[ np.argmax(np.array(score)[:,0]) ] 


    '''

    comp_cluster_louvain_pd, n = metacomp.cluster_louvain( rcp_concat_z_topg.T, resolution = 1.25, min_cluster_size = min_cluster_size )
    print(comp_cluster_louvain_pd) #[testprint]
    comp_cluster_louvain_pd = metacomp.trim_comp_cluster( comp_cluster_louvain_pd, min_cluster_size = min_cluster_size )
    print(comp_cluster_louvain_pd) #[testprint]


    # ### Plotting the MeC heatmap


    # value: cosine similarity
    sim = sklearn.metrics.pairwise.cosine_similarity(rcp_concat_z_topg.T)
    sim = pd.DataFrame(sim, index  = rcp_concat_z_topg.columns, columns = rcp_concat_z_topg.columns )


    # In[29]:


    # color clustesr and dataset and cancerType
    comp_order = comp_cluster_louvain_pd.sort_values(by='cluster').index
    comp_cluster_louvain_pd = comp_cluster_louvain_pd.loc[comp_order]

    comp_cluster_louvain_pd, cluster_color_dict = complexClustermap.colorCluster( comp_cluster_louvain_pd )
    comp_cluster_louvain_pd['dataset'] = comp_cluster_louvain_pd.index.str.split('.').str.get(0)
    comp_cluster_louvain_pd, dataset_color_dict = complexClustermap.colorCluster( comp_cluster_louvain_pd , colorcol = 'dataset')
    comp_cluster_louvain_pd['type'] = comp_cluster_louvain_pd['dataset'].str.split('_').str.get(0)
    comp_cluster_louvain_pd, cancertype_color_dict = complexClustermap.colorCluster( comp_cluster_louvain_pd , colorcol = 'type')

#    print(comp_cluster_louvain_pd[:2]) #[testprint]
#    print(comp_cluster_louvain_pd.shape)
    # In[30]:


    # reorder
    dat = sim.loc[comp_order, comp_order[::-1]]
    cfig = sns.clustermap(dat, figsize = (18,18), row_cluster=False, col_cluster=False,
                    dendrogram_ratio = 0.2, cmap = 'RdBu_r',
                        vmin = -1, vmax = 1, center = 0,
                        xticklabels = False,yticklabels=False,
                        row_colors = \
                            comp_cluster_louvain_pd[[ 'cluster_color','dataset_color','type_color']].rename(columns={'cluster_color': 'MeC','dataset_color': 'DataSet', 'type_color':'Type'}),
                        )
    plt.savefig(os.path.join(OUTDIR,'MeC_louvain_cluster_1.25.pdf'),dpi = 300)


    # In[869]:

    # ## save best model

    # In[31]:


    comp_cluster_label_file = os.path.join(OUTDIR, 'comp_cluster_labels.txt')
    comp_cluster_louvain_pd.to_csv( comp_cluster_label_file, sep='\t')

    '''
    # ## Check Top genes

    # In[32]:


    for i in range(len(comp_cluster_louvain_pd.drop_duplicates())):
        
        #if(i==0):
        #    continue
        icfun.gene_clusteri( cluster_i= i, 
                        comp_cluster_labels_pd=comp_cluster_louvain_pd,
                        component_concat=rcp_concat_z_topg, 
                topn=30, polar='positive',
                #savefile = '../analysis/ic100_corr0.3_cluster80/c'+str(i)+'.sig.txt',
                # saverowf = open('../analysis/ic100_corr0.3_cluster80/sig.20.txt','a')
                )
        
    '''

    # ## Save MeC genes for interpretmec


    # In[37]:


    metacomp_rank, metacomp_weightaverage = metacomp.mec_gene_full_table( comp_cluster_louvain_pd , rcp_concat_z_topg )
    mecmodule = metacomp.mec_gene_top( metacomp_weightaverage  )
    mec_topg_weightaverage = metacomp.mec_score_topg( metacomp_weightaverage , mecmodule, meclst_gene_col = 'TopGeneSig5' ) 
    metacomp_rank.to_csv( os.path.join( OUTDIR, 'MeC_allgene_ranksum-ranks.txt') ,sep='\t')
    metacomp_weightaverage.to_csv( os.path.join( OUTDIR, 'MeC_allgene_average-weights.txt') ,sep='\t')
    mec_topg_weightaverage.to_csv( os.path.join( OUTDIR, 'MeC_TopGeneSig5_average-weights.txt' ), sep='\t')
    mecmodule.to_csv( os.path.join( OUTDIR, 'MeC_topgene.txt') ,sep='\t')

    metacomp.mec_gene_lstfiles( mecmodule, OUTDIR = OUTDIR )

# loading
"""
mecres = metacomp.load_mec_precomputed()
mec_topg = mecres['mec_topg']
mec_rank = mecres['mec_rank']
mec_score = mecres['mec_score']
mec_cluster = mecres['mec_cluster']
"""



__doc__ = """"""
if __name__ == '__main__':
    # python callmec.py -c ../test/analysis/decompose_pull/ -o ../test/analysis/MeC/ -u True -s 3
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-c','--compconcatdir', help="input dir with concatenated components")
    parser.add_argument('-t','--threads',default=1,type=int,help='Number of threads')
    parser.add_argument('-o', '--outmecdir', default='../',type=str, help="output directory")
    parser.add_argument('-s', '--minmecsize', default=5,type=int, help="minimum mec cluster size.")
    parser.add_argument('-u', '--unittest', default=False,type=bool, help="unit test mode or not. affect minmecsize.")

    args = parser.parse_args()

    if not os.path.isdir(args.outmecdir): #generate output folder if not there
        os.makedirs(args.outmecdir)
    
    #OUTDIR = config.SCMECDIR#
    if args.unittest == True:
        args.minmecsize = 3
    callmec( args.compconcatdir , args.outmecdir, args.minmecsize)
    
