#!/usr/bin/env python
# coding: utf-8

## function for meta components
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sklearn
import config
from src import getgene

#### Clustering of components into meta-components ####
# def comp_cluster_agg_pl( data, N_CLUSTERS = 72):
#     """ Hierachical clustering . euclidean, ward.
#     ---
#     Input: 
#     data : pd.DataFrame. pairwise simliarity between all components. s
#     ---
#     Return:
#     gene assignment to metacomponent, result need for re-plotting. 
#     ----
#     Usage
#     X = sklearn.metrics.pairwise.cosine_similarity(rcp_concat_z_topg.T)
#     x = pd.DataFrame(X, index = rcp_concat_z_topg.columns, columns = rcp_concat_z_topg.columns)
#     data = x.T
#     """
#     from scipy.spatial import distance
#     from sklearn.cluster import AgglomerativeClustering
#     from scipy.cluster.hierarchy import linkage
# 
#     gene_cluster_labels = AgglomerativeClustering( n_clusters = N_CLUSTERS, 
#                                                   #affinity = 'euclidean', 
#                                                   #linkage = 'ward',
#                                                   affinity = 'precomputed', # precomputed: input distance
#                                                   )\
#                         .fit(data).labels_
#     gene_cluster_colors = dict(zip(np.unique(gene_cluster_labels), sns.color_palette( n_colors = N_CLUSTERS)))
#     linkage = linkage(data, method = 'ward', optimal_ordering=True)
#     cohort_labels = data.index.str.split('.').str.get(0)
#     cohort_colors = dict(zip(np.unique(cohort_labels), sns.color_palette( n_colors = len( np.unique(cohort_labels)) )))
#     cfig = sns.clustermap(data, figsize = (18,18),
#                   dendrogram_ratio = 0.2, cmap = 'RdBu_r', cbar_pos = None,
#                    center = 0, 
#                    row_colors = [[gene_cluster_colors[cl] for cl in gene_cluster_labels],
#                                 [cohort_colors[cl] for cl in cohort_labels]
#                                 ],
#                    row_linkage = linkage,
#                   col_linkage = linkage)
#     gene_cluster_labels_pd = pd.DataFrame( gene_cluster_labels, 
#                  index = data.index[cfig.dendrogram_col.reordered_ind], 
#                 columns = ['cluster'])
#     #gene_cluster_labels_pd_icagg80 = gene_cluster_labels_pd
#     res = {' cfig' : cfig ,
#            'gene_cluster_colors': gene_cluster_colors,
#            'gene_cluster_labels': gene_cluster_labels,
#            'cohort_colors': cohort_colors,
#            'cohort_labels' : cohort_labels,
#            'linkage': 'linkage',
#     }
#     return( gene_cluster_labels_pd , res)

def comp_cluster_agg( data, N_CLUSTERS = 72):
    """ Hierachical clustering . euclidean, ward.
    ---
    Input: 
    data : pd.DataFrame. pairwise simliarity between all components. s
    ---
    Return:
    gene assignment to metacomponent, result need for re-plotting. 
    ----
    Usage
    X = sklearn.metrics.pairwise.cosine_similarity(rcp_concat_z_topg.T)
    x = pd.DataFrame(X, index = rcp_concat_z_topg.columns, columns = rcp_concat_z_topg.columns)
    data = x.T
    """
    from scipy.spatial import distance
    from sklearn.cluster import AgglomerativeClustering
    from scipy.cluster.hierarchy import linkage

    gene_cluster_labels = AgglomerativeClustering( n_clusters = N_CLUSTERS, 
                                                  affinity = 'precomputed', # precomputed: input distance 
                                                  linkage = 'average')\
                        .fit(data).labels_
    gene_cluster_labels_pd = pd.DataFrame( gene_cluster_labels, 
                 index = data.index, 
                columns = ['cluster'])
    return(gene_cluster_labels_pd, 1)    
'''
def eval_comp_cluster( data, gene_cluster_labels_pd ):
    """ Evaluator that helps choosing a clustering number.
    Method: compute 2 metrics. 1. intra-clustersimilarity vs. intra+inter. 2. intra-sim vs intra+similarity with second closet cluster
    ------
    Input:
    data: pairwise cosine similarity for all components.
    gene_cluster_labels_pd: df with gene as index and 'cluster' column as cluster index.
    ------
    Return:
    mean_dist_ratio
    mean_close2_ratio
    ------
    Usage:
       ## Given rcp_concat_z_topg
       X = sklearn.metrics.pairwise.cosine_similarity(rcp_concat_z_topg.T)
       x = pd.DataFrame(X, index = rcp_concat_z_topg.columns, columns = rcp_concat_z_topg.columns)
       data = x.T
       eval_comp_cluster( data, gene_cluster_labels_pd )     
    """
    complist = gene_cluster_labels_pd.index
    data = data.loc[ complist,complist ]
    dist_ratio_list = []
    close2_ratio_list = []
    # remove the largest cuz that would be noise
    #i_biggist_cluster = gene_cluster_labels_pd.reset_index().groupby('cluster').count().sort_values(by = 'index', ascending=False).index[0]
    #gene_cluster_labels_pd = gene_cluster_labels_pd[ gene_cluster_labels_pd['cluster']!=i_biggist_cluster]
    
    for cluster_i in sorted(gene_cluster_labels_pd['cluster'].drop_duplicates()):
        comp_incluster = gene_cluster_labels_pd[gene_cluster_labels_pd[ 'cluster']== cluster_i].index.tolist()
        if(len(comp_incluster)<=3): # skip cluster with less than 3 components
            continue
        intrasim = data.loc[comp_incluster, comp_incluster].values.flatten().mean()
        intersim = data.loc[data.index.difference(comp_incluster), comp_incluster].values.flatten().mean()
        dist_ratio = intrasim/(intrasim + intersim )
        dist_ratio_list.append(dist_ratio)
        ## the second closest cluster
        for cluster_j in sorted(gene_cluster_labels_pd['cluster'].drop_duplicates()):
            close2sim=100 
            if(cluster_i == cluster_j):
                continue
            comp_incluster_j = gene_cluster_labels_pd[gene_cluster_labels_pd[ 'cluster']== cluster_j].index.tolist()
            t = np.mean( data.loc[comp_incluster_j, comp_incluster].values.flatten() )
            if(t < close2sim):
                close2sim = t
                # 
        close2_ratio = intrasim/(intrasim + close2sim  )
        close2_ratio_list.append(close2_ratio)
    mean_dist_ratio = np.mean(dist_ratio_list)
    mean_close2_ratio = np.mean( close2_ratio_list)
    return(mean_dist_ratio,mean_close2_ratio, len(dist_ratio_list))
'''
'''
def eval_comp_cluster( data, gene_cluster_labels_pd ):
    """ Evaluator that helps choosing a clustering number.
    Method: compute 2 metrics. 1. intra-clustersimilarity vs. intra+inter. 2. intra-sim vs intra+similarity with second closet cluster
    ------
    Input:
    data: pairwise cosine similarity for all components.
    gene_cluster_labels_pd: df with gene as index and 'cluster' column as cluster index.
    ------
    Return:
    mean_dist_ratio
    mean_close2_ratio
    ------
    Usage:
       ## Given rcp_concat_z_topg
       X = sklearn.metrics.pairwise.cosine_similarity(rcp_concat_z_topg.T)
       x = pd.DataFrame(X, index = rcp_concat_z_topg.columns, columns = rcp_concat_z_topg.columns)
       data = x.T
       eval_comp_cluster( data, gene_cluster_labels_pd )     
    """
    complist = gene_cluster_labels_pd.index
    data = data.loc[ complist,complist ]
    dist_ratio_list = []
    close2_ratio_list = []
    # remove the largest cuz that would be noise
    #i_biggist_cluster = gene_cluster_labels_pd.reset_index().groupby('cluster').count().sort_values(by = 'index', ascending=False).index[0]
    #gene_cluster_labels_pd = gene_cluster_labels_pd[ gene_cluster_labels_pd['cluster']!=i_biggist_cluster]
    
    for cluster_i in sorted(gene_cluster_labels_pd['cluster'].drop_duplicates()):
        comp_incluster = gene_cluster_labels_pd[gene_cluster_labels_pd[ 'cluster']== cluster_i].index.tolist()
        if(len(comp_incluster)<=3): # skip cluster with less than 3 components
            continue
        intrasim = data.loc[comp_incluster, comp_incluster].values.flatten().mean()
        intersim = data.loc[data.index.difference(comp_incluster), comp_incluster].values.flatten().mean()
        dist_ratio = intrasim/ intersim 
        dist_ratio_list.append(dist_ratio)
        ## the second closest cluster
        for cluster_j in sorted(gene_cluster_labels_pd['cluster'].drop_duplicates()):
            close2sim= 0 
            if(cluster_i == cluster_j):
                continue
            comp_incluster_j = gene_cluster_labels_pd[gene_cluster_labels_pd[ 'cluster']== cluster_j].index.tolist()
            t = np.mean( data.loc[comp_incluster_j, comp_incluster].values.flatten() )
            if(t > close2sim):
                close2sim = t
                # 
        close2_ratio = intrasim/ (intrasim+close2sim  )
        close2_ratio_list.append(close2_ratio)
    mean_dist_ratio = np.mean(dist_ratio_list)
    mean_close2_ratio = np.mean( close2_ratio_list)
    return(mean_dist_ratio,mean_close2_ratio, len(dist_ratio_list))
'''

def eval_comp_cluster( data, comp_cluster_labels_pd, metric='precomputed' ):
    """ Evaluator that helps choosing a clustering number.
    data: distance matrix cuz metric='precomputed'
    """
    from  sklearn.metrics import  silhouette_score
    ## remove cluster with <=3 compoentns
    #t = comp_cluster_labels_pd.reset_index().groupby('cluster').count()
    #clusters_kept = t[t['index']>3].index.tolist()
    #comp_cluster_labels_pd = comp_cluster_labels_pd[comp_cluster_labels_pd['cluster'].isin(clusters_kept)]

    labels = comp_cluster_labels_pd[ comp_cluster_labels_pd.columns[0] ].values
    if(metric=='precomputed'):
        data = data.loc[comp_cluster_labels_pd.index, 
                                    comp_cluster_labels_pd.index]
    res = silhouette_score( data ,  labels, metric= metric )

    return(res)

def jaccard_dist( x, y, topn = 50):
    """ compute jaccard index between two pd series"""
    if(isinstance(y, pd.Series)):
        y=y[y>0]
        yset = y.sort_values(ascending=False)[:topn].index.tolist()
    elif( isinstance(y, list)):
        yset = y[:topn]
    x= x[x>0]; 
    xset = x.sort_values(ascending=False)[:topn].index.tolist()
    
    b = set(xset).union(set(yset))
    a = ( set(xset) - set(yset) ).union( set(yset) - set(xset))
    #    print(xset,yset)
    return( len(a) /len(b)  )

def trim_comp_cluster( comp_cluster_labels_pd, min_cluster_size = 5):
    t = comp_cluster_labels_pd.reset_index().groupby('cluster').count()
    clusters_kept = t[t['index']>=min_cluster_size].index.tolist()
    comp_cluster_labels_pd = comp_cluster_labels_pd[comp_cluster_labels_pd['cluster'].isin(clusters_kept)]
    # print(len(comp_cluster_labels_pd['cluster'].drop_duplicates()))
    return(comp_cluster_labels_pd)


def cluster_louvain( rcp_concat_z_topg, resolution = 1.25,min_cluster_size = 5 ):
    """ Graph clustering for either gene or components
    defalt: cosine similarity"""

    from sknetwork.clustering import Louvain
    from sklearn.metrics.pairwise import cosine_similarity
    gg_adjacency = cosine_similarity( rcp_concat_z_topg.fillna(0) )
    df = pd.DataFrame( np.exp(gg_adjacency), index = rcp_concat_z_topg.index, columns= rcp_concat_z_topg.index)

    # louvain clustering
    louvain = Louvain( resolution = resolution, random_state = 2)
    labels= louvain.fit_transform( df.values )
    cluster_df = pd.DataFrame( {'cluster' :labels}, index = df.index)

    # remove cluster size <=5
    t = cluster_df.reset_index().groupby('cluster').count()
    t= t[t['index']>= min_cluster_size ]
    n_cluster = t.shape[0]
    cluster_df = cluster_df[ cluster_df['cluster'].isin( t.index  )]

    # create graph for plot.
    #G, links_filtered = graphfun._construct_graph( df )
    #[colors_nodes, colors_nodes_map] = graphfun.List2Colors( labels, cmap = matplotlib.cm.tab20c, return_valmap=True, colorexpand=True, alpha=0.7) # best colormap! rainbow!
    
    return(cluster_df, n_cluster)




######### Save Meta-component matrix and ordered list ########


def mec_gene_full_table( comp_cluster_labels_pd , rcp_concat_z_topg):
    metacomp_rank_lst = []
    metacomp_weightaverage_lst = []
    for cluster_i in sorted( comp_cluster_labels_pd['cluster'].drop_duplicates()) :

        ics = comp_cluster_labels_pd[comp_cluster_labels_pd['cluster'] == cluster_i].index.tolist()
        components_cluster = rcp_concat_z_topg[ ics ]
        ranksum = components_cluster.rank(ascending=False).applymap(lambda t: np.log(t+1)).sum(axis=1).sort_values()
        # Gene rank of the mega-component
        metacomp_rank = pd.DataFrame( ranksum.rank(), columns = [cluster_i] )
        metacomp_weightaverage = pd.DataFrame( components_cluster.mean(axis=1), columns = [cluster_i] )
        metacomp_rank_lst.append(metacomp_rank)
        metacomp_weightaverage_lst.append(metacomp_weightaverage)
    metacomp_rank = pd.concat( metacomp_rank_lst, axis=1)
    metacomp_weightaverage = pd.concat( metacomp_weightaverage_lst, axis=1)

    return( metacomp_rank, metacomp_weightaverage)

def mec_gene_top( metacomp_weightaverage  ):
    """
    Save the table for results
    Each row is a MeC, Top gene and top tf (calls gene module)
    """

    module = pd.DataFrame(index = metacomp_weightaverage.columns)

    alltfs = getgene.getalltfs() #[todo] make TF work for mouse
    for cluster_i in metacomp_weightaverage.columns :
        tmp = metacomp_weightaverage[cluster_i].sort_values(ascending=False)
        tfs = tmp[:100].index.intersection(alltfs)
        #print( cluster_i,'\t',','.join( tfs[:10].tolist()),'\t',','.join( tmp[:10].index.tolist()) )
        module.loc[cluster_i, 'TopTF'] = ','.join( tfs[:10].tolist())
        module.loc[cluster_i, 'TopGene20'] = ','.join( tmp[:20].index.tolist())
        module.loc[cluster_i, 'TopGene100'] = ','.join( tmp[:100].index.tolist())
        module.loc[cluster_i, 'TopGeneSig5'] = ','.join( tmp[tmp>5].index.tolist())
    return(module)

def mec_score_topg( mecscore, meclst, meclst_gene_col = 'TopGeneSig5'  ):
    """
    only keep topgene in the score matrix. masking useful for tsmec projection.
    mecscore columns are string. meclst index are integers.
    """
    g = [] # gene order
    mecscore_topg_mask = pd.DataFrame(False, index= mecscore.index, columns = mecscore.columns )
    for i in meclst.index:
        colname = str(i)
        gs = meclst.loc[i, meclst_gene_col].split(',')
        # gene order
        g = g + [ t for t in gs if t not in g ]
        mecscore_topg_mask.loc[ gs, colname ] = True
    mecscore_topg = mecscore[ mecscore_topg_mask ].loc[g].fillna(0)
    return(mecscore_topg)

def mec_maxnorm( dat ):
    """ normalization using max value, independently for + and - weights. 
    ----
    E.g.
    mec_1score = mec_maxnorm( mec_score )"""
    dat1 = dat[dat>0]/ dat.max()
    dat2 = -1 *dat[dat<0]/ dat.min()
    dat = dat1.fillna(0)+dat2.fillna(0)
    return(dat)


def mec_gene_lstfiles( mec_topg , OUTDIR, genecol = 'TopGene100', topn = 30):
    """ save seperate gene lsit for interpretation"""
    for i in mec_topg.index:
        filename = 'mec.'+str(i)+'.genelist.txt' 
        glst = pd.DataFrame( {'g': mec_topg.loc[i, genecol ].split(',')[:topn]  })
        glst.to_csv( os.path.join( OUTDIR, filename ) , sep='\t', index=False,header=False)

def load_mec_precomputed():
    """Load pre-computed Meta-component matrix and ordered list
    #Look for precomputed files in config.SCMECDIR
    MeC_allgene_ranksum-ranks.txt
    MeC_allgene_average-weights.txt
    MeC_topgene.txt
    MeC_TopGeneSig5_average-weights.txt
    -----
    Usage:
    mecres = metacomp.load_mec_precomputed()
    mec_topg = mecres['mec_topg']
    mec_rank = mecres['mec_rank']
    mec_score = mecres['mec_score']
    mec_score_topg = mecres['mec_score_topg']
    mec_cluster = mecres['mec_cluster']
    """
    import config
    mec_rank = pd.read_table( os.path.join( config.SCMECDIR, 'MeC_allgene_ranksum-ranks.txt') , index_col = 0)
    mec_score = pd.read_table( os.path.join( config.SCMECDIR, 'MeC_allgene_average-weights.txt') , index_col = 0)
    mec_score_topg = pd.read_table( os.path.join( config.SCMECDIR, 'MeC_TopGeneSig5_average-weights.txt') , index_col = 0)
    mec_topg = pd.read_table( os.path.join( config.SCMECDIR, 'MeC_topgene.txt') , index_col = 0)
    mec_anno = pd.read_table( os.path.join( config.SCMECDIR, 'MeC_anno.txt') , index_col = 0)
    mec_cluster = pd.read_table( os.path.join( config.SCMECDIR, 'comp_cluster_labels.txt' ), index_col = 0  )
    return({'mec_rank': mec_rank, 
            'mec_score': mec_score, 
            'mec_score_topg': mec_score, 
            'mec_topg' : mec_topg,  # can also obtain from mecmapper.load_mec(mode='list)
            'mec_cluster': mec_cluster,
            'mec_anno': mec_anno,
           })



def load_mecname( mode='mecnamedict' ):
    """
        ----Can use to generate config scorectcols
        mecname = metacomp.load_mecname( mode ='table')
        t = mecname[mecname['UseForCellTypeAnno']==1]
        t['score_id'] = 'score_'+ t['MeC_id'].str.split('_').str.get(1)
        scorectcols = t['score_id'].values.tolist()
    """
    import config
    mecname = pd.read_table(os.path.join(config.SCMECDIR, 'MeC_anno_name.txt'))
    if(mode=='table'):
        mectable=mecname.copy()
        #mectable['color_level0'] = mectable['MajorLineage_level0'].apply(lambda t: config.level0colordict[t])
        return(mectable)        
    else:
        mecname['scorecol']=config.scorecols
        mask = (mecname['Annotation'].isna())
        mecname.loc[mask,'Annotation']=mecname.loc[mask,'MeC_id' ]
        mecnamedict = mecname.set_index('scorecol')[['Annotation']].to_dict()['Annotation']

    if(mode=='mecnamedict'):
        return(mecnamedict)
    elif(mode =='meciddict'):
        meciddict={}
        _=[meciddict.update({mecnamedict[t]:t}) for t in mecnamedict.keys()]
        return(meciddict)
    print('modes: mecnamedict, table , meciddict')
        

def load_meclineage():
    import config
    """ DNU.
    Pre-assigned , per mec, lineage assignment
    """
    meclineage = pd.read_table(os.path.join(config.SCMECDIR, 'MeCname_lineage_types.txt'))
    # mec with no lineage TODO: can return it
    mec_nolineage = meclineage[meclineage['MajorLineage'].isin([None,'Signaling'])]['name'].values.tolist()
    mec_islineage = meclineage[meclineage['MajorLineage'].isin(config.ct_lineages)]['name'].values.tolist()
    return(meclineage ,mec_nolineage,mec_islineage)


def print_gene_marker3():
    """ automatically generated per signature marker , for manual pick. 
    import json
    with open('../analysis/20211001_scICA/gene_marker_selected_for_plot_known.json','r' ) as rp:
        markerjson = json.load(rp)
    """
    import json
    mecnamedict = load_mecname()
    mecres = load_mec_precomputed()
    mec_topg = mecres['mec_topg']
    mecname = load_mecname(mode='table')
    # follow major lineage order
    mecname = mecname[mecname['Type']!='u']
    mecname = mecname.sort_values(by=['MajorLineage_level0','Type'])
    markerd={}
    for i in mecname.index:
        k = mecname.loc[i,'Annotation']
        v = mec_topg.loc[i,'TopGene20'].split(',')[:1]
        markerd.update({k:v})
    with open('../analysis/20211001_scICA/gene_marker1_selected_for_plot.json','w') as fp:
        json.dump(markerd, fp, indent=4, )
    print('Output to ../analysis/20211001_scICA/gene_marker_selected_for_plot.json ')
    
"""
"""