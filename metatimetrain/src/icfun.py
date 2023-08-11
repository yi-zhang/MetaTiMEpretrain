import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sklearn
import sklearn.decomposition as decomp
import scipy.stats as ss
import config

########## run low-dim decomposition  ##########


def runIC(X, n_components, random_state = 42 ):
    X -= X.mean(axis=0)
    #[original: fixed with FastICA_correction. research into this more.
    #]ica_correction = decomp.FastICA_correction( n_components = n_components, whiten = True, random_state = random_state , max_iter = 800)
    ica_correction = decomp.FastICA( n_components = n_components, whiten = True, random_state = random_state , max_iter = 800)
    ica_correction.fit( X )
    components = pd.DataFrame( ica_correction.components_.T , index = X.columns )
    #components.to_csv( 'SCC.sc.CD8.c'+str(n_components)+'.txt', sep='\t')
    return(ica_correction, components)

def runNMF( X, n_components, random_state = 42):
    nmf = decomp.NMF( n_components = n_components ,random_state = random_state )
    nmf.fit(X)
    components = pd.DataFrame( nmf.components_.T , index = X.columns )
    return(nmf, components)

########## IC rotation ##########


def skew_outlier(tmp):
    
    ztmp  = (tmp-tmp.mean())/tmp.std()
    return( ss.skewtest(ztmp))

def ic_rotate( components, n_components ):
    """ Rotate components based on purely weights 
    ====
    Calls: skew_outlier
    """
    skewtests = [skew_outlier( components[ col ].values) for col in components.columns ]  
    skews = np.array( [t[0] for t in skewtests ] )
    components = components * (  ( skews >0 )*2-1  )
    return( components, skewtests )

########## IC reproducible ##########

# concatenate ics and QC filtering.

def qc_comp_ztrans(components_all_concat ):
    """ z-transform and clean the compoennts. 
    # z-transform, an empirical gene cutoff strategy.
    ## processing: scale weight using std. keep those out of 1std, remove genes that do not have any <1std, subset those genes. 
    
    """
    components_all_concat_z = (components_all_concat.fillna(0) ) / components_all_concat.fillna(0).std()
    components_all_concat_z_filtered = components_all_concat_z[components_all_concat_z.abs()>1]
    components_all_concat_z_filtered = components_all_concat_z_filtered[ components_all_concat_z_filtered.isna().sum(axis=1) < components_all_concat_z_filtered.shape[1] ]

    components_all_concat_z = components_all_concat_z.loc[ components_all_concat_z_filtered.index].fillna(0)

    components_all_concat_z_filtered = components_all_concat_z[components_all_concat_z.abs()>1]
    components_all_concat_z_filtered = components_all_concat_z_filtered.fillna(0)
    return(components_all_concat_z_filtered)

def qc_comp_topg( components_all_concat_z_filtered , topn = 200, polar='both'):  
    """ keep top genes with extreme Component Score only.  An empirical gene cutoff strategy.
    """
    ## keep genes within rank 200 for ic clustering. keep genes showing up in> 3 components.
    if( polar == 'both'):
        components_all_concat_z_filtered_topg= components_all_concat_z_filtered[ (components_all_concat_z_filtered.rank(ascending=False)< topn ) | \
                                    (components_all_concat_z_filtered.rank(ascending=True)< topn )]
    elif( polar == 'pos'):
        components_all_concat_z_filtered_topg= components_all_concat_z_filtered[ (components_all_concat_z_filtered.rank(ascending=False)< topn ) ]
#    elif( polar == 'neg'):
#        components_all_concat_z_filtered_topg= components_all_concat_z_filtered[  (components_all_concat_z_filtered.rank(ascending=True)< topn )]
       
    components_all_concat_z_filtered_topg = components_all_concat_z_filtered_topg[(components_all_concat_z_filtered_topg.isna().sum(axis=1))< components_all_concat_z_filtered_topg.shape[1]-3]
    components_all_concat_z_filtered_topg = components_all_concat_z_filtered_topg.fillna(0)
    components_all_concat_z_filtered = components_all_concat_z_filtered_topg.copy()

    return( components_all_concat_z_filtered )

def qc_comp_maxcor( components, ic_sims, CORRCUT = 0.3 ):
    """ Only keep components with correlation >=0.3 with at least one compoennt from another cohort. """

    #simarray = np.random.choice( ic_sims.values.flatten(), size = 10000 )
    ##sns.histplot( simarray)
    #cut = 2*simarray.std() 
    ##np.quantile(simarray, q = 0.999)
    #print('cut: ', cut)

    ickept = []

    for c in ic_sims.columns:
        cohort = c.split('.')[0]
        index_othercohort = ic_sims.index[ ic_sims.index.str.split('.').str.get(0)!=cohort ] 
        maxcor = ic_sims.loc[ index_othercohort, c].sort_values(ascending=False)[1]
        if(maxcor>= CORRCUT ):
            ickept.append(c)
    # use these filtered components: ickept
    components = components[ ickept ]
    return(components)
    

def concat_ic( compDIR = config.SCICDIR , OUTDIR = config.SCICDIR, 
                readingsuffix = '.sc.rotate.c100.txt', 
                outfile = '', #'components_all_concat.txt' ,
                filter = False ):
    """ Concatenate all components after rotation.
    compDIR: place to save individual gene by comp matrix.
    suffix: file name fixed suffix
    ----
    return: concatenated comp table
    ----
    Calls: qc_comp_ztrans
    
    """

    files = os.listdir( compDIR )
    files = [t for t in files if  readingsuffix in t]
    pdlst = []
    for file in files:
        name = file.split('.')[0] 
        dat = pd.read_table( os.path.join( compDIR , file) , index_col = 0 ) 
        dat.columns = name + '.'+ dat.columns
        pdlst.append( dat)
        
    components_all_concat = pd.concat(pdlst, axis=1 )
    components_all_concat.loc[components_all_concat.isna().sum(axis=1)<components_all_concat.shape[1]].shape
    # filter
    if(filter):
        print( "[Z score tranform and filtering]")
        components_all_concat =  qc_comp_ztrans(components_all_concat )
    # save
    #compall_OUTfile = os.path.join(OUTDIR, outfile )
    #components_all_concat.to_csv( compall_OUTfile ,sep='\t')
    #print(f'[Log] Saved concatenated component matrix in :\n {compall_OUTfile}')
    return(components_all_concat)

def comps_similarity( components_all_concat_z_filtered_topg , file = 'sim_pairwise_components.txt', OUTDIR = config.SCICDIR):
    """ 
    ---
    generating filtered similarity:
    components_all_concat_z_filtered_topg = qc_comp_topg(components_all_concat )
    ----
    Calls: qc_comp_topg
    """

    X = components_all_concat_z_filtered_topg.T.values
    ic_sims = pd.DataFrame(X, index = components_all_concat_z_filtered_topg.columns.tolist() ).T.corr(method='pearson')
    if(file):
        print(f'[Save file] Saving to {file} in {OUTDIR}')
        ic_sims.to_csv( os.path.join(OUTDIR, file),sep='\t')
    return( ic_sims )
        
def qc_comp_highsim( components, ic_sims ):  
    """ keep top components with high similarty (correlation > 0.1)
    """
    
    

    return( components )



########## [pre-computed] pre-computed component objects #########
def load_components_precomputed( compconcatDIR = config.SCICCONCATDIR, compDIR = config.SCICDIR ):
    """
    #Concatenate all avaiable components, and QC.
    #Look for precomputed files in config.SCICDIR
    components_all_concat.txt
    components_topg_concat.txt
    components_similarity.txt
    rcp_z.txt
    rcp_topg.txt
    -----
    Usage:
    ic100res = icfun.load_components_precomputed()
    comp_concat_z = ic100res['comp_concat_z']
    comp_concat_z_topg = ic100res['comp_concat_z_topg']
    ic_sims = ic100res['ic_sims']
    rcp_concat_z = ic100res['rcp_concat_z']
    rcp_concat_z_topg = ic100res['rcp_concat_z_topg']
    """
    ## compute
    comp_concat_z_file = os.path.join(compconcatDIR, 'components_all_concat.txt' )
    comp_concat_z_topg_file = os.path.join(compconcatDIR, 'components_topg_concat.txt' )
    ic_sims_file = os.path.join(compconcatDIR, 'components_similarity.txt' )
    rcp_z_file = os.path.join(compconcatDIR, 'rcp_z.txt' )
    rcp_topg_file = os.path.join(compconcatDIR, 'rcp_topg.txt' )

    if( (not os.path.isfile( comp_concat_z_file) ) | \
        (not os.path.isfile( comp_concat_z_topg_file) ) | \
        (not os.path.isfile( ic_sims_file) ) 
            ):
        print('[Log] Concatenating all components')
        comp_concat = concat_ic( compDIR = compDIR , OUTDIR = compconcatDIR )
        print('[Log] ztransform all components')
        comp_concat_z = qc_comp_ztrans(comp_concat)
        print('[Log] Getting top genes on all components')
        comp_concat_z_topg = qc_comp_topg(comp_concat_z )
        print('[Log] Eval similarity of all components')
        ic_sims = comps_similarity(comp_concat_z_topg)
        print('[Log] Only keep components with correlation >=0.3 with at least one compoennt from another cohort.  ')
        rcp_z = qc_comp_maxcor(comp_concat_z_topg, ic_sims)
        print('[Log] Getting top genes on selected components')
        rcp_topg = qc_comp_topg( rcp_z)

        comp_concat_z.to_csv( comp_concat_z_file , sep='\t')
        comp_concat_z_topg.to_csv(comp_concat_z_topg_file , sep='\t')
        ic_sims.to_csv( ic_sims_file, sep='\t')
        rcp_z.to_csv( rcp_z_file , sep='\t')
        rcp_topg.to_csv( rcp_topg_file , sep='\t')

        print('comp_concat.shape', comp_concat.shape)
        print('comp_concat_z.shape', comp_concat_z.shape)
        print('comp_concat_z_topg.shape', comp_concat_z_topg.shape)
        print('rcp_z.shape', rcp_z.shape)
        print('rcp_topg.shape', rcp_topg.shape)

    else:
        comp_concat_z = pd.read_table( comp_concat_z_file , index_col = 0)
        comp_concat_z_topg = pd.read_table( comp_concat_z_topg_file , index_col = 0)
        ic_sims = pd.read_table( ic_sims_file , index_col = 0)
        rcp_z = pd.read_table( rcp_z_file , index_col = 0)
        rcp_topg = pd.read_table( rcp_topg_file , index_col = 0)

    return({'comp_concat_z': comp_concat_z, 
            'comp_concat_z_topg': comp_concat_z_topg, 
            'ic_sims' : ic_sims, 
            'rcp_concat_z':rcp_z,
            'rcp_concat_z_topg': rcp_topg,
           })




########## IC interpreting 

def comp2dict(complist):
    """ reformatting component list based on dict"""
    #complist = component_concat_corr.index[np.where(gene_cluster_labels == 13)].tolist()
    res = {}
    _=[ res.update({ t.split('.')[0] : t.split('.')[1] }) for t in complist ]
    return(res)

def gene_clusteri( cluster_i, comp_cluster_labels_pd, component_concat, topn=20, savefile =None, saveweighted = None, saverowf = None, polar= 'positive', printing = True, geneonly = False):
    """ For inner check. Simple printing of gene members of clustered programs"""
    #cluster_i = 5
    clucol = comp_cluster_labels_pd.columns[0]
    complist = comp_cluster_labels_pd[comp_cluster_labels_pd[clucol] == cluster_i ]
    #complistdict = comp2dict(complist)
    comp_clustered = complist.index.tolist()
    if(len(comp_clustered) ==0 ):
        return(None)
    if( polar == 'positive'):
        #t = component_concat[comp_clustered].rank(ascending=False).mean(axis=1).sort_values(ascending=True)
        t = component_concat[comp_clustered].mean(axis=1).sort_values(ascending=False)
        genes = t[: topn].index.tolist()
        #weightedgenes = component_concat[comp_clustered].mean(axis=1).sort_values(ascending=False)[: topn]
        weightedgenes = component_concat[comp_clustered].mean(axis=1).loc[genes]
    elif( polar == 'negative'):
        genes = component_concat[comp_clustered].mean(axis=1).sort_values(ascending=True)[: topn].index.tolist()
        weightedgenes = component_concat[comp_clustered].mean(axis=1).sort_values(ascending=True)[: topn]
    #print(','.join(comp_clustered))
    if(printing == True):
        print('Cluster',cluster_i, '\t', len(comp_clustered),'\t',
        ','.join(genes))
    if(savefile):
        pd.DataFrame({ 'gene': genes} ).to_csv(savefile, header=False, index=False)
    if(saverowf):
        if(geneonly == True):
            saverowf.write( 
             ','.join(genes)+ '\n'
            ) 
        else:
            saverowf.write( 
             'Cluster' +str(cluster_i)+ '\t'+str( len(comp_clustered)) + '\t' + ','.join(genes)+ '\n'
            )
    if(saveweighted):
        pd.DataFrame( weightedgenes, columns = ['sig'] ).to_csv(saveweighted, sep='\t', header=False, index=True)
    return( weightedgenes )

def genelistrowwise( genelistfiles, genelistrowfile ):
    """ Reformat gene list files into one file where each row is a gene list"""
    wf = open(genelistrowfile,'w')

    for file in genelistfiles:
        dat = pd.read_csv(file, names=['gene'])
        wf.write( ','.join( list(dat['gene'].values) ) + '\n' )
    wf.close()
    
    
########## plotting ##########
def markscatter(dat, xcol = None, ycol = None):
    import adjustText
    xcol = dat.columns[0]
    ycol = dat.columns[1]
    dat = dat.loc[ dat.rank(ascending=False).sum(axis=1).sort_values(ascending=True).index ].fillna(0)


    markgenes = dat[:40].index.tolist()
    markgenes = markgenes + dat[ xcol ].sort_values(ascending=False)[:10].index.tolist()
    markgenes = markgenes + dat[ xcol ].sort_values(ascending=True)[:5].index.tolist()
    markgenes = markgenes + dat[ ycol ].sort_values(ascending=False)[:5].index.tolist()
    markgenes = markgenes + dat[ ycol ].sort_values(ascending=True)[:5].index.tolist()
    print(' '.join( markgenes))
    markgenes = list(set(markgenes))

    sns.scatterplot(data=dat, x = xcol , y = ycol)
    #plt.gca().set_ylim(bottom=-2)
    texts = []
    for gene in markgenes:
        texts.append(plt.text(x = dat.loc[gene][ xcol ], y = dat.loc[gene][ ycol ], s = gene))
    adjustText.adjust_text( texts , arrowprops=dict(arrowstyle='-', color='pink'))



def negative_entropy( w ):
    import numpy
    '''
    Calculate negative entropy of the component w
    '''
    w_mean = np.mean(w)
    w_std = np.std(w)
    y = (w-w_mean) / w_std
    Egv = -1 / 2**0.5
    gy = - np.exp(-y**2/2)
    Egy = np.mean(gy)
    negentropy = np.mean( (Egy-Egv)**2 )
    return( negentropy )

## aggregate a list of component tables
def agg_ic_files(cohort_ic_files, icDIR, suffix ='.ic50.txt'):
    ic_dict = {}
    ic_all_list = []
    for file in cohort_ic_files:
        thispat = file.replace( suffix ,'')
        ic = pd.read_table(os.path.join( icDIR, file ), index_col = 0)
        ic_dict.update( { 
            thispat: ic
        })
        ic.columns = [thispat+':'+str(t) for t in ic.columns]
        ic_all_list.append(ic)
    ic_all = pd.concat(ic_all_list, axis=1)
    return(ic_all)


