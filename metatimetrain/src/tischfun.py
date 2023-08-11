import numpy as np
import scanpy as sc
import pandas as pd
import os

### pre-processing adata. tisch style. fun done.

def adatapp(
    adata_input, 
    mode='pp', # or  'umap' for just computing umap visualization embedding
    random_state = 42,
    MAX_MITO_PERCENT=5,
    MIN_GENES = 500,
    MIN_COUNTS = 1000,
    MIN_CELLS = 5
):
    
    """ 
    Pre-processing. if adata.X has integer record, or adata has 'counts' layer, go through standard preprocessing.
    remove cells with <500 genes and <1000 counts. remove genes with min cells < 5. 
    mitocondrial <5%. 
    normalize to 1e6. log. pca, umap. normalized value savedt to 'norm_data' layer
    ----
    Notes: human only
    -----
    Usage:
    adata = adatapp(adata )
    """

    adata = adata_input.copy()
    #print( ' preprocessing ')
    if(mode=='pp'):
        
        if('counts' in list(adata.layers)):
            # if 'counts' in layers, redo log,pca,neibours
            adata.X = adata.layers['counts']
            print( '[Log] use count layer ')
            adata.var_names_make_unique()
            adata.obs['n_counts'] = adata.X.sum(1).A1
            adata.obs['log_counts'] = np.log(adata.obs['n_counts'])
            adata.obs['n_genes'] = (adata.X > 0).sum(1).A1

            sc.pp.filter_cells(adata, min_genes = MIN_GENES)  
            sc.pp.filter_cells(adata, min_counts = MIN_COUNTS )   
            sc.pp.filter_genes(adata, min_cells = MIN_CELLS) 

            adata.var['mt'] = adata.var_names.str.startswith('MT-') | adata.var_names.str.startswith('mt-')
            sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
            adata = adata[adata.obs.pct_counts_mt < MAX_MITO_PERCENT , :]
            sc.pp.normalize_total(adata, target_sum = 1e4)  # match NormalizeData in Seurat
            sc.pp.log1p(adata) # checked the values are logged.
            sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=4000 ) 
            sc.tl.pca(adata,  svd_solver='arpack')
            sc.pp.neighbors(adata )
            sc.tl.umap(adata, random_state = random_state )

        elif( np.sum( np.sum(adata.X[:,0:100] , axis=0) %1 ) == 0 ): # if top 100 cells all mod 0, then assume the .X has integer
            # if .X has interger, redo log, pca, neighbours
            adata.layers["counts"] = adata.X.copy()
            print( '[Log] adata.X is count ')

            adata.var_names_make_unique()
            adata.obs['n_counts'] = adata.X.sum(1).A1
            adata.obs['log_counts'] = np.log(adata.obs['n_counts'])
            adata.obs['n_genes'] = (adata.X > 0).sum(1).A1

            sc.pp.filter_cells(adata, min_genes = MIN_GENES) 
            sc.pp.filter_cells(adata, min_counts = MIN_COUNTS )
            sc.pp.filter_genes(adata, min_cells = MIN_CELLS) 

            adata.var['mt'] = adata.var_names.str.startswith('MT-') | adata.var_names.str.startswith('mt-')
            sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
            adata = adata[adata.obs.pct_counts_mt < MAX_MITO_PERCENT, :]
            sc.pp.normalize_total(adata, target_sum = 1e4)  # match NormalizeData in Seurat
            sc.pp.log1p(adata) # checked the values are logged.
            sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=4000 )
            sc.tl.pca(adata,  svd_solver='arpack')
            sc.pp.neighbors(adata )
            sc.tl.umap(adata, random_state = random_state)
            
        else: # directly use the values if X not continuous: already normalized
            # if .X is continuous, use as is
            print( '[Log] adata.X is continuous ')
            adata.var_names_make_unique()
            adata.layers["norm_data"] = adata.X
            sc.pp.filter_cells(adata, min_genes = MIN_GENES) 
            sc.pp.filter_genes(adata, min_cells = MIN_CELLS) 
            sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=4000 )
            sc.tl.pca(adata,  svd_solver='arpack')
            sc.pp.neighbors(adata )
            sc.tl.umap(adata, random_state = random_state)


    elif(mode=='umap'):
        """ just for computing umap visulization embedding"""
        sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=4000) 
        sc.tl.pca(adata)
        sc.pp.neighbors(adata)
        sc.tl.umap(adata, random_state = random_state)

    # extra 
    adata.obs.columns = [t.replace('assign.','assign_') for t in adata.obs.columns]
    
    return(adata)


def readTischAnn(
    name, 
    preprocessing = False,
    suffix = '.h5ad',
    DATADIR = '',  
    MAX_MITO_PERCENT=5,
    MIN_GENES = 500,
    MIN_COUNTS = 1000,
    MIN_CELLS = 5
):
    """ 
    Read in a cohort as scanpy object.
    ----
    Notes: preprocessing: call adatapp.
    -----
    Usage:
    adata = readTischAnn(name , preprocessing = False, 
                               DATADIR = '/liulab/yiz/DATA/scRNA_processed/scanpy/pp/',
                                suffix = '_res.pp.h5ad')
    """
    #scdata_file = '/liulab/yiz/DATA/scRNA_processed/scanpy/'+name+'_res.h5ad'
    scdata_file = os.path.join( DATADIR , name + suffix )
    adata = sc.read(scdata_file)

    if ( preprocessing == True ):
        #if 'norm_data' not in adata.layers:
        ## do preprocessing from raw
        print( ' Redo preprocessing '+name )
        adata = adatapp(adata, 
                        mode='pp', # or  'umap' for just computing umap visualization embedding
                        random_state = 42,
                        MAX_MITO_PERCENT=MAX_MITO_PERCENT,
                        MIN_GENES = MIN_GENES,
                        MIN_COUNTS = MIN_COUNTS,
                        MIN_CELLS = MIN_CELLS
                        )
    else:
        print( ' Loaded '+ name)
    return(adata)


###########

### tischdata function  ###

import config

def allTischCohorts(
    mode = 'h5', 
    TISCHWEBDIR = config.TISCHWEBDIR, 
    SCPPDIR = config.SCPPDIR,
    SCDIR = config.SCDIR,
):
    
    #Get a list of human  TISCH cohort names 
    #mode in ('pp', 'web' ): 'web' has _Data meta , 'pp' from post-processed h5ad
    
    if(mode == 'web'):
        allnames = [ t.replace('_Data', '') for t in os.listdir( TISCHWEBDIR ) if ( t[-5:] == '_Data'  and  'mouse' not in t ) ]
    elif( mode == 'pp' ):
        allnames = [ t.replace('_res.pp.h5ad', '') for t in os.listdir( SCPPDIR ) if ( t.split('.')[-1] == 'h5ad'  and  'mouse' not in t ) ]
    elif( mode == 'h5' ):
        allnames = [ t.replace('_res.h5ad', '') for t in os.listdir( SCDIR ) if ( t.split('.')[-1] == 'h5ad'  and  'mouse' not in t ) ]
    ## load public-only data
    private_only = ['inhouse_bz_grmix_2h48h','inhouse_bz_grmix_2h48h_res_integrated.h5ad']
    allnames = list(set(allnames)-set(private_only))
    ## human only
    allnames = [t for t in allnames if 'mouse' not in t]
    return(allnames)

### Single cohort tisch fun ###

def readmeta(
    cohort, 
    TISCHWEBDIR  = config.TISCHWEBDIR
):
    cohortDIR = cohort + '_Data'
    filename = cohort + '_CellMetainfo_table.tsv'
    filepath = os.path.join( TISCHWEBDIR , cohortDIR, filename)
    data = pd.read_csv(filepath, sep = '\t', index_col = 0)
    return(data)

### special: fix some real adata

def fixadataobs(
    adata, name
):
    # Special treatment for the clinical info in adata
    # Return: cleaned adata.obs
    # Calls: readmeta
    # ----
    if ( name == 'SKCM_GSE120575_aPD1aCTLA4' ): # patient id has Pre_ or Post_.
        ## this meta misses pre and post information.
        adata.obs['Timing'] = adata.obs['patient'].str.split('_').str.get(0) 
    if (name == 'BCC_GSE123813_aPD1'): # extra clinical info
        # fix meta data
        meta = readmeta(name, config.TISCHWEBDIR )
        adata.obs = adata.obs.merge( meta[ ['Celltype (malignancy)','Celltype (minor-lineage)', 'Celltype (original)', 'Patient','Treatment'] ] , 
                left_index = True, right_index = True , how = 'left')

    return(adata)

def isRawCount( ds ):
    # Tell whether the matrix is raw  count or float 
    # ------
    # ds: loom, or just a matrix 2d. 
    
    if( ds[:100,:100].sum() - int( ds[:100,:100].sum() ) == 0 ) :
        datatype = 'rawcount'
        datatype = 1
    else:
        datatype = 'notcount'
        datatype = 0
    return(datatype)

