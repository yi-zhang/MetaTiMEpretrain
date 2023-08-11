#!/usr/bin/env python
# coding: utf-8

from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

######### Generating colored df ########

def colorCluster( cluster_labels_pd, colorcol = 'cluster', comp_order = None,
palette=None, ):
    """ color cluster for a column in dataframe used for clustermap
    ----
    Usage:
    cluster_labels_pd, color_dict = colorCluster( cluster_labels_pd )
    # then use in clustermap:
    #row_colors = cluster_labels_pd[ colorcol + '_color' ]
    """
    if not colorcol:
        colorcol = cluster_labels_pd.columns[0]
    if comp_order == 'reorder':
        comp_order = cluster_labels_pd.sort_values( by = colorcol ).index
        cluster_labels_pd = cluster_labels_pd.loc[ comp_order ]
        
    N_CLUSTERS = cluster_labels_pd[ colorcol].drop_duplicates().shape[0]
    
    cluster_labels = cluster_labels_pd[ colorcol ].values
    color_dict = dict(zip(np.unique(cluster_labels), sns.color_palette( palette=palette, n_colors = N_CLUSTERS)))
    row_colors = [color_dict[cl] for cl in cluster_labels] 
    cluster_labels_pd[colorcol+'_color'] = row_colors
    return( cluster_labels_pd, color_dict  )

######### Generate heatmap line positions ########

def get_line_pos(gene_meta_pd, comp_meta_pd ):
    gene_line_j = gene_meta_pd.groupby( 'mec_i').max()['j'].values
    gene_line_j=gene_line_j+1
    comp_line_j = comp_meta_pd.groupby( 'mec_i').max()['j'].values
    comp_line_j=comp_line_j+1 
    return( gene_line_j,comp_line_j  )


######### Formating text blocks ########

def mecFormatGeneBlock( mec_order, gene_labels_pd, mec_anno_name_file = '../analysis/20211001_scICA/MeC/MeC_anno_name.txt' ,
                        genes_newline_str_len = 16, 
                        annotext_newline_str_len = 16,
                        Ntopgenes = 2 ):
    """ formatting the marked genes
    Usage
    -----
    mec_gene_block = mecFormatGeneBlock( mec_order, gene_labels_pd )
    ----
    mec_order: e.g. ['0','1']
    gene_labels_pd : df, [index]gene, columns[mec_i (cluster id), mec_i_color, j (gene position integer)] 
    mec_anno_name_file
    ----
    Return:
    Dict with color, text, block name
    """
    i = 0
    mec_block = dict()
    mec_name = pd.read_table(mec_anno_name_file, index_col = 0 )
    #mec_name.loc[mec_name['Annotation'].isna(), 'Annotation'] = mec_name.index
    for mec_i in mec_order:
        # the fillcolor in between
        i=i+1;fill_01 = i%2; colorfill = 'white' if fill_01 else 'lightgray'; 
        # format the text block
        t = gene_labels_pd.loc[gene_labels_pd['mec_i']==int(mec_i), 'j']
        mec_block_key = ( t.min(), t.max()+1 )
        # [TODO] printing more genes with new line
        #gs_block = gene_labels_pd.loc[gs].loc[gene_labels_pd['mec_i'] == int(mec_i)].index.tolist()
        gs_block = gene_labels_pd.loc[gene_labels_pd['mec_i'] == int(mec_i)].index.tolist()[: Ntopgenes] #2: N genes marked
        #mec_block_text = ' ' + ', '.join(gs_block)  # need the first space
        my_str = ' ' + ', '.join(gs_block)
        mec_block_text = '\n'.join(my_str[i:i+ genes_newline_str_len ] for i in range(0, len(my_str), genes_newline_str_len )) #[:15] # limit character per line
        
        mec_block_nametext = ' ' + str( mec_name.loc[ int(mec_i), 'Annotation' ]  )
        my_str = mec_block_nametext
        mec_block_nametext = '\n'.join(my_str[i:i+ annotext_newline_str_len ] for i in range(0, len(my_str), annotext_newline_str_len ))
        mec_block.update( {
            mec_block_key : 
            { 'color' : colorfill,
             'text': mec_block_text,
             'textname': mec_block_nametext,
            }})
    return(mec_block)


######### Plotting Meta-component heatmap ########

def complex_clustermap(
    dat, 
    row_colors, 
    col_colors, 
    hlines_loc, 
    vlines_loc, 
    figsize = (16,16),
    annotation_texts=None, 
    annotation_texts_block=None,
    annotation_texts_block1_fontsize = 8,
    annotation_texts_block1_fontweight= 'regular',
    annotation_texts_block2_fontsize = 10,
    annotation_texts_block2_fontweight= 'bold',
    annotation_texts_block2_fontcolor = 'black',
    block2_text_yoffset = 0.15,
    fillcoloralpha=1,
    row_cluster_fillcolor_df=None,
    width0=1,
    width1=5,
    width2=6,
    colors_ratio = (0.015, 0.015), # 
    annoname_rotation = 0,
    annotation_panel_ratio ="35%"  ,
    ):
    """Drawing the complex cluster map
    
    Args:
        dat: input data matrix (n_row by m_col) to draw
        row_colors: row color definition. row_colors[i] is corresponding to row i.
        col_colors: col color definition. col_colors[j] is corresponding to col j.
        hlines_loc: positions of horizontal lines to draw. The hlines_loc[p] indicates 
            the p'th horizontal line's y coordinate (y coordinate = row index).
        vlines_loc: positions of vertical lines to draw. The vlines_loc[q] indicates
            the q'th vertical line's x coordinate (x coordinate = col index)
        annotation_texts: a dictionary that stores text annotations.
            e.g {"a": 0, "b":14} means row 0's annotation text is "a".
            and row 14's annotation text is "b".
        annotation_texts_block: a dictionary that stores block text annotations. Key is a
            two_tuple (y1, y2) denoting the start and ending position of the data matrix's 
            row index. value is a dictionary: 
            {"text": list of strings to annotate, "color": the color of that block}.
        width0, width1,width2: width for each annotation block. need to tue with figsize.
    
    Note:
        The coordination system for drawing:
            (0,0) is at the uppler left corner.
            x coordinate corresponding to col index
            y coordinate corresponding to row index.
        Partly borrowed from 
            https://stackoverflow.com/questions/49582178/add-1-histogram-to-side-of-clustermap
    Example:


        # new annotation text style
        gs_block = {
            (2,4): {
                "color":'white',
                "text": " gje134\n rib098\n lti123"
            },
            (4,9): {
                "color":'gray',
                "text": " oqr053\n kvw102\n boq094\n oqr053\n kvw102\n boq094\n oqr053\n kvw102\n boq094\n oqr053"
            },
            (9,14): {
                "color":'white',
                "text": " oqr053\n kvw102\n boq094"
                "textname": " MeC_0"
            }
            }

    """

    g = sns.clustermap(data = dat, row_cluster=False, col_cluster=False,
                       dendrogram_ratio = 0.2, cmap = 'RdBu_r',
                          center = 0,
                   xticklabels = False,yticklabels = False,
                   row_colors = row_colors,
                   col_colors = col_colors,
                   colors_ratio = colors_ratio, #https://github.com/mwaskom/seaborn/issues/437
        
                   figsize = figsize )

    # not show colorbar
    g.cax.set_visible(False)
    g.ax_heatmap.set_ylabel('')
    # draw vertical and horizontal lines
    ax = g.ax_heatmap
    ax.hlines(hlines_loc , *ax.get_xlim(), linewidths = 0.5 , colors='gray', linestyle='--')
    ax.vlines(vlines_loc, *ax.get_ylim(), linewidths = 0.5 , colors='gray', linestyle='--')

    # shift dendrogram and col color plots left to make room for annotation plot
    divider_dendrogram = make_axes_locatable(g.ax_col_dendrogram)
    divider_colcolor = make_axes_locatable(g.ax_col_colors)
    ax_dendrogram_shift = divider_dendrogram.new_horizontal(size= annotation_panel_ratio , pad=0)
    ax_colcolor_shift = divider_colcolor.new_horizontal(size= annotation_panel_ratio , pad=0)

    ## create annotation plot axis
    divider_heatmap = make_axes_locatable(g.ax_heatmap)
    ax_annotate = divider_heatmap.append_axes("right", size= annotation_panel_ratio , pad=0)

    if annotation_texts and annotation_texts_block:
        raise ValueError("Two annotation methods are provided. You should only provide one.")

    # Draw annotation plot, simply pointing to some genes
    if annotation_texts:
        ## detemine text and annotatoin line locations
        texts = sorted(list(annotation_texts.items()), key=lambda x: x[1]) # convert dict to list, order by their row_name_index
        text_left = [i[1]+0.5 for i in texts]
        text_right = list(np.arange(0, len(dat), len(dat) / (len(texts)-1)))
        text_right.append(len(dat)-1)
        ## drawing
        for i, t in enumerate(texts):
            ax_annotate.plot([0,1],[text_left[i],text_right[i]], color='gray', alpha=0.5, linewidth = 0.5)
            ax_annotate.text(1, text_right[i]+1,t[0], fontsize= annotation_texts_block1_fontsize )
        ax_annotate.set_ylim([dat.shape[0],0]) # make sure the y coordinate aligns with heatmap's y coordinate
    
    # Draw annotation text block plot
    if annotation_texts_block and row_cluster_fillcolor_df is None:
        # sort annotation text block according to their y1
        box_height = len(dat) / len(annotation_texts_block)
        box_y_boundary = [i*box_height for i in range(len(annotation_texts_block)+1)]
        for i, block in enumerate(sorted(annotation_texts_block, key=lambda x:x[0])):
            # 0-1: arrow space. 1-5: gene block. 5-9: mec name block
            line1_left = block[0]
            line2_left = block[1]
            line1_right = box_y_boundary[i]
            line2_right = box_y_boundary[i+1]
            ax_annotate.fill_between(
                [0,width0 ],  
                [line1_left, line1_right],
                [line2_left, line2_right],
                #facecolor = annotation_texts_block[block]['color'],
                #edgecolor = 'gray'
                facecolor = row_colors[i],
                edgecolor = 'white'
                
            )
            ax_annotate.fill_between(
                [width0,width0+width1 ],
                [line1_right, line1_right],
                [line2_right, line2_right],
                edgecolor = "gray",
                facecolor = 'white'
            )
            ax_annotate.fill_between(
                [width0+width1, width0+width1+width2],
                [line1_right, line1_right],
                [line2_right, line2_right],
                edgecolor = "gray",
                facecolor = 'white'
            )
            annotation_n_line = len(
                annotation_texts_block[block]['text'].split('\n'))
            ax_annotate.text(
                width0+0.15, 
                #line1_right + figsize[1]*0.2*annotation_n_line  ,
                line2_right - 1 ,
                annotation_texts_block[block]['text'],
                fontsize = annotation_texts_block1_fontsize,
                 weight= annotation_texts_block1_fontweight,
                )
            ax_annotate.text(
                width0+width1+0.15, 
                #line1_right + figsize[1]*0.2*annotation_n_line  ,
                line2_right - 1 ,
                annotation_texts_block[block]['textname'],
                fontsize = annotation_texts_block2_fontsize,
                 weight= annotation_texts_block2_fontweight,
                )
        ax_annotate.set_ylim([dat.shape[0],0]) # make sure the y coordinate aligns with heatmap's y coordinate

    # color fill text block
    if annotation_texts_block and row_cluster_fillcolor_df is not None:
        # row_cluster_fillcolor_df need to have same order with annotation_texts_block
        facecolorcol = row_cluster_fillcolor_df.columns[0]

        # sort annotation text block according to their y1
        box_height = len(dat) / len(annotation_texts_block)
        box_y_boundary = [i*box_height for i in range(len(annotation_texts_block)+1)]
        
        for i, block in enumerate(sorted(annotation_texts_block, key=lambda x:x[0])):
            # 0-1: arrow space. 1-5: gene block. 5-9: mec name block
            line1_left = block[0]
            line2_left = block[1]
            line1_right = box_y_boundary[i]
            line2_right = box_y_boundary[i+1]
            ax_annotate.fill_between(
                [0,width0 ],  
                [line1_left, line1_right],
                [line2_left, line2_right],
                facecolor = row_cluster_fillcolor_df.iloc[i][facecolorcol ],
                edgecolor = 'gray',
                alpha=fillcoloralpha,
                
            )
            ax_annotate.fill_between(
                [width0,width0+width1 ],
                [line1_right, line1_right],
                [line2_right, line2_right],
                edgecolor = "gray",
                facecolor = 'white'
            )
            ax_annotate.fill_between(
                [width0+width1, width0+width1+width2],
                [line1_right, line1_right],
                [line2_right, line2_right],
                facecolor = row_cluster_fillcolor_df.iloc[i][facecolorcol ],
                edgecolor = 'gray',
                alpha=fillcoloralpha,
            )
            annotation_n_line = len(
                annotation_texts_block[block]['text'].split('\n'))
            ax_annotate.text(
                width0+0.15, 
                #line1_right + figsize[1]*0.2*annotation_n_line  ,
                line2_right - 1 ,
                annotation_texts_block[block]['text'],
                fontsize = annotation_texts_block1_fontsize,
                weight = annotation_texts_block1_fontweight
                )
            ax_annotate.text(
                width0+width1+0.15, 
                #line1_right + figsize[1]*0.2*annotation_n_line  ,
                line2_right - 1 - block2_text_yoffset ,
                annotation_texts_block[block]['textname'],
                fontsize = annotation_texts_block2_fontsize,
                color= annotation_texts_block2_fontcolor,
                weight= annotation_texts_block2_fontweight,
                rotation=annoname_rotation,
                )
        ax_annotate.set_ylim([dat.shape[0],0]) # make sure the y coordinate aligns with heatmap's y coordinate

    ax_annotate.set_axis_off()
    return g




################## Plot mecname and top genes only. ###################

def annotation_plot( mectable ,mec_order ,mec_topg, mec_gene_block,
               width0 = 1.5,
                width1 = 1.5,
                height = 15,
                block1_fontsize = 5,
                block2_fontsize = 4,
                block1_weight='bold',
                Ntopgenes=4,
                annotext_newline_str_len = 30 ):


    fig,ax = plt.subplots(1,figsize=( width0 + width1, height))
    fig.gca().invert_yaxis()
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    for i, ( meci, block)  in enumerate( zip( mec_order,  mec_gene_block ) ):
        block1_text = mec_gene_block[block]['textname']
        block1_color = mectable[mectable['MeCi']== meci ]['lineagecolor'].values[0]
        #block2_text = mec_gene_block[block]['text']
        block2_genes =  mec_topg['TopGene20'].loc[ meci ].split(',')[:Ntopgenes] 
        my_str = ', '.join( block2_genes )
        block2_text = '\n'.join(my_str[i:i+ annotext_newline_str_len ] for i in range(0, len(my_str), annotext_newline_str_len ))

        ax.fill_between(
                [0,width0], [i, i], [i+1,i+1], 
            facecolor = block1_color ,
            edgecolor='gray',
            alpha=1
        )
        ax.fill_between(
                [width0,width0+width1], [i, i], [i+1,i+1], 
            facecolor = 'white',
            edgecolor='gray',
            alpha=1
        )
        ax.text(
                    0.1, 
                    i+1-0.2,
                    s = block1_text,
                    fontsize = block1_fontsize,
                    weight= block1_weight,
        )
        ax.text(
                    0.05 + width0, 
                    i+1-0.1,
                    s = block2_text,
                    fontsize = block2_fontsize,
                    weight= 'regular',
        )
    _=ax.set_yticks([])
    _=ax.set_xticks([])

    return( fig )
