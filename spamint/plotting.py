# import anndata
# import os
# import pandas as pd
# from scipy.sparse import csr_matrix
# import scanpy as sc

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
# from matplotlib_venn import venn2
# from matplotlib_venn import venn3
# import plotly.express as px
# import math
# from . import utils

def sc_celltype(adata, color_map = None, tp_key = 'celltype', 
            subset_idx = None, legend = False, figsize = (4,4),
            save_path = None, size = 10, alpha = 0.8,title = None,theme = 'white'):
    '''
    @ Wang Jingwan 0314
    This function is used to draw scatter plot of single cell data.
    '''
    sc_agg_meta = adata.obs.copy()
    sc_agg_meta['pivot'] = 1
    #sort celltype number from large to small from sc_agg_meta
    sc_agg_meta['celltype_num'] = sc_agg_meta.groupby(tp_key)['pivot'].transform('count')
    #draw scater plot of each celltype in a order of celltype_num, large to small
    sc_agg_meta = sc_agg_meta.sort_values(by=['celltype_num'],ascending=False)
    if 'adj_spex_UMAP1' in sc_agg_meta.columns:
        cols = ['adj_spex_UMAP1','adj_spex_UMAP2']
    elif 'adj_UMAP1' in sc_agg_meta.columns:
        cols = ['adj_UMAP1','adj_UMAP2']
    elif 'st_x' in sc_agg_meta.columns:
        cols = ['st_x','st_y']
    elif 'col' in sc_agg_meta.columns:
        cols = ['row','col']
    else:
        cols = ['x','y']
    plt.figure(figsize=figsize)
    with sns.axes_style("white"):
        if subset_idx is None:
            sns.scatterplot(data=sc_agg_meta, x=cols[0], y=cols[1],hue = tp_key, 
                            s = size,alpha = alpha,palette=color_map,edgecolor = None)
        else:
            subset_data = sc_agg_meta.copy()
            subset_data['subset'] = 'Other cells'
            subset_data.loc[subset_idx,'subset'] = subset_data.loc[subset_idx,tp_key]
            sns.scatterplot(data=subset_data[subset_data['subset'] == 'Other cells'], x=cols[0], y=cols[1],hue = 'subset', 
                            s = 5,alpha = 1,palette=['#cccccc'],edgecolor = None)
            sns.scatterplot(data=subset_data.loc[subset_idx], x=cols[0], y=cols[1],hue = 'subset', 
                            s = size+5,alpha = alpha,palette=color_map,edgecolor = None)
        if legend:
            if sc_agg_meta[tp_key].nunique() > 30:
                ncol = 2
            else:
                ncol = 1
            handles, labels = plt.gca().get_legend_handles_labels()
            # Sort the handles and labels based on the labels
            handles_labels_sorted = sorted(zip(handles, labels), key=lambda x: (x[1] == "Other cells", x[1]))
            handles_sorted, labels_sorted = zip(*handles_labels_sorted)
            # Create a new legend with the sorted handles and labels
            leg = plt.legend(handles_sorted, labels_sorted,loc='center left', bbox_to_anchor=(0.95, 0.5),
                             ncol=ncol, handletextpad=0.5,columnspacing=0.4,labelspacing=0.1,
                             fontsize = 16,markerscale = 2,handlelength = 0.5)
            leg.get_frame().set_linewidth(0.0)  # Remove legend frame
        else:
            plt.legend([],[], frameon=False)  
        plt.title(title,fontsize=21)
        plt.axis('equal')
    # plt.tight_layout()
    if theme == 'white':
        # sns.despine(left=True, bottom=True)
        plt.xlabel('',fontsize=16)
        plt.ylabel('',fontsize=16)
        plt.xticks([],fontsize=14)
        plt.yticks([],fontsize=14)
        # plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9)
    
    if save_path:
        plt.savefig(f'{save_path}')
    plt.show()
    plt.clf()



def decide_figsize(x = None, y = None,unit = 4):
    # small one be the unit
    xrange = x.max() - x.min()
    yrange = y.max() - y.min()
    # col_L is x
    # row_L is y
    if xrange < yrange:
        ROW_L = round(unit * (yrange/xrange))
        COL_L = unit
    else:
        COL_L = round(unit * (xrange/yrange))
        ROW_L = unit
    return ROW_L,COL_L


def sc_subtype(adata,color_map = None,tp_key = 'celltype', target_tp = None, 
                   theme_white = True, save_path = None, COL = None, hue = None):
    sc_agg_meta = adata.obs.copy()
    if hue is None:
        sc_agg_meta['pivot'] = 1
        hue = 'pivot'

    if target_tp is None:
        target_tp = sc_agg_meta[tp_key].unique().tolist()
        name = 'all'
    else:
        name = target_tp[0]
    
    if COL is None:
        # e.g. sqrt(9) = 3x3
        # e.g. sqrt(12) = 3x4 (lesser row more col, thus np.floor)
        ROW = int(np.floor(np.sqrt(len(target_tp))))
        COL = int(np.ceil(len(target_tp) / ROW))
    else:
        ROW = int(np.ceil(len(target_tp) / COL))

    if 'adj_spex_UMAP1' in sc_agg_meta.columns:
        cols = ['adj_spex_UMAP1','adj_spex_UMAP2']
    elif 'adj_UMAP1' in sc_agg_meta.columns:
        cols = ['adj_UMAP1','adj_UMAP2']
    else:
        cols = ['x','y']

    if 'st_x' in sc_agg_meta.columns:
        st_cols = ['st_x','st_y']
    elif 'x' in sc_agg_meta.columns:
        st_cols = ['x','y']
    else:
        st_cols = cols
    i = 0
    ROW_L,COL_L = decide_figsize(x = sc_agg_meta['st_x'], y = sc_agg_meta['st_y'])
    plt.figure(figsize=(COL_L*COL, ROW_L* ROW))
    for tp in target_tp:
        with sns.axes_style("white"):
            plt.subplot(ROW, COL, i + 1)
            sns.scatterplot(data=sc_agg_meta[[st_cols[0],st_cols[1],'pivot']].drop_duplicates(), x=st_cols[0], y=st_cols[1],hue = hue, 
                        s = 20,alpha = 0.8,palette=['#ccc'],edgecolor = None
                        )
            sns.scatterplot(data=sc_agg_meta[sc_agg_meta[tp_key] == tp], x=cols[0], y=cols[1],hue = hue, 
                            s = 20,alpha = 0.8,palette=[color_map[tp]],edgecolor = None
                            )
            if theme_white:
                plt.legend([],[], frameon=False)  
                # sns.despine(left=True, bottom=True)
                plt.tick_params(left=False, bottom=False, top = False, right = False)
                plt.tight_layout()
                plt.xlabel('',fontsize=16)
                plt.ylabel('',fontsize=16)
                plt.xticks([],fontsize=14)
                plt.yticks([],fontsize=14)
            else:
                plt.legend([],[], frameon=False)  
                plt.axis('equal')
                plt.xlabel('x',fontsize=16)
                plt.ylabel('y',fontsize=16)
                plt.subplots_adjust(wspace=0.4,hspace=0.4)
        plt.title(f'{tp}',fontsize=22)
        i+=1
    if save_path:
        plt.savefig(f'{save_path}')
    plt.show()
    plt.clf()



def boxplot(adata, metric = 'spot_cor', palette_dict = None,
                 x='method', y='pair_num', hue='method',figsize = (2.4, 3),
                 ylabel='LRI count', dodge=False,legend = False,
                 test = 't-test_ind',rotate_x = False,
                 save_path = None, title = None):
    '''
    @ Wang Jingwan 0314
    test method can choose from 
    't-test_welch', 't-test_paired', 'Mann-Whitney', 'Mann-Whitney-gt', 
    'Mann-Whitney-ls', 'Levene', 'Wilcoxon', 'Kruskal', 'Brunner-Munzel'
    '''
    if metric == 'spot_cor':
        draw_df = adata.uns['spot_cor']

    from statannotations.Annotator import Annotator

    plt.figure(figsize=figsize)
    ax = sns.boxplot(x=x, y=y, data=draw_df, hue=hue, 
                     dodge=dodge, palette=palette_dict,
                     width=.8)
    pairs = [tuple(draw_df['map'].unique())]
    annot = Annotator(ax, pairs, x=x, y=y, data=draw_df)
    annot.configure(test=test, comparisons_correction="BH",correction_format="replace")
    annot.apply_test()
    annot.annotate()
    plt.tight_layout()
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel('', size=16)
    plt.ylabel(ylabel, size=16)
    plt.title(title,size = 16)
    if rotate_x:
        plt.xticks(rotation=90)

    if legend:
        # plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
        leg = ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.16), ncol=2,
                        handletextpad=0.3,columnspacing=0.3,fontsize = 14)
        leg.get_frame().set_linewidth(0.0)  # Remove legend frame
    else:
        plt.legend([],[], frameon=False)
    if save_path:
        plt.savefig(f'{save_path}')
