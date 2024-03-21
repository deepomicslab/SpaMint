# import anndata
# import os

# from scipy.sparse import csr_matrix
# import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
# from matplotlib_venn import venn2
# from matplotlib_venn import venn3
# import plotly.express as px
# import math
# from . import utils

def sc_celltype(adata, color_map = None, tp_key = 'celltype', 
            subset_idx = None, legend = False, figsize = (4,4),
            savefig = False, size = 10, alpha = 0.8,title = None,theme = 'white'):
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
    
    if savefig:
        if isinstance(savefig,str):
            if '/' in savefig:
                savefig = f'{savefig}'
            else:
                save_path = adata.uns['figpath']
                savefig = f'{save_path}/{savefig}'
        else:
            save_path = adata.uns['figpath']
            savefig = f'{save_path}/celltype_sc.pdf'
        plt.savefig(f'{savefig}')
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
                   theme_white = True, savefig = False, COL = None, hue = None):
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

    if savefig:
        save_path = adata.uns['figpath']
        savefig = f'{save_path}/{name}_sc_solo.pdf'
        plt.savefig(f'{savefig}')
    plt.show()
    plt.clf()



def boxplot(adata, metric = 'spot_cor', palette_dict = None,
                 x='method', y='pair_num', hue='method',figsize = (2.4, 3),
                 ylabel='LRI count', dodge=False,legend = False,
                 test = 't-test_ind',rotate_x = False,
                 savefig = False, title = None):
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

    if savefig:
        save_path = adata.uns['figpath']
        savefig = f'{save_path}/{metric}_box.pdf'
        plt.savefig(f'{savefig}')


def exp_violin(adata, gene=None, tp_key=None, types=None,
               palette_dict=None, test='t-test_ind',
               highlight=[], log=True, scale=True,
               figsize=(6, 4), x_rotate=False,
               savefig = False, title=''):
    from statannotations.Annotator import Annotator
    from itertools import combinations
    df = adata[:, gene].to_df()

    if log:
        df = np.log(df + 1)

    if scale:
        df = (df - df.mean()) / df.std()

    df = pd.concat((df, adata.obs[tp_key]), axis=1)
    df = df[df[tp_key].isin(types)]

    plt.figure(figsize=figsize)
    ax = sns.violinplot(data=df, x=tp_key, y=gene,
                        order=types, palette=palette_dict,
                        linewidth=.9, cut = 0, scale = 'width')
    if highlight:
        pairs = [(cell_type1, cell_type2) for cell_type1 in highlight for cell_type2 in types if cell_type1 != cell_type2]
    else:
        pairs = list(combinations(types, 2))

    annot = Annotator(ax, pairs, x=tp_key, y=gene, data=df)
    annot.configure(test=test, comparisons_correction="BH", correction_format="replace")
    annot.apply_test()
    annot.annotate()

    plt.tight_layout()
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel('')
    plt.ylabel(gene, size=21)
    plt.title(title, size=21)

    if x_rotate:
        plt.xticks(rotation=90)

    if savefig:
        save_path = adata.uns['figpath']
        savefig = f'{save_path}/{title}_exp_violin.pdf'
        plt.savefig(f'{savefig}')



def highlight_y(draw_df, y, y_highlight):
    arr = draw_df[y].unique()
    highlighted_ytick = np.where(np.isin(arr, y_highlight))[0]
    # print(highlighted_ytick) # The ytick to be highlighted
    # Customize the ytick color
    yticks, _ = plt.yticks()
    ytick_labels = plt.gca().get_yticklabels()
    for ytick, ytick_label in zip(yticks, ytick_labels):
        # print(ytick)
        if ytick in highlighted_ytick:
            ytick_label.set_color('red')


def highlight_x(draw_df, x, x_highlight):
    arr = draw_df[x].unique()
    highlighted_xtick = np.where(np.isin(arr, x_highlight))[0]
    # print(highlighted_xtick) # The ytick to be highlighted
    # Customize the ytick color
    xticks, _ = plt.xticks()
    xtick_labels = plt.gca().get_xticklabels()
    for xtick, xtick_label in zip(xticks, xtick_labels):
        # print(ytick)
        if xtick in highlighted_xtick:
            xtick_label.set_color('red')


def draw_bubble(draw_df, x=None, y="Description", x_highlight = None, y_highlight=None, 
                savefig=None, figsize=None, legend=False, xrotate=False,
                showlabel=False, xlabel = None, ylabel = None, title=None, 
                  cmap='viridis',hue='Count', size='GeneRatio'):
    '''
    Red label item in x, y axis
    highlight = 'T_cells_CD4+'
    or
    highlight = ['T_cells_CD4+','B cell]
    '''
    if figsize is None:
        width = np.max([int(np.ceil(len(draw_df[x].unique()) / 3.5)), 2])+0.2
        height = int(np.ceil(len(draw_df[y].unique())/4)+1)
        if xlabel:
            height += 0.5
        if ylabel:
            width += 0.5
        print(f'Auto figsize {(width,height)}')
    else:
        width = figsize[0]
        height = figsize[1]

    # print(width,height)
    with sns.axes_style("whitegrid"):
        plt.figure(figsize=(width, height))
        sns.scatterplot(data=draw_df, x=x, y=y, hue = hue, 
                        legend=legend,
                        palette=cmap, 
                        size=size, sizes=(50, 300)
                        )
        
        highlight_x(draw_df, x, x_highlight)
        highlight_y(draw_df, y, y_highlight)

        plt.yticks(fontsize=16)
        plt.xticks(fontsize=16)
        if showlabel == False:
            plt.xlabel('')
            plt.ylabel('')
        else:
            if xlabel:
                plt.xlabel(xlabel,fontsize = 22)
            else:
                plt.xlabel('',fontsize = 22)

            if ylabel:
                plt.ylabel(ylabel,fontsize = 22)
            else:
                plt.ylabel('',fontsize = 22)

            if title:
                plt.title(title,fontsize=22)

        if xrotate:
            plt.xticks(rotation=90)

        if legend:
            plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5),fontsize=14)

        plt.xlim(-0.5, len(draw_df[x].unique())-0.5)    
        # if len(draw_df[x].unique()) == 2:
        #     plt.xlim(-0.5, 1.5)
        #     # plt.ylim(10, -0.8)
        # elif len(draw_df[x].unique()) == 3:
        #     plt.xlim(-0.5, 2.5)
        if savefig:
            plt.savefig(f'{savefig}')
    plt.show()
    plt.clf()



def celltype_cci(adata, is_sender_per = True, figsize = (3.6,3.6), y_highlight='', x_highlight='',
                cmap = 'Reds', hue = 'LRI',showlabel = True,
                size = 'CCI',xrotate = True,legend = True,
                savefig = False, title = ''):
    '''
    is_sender_per: if True, use the cci percentage file that sender row sum is 1 [send_per].
    '''
    if is_sender_per:
       df = adata.uns['send_per']
       name = 'sender'
    else:
        df = adata.uns['rec_per']
        name = 'receiver'
    
    if savefig:
        save_path = adata.uns['figpath']
        savefig = f'{save_path}/{name}_CCI.pdf'
    
    draw_bubble(df, x='Receiver', y='Sender', x_highlight = x_highlight, y_highlight=y_highlight,
                    savefig = savefig, title = title, figsize = figsize, legend = legend,
                    showlabel = showlabel, xrotate = xrotate, cmap=cmap, hue=hue, size=size)
    

def clustermap(df, index = 'ligand', col = 'receptor', value = 'lr_co_exp_num',
                 aggfunc = 'sum', log = True, row_cluster = False, col_cluster = False,
                 highlight = None,  cmap = 'coolwarm', title = '',
                 xticks = False, yticks = True, figsize = (5,5),
                 savefig = False):
    '''
    df: dataframe
    if df has there columns and required to be pivot, then use the following parameters
    otherwise set aggfunc as None
        index: row index
        col: column index
        value: value to fill the pivot table
        aggfunc: how to aggregate the value
    log: log the value or not
    row_cluster: cluster row or not
    col_cluster: cluster column or not
    highlight: highlight the rows
    title: title of the plot
    cmap: color map of heatmap
    col_color: color of the column [sorted by the df row's order]
    xticks: show xticks or not
    yticks: show yticks or not
    savefig: output directory
    '''
    if aggfunc is None:
        n_tp_lri = df.copy()
    else:
        if aggfunc == 'sum':
            n_tp_lri = pd.crosstab(index=df[index], columns = df[col],values=df[value],aggfunc = sum)
    n_tp_lri = n_tp_lri.fillna(0)
    if log:
        n_tp_lri = np.log(n_tp_lri + 1)
        legend_label = f'{value} (log)'
    else:
        legend_label = f'{value}'

    clustermap = sns.clustermap(n_tp_lri,row_cluster=row_cluster,col_cluster=col_cluster,
               standard_scale=None,dendrogram_ratio=0.001,square = True,cmap = cmap,
               figsize=figsize, cbar_pos=(1, 0.5, 0.02, 0.2),cbar_kws={'orientation': 'vertical','label':legend_label})
    if not xticks:
        clustermap.ax_heatmap.set_xticklabels([])
        clustermap.ax_heatmap.set_xticks([])
    
    if not yticks:
        clustermap.ax_heatmap.set_yticklabels([])
        clustermap.ax_heatmap.set_yticks([])

    if highlight is not None:
        arr = n_tp_lri.index.to_list()
        # Reorganize the index labels based on the cluster order
        reordered_index = [arr[i] for i in clustermap.dendrogram_row.reordered_ind]
        # Customize the ytick color and font weight
        yticks, _ = plt.yticks()
        ytick_labels = clustermap.ax_heatmap.get_yticklabels()
        for index, ytick_label in enumerate(ytick_labels):
            if reordered_index[index] in highlight:
                ytick_label.set_color('red')
                ytick_label.set_weight('bold')
    
    if title:
        clustermap.ax_heatmap.set_title(title,fontsize=22) 

    if savefig:
        plt.savefig(f'{savefig}')
    

    plt.show()


def lri_heatmap(adata, target_sender = [], target_receiver = [], title = None, 
                figsize = (8,3), log = True, row_cluster = True, col_cluster = False,
                highlight = None, cmap = 'Reds',cellTypeAsCol = True,
                xticks = True, yticks = True,savefig = False):
    
    lri_df = adata.uns['spatalk']
    if len(target_sender) == 0:
        target_sender = lri_df['celltype_sender'].unique()
    
    if len(target_receiver) == 0:
        target_receiver = lri_df['celltype_receiver'].unique()
        
    lri_df = lri_df[lri_df['celltype_sender'].isin(target_sender) & lri_df['celltype_receiver'].isin(target_receiver)].copy()
    
    if savefig:
        save_path = adata.uns['figpath']
        savefig = f'{save_path}/{title}_LRI_co_exp_heatmap.pdf'
    if cellTypeAsCol:
        clustermap(lri_df,index = 'LRI', col = 'celltype_receiver', value = 'lr_co_exp_num',
                        aggfunc = 'sum', log = log, row_cluster = row_cluster, col_cluster = col_cluster,
                        figsize = figsize, cmap = cmap, title = title, savefig = savefig,
                        highlight = highlight, xticks = xticks, yticks = yticks)        
    else:
        clustermap(lri_df,index = 'celltype_receiver', col = 'LRI', value = 'lr_co_exp_num',
                        aggfunc = 'sum', log = log, row_cluster = row_cluster, col_cluster = col_cluster,
                        figsize = figsize, cmap = cmap, title = title, savefig = savefig,
                        highlight = highlight, xticks = xticks, yticks = yticks)
    return lri_df



def top_kegg_enrichment(adata, top_n = 10, groupby_sender = False,
                        target_sender = [], target_receiver = [],target_path = [], 
                        use_lig_gene = True, use_rec_gene = True,
                        cmap = 'RdPu', hue = '-log10 pvalue', size = 'GeneRatio',
                        legend = True, figsize = None, savefig = False):
    # TODO discard cancer related
    kegg_res = adata.uns['kegg_enrichment']

    lig_gene = 'L' if use_lig_gene else 'n'
    rec_gene = 'R' if use_rec_gene else 'n'
    used_genes = f'{lig_gene}_{rec_gene}'

    kegg_res = kegg_res[kegg_res['celltype_sender'].isin(target_sender) & 
                        kegg_res['celltype_receiver'].isin(target_receiver) & 
                        (kegg_res['used_genes'] == used_genes)]
    # print(kegg_res.head(5))

    if groupby_sender:
        show_row = 'celltype_sender'
        xlabel = 'Sender'
        target_row = 'celltype_receiver'
        targets = target_receiver
    else:
        show_row = 'celltype_receiver'
        xlabel = 'Receiver'
        target_row = 'celltype_sender'
        targets = target_sender        

    # TODO subset by sender, plot two with title why plot together
    top_keggs = pd.DataFrame()
    for target in targets:
        title = f'{target}'
        top_kegg = kegg_res[kegg_res[target_row] == target].copy()
        top_kegg = top_kegg.groupby([show_row]).apply(lambda x: x.nlargest(top_n, 'Count')).reset_index(drop=True)
        uniq_path = list(set(top_kegg['Description']))
        # incase some path shows but not in top
        top_kegg = kegg_res[kegg_res['Description'].isin(uniq_path)].copy()
        # top_kegg = top_kegg[~top_kegg['Description'].isin(discard)]
        top_kegg = top_kegg.sort_values(by = [show_row,'Count'],ascending = [True,False])
        top_keggs = pd.concat([top_keggs,top_kegg]) 
        draw_bubble(top_kegg, x = show_row, y_highlight = target_path,
                    y = "Description",cmap = cmap, hue = hue,
                    size = size,xrotate = True,legend = True, 
                    title = title, showlabel = True,
                    xlabel = xlabel,figsize = figsize
                    )
        plt.clf()
    return top_keggs



def keggPath_lri(adata, groupby_sender = False, value = 'lr_co_exp_num', thred = 5,
                target_sender = [], target_receiver = [], target_path = '', 
                use_lig_gene = True, use_rec_gene = True,
                savefig = False, show_title = True,
                figsize = None, log = True, row_cluster = True, col_cluster = False,
                highlight = None, cmap = 'Reds',cellTypeAsCol = True,
                xticks = True, yticks = True):
    kegg_res = adata.uns['kegg_enrichment']
    lri_df = adata.uns['spatalk']

    lig_gene = 'L' if use_lig_gene else 'n'
    rec_gene = 'R' if use_rec_gene else 'n'
    used_genes = f'{lig_gene}_{rec_gene}'

    kegg_res = kegg_res[kegg_res['celltype_sender'].isin(target_sender) & 
                        kegg_res['celltype_receiver'].isin(target_receiver) & 
                        (kegg_res['used_genes'] == used_genes) &
                        (kegg_res['Description'] == target_path)]
    if len(kegg_res) == 0:
        raise ValueError(f'No {target_path} found in the enrichment result.')
    
    lri_df = lri_df[lri_df['celltype_sender'].isin(target_sender) &
                    lri_df['celltype_receiver'].isin(target_receiver)]
    
    path_genes = []
    for _, row in kegg_res.iterrows():
        tmp_genes = row['geneSymbol'].split('/')
        # print(len(tmp_genes))
        path_genes.extend(tmp_genes)
    path_genes = list(set(path_genes))
    if (lig_gene == 'L') and (rec_gene == 'R'):
        target_df = lri_df[lri_df['ligand'].isin(path_genes) | lri_df['receptor'].isin(path_genes)]
    elif (lig_gene == 'L') and (rec_gene == 'n'):
        target_df = lri_df[lri_df['ligand'].isin(path_genes)]
    elif (lig_gene == 'n') and (rec_gene == 'R'):
        target_df = lri_df[lri_df['receptor'].isin(path_genes)]
    else:
        raise ValueError('Both use_lig_gene and use_rec_gene are False. No gene used.')
    
    target_df = target_df[target_df[value]>thred].copy()
    if savefig:
        save_path = adata.uns['figpath']
        savefig = f'{save_path}/LRI_{target_path}_heatmap.pdf'

    if groupby_sender:
        show_row = 'celltype_sender'
    else:
        show_row = 'celltype_receiver'

    if show_title:
        title = target_path

    if cellTypeAsCol:
        clustermap(target_df,index = 'LRI', col = show_row, value = value,
                        aggfunc = 'sum', log = log, row_cluster = row_cluster, col_cluster = col_cluster,
                        figsize = figsize, cmap = cmap, title = title, savefig = savefig,
                        highlight = highlight, xticks = xticks, yticks = yticks)        
    else:
        clustermap(lri_df,index = 'celltype_receiver', col = 'LRI', value = value,
                        aggfunc = 'sum', log = log, row_cluster = row_cluster, col_cluster = col_cluster,
                        figsize = figsize, cmap = cmap, title = title, savefig = savefig,
                        highlight = highlight, xticks = xticks, yticks = yticks)
    return kegg_res



def draw_lr_flow2(df, left_panel = 'ligand', right_panel = 'receptor',
                 figsize = (10,10)):
    import plotly.graph_objects as go
    ligs = list(df[left_panel].unique())
    recs = list(df[right_panel].unique())
    labels = list(df[left_panel].unique())
    label_rec = list(df[right_panel].unique())
    labels.extend(label_rec)
    # define source and target indices for each link
    source = df[left_panel].astype('category').cat.codes.tolist()
    target = df[right_panel].astype('category').cat.codes.tolist()
    target = [x + len(set(source)) for x in target]
    value = [np.random.randint(1, 2) for _ in range((len(df)))]
    trace = go.Sankey(
        node=dict(
            pad=5,
            thickness=20,
            line=dict(color='black', width=0.1),
            label=labels,
            color=['#EFA2B5']*len(ligs) + ['#9DCD82']*len(recs)
        ),
        link=dict(
            source=source, # indices correspond to labels, eg A1, A2, A1, B1 
            target=target, 
            value=value,
        )
    )
    # create layout
    layout = go.Layout(
        title='',
        font=dict(size=18)
    )
    # create figure
    fig = go.Figure(data=[trace], layout=layout)
    width = figsize[0]*100
    height = figsize[1]*100
    fig.update_layout(width=width, height=height)
    # fig.write_image(f'./8.Ref/figures/main/lr_network.pdf') 
    fig.show()
    return


def draw_lr_flow3(df, left_panel = 'ligand', mid_panel = 'receptor', right_panel = 'pathway',
                 figsize = (10,10)):
    import plotly.graph_objects as go
    ligs = list(df[left_panel].unique())
    recs = list(df[mid_panel].unique())
    paths = list(df[right_panel].unique())

    labels = list(df[left_panel].unique())
    label_rec = list(df[mid_panel].unique())
    label_path = list(df[right_panel].unique())
    labels.extend(label_rec)
    labels.extend(label_path)
    print('label',labels)
    label_num_map = dict(zip(labels,range(len(labels))))
    source1 = list(df[left_panel].map(label_num_map))
    target1 = list(df[mid_panel].map(label_num_map))
    target2 = list(df[right_panel].map(label_num_map))
    source = source1 + target1
    target = target1 + target2
    print('source1',source1)
    print('target1',target1)
    print('target2',target2)
    value = [np.random.randint(1, 2) for _ in range((len(source)))]
    trace = go.Sankey(
        node=dict(
            pad=5,
            thickness=20,
            line=dict(color='black', width=0.1),
            label=labels,
            color=["#F56867"]*len(ligs) + ["#3A84E6"]*len(recs) + ["#59BE86"] *len(paths)
        ),
        link=dict(
            source=source, # indices correspond to labels, eg A1, A2, A1, B1 
            target=target, 
            value=value,
        )
    )
    # create layout
    layout = go.Layout(
        title='',
        font=dict(size=18)
    )
    # create figure
    fig = go.Figure(data=[trace], layout=layout)
    width = figsize[0]*100
    height = figsize[1]*100
    fig.update_layout(width=width, height=height)
    # fig.write_image(f'./8.Ref/figures/main/lr_network.pdf') 
    fig.show()
    return