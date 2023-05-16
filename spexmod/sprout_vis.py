import anndata
import os
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
import scanpy as sc

import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from matplotlib_venn import venn3
import plotly.express as px
import math
from . import utils


def creat_vis_dir(save_path):
    vis_path = f'{save_path}/vis/'
    if not os.path.exists(vis_path):
        os.makedirs(vis_path)
    print(f'Saving plots in {vis_path}...')

def findMarkerGenes(adata_in, tp_columns):
    adata = adata_in.copy()
    # 1.prep
    sc.pp.filter_genes(adata, min_cells=30)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    # 2.test
    sc.tl.rank_genes_groups(adata, tp_columns, method='t-test')
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
    marker_genes_df = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(20)
    marker_genes_df
    return adata,marker_genes_df

def adata_umap(adata):
    sc.pp.scale(adata, max_value=10000)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=5, n_pcs=20)
    sc.tl.umap(adata)
    sc.tl.leiden(adata)

def loss_plot(save_path):
    vis_path = f'{save_path}/vis/'
    # read file
    loss_df = pd.read_csv(f'{save_path}/loss.csv',sep = ',',header=0,index_col=0)
    # print loss values    
    n_epoch =  loss_df.shape[0]    
    start_loss = loss_df.loc[0,'total']
    end_loss = loss_df.loc[n_epoch - 1,'total']
    #de_level = (start_loss - end_loss)/start_loss
    print(f'Loss drop from {start_loss} to {end_loss} in {n_epoch} epochs.')
    loss_df.plot()
    plt.savefig(f'{vis_path}/loss.pdf')

def alter_umap_plot(alter_adata,cell_tp,save_path):
    vis_path = f'{save_path}/vis/'
    adata_umap(alter_adata)
    sc.pl.umap(alter_adata, color=[cell_tp],show=False)
    plt.savefig(f'{vis_path}/alter_umap.pdf')


def venn_plot3(save_path,svg,hvg_sc,hvg_ref):
    vis_path = f'{save_path}/vis/'
    plt.figure(figsize=(6,6))
    venn3([svg, hvg_sc, hvg_ref], ('ST','SC','SC ref'))
    plt.savefig(f'{vis_path}/venn3.pdf')


def venn_plot(save_path,marker_genes_df,alter_marker_genes_df):
    vis_path = f'{save_path}/vis/'
    # 1.paras
    num_inter_df = pd.DataFrame()
    tp_num = len(marker_genes_df.columns)
    COL_NUM = 5
    ROW_NUM = int(np.ceil(tp_num/COL_NUM))
    i=0
    # 2.venn plot
    plt.figure(figsize=(3*COL_NUM,3*ROW_NUM))
    for tp in marker_genes_df.columns:
        plt.subplot(ROW_NUM, COL_NUM, i + 1)
        set2 = set(alter_marker_genes_df[tp].tolist())
        set1 = set(marker_genes_df[tp].tolist())
        inter_genes = set1.intersection(set2)
        num_inter_df.loc[tp,'num'] = len(inter_genes)/20
        plt.title(tp)
        venn2([set1, set2], ('Original', 'Imputed'))
        #plt.show()
        i+=1
        #break
    print(f'{num_inter_df["num"].mean():.0%} top 20 marker genes are preserved.')
    plt.savefig(f'{vis_path}/marker_venn.pdf')
    ## 3.boxplot
    fig = px.box(num_inter_df, y = 'num', color_discrete_sequence=['#FB8E68','#E88AC2','#63C1A3', '#9DB0D4'],points = False)
    fig.update_layout(
        autosize=False,
        width=250,
        height=400,
        margin=dict(l=50,r=50,b=100,t=100,pad=4),
        plot_bgcolor="White", paper_bgcolor="White",
        yaxis=dict(title_text="Num of kept top 20 marker genes",titlefont=dict(size=12),)
    )
    fig.update_xaxes(showline=True, linewidth=1.3, linecolor='#313F5F')
    fig.update_yaxes(showline=True, linewidth=1.3, linecolor='#313F5F')
    fig.update_layout(showlegend=False)
    #fig.show()
    fig.write_image(f'{vis_path}/marker_boxplot.pdf') 

def sc_coord_plot(save_path,meta,tp):
    vis_path = f'{save_path}/vis/'
    sc_coord = pd.read_csv(f'{save_path}/coord_best_in_shape.csv',sep = ',',header=None,index_col=None)
    coord = pd.DataFrame(sc_coord.values,columns=['UMAP1','UMAP2'])
    coord[tp] = meta[tp].values
    fig = px.scatter(coord, x='UMAP1', y='UMAP2',color=tp,
                    opacity = 1,
                    width=600,height=500,
                    template = "simple_white",
                    title = 'SPROUT SC raw')
    fig.update_traces(marker_size=4)
    fig.update_layout(coloraxis_showscale=False)
    fig.show()
    fig.write_image(f'{vis_path}/sc_coord.pdf') 

def rotate_via_numpy(xy, radians):
    """Use numpy to build a rotation matrix and take the dot product."""
    x, y = xy
    c, s = np.cos(radians), np.sin(radians)
    j = np.matrix([[c, s], [-s, c]])
    m = np.dot(j, [x, y])

    return float(m.T[0]), float(m.T[1])


def rotate_sc_coord(save_path,meta,tp,angle):
    '''This function rotates the coordinates clockwise'''
    vis_path = f'{save_path}/vis/'
    coord = pd.read_csv(f'{save_path}/coord_best_in_shape.csv',sep = ',',header=None,index_col=None)
    new_coord = pd.DataFrame()
    for po in np.array(coord)[:,0:2]:
        theta = math.radians(angle)
        #theta = math.radians(90)
        tmp = rotate_via_numpy(tuple(po), theta)
        new_coord= new_coord.append(pd.DataFrame(tmp).T)
    new_coord.columns = ['UMAP1','UMAP2']
    new_coord[tp] = meta[tp].values
    fig = px.scatter(new_coord, x='UMAP1', y='UMAP2',color=tp,
                    opacity = 1,
                    width=600,height=500,
                    template = "simple_white",
                    title = 'SPROUT SC rotated')
    fig.update_traces(marker_size=4)
    fig.update_layout(coloraxis_showscale=False)
    fig.show()
    fig.write_image(f'{vis_path}/sc_coord_rotate.pdf')
    del new_coord[tp]
    return new_coord

def rotate_st_coord(coord,angle):
    '''This function rotates the coordinates clockwise'''
    new_coord = pd.DataFrame()
    for po in np.array(coord)[:,0:2]:
        theta = math.radians(angle)
        #theta = math.radians(90)
        tmp = rotate_via_numpy(tuple(po), theta)
        new_coord= new_coord.append(pd.DataFrame(tmp).T)
    new_coord.columns = ['x','y']
    fig = px.scatter(new_coord, x='x', y='y',
                    opacity = 1,
                    width=600,height=500,
                    template = "simple_white",
                    title = 'ST rotated')
    fig.update_traces(marker_size=4)
    fig.update_layout(coloraxis_showscale=False)
    fig.show()
    return new_coord

def get_mutual_marker(marker_genes_df,alter_marker_genes_df):
    tp_inter_genes = {}
    for tp in marker_genes_df.columns:
        set2 = set(alter_marker_genes_df[tp].tolist())
        set1 = set(marker_genes_df[tp].tolist())
        inter_genes = set1.intersection(set2)
        tp_inter_genes[tp] = inter_genes
    return tp_inter_genes

def marker_plot(save_path, st_adata, sc_adata, alter_adata, tp_inter_genes):
    sc_agg_adata = sc_adata[alter_adata.obs['sc_id']]
    sc_coord = alter_adata.obs[['UMAP1','UMAP2']]
    st_coord = st_adata.obs[['x','y']]
    vis_path = f'{save_path}/vis/'
    cmap = 'OrRd'
    cmap = 'viridis'
    #cmap = 'plasma'
    #cmap = 'RdPu'
    MIN = 0
    MAX = 6
    for key, values in tp_inter_genes.items():
        print(key)
        COL_NUM = 5
        ROW_NUM = len(values)
        i=0
        a_tmp = alter_adata[:,list(values)].to_df()
        b_tmp = sc_agg_adata[:,list(values)].to_df()
        c_tmp = st_adata[:,list(values)].to_df()

        a_tmp['spot'] = alter_adata.obs.spot
        d_tmp = a_tmp.groupby('spot').sum()
        d_coord = st_coord.loc[d_tmp.index]

        b_tmp['spot'] = alter_adata.obs.spot.values
        e_tmp = b_tmp.groupby('spot').sum()
        e_coord = st_coord.loc[e_tmp.index]
        plt.figure(figsize=(4*COL_NUM,4*ROW_NUM))
        for gene in values:
            plt.subplot(ROW_NUM, COL_NUM, i + 1)
            t = np.array(utils.scale_01(np.log(c_tmp[gene]+1),MIN,MAX))
            #t = np.array(np.log(c_tmp[gene]+1))
            plt.scatter(st_coord['x'], st_coord['y'],c = t, s = 20,cmap=cmap)
            plt.title(f'{gene} ST', fontsize=12)
            plt.xticks([])
            plt.yticks([])

            plt.subplot(ROW_NUM, COL_NUM, i + 4)
            t = np.array(utils.scale_01(np.log(b_tmp[gene]+1),MIN,MAX))
            #t = np.array(np.log(b_tmp[gene]+1))
            plt.scatter(sc_coord['UMAP1'], sc_coord['UMAP2'],c = t, s = 5,cmap=cmap)
            plt.title(f'{gene} Original SC', fontsize=12)
            plt.xticks([])
            plt.yticks([])

            plt.subplot(ROW_NUM, COL_NUM, i + 5)
            t = np.array(utils.scale_01(np.log(a_tmp[gene]+1),MIN,MAX))
            #t = np.array(np.log(a_tmp[gene]+1))
            plt.scatter(sc_coord['UMAP1'], sc_coord['UMAP2'],c = t, s = 5,cmap=cmap)
            plt.title(f'{gene} Altered SC', fontsize=12)
            plt.xticks([])
            plt.yticks([])

            plt.subplot(ROW_NUM, COL_NUM, i + 2)
            t = np.array(utils.scale_01(np.log(d_tmp[gene]+1),MIN,MAX))
            #t = np.array(np.log(d_tmp[gene]+1))
            plt.scatter(d_coord['x'], d_coord['y'],c = t, s = 20,cmap=cmap)
            plt.title(f'{gene} Altered summed ST', fontsize=12)
            plt.xticks([])
            plt.yticks([])

            plt.subplot(ROW_NUM, COL_NUM, i + 3)
            t = np.array(utils.scale_01(np.log(e_tmp[gene]+1),MIN,MAX))
            #t = np.array(np.log(e_tmp[gene]+1))
            plt.scatter(e_coord['x'], e_coord['y'],c = t, s = 20,cmap=cmap)
            plt.title(f'{gene} Original summed ST', fontsize=12)
            plt.xticks([])
            plt.yticks([])
            i+=5
            #break
            plt.savefig(f'{vis_path}/maker_exp_{key}.pdf')
        # break
    #


