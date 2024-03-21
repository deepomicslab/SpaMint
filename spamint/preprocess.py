import pickle
import numpy as np
import pandas as pd
import anndata
# import time
# import logging
import os


def read_csv_tsv(filename):
    if ('csv' in filename) or ('.log' in filename):
        tmp = pd.read_csv(filename, sep = ',',header = 0,index_col=0)
    else:
        tmp = pd.read_csv(filename, sep = '\t',header = 0,index_col=0)
    return tmp



def load_lr_df(species = 'Human',lr_dir = None):
    if lr_dir:
        lr_df = read_csv_tsv(lr_dir)
    else:
        lr_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__))) + '/LR/'
        # print(lr_path)
        if species in ['Human','Mouse']:
            # to lowercase
            species = species.lower()
            lr_df = pd.read_csv(f'{lr_path}/{species}_LR_pairs.txt',sep='\t',header=None)
        else:
            raise ValueError(f'Currently only support Human and Mouse, get {species}')
    return lr_df


def make_adata(mat,meta,species,save_path = None, save_adata = True ):
    # mat: exp matrix, should be cells x genes
    # index should be strictly set as strings
    meta.index = meta.index.map(str)
    mat.index = mat.index.map(str)
    mat = mat.loc[meta.index]
    adata = anndata.AnnData(mat,dtype=np.float32)
    adata.obs = meta
    adata.var = pd.DataFrame(mat.columns.tolist(), columns=['symbol'])
    adata.var_names = adata.var['symbol'].copy()
    #sc.pp.filter_cells(adata, min_genes=200)
    #sc.pp.filter_genes(adata, min_cells=3)
    # remove MT genes for spatial mapping (keeping their counts in the object)
    if species == 'Mouse':
        adata.var['MT_gene'] = [gene.startswith('mt-') for gene in adata.var['symbol']]
    if species == 'Human':
        adata.var['MT_gene'] = [gene.startswith('MT-') for gene in adata.var['symbol']]
    adata.obsm['MT'] = adata[:, adata.var['MT_gene'].values].X.toarray()
    adata = adata[:, ~adata.var['MT_gene'].values]
    if save_path:
        if not os.path.exists(save_path):
            os.makedirs(save_path)
    else:
        save_path = os.getcwd()
    
    figpath = save_path + '/figures/'
    
    if not os.path.exists(figpath):
        os.makedirs(figpath)

    if save_adata:
        adata.write(f'{save_path}/adata.h5ad')
    adata.uns['save_path'] = save_path
    adata.uns['species'] = species
    adata.uns['figpath'] = figpath
    return adata



def lr2kegg(lri_df, use_lig_gene = True, use_rec_gene = True):
    '''
    Use both ligand and receptor
    '''
    if use_lig_gene:
        a = lri_df[['ligand','lr_co_exp_num','lr_co_ratio_pvalue']]
        a.columns = ['gene','lr_co_exp_num','lr_co_ratio_pvalue']
    else:
        a = pd.DataFrame(columns = ['gene','lr_co_exp_num','lr_co_ratio_pvalue'])

    if use_rec_gene:
        b = lri_df[['receptor','lr_co_exp_num','lr_co_ratio_pvalue']]
        b.columns = ['gene','lr_co_exp_num','lr_co_ratio_pvalue']
    else:
        b = pd.DataFrame(columns = ['gene','lr_co_exp_num','lr_co_ratio_pvalue'])
    c = pd.concat((a,b))
    c = c.groupby('gene').mean().reset_index()
    return c


def filter_kegg(df, pval_thred = 0.05):
    tmp = df.copy()
    tmp = tmp[tmp['pvalue'] < pval_thred].copy()
    tmp['-log10 pvalue'] = np.log10(tmp['pvalue']) * (-1)
    tmp[['tmp1','tmp2']] = tmp['GeneRatio'].str.split('/',expand=True)
    tmp['GeneRatio'] = tmp['tmp1'].astype(int) / tmp['tmp2'].astype(int)
    tmp['Count'] = tmp['Count'].astype(int)
    tmp['-log10 pvalue'] = tmp['-log10 pvalue'].astype(float)
    tmp = tmp.sort_values('GeneRatio',ascending=False)
    # remove suffix
    tmp['Description'] = tmp['Description'].str.split(' - ', expand=True)[0]
    return tmp