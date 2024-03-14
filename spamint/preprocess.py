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



def load_lr_df(species = 'human',lr_dir = None):
    if lr_dir:
        lr_df = read_csv_tsv(lr_dir)
    else:
        lr_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__))) + '/LR/'
        # print(lr_path)
        if species in ['human','mouse']:
            lr_df = pd.read_csv(f'{lr_path}/{species}_LR_pairs.txt',sep='\t',header=None)
        else:
            raise ValueError(f'Currently only support human and mouse, get {species}')
    return lr_df


def make_adata(mat,meta,species,save_path = None):
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
    if species == 'mouse':
        adata.var['MT_gene'] = [gene.startswith('mt-') for gene in adata.var['symbol']]
    if species == 'human':
        adata.var['MT_gene'] = [gene.startswith('MT-') for gene in adata.var['symbol']]
    adata.obsm['MT'] = adata[:, adata.var['MT_gene'].values].X.toarray()
    adata = adata[:, ~adata.var['MT_gene'].values]
    if save_path:
        adata.write(f'{save_path}/adata.h5ad')
    return adata
