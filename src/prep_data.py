import pandas as pd
import numpy as np
# from scipy.sparse import csr_matrix
import scanpy as sc
import anndata


def scale_sum(x,SUM):
    res = x.divide(x.sum(axis = 1),axis=0)
    return res*SUM

def prep_adata(mat,meta,species):
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
    return adata

def data_clean(sc_exp, st_exp):
    # cell x genes
    # 1. remove unexpressed genes
    filtered_sc = sc_exp.loc[:,(sc_exp != 0).any(axis=0)]
    filtered_st = st_exp.loc[:,(st_exp != 0).any(axis=0)]
    st_gene = set(filtered_st.columns)
    sc_gene = set(filtered_sc.columns)
    shared_genes = list(st_gene.intersection(sc_gene))
    filtered_sc1 = filtered_sc.loc[:,shared_genes]
    filtered_st1 = filtered_st.loc[:,shared_genes]
    return filtered_sc1, filtered_st1 

def denoise_genes(sc_exp, st_exp, sc_distribution,species):
    sc_genes = sc_exp.columns.tolist()
    st_genes = st_exp.columns.tolist()
    genes = list(set(sc_genes).intersection(set(st_genes)))
    genes = list(set(genes).intersection(set(sc_distribution.columns)))

    if species == 'mouse':
        mt = [gene for gene in genes if gene.startswith('mt-')]
    if species == 'human':
        mt = [gene for gene in genes if gene.startswith('MT-')]
    genes = list(set(genes).difference(set(mt)))
    genes.sort()
    return genes


def prep_all_adata(sc_exp = None, st_exp = None, sc_distribution = None, 
                   sc_meta = None, st_coord = None, lr_df = None, SP = 'human'):
    '''
    1. remove unexpressed genes
    2. select shared genes
    3. transform to adata format
    '''
    # scale all genes to [0,10]
    # v5 
    # SUM = st_exp.sum(axis = 1).mean()
    # v6 from st sum to 1e4
    if (SP != 'human') and (SP != 'mouse'):
        raise ValueError(
            f'Species should be choose among either human or mouse.')
    SUM = 1e4
    # Data Clean
    sc_exp, st_exp = data_clean(sc_exp, st_exp)
    genes = denoise_genes(sc_exp, st_exp, sc_distribution, SP)
    sc_exp = sc_exp[genes]
    st_exp = st_exp[genes]
    sc_distribution = sc_distribution[genes]
    lr_df = lr_df[lr_df[0].isin(genes) & lr_df[1].isin(genes)]
    # Adata Preparation
    # 1. SC to adata
    scale_sc_exp = scale_sum(sc_exp,SUM)
    sc_adata = prep_adata(scale_sc_exp,sc_meta,SP)
    # 2. ST to adata
    scale_st_exp = scale_sum(st_exp,SUM)
    st_adata = prep_adata(scale_st_exp,st_coord,SP)
    # 3. distribution to adata
    sc_ref = scale_sum(sc_distribution,SUM)
    # v6 canceled ref adata
    # sc_ref = prep_adata(scale_poisson_spot,sc_ref_meta,SP)
    if sc_adata.shape[1] == st_adata.shape[1] and st_adata.shape[1] == sc_ref.shape[1]:
        print(f'Data clean is done! Using {st_adata.shape[1]} shared genes .')
    return sc_adata, st_adata, sc_ref, lr_df

