import pickle
import numpy as np
import pandas as pd

def scale_01(x,a,b):
    MAX = np.max(x)
    MIN = np.min(x)
    res = (b-a)*(x-MIN)/(MAX-MIN)+a
    return res

def scale_global_MIN_MAX(df,MIN,MAX):
    df_max = df.max().max()
    df_min = df.min().min()
    df_01 = (df - df_min)/(df_max - df_min)
    res = df_01*(MAX - MIN) - MIN
    return res
    
def save_object(obj, filename):
    with open(filename, 'wb') as outp:  # Overwrites any existing file.
        pickle.dump(obj, outp, pickle.HIGHEST_PROTOCOL)

def load_object(filename):
    pkl_file = open(filename, 'rb')
    data = pickle.load(pkl_file)
    return data

def check_positive(**params):
    """Check that parameters are positive as expected
    Raises
    ------
    ValueError : unacceptable choice of parameters
    """
    for p in params:
        if params[p] <= 0:
            raise ValueError(
                "Expected {} > 0, got {}".format(p, params[p]))


def check_int(**params):
    """Check that parameters are integers as expected
    Raises
    ------
    ValueError : unacceptable choice of parameters
    """
    for p in params:
        if not isinstance(params[p], int):
            raise ValueError(
                "Expected {} integer, got {}".format(p, params[p]))

def check_num(**params):
    """Check that parameters are numeric as expected
    Raises
    ------
    ValueError : unacceptable choice of parameters
    """
    for p in params:
        if not params[p].isnumeric():
            raise ValueError(
                "Expected {} numeric, got {}".format(p, params[p]))

def check_st_coord(st_adata):
    st_meta = st_adata.obs
    if 'x' and 'y' in st_meta.columns:
        st_coord = st_meta[['x','y']]
    elif 'row' and 'col' in st_meta.columns:
        st_coord = st_meta[['row','col']]
    else:
        raise ValueError(
            f'st_adata expected two columns either name x y or row col to represent spatial coordinates, but got None in {st_meta.shape[1]} columns.')
    st_coord.columns = ['x','y']
    st_coord.index = st_coord.index.map(str)
    return st_coord

def check_empty_dict(mydict):
    if not any(mydict.values()):
        raise ValueError(
            "No cell has neighbor, check parameter st_tp")


def convert_str_int_tp(sc_meta,tp_columns):
    sc_meta1 = sc_meta.copy()
    tp_array = sc_meta[tp_columns]
    tp_uniq = pd.DataFrame(set(tp_array))[0]
    idx_dict = {k: v for v, k in enumerate(tp_uniq)}
    sc_meta1['num_tp'] = sc_meta1[tp_columns].map(idx_dict)
    return sc_meta1['num_tp']

def align_lr_gene(self):
    lr_df = self.lr_df
    genes = list(self.sc_adata.var.index)
    lr_df = lr_df[lr_df[0].isin(genes) & lr_df[1].isin(genes)]
    return lr_df

def check_sc_meta(sc_adata):
    '''
    sc_meta should have 'sc_id' and 'spot' columns
    if have x and y return false for running embedding as input coordinates
    '''
    col = sc_adata.obs.columns
    if 'sc_id' not in col:
        raise ValueError(f"Expected sc_id as cell index column in sc meta. Not found.")
    if 'spot' not in col:
        raise ValueError(f"Expected spot as cell's spot column in sc meta. Not found.")
    if sc_adata.obs['spot'].dtype != 'object':
        sc_adata.obs['spot'] = sc_adata.obs['spot'].astype('str')

def check_sc_coord(init_sc_embed):
    if not isinstance(init_sc_embed, pd.DataFrame):
        init_sc_embed = pd.DataFrame(init_sc_embed)
    
    if init_sc_embed.shape[1] == 2:
        sc_coord = init_sc_embed
    else:
        # subset only coordinates and rename
        if 'x' and 'y' in init_sc_embed.columns:
            sc_coord = init_sc_embed[['x','y']]
        elif 'row' and 'col' in init_sc_embed.columns:
            sc_coord = init_sc_embed[['row','col']]
        else:
            raise ValueError(
                f'st_adata expected two columns either name x y or row col to represent spatial coordinates, but got None in {init_sc_embed.shape[1]} columns.')
    sc_coord.columns = ['x','y']
    return sc_coord

def check_st_tp(st_tp):
    if (st_tp != 'visum') and (st_tp != 'st') and (st_tp != 'slide-seq'):
        raise ValueError(
            f'st_tp should be choose among either visum or st or slide-seq, get {st_tp}')

