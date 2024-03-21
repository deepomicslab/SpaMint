import subprocess
import os
import pandas as pd
from . import preprocess as pp


def load_spatalk(result, pvalue_thred = 0.005, tp_map = None):
    result['celltype_sender'] = result['celltype_sender'].map(tp_map)
    result['celltype_receiver'] = result['celltype_receiver'].map(tp_map)
    result = result[result['lr_co_ratio_pvalue'] < pvalue_thred].copy()
    result["name"] = result[['ligand','receptor','celltype_sender','celltype_receiver']].apply("-".join, axis=1)
    # some are repeat with different score somehow
    result = result.groupby(['ligand','receptor','celltype_sender','celltype_receiver',"name"]).mean(numeric_only = True).reset_index()
    result["CCI"] = result[['celltype_sender','celltype_receiver']].apply("-".join, axis=1)
    result["LRI"] = result[['ligand','receptor']].apply("-".join, axis=1)
    return result



def runSpaTalk(adata, rscript_executable = '/apps/software/R/4.2.0-foss-2021b/bin/Rscript',
               meta_key = 'celltype',species = 'Human',overwrite = False):
    '''
    This function is to run SpaTalk and add its results on the adata object.
    '''
    save_path = adata.uns['save_path']
    out_f = f'{save_path}/spatalk/'
    tp_key = adata.uns['tp_key']
    if overwrite or not os.path.exists(f'{out_f}/lr_pair.csv'):
        script_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__))) + '/pipelines/'
        r_script_file = f'{script_path}/run_spatalk_lr.R'
        # TODO change name to sc_count.tsv and sc_meta.tsv
        st_dir = f'{save_path}/alter_sc_exp.tsv'
        st_meta_dir = f'{save_path}/spexmod_sc_meta.tsv'
        args = [st_dir,st_meta_dir,st_meta_dir,meta_key,species,out_f]
        subprocess.run([rscript_executable, "--vanilla", r_script_file]+ args)

    if not os.path.exists(f'{out_f}/lr_pair.csv'):
        raise ValueError(f'Error in running SpaTalk, excuation halted.')
    else:
        # spatalk changes '-' to '_' in celltype.
        tp4spatalk = adata.obs['celltype'].str.replace('-','_')
        tp_map = dict(zip(tp4spatalk,adata.obs[tp_key]))
        adata.uns['tp_map_spatalk'] = tp_map

        df = pd.read_csv(f'{out_f}/lr_pair.csv',sep = ',',header=0,index_col=0)
        adata.uns['spatalk'] = load_spatalk(df, pvalue_thred = 0.005,tp_map = tp_map)
    # no need for return adata


def generate_tp_lri(adata,col4Rec,sender_order,receiver_order):
    '''
    This function is used to generate the LRI for each cell type pair.
    '''
    draw_lr = adata.uns['spatalk']
    if sender_order == None:
        sender_order = col4Rec.index.tolist()
    if receiver_order == None:
        receiver_order = col4Rec.columns.tolist()
        
    col4Rec_melt = col4Rec.reset_index()
    col4Rec_melt = col4Rec_melt.melt(id_vars = 'sender_tp',value_vars = col4Rec_melt.columns,var_name = 'receiver_tp',value_name = 'CCI')
    draw_target_pattern = col4Rec_melt.copy()
    draw_target_pattern['LRI'] = 0

    for idx,row in draw_target_pattern.iterrows():
        sender = row['sender_tp']
        receiver = row['receiver_tp']
        count = row['LRI']
        tmp = draw_lr[(draw_lr['celltype_sender'] == sender)&(draw_lr['celltype_receiver'] == receiver)]
        if tmp.shape[0] != 0:
            draw_target_pattern.loc[idx,'LRI'] = int(len(tmp['LRI']))
    draw_target_pattern.columns = ['Sender','Receiver','CCI','LRI']
    # print(draw_target_pattern.head())
    draw_target_pattern = draw_target_pattern[draw_target_pattern['Sender'].isin(sender_order) & 
                                            draw_target_pattern['Receiver'].isin(receiver_order)]

    draw_target_pattern['Sender'] = pd.Categorical(draw_target_pattern['Sender'], categories=sender_order, ordered=True)
    draw_target_pattern['Receiver'] = pd.Categorical(draw_target_pattern['Receiver'], categories=receiver_order, ordered=True)
    draw_target_pattern = draw_target_pattern.sort_values(by=['Sender', 'Receiver'])
    return draw_target_pattern


def generate_cci(adata, target_sender = None, target_receiver = None, return_df = False):
    save_path = adata.uns['save_path']+'/spatalk/'
    tp_key = adata.uns['tp_key']
    tp_map = adata.uns['tp_map_spatalk']
    # print(save_path)

    cellpair = pp.read_csv_tsv(f'{save_path}/cellpair.csv')
    cellpair[['sender_tp','receiver_tp']] = cellpair['Name'].str.split(' -- ',expand = True)
    cellpair['sender_tp'] = cellpair['sender_tp'].map(tp_map)
    cellpair['receiver_tp'] = cellpair['receiver_tp'].map(tp_map)
    inter_tp = set(cellpair['sender_tp'].unique()).intersection(set(adata.obs[tp_key].unique()))
    if len(inter_tp) < len(adata.obs[tp_key].unique())/2:
        print(f'The cell type in the spatalk cellpair file is not consistent with the {tp_key} cell type in the adata object.')
    adata.uns['cellpair'] = cellpair


    nn_df = cellpair.groupby(['sender_tp','receiver_tp']).count().reset_index()
    nn_df = nn_df.pivot(index = 'sender_tp',columns = 'receiver_tp',values = 'Name')
    nn_df.fillna(0,inplace = True)
    nn_df = nn_df.astype(int)
    celltype_nn = nn_df.reset_index()
    celltype_nn = celltype_nn.melt(id_vars = 'sender_tp',value_vars = celltype_nn.columns[1:],var_name = 'receiver_tp',value_name = 'CCI')
    adata.uns['celltype_nn'] = celltype_nn
    # cell type x cell type, CCI percentage
    col4Rec = nn_df/nn_df.sum(axis = 0)
    row4Send = (nn_df.T/nn_df.sum(axis = 1)).T

    draw_rec = generate_tp_lri(adata,col4Rec,target_sender,target_receiver)
    draw_send = generate_tp_lri(adata,row4Send,target_receiver,target_sender)
    adata.uns['rec_per'] = draw_rec
    adata.uns['send_per'] = draw_send
    if return_df:
        return draw_rec,draw_send


def runKEGG(adata, rscript_executable = '/apps/software/R/4.2.0-foss-2021b/bin/Rscript', input_fn = None):
    save_path = adata.uns['save_path']
    out_f = f'{save_path}/kegg/'
    if not os.path.exists(out_f):
        os.makedirs(out_f)
        
    script_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__))) + '/pipelines/'
    r_script_file = f'{script_path}/kegg.R'
    if input_fn:
         args = [out_f,'mouse',input_fn]
    else:
        # no file specified, run all kegg file under out_F
        args = [out_f,'mouse']
    subprocess.run([rscript_executable, "--vanilla", r_script_file]+ args)   


def lri_kegg_enrichment(adata, target_sender = [], target_receiver = [], 
             use_lig_gene = True, use_rec_gene = True,overwrite = False):
    lri_df = adata.uns['spatalk']
    save_path = adata.uns['save_path']
    out_f = f'{save_path}/kegg/'
    if not os.path.exists(out_f):
        os.makedirs(out_f)

    if len(target_sender) == 0:
        target_sender = lri_df['celltype_sender'].unique()
    
    if len(target_receiver) == 0:
        target_receiver = lri_df['celltype_receiver'].unique()
    
    lig_gene = 'L' if use_lig_gene else 'n'
    rec_gene = 'R' if use_rec_gene else 'n'
    
    for sender in target_sender:
        for rec in target_receiver:
            target_lri = lri_df[(lri_df['celltype_sender']== sender) & (lri_df['celltype_receiver'] == rec)].copy()
            # incase of invalid file name
            sender = sender.replace('/','_')
            rec = rec.replace('/','_')
            tmp = pp.lr2kegg(target_lri, use_lig_gene = use_lig_gene, use_rec_gene = use_rec_gene).reset_index()
            input_fn = f'{out_f}/{sender}_{rec}_{lig_gene}_{rec_gene}_kegg.tsv'
            out_fn = f'{out_f}/{sender}_{rec}_{lig_gene}_{rec_gene}_kegg_enrichment.tsv'
            if not os.path.exists(out_fn) or overwrite:
                tmp.to_csv(input_fn,index = True,sep = '\t',header = True)
                runKEGG(adata, rscript_executable = '/apps/software/R/4.2.0-foss-2021b/bin/Rscript', input_fn = input_fn)

    # load
    if 'kegg_enrichment' in adata.uns.keys():
        kegg_res = adata.uns['kegg_enrichment']
        gene_dict = adata.uns['gene_dict']
    else:
        kegg_res = pd.DataFrame()
    geneid = pd.DataFrame()
    # for query_tp in glia:
    for sender in target_sender:
        for rec in target_receiver:
            sender = sender.replace('/','_')
            rec = rec.replace('/','_')
            tmp = pd.read_csv(f"{out_f}/{sender}_{rec}_{lig_gene}_{rec_gene}_kegg_enrichment.tsv",sep = '\t',header=0,index_col=0)
            tmp = pp.filter_kegg(tmp,pval_thred = 0.05)
            tmp['celltype_sender'] = sender
            tmp['celltype_receiver'] = rec
            tmp['used_genes'] = f'{lig_gene}_{rec_gene}'
            kegg_res = pd.concat((kegg_res,tmp))

            tmp = pd.read_csv(f"{out_f}/{sender}_{rec}_{lig_gene}_{rec_gene}_kegg_geneID.tsv",sep = '\t',header=0,index_col=0)
            geneid = pd.concat((geneid,tmp))
    tmp_gene_dict = dict(zip(geneid['ENTREZID'],geneid['SYMBOL']))
    gene_dict.update(tmp_gene_dict)
    kegg_res.index = range(len(kegg_res))
    for index, row in kegg_res.iterrows():
        gene_ids = row['geneID'].split('/')
        symbols = [gene_dict[int(gene_id)] for gene_id in gene_ids]
        new_gene_ids = '/'.join(symbols)
        kegg_res.loc[index, 'geneSymbol'] = new_gene_ids
    adata.uns['kegg_enrichment'] = kegg_res  
    adata.uns['gene_dict'] = gene_dict
    return kegg_res, gene_dict