from scipy.spatial import distance_matrix
from . import optimizers
from . import utils

import pandas as pd
import numpy as np
import os

class SproutImpute:
    '''
    Created on 2022/11/09
    Updated on 2022/12/20
    SPROUT_impute_v2
        1. Gradient <0 =0
    SPROUT_impute_v3 2022/11/22
        1. Term1 was equalized by dividing sc_spot_sum by cell number per spot
        2. Term3 was equalized by: 
            1) np.sqrt(aff/2); 
            2) summed expression of pair RL to the mean when calculating rl_agg
    SPROUT_impute_v4 2022/11/28
        1. Term3 LRagg calculation
            Given that a gene g1 can be consider as both L and R
            Add the sum of L(g1) and R(g1)
        2. Add regulization term \sum X^2 in gradient descent
        3. Change Learning rate parameter from gamm to ETA
        4. Add hyperparameter DELTA for the newly add regulization term
        5. Change the naming format of the hyperparameters
    SPROUT_impute_v5 2022/12/20
        1. Added weight on hvg genes
        2. Updated the calculation of term3 with neighobring-indicator mat, accelerate from 22s to 2s.
        3. Remove parameter p_dist, change to auto scale.
        4. Correct the scale method from same max each cell to same sum each cell.
        5. Added inital adjustment of each term's weight
    SPROUT_impute_v6 2023/03/14
        1. Added center shift embedding
        2. Add a new neighbor calculation method to adapt to slide-seq data
        3. Assign coordinates if is slide-seq data 
        4. Adapt to spatalk input
        5. Enable user-specific init embedding
        6. Deleted st_coord parameters, subset from st_adata

    st_tp: choose among either visum or st or slide-seq
    '''
    def __init__(self, save_path, st_adata, sc_ref, sc_adata, lr_df, 
                 alpha, beta, gamma, delta, eta, st_tp = 'st', init_sc_embed = False,
                 iteration = 3, k= 2, W_HVG = 2,
                 left_range = 1, right_range = 2, steps = 1, dim = 2):
        self.save_path = save_path +'/'
        self.st_adata = st_adata
        self.lr_df = lr_df
        self.st_tp = st_tp
        self.sc_adata = sc_adata
        self.sc_ref = sc_ref
        # v6 init sc embedding
        self.init_sc_embed = init_sc_embed

        self.ALPHA = alpha
        self.BETA = beta
        self.GAMMA = gamma
        self.DELTA = delta
        self.ETA = eta
        self.K = k
        # v5 weight of hvg genes
        self.W_HVG = W_HVG
        # embedding
        self.iteration = iteration
        self.left_range = left_range
        self.right_range = right_range
        self.steps = steps
        self.dim = dim
    
    def _check_params(self):
        if not os.path.exists(self.save_path):
            os.makedirs(self.save_path)
        
        # utils.check_int(self.ALPHA, self.BETA, self.DELTA, self.ETA, self.iteration, self.K,self.left_range, self.right_range, self.steps, self.dim)
        # utils.check_int(self.iteration, self.K, self.left_range, self.right_range, self.steps, self.dim)
        # utils.check_num(self.ALPHA, self.BETA, self.DELTA, self.ETA)
        utils.check_st_tp(self.st_tp)
        utils.check_sc_meta(self.sc_adata)
        self.st_coord = utils.check_st_coord(self.st_adata)
        self.lr_df = utils.align_lr_gene(self)
        print('Parameters checked!')
        
    
    def run_gradient(self):
        # v5
        self.term1_df,self.loss1 = optimizers.cal_term1(self.alter_sc_exp,self.sc_meta,self.st_exp,self.svg,self.W_HVG)
        print('First-term calculation done!')
        self.term2_df,self.loss2 = optimizers.cal_term2(self.alter_sc_exp,self.sc_ref)
        print('Second-term calculation done!')
        '''
        3. Third term, closer cells have larger affinity
        '''
        if not (self.st_tp == 'slide-seq' and hasattr(self, 'sc_knn')):
            # if slide-seq and already have found sc_knn 
            # dont do it again
            self.sc_dist = pd.DataFrame(distance_matrix(self.sc_coord,self.sc_coord), index = self.alter_sc_exp.index, columns = self.alter_sc_exp.index)
            # 3.2 get c' = N(c)
            self.sc_knn = optimizers.findCellKNN(self.st_coord,self.st_tp,self.sc_meta,self.sc_coord,self.K)
            utils.check_empty_dict(self.sc_knn)
        # 3.3 get the paring genes (g') of gene g for each cells
        self.rl_agg = optimizers.generate_LR_agg(self.alter_sc_exp,self.lr_df)
        # 3.4 get the affinity
        self.aff = optimizers.calculate_affinity_mat(self.lr_df, self.alter_sc_exp.T)
        np.fill_diagonal(self.aff,0)
        self.aff = pd.DataFrame(self.aff, index = self.sc_meta.index, columns=self.sc_meta.index)
        # 3.5 Calculate the derivative
        self.term3_df,self.loss3 = optimizers.cal_term3(self.alter_sc_exp,self.sc_knn,self.aff,self.sc_dist,self.rl_agg)
        print('Third term calculation done!')
        '''
        4. Fourth term, towards spot-spot affinity profile
        '''
        self.term4_df,self.loss4 = optimizers.cal_term4(self.st_exp,self.sc_knn,self.st_aff_profile_df,self.alter_sc_exp,self.sc_meta,self.spot_cell_dict,self.lr_df)
        print('Fourth term calculation done!')
        '''
        5. Fifth term, norm2 regulization
        '''
        self.term5_df,self.loss5 = optimizers.cal_term5(self.alter_sc_exp)
        return


    def init_grad(self):
        if isinstance(self.init_sc_embed, pd.DataFrame):
            self.sc_coord = utils.check_sc_coord(self.init_sc_embed)
            print('Using user provided init sc_coord.')
        else:
            print('Init sc_coord by affinity embedding...')
            # TODO 减少aff计算次数；使用sparse array
            self.sc_coord,_,_,_ = optimizers.aff_embedding(self.alter_sc_exp,self.st_coord,self.sc_meta,self.lr_df,
                            self.save_path,self.left_range,self.right_range,self.steps,self.dim)
        self.run_gradient()
        # v5 calculte the initial loss of each term to balance their force.
        adj2,adj3,adj4,adj5 = optimizers.loss_adj(self.loss1,self.loss2,self.loss3,self.loss4,self.loss5)
        self.ALPHA,self.BETA,self.GAMMA,self.DELTA = self.ALPHA*adj2,self.BETA*adj3,self.GAMMA*adj4,self.DELTA*adj5
        print('Hyperparameters adjusted.')


    def gradient_descent(self):
        print('Running v7 now...')
        self._check_params()
        ######### init ############
        # v5
        self.svg = optimizers.get_hvg(self.st_adata)
        self.sc_ref = np.array(self.sc_ref.loc[self.sc_adata.obs['sc_id']])
        self.alter_sc_exp = self.sc_adata.to_df()
        self.st_exp = self.st_adata.to_df()
        self.sc_meta = self.sc_adata.obs.copy()
        self.spot_cell_dict = self.sc_meta.groupby('spot').apply(optimizers.apply_spot_cell).to_dict()
        self.st_aff_profile_df = optimizers.cal_aff_profile(self.st_exp, self.lr_df)
        res_col = ['loss1','loss2','loss3','loss4','loss5','total']
        result = pd.DataFrame(columns=res_col)
        self.init_grad()
        # loss = self.loss1 + self.ALPHA*self.loss2 + self.BETA*self.loss3 + self.GAMMA*self.loss4 + self.DELTA*self.loss5
        # tmp = pd.DataFrame(np.array([[self.loss1,self.ALPHA*self.loss2,self.BETA*self.loss3,self.GAMMA*self.loss4,self.DELTA*self.loss5,loss]]),columns = res_col, index = [0])
        # result = pd.concat((result,tmp),axis=0)
        ######### init done ############
        for ite in range(self.iteration):
            print(f'-----Start iteration {ite} -----')
            # TODO 减少aff计算次数；使用sparse array
            if self.st_tp != 'slide-seq':
                self.sc_coord,_,_,_ = optimizers.aff_embedding(self.alter_sc_exp,self.st_coord,self.sc_meta,self.lr_df,
                                self.save_path,self.left_range,self.right_range,self.steps,self.dim)
            self.run_gradient()
            gradient = self.term1_df - self.ALPHA*self.term2_df + self.BETA*self.term3_df + self.GAMMA*self.term4_df + self.DELTA*self.term5_df
            self.alter_sc_exp = self.alter_sc_exp - self.ETA * gradient
            self.alter_sc_exp[self.alter_sc_exp<0] = 0
            # v2 added 
            print(f'---{ite} self.loss4 {self.loss4} self.GAMMA {self.GAMMA} self.GAMMA*self.loss4 {self.GAMMA*self.loss4}')
            loss = self.loss1 + self.ALPHA*self.loss2 + self.BETA*self.loss3 + self.GAMMA*self.loss4 + self.DELTA*self.loss5
            tmp = pd.DataFrame(np.array([[self.loss1,self.ALPHA*self.loss2,self.BETA*self.loss3,self.GAMMA*self.loss4,self.DELTA*self.loss5,loss]]),columns = res_col, index = [ite])
            result = pd.concat((result,tmp),axis=0)
            print(f'---In iteration {ite}, the loss is:loss1:{self.loss1:.5f},loss2:{self.loss2:.5f},loss3:{self.loss3:.5f},', end="")
            print(f'loss4:{self.loss4:.5f},loss5:{self.loss5:.5f}.')
            print(f'---In iteration {ite}, the loss is:loss1:{self.loss1:.5f},loss2:{self.ALPHA*self.loss2:.5f},loss3:{self.BETA*self.loss3:.5f},', end="")
            print(f'loss4:{self.GAMMA*self.loss4:.5f},loss5:{self.DELTA*self.loss5:.5f}.')
            print(f'The total loss after iteration {ite} is {loss:.5f}.')
        ### v5 add because spatalk demo
        self.alter_sc_exp[self.alter_sc_exp < 1] = 0  
        self.alter_sc_exp.to_csv(f'{self.save_path}/alter_sc_exp.tsv',sep = '\t',header=True,index=True)
        result.to_csv(f'{self.save_path}/loss.tsv',sep = '\t',header=True,index=True)
        self.result = result
        ### v6 add for center shift  
        if self.st_tp != 'slide-seq':
            self.sc_coord,_,_,_ = optimizers.aff_embedding(self.alter_sc_exp,self.st_coord,self.sc_meta,self.lr_df,
                                self.save_path,self.left_range,self.right_range,self.steps,self.dim)
            _, sc_spot_center = optimizers.sc_prep(self.st_coord, self.sc_meta)
            self.sc_meta[['st_x','st_y']] = sc_spot_center
            self.sc_meta = optimizers.center_shift_embedding(self.sc_coord, self.sc_meta, max_dist = 1)
        
        return self.alter_sc_exp,self.sc_meta