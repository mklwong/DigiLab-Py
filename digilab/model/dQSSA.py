# -*- coding: utf-8 -*-
"""
Created on Sun Aug 19 22:12:58 2018

@author: Martin
"""

import pandas as pd
from digilab.model import DigilabModel

class DigilabModeldQSSA(DigilabModel):
    def sigma(self,t):
        return np.asarray([0,0,0,0])
    
    def reset_params(self):
        raw_params = {}
        raw_params['k0'] = pd.DataFrame(columns=['y','k','p'])
        raw_params['k1'] = pd.DataFrame(columns=['y','i_x','k','p'])
        raw_params['k2'] = pd.DataFrame(columns=['y','i_x','j_x','k','p'])
        raw_params['Km'] = pd.DataFrame(columns=['y','i_x','j_x','Km','p'])
        raw_params['H'] = pd.DataFrame(columns=['y','i_x','E','k','K','n','p'])
        self.raw_params = raw_params
        self.consolidate_params()
        return self
    
    def classify_reaction(self,reaction):
        reaction.compile = None
        # Classify the reactions
        if len(reaction.reactants) == 0:
            if len(reaction.enzymes) == 0:
                reaction.compile = synthesis
            elif len(reaction.enzymes) == 1:
                reaction.compile = enzymatic_synthesis
        elif len(reaction.reactants) == 1:
            if len(reaction.enzymes) == 0:
                reaction.compile = unimolecular
            elif len(reaction.enzymes) == 1:
                reaction.compile = dQSSA_enzymatic
        elif len(reaction.reactants) == 2:
            if len(reaction.enzymes) == 0:
                reaction.compile = bimolecular
        # Raise error if no rule function applied to reaction
        if reaction.compile is None:
            raise(ValueError('Reaction {} does not fit into any implemented reaction patterns. Please review this reaction.'))
        return reaction
    
    def consolidate_params(self):
        self.params = deepcopy(self.raw_params)
        tmp_df = self.params['k1']
        tmp_df.loc[:,'val'] = tmp_df.loc[:,'k']*tmp_df.loc[:,'p']
        tmp_df = self.params['k2']
        comp = self.get_comp_two_species(tmp_df)
        tmp_df.loc[:,'val'] = comp*tmp_df.loc[:,'k']*tmp_df.loc[:,'p']
        tmp_df = self.params['Km']
        comp = self.get_comp_two_species(tmp_df)
        tmp_df.loc[:,'val'] = comp*tmp_df.loc[:,'p']/(tmp_df.loc[:,'Km']*self.compartment[tmp_df.loc[:,'y'].values.astype(int)])
        tmp_df = self.params['H']
        comp = self.get_comp_two_species(tmp_df)
        tmp_df.loc[:,'val'] = comp*tmp_df.loc[:,'k']*tmp_df.loc[:,'p']
        return self
        
    def equation(self,t,y):
        M_template = np.zeros([y.shape[0]]*2)
        
        # k1 matrix
        k1 = self.params['k1'].copy()
        k1 = make_matrix(k1,M_template)
        
        # k2 matrix
        k2 = self.params['k2'].copy()
        k2.loc[:,'val'] = y[k2.loc[:,'j_x']]*k2.loc[:,'val'].values
        k2 = make_matrix(k2,M_template)
        
        # dQSSA term matrix     
        Km = self.params['Km'].copy()
        Km.loc[:,'val'] = y[Km.loc[:,'j_x']]*Km.loc[:,'val'].values
        Km = make_matrix(Km,M_template)
         
        # Hill term matrix     
        H = self.params['H'].copy()
        H_ratio = np.power(y[H.loc[:,'E']]/H.loc[:,'K'],H.loc[:,'n'])
        H.loc[:,'val'] = H.loc[:,'val']*H_ratio/(H_ratio+1)
        H = make_matrix(H,M_template)
        top = np.dot(k1,y)+np.dot(k2+H,y)/model.compartment + model.sigma(t)
        dydt = np.dot(np.linalg.inv(np.eye(y.shape[0])+Km),top)
        return dydt

def synthesis(model):
    pass
    
def enzymatic_synthesis(model):
    pass

def dQSSA_enzymatic(model):
    pass

def unimolecular(model):
    y = [reaction.reactant[0],reaction.product[0]]
    x = [reaction.reactant[0],reaction.reactant[0]]
    df = [[pd.DataFrame({'y':[reaction.reactant[0]],'i_x':[y],'k':[],'p':[reaction.]});

def bimolecular(model):
       
def make_matrix(df,template):
    M = template.copy()
    if df.shape[0]>0:
        df = df.pivot_table(index='y',columns='i_x',values='val',aggfunc='sum')
        M_t = df.reindex(list(range(df.columns.max()+1)),axis=1).reindex(list(range(df.index.max()+1))).fillna(0).values
        M[:M_t.shape[0],:M_t.shape[1]] = M_t
    return M