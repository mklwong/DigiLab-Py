# -*- coding: utf-8 -*-
"""
Created on Sun Aug 19 22:12:58 2018

@author: Martin
"""

import pandas as pd
from digilab.model import DigilabModel
from digilab.model.elements import Parameter
import numpy as np
from copy import deepcopy

class DigilabModeldQSSA(DigilabModel):
    def sigma(self,t):
        return np.asarray([0,0,0,0])
    
    def reset_params(self):
        raw_params = {}
        raw_params['k0'] = pd.DataFrame(columns=['y','k'])
        raw_params['k1'] = pd.DataFrame(columns=['y','i_x','k'])
        raw_params['k2'] = pd.DataFrame(columns=['y','i_x','j_x','k','r'])
        raw_params['Km'] = pd.DataFrame(columns=['y','i_x','j_x','K','r'])
        raw_params['H']  = pd.DataFrame(columns=['y','i_x','j_x','k','K','n','r'])
        self.raw_params = raw_params
        self.consolidate_params()
        return self
    
    def classify_reaction(self,reaction):
        # Classify the reactions
        if reaction.kinetic_law == 'infer':
            if len(reaction.reactants) == 0:
                if len(reaction.enzymes) == 0:
                    reaction.kinetic_law = 'synthesis'
                elif len(reaction.enzymes) == 1:
                    reaction.kinetic_law = 'enzymatic_synthesis'
            elif len(reaction.reactants) == 1:
                if len(reaction.enzymes) == 0:
                    reaction.kinetic_law = 'unimolecular'
                elif len(reaction.enzymes) == 1:
                    reaction.kinetic_law = 'enzymatic_unimolecular'
            elif len(reaction.reactants) == 2:
                if len(reaction.enzymes) == 0:
                    reaction.kinetic_law = 'bimolecular'
        
        # Apply the relevant compile reactions
        if reaction.kinetic_law == 'synthesis':
            reaction.apply_rule(synthesis)
        elif reaction.kinetic_law == 'enzymatic_synthesis':
            reaction.apply_rule(enzymatic_synthesis)
        elif reaction.kinetic_law == 'unimolecular':
            reaction.apply_rule(unimolecular)
        elif reaction.kinetic_law == 'enzymatic_unimolecular':
            reaction.apply_rule(enzymatic_unimolecular)
        elif reaction.kinetic_law == 'bimolecular':
            reaction.apply_rule(bimolecular)
        elif reaction.kinetic_law == 'hill_unimolecular':
            reaction.apply_rule(hill_unimolecular)
            
        # Raise error if no rule function applied to reaction
        if reaction.compile is None:
            raise(ValueError('Reaction {} does not fit into any implemented reaction patterns. Please review this reaction.'))
        return reaction
    
    def consolidate_params(self):
        self.params = deepcopy(self.raw_params)
        tmp_df = self.params['k1']
        tmp_df.loc[:,'val'] = tmp_df.loc[:,'k']*tmp_df.loc[:,'r']
        tmp_df = self.params['k2']
        comp = self.get_comp_two_species(tmp_df)
        tmp_df.loc[:,'val'] = comp*tmp_df.loc[:,'k']*tmp_df.loc[:,'r']
        tmp_df = self.params['Km']
        comp = self.get_comp_two_species(tmp_df)
        tmp_df.loc[:,'val'] = comp*tmp_df.loc[:,'r']/(tmp_df.loc[:,'m']*self.compartment[tmp_df.loc[:,'y'].values.astype(int)])
        tmp_df = self.params['H']
        comp = self.get_comp_two_species(tmp_df)
        tmp_df.loc[:,'val'] = comp*tmp_df.loc[:,'k']*tmp_df.loc[:,'r']
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
        H_ratio = np.power(y[H.loc[:,'j_x']]/H.loc[:,'K'],H.loc[:,'n'])
        H.loc[:,'val'] = H.loc[:,'val']*H_ratio/(H_ratio+1)
        H = make_matrix(H,M_template)
        top = np.dot(k1,y)+np.dot(k2+H,y)/self.compartment + self.sigma(t)
        dydt = np.dot(np.linalg.inv(np.eye(y.shape[0])+Km),top)
        return dydt

def synthesis(self,model):
    k0 = Parameter(name='k0|{}'.format(rxn=self))
    df = pd.DataFrame({'y':self.products,'k':[k0]*len(self.products)})
    self.parameters = [k0]
    model.parameters += self.parameters
    model.raw_params['k0'] = pd.concat([model.raw_params['k0'],df],axis=0)
    return self

def enzymatic_synthesis(self,model):
    k0 = Parameter(name='k0|{}'.format(rxn=self))
    df = pd.DataFrame({'y':self.products,\
                       'i_x':self.enzymes*len(self.products),\
                       'k':[k0]*len(self.products)})
    self.parameters = [k0]
    model.parameters += self.parameters
    model.raw_params['k1'] = pd.concat([model.raw_params['k1'],df],axis=0)
    return self

def unimolecular(self,model):
    k1 = Parameter(name='k1|{}'.format(rxn=self))
    df1 = pd.DataFrame({'y':self.products,\
                       'i_x':self.reactant*len(self.products),\
                       'k':[k1]*len(self.products)})
    df2 = pd.DataFrame({'y':self.reactant,\
                       'i_x':self.reactant,\
                       'k':[-k1]})
    df = pd.concat([df1,df2],axis=0)
    self.parameters = [k1]
    model.parameters += self.parameters
    model.raw_params['k1'] = pd.concat([model.raw_params['k1'],df],axis=0)
    return self

def enzymatic_unimolecular(self,model):
    pass

def hill_unimolecular(self,model):
    pass

def bimolecular(self,model):
    k2 = Parameter(name='k2|{}'.format(rxn=self))
    df1 = pd.DataFrame({'y':self.products,\
                       'i_x':self.reactant*len(self.products),\
                       'k':[k2]*len(self.products)})
    df2 = pd.DataFrame({'y':self.reactant,\
                       'i_x':self.reactant,\
                       'k':[-k2]})
    df = pd.concat([df1,df2],axis=0)
    self.parameters = [k2]
    model.parameters += self.parameters
    model.raw_params['k1'] = pd.concat([model.raw_params['k1'],df],axis=0)
    return self

def make_matrix(df,template):
    M = template.copy()
    if df.shape[0]>0:
        df = df.pivot_table(index='y',columns='i_x',values='val',aggfunc='sum')
        M_t = df.reindex(list(range(df.columns.max()+1)),axis=1).reindex(list(range(df.index.max()+1))).fillna(0).values
        M[:M_t.shape[0],:M_t.shape[1]] = M_t
    return M