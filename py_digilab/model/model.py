# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 12:57:28 2018

@author: Martin
"""

import numpy as np
import pandas as pd
from py_digilab.solver.equation import dQSSA
from copy import deepcopy

class DigilabModelTemplate(object):
    def __init__(self):
        self.reset_params()
    
    def get_odefun(self):
        raise NotImplemented('get_odefun not implemented for this model. Must return a function that takes in t, y and model as inputs.')
    
    def reset_params(self):
        raise NotImplemented('parameter reset not implemented.')
    
    def copy(self):
        raise NotImplemented('copy for this model not implemented.')
    
class DigilabModel(DigilabModelTemplate):
    def __init__(self,species=[],reactions=[],compartments=[]):
        super(DigilabModel,self).__init__()
        self.species = species
        self.reactions = reactions
        self.compartments = compartments
    
    def add_reaction(self,reaction):
        return self
    
    def sigma(self,t):
        return np.asarray([0,0,0,0])
    
    def copy(self):
        new_model = DigilabModel()
        new_model.param = deepcopy(self.params)
        new_model.compartment = deepcopy(self.compartment)
        return new_model
    
    def reset_params(self):
        raw_params = {}
        raw_params['k1'] = pd.DataFrame(columns=['y','i_x','k','p'])
        raw_params['k2'] = pd.DataFrame(columns=['y','i_x','j_x','k','p'])
        raw_params['Km'] = pd.DataFrame(columns=['y','i_x','j_x','Km','p'])
        raw_params['H'] = pd.DataFrame(columns=['y','i_x','E','k','K','n','p'])
        self.raw_params = raw_params
        self.consolidate_params()
        return self
    
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
        
    def get_comp_two_species(self,df):
        if df.shape[0] == 0:
            return []
        comp1 = self.compartment[df.loc[:,'i_x'].values.astype(int)]
        comp2 = self.compartment[df.loc[:,'i_x'].values.astype(int)]
        comp = np.min((comp1,comp2),axis=0)
        return comp
        
    def get_odefun(self):
        return dQSSA
    
class DigilabModel_test(DigilabModel):
    def set_test(self):
        self.raw_params['k1'] = pd.DataFrame({'y':[0,3,0,3],'i_x':[3,3,2,2],'k':[1,-1,-1,1],'p':[1,1,1,1]})
        self.raw_params['Km'] = pd.DataFrame({'y':[2,2,0,0,1,1],'i_x':[0,1,0,1,0,1],'j_x':[1,0,1,0,1,0],'Km':[-2,-2,2,2,2,2],'p':[1,1,1,1,1,1]})
        self.compartment = np.asarray([1,1,1,1])
        self.consolidate_params()
        return self
    
def sigma(t):
    return t