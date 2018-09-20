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
    def __init__(self,reactions=[]):
        self.reset_params()
        self._species = {}
        self._compartments = {}
        for reaction in reactions:
            self.add_reaction(reaction)
        
    def copy(self):
        return deepcopy(self)
    
    def add_reaction(self,reaction):
        self.add_species(reaction.reactants+reaction.products+reaction.enzymes)
        self.reaction += [reaction]
        return self
    
    def add_species(self,_species):
        staging = {}
        for species in _species:
            if species in staging and not staging[species]._like_(species):
                raise(ValueError("New species '{}' ID same as added species '{}', but names or compartment location are different. Check species input and try again.".format(species,staging[species])))
            elif species in self._species and not self._species[species]._like_(species):
                raise(ValueError("New species '{}' ID same as existing model species '{}', but names compartment location are different. Check species input and try again.".format(species,self._species[species])))
            else:
                staging[species] = species
        self._species.update(staging)
        self.species = self._species.keys()
    
    def add_compartment(self,compartments):
        staging = {}
        for compartment in compartments:
            if compartment in staging and not staging[compartment]._like_(compartment):
                raise(ValueError("New compartment '{}' ID same a added compartment '{}', but names different. Check compartment input and try again.".format(compartment,staging[compartment])))
            elif compartment in self._compartments and not self._compartments[compartment]._like_(compartment):
                raise(ValueError("New compartment '{}' ID same as existing model compartment '{}', but names different. Check compartment input and try again.".format(compartment,self._species[compartment])))
            else:
                staging[compartment] = compartment
        self._compartments.update(staging)
        self.compartments = self._compartments.keys()
    
    def reset_params(self):
        raise NotImplemented('reset_parameter not implemented.')
    
    def consolidate_params(self):
        raise NotImplemented('consolidate_parameter not implemented.')
    
    def get_odefun(self):
        raise NotImplemented('get_odefun not implemented for this model. Must return a function that takes in t, y and model as inputs.')
    
class DigilabModel(DigilabModelTemplate):
    def sigma(self,t):
        return np.asarray([0,0,0,0])
    
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
    
def sigma(t):
    return t