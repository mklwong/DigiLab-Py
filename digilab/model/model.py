# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 12:57:28 2018

@author: Martin
"""

import numpy as np
import pandas as pd
from copy import deepcopy

class DigilabModel(object):
    def __init__(self,reactions=[]):
        self.reset_params()
        self._species = {}
        self._compartments = {}
        self.reactions = []
        for reaction in reactions:
            self.add_reaction(reaction)
        
    def copy(self):
        return deepcopy(self)
    
    def add_reaction(self,reaction):
        self.add_species(reaction.reactants+reaction.products+reaction.enzymes)
        self.reactions += [reaction]
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
        return self
    
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
        return self
    
    def get_comp_two_species(self,df):
        if df.shape[0] == 0:
            return []
        comp1 = self.compartment[df.loc[:,'i_x'].values.astype(int)]
        comp2 = self.compartment[df.loc[:,'i_x'].values.astype(int)]
        comp = np.min((comp1,comp2),axis=0)
        return comp
    
    def classify_reaction(self,reaction):
        raise NotImplemented('classify_reaction not implemented.')
    
    def compile_model(self):
        self.reactions = [self.classify_reaction(reaction) for reaction in self.reactions]
        self.parameters = []
        for reaction in self.reactions:
            reaction.compile(self)
            
    def reset_params(self):
        raise NotImplemented('reset_parameter not implemented.')
    
    def consolidate_params(self):
        raise NotImplemented('consolidate_parameter not implemented.')
    
    def get_odefun(self):
        raise NotImplemented('get_odefun not implemented for this model. Must return a function that takes in t, y and model as inputs.')
    