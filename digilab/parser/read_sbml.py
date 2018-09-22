# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 21:41:20 2018

@author: Martin
"""

from xml.etree import ElementTree as ET
from digilab.model.elements import Species, Reaction, Compartment
import re

def read_sbml(self,file):
    out = ET.parse(file)
    root = out.getroot()
    models = [self.parse_model(x) for x in root if 'model' in x.tag]
    return models

def parse_model(self,model):
    model_elements = {re.sub(r"\{.*?\}","",x.tag):x for x in model}
    compartments = [self.parse_compartment(x.attrib) for x in model_elements['listOfCompartments']]
    species = [self.parse_species(x,compartments) for x in model_elements['listOfSpecies']]
    reactions = [self.parse_reaction(x,species) for x in model_elements['listOfReactions']]
    return {'compartment':compartments,'species':species,'reactions':reactions}

def parse_compartment(self,compartment):
    compartment = dict(compartment)
    compartment['name'] = compartment['id']
    compartment['size'] = float(compartment['size'])
    return Compartment(**compartment)
    
def parse_species(self,species,compartments):
    compartments = {x['name']:x for x in compartments}
    species = dict(species.attrib)
    species['name'] = species.get('name',species['id'])
    species['compartment'] = compartments[species['compartment']]
    return Species(**species)

def parse_reaction(self,reaction,species_list):
    species_list = {x['id']:x for x in species_list}
    reaction = {re.sub(r"\{.*?\}","",x.tag):x for x in reaction}
    reaction['reactants'] = [species_list[dict(x.attrib)['species']] for x in reaction['listOfReactants']]
    reaction['products'] = [species_list[dict(x.attrib)['species']] for x in reaction['listOfProducts']]
    if 'ListOfModifierSpeciesReferences' in reaction.keys():
        reaction['enzymes'] = [species_list[dict(x.attrib)['species']] for x in reaction['ListOfModifierSpeciesReferences']]
        reaction.pop('ListOfModifierSpeciesReferences')
    reaction.pop('listOfReactants')
    reaction.pop('listOfProducts')
    return Reaction(**reaction)

def __repr__(self):
    return 'digilab:SBMLParser()'