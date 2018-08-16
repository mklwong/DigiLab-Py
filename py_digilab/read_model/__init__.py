# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import xml.etree.ElementTree as ET
import warnings

# Exceptions
class InputFileError(Exception):
    pass
class DuplicateID(Exception):
    pass

# Warnings
class DuplicateNames(Warning):
    pass
class UnnamedItem(Warning):
    pass


def read_sbml(fname):
    tree = ET.parse(fname)
    root = tree.getroot()
    if 'sbml' not in root.tag:
        raise(InputFileError('Input model must be an SBML file.'))
    #Extract components
    models = [_parseModel(x) for x in root]
    return(models)
        
def _parseModel(model):
    modelParts = [x for x in model]
    compartInd = [ii for ii,x in enumerate(modelParts) if 'Compartments' in x.tag][0]
    speciesInd = [ii for ii,x in enumerate(modelParts) if 'Species' in x.tag][0]
    reactxnInd = [ii for ii,x in enumerate(modelParts) if 'Reaction' in x.tag][0]
    comp = _parseCompartments(modelParts[compartInd])
    specs = _parseSpecies(modelParts[speciesInd])
    rxns = _parseReactions(modelParts[reactxnInd])
    return(comp,specs,rxns)
    
def _parseCompartments(sbmlComprts):
    itemName = 'compartment'
    keepKeys = ['name','size','units']
    return(_parseElement(sbmlComprts,keepKeys,itemName))

def _parseSpecies(sbmlSpecies):
    itemName = 'species'
    keepKeys = ['name','compartment','initialConcentration']
    return(_parseElement(sbmlSpecies,keepKeys,itemName))

def _parseElement(element,keepKeys,itemName):
    attrib= {}
    name2id = {}
    undefName = 1
    for ii,item in enumerate(element):
        if item.attrib['id'] in item.keys():
            raise(DuplicateID('Model ' + itemName + ' have duplicate IDs (ID: ' + str(item.attrib['id']) + '). IDs must be unique.'))
        if 'name' not in item.attrib.keys():
            warnings.warn('Unnamed ' + itemName + ' found with ID: ' + str(item.attrib['id']) + '. Renamed ' + itemName + str(undefName),UnnamedItem)
            item.attrib['name'] = itemName + str(undefName)
            undefName += 1
        if item.attrib['name'] not in name2id.keys():
            name2id[item.attrib['name']] = item.attrib['id']
        else:
            warnings.warn('There are duplicated ' + itemName + ' names in this model. Model will be more difficult to manipulate.',DuplicateNames)
            if type(name2id[item.attrib['name']]) in [int,str] :
                name2id[item.attrib['name']] = [name2id[item.attrib['name']]]
            name2id[item.attrib['name']] += [item.attrib['id']]
        attrib[item.attrib['id']] = {x:item.attrib[x] for x in item.attrib.keys() if x in keepKeys}
    return({'specs':attrib,'dict':name2id})

def _parseReactions(sbmlRxns):
    rxnDict = {}
    attrib= {}
    name2id = {}
    undefName = 1
    for ii,rxn in enumerate(sbmlRxns):
        if rxn.attrib['id'] in rxn.keys():
            raise(DuplicateID('Model reactions have duplicate IDs (ID:' + str(rxn.attrib['id']) + '). IDs must be unique.'))
        rxnDict = {}
        rxnDict['reactants'] = []
        rxnDict['products'] = []
        for jj,item in enumerate(rxn):
            if 'listOfReactants' in item.tag:
                rxnDict['reactants'] = [x.attrib['species'] for x in item]
            elif 'listOfProducts' in item.tag:
                rxnDict['products'] = [x.attrib['species'] for x in item]
            elif 'listOfModifiers' in item.tag:
                rxnDict['modifier'] = [x.attrib['species'] for x in item]
            elif 'kineticLaw' in item.tag:
                for kk,subItem in enumerate(item):
                    if 'parameters' in subItem.tag:
                        print('hi')
                # Custom maths embedded in the SBML file is not implemented yet
                pass
        
        if 'name' not in rxn.attrib.keys():
            rxnName = _mkRxnName(rxnDict,undefName)
            undefName += 1
            warnings.warn('Unnamed reaction found with ID:' + str(rxn.attrib['id']) + '. Renamed ' + rxnName)
            rxn.attrib['name'] = rxnName
        if rxn.attrib['name'] not in name2id.keys():
            name2id[rxn.attrib['name']] = rxn.attrib['id']
        else:
            warnings.warn('There are duplicated reaction names in this model. Model will be more difficult to manipulate.',DuplicateNames)
            if type(name2id[rxn.attrib['name']]) in [int,str] :
                name2id[rxn.attrib['name']] = [name2id[rxn.attrib['name']]]
            name2id[rxn.attrib['name']] += [rxn.attrib['id']]
        rxnDict['name'] = rxn.attrib['name']
        attrib[rxn.attrib['id']] = rxnDict
    return({'specs':attrib,'dict':name2id})

def _mkRxnName(rxnDict,num):
    rxnName = 'rxn' + str(num) + ': ' + '+'.join(rxnDict['reactants']) + ' -> ' +  '+'.join(rxnDict['products'])
    if 'modifier' in rxnDict.keys():
        rxnName += ' | ' + ','.join(rxnDict['modifier'])
    return(rxnName)