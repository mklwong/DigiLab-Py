# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 12:57:28 2018

@author: Martin
"""

from copy import deepcopy
import types
from future.utils import string_types

class DigilabElement(object):
    """
    Base class for DigilabElements: model elements in a Digilab model
    
    Describes how equality is defined for DigilabElements that inherit from 
    dict.
    """
    def __hash__(self):
        return hash(self['id'])
    
    def __eq__(self,other):
        if isinstance(other,type(self)):
            return self['id']==other['id']
        return False
    
    def __ne__(self,other):
        return not self.__eq__(other)
    
    def __call__(self):
        if 'value' not in self:
            raise(KeyError(self._call_err_str))
        return self['value']

class CellElements(DigilabElement,dict):
    """
    Base class for DigilabElements that need parameter of value substitution
    """
    def __init__(self):
        super(CellElements,self).__init__()
        self['sign'] = 1
        
    def copy(self):
        return deepcopy(self)
    
    def _like_(self,other):
        return self['name'] == other['name'] and type(self)==type(other)
    
    def __neg__(self):
        out = self.copy()
        out['sign'] = -out['sign']
        return out
        
    def __mul__(self,other):
        return self['sign']*other
    
class Compartment(CellElements):
    """
    DigilabElements that defines a compartment.
    
    Inputs:
        name: Name of the compartment as a string.
        size: Optional. Default size of the compartment. This value is used in 
              a model if no value provided at model runtime.
        id  : Optional. Identifier for the compartment. If two compartments 
              have the same id they are considered the same. If not passed, 
              compartment id will be generated from name.
    """
    _call_err_str = 'Compartment size not defined'
    
    def __init__(self,name,size=None,**kwargs):
        super(Compartment,self).__init__()
        self['name'] = name
        self['value'] = size
        self.update(kwargs)
        self['id'] = kwargs.get('id',self.reset_id()['id'])
        
    def reset_id(self):
        self['id'] = self['name'].replace(' ','_')
        return self
    
    def __neg__(self):
        return self
    
    def __repr__(self):
        x = {x:y for x,y in self.items() if x in ['name','size']}
        return "Compartment({})".format(_pprint(x))
    
    
class Species(CellElements):
    """
    DigilabElements that defines a species.
    
    Inputs:
        name         : Name of the species as a string.
        compartment  : A Digilab Compartment object
        id           : Optional. Identifier for the species. If two species 
                       have the same id they are considered the same. If not 
                       passed, species id will be generated from species and 
                       compartment names.
        initial_conc : Optional. Default concentration 
    """
    
    _call_err_str = 'Species index not defined. Is your model initialised correctly?'
    
    def __init__(self,name,compartment,**kwargs):
        super(Species,self).__init__()
        self['name'] = name
        assert isinstance(compartment,Compartment)
        self['compartment'] = compartment
        self.update(kwargs)
        self['id'] = kwargs.get('id',self.reset_id()['id'])
    
    def reset_id(self):
        self['id'] = self['name'].replace(' ','_') + '_' + self['compartment']['id'].replace(' ','_')
        return self
    
    def _like_(self,other):
        return self['name'] == other['name'] and type(self)==type(other) and self['compartment'] == other['compartment']
    
    def __repr__(self):
        x = {x:y for x,y in self.items() if x in ['name','compartment','id']}
        if self['sign']==1:
            sign_str = ''
        else:
            sign_str = '-'
        return "{}Species({})".format(sign_str,_pprint(x))

class Parameter(CellElements):
    """
    DigilabElements that defines a Parameter.
    
    Inputs:
        name  : Name of the parameter as a string.
        value : Optional. Default value for this parameter.
    """
    
    _call_err_str = 'Parameter value not defined. Ensure parameter values are passed correctly.'
    
    def __init__(self,name,value=None,**kwargs):
        super(Parameter,self).__init__()
        self['name'] = name
        self['value'] = value
        self.update(kwargs)
        self['id'] = kwargs.get('id',self.reset_id()['id'])
        self._is_copy = False
        self._sign = 1
    
    def reset_id(self):
        self['id'] = self['name'].replace(' ','_')
        return self
    
    def __repr__(self):
        x = {x:y for x,y in self.items() if x in ['name','value']}
        if self['sign']==1:
            sign_str = ''
        else:
            sign_str = '-'
        return "{}Parameter({})".format(sign_str,_pprint(x))
    
class Reaction(DigilabElement):
    compile = None
    def __init__(self,reactants,products,kinetic_law='infer',enzymes=[],xsection_ratio=1,**kwargs):
        for x in reactants + products + enzymes:
            assert(isinstance(x,Species))
        self.reactants = reactants
        self.products = products
        self.enzymes = enzymes
        assert(isinstance(kinetic_law,string_types))
        self.kinetic_law = kinetic_law
        self.extras = kwargs
        self.xsection_ratio = xsection_ratio
    
    def apply_rule(self,func):
        self.compile = types.MethodType(func,self)
    
    def _like_(self,other):
        return self.__eq__(other)
    
    def __eq__(self,other):
        pass
    
    def __ne__(self,other):
        return not self.__eq__(other)
    
    def __repr__(self):
        r = '['+','.join([x['name'] for x in self.reactants])+']'
        p = '['+','.join([x['name'] for x in self.products])+']'
        e = '['+','.join([x['name'] for x in self.enzyme])+']'
        outstr = "Reaction(reactant={},product={},enzyme={},kinetic_law={},...)"\
                 .format(r,p,e,self.kinetic_law)
        return outstr
    
    def __call__(self):
        raise(TypeError("'Reaction' object is not callable"))
    
def _pprint(d):
    out = []
    for key,val in d.items():
        if isinstance(val,string_types):
            out += ["{}='{}'".format(key,val)]
        else:
            out += ["{}={}".format(key,val)]
    return ','.join(out)