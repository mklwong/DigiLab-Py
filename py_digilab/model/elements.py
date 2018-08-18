# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 12:57:28 2018

@author: Martin
"""

from future.utils import string_types

class Species(dict):
    def __init__(self,name,compartment,**kwargs):
        super(Species,self).__init__()
        self['name'] = name
        self['compartment'] = compartment
        self.update(kwargs)
    def __repr__(self):
        x = {x:y for x,y in self.items() if x in ['name','compartment']}
        return "Species({})".format(_pprint(x))
    
class Parameter(dict):
    def __init__(self,name,value=None,**kwargs):
        super(Parameter,self).__init__()
        self['name'] = name
        self['value'] = value
        self.update(kwargs)
    def __repr__(self):
        x = {x:y for x,y in self.items() if x in ['name','value']}
        return "Parameter({})".format(_pprint(x))
                         
class Compartment(dict):
    def __init__(self,name,size=None,**kwargs):
        super(Compartment,self).__init__()
        self['name'] = name
        self['size'] = size
        self.update(kwargs)
    def __repr__(self):
        x = {x:y for x,y in self.items() if x in ['name','size']}
        return "Compartment({})".format(_pprint(x))
    
class Reaction(object):
    def __init__(self,reactants,products,kinetic_law='infer',enzyme=[],**kwargs):
        for x in reactants + products + enzyme:
            assert(isinstance(x,Species))
        self.reactants = reactants
        self.products = products
        self.enzyme = enzyme
        assert(isinstance(kinetic_law,str))
        self.kinetic_law = kinetic_law
        self.extras = kwargs
        
    def __repr__(self):
        r = '['+','.join([x['name'] for x in self.reactants])+']'
        p = '['+','.join([x['name'] for x in self.products])+']'
        e = '['+','.join([x['name'] for x in self.enzyme])+']'
        outstr = "Reaction(reactant={},product={},enzyme={},kinetic_law={},...)"\
                 .format(r,p,e,self.kinetic_law)
        return outstr
    
def _pprint(d):
    out = []
    for key,val in d.items():
        if isinstance(val,string_types):
            out += ["{}='{}'".format(key,val)]
        else:
            out += ["{}={}".format(key,val)]
    return ','.join(out)