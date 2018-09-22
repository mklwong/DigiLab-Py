# -*- coding: utf-8 -*-
"""
Created on Sun Aug 19 16:00:44 2018

@author: Martin
"""

import warnings
from digilab.model.elements import DigilabElement

def like(obj1,obj2):
    """
    This logic function checks if py_digilab objects are like each other. This
    exists because py_digilab model elements equality is determined using the 
    object id. However, two elements may actually be the same based on their 
    definitions (e.g. name and location of a species) while having different
    ids (e.g. if the model is an imported CellDesigner model.)
    
    "Like" is used to identify these similar cases which can identify 
    opportunities to streamline model definitions, identify poorly named 
    species that should be renamed, or identify poorly given element ids.
    
    If elements have poorly named ids, they can be automatically regenerated
    using the element's reset_id method.
    
    This function uses obj1's _like_ method of 
    """
    if not isinstance(obj1,DigilabElement):
        raise TypeError('{} is not of type DigilabElement'.format(obj1))
    if not (obj2,DigilabElement):
        raise TypeError('{} is not of type DigilabElement'.format(obj2))
    if obj1._like_(obj2):
        warn_str = 'Definition of {} is very similary to {} but do not share'\
                   'the same id. Consider making their names distinct or '\
                   'collapsing into one {}'
        warnings.warn(warn_str.format(obj1,obj2,type(obj1)))