# -*- coding: utf-8 -*-
"""
Created on Sun Aug 19 22:12:58 2018

@author: Martin
"""

import pandas as pd

def unimolecular_reaction(self):
    y = [reaction.reactant[0],reaction.product[0]]
    x = [reaction.reactant[0],reaction.reactant[0]]
    df = [[pd.DataFrame({'y':[reaction.reactant[0]],'i_x':[y],'k':[],'p':[reaction.]});
           
           
