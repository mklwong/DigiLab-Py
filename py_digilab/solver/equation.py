# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 21:46:36 2018

@author: Martin
"""

import numpy as np
        
def dQSSA(t,y,model):
    M_template = np.zeros([y.shape[0]]*2)
    
    # k1 matrix
    k1 = model.params['k1'].copy()
    k1 = make_matrix(k1,M_template)
    
    # k2 matrix
    k2 = model.params['k2'].copy()
    k2.loc[:,'val'] = y[k2.loc[:,'j_x']]*k2.loc[:,'val'].values
    k2 = make_matrix(k2,M_template)
    
    # dQSSA term matrix     
    Km = model.params['Km'].copy()
    Km.loc[:,'val'] = y[Km.loc[:,'j_x']]*Km.loc[:,'val'].values
    Km = make_matrix(Km,M_template)
     
    # Hill term matrix     
    H = model.params['H'].copy()
    H_ratio = np.power(y[H.loc[:,'E']]/H.loc[:,'K'],H.loc[:,'n'])
    H.loc[:,'val'] = H.loc[:,'val']*H_ratio/(H_ratio+1)
    H = make_matrix(H,M_template)
    top = np.dot(k1,y)+np.dot(k2+H,y)/model.compartment + model.sigma(t)
    dydt = np.dot(np.linalg.inv(np.eye(y.shape[0])+Km),top)
    return dydt

def make_matrix(df,template):
    M = template.copy()
    if df.shape[0]>0:
        df = df.pivot_table(index='y',columns='i_x',values='val',aggfunc='sum')
        M_t = df.reindex(list(range(df.columns.max()+1)),axis=1).reindex(list(range(df.index.max()+1))).fillna(0).values
        M[:M_t.shape[0],:M_t.shape[1]] = M_t
    return M