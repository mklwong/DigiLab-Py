# -*- coding: utf-8 -*-
"""
Created on Sun Apr 16 11:23:56 2017

@author: Marty
"""

import numpy as np

def odeFun(t,y):
    k = np.ones((2))
    dy_dt = np.zeros(y.shape)
    dy_dt[0]= k[1]*y[1]-k[0]*y[0]
    dy_dt[1]=-dy_dt[0]
    return(dy_dt)

t = np.linspace(0,10,1000)
yOut = int.odeint(odeFun,[1,0],[0,10])
