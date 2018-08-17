# -*- coding: utf-8 -*-
"""
Created on Wed Aug 15 21:29:31 2018

@author: Martin
"""

from scipy.integrate import ode
import numpy as np

def make_solver(func):
    def _make_solver(odefunc,tspan,y0,options=None):
        y_out = []
        t_out = []
        def solution_getter(t,y):
           t_out.append(t)
           y_out.append(y.copy()) 
        ode_solver = func(odefunc)
        if len(tspan) > 2:
            ode_solver.set_initial_value(y=y0, t=tspan[0])
            for t in tspan:
                solution_getter(t,ode_solver.integrate(t))
        else:
            ode_solver.set_solout(solution_getter)
            ode_solver.set_initial_value(y=y0, t=tspan[0])
            ode_solver.integrate(tspan[1])
        return np.asarray(t_out), np.asarray(y_out)
    return _make_solver

def ode15s(odefunc,tspan,y0,options=None,tol=0.01):
    ode_solver = ode(odefunc).set_integrator('vode',method='bdf',order=15)
    ode_solver.set_initial_value(y=y0, t=tspan[0])
    
    t_out = [tspan[0]]
    y_out = [np.asarray(y0).copy()]
    
    if len(tspan) > 2:
        for t in tspan[1:]:
            y_out.append(ode_solver.integrate(t))
        return tspan, np.asarray(y_out)
    
    # Initialise
    h = float(tspan[1])/5
    te = tspan[0]
    ye = np.asarray(y0)
    
    while te < tspan[1]:
        max_diff = 1
        ta = te
        h = h*4
        while not (max_diff < tol or ta >= tspan[1]):
            h = h/2
            # Predict based on Newtons
            yp = ye+odefunc(te,ye)*h
            
            # Use solver to get actual
            ta = min(tspan[1],te+h)
            ode_solver.set_initial_value(y=ye, t=te)
            ya = ode_solver.integrate(ta)
            
            # Compare newton and vode prediction
            mask = (ya!=0)&(yp!=0)
            max_diff = np.max((yp[mask]-ya[mask])/ya[mask])
        te = ta
        ye = ya
        t_out.append(te)
        y_out.append(ye)
    return np.asarray(t_out), np.asarray(y_out)

@make_solver
def ode45(odefunc):
   return ode(odefunc).set_integrator('dopri5')

def odefun(t,y):
    return np.asarray([-0.5*y[0]+0.1*y[1],0.5*y[0]-0.1*y[1]])