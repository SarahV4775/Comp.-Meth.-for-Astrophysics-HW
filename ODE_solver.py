# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 20:02:46 2020

@author: Sarah
"""
import matplotlib.pyplot as plt
import numpy as np


def Euler(func, y0, t_array, h):
    '''
    ODE Eulers method
    y[i+1] = y[i] + h * func(t[i], y[i])
    ---------------------
    Parameters:
        func = function that we are trying to
        
        y = initial guess for the solution
        
        t = array of vlaues
            (t = np.linspace(start,stop,h))
        
        h = time step
    ----------------------
    Returns:
        an array of solutions for y (ysol)
    '''
    sub_intervals = int((t_array[-1] - t_array[0])/h) # Number of sub-intervals
    numOfODEs = y0.shape[0] # Number of ODEs
    
    t = t_array[0]
    y = y0
    
    # Creates arrays for solutions
    tsol = [t]
    ysol = [y[0]]
    
    for i in range(sub_intervals):
        yprime = func(t, y) # Evaluates dy/dt
        
        for j in range(numOfODEs):
            y[j] = y[j] + (h*yprime[j])
            try:
                ysol = np.append(ysol, y[j]) # store updated y's
            except NameError:
                ysol = y[j]
        t += h # increment by time step
        try:
            tsol = np.append(tsol, t) # store evaluation points
        except NameError:
            tsol = t
     
    return tsol, ysol

def Heun(func, y0, t_array, h):
    '''
    ODE Heuns Method
    y[i+1] = y[i] + h/2 *( f(t[i], y[i]) + f(t[i+1], y[i+1]) )
    ---------------------
    Parameters:
        func = function that we are trying to
        
        y = initial guess for the solution
        
        t = array of vlaues
            (t = np.linspace(start,stop,h))
        
        h = time step
    ----------------------
    Returns:
        an array of solutions for y (ysol)
    '''
    sub_intervals = int((t_array[-1] - t_array[0])/h) # Number of sub-intervals
    numOfODEs = y0.shape[0] # Number of ODEs
    
    t = t_array[0] # Initializes t values
    y = y0 # Initializes variables y
    
    # Creates arrays for solutions
    tsol = [t]
    ysol = [y[0]]
    
    for i in range(sub_intervals):
        y1prime = func(t, y) # Evaluates dy/dt
        
        ypred = y + (h*y1prime)
        
        y2prime = func((t+h), ypred)
        
        for j in range(numOfODEs):
            y[j] = y[j] + (h/2)*(y1prime[j] + y2prime[j])
            try:
                ysol = np.append(ysol, y[j]) # store updated y's
            except NameError:
                ysol = y[j]
        t += h # increment by time step
        try:
            tsol = np.append(tsol, t) # store evaluation points
        except NameError:
            tsol = t
           
    return tsol, ysol


def RK4(func, y0, t_array, h):
    '''
    ODE 4th order Runge-Kutta
    y[i+1] = y[i] + (h/6)*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i])
    k1 = h*f(x[i], y[i])
    k2 = h*f(x[i]+(h/2), y[i]+(k1/2))
    k3 = h*f(x[i]+(h/2), y[i]+(k2/2))
    k4 = h*f(x[i]+h, y[i]+k3)
    ---------------------
    Parameters:
        func = function that we are trying to
        
        y0 = initial guess for the solution
        
        t = array of vlaues
            (t = np.linspace(start,stop,h))
        
        h = time step
    ----------------------
    Returns:
        an array of solutions for y (ysol)
    '''
    sub_intervals = int((t_array[-1] - t_array[0])/h) # Number of sub-intervals
    numOfODEs = y0.shape[0] # Number of ODEs
    
    t = t_array[0] # Initializes t values
    y = y0 # Initializes variables y
    
    # Creates arrays for solutions
    tsol = [t]
    ysol = [y[0]]
    
    for i in range(sub_intervals):
        k1 = func(t, y)
        
        yp2 = y + k1*(h/2)
        
        k2 = func(t+h/2, yp2)
        
        yp3 = y + k2*(h/2)
        
        k3 = func(((t+h)/2), yp3)
        
        yp4 = y + k3*h
        
        k4 = func((t+h), yp4)
        
        for j in range(numOfODEs):
            y[j] = y[j] + (h/6)*(k1[j] + 2*k2[j] + 2*k3[j] + k4[j])
            try:
                ysol = np.append(ysol, y[j]) # store updated y's
            except NameError:
                ysol = y[j]
       
        t += h # increment by time step
        try:
            tsol = np.append(tsol, t) # store evaluation points
        except NameError:
            tsol = t 
        
    return (tsol, ysol)

