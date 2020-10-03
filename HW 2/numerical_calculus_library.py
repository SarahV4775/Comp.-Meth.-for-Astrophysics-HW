# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 16:42:09 2020

@author: Sarah
"""
#func = lambda x: x**2 - 45
import numpy as np

def derivative(func):
    '''
    Symmetric derivative of a function initialized as func
    returns dfunc: 
    -----
    Parameters:
        func = function of a single variable

        a, b = start and stop of the interval that the function is plotted
            [a,b]
    '''
    h = 0.001 #value commonly used in calculators (smaller == more effective)
    dfunc = (func(x+h)-func(x-h))/(2*h)
    return (dfunc)

def midpoint_rule(func, a, b, n):
    '''
    Midpoint rule (rectangle rule) for estimating a definate integral using a
        Riemann sum with midpoints of subintervals
    returns M: the value of the integration with the given parameters
    -----
    Parameters:
        a, b = numbers Interval of integration [a,b]
        
        n = Number of subintervals of [a,b]
        
        func = function of a single variable
    '''
    deltx = (b-a)/n #step size
    M = 0
    for i in range(n):
        M += func( (a + (deltx/2)) + (i*deltx) )
    M = M*deltx
    return (M)

def trapezoidal_rule(func, a, b, n):
    '''
    Trapezoidal rule for estimating a definate integral using trapizoids
        instead of rectangles like midpoint rule
    returns T: the value of the integration with the given parameters
    -----
    Parameters:
        a, b = numbers Interval of integration [a,b]
        
        n = Number of subintervals of [a,b]
        
        func = function of a single variable
    '''
    deltx = (b-a)/n #step size
    T = 0
    for i in range(n):
        T += ( deltx/2 )*( func(a + (i*deltx)) + func(a + ((i + 1)*deltx)) )
    return (T)

def simpsons_rule(func, a, b, n):
    '''
    Simpsons rule for estimating a definate integral using 
        piecewise quadratic functions
    returns S: the value of the integration with the given parameters
    -----
    Parameters:
        a, b = numbers Interval of integration [a,b]
        
        n = Number of subintervals of [a,b]
        
        func = function of a single variable
    '''
    deltx = (b-a)/n #step size
    c = 0
    x = a + deltx
    n1 = (n//2) + 1
    for i in range(1, n1):
        c += 4*func(x)
        x += 2*deltx

    x = a + 2*deltx
    n2 = (n//2)
    for i in range(1, n2):
        c += 2*func(x)
        x += 2*deltx
    S = ( deltx/3 )*( func(a) + func(b) + c )
    
    return (S)

