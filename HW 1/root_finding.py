# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 14:59:21 2020

Library of root-finding algorithms 

@author: Sarah
"""

import numpy as np


def Bisection(f,a,b,thresh):
    '''
    bisection method of root-finding
    -----
    Parameters:
        f = function that we are trying to approximate
        
        a,b = initial numbers of the interval in which to search for the solution
        
        thresh = the threshold of the final value. when the function of the new guess
            is less then the threshold then the current guess is the solution
    '''
    if f(a)*f(b) > 0:
        print("Bisection method fails. No root.")
        return None

    a_i = a
    b_i = b
    #counter = iteration
    counter = 0
    c = (a_i+b_i)/2
    
    while abs(f(c)) > thresh:
        if f(a)*f(c) < 0:
            b = c
        elif f(b)*f(c) < 0:
            a = c
        elif f(c) == 0:
            return c, counter
        else:
            print("Bisection method fails.")
            #return None
        c = (a+b)/2
        counter += 1
    return c, counter


        
def Newton(f,df,x0,thresh):
    '''
    Newton method of root-finding
    -----
    Parameters:
        f = function that we are trying to approximate
        
        df = derivative of function
        
        x0 = initial guess for the solution
        
        thresh = the threshold of the final value. when the function of the new guess
            is less then the threshold then the current guess is the solution
    '''
    #diff = 1
    counter = 0
    xi = x0
    while abs(f(xi)) > thresh:
        #diff = np.abs(f(xi)/df(xi)-xi)
        x1 = xi - (f(xi)/df(xi))
        xi = x1
        counter += 1
    return xi, counter
    
    
def Secant(f,x0,x1,thresh):
    '''
    Secant method of root-finding
    -----
    Parameters:
        f = function that we are trying to approximate
        
        x0,x1 = initial numbers of the interval in which to search for the solution
        
        thresh = the threshold of the final value. when the function of the new guess
            is less then the threshold then the current guess is the solution
    '''
    counter = 0
    while abs(f(x1)) > thresh:
        x = x1 - f(x1)*((x1 - x0)/(f(x1) - f(x0)))
        x0 = x1
        x1 = x
        counter += 1
    return x, counter

#func = lambda x: x**2 - 45
#
#test = Bisection(func, 6, 7, 0.01)
#print(test)
