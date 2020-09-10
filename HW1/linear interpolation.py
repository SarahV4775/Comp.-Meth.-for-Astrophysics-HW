# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 20:15:34 2020

@author: Sarah
"""

import numpy as np
import matplotlib.pyplot as plt

def lin_interpolation(x_new, x_list, y_list):
    '''
    Function that will interpolate the new x value between the two points that 
    are initiated in the begining. will return a new value for y to corispond
    to the new x value.
    -----
    Parameters:
        x_new = New Value of X
        
        x_list = list of initial x values from the text file
        
        y_list = list of initial y values from the text file
    '''
    i = 0
    n=0
    
    while x_list[0] <= x_new[i] <= x_list[len(x_list)-1]:
        y_newlist = []
        if x_list[n] <= x_new[i] <= x_list[n+1]:
            xd = (x_new[i] - x_list[n])/(x_list[n+1] - x_list[n])
            y_new = y_list[n]*(1-xd) + (y_list[n]*xd)
            #y_new = y_list[0]+(x_new[i] - x_list[n])*((y_list[n+1]-y_list[n])/(x_list[n+1]-x_list[n]))
            y_newlist.append(y_new)
            i =+ 1
        else:
            n =+ 1
        
        return y_newlist

#read in text file
data = np.loadtxt('lens_density.txt', skiprows=1)
# X values
X = data[:, 0]
# Y values
N_e = data[:, 1]

#list of new x values
newX = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 
        13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5]

# new y interpolated from the new x values
newY = np.interp(newX, X, N_e)

#plot points
plt.plot(X,N_e, 'ro', label = 'Original points')
plt.plot(newX,newY,'bs', label = 'Interpolated points')
plt.title('Piecewise Linear Interpolation',fontsize=13,color='black')
plt.legend()
plt.grid(True)

#plt.savefig('Piecewise_Linear_Interpolation.png')
plt.show()
