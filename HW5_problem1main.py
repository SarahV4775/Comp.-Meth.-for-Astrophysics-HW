# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 14:49:49 2020

@author: Sarah
"""

import numpy as np
import matplotlib.pyplot as plt
from NodeClass_BarnesHut import Accel
from NodeClass_BarnesHut import Node
from NodeClass_BarnesHut import tree_contribution


def verlet(points0, points1, t, Accel_Node):
    '''
    symplectic integration method to calculate the next position of the 
    galaxy as it evolves
    ---------
    points0 = list of points from the 1st list of points gal0
    points1 = list of points from the 1st list of points gal1 
    t: time step
    Accel_Node: list of Acceleration values from the Accel func calculated in 
            NodeClass_BarnesHut
    '''
    #x and y values seperated
    x0 = points0[:,0]
    x1 = points1[:,0]
    y0 = points0[:,1]
    y1 = points1[:,1]
    #empty list of new x and y values for the next galaxies
    next_x = np.zeros_like(x0)
    next_y = np.zeros_like(y0)
    
    for i in range(len(x0)):
        n_x = 2*x1[i] - x0[i] + Accel_Node[i]*t**2
        n_y = 2*y1[i] - y0[i] + Accel_Node[i]*t**2
        next_x = np.append(next_x, n_x)
        next_y = np.append(next_y, n_y)
    return next_x, next_y
    
    


gal0 = np.load('galaxies0.npy')
gal1 = np.load('galaxies1.npy')
time = 1000 #years
#list of mass
mass = np.full(len(gal0), 1e12) #solar mass

#number of iterations set for the evolution of the cluster before it prints it
iterations = 10
#initial two cluster lists
prev = gal0
current = gal1
for i in range(iterations):
    A = Accel(current,mass) 
    #calculate the new x and y values
    x, y = verlet(prev, current, time, A)
    New = [x,y]
    #update the prev and current variables to evolve the loop
    prev = current
    current = New

plt.scatter(current[:,0], current[:,1])
#plt.scatter(gal0[:,0],gal0[:,1])


