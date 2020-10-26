# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 16:53:10 2020

@author: Sarah
"""

import numpy as np
import matplotlib.pyplot as plt

def potential(x,y,numpoints, cenx, ceny, rad):
    M = 1e12 * numpoints #* (4/3)*np.pi*rad**3 #galaxies all with 10^12 solar mass
    rx = (x - cenx)
    ry = (y - ceny)
    a = 1 #scale parameter Plummer Radius
    G = 4.301e-9 #km**2 Mpc MSun**-1 s
    #Plummer model for 
    phix = -G*M/(np.sqrt((rx**2) + (a**2)))
    phiy = -G*M/(np.sqrt((ry**2) + (a**2)))
    phi = abs(phix + phiy)
    return(phi)

x = np.linspace(0, 10, 10)
y = np.linspace(0, 10, 10)
[X, Y] = np.meshgrid(x, y)
Z = np.zeros_like(X) #empty list of potentials

numPoints1 = 400 #400 galaxies
numPoints2 = 255 #255 galaxies
cen1x, cen1y = 7,7 #1st cluster centered at (7,7)
cen2x,cen2y = 2, 8  #2nd cluster centered at (8,2)
rad1 = 3
rad2 = 2

for i in range(0, 10):
    for j in range(0, 10):
        Phi1 = potential(i, j, numPoints1, cen1x, cen1y, rad1)
        Phi2 = potential(i, j, numPoints2, cen2x, cen2y, rad2)
        Phitot = Phi1+ Phi2 #adding both potentials
        Z[i,j] = Phitot

plt.contourf(X,Y,Z,1000)
plt.xlabel('x Mpc')
plt.ylabel('y Mpc')
plt.title('Gravitational Potential')
plt.colorbar()
#plt.savefig('Grav_potential_Q2.png')