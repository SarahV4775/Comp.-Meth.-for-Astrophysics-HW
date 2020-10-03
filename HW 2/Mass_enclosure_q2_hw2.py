# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 13:49:29 2020

@author: Sarah
"""

import matplotlib.pyplot as plt
import numpy as np

from numerical_calculus_library import derivative

c1 = 15 
c2 = 5

v200_1 = 160 #[km/s]
v200_2 = 200

r200 = 230000 #[pc]
G = 4.3*10**-3 #[pc*Mo^-1*(km/s)^2]
x = np.arange(0, 3*10**5, 1)
def M_enc(r,c,v200):
    x = r/r200
    vc = v200*((1/x)*((np.log(1 + c*x) - ((c*x)/(1+c*x)))/(np.log(1 + c) - ((c)/(1+c)))))**(1/2)
    Menc = (r*(vc**2))/G
    return (Menc)

def M_fun(r):
    ro = 1/(((r*c1)/r200)*(1+((r*c1)/r200))**2)
    m = ro*(4/3)*np.pi(r**3)
    return m
der = derivative(M_fun(x))

#print(M_enc(2.3*10**5))
#print(der)

#plot points
plt.plot(x, M_enc(x,c1,v200_1), 'r-', label = 'c = 15, V_200 = 160')
plt.plot(x, M_enc(x,c1,v200_2), 'b-', label = 'c = 15, V_200 = 200 ')
plt.plot(x, M_enc(x,c2,v200_1), 'y-', label = 'c = 5, V_200 = 160')


plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.title('Mass Profile',fontsize=13,color='black')
plt.legend()
plt.grid(True)
plt.ylabel('Mass [$M_{\odot}$]')
plt.xlabel('Radius [pc]')
#plt.savefig('Mass_Profile.png')
plt.show()
