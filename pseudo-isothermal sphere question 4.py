# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 19:51:18 2020

@author: Sarah
"""

import matplotlib.pyplot as plt
import numpy as np

from root_finding import Bisection
from root_finding import Newton
from root_finding import Secant

#width (x) of pseudo-isothermal sphere in terms of rc
# 
func = lambda x: ((1+(x/(1)))**(-1/2))-(1/2)
#dfunc = lambda x: (x/1.496*10**13)/((1 + (x/1.496*10**13)**2)**(3/2))

#bi = Bisection(func, 1, 5, 0.001)
#Newt = Newton(func, dfunc, 1, 0.1)
Sec = Secant(func, 0, 5, 0.001)

#print("using the Bisection method, the root value is: ", bi[0])
#print("in ", bi[1], "iterations.")
#print("using the Newtons method, the root value is: ", Newt[0])
#print("in ", Newt[1], "iterations.")
print("using the Secant method, the root value is: ", Sec[0])
print("in ", Sec[1], "iterations.")

root = Sec[0]
t = np.arange(0, 5, 0.1)
#plot points
plt.plot(t,func(t), 'r.-', label = 'Function of a pseudo-isothermal sphere ')
plt.plot(root,func(root),'bv', label = 'Root of the function')
plt.title('Pseudo-isothermal Sphere',fontsize=13,color='black')
plt.legend()
plt.grid(True)
plt.ylabel('[AU]')
plt.xlabel('[AU]')
#plt.savefig('pseudo-isothermal_sphereQuestion4.png')
plt.show()