# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 21:55:38 2020

@author: Sarah
"""
import matplotlib.pyplot as plt
import numpy as np

from root_finding import Bisection
from root_finding import Newton
from root_finding import Secant

#width (x) of pseudo-isothermal sphere in terms of rc
# 
func = lambda u: ((1+(u**2))**(-1/2))-(1/2)
dfunc = lambda u: u/((1 + u**2)**(3/2))

bi = Bisection(func, 1, 5, 0.001)
Newt = Newton(func, dfunc, 1, 0.1)
Sec = Secant(func, 1, 5, 0.001)



print("using the Bisection method, the root value is: ", bi[0])
print("in ", bi[1], "iterations.")
print("using the Newtons method, the root value is: ", Newt[0])
print("in ", Newt[1], "iterations.")
print("using the Secant method, the root value is: ", Sec[0])
print("in ", Sec[1], "iterations.")

Broot = bi[0]
Nroot = Newt[0]
Sroot = Sec[0]
t = np.arange(0, 5, 0.1)
#t1 = np.arange(0,1, 0.1)

plt.plot(t, func(t))
plt.plot(Broot, func(Broot), 'k|', markeredgewidth = 3, markersize = 10, label="Bisection method root")
plt.plot(Nroot, func(Nroot), 'r|', markeredgewidth = 3, label="Newtons method root")
plt.plot(Sroot, func(Sroot), 'b|', markeredgewidth = 3, label="Secant method root")
plt.title('Piecewise Linear Interpolation',fontsize=13,color='black')
plt.ylabel('Y')
plt.xlabel('X')

plt.legend()

#plt.plot(t1,x(t1))
plt.grid(True)
#plt.savefig('Pseudo-Isothermal Sphere.png')

plt.show()
