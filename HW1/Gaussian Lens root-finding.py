# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 13:11:58 2020

@author: Sarah
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText


from root_finding import Bisection
from root_finding import Newton
from root_finding import Secant

def f(x):
    lam = 1.403*10**-12
    N0 = 2.238*10**24
    D = 2.603*10**8
    a = 1
    r = 1.884*10**-26
    dx = x*(1 + (((lam**2)*r*N0*D)/(np.pi*a**2))*np.exp(-(x/a)**2))-2
    #Theta = -(((lam**2)*r*N0)/(np.pi*a**2))*x*np.exp(-(x/a)**2)

    return dx

R = Secant(f, 0 ,5,.001)

print(R)

xprime = np.arange(0, 2.1, 0.1)
#xline = np.arange(0, R, 0.1)

plt.plot((2, 2), (0, 2),'b-', label = 'ray before passing through lens plane ')
plt.plot(xprime,f(xprime), 'r--', label = 'ray after passing through lens plane ')

plt.axhline(0)

plt.legend()
plt.annotate("Lens Plane", (0.25,0.1))
plt.title('Gaussian Lens root-finding',fontsize=13,color='black')
plt.ylabel('[AU]')
plt.xlabel('[AU]')
#plt.savefig('Gaussian_Lens_root-finding.png')
plt.show()
