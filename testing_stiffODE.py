# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 21:10:08 2020

@author: Sarah
"""

import numpy as np
import matplotlib.pyplot as plt
from ODE_solver import Euler
from ODE_solver import Heun
from ODE_solver import RK4


def stiffODE(Lamda):
    def eq(t,y):
        dydt = -Lamda*(y - np.cos(t))
        return dydt
    return eq


# Parameters
h = 0.01
yinit = np.array([0.0])
t = np.linspace(0, 10, 1000)
L = 10

Ets, Eys = Euler(stiffODE(L), yinit, t, h)
Hts, Hys = Heun(stiffODE(L), yinit, t, h)
RKts, RKys = RK4(stiffODE(L), yinit, t, h)

# Exact solution, for comparison

def stiff_ODE_exact(lam, x):
    return lam/(1 + lam**2.0) * np.sin(x) + lam**2/(1 + lam**2) * np.cos(x) - lam**2/(1 + lam**2) * np.exp(-lam * x)

yexact = stiff_ODE_exact(L, t)

plt.plot(Ets, Eys, 'r')
plt.plot(Hts, Hys, 'y')
#plt.plot(RKts, RKys, 'k')
plt.plot(t, yexact, 'b')
plt.xlim([0, 2])
plt.ylim([-1, 1])
plt.legend(["Forward Euler method",
            "Heuns Method",
            "Runge-Kutta Method",
            "Exact solution"])
plt.xlabel('Time')
plt.ylabel('y')
plt.title('Stiff ODE')
plt.tight_layout()
#plt.savefig('testing_stiffODE.png')
plt.show()