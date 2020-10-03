# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 20:25:55 2020

@author: Sarah
"""
import matplotlib.pyplot as plt
import numpy as np
from ODE_solver import Euler
from ODE_solver import Heun
from ODE_solver import RK4

from scipy.integrate import odeint



def ODEpend(b, c):
    def eq(t, y):
        theta, omega = y
        dydt =np.array([omega, -b*omega - c*np.sin(theta)])
        return dydt
    return eq
    
    

def Odeint_pend(y, t, b, c):
    theta, omega = y
    dydt = [omega, -b*omega - c*np.sin(theta)]
    return dydt

# Parameters
b = 0.25
c = 5.0
h = 0.001
y0 = np.array([(np.pi - 0.1), 0.0])
t = np.linspace(0, 10, 1000)

sol = odeint(Odeint_pend, y0, t, args=(0.25, 5.0))

Ets, Eys = Euler(ODEpend(b, c), y0, t, h)
Hts, Hys = Heun(ODEpend(b, c), y0, t, h)
RKts, RKys = RK4(ODEpend(b, c), y0, t, h)

#Plotting Eulers Method
plt.plot(Ets, Eys[np.linspace(0, 19998, 10001).astype(int)], label='omega euler') # odd
plt.plot(Ets, Eys[np.linspace(1, 19999, 10001).astype(int)], label='theta euler') # even

#Plotting Heuns Method
plt.plot(Hts, Hys[np.linspace(0, 19998, 10001).astype(int)], label='omega heun')
plt.plot(Hts, Hys[np.linspace(1, 19999, 10001).astype(int)], label='theta heun')

#Plotting Runga-Kutta 4 Method
plt.plot(RKts, RKys[np.linspace(0, 19998, 10001).astype(int)], label='omega rk4')
plt.plot(RKts, RKys[np.linspace(1, 19999, 10001).astype(int)], label='theta rk4')

plt.plot(t, sol[:, 0], label='theta Odeint')
plt.plot(t, sol[:, 1], label='omega Odeint')

plt.legend(loc = 'lower right', ncol=2)
plt.xlabel('Time')
plt.grid()
#plt.savefig('Odeint_vs_ODEsolver.png')
plt.show()