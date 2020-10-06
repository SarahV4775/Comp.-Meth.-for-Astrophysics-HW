# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 03:50:41 2020

@author: Sarah
"""

import matplotlib.pyplot as plt
import numpy as np
from ODE_solver import RK4


#hydrostatic equations
def hydrostatic(rho):
    def equations(X, r):
        G = 6.67e-12    #cm^3 g^-1 s^-2
        M = X
        dMdr = 4*np.pi*(r**2)*rho
        dPdr = -(G*M*rho)/(r**2)
        dydr = [dMdr, dPdr]
        return dydr
    return equations

#Relativistic Hydrostatic equilibrium equations
def rel_hydrostatic(rho):
    def rel_equations(X, r):
        G = 6.67e-12    #cm^3 g^-1 s^-2
        c = 3*10e10    #cm s^-1
        M, P = X
        dMdr = 4*np.pi*(r**2)*rho
        dPdr = -(G*M/(r**2))*rho*(1+(P/(rho*(c**2))))*(1+((4*np.pi*(r**3)*P)/(M*(c**2))))*(1-(2*G*M/(r*(c**2))))**-1
        dydr = [dMdr, dPdr]
        return dydr
    return rel_equations

rho_c_wd = np.linspace(1e4, 1e6, 3)#white dwarf central densities 
rho_c_ns = np.linspace(1e14, 1e16, 3)#Neutron star central densities
h = 1000

#mass-radius curve for White dwarf 
for i in range(len(rho_c_wd)):
    r = np.linspace(1, 3e9, 5) #cm
    P_c = (1e13)*(rho_c_wd[i]/2)**(5/3) #eq 4
    y0 = np.array([1e-4, P_c])
    Msol, Rsol = RK4(hydrostatic(rho_c_wd[i]), y0, r, h)
    plt.plot(r, Msol)
plt.yscale('log')
plt.xlabel('Radius (cm)')
plt.ylabel('Mass (g)')
plt.title('Mass-Radius Curve for White Dwarf ')


#mass-radius curve for Neutron Stars 
for i in range(len(rho_c_ns)):
    r = np.linspace(1, 3e9, 5) #cm
    P_c = (5.4e9)*(rho_c_ns[i])**(5/3) #eq 5
    y0 = np.array([1e-4, P_c])
    Msol, Rsol = RK4(rel_hydrostatic(rho_c_ns[i]), y0, r, h)
    plt.plot(r, Msol)
plt.yscale('log')
plt.xlabel('Radius (cm)')
plt.ylabel('Mass (g)')
plt.title('Mass-Radius Curve for Neutron Star')


#calculated Mass of J0030+0451 where R = 1.302e6 [cm] using TOV equations
#rho values == Neutron Star rhos
for i in range(len(rho_c_ns)):
    r = np.linspace(1, 1.302e6, 5) #cm
    P_c = (5.4e9)*(rho_c_ns[i])**(5/3) #eq 5
    y0 = np.array([1e-4, P_c])
    Msol, Rsol = RK4(rel_hydrostatic(rho_c_ns[i]), y0, r, h)
    print("Mass for J0030+0451 at radius 13.02 [km] and centrl density ",
          rho_c_ns[i], " is: ", Msol[i], "Solar Masses")
    