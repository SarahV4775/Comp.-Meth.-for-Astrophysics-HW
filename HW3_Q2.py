# -*- coding: utf-8 -*-
"""
Created on Sat Sep 23 18:08:38 2020

@author: Sarah
"""

import numpy as np
import matplotlib.pyplot as plt
from Matrix_class import Matrix as M

#read in the data from the text file
l, u, As = np.loadtxt('A_coefficients.dat', unpack=True, delimiter=",",
                          dtype={"names": ('l', 'u', 'A_ul'),
                                 "formats": (np.int, np.int, np.float)})

#variables
temp = np.logspace(1, 10, 10)
k = 1.3807e-16 # Boltzmann constant [cm^2*g*s^-2*K^-1]
h = 6.626e-27 # Plancks constant [cm^2*g*s^-1]
c = 3e10 # speed of light [cm/s]
nu = [] # frequencies
for i in range(0,9):  #matrix
    for j in range(0,9):
        if i >= j:
            nu= np.append(nu, 0) 
        else:
            dE = -13.6*((1/((j+1)**2))-(1/((i+1)**2)))
            v = h/dE
            nu= np.append(nu, c/(v*10e-7))


#new matrices initialized to be zeros
A_ul = M(9, 9, 0)
B_lu = M(8, 8, 0)
B_ul = M(8, 8, 0)
J = M(8, 8, 0)

#fill the matrices with the values according to the text file
# A_ul values in text file, B_lu, B_ul and J found through equations relating to A_ul
for i in range(As.shape[0]):
	A_ul.rows[l[i]-1, u[i]-1] = As[i]
A_ul.rows = A_ul.rows[:-1, 1:]

#for the value in the temp range calculate the rest of the matices
for t in temp:
    
    for row in range(A_ul.rows.shape[0]):
        for col in range(A_ul.rows.shape[1]):
            J.rows[row, col] = ((2*h*nu**3)/c**2)(1/( np.exp((h*nu)/(k*t))-1))
            B_ul.rows[row, col] = A_ul.rows[row, col]*((c**2)/(2*h*nu**3))
            B_lu.rows[col, row] = B_ul.rows[row, col]*((col)**2/(row)**2)
    
    #n-state system equations in matrix form
    # system_mat * n = zero_matrix    
    system_mat = M(8,8) #new empty matrix 
    
    for row in range(system_mat.width):
        for col in range(system_mat.height):
            # u = col, l = row
            if row > col:
                system_mat.rows[row, col] = -(A_ul.rows[col, row] + M.__mult__(B_ul.rows[col, row], J.rows[row, col]))
            if row < col:
                system_mat.rows[row, col] = -(M.__mult__(B_ul.rows[col, row], J.rows[row, col]))
            if row == col:
                system_mat.rows[row, col] = (np.sum(B_ul.rows[col, :row != row]))*J.rows[row, col] + (np.sum(A_ul.rows[col, :row != row]))
    
    # Ax = b => x = b*A^-1
    #where A is the system_mat and b is and 1x8 matrix of zeros solve for x
    b = M.__init__(1,8,0)
    A_invert = M.__invert__(system_mat)
    X = M.__mult__(b, A_invert)
    
    N = np.sum(X.rows) # should equal 1
    
print(t, X[0], label= '0')
print(t, X[1], label= '1')
print(t, X[2], label= '2')
print(t, X[3], label= '3')
print(t, X[4], label= '4')
print(t, X[5], label= '5')
print(t, X[6], label= '6')
print(t, X[7], label= '7')

plt.legend(loc='best')
plt.xlabel('Time')
plt.grid()
plt.show()

    
    