# -*- coding: utf-8 -*-
"""
Created on Sun Nov  1 18:15:14 2020

@author: Sarah
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import f
from scipy.stats import chisquare

def Weighted (X, V, Y):
    '''
    Using Matrix manipulation to solve for the values of Alpha, Beta, Gamma
    while taking into account the uncertainty (eq 37)
        Theta = (X.transpose * V^-1 * X)^-1 * X.transpose*V^-1*Y
    ----------
    Parameters:
        X: Design matrix
        V: Identity matrix of the uncertainty
        Y: matrix of initial observations
    ----------
    returns:
        Theta: parameter vector, the best fit parameters for the matrix
    '''
    Xdag = np.transpose(X)
    Vinv = np.linalg.inv(V)
    XdagV = np.matmul(Xdag, Vinv)
    XdagVX = np.matmul(XdagV, X)
    invXdagVX = (XdagVX)**(-1)
    invXdagXdagVX = np.matmul(invXdagVX, Xdag)
    invXdagXdagVXV = np.matmul(invXdagXdagVX, Vinv)
    theta = np.matmul(invXdagXdagVXV, Y)
    
    return(theta)

def WeightedUncertainties(X,V):
    '''
    uncertainties of the coefficients (eq 38)
    sigma^2 = ((X.transpose* V^-1 * X)^-1) along the diagnal
    -------------
    Parameter:
       X: Design matrix
       V: Identity matrix of the uncertainty
    ------------
    return:
        theta_sig: the weighted uncertainty of each variable
    '''
    Xdag = np.transpose(X)
    Vinv = np.linalg.inv(V)
    XdagV = np.matmul(Xdag, Vinv)
    XdagVX = np.matmul(XdagV, X)
    invXdagVX = (XdagVX)**(-1)
    diagnal = invXdagVX.diagonal()
    theta_sig = np.sqrt(diagnal)
    
    return (theta_sig)


data = pd.read_csv('cepheid_data.txt',skiprows=1, header=None,)
Period = data.iloc[:,1] #days
Distance = data.iloc[:,2] #kpc
K_mag = data.iloc[:,6]
Excess = data.iloc[:,7]
Z = data.iloc[:,8] #Metallicity [Fe/H]

#from Homework pdf
R_v = 3.1
#from table
A_k_lam = 0.117

#Design Matrix for nested model
#Colums: 1, log10(Period)
array_nest = np.zeros((Period.shape[0], 2))
array_nest[:, 0] = 1 # column 1
array_nest[:, 1] = np.log10(Period) # column 2
#array to matrix
X_nest = np.matrix(array_nest)

#Design Matrix for Full Model
#Colums: 1, log10(Period), Z
array_Full = np.zeros((Period.shape[0], 3))
array_Full[:, 0] = 1 # column 1
array_Full[:, 1] = np.log10(Period) # column 2
array_Full[:, 2] = Z # column 3
#array to matrix
X_Full = np.matrix(array_Full)

#absolute magnitude
A_lam = R_v * Excess * A_k_lam
M_K = K_mag - 5*np.log10(Distance*1000) + 5 - A_lam

#identity matrix of sig**2=0.1**2   
V = np.identity(M_K.shape[0])*(0.1**2)

NestweightedTheta= Weighted(X_nest,V,M_K)
NesttransweightedTheta = np.transpose(NestweightedTheta)
Nestweighted_Theta_k_sig = WeightedUncertainties(X_nest,V)

FullweightedTheta= Weighted(X_Full,V,M_K)
FulltransweightedTheta = np.transpose(FullweightedTheta)
Fullweighted_Theta_k_sig = WeightedUncertainties(X_Full,V)


Nested_model = NesttransweightedTheta[0] + NesttransweightedTheta[1]*np.log10(np.sort(Period))
Full_model = FulltransweightedTheta[0] + FulltransweightedTheta[1]*np.log10(np.sort(Period)) + FulltransweightedTheta[2]*np.sort(Z)

#printing the estimated values for Alpha, Beta, and Gamma
print(NestweightedTheta)
print(Nestweighted_Theta_k_sig)

#print plots (problem 3)
plt.plot(np.log10(Period), M_K,'o',markersize=4.0)
plt.errorbar(np.log10(Period), M_K, yerr=0.1, fmt=' ')
plt.plot(np.log10(np.sort(Period)), np.transpose(Nested_model), color='r')
plt.ylim([5,-15])
plt.xlabel('log(Period) [days]')
plt.ylabel('Absolute Magnitude')
plt.title('The Cepheid Period-Luminosity-Metallicty Relation\n Nested Model  \n K-band, Apparent Magnitude error: 0.1')
#plt.savefig('HW6CepheidNestedModel.png')

Nested_Model_Array = []
Full_Model_Array = []
i = 0
for i in range(len(Nested_model.T)):
    Nested_Model_Array = np.append(Nested_Model_Array, Nested_model.T[i])

for i in range(len(Full_model.T)):
    Full_Model_Array = np.append(Full_Model_Array, Full_model.T[i])

Chi_sqr_Full = chisquare(Full_Model_Array, f_exp=M_K).statistic
Chi_sqr_Nested = chisquare(Nested_Model_Array, f_exp=M_K).statistic

full_nu = M_K.shape[0] - 3
nest_nu = M_K.shape[0] - 2

top = (Chi_sqr_Nested - Chi_sqr_Full) / (full_nu - nest_nu)
bottom = (Chi_sqr_Full/full_nu)

# F calculated using equations
F = top/bottom

# F calculated using scipy.stats.f.cdf()
F_cdf = f.cdf(F, np.abs(full_nu - nest_nu), full_nu)

print(F, F_cdf)