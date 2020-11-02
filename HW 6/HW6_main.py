# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 09:40:02 2020

@author: Sarah
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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

def LeastSquare(X, Y):
    '''
    least square will calculate the values for the coefficients (eq 18)
    by calculating the equation:
        Theta = (X.transpose * X)^-1 * X.transpose*Y
    --------------
    X: Design matrix
    Y: matrix of initial observations
    --------------
    returns:
        Theta: parameter vector, the best fit parameters for the matrix
    '''
    Xdag = np.transpose(X)
    XdagX = np.matmul(Xdag, X)
    invXdagX = (XdagX)**-1
    invXdagXdagX = np.matmul(invXdagX, Xdag)
    theta = np.matmul(invXdagXdagX, Y)

    return (theta)

def Uncertainties(X):
    '''
    uncertainties of the coefficients (eq 33)
    sigma^2 = ((X.transpose * X)^-1) along the diagnal
    -------------
    Parameter:
       X: Design matrix
    ------------
    return:
        theta_sig: the uncertainty of each variable
    '''
    Xdag = np.transpose(X)
    XdagX = np.matmul(Xdag, X)
    invXdagX = (XdagX)**-1
    diagnal = invXdagX.diagonal()
    theta_sig = np.sqrt(diagnal)
    
    return (theta_sig)
    
'''
problem 1 and 2
'''
#Design Matrix
#Colums: 1, log10(Period), Z
array = np.zeros((Period.shape[0], 3))
array[:, 0] = 1 # column 1
array[:, 1] = np.log10(Period) # column 2
array[:, 2] = Z # column 3
#array to matrix
X = np.matrix(array)

#absolute magnitude
A_lam = R_v * Excess * A_k_lam
M_K = K_mag - 5*np.log10(Distance*1000) + 5 - A_lam

#calculate the K parameters
Theta_K = LeastSquare(X, M_K)
trans_K = np.transpose(Theta_K)
Theta_K_sig = Uncertainties(X)

#equation 3 in homework
model = trans_K[0] + trans_K[1]*np.log10(np.sort(Period)) + trans_K[2]*np.sort(Z)

#printing the estimated values for Alpha, Beta, and Gamma
#print(Theta_K)#uncomment for problem 1
#print(Theta_K_sig)

#print plots (problem 2)
#plt.plot(np.log10(Period), M_K,'o',markersize=4.0)
#plt.plot(np.log10(np.sort(Period)), np.transpose(model), color='r')
#plt.ylim([5,-15])
#plt.xlabel('log(Period) [days]')
#plt.ylabel('Absolute Magnitude')
#plt.title('The Cepheid Period-Luminosity-Metallicty Relation \n K-band, No Error')

'''
Problem 3
'''
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
    invXdagVX = (XdagVX)**-1
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
    invXdagVX = (XdagVX)**-1
    diagnal = invXdagVX.diagonal()
    theta_sig = np.sqrt(diagnal)
    
    return (theta_sig)

#identity matrix of sig**2=0.1**2   
V = np.identity(M_K.shape[0])*(0.1**2)

weightedTheta= Weighted(X,V,M_K)
transweightedTheta = np.transpose(weightedTheta)
weighted_Theta_k_sig = WeightedUncertainties(X,V)

#equation 3 in homework
model2 = transweightedTheta[0] + transweightedTheta[1]*np.log10(np.sort(Period)) + transweightedTheta[2]*np.sort(Z)

#printing the estimated values for Alpha, Beta, and Gamma
print(weightedTheta)
print(weighted_Theta_k_sig)

#print plots (problem 3)
plt.plot(np.log10(Period), M_K,'o',markersize=4.0)
plt.errorbar(np.log10(Period), M_K, yerr=0.1, fmt=' ')
plt.plot(np.log10(np.sort(Period)), np.transpose(model2), color='r')
plt.ylim([5,-15])
plt.xlabel('log(Period) [days]')
plt.ylabel('Absolute Magnitude')
plt.title('The Cepheid Period-Luminosity-Metallicty Relation\n K-band, Apparent Magnitude error: 0.1')
#plt.savefig('HW6CepheidAppMagErr.png')


