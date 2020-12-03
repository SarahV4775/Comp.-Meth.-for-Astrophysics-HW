# -*- coding: utf-8 -*-
"""
Non-Linear Least Squares

Use a Non-Linear least Squares method, Levenberg-Marquardt, to fit the curve of 
the lightcurve in UM234_lightcurve.txt

Sarah Vaughn
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def sine_function(param, mag,time):
    """
    General Sine function to fit to curve data.
    -----------
    Parameters:
        param : list of parameters 
            [Amplitude, Frequency, Phase, Mean]
        mag : the magnitude array from the light curve data
        time : time array
    """
    sf = param[0]*np.sin(param[1]*time+param[2]) + param[3] - mag
    return(sf)
    
def LM(X, param, nmax, step, uncertainty, time):
    """
    Levenberg-Marquardt Method to solve the non-linear least squares problem
    of finding the best fit parameters for the model that is fitted to the 
    light curve
    
    Simplified Equation for the Levenberg-Marquardt Method:
    Gradient_deltheta sf(theta) = Gradient_theta sf(theta) + K*deltheta
    solving for deltheta, then using:
    theta_best = deltheta + theta
    ------
    Parameters:
        X : arguments from the lightcurve
        param : list of parameters 
            [Amplitude, Frequency, Phase, Mean]
        nmax : max number of iterations
        step : the step between the guess param to new param
        uncertainty = uncertainty in data
        time : time array
    ------
    Return:
        the best values for the parameters
    """
    #lsit of all the parameters and thier reduced chi squared value
    #best value is the one thats closest to 1
    ParamList = pd.DataFrame([], colums = ["Rchisquare","theta[0]","theta[1]","theta[2]","theta[3]"])
    n = 0
    while n < nmax:
        
        theta_current = param[:] # Inital guess parameters
        deltheta = theta_current + step # incrimented parameter
        G_dt = np.gradient(sine_function(deltheta, X, time)) #gradient with respect to delta theta
        G_t = np.gradient(sine_function(theta_current, X, time)) #gradient with respect to theta
        # curvature matrix from notes 15 eq(8)
        #K = ?
        
        #solve for delta theta
        DeltaTheta = (G_dt - G_t) * K**-1
        
        theta_best = DeltaTheta + theta_current
        
        #check how close the current value is 
        #Chi Squared = (data-model/uncertainty)^2 
        ChiSquare = np.sum(((X-sine_function(theta_best, X, time))/uncertainty)**2)
        #Reduced Chi Square 
        numParam = len(param) #number of paramters in theta             
        # Number of degrees of freedom (N) = data size - number of parameters  
        N = X.size- numParam
        RChiSq = ChiSquare/N
        
        #update list
        NewChiandTheta = [RChiSq,theta_best[0],theta_best[1],theta_best[2],theta_best[3]]
        lengthList = len(ParamList)
        ParamList.loc[lengthList] = NewChiandTheta
        
        #update the parameters
        theta_current = theta_best
        n =+ 1
    
    #the index of the the Chi Squared value that is the closest to 1
    Bestindex = min(range(len(ParamList)), key=lambda i: abs(ParamList[i]-1))
    # the best param
    BestParam = ParamList.iloc[Bestindex,1:4] #[row,col]
    return(BestParam)
    

#read in the lightcurve data
data = pd.read_csv('UM234_lightcurve.txt',sep='\t',skiprows=0, header=None,)
Mag = data.iloc[:,2]
MJD = data.iloc[:,6]
magerr = data.iloc[:,3]

#guess Parameters
guess_mean = np.mean(Mag)
guess_phase = 1
guess_freq = (2*np.pi)/1800
guess_amp = 1
guess = [guess_amp, guess_freq, guess_phase, guess_mean]

#Least Squares to get the best parameters
# number of iteration = 500
# step = 0.01
bestparam = LM(Mag, guess, 500, 0.01, magerr)

#sine function with the bestparam values
data_fit = sine_function(bestparam, Mag, MJD)

#plot of the data
plt.plot(MJD, Mag,'o',markersize=4.0)
plt.errorbar(MJD,Mag,yerr=magerr,fmt=' ')

#fitted plot 
plt.plot(MJD, data_fit, 'g', label='after fitting')

#Y axis limits
plt.ylim([max(Mag),min(Mag)])

# Give x axis label for the sine wave plot
plt.xlabel('Time[days]')

# Give y axis label for the sine wave plot
plt.ylabel('Magnitude')

plt.grid(True, which='both')
plt.axhline(y=0, color='k')
plt.title("Curve Fit of UM 234")
plt.legend()