# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data = np.loadtxt('lightcurve_data.txt',skiprows=18, usecols=(0,1))
HJD = data[:,0] #Heliocentric Julian Date
DN_Flux = data[:,1] #detrended normalized flux in Kepler bandpass


#folded lightcurve
Period = 3.5485 #days

folded_time = np.zeros_like(HJD)
t_0 = HJD[0]

for i in range(len(HJD)):
    folded_time[i] = (HJD[i] - t_0) % Period


def model(folded_time, tau, t_ref, deltaI):
    '''
    Model to fit the data
    --------------------
    parameter:
        folded_time: the time from the folded light curve
        tau: width of the transit
        t_ref: time at the start of the dip in the lightcurve
        deltaI: difference in the intensity of the dip
        
    return:
        model: the model that should fit the curve
    '''
    model = []
    for i in range(len(folded_time)):
        if folded_time[i] < t_ref:
            model.append(1)
        elif folded_time[i] > (t_ref + tau):
            model.append(1)
        else:
            model.append(1-deltaI)
    return model

def radius(deltaI,I,Radstar):
    '''
    calculates the planets radius
    ---------------
    parameters:
        deltaI: difference in the intensity of the dip
        I: inital intensity (seems to be 1)
        Radstar: radius of the star (given)

    '''
    Radplanet = np.sqrt((deltaI/I)*Radstar**2)
    return Radplanet

def mcmc_alg(iterations, guesses, Period, time, Flux):
    '''
    Performs The Metropolis-Hastings Algorithm
    -----------
    Parameters:
        iterations: number of time to run through the algorithtm
        guesses: list of initial guesses
        Period: the period of hte transit(given)
        time: the folded time
        Flux: detrended normalized flux in Kepler bandpass
    returns:
        best parameters for the model to fit
    '''
    #list of parameters starting with the guess and ending with the best
    tau = [guesses[0]]
    t_ref = [guesses[1]]
    deltaI = [guesses[2]]
    for i in range(iterations):
        #current param is the last param in list until best value is found
        Taucurr = tau[-1]
        trefcurr = t_ref[-1]
        deltaIcurr = deltaI[-1]
        #initital likelihood with the first guesses
        model_init = model(time,Taucurr,trefcurr,deltaIcurr)
        for w in range(len(Flux)):
            Lik_init = (1/np.sqrt(2*np.pi))*((-1/2)*(Flux[i]-model_init[i])**2)
        
        #random change of guess to find optimal value
        rand = np.random.uniform(-0.0001,0.0001,1)
        newTau = Taucurr + rand
        newtref = trefcurr + rand
        newdelI = deltaIcurr + rand
        
        #finding best Tau
        model_1 = model(time,newTau,trefcurr,deltaIcurr)
        for x in range(len(Flux)):
            Lik_1 = (1/np.sqrt(2*np.pi))*((-1/2)*(Flux[i]-model_1[i])**2)
        #Metropolis ratio
        r1 =  Lik_1/Lik_init
        if r1 >= 1:
            tau.append(newTau)
        else:
            #compair the ration to a random variable between 0 and 1
            U = np.random.uniform(0,1,1)
            if U <= r1:
                tau.append(newTau)
        
        #finding best t_ref
        model_2 = model(time,Taucurr,newtref,deltaIcurr)
        for y in range(len(Flux)):
            Lik_2 = (1/np.sqrt(2*np.pi))*((-1/2)*(Flux[i]-model_2[i])**2)
        #Metropolis ratio
        r2 =  Lik_2/Lik_init
        if r2 >= 1:
            t_ref.append(newtref)
        else:
            #compair the ration to a random variable between 0 and 1
            U = np.random.uniform(0,1,1)
            if U <= r2:
                t_ref.append(newtref)
        
        #finding best DeltaI
        model_3 = model(time,Taucurr,trefcurr,newdelI)
        for z in range(len(Flux)):
            Lik_3 = (1/np.sqrt(2*np.pi))*((-1/2)*(Flux[i]-model_3[i])**2)
        #Metropolis ratio
        r3 =  Lik_3/Lik_init
        if r3 >= 1:
            deltaI.append(newdelI)
        else:
            #compair the ration to a random variable between 0 and 1
            U = np.random.uniform(0,1,1)
            if U <= r3:
                deltaI.append(newdelI)
    return(tau, t_ref, deltaI)
        

#initialguesses
guess_tau = 0.2
guess_tref = 2.26
guess_deltaI = 0.007
guesses = np.array([guess_tau, guess_tref, guess_deltaI])

optTau, optTref, optDeltaI = mcmc_alg(1000,guesses,Period,folded_time,DN_Flux)

I = 1.0
star_radius = 1.79 #solar masses
planet_radius = radius(optDeltaI[-1], I, star_radius)

print(planet_radius)
#plt.plot(folded_time, DN_Flux, '.')
##plt.plot(sorted(folded_time), model(sorted(folded_time), guess_tau, guess_tref, guess_deltaI))
#plt.plot(sorted(folded_time), model(sorted(folded_time), optTau[-1], optTref[-1], optDeltaI[-1]))
##Zoom in
#plt.ylim(0.990, 1.005)
#plt.xlim(2.2, 2.55)
#plt.xlabel('Time HJD')
#plt.ylabel('Folded Relative Intensity')
#plt.title('Kepler Folded Lightcurve data')
#plt.savefig('zoombestfitLightcurve1.png')
#
#plt.hist(optDeltaI, alpha=0.5, bins=50)
#plt.xlabel('Relative Intensity')
#plt.ylabel('Frequency')
#plt.title('Distribution of values of Intensity Difference')
#plt.savefig('DeltaIhistogram1.png')