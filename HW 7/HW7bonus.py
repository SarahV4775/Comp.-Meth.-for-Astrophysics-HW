# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data = np.loadtxt('RV_data.txt',skiprows=17, usecols=(0,1,2,3,4))
BJD = data[:,0] #Barycentric Julian Date
Rad_velocity = data[:,1] #Radial velocity
RV_uncertainty = data[:,2] #Uncertainty in target radial velocity
Bisector_Span = data[:,3] #bisector span
BS_uncertainty = data[:,4] #uncertainty in bisector span

Period = 3.5485 #days
t_0 = BJD[0]
folded_time = np.zeros_like(BJD)

for i in range(len(BJD)):
    folded_time[i] = (BJD[i] - t_0) % Period
    
freq = (2*np.pi)/Period
phase = 1.5
amp = Rad_velocity[0]

def sine(amp,freq,time,phase):
    '''
    sine model to fit to the points
    -----
    parameters:
        amp: amplitude of sine wave
        freq: frequency of sine wave
        time: list of time of the data
        phase: phase of sine wave
    '''
    model = np.zeros_like(time)
    for i in range(len(time)):
        model[i] = amp*np.sin(freq*time[i] + phase)
    return model

plt.plot(folded_time, Rad_velocity, '.')
#fitted plot 
plt.plot(sorted(folded_time), sine(amp,freq,sorted(folded_time),phase), 'g', label='after fitting')
plt.ylim(-300, 200)
plt.xlabel('Time BJD')
plt.ylabel('Radial Velocity')
plt.title('Keck-HIRES Radial Velocity Data')
plt.savefig('hw7bonusplot.png')