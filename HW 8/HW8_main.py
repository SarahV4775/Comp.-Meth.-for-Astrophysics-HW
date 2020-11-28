# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 15:04:22 2020

@author: Sarah
"""

import numpy as np
import matplotlib.pyplot as plt

strain = np.load('strain.npy')
time_min = np.linspace(1,len(strain) + 1, len(strain))
time_hour = time_min/60
time_days = time_hour/24
time_sec = time_min*60 #for later plotting FFT

def DFT(data):
    '''
    1 dimentional Discrete Fourier Transform using matrix multiplication
    --------
    Parameter:
    data = data that is pulled in (the strain)
    -------
    return:
        H = dot product of W_kn matrix and the h-array
    '''
    h = np.array(data)
    N = h.shape[0]
    n = np.arange(N)
    k = n.reshape((N, 1))
    W_kn = np.exp((-2j * np.pi * k * n)/N)
    H = np.dot(W_kn, h)
    return(H)
    
    
def FFT(data):
    '''
    Cooley-Tukey algorithm to perform the Fast fourier transform
    -------
    Parameter:
    data = data that is pulled in (the strain)
    '''
    h = np.array(data)
    N = h.shape[0]
    
    if N % 2 > 0:
        raise ValueError("data size must be a power of 2")
    elif N <= 32: #cutoff
        return DFT(data)
    else:
        #break into odd and even
        even = FFT(h[::2]) #start at 1st and jump in interval of 2 (even)
        odd = FFT(h[1::2]) #start at 2nd and jump in interval of 2 (odd)
        k = np.arange(N)
        wk = np.exp((-2j * np.pi * k)/ N)
        #join the arrays with the two different factors (w)
        #wk[:N//2] : slice of first half of the list
        #wk[N//2:] : slice of last half of the list
        H_k= np.concatenate([even + wk[:N//2] * odd, even + wk[N//2:] * odd]) 
        return H_k

#FFT of the strain data
fourier_transform = FFT(strain)
real_tran = np.real(fourier_transform)*2 #real values
N = len(time_sec)
#print(real_tran)
#freq from time
time_samp =np.mean(np.diff(time_sec))
freq_samp = 1/(time_samp)
Nyfreq = freq_samp/2 #Nyquist Frequency
fstep = Nyfreq/len(time_sec)
freq = np.array([i*fstep for i in range(len(time_sec))])

#Values of the spike in the data
#frequency
freq_range = freq[(freq > 10**(-3.5)) & (freq < 10**(-2.5))]
peak_freq = np.max(freq_range)
#Amplitude
amp_range = real_tran[(freq > 10**(-3.5)) & (freq < 10**(-2.5))]
amp_at_peak = np.max(amp_range) #np.interp(peak_freq, freq, fourier_transform)
print('Amplitude: ' + str(amp_at_peak))
print('Frequency: ' + str(peak_freq))

plt.plot(np.log10(freq[:N//2]), np.log10(fourier_transform[:N//2]))
plt.xlabel('log(frequency [Hz])')
plt.ylabel('log(Amplitude Spectrum)')
plt.title('Log-Log Amplitude Spectrum vs Frequency')
#plt.savefig('loglogAmpplot.png')

# Calculations for the total mass M and the seperation R
h = amp_at_peak
f_gw = peak_freq
con1 = 2.6*10**(-21) # constants G/4pi^2c^4
con2 = 10**(-4) # constants G^1/2/2pi
D = 12 #pc
R = ((D*h*con2)/(con1*f_gw))**(1/2)
M = ((f_gw*(R**(3/2)))/con2)**2

print("total mass: " + str(M)) # Solar Mass
print("Separation: " + str(R)) # Solar Radius