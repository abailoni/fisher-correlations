import numpy as np

#---------------------------------------
# RADIAL CONVOLUTION using FFT:
#---------------------------------------

from numpy.fft import rfft, irfft, fftfreq,fftshift
from numpy.fft import fft, ifft

def sine_fft(array):
    N = array.shape[0]
    array_mod = np.concatenate( (array, np.array([0.]), -array[::-1][:-1] ) )
    return - np.imag(fft(array_mod))[:N]

def radial_fft(c_n, R):
    N = c_n.shape[0]
    indices = np.arange(N)
    c_k =  R/(np.pi*indices[1:]) * sine_fft(c_n * indices*R / N)[1:]
    c_k = np.insert(c_k, 0, 2* np.sum(indices*indices * (R/N)*(R/N) * c_n) )
    # return  c_k
    return  2*np.pi*R/N* c_k

def sine_ifft(array):
    N = array.shape[0]
    array_mod = np.concatenate( (array, np.array([0.]), -array[::-1][:-1] ) )
    return - 1./(2.*N)*np.imag(fft(array_mod))[:N]

def radial_ifft(c_k, R):
    N = c_k.shape[0]
    indices = np.arange(N)
    c_n =  N/(indices[1:]*R) * sine_ifft(c_k * np.pi*indices / R)[1:]
    c_n = np.insert(c_n, 0, 1./N * np.sum( indices*indices * (np.pi/R)*(np.pi/R) * c_k ) )
    # return c_n
    return N/(2*np.pi*R) * c_n

def radial_convolution(array1,array2,R):
    F1 = radial_fft(array1,R)
    F2 = radial_fft(array2,R)
    return radial_ifft(F1 * F2, R)
