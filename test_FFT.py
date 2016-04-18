import modules.uFM as uFM
uFM.all()

import numpy as np
import sys
import sympy as sym
import scipy as scy
import matplotlib.pyplot as pl
import os
import time


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

'''
#---------------------------------------
# 3D FFT: (for intDER)
#---------------------------------------
from numpy.fft import fftn, ifftn

bin1, bin2 = 0, 0

# 3D FFT:
N_k = 8000 # for every dimension, so 8.000.000 in total
vect_k = np.linspace(0.0001,0.5,N_k)

samples=np.empty(N_k)
for i in range(N_k):
    samples[i] = uFM.wrapper_derW(vect_k[i],bin1,bin2)

startFFT=time.clock()
sample_derK = uFM.W2_x_derW1_FFT_3D(vect_k,bin1,bin2)
[kx, ky, kz] = np.meshgrid(vect_k,vect_k,vect_k)
mod_k = np.sqrt(kx*kx + ky*ky + kz*kz) #Can be optimised, because it's sym.
sample_P = uFM.zero_spectrum_FFT_3D(mod_k)

# sample_derK, sample_P = np.empty((N_k,N_k,N_k)), np.empty((N_k,N_k,N_k))
# # This is really shitty python:
# for i in range(N_k):
#     for j in range(N_k):
#         for k in range(N_k):
#             modul_k = mod_k[i,j,k]
#             sample_derK[i,j,k] = uFM.W2_x_derW1_FFT(modul_k,bin1,bin2)*kz[i,j,k]/modul_k
#             sample_P[i,j,k] = uFM.zero_spectrum_py(modul_k)
derK_FFT, P_FFT = fftn(sample_derK), fftn(sample_P)
print "3D FFT: %g sec" %(time.clock()-startFFT)

# 3D inverse FFT:
bohx =  (0.5/N_k)**3 * ifftn(np.multiply(derK_FFT,P_FFT))
bohx_squared = bohx[0,:,0]

fig2=pl.figure()
ax1=fig2.add_subplot(111)
ax1.plot(vect_k[-200:-2],samples[-200:-2],'r-',label="derW")
# ax1.plot(vect_k[:200],naive_convol[:200],'b-',label="Naivily convolved")
# ax1.plot(vect_k,bohx_squared,'g-',label="Shitty convolved")
ax1.grid(True)
ax1.legend(loc='best')
# ax1.set_yscale('log')
# ax1.set_xlabel("$k$ [$h$/Mpc]")
# ax1.set_ylabel("$P(x)$ [(Mpc/$h$)$^3$]")
fig2.savefig('plots/is_FFT_so_magic.pdf')


quit()
'''

#########################################
# Int1 solved with radial FFT:
#########################################

bin1, bin2 = 0, 2
N_k = 10000
R = 1.4
vect_k = np.linspace(0.0001,R,N_k)
K_samples, P_samples = np.empty(N_k), np.empty(N_k)
for i in range(N_k):
    K_samples[i] = uFM.K_py(vect_k[i],bin1,bin2)
    P_samples[i] = uFM.zero_spectrum_py(vect_k[i])
    # P_samples[i] = uFM.CLASS_der(vect_k[i],1)

# # Naive:
# K_fft, P_fft = rfft(K_samples), rfft(P_samples)
# convol_fft = K_fft*P_fft
# naive_convol = 0.5/5000 * irfft(convol_fft)

# Radial:
start = time.clock()
convolved_spectrum = np.sqrt(uFM.vol_shell_original_py(bin1)*uFM.vol_shell_original_py(bin2))/(2*np.pi)**3 *radial_convolution(P_samples,K_samples,R)
stop = time.clock()
print stop-start

# log_P = np.log(P_samples)
# log_convolved = np.empty(convolved_spectrum.shape[0])
# for i in range(convolved_spectrum.shape[0]):
#     if convolved_spectrum[i]>=0:
#         log_convolved[i]=np.log(convolved_spectrum[i])
#     else:
#         log_convolved[i]=-np.log(-convolved_spectrum[i])

fig2=pl.figure()
ax1=fig2.add_subplot(111)
ax1.plot(vect_k[:-250],P_samples[:-250],'r-',label="Spectrum")
# ax1.plot(vect_k[:200],naive_convol[:200],'b-',label="Naivily convolved")
ax1.plot(vect_k[:-250],abs(convolved_spectrum[:-250]),'g-',label="Radially convolved")
ax1.grid(True)
ax1.legend(loc='best')
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_xlabel("$k$ [$h$/Mpc]")
ax1.set_ylabel("$P(x)$ [(Mpc/$h$)$^3$]")
ax1.set_xlim(1e-4,1.8)
fig2.savefig('plots/spectrum_%d%d.pdf'%(bin1,bin2))

'''

###################################################
# Try to compute radial FT manually and compare:
###################################################

def x(t):
    alpha = 4.
    return np.exp(-(t-alpha)**2 )
def g(t):
    if (t>1. and t<2.):
        return 1.
    else:
        return 0.

from scipy.integrate import quad
def Fx_manual(om):
    return 4*np.pi/om * quad(lambda r:x(r)*r*np.sin(r*om),0,100 )[0]

R, N_points = 20., 200
samples = np.linspace(0.,R,N_points)
x_sam = x(samples)

frequencies = fftfreq(N_points,R)*2*np.pi

Fx_man, g_sam = np.empty(N_points), np.empty(N_points)
for i in range(N_points):
    # Fx_man[i]=Fx_manual(frequencies[i])
    g_sam[i]=g(samples[i])

# Fx = radial_fft(x_sam,R)
# Fg = radial_fft(g_sam,R)

# convol_radial = radial_ifft(Fx * Fg, R)

fig2=pl.figure()
ax1=fig2.add_subplot(111)
ax1.plot(samples,x_sam,'r-',label="x(t)")
ax1.plot(samples,radial_convolution(x_sam,g_sam,R),'b-',label="Convolution")
# ax1.plot(frequencies,Fx,'b-',label="F. Transform")
# ax1.plot(frequencies,Fx_man,'g-',label="Manual")
ax1.grid(True)
ax1.legend(loc='best')
# ax1.set_yscale('log')
# ax1.set_xlabel("$k$ [$h$/Mpc]")
# ax1.set_ylabel("$P(x)$ [(Mpc/$h$)$^3$]")
fig2.savefig('plots/test_convol_radial.pdf')


# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# # ########################
# # Test FFT normale:
# # ########################
# #
# def x(t):
#     beta = 1.
#     return np.exp(-beta*abs(t))
# def h(t):
#     delta, A = 10., 2.
#     if (t>0. and t<delta):
#         return A
#     else:
#         return 0.


# N_k = 5000
# samples = np.concatenate((np.linspace(0.,30.,N_k/2),np.linspace(-30.,0.,N_k/2)))
# x_sam, h_sam = np.empty(N_k), np.empty(N_k)
# for i in range(N_k):
#     x_sam[i], h_sam[i] = x(samples[i]), h(samples[i])

# x_FFT, h_FFT = rfft(x_sam), rfft(h_sam)
# convol = 60./N_k * irfft(x_FFT*h_FFT)

# fig2=pl.figure()
# ax1=fig2.add_subplot(111)
# ax1.plot(samples,x_sam,'r-',label="x")
# ax1.plot(samples,h_sam,'b-',label="h")
# ax1.plot(samples,convol,'g-',label="Convol")
# ax1.grid(True)
# ax1.legend(loc='best')
# # ax1.set_yscale('log')
# # ax1.set_xlabel("$k$ [$h$/Mpc]")
# # ax1.set_ylabel("$P(x)$ [(Mpc/$h$)$^3$]")
# fig2.savefig('plots/test_convol.pdf')








'''








