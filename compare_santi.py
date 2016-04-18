import sys
import numpy as np
import sympy as sym
import scipy as scy
import matplotlib.pyplot as pl
import os
import time

SANTI_FOLDER = "santiRESULTS/"

FM_NOAP = np.loadtxt(open(SANTI_FOLDER+"FisherMatrix--no-AP14bins.txt","rb"),delimiter="\t",skiprows=0)
FM_YESAP = np.loadtxt(open(SANTI_FOLDER+"FisherMatrix--yes-AP14bins.txt","rb"),delimiter="\t",skiprows=0)

own_FM = np.loadtxt(open("OUTPUT/FM_uncorr_K0_kmax0.2.csv","rb"),delimiter=" ",skiprows=0)
own_FM = np.delete(np.delete(own_FM,5,0),5,1)

percNO = 100*own_FM/FM_NOAP-100
percYES = 100*own_FM/FM_YESAP-100
bin1_N, bin1_Y = percNO[0,:], percYES[0,:]

np.savetxt(SANTI_FOLDER+"comparisonNO.csv", percNO)
np.savetxt(SANTI_FOLDER+"comparisonYES.csv", percYES)
np.savetxt(SANTI_FOLDER+"bin1_compare.csv", np.column_stack(bin1_Y,bin1_N))
