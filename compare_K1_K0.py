import sys
import numpy as np
import sympy as sym
import scipy as scy
import matplotlib.pyplot as pl
import os
import time

INPUT_FOLDER = "OUTPUT/FM_7bins/"
OUTPUT_FOLDER = "plots/marginalisation/comparison/7bins/"
name = "K0_kmax0"
cFM0 = np.loadtxt(open(INPUT_FOLDER+"FM_corr_"+name+".csv","rb"),delimiter=" ",skiprows=0)
uFM0 = np.loadtxt(open(INPUT_FOLDER+"FM_uncorr_"+name+".csv","rb"),delimiter=" ",skiprows=0)
name = "K1_kmax0"
cFM1 = np.loadtxt(open(INPUT_FOLDER+"FM_corr_"+name+".csv","rb"),delimiter=" ",skiprows=0)
uFM1 = np.loadtxt(open(INPUT_FOLDER+"FM_uncorr_"+name+".csv","rb"),delimiter=" ",skiprows=0)



# uFM_K1 = np.loadtxt(open("OUTPUT/FM_uncorr_K1.csv","rb"),delimiter=" ",skiprows=0)

# Correction factor volume:
correction_factor = np.pi * 15000. / (180.**2 * 4.)
cFM0, uFM0 = cFM0*correction_factor, uFM0*correction_factor
cFM1, uFM1 = cFM1*correction_factor, uFM1*correction_factor


# Adding k_derivative term: (only non-correlated at the moment)
# uFM = uFM_K1

# ------------------------------------
# Checking dependence on k_max:
# ------------------------------------
# uFM_04 = np.loadtxt(open(INPUT_FOLDER+"../OUTPUT/FM_uncorr_0.45.csv","rb"),delimiter=" ",skiprows=0)
# uFM_02 = np.loadtxt(open(INPUT_FOLDER+"../OUTPUT/FM_uncorr_0.2.csv","rb"),delimiter=" ",skiprows=0)

# difference = cFM / uFM

# print abs((uFM_04-uFM_02) / uFM_02)

# np.savetxt(INPUT_FOLDER+"FM_difference.csv", difference )


# ------------------------------------
# Computing uncertaintes: (finally)
# ------------------------------------

# inv_cFM = np.linalg.inv(cFM)
# inv_uFM = np.linalg.inv(uFM)

n_var = {"h": 0, "\\Omega_b": 1, "n_s": 2, "\\Omega_m": 3, "w_p": 4, "w_1": 5, "\\gamma": 6, "b_0": 7,"b_1": 8,"b_2": 9,"b_3": 10,"b_4": 11,"b_5": 12,"b_6": 13, "b_7": 14}
var_names = ["h", "\\Omega_b", "n_s", "\\Omega_m", "w_p", "w_1", "\\gamma", "b_0","b_1","b_2","b_3","b_4","b_5","b_6", "b_7"]
ref_values = {'\\Omega_m': 0.25, 'h': 0.7, '\\Omega_b': 0.0445, 'n_s': 0.967, "w_p": -0.9, "w_1": 0., "\\gamma": 0.545, "b_0": 1.,"b_1": 1.,"b_2": 1.,"b_3": 1.,"b_4": 1.,"b_5": 1.,"b_6": 1., "b_7": 1.}

# Fully marginalised uncertaintes:
# for i in range(cFM.shape[0]):
#     sigmaC, sigmaU = np.sqrt(inv_cFM[i,i]), np.sqrt(inv_uFM[i,i])
#     print "%s\t\tFM_corr=%g \tFM_uncorr=%g  \t%g%%" %(var_names[i],sigmaC, sigmaU, (sigmaC-sigmaU)/sigmaC*100)


# Plot ellipse for two variables:

def marginalise(inv_matrix,var1,var2):
    temp_matrix = np.row_stack((inv_matrix[var1,:],inv_matrix[var2,:]))
    return np.column_stack((temp_matrix[:,var1], temp_matrix[:,var2]))

# sub_cFM = marginalise(inv_cFM,n_va1,n_var2)

# Return axes and angle:
def sigma1_ellipse(FM_matrix, var1, var2):
    inv_matrix = np.linalg.inv(FM_matrix)

    sub_matrix = marginalise(inv_matrix,var1,var2)
    # Find axes:
    term1 = (sub_matrix[0,0]+sub_matrix[1,1])/2.
    term2 = np.sqrt((sub_matrix[0,0]-sub_matrix[1,1])**2/4. +sub_matrix[0,1]**2)
    a = np.sqrt( term1+term2 )
    b = np.sqrt( term1-term2 )
    # eigenvalues = np.linalg.eig(inv_matrix)[0]

    angle = np.arctan( 2*sub_matrix[0,1]/(sub_matrix[0,0]-sub_matrix[1,1]) ) / 2.
    return 2*1.52*a, 2*1.52*b, angle



from math import degrees
import matplotlib.patches as patch

def compare(uFM0,cFM0,uFM1,cFM1,var1,var2): #var names!
    res_uFM0 = sigma1_ellipse(uFM0,n_var[var1],n_var[var2])
    res_cFM0 = sigma1_ellipse(cFM0,n_var[var1],n_var[var2])
    res_uFM1 = sigma1_ellipse(uFM1,n_var[var1],n_var[var2])
    res_cFM1 = sigma1_ellipse(cFM1,n_var[var1],n_var[var2])
    print res_uFM0[0], res_uFM0[1]
    print degrees(res_uFM0[2]), degrees(res_cFM0[2])
    ellipse1 = patch.Ellipse(xy=(ref_values[var1], ref_values[var2]), width=res_uFM0[0], height=res_uFM0[1], angle=degrees(res_uFM0[2]), alpha= 0.7, color="r", label="K0 w/o corr")
    ellipse2 = patch.Ellipse(xy=(ref_values[var1], ref_values[var2]), width=res_cFM0[0], height=res_cFM0[1], angle=degrees(res_cFM0[2]), alpha= 0.7, color="b", label="K0 with corr")
    ellipse3 = patch.Ellipse(xy=(ref_values[var1], ref_values[var2]), width=res_uFM1[0], height=res_uFM1[1], angle=degrees(res_uFM1[2]), alpha= 0.7, color="g", label="K1 w/o corr")
    ellipse4 = patch.Ellipse(xy=(ref_values[var1], ref_values[var2]), width=res_cFM1[0], height=res_cFM1[1], angle=degrees(res_cFM1[2]), alpha= 0.7, color="y", label="K1 with corr")
    fig=pl.figure()
    ax1=fig.add_subplot(111)
    ax1.add_patch(ellipse2)
    ax1.add_patch(ellipse1)
    ax1.add_patch(ellipse4)
    ax1.add_patch(ellipse3)

    # ax1.set_aspect('equal')
    # Set limits plot:
    if abs(res_cFM0[0]/2*np.cos(res_cFM0[2]))>abs(res_cFM0[1]/2*np.sin(res_cFM0[2])):
        x_range = abs(res_cFM0[0]/2*np.cos(res_cFM0[2]))
    else:
        x_range = abs(res_cFM0[1]/2*np.sin(res_cFM0[2]))
    if abs(res_cFM0[1]/2*np.cos(res_cFM0[2]))>abs(res_cFM0[0]/2*np.sin(res_cFM0[2])):
        y_range = abs(res_cFM0[1]/2*np.cos(res_cFM0[2]))
    else:
        y_range = abs(res_cFM0[0]/2*np.sin(res_cFM0[2]))


    # if results2[0]>results2[1]:
    #     x_range = results2[0]/2*np.cos(results2[2])
    #     y_range = results2[0]/2*np.sin(results2[2])
    # else:
    #     x_range = -results2[1]/2*np.sin(results2[2])
    #     y_range = results2[1]/2*np.cos(results2[2])
    # ax1.set_xlim(-1,1)
    # ax1.set_ylim(-1,1)
    ax1.set_xlim( ref_values[var1]-x_range*1.5,ref_values[var1]+x_range*1.5 )
    ax1.set_ylim( ref_values[var2]-y_range*1.5, ref_values[var2]+y_range*1.5 )
    ax1.set_xlabel("$%s$" %(var1))
    ax1.set_ylabel("$%s$" %(var2))
    ax1.grid(True)
    ax1.legend(loc='lower right')
    # ax1.legend().set_visible(False)
    fig.savefig(OUTPUT_FOLDER+'compare_vars_%d_%d.pdf' %(n_var[var1],n_var[var2]))



# var1, var2 = "h", "\\Omega_m"
# var1, var2 = "h", "\\Omega_b"
# var1, var2 = "h", "w_p"
var1, var2 = "w_p", "b_3"
n_var1, n_var2 = n_var[var1], n_var[var2]
compare(uFM0,cFM0,uFM1,cFM1,var1,var2)

# for i in range(15):
#     for j in range(i+1,15):
#         print_ellipse(uFM,cFM,var_names[i],var_names[j])
#
