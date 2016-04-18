import modules.uFM as uFM

binList = [0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05]


uFM.set_survey(bins_list=binList, dens_list=[1.]*(len(binList)-1) )
uFM.import_CLASS_data()
uFM.compute_survey_DATA()
uFM_14_0 = uFM.FM()
uFM_14_1 = uFM.FM(1)



binList = [0.65, 0.85, 1.05, 1.25, 1.45, 1.65, 1.85, 2.05]
uFM.set_survey(bins_list=binList, dens_list=[1.]*(len(binList)-1) )
uFM.import_CLASS_data()
uFM.compute_survey_DATA()
uFM_7_0 = uFM.FM()
uFM_7_1 = uFM.FM(1)


import sys
import numpy as np
import sympy as sym
import scipy as scy
import matplotlib.pyplot as pl
import os
import time

OUTPUT_FOLDER = "plots/marginalisation/comparison/test_density/"


# Correction factor volume:
correction_factor = np.pi * 15000. / (180.**2 * 4.)
uFM_7_0, uFM_14_0 = uFM_7_0*correction_factor, uFM_14_0*correction_factor
uFM_7_1, uFM_14_1 = uFM_7_1*correction_factor, uFM_14_1*correction_factor


# ------------------------------------
# Computing uncertaintes: (finally)
# ------------------------------------

n_var = {"h": 0, "\\Omega_b": 1, "n_s": 2, "\\Omega_m": 3, "w_p": 4, "w_1": 5, "\\gamma": 6, "b_0": 7,"b_1": 8,"b_2": 9,"b_3": 10,"b_4": 11,"b_5": 12,"b_6": 13, "b_7": 14}
var_names = ["h", "\\Omega_b", "n_s", "\\Omega_m", "w_p", "w_1", "\\gamma", "b_0","b_1","b_2","b_3","b_4","b_5","b_6", "b_7"]
ref_values = {'\\Omega_m': 0.25, 'h': 0.7, '\\Omega_b': 0.0445, 'n_s': 0.967, "w_p": -0.9, "w_1": 0., "\\gamma": 0.545, "b_0": 1.,"b_1": 1.,"b_2": 1.,"b_3": 1.,"b_4": 1.,"b_5": 1.,"b_6": 1., "b_7": 1.}

uFM_7_0_inv, uFM_14_0_inv = np.linalg.inv(uFM_7_0), np.linalg.inv(uFM_14_0)
uFM_7_1_inv, uFM_14_1_inv = np.linalg.inv(uFM_7_1), np.linalg.inv(uFM_14_1)

np.savetxt(OUTPUT_FOLDER+"7.txt",uFM_7_0_inv)
np.savetxt(OUTPUT_FOLDER+"14.txt",uFM_14_0_inv)
np.savetxt(OUTPUT_FOLDER+"7_1.txt",uFM_7_1_inv)
np.savetxt(OUTPUT_FOLDER+"14_1.txt",uFM_14_1_inv)

# for i in range(7):
#     sigma7, sigma14 = np.sqrt(uFM_7_0_inv[i,i]), np.sqrt(uFM_14_0_inv[i,i])
#     print "%s\t\tFM_corr=%g \tFM_uncorr=%g  \t%g%%" %(var_names[i],sigma7, sigma14, (sigma7-sigma14)/sigma7*100)

# # Fully marginalised uncertaintes:
# for i in range(cFM.shape[0]):
#     sigmaC, sigmaU = np.sqrt(inv_cFM[i,i]), np.sqrt(inv_uFM[i,i])
#     print "%s\t\tFM_corr=%g \tFM_uncorr=%g  \t%g%%" %(var_names[i],sigmaC, sigmaU, (sigmaC-sigmaU)/sigmaC*100)


# Plot ellipse for two variables:

def marginalise(inv_matrix,var1,var2):
    temp_matrix = np.row_stack((inv_matrix[var1,:],inv_matrix[var2,:]))
    return np.column_stack((temp_matrix[:,var1], temp_matrix[:,var2]))

# Return axes and angle:
def sigma1_ellipse(FM_matrix, var1, var2):
    inv_matrix = np.linalg.inv(FM_matrix)
    sub_matrix = marginalise(inv_matrix,var1,var2)
    # Find axes:
    term1 = (sub_matrix[0,0]+sub_matrix[1,1])/2.
    term2 = np.sqrt((sub_matrix[0,0]-sub_matrix[1,1])**2/4. +sub_matrix[0,1]**2)
    a = np.sqrt( term1+term2 )
    b = np.sqrt( term1-term2 )
    angle = np.arctan( 2*sub_matrix[0,1]/(sub_matrix[0,0]-sub_matrix[1,1]) ) / 2.
    return 2*1.52*a, 2*1.52*b, angle



from math import degrees
import matplotlib.patches as patch

def compare(uFM_7_0,uFM_7_1,uFM_14_0,uFM_14_1,var1,var2): #var names!
    res_uFM0 = sigma1_ellipse(uFM_7_0,n_var[var1],n_var[var2])
    res_cFM0 = sigma1_ellipse(uFM_14_0,n_var[var1],n_var[var2])
    res_uFM1 = sigma1_ellipse(uFM_7_1,n_var[var1],n_var[var2])
    res_cFM1 = sigma1_ellipse(uFM_14_1,n_var[var1],n_var[var2])
    print res_uFM0[0], res_uFM0[1]
    print degrees(res_uFM0[2]), degrees(res_cFM0[2])
    ellipse1 = patch.Ellipse(xy=(ref_values[var1], ref_values[var2]), width=res_uFM0[0], height=res_uFM0[1], angle=degrees(res_uFM0[2]), alpha= 0.7, color="r", label="7 bins - K=0")
    ellipse2 = patch.Ellipse(xy=(ref_values[var1], ref_values[var2]), width=res_cFM0[0], height=res_cFM0[1], angle=degrees(res_cFM0[2]), alpha= 0.7, color="b", label="14 bins - K=0")
    ellipse3 = patch.Ellipse(xy=(ref_values[var1], ref_values[var2]), width=res_uFM1[0], height=res_uFM1[1], angle=degrees(res_uFM1[2]), alpha= 0.7, color="g", label="7 bins - K=1")
    ellipse4 = patch.Ellipse(xy=(ref_values[var1], ref_values[var2]), width=res_cFM1[0], height=res_cFM1[1], angle=degrees(res_cFM1[2]), alpha= 0.7, color="y", label="14 bins - K=0")
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
    ax1.set_xlim( ref_values[var1]-x_range*1.5,ref_values[var1]+x_range*1.5 )
    ax1.set_ylim( ref_values[var2]-y_range*1.5, ref_values[var2]+y_range*1.5 )
    ax1.set_xlabel("$%s$" %(var1))
    ax1.set_ylabel("$%s$" %(var2))
    ax1.grid(True)
    ax1.legend(loc='lower right')
    # ax1.legend().set_visible(False)
    fig.savefig(OUTPUT_FOLDER+'compare_vars_%d_%d.pdf' %(n_var[var1],n_var[var2]))



# var1, var2 = "h", "\\Omega_m"
var1, var2 = "h", "\\Omega_b"
# var1, var2 = "h", "w_p"
# var1, var2 = "w_p", "b_3"
n_var1, n_var2 = n_var[var1], n_var[var2]
# compare(uFM_7_0,uFM_7_1,uFM_14_0,uFM_14_1,var1,var2)

# for i in range(15):
#     for j in range(i+1,15):
#         print_ellipse(uFM,cFM,var_names[i],var_names[j])
#

