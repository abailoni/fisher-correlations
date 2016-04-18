# - CHANGE ALL VARIABLE NAMES!
# - set all numbers of variable from 1 to 7 (not zero...) or put spectrum in another place

import numpy as np


ref_values = {'Om_m': 0.25, 'h': 0.7, 'Om_b': 0.0445, 'n_s': 0.967, 'Om_k':0., 'gamma': 0.545, 'w_p':-0.9, 'w_1':0.}


# Options:
#   - the names of the variables are: 'h', 'Om_m', 'Om_b', 'n_s', 'gamma', 'Om_k', 'w_0', 'w_1'
def set_ref_values(**args):
    var_names = ['h', 'Om_m', 'Om_b', 'n_s', 'gamma', 'Om_k', 'w_0', 'w_1']
    for var_name in args:
        if var_name in var_names:
            ref_values[var_name] = args[var_name]
        else:
            print "Not valid option inserted!!"
            break



z_in = np.array( [0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1] )
N_bins = z_in.shape[0]-1
z_avg, dz = np.zeros(N_bins), np.zeros(N_bins)
survey_area = 15000.


#
# After updating these data it is necessary to compute again distancese and derivatives!!
#
#
# Options:
#   - "bins_list": the python list should include both extremes. Thus the last number is the border of the last bin
#
#   - "survey_area": area in sq deg [default 15000]

def set_survey(**options):
    options_val = ["bins_list", "survey_area"]
    for option in options:
        if option not in options_val:
            print "Not valid option inserted!!"
            break
        if option=="bins_list":
            global N_bins, z_in
            N_bins = len(option["bins_list"])-1
            z_in = np.array(option["bins_list"])
        if option=="survey_area":
            global survey_area
            survey_area = option["survey_area"]


def compute_survey_DATA():
    # Setting bins:
    global z_avg, dz
    z_avg = np.array([ (z_in[i]+z_in[i+1])/2. for i in range(N_bins)])
    dz = np.array([z_in[i+1]-z_in[i] for i in range(N_bins)])
    # Distances:
    global com_zbin, com_zbin_avg
    com_zbin = np.array([ comov_dist(zbn) for zbn in z_in])
    com_zbin_avg = np.array([ (com_zbin[i]+com_zbin[i+1])/2. i in range(N_bins)])
    # Derivatives and funtions:
    global bias_bins, lnG_der_data, Beta_der_data, lnH_der_data, lnD_der_data, Growth_bins, beta_bins
    bias_bins = np.array([bias(zx) for zx in z_avg])
    Growth_bins = np.array([Growth(zx) for zx in z_avg])
    beta_bins = np.array([ beta(bin) for bin in range(N_bins)])
    four_parameters = ['Omega_m', 'w_p', 'w_1','gamma']
    for var in four_parameters: # num_var = [3-5] + gamma
        for bin in range(N_bins):
            # for the first two also add bias:
            lnG_der_data[n_var[var]][bin] = lnG_der[var](z_avg[bin])
            Beta_der_data[n_var[var]][bin] = 1./bias_bins[bin] * Beta_der[var](z_avg[bin])
            lnH_der_data[n_var[var]][bin] = lnH_der[var](z_avg[bin])
            lnD_der_data[n_var[var]][bin] = lnD_der[var](z_avg[zbn])




















