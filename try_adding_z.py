import numpy as np
import modules.cFM as cFM
cFM.compute_CAMB_spectra()

N_bins = 7
z_in = np.linspace(0.65,2.05,N_bins+1)
dens = cFM.density_advanced_py(z_in)
bias = cFM.bias_advanced_py(z_in)
cFM.set_survey(bins_list=z_in, dens_list=dens,bias_list=bias)

cFM.compute_survey_DATA()
cFM.compute_data_adding_z()

# cFM.store_conv_spectra()
# print cFM.obsSpectr_redshift_der(0,0,0.1,1.,0)
# cFM.set_typeFM("correlations+windowFun")
# print cFM.trace_adding_z(0.2, 0.5, 0, 14)
cFM.FM_plusRedshift(0,"test_redshift","uncorr")

# print cFM.trace_adding_z(0.1,1.,1,20)
