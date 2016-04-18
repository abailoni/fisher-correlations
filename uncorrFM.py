import sys
import numpy as np
check_AP, name, N_bins, check_W = int(sys.argv[1]), sys.argv[2], int(sys.argv[3]), int(sys.argv[4])

import modules.uFM as uFM
if N_bins==14:
    uFM.init()
    uFM.FM(check_AP, check_W, name)
elif N_bins==8:
    # Amendola's settings:
    dens = np.array([3.56, 2.42, 1.81, 1.44, 0.99, 0.55, 0.29, 0.15])*0.5
    uFM.set_survey(bins_list=[0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1], dens_list=[ds for ds in dens],bias_list=[1.]*8)
    uFM.init()
    uFM.FM(check_AP, check_W, name)
elif N_bins==7:
    # Half EUCLID-2012 bins:
    z_in = [0.65, 0.85, 1.05, 1.25, 1.45, 1.65, 1.85, 2.05]
    z_avg = [ (z_in[i]+z_in[i+1])/2. for i in range(N_bins)]
    dens = uFM.density_py(np.array(z_avg))
    bias = uFM.bias_py(z_avg)
    uFM.set_survey(bins_list=z_in, dens_list=dens,bias_list=bias)
    uFM.init()
    uFM.FM(check_AP, check_W, name)
