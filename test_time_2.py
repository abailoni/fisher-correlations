import time
import numpy as np

start = time.clock()
import modules.uFM as uFM
uFM.set_survey(bins_list=[0.3,0.42,0.54,0.66,0.78,0.9,1.02,1.14,1.26,1.38,1.5,1.62,1.74,1.86,1.98,2.1], dens_list=[4.,3.725,3.45,3.175,2.9,2.625,2.35,2.075,1.8,1.525,1.25,0.975,0.7,0.425,0.15])
uFM.import_CLASS_data()
uFM.compute_survey_DATA()
uFM.FM()
stop = time.clock()
print "Time: %g seconds" %(stop-start)
