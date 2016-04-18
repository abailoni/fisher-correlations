import time

# start = time.clock()
import modules.uFM as uFM
uFM.compute_CAMB_spectra()
uFM.compute_survey_DATA()
print uFM.vol_shell_py(0)
# uFM.FM(1,0.2)

# stop = time.clock()
# print "Time: %g seconds" %(stop-start)
