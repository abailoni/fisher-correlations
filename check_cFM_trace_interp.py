import modules.cFM as cFM
import matplotlib.pyplot as pl
from scipy.interpolate import spline, interp1d, Akima1DInterpolator
from scipy.signal import resample
import numpy as np

import time

cFM.init()
var1, var2 = 0,0
# start=time.time()
# cFM.init_Trace_term(var1,var2,50,10)
# print "First method: %f sec." %(time.time()-start)
# print cFM.fisher_matrix_element(var1,var2,0,True)
# start=time.time()
# cFM.adapt_trace_tot(var1,var2)
# print "Second method: %f sec." %(time.time()-start)
# print cFM.fisher_matrix_element(var1,var2,0,True)

# quit()
x, y = cFM.adapt_trace_term(var1, var2, 50, -1., 0.5, 10, 20)
x_smooth = np.linspace(1e-3,0.2,3000)
y_smooth = Akima1DInterpolator(x,y)(x_smooth)
print x.shape[0]

# x1, y1 = cFM.adapt_trace_term(var1, var2, 20, -1., 0.1)
# print x1.shape[0]

x2, y2 = cFM.adapt_trace_term(var1, var2, 10, -1., 0.00000000005, 10, 20)
x2_smooth = np.linspace(1e-3,0.2,3000)
# y2_smooth = interp1d(x2,y2,kind="slinear")(x2_smooth)
y2_smooth2 = Akima1DInterpolator(x2,y2)(x2_smooth)
# y2_smooth2, x2_smooth2 = resample(y2,3000,x2)
y2_smooth3 = interp1d(x2,y2,kind="cubic")(x2_smooth)
print x2.shape[0]


fig2=pl.figure()
ax1=fig2.add_subplot(111)
ax1.plot(x2_smooth,cFM.prova(x2_smooth,y2_smooth2),'b-',label="Trace along $\\mu=-1$")
ax1.plot(x2,cFM.prova(x2,y2),'b.',label="Trace along $\\mu=-1$")
# ax1.plot(x,y,'b.')

# ax1.plot(x2_smooth,y2_smooth2,'y-')
# ax1.plot(x2,y2,'g.')

# ax1.plot(k_vect_2,trace_int,'b-',label="Trace interp $\\mu=-1$")
# ax1.plot(vect_k,class_fct['P_0'](vect_k),'b-',label="Class Spectrum")
ax1.grid(True)
ax1.legend(loc='best')
# ax1.set_yscale('symlog')
# ax1.set_xlabel("$k$ [$h$/Mpc]")
# ax1.set_ylabel("$P(x)$ [(Mpc/$h$)$^3$]")
fig2.savefig('plots/trace/vabbe_advanced_%d_%d.pdf' %(var1,var2))


'''

cFM.init_Trace(50,10)
cFM.FM(0,"EUCLID2012-interpolated")


computed =  234276.21122811467

# cFM.init_Trace_term(var1,var2,80,4)
# this = cFM.fisher_matrix_element(var1,var2,0,True)
# print "80,4: %g --> Difference: %g" %(this, computed-this)

# cFM.init_Trace_term(var1,var2,50,8)
# this = cFM.fisher_matrix_element(var1,var2,0,True)
# print "50,8: %g --> Difference: %g" %(this, computed-this)

# cFM.init_Trace_term(var1,var2,80,8)
# this = cFM.fisher_matrix_element(var1,var2,0,True)
# print "80,8: %g --> Difference: %g" %(this, computed-this)

cFM.init_Trace_term(var1,var2,25,10)
this = cFM.fisher_matrix_element(var1,var2,0,True)
print "25,10: %g --> Difference: %g" %(this, computed-this)

cFM.init_Trace_term(var1,var2,35,10)
this = cFM.fisher_matrix_element(var1,var2,0,True)
print "35,10: %g --> Difference: %g" %(this, computed-this)

cFM.init_Trace_term(var1,var2,50,10)
this = cFM.fisher_matrix_element(var1,var2,0,True)
print "50,10: %g --> Difference: %g" %(this, computed-this)

cFM.init_Trace_term(var1,var2,80,10)
this = cFM.fisher_matrix_element(var1,var2,0,True)
print "80,10: %g --> Difference: %g" %(this, computed-this)

cFM.init_Trace_term(var1,var2,100,10)
this = cFM.fisher_matrix_element(var1,var2,0,True)
print "100,10: %g --> Difference: %g" %(this, computed-this)
'''
