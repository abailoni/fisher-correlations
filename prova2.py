import sys
import numpy as np
import sympy as sym
import scipy as scy
import matplotlib.pyplot as pl
import os
import time

from scipy.integrate import quad
from scipy.optimize import newton
from os import walk
from re import match
from sympy import srepr
from re import findall

# - Return variable of an expression in order
def varExpr(expression):
    variables_Str = list(set(findall(r"Symbol\(\'(\w+)\'\)", srepr(expression))))
    return [ sym.symbols(var_name) for var_name in variables_Str]

# Return just strings: (instead of sympy variables)
def varExpr2(expression):
    variables_Str = list(set(findall(r"Symbol\(\'(\w+)\'\)", srepr(expression))))
    return [var_name for var_name in variables_Str]

def SymToPy(expression):
    return [sym.lambdify(varExpr(expression), expression), varExpr2(expression)]

def fnExpr(expression):
    return sym.lambdify(varExpr(expression), expression)

def fnEv(data,**args):
    arguments = []
    for option in data[1]:
        if option in args:
            arguments.append(args[option])
        else:
            print "Error! Wrong or missing arguments! Accepted:"
            print data[1]
            break
    if len(data[1])==len(arguments):
        return data[0](*arguments)

# 'var' must be a sympy variable, not a string
# The slow parts are given by the two lambda functions and the lists joining
def NInt(expression, var, start, stop, **args):
    data = SymToPy(expression)
    arguments = []
    for option in data[1]:
        if option in args:
            arguments.append(args[option])
        elif option!=sym.symbol(var):
            print "Error! Wrong or missing arguments! Accepted:"
            print data[1]
            return
    var_array = varExpr(expression)
    if var in var_array:
        index_var = var_array.index(var)
        arguments_before = arguments[0:index_var-1]
        arguments_after = arguments[index_var:]
    else:
        print "Not valid integration variable inserted"
        return
    return quad(lambda x: data[0](*[arguments_before]+[x]+[arguments_after]),start,stop,epsrel=INT_PREC)[0]


Omega_m, h, Omega_b, bias, w_p, w_1, n_s, gamma, sigma8 = sym.symbols('Omega_m h Omega_b bias w_p w_1 n_s gamma sigma8')

# Other variables:
z = sym.symbols('z')
# Toy variables:
zx = sym.symbols('zx') #Actually this is always a number...

# Pivot redshift:
z_p = 0
# Cosmological constant w(z):
w = w_p + w_1*(z-z_p)

# Dark Energy density (flat universe k=0):
Omega_DE = 1 - Omega_m

Hub = sym.sqrt( Omega_m*(1+z)**3 + Omega_DE*sym.exp(3* sym.integrate( (1+w.subs(z,zx))/(1+zx), (zx,0,z)) ) ) #[z, w_1, w_p, Omega_m]

Omega_m_z = Omega_m * (1+z)**3 / Hub**2
Om_DE = SymToPy(Hub)
Om_DE_2 = fnExpr(Hub) #['z', 'w_1', 'w_p', 'Omega_m']

# NIntegrate( Om_m_z**gamma/(1+z), z, 0., zx , [ref_values['gamma'], ref_values['w_1'], ref_values['w_p'], ref_values['Om_m']] )


start=time.clock()
for i in range(15000):
    fnEv(Om_DE,Omega_m=0.7,z=1.,w_1=0.3,w_p=0.4)
stop=time.clock()
print stop-start
start=time.clock()
for i in range(15000):
    Om_DE_2(1.,0.3,0.4,0.7)
stop=time.clock()
print stop-start

ciaone =[2,3,4,51,3]
print ciaone[:2]+[3.]+ciaone[2:]









import StringIO
from re import sub

def SympyCythonization(SympyExpression, fnName):
    #---------------------------------
    # CAPTURE EXPRESSION in a string:
    #---------------------------------
    codeOut = StringIO.StringIO()
    # Capturing output expression:
    sys.stdout = codeOut
    print SympyExpression
    sys.stdout = sys.__stdout__
    # Check:
    stringExpression = codeOut.getvalue()
    print stringExpression
    codeOut.close()
    #---------------------------------
    # Modifing the string:
    #--------------------------------- (sdsds)**2
    test_string = "(-Piecewise((0, sqrt(k**2 - 2*k*kp*z + kp**2) == 0), (-r3*cos(r3*sqrt(k**2 - 2*k*kp*z + kp**2))/sqrt(k**2 - 2*k*kp*z + kp**2) + sin(r3*sqrt(k**2 - 2*k*kp*z + kp**2))/(k**2 - 2*k*kp*z + kp**2), True))/sqrt(k**2 - 2*k*kp*z + kp**2) + Piecewise((0, sqrt(k**2 - 2*k*kp*z + kp**2) == 0), (-r4*cos(r4*sqrt(k**2 - 2*k*kp*z + kp**2))/sqrt(k**2 - 2*k*kp*z + kp**2) + sin(r4*sqrt(k**2 - 2*k*kp*z + kp**2))/(k**2 - 2*k*kp*z + kp**2), True))/sqrt(k**2 - 2*k*kp*z + kp**2))*(-(-k + kp*z)*Piecewise((0, sqrt(k**2 - 2*k*kp*z + kp**2) == 0), (-r1*cos(r1*sqrt(k**2 - 2*k*kp*z + kp**2))/sqrt(k**2 - 2*k*kp*z + kp**2) + sin(r1*sqrt(k**2 - 2*k*kp*z + kp**2))/(k**2 - 2*k*kp*z + kp**2), True))/(k**2 - 2*k*kp*z + kp**2)**(3/2) + (-k + kp*z)*Piecewise((0, sqrt(k**2 - 2*k*kp*z + kp**2) == 0), (-r2*cos(r2*sqrt(k**2 - 2*k*kp*z + kp**2))/sqrt(k**2 - 2*k*kp*z + kp**2) + sin(r2*sqrt(k**2 - 2*k*kp*z + kp**2))/(k**2 - 2*k*kp*z + kp**2), True))/(k**2 - 2*k*kp*z + kp**2)**(3/2) - Piecewise((0, sqrt(k**2 - 2*k*kp*z + kp**2) == 0), (r1**2*(k - kp*z)*sin(r1*sqrt(k**2 - 2*k*kp*z + kp**2))/(k**2 - 2*k*kp*z + kp**2) - r1*(-k + kp*z)*cos(r1*sqrt(k**2 - 2*k*kp*z + kp**2))/(k**2 - 2*k*kp*z + kp**2)**(3/2) + r1*(k - kp*z)*cos(r1*sqrt(k**2 - 2*k*kp*z + kp**2))/(k**2 - 2*k*kp*z + kp**2)**(3/2) + (-2*k + 2*kp*z)*sin(r1*sqrt(k**2 - 2*k*kp*z + kp**2))/(k**2 - 2*k*kp*z + kp**2)**2, True))/sqrt(k**2 - 2*k*kp*z + kp**2) + Piecewise((0, sqrt(k**2 - 2*k*kp*z + kp**2) == 0), (r2**2*(k - kp*z)*sin(r2*sqrt(k**2 - 2*k*kp*z + kp**2))/(k**2 - 2*k*kp*z + kp**2) - r2*(-k + kp*z)*cos(r2*sqrt(k**2 - 2*k*kp*z + kp**2))/(k**2 - 2*k*kp*z + kp**2)**(3/2) + r2*(k - kp*z)*cos(r2*sqrt(k**2 - 2*k*kp*z + kp**2))/(k**2 - 2*k*kp*z + kp**2)**(3/2) + (-2*k + 2*kp*z)*sin(r2*sqrt(k**2 - 2*k*kp*z + kp**2))/(k**2 - 2*k*kp*z + kp**2)**2, True))/sqrt(k**2 - 2*k*kp*z + kp**2))"
    test_string = sub(r'Piecewise.+?(?=\), \()\), \((.+?(?=,)), True\)\)', r'(\1)', test_string.rstrip())
    test_string = sub(r'\(([^()]|(?R))*\)\*\*', r'(\1)', test_string.rstrip())

    print test_string
    return 0.

SympyCythonization(Omega_DE, "fdfd")

'''
code = """
def f(%s):
    x = x + 1
    return x

print 'This is my output.'
""" %("x")

# capture output and errors
# sys.stdout = codeOut
# sys.stderr = codeErr

exec code
print f(3)

# restore stdout and stderr
# sys.stdout = sys.__stdout__
# sys.stderr = sys.__stderr__

# s = codeErr.getvalue()

# print "error:\n%s\n" % s

# s = codeOut.getvalue()

# print "output:\n%s" % s

# codeOut.close()
# codeErr.close()


'''
