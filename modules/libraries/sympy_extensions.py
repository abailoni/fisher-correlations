###################################
# Some routine useful for Sympy:
###################################

import sys
import numpy as np
import sympy as sym
import scipy as scy

from sympy.utilities.autowrap import ufuncify
from scipy.integrate import quad

import sympy as sym

from sympy import srepr
from re import findall

INT_PREC = 1e-3

# - Return variable of an expression in order
def varExpr(expression):
    variables_Str = list(set(findall(r"Symbol\(\'(\w+)\'\)", srepr(expression))))
    return [ sym.symbols(var_name) for var_name in variables_Str]

# Convert to a function lambda:
def fnExpr(expression):
    return sym.lambdify(varExpr(expression), expression)


# Convert function to lamda without requiring order of arguments,
# but 3 times slower wrt fnExpr:
def varExpr2(expression):
    variables_Str = list(set(findall(r"Symbol\(\'(\w+)\'\)", srepr(expression))))
    return [var_name for var_name in variables_Str]
def SymToPy(expression):
    return [sym.lambdify(varExpr(expression), expression), varExpr2(expression)]
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
        elif sym.symbols(option)!=var:
            print "Error! Wrong or missing arguments! Accepted:"
            print data[1]
            return -1
    var_array = varExpr(expression)
    if var in var_array:
        index_var = var_array.index(var)
        arguments_before = arguments[:index_var]
        arguments_after = arguments[index_var:]
    else:
        print "Not valid integration variable inserted"
        return -1
    return quad(lambda x: data[0](*arguments_before+[x]+arguments_after),start,stop,epsrel=INT_PREC)[0]


# BOTH HORRIBLE PREFORMANCE because the function is lambdified every time the integrand is called or the sub method is used.
# Furthermore the first version is NOT WORKING! If the parameter is not the first one, it screw up...
# - NIntegrate of a sympy expression :D --> Brao ti...
def NIntegrate(expression, var, start, stop, par_array=0):
    if par_array==0:
        return quad(fnExpr(expression),start,stop,epsrel=INT_PREC)[0]
    else:
        var_array = varExpr(expression)
        # Put the integrated variable in the front of the array:
        index = var_array.index(var)
        if index!=0:
            var_array.insert(0,var_array.pop(index))
        # Integrate:
        if par_array==-1:
            print "## Inputs: ",var_array[1:]," ##"
            return 0
        else:
            return quad(fnExpr(expression),start,stop,args=tuple(par_array),epsrel=INT_PREC)[0]


# This version does not care of the order of the parameters
# Format par_array: [('z',3),('Omega_m',3),...]
def NIntegrate2(expression, var, start, stop, par_array=0):
    if par_array==0:
        return quad(fnExpr(expression),start,stop,epsrel=INT_PREC)[0]
    else:
        var_array = varExpr(expression)
        # Put the integrated variable in the front of the array:
        index = var_array.index(var)
        if index!=0:
            var_array.insert(0,var_array.pop(index))
        # Integrate:
        if par_array==-1:
            print "## Inputs: ",var_array[1:]," ##"
            return 0
        else:
            return quad(lambda x: expression.subs([(var,x)]+par_array),start,stop,epsrel=INT_PREC)[0]

# - substitute ;)
# Changing info prints the list of required variables
def sub(expression, values, info=0):
    expression2 = expression
    if info!=0:
        print "## Order inputs: ",varExpr(expression)," ##"
        return 0
    else:
        var_array = varExpr(expression)
        return expression.subs( [ (var_array[i],values[i]) for i in range(0,len(var_array)) ] )
