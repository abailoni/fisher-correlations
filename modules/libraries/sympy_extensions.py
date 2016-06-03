###################################
# Some routine useful for Sympy:
###################################

import sys
import numpy as np
import sympy as sym

# from theano.compile.io import In
# from sympy.printing.theanocode import theano_function
# from sympy.printing.theanocode import *

from scipy.integrate import quad

from sympy.utilities.autowrap import ufuncify

from sympy import srepr
from re import findall

INT_PREC = 1e-5


# Check if we have a numpy array:
def numpy_check(data):
    return True if isinstance(data, np.ndarray) else False


# - Return variables (symbols) of an expression in order
def varExpr(expression):
    variables_Str = list(set(findall(r"Symbol\(\'(\w+)\'\)", srepr(expression))))
    return [ sym.symbols(var_name) for var_name in variables_Str]
# - Return variables (names) of an expression in order
def varExprNames(expression):
    variables_Str = list(set(findall(r"Symbol\(\'(\w+)\'\)", srepr(expression))))
    return [var_name for var_name in variables_Str]


# --------------------------------------
# THE UFUNCIFY EVALUATION FUNCTION:
# --------------------------------------
#
#
# The following routine converts a sympy expression to a compiled and
# optimized Theano function. It accepts default values or input-ordering
# as optional inputs.
#
# INPUTS:
#   - expression: sympy expression
#   - optional arguments:
#        ----> default values for the values at which the expression should be
#              evaluated. An example of arguments could be Omega_m=0.3, h=0.7
#              and so on. The names of the arguments should match the names of
#              the sympy variables.
#              Passed arguments that are not necessary for the evaluation of
#              the expression are ignored.
#        ----> "order_vars": a list containing all the variable names required
#              for the evaluation of the expression which were not passed as
#              default values. The Theano function in output will follow this
#              ordering of the variables.
#              If some variables are missing an error is returned. If too many
#              variables are passed (or some of them are included in the
#              default ones) they are ignored. The list returned as second
#              output will show just the required ones.
#        ----> "NIntegrate": if set to something then it return a function with the default values fixed (can not be changed)
#
# OUTPUTS:
#   - Theano function: it will be of the form fun(var1, var2, def_var1=val1,
#     def_var2=val2) so the default values passed can be modified by calling
#     the argument names. The inputs without a default values are mandatory
#     (and if an ordering was passed, they will follow it)
#   - list with the names of the required inputs in the default or
#     passed order
#   - list with the names of the arguments to call to modify the default
#     values. Without default values passed it is just an empty list. Of
#     course there is no order to follow for these optional arguments.
#
# MISSING:
#   - Missing optimisation of the input types (float64, etc...)


def SymToUfuncify(expression,**args):
    # expression = sym.simplify(expression)
    required_vars = varExprNames(expression)
    if len(args)!=0:
        def_varNames, def_values = [], []

        # Structure the default values:
        for var in args:
            if var in required_vars:
                def_varNames.append(var)
                def_values.append(args[var])
        def_symVar = [sym.symbols(var_name) for var_name in def_varNames]

        # Find the missing ones and order them if an order is passed:
        missing_varNames = []
        for var in required_vars:
            if var not in def_varNames:
                missing_varNames.append(var)
        if "order_vars" in args:
            # Check not to have too many variables and if so delete them:
            for i, var in enumerate(args["order_vars"]):
                if var not in missing_varNames:
                    args["order_vars"].pop(i)
            # Check to have all of them:
            if len(args["order_vars"])!=len(missing_varNames):
                print "Not all vars passed! Required:"
                print missing_varNames
                return -1
            missing_varNames = args["order_vars"]
        missing_symVar = [sym.symbols(var_name) for var_name in missing_varNames]

        # Compile the ufuncify function:
        if "NIntegrate" in args:
            expression.subs( [ (symVar,val) for symVar, val in zip(def_symVar,def_values) ] )
            # Output a one-only input function:
            return ufuncify(missing_symVar, expression)
        else:
            return [ufuncify(missing_symVar+def_symVar, expression), missing_varNames, def_varNames, def_values]
    else:
        # Just use the default sympy ordering:
        return [ufuncify(varExpr(expression), expression), required_vars]


# Evaluation:
# The arguments are used only if it is necessary to change some default data (slow). Not necessary arguments are ignored.
# inputs are the ordered mandatory inputs (they can be numpy arrays).
def unfun_Ev(SymToUfuncify_data,*inputs,**kargs):
    # Change eventual default data:
    defValues = list(SymToUfuncify_data[3])
    if len(kargs)!=0:
        for i, def_varName in enumerate(SymToUfuncify_data[2]):
            if def_varName in kargs:
                defValues[i] = kargs[def_varName]
    # With numpy arrays Vectorize default inputs:
    if numpy_check(inputs[0]):
        if len(defValues)!=0:
            defValues = [np.ones(inputs[0].shape)*value for value in defValues]
    return SymToUfuncify_data[0](*inputs+tuple(defValues))


#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------

# Accept option "numpy"=True for a faster vectorised version. Default is False
#
def SymToLambda(expression,**args):
    numpy_option = []
    if "numpy" in args:
        numpy_option = ['numpy']
    # expression = sym.simplify(expression)
    required_vars = varExprNames(expression)
    if len(args)!=0:
        def_varNames, def_values = [], []

        # Structure the default values:
        for var in args:
            if var in required_vars:
                def_varNames.append(var)
                def_values.append(args[var])
        def_symVar = [sym.symbols(var_name) for var_name in def_varNames]

        # Find the missing ones and order them if an order is passed:
        missing_varNames = []
        for var in required_vars:
            if var not in def_varNames:
                missing_varNames.append(var)
        if "order_vars" in args:
            # Check not to have too many variables and if so delete them:
            for i, var in enumerate(args["order_vars"]):
                if var not in missing_varNames:
                    args["order_vars"].pop(i)
            # Check to have all of them:
            if len(args["order_vars"])!=len(missing_varNames):
                print "Not all vars passed! Required:"
                print missing_varNames
                return -1
            missing_varNames = args["order_vars"]
        missing_symVar = [sym.symbols(var_name) for var_name in missing_varNames]

        # Compile the ufuncify function:
        if "NIntegrate" in args:
            if len(missing_varNames)>1:
                print "Attention: not all the necessary arguments passed!"
                print "Missing: (ignore integration variable)"
                print missing_varNames
                return -1
            expression = expression.subs( [ (symVar,val) for symVar, val in zip(def_symVar,def_values) ] )
            # Output a one-only input function:
            return sym.lambdify(missing_symVar, expression,*numpy_option)
        else:
            return [sym.lambdify(missing_symVar+def_symVar, expression,*numpy_option), missing_varNames, def_varNames, def_values]
    else:
        # Just use the default sympy ordering:
        return [sym.lambdify(varExpr(expression), expression,*numpy_option), required_vars]



# Evaluation:
# The arguments are used only if it is necessary to change some default data (slow). Not necessary arguments are ignored.
# inputs are the ordered mandatory inputs (they can be numpy arrays).
# Now also matrices are supported, so the kargs can be arrays of parameters
def Lambda_Ev(SymToLambda_data,*inputs,**kargs):
    # Change eventual default data:
    defValues = list(SymToLambda_data[3])
    if len(kargs)!=0:
        for i, def_varName in enumerate(SymToLambda_data[2]):
            if def_varName in kargs:
                defValues[i] = kargs[def_varName]
    # With numpy arrays Vectorize default inputs:
    if numpy_check(inputs[0]):
        if len(defValues)!=0:
            # Check if kargs are also arrays:
            if numpy_check(defValues[0]):
                defValues = [np.ones(inputs[0].shape)*np.reshape(value,(-1,1)) for value in defValues]
            else:
                defValues = [np.ones(inputs[0].shape)*value for value in defValues]
    return SymToLambda_data[0](*inputs+tuple(defValues))


# --------------------------------------------
# The definitive NIntegrate tools:
# --------------------------------------------
# The following routines transform slow sympy expressions in Theano compiled one and then performe an integration with scipy quad. All the parameters are passed as arguments. No order required.
#
# MISSING POSSIBLE IMPROVEMENTS:
#   - Missing optimisation of the input types (float64, etc...)
#   - GSL integration...?

# ---------------
# MAIN ROUTINE:
# ---------------
#
# Here start and stop can be numpy arrays. In this way the sympy expression is "Theanized" only once.
#
#
def NIntegrate(expression, integrated_varName, start, stop, INT_PREC=1e-6, **args):
    # By setting the "ordering_vars" we check if some arg is missing [if so SymToLambda gives an error]
    lamda_fun = SymToLambda(expression,order_vars=[integrated_varName],NIntegrate=True,**args)
    # Check if start and stop are numpy arrays:
    if numpy_check(start) and numpy_check(stop):
        if start.shape==stop.shape:
            results = np.empty(start.shape)
            for i in range(start.shape[0]):
                results[i] = quad(lamda_fun,start[i],stop[i],epsrel=INT_PREC)[0]
            return results
    return quad(lamda_fun,start,stop,epsrel=INT_PREC)[0]

# --------------------
# AUXILIARY ROUTINES:
# --------------------
#
# The following two auxiliary routines do respectively:
#    - return the compiled Theano fct required for a quad integration
#    - perform the quad integration with an already compiled Theano fct
#
def integrationFun(expression, integrated_varName, **args):
    return SymToLambda(expression,order_vars=[integrated_varName],NIntegrate=True,**args)

def NIntegrate_fun(integration_fun,start,stop,INT_PREC=1e-6): #start and stop can be arrays
    if numpy_check(start) and numpy_check(stop):
        if start.shape==stop.shape:
            results = np.empty(start.shape)
            for i in range(start.shape[0]):
                results[i] = quad(integration_fun,start[i],stop[i],epsrel=INT_PREC)[0]
            return results
    return quad(integration_fun,start,stop,epsrel=INT_PREC)[0]

#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------


# OLD THEANO CODE NOT WORKING

# # --------------------------------------
# # THE DEFINITIVE EVALUATION FUNCTION:
# # --------------------------------------
# # To good to be true... Probably not that fast and not working...
# #
# #
# # The following routine converts a sympy expression to a compiled and
# # optimized Theano function. It accepts default values or input-ordering
# # as optional inputs.
# #
# # INPUTS:
# #   - expression: sympy expression
# #   - optional arguments:
# #        ----> default values for the values at which the expression should be
# #              evaluated. An example of arguments could be Omega_m=0.3, h=0.7
# #              and so on. The names of the arguments should match the names of
# #              the sympy variables.
# #              Passed arguments that are not necessary for the evaluation of
# #              the expression are ignored.
# #        ----> "order_vars": a list containing all the variable names required
# #              for the evaluation of the expression which were not passed as
# #              default values. The Theano function in output will follow this
# #              ordering of the variables.
# #              If some variables are missing an error is returned. If too many
# #              variables are passed (or some of them are included in the
# #              default ones) they are ignored. The list returned as second
# #              output will show just the required ones.
# #
# # OUTPUTS:
# #   - Theano function: it will be of the form fun(var1, var2, def_var1=val1,
# #     def_var2=val2) so the default values passed can be modified by calling
# #     the argument names. The inputs without a default values are mandatory
# #     (and if an ordering was passed, they will follow it)
# #   - list with the names of the required inputs in the default or
# #     passed order
# #   - list with the names of the arguments to call to modify the default
# #     values. Without default values passed it is just an empty list. Of
# #     course there is no order to follow for these optional arguments.
# #
# # MISSING:
# #   - Missing optimisation of the input types (float64, etc...)


# def SymToTheano(expression,**args):
#     expression = sym.simplify(expression)
#     required_vars = varExprNames(expression)
#     if len(args)!=0:
#         def_varNames, def_values = [], []

#         # Structure the default values:
#         for var in args:
#             if var in required_vars:
#                 def_varNames.append(var)
#                 def_values.append(args[var])
#         def_symVar = [sym.symbols(var_name) for var_name in def_varNames]
#         theano_default_inputs = [In(theano_code(sym_var), value=def_val, name=varName) for sym_var, def_val, varName in zip(def_symVar,def_values,def_varNames)]

#         # Find the missing ones and order them if an order is passed:
#         missing_varNames = []
#         for var in required_vars:
#             if var not in def_varNames:
#                 missing_varNames.append(var)
#         if "order_vars" in args:
#             # Check not to have too many variables and if so delete them:
#             for i, var in enumerate(args["order_vars"]):
#                 if var not in missing_varNames:
#                     args["order_vars"].pop(i)
#             # Check to have all of them:
#             if len(args["order_vars"])!=len(missing_varNames):
#                 print "Not all vars passed! Required:"
#                 print missing_varNames
#                 return -1
#             missing_varNames = args["order_vars"]
#         missing_symVar = [sym.symbols(var_name) for var_name in missing_varNames]
#         print def_varNames, missing_varNames
#         # Compile the Theano function:
#         return theano_function(missing_symVar+theano_default_inputs, [expression]), missing_varNames, def_varNames
#     else:
#         # Just use the default sympy ordering:
#         return theano_function(varExpr(expression), [expression]), required_vars, []


# # Slow (obsolete) evaluation:
# def theano_simpleEv(symToTheano_data,**vars):
#     arguments = []
#     for option in symToTheano_data[1]:
#         if option in vars:
#             arguments.append(vars[option])
#         else:
#             print "Error! Wrong or missing arguments! Required:"
#             print symToTheano_data[1]
#             break
#     return symToTheano_data[0](*arguments)


# # --------------------------------------------
# # The definitive NIntegrate tools:
# # --------------------------------------------
# # The following routines transform slow sympy expressions in Theano compiled one and then performe an integration with scipy quad. All the parameters are passed as arguments. No order required.
# #
# # MISSING POSSIBLE IMPROVEMENTS:
# #   - Missing optimisation of the input types (float64, etc...)
# #   - GSL integration...?

# # ---------------
# # MAIN ROUTINE:
# # ---------------
# #
# # Here start and stop can be numpy arrays. In this way the sympy expression is "Theanized" only once.
# #
# def NIntegrate(expression, integrated_varName, start, stop, **args):
#     # By setting the "ordering_vars" we check if some arg is missing [if so SymToTheano gives an error]
#     theano_function,_,_ = SymToTheano(expression,order_vars=[integrated_varName],**args)
#     # Check if start and stop are numpy arrays:
#     if (type(start).__module__ == np.__name__) and (type(stop).__module__ == np.__name__):
#         if start.shape==stop.shape:
#             # results = np.empty(start.shape)
#             # for i in range(start.shape):
#             #     results[i] = quad(theano_function,start[i],stop[i],epsrel=INT_PREC)[0]
#             results = quad(theano_function,start,stop,epsrel=INT_PREC)[0]
#             return results
#     return quad(theano_function,start,stop,epsrel=INT_PREC)[0]

# # --------------------
# # AUXILIARY ROUTINES:
# # --------------------
# #
# # The following two auxiliary routines do respectively:
# #    - return the compiled Theano fct required for a quad integration
# #    - perform the quad integration with an already compiled Theano fct
# #
# def Theano_integrationFun(expression, integrated_varName, **args):
#     # By setting the "ordering_vars" we check if some arg is missing [if so SymToTheano gives an error]
#     theano_function,_,_ = SymToTheano(expression,order_vars=[integrated_varName],**args)
#     return theano_function

# def NInt_TheanoFun(theano_fct,start,stop): #start and stop can be arrays
#     if (type(start).__module__ == np.__name__) and (type(stop).__module__ == np.__name__):
#         if start.shape==stop.shape:
#             results = np.empty(start.shape)
#             for i in range(start.shape):
#                 results[i] = quad(theano_function,start[i],stop[i],epsrel=INT_PREC)[0]
#             return results
#     return quad(theano_function,start,stop,epsrel=INT_PREC)[0]

# #-------------------------------------------------------------------------
# #-------------------------------------------------------------------------
# #-------------------------------------------------------------------------


# --------------------------
# Old shitty things:
# --------------------------

# CONVERSION TO FASTER COMPILED VERSIONS:
# Convert to a function lambda:
def fnExpr(expression):
    return sym.lambdify(varExpr(expression), expression)


# Convert function to lamda without requiring order of arguments,
# but 3 times slower wrt fnExpr:
def SymToPy(expression):
    return [sym.lambdify(varExpr(expression), expression), varExprNames(expression)]
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
# The slow parts are given by the evaluation of the two lambda functions and the lists joining (but that perhaps happens only once..)
def NInt(expression, var, start, stop, **args):
    SymToPydata = SymToPy(expression)
    arguments = []
    for option in SymToPydata[1]:
        if option in args:
            arguments.append(args[option])
        elif sym.symbols(option)!=var:
            print "Error! Wrong or missing arguments! Accepted:"
            print SymToPydata[1]
            return -1
    var_array = varExpr(expression)
    if var in var_array:
        index_var = var_array.index(var)
        arguments_before = arguments[:index_var]
        arguments_after = arguments[index_var:]
    else:
        print "Not valid integration variable inserted"
        return -1
    return quad(lambda x: SymToPydata[0](*arguments_before+[x]+arguments_after),start,stop,epsrel=INT_PREC)[0]



# BOTH HORRIBLE PREFORMANCE because the function is lambdified every time the integrand is called or the sub method is used.
# Furthermore the first version is NOT WORKING! If the parameter is not the first one, it screws up...
# - NIntegrate of a sympy expression :D --> Brao ti...
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
            return quad(fnExpr(expression),start,stop,args=tuple(par_array),epsrel=INT_PREC)[0]


# This version does not care of the order of the parameters
# Format par_array: [('z',3),('Omega_m',3),...]
def NIntegrate3(expression, var, start, stop, par_array=0):
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
