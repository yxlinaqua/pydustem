#-*- coding: utf-8 -*-
'''functions for fitting dust mass'''

from math import log10
import numpy as np
import scipy.optimize
import scipy.interpolate
import run

def counter(start=1):
    '''counter for repeatedly-executed DustEM'''
    i = start
    while True:
        yield i
        i += 1

def fitting_mass(bands, values, grain_composition, cache,
                 var_g0=True, var_gamma=False, var_umax=False, var_alpha=False,
                 variables=None, absolute=True, approx_g0=True):
    '''fitting dust mass (and other parameters) using scipy.optimize.fmin_slsqp'''
    if len(bands) != len(values): raise ValueError
    if not absolute and None not in variables.values():
        raise ValueError, "Need value 'None' for at least one dust component when calculating relative abundances."
    if not variables:
        variables = {}
        for i, grain in enumerate(grain_composition.grains.keys()):
            if absolute:
                variables[grain] = i+1
            else:
                if i == 0: variables[grain] = None
                else: variables[grain] = i

    cnt = counter()

    def func(x):
        print x
        ans = np.array(values)
        grain_input = grain_composition.copy()

        i = 0
        if var_g0:
            grain_input.g0 = grain_input.g0*0.1*10**x[i]
            i += 1
        if var_gamma:
            grain_input.gamma = grain_input.gamma*x[i]
            i += 1
        if var_umax:
            grain_input.umax = grain_input.umax*0.1*10**x[i]
            i += 1
        if var_alpha:
            grain_input.alpha = grain_input.alpha*x[i]
            i += 1

        for k, v in variables.items():
            if v: grain_input[k].mass *= x[v+i-1]
        
        x, y = run.calc_sed(grain_input, cache=cache, msg='running DUSTEM... '+str(cnt.next()))
        res = adjust_sed(bands, x, y)

        diff = ((ans-res)**2/ans).sum()/ans[0]
        print diff
        return diff

    def cons(x): return x[0]

    gnum = max(variables.values())
    x0 = np.ones((gnum+var_g0+var_gamma+var_umax+var_alpha,))
    bounds = []
    if var_g0: bounds.append(tuple(np.array([-4,6])-log10(grain_composition.g0)))
    if var_gamma: bounds.append((0,1./grain_composition.gamma))
    if var_umax: bounds.append(tuple(np.array([4,10])-log10(grain_composition.umax)))
    if var_alpha: bounds.append((1.01/grain_composition.alpha,4./grain_composition.alpha))

    x = scipy.optimize.fmin_slsqp(func, x0, bounds=bounds+[(1e-5,100)]*gnum, epsilon=1e-5)

    newgc = grain_composition.copy()
    
    i = 0
    if var_g0:
        print 'G0 = %f'%(grain_composition.g0*10**x[i]*0.1)
        newgc.g0 = grain_composition.g0*10**x[i]*0.1
        i += 1
    if var_gamma:
        print 'gamma = %f'%(grain_composition.gamma*x[i])
        newgc.gamma = grain_composition.gamma*x[i]
        i += 1
    if var_umax:
        print 'Umax = %f'%(grain_composition.umax*10**x[i]*1)
        newgc.umax = grain_composition.umax*10**x[i]*0.1
        i += 1
    if var_alpha:
        print 'alpha = %f'%(grain_composition.alpha*x[i])
        newgc.alpha = grain_composition.alpha*x[i]
        i += 1
    for k, v in variables.items():
        if v:
            print '%s: %E (x%f)'%(k,grain_composition[k].mass*x[v+i-1],x[v+i-1])
            newgc[k].mass = grain_composition[k].mass*x[v+i-1]
        else:
            print '%s: %E'%(k,grain_composition[k].mass)

    return newgc

def display_fitting_result(): pass

def adjust_sed(bands, x, y):
    '''calculate averaged/interpolated intensity'''
    if isinstance(bands[0], str):
        return filter.weighted_mean(bands, x, y)
    else:
        tck = scipy.interpolate.splrep(x, y)
        return scipy.interpolate.splev(bands, tck)
