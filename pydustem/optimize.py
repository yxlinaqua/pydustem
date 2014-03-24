#-*- coding: utf-8 -*-

from math import log10
import numpy as np
import scipy.optimize
import scipy.interpolate
import run
import grainmodel

def counter(start=1):
    i = start
    while True:
        yield i
        i += 1

def fitting_mass(bands, values, grain_composition, cache,
                 var_g0=True, var_gamma=False, var_umax=False, var_alpha=False,
                 variables=None, absolute=True, approx_g0=True):
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

#        print (np.log10(ans/res)**2).sum(), np.log10(ans/res)     
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
#    x = scipy.optimize.leastsq(func, x0, epsfcn=0.001)[0]
#    x = scipy.optimize.fmin_cobyla(func, x0, [lambda x, i=i: x[i] for i in xrange(len(x0))])
#    x = scipy.optimize.fmin_l_bfgs_b(func, x0, epsilon=1e-4, approx_grad=True,
#                                     bounds=bounds+[(1e-5,100)]*gnum)[0]
#    x = scipy.optimize.fmin_tnc(func, x0, epsilon=1e-4, approx_grad=True, maxfun=500,
#                                bounds=bounds+[(1e-5,100)]*gnum)
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
    if isinstance(bands[0], str):
        return filter.weighted_mean(bands, x, y)
    else:
        tck = scipy.interpolate.splrep(x, y)
        return scipy.interpolate.splev(bands, tck)

if __name__ == '__main__':
    cache = {}

#    bands = ['S11', 'S7', 'L15', 'L24', 'M70', 'M160']
#    values = np.array([0.369/11., 0.390/7., 0.619/15., 1.349/24., 10.3/70., 11.5/160.])
#    values /= 6e21

    aryS = np.loadtxt('/home/tnakamura/analysis/NGC2782/Spitzer-IRS/SL.dat')
    aryL = np.loadtxt('/home/tnakamura/analysis/NGC2782/Spitzer-IRS/LL.dat')
    aryM = np.loadtxt('/home/tnakamura/analysis/NGC2782/Spitzer-IRS/MIPS.dat')
    bands = np.r_[aryS[:,0],aryL[:,0],aryM[:,0]]
    values = np.r_[aryS[:,1]/aryS[:,0]/1.1,aryL[:,1]/aryL[:,0]/1.1,aryM[:,1]/aryM[:,0]]
    values /= 6e21
    print bands, values
    
    grain_composition = grainmodel.MC10(g0=10.0, gamma=0.01, alpha=2.0, umax=1e6,
                                  factor={'PAH0':1.0,'PAH1':1.0,'SamC':5.0,'LamC':1.0,'aSil':1.0})
#    grain_composition = dust.MC10(g0=1.0, gamma=0.01)
    var_g0 = True
    var_gamma = True
    var_umax = False
    var_alpha = False
#    variables = {'PAH0': 1, 'PAH1': 1, 'SamC': 2, 'LamC': None, 'aSil': None}
#    variables = {'PAH0': 1, 'PAH1': 1, 'SamC': 1, 'LamC': 2, 'aSil': 2}
    variables = {'PAH0': 1, 'PAH1': 1, 'SamC': 2, 'LamC': 3, 'aSil': 3}

    newgc = fitting_mass(bands, values, grain_composition, cache, var_g0, var_gamma, var_umax, var_alpha,
                         variables=variables, absolute=True)

#    x, y = run.calc_sed(newgc, cache=cache)
#    np.savetxt('test.dat', np.transpose(np.array([x, y])))https://github.com/dpforest/pydustem

    res = run.calc_sed(newgc, cache=cache, readraw=True)
    np.savetxt('test.dat', np.transpose(np.vstack((res['lambda'], res['PAH0']+res['PAH1'], res['SamC'], res['LamC']+res['aSil'], res['Total']))))