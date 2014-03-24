#-*- coding: utf-8 -*-

import os
import sys
import copy
from collections import OrderedDict
from math import log10
import numpy as np
from config import RES_FILE, GRAIN, DUSTEM, RES_HEADER

def run_dustem(grain_composition, output='SED', msg=None, silent=True):
    if not msg: msg = 'running DUSTEM...'
    grain_composition.correct_params()
    grain_composition.savefile(GRAIN)
    print >> sys.stderr, msg
    if silent:
        os.system('%s > /dev/null'%DUSTEM)
    else:
        os.system(DUSTEM)

def calc_sed(gc, du=0.5, msg=None, silent=True, readraw=False, cache=None):
    if gc.gamma:
        # ulist = [g0, 10**(k_1), ..., 10**(k_n), umax]
        ulist = 10**np.arange(int(log10(gc.g0)), int(log10(gc.umax))+1, du)
        ulist = np.r_[gc.g0, ulist[(ulist>gc.g0)&(ulist<gc.umax)], gc.umax]
    else:
        ulist = [gc.g0]
    
    gclist = [gc.dup(g0=i, gamma=0.0) for i in ulist]

    if cache is None:
        jlist = map(calc_single_sed, gclist)
    else:
        jlist = []
        for gcsub in gclist:
            if '%E'%gcsub.g0 in cache:
                jlist.append(calc_sed_from_cache(gcsub, cache))
            else:
                res = calc_single_sed(gcsub)
                jlist.append(copy.deepcopy(res))
                for k in res.keys():
                    if k in gcsub.grains.keys(): res[k] /= gcsub.grains[k].mass
                cache['%E'%gcsub.g0] = res

    res = integ_sed_gamma(gc, ulist, jlist)
    return res if readraw else (res['lambda'], res['Total'])

def calc_single_sed(gc, msg=None, silent=True, readraw=True):
    run_dustem(gc, output='SED', msg=msg, silent=silent)
    ary = read_output_raw(RES_FILE, skiprows=RES_HEADER, unpack=True)
    ylist = ['lambda'] + gc.grains.keys() + ['Total']
    res = OrderedDict(zip(ylist, ary[:]))

    return res if readraw else (res['lambda'], res['Total'])

def calc_sed_from_cache(grain_composition, cache):
    res = copy.deepcopy(cache['%E'%grain_composition.g0])
    for k in grain_composition.grains.keys():
        res[k] *= grain_composition.grains[k].mass
    return res

def integ_sed_gamma(gc, ulist, jlist):
    res = OrderedDict({'lambda': jlist[0]['lambda']})
    res['Total'] = np.zeros(jlist[0]['lambda'].shape)
    for k in gc.grains.keys():
        res[k] = jlist[0][k]*(1-gc.gamma)

    for i in xrange(len(jlist)-1):
        coef = gc.gamma*(gc.alpha-1)/(gc.g0**(1-gc.alpha)-gc.umax**(1-gc.alpha))

        for k in gc.grains.keys():
            j0 = jlist[i][k]*ulist[i]**(-gc.alpha)*coef
            j1 = jlist[i+1][k]*ulist[i+1]**(-gc.alpha)*coef
            res[k] += np.sqrt(j0*j1)*(ulist[i+1]-ulist[i])

    for k in gc.grains.keys():
        res['Total'] += res[k]

    return res

def read_output(filename, skiprows=0, xcol=0, ycol=-1, unpack=True, reverse=False):
    x, y =  np.loadtxt(filename, skiprows=skiprows, usecols=(xcol,ycol), unpack=unpack)
    if reverse:
        x = x[::-1]
        y = y[::-1]
    return (x,y)

def read_output_raw(filename, skiprows=0, unpack=False):
    return np.loadtxt(filename, skiprows=skiprows, unpack=unpack)

