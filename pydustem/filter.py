#-*- coding: utf-8 -*-
'''functions for filter convolution'''

import os
import numpy as np
from scipy import interpolate
from config import FILTER_DIR, FILTER_PATH

def response_curve(filter_name, dtype='array', unpack=True, **kwargs):
    '''Get a response curve of a filter'''
    ary = np.loadtxt(get_path(filter_name))
    if dtype == 'spline':
        return interpolate.splrep(ary[:,0], ary[:,1], **kwargs)
    elif dtype == 'array':
        if unpack:
            return np.transpose(ary)
        else:
            return ary
        
    else:
        raise TypeError, 'Returned type should be "spline" or "array".'
    
def get_path(filter_name):
    '''Get a filter path'''
    return os.path.join(FILTER_DIR,FILTER_PATH[filter_name])

def weighted_mean(filter_list, x, y):
    ''' Calculate weighted mean values of the spectrum for the selected filters'''
    res = []
    for filter_name in filter_list:
        lmb, response = response_curve(filter_name)
        tck = interpolate.splrep(x, y)
        y2 = interpolate.splev(lmb, tck)
        y_weighted = y2 * response
        res.append(np.trapz(y_weighted, lmb)/np.trapz(response, lmb))
    return np.array(res)
