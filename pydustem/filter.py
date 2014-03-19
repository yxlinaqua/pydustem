#-*- coding: utf-8 -*-

import os
import numpy as np
from scipy import interpolate

PATH = os.path.join(os.path.abspath(os.path.dirname(__file__)),'filterdata')
FILE = {
         'N2'  : 'n2.fad',  'N3'  : 'n3.fad',   'N4'  : 'n4.fad',    # NIR
         'S7'  : 's7.fad',  'S9W' : 's9w.fad',  'S11' : 's11.fad',   # MIR-S
         'L15' : 'l15.fad', 'L18W': 'l18w.fad', 'L24' : 'l24.fad',   # MIR-L
         'M24' : 'm24.fad', 'M70' : 'm70.fad',  'M160': 'm160.fad',  # MIPS
       }


def response_curve(filter_name, dtype='array', unpack=True, **kwargs):
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
    ''' Get the filter path'''
    return os.path.join(PATH,FILE[filter_name])

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
