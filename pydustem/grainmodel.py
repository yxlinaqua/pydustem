#-*- coding: utf-8 -*-

from collections import OrderedDict
from grain import Grain, GrainComposition

def MC10_PAH0(mass=7.8000E-04):
    return Grain('PAH0', 10, ('mix','logn'), mass, 2.2400E+00, 3.5000E-08, 1.2000E-07,
                 6.4000E-08, 1.0000E-01)
def MC10_PAH1(mass=7.8000E-04):
    return Grain('PAH1', 10, ('mix','logn'), mass, 2.2400E+00, 3.5000E-08, 1.2000E-07,
                 6.4000E-08, 1.0000E-01)
def MC10_SamC(mass=1.6500E-04):
    return Grain('amCBEx', 15, ('logn',), mass, 1.8100E+00, 6.0000E-08, 2.0000E-06,
                 2.0000E-07, 3.5000E-01)
def MC10_LamC(mass=1.4500E-03):
    return Grain('amCBEx', 25, ('plaw','ed'), mass, 1.8100E+00, 4.0000E-07, 2.0000E-04,
                 -2.8000E+00, 1.5000E-05, 1.5000E-05, 2.0000E+00)
def MC10_aSil(mass=7.8000E-03):
    return Grain('aSil', 25, ('plaw','ed'), mass, 3.5000E+00, 4.0000E-07, 2.0000E-04,
                 -3.4000E+00, 2.0000E-05, 2.0000E-05, 2.0000E+00)

def MC10(g0=1.000, factor=None, umax=1e6, gamma=0.0, alpha=2.0):
    mass = {'PAH0': 7.8000E-04, 'PAH1': 7.8000E-04, 'SamC': 1.6500E-04, 'LamC': 1.4500E-03, 'aSil': 7.8000E-03}

    if factor:
        for item in mass.keys():
            if item in factor: mass[item] *= factor[item]

    return GrainComposition('sdist', g0,
                            OrderedDict({
                              'PAH0': MC10_PAH0(mass['PAH0']),
                              'PAH1': MC10_PAH1(mass['PAH1']),
                              'SamC': MC10_SamC(mass['SamC']),
                              'LamC': MC10_LamC(mass['LamC']),
                              'aSil': MC10_aSil(mass['aSil'])
                            }),
                            umax, gamma, alpha)

def DL07_PAH0(mass=5.4000E-04):
    return Grain('PAH0_DL07', 10, ('mix','logn'), mass, 2.2400E+00, 3.1000E-08, 1.2000E-07,
                 3.5000E-08, 4.0000E-01)
def DL07_PAH1(mass=5.4000E-04):
    return Grain('PAH1_DL07', 10, ('mix','logn'), mass, 2.2400E+00, 3.1000E-08, 1.2000E-07,
                 3.5000E-08, 4.0000E-01)
def DL07_SGra(mass=1.8000E-04):
    return Grain('Gra', 30, ('logn',), mass, 2.2400E+00, 3.1000E-08, 4.0000E-06,
                 3.0000E-07, 4.0000E-01)
def DL07_LGra(mass=2.3300E-03):
    return Grain('Gra', 70, ('plaw','ed', 'cv'), mass, 2.2400E+00, 3.1000E-08, 2.0000E-04,
                 -2.5400E+00, 1.0700E-06, 4.2800E-05, 3.0000E+00, 1.0700E-06, -1.6500E-01, 1.0000E+00)
def DL07_aSil(mass=8.2700E-03):
    return Grain('aSil', 70, ('plaw','ed', 'cv'), mass, 3.5000E+00, 3.1000E-08, 2.0000E-04,
                 -3.2100E+00, 1.6400E-05, 1.0000E-05, 3.0000E+00, 1.6400E-05, 3.0000E-01, 1.0000E+00)

def DL07(g0=1.000, factor=None, umax=1e6, gamma=0.0, alpha=2.0):
    mass = {'PAH0': 5.4000E-04, 'PAH1': 5.4000E-04, 'SGra': 1.8000E-04, 'LGra': 2.3300E-03, 'aSil': 8.2700E-03}

    if factor:
        for item in mass.keys():
            if item in factor: mass[item] *= factor[item]

    return GrainComposition('sdist', g0,
                            OrderedDict({
                              'PAH0': DL07_PAH0(mass['PAH0']),
                              'PAH1': DL07_PAH1(mass['PAH1']),
                              'SGra': DL07_SGra(mass['SGra']),
                              'LGra': DL07_LGra(mass['LGra']),
                              'aSil': DL07_aSil(mass['aSil'])
                            }),
                            umax, gamma, alpha)
