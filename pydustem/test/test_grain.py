import os
import numpy as np
from pydustem import run, grainmodel

def test_MC10():
    MC10_REF = os.path.join(os.path.dirname(__file__), 'data/SED_MC10.RES')
    with open(MC10_REF) as f:
        names = f.readline().split()
    grain_composition = grainmodel.MC10()
    result = run.calc_single_sed(grain_composition, readraw=True)
    answer = np.loadtxt(MC10_REF,skiprows=1)
    for i, name in enumerate(names):
        assert ((result[name]-answer[:,i]).sum()<1e-15).all()

def test_DL07():
    DL07_REF = os.path.join(os.path.dirname(__file__), 'data/SED_DL07.RES')
    with open(DL07_REF) as f:
        names = f.readline().split()
    grain_composition = grainmodel.DL07()
    result = run.calc_single_sed(grain_composition, readraw=True)
    answer = np.loadtxt(DL07_REF,skiprows=1)
    for i, name in enumerate(names):
        assert ((result[name]-answer[:,i]).sum()<1e-15).all()
