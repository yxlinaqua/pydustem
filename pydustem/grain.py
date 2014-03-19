#-*- coding: utf-8 -*-

import copy
from collections import OrderedDict

DEFAULT_UMAX = 1e6
DEFAULT_GAMMA = 0.0
DEFAULT_ALPHA = 2.0
MASS_MIN = 1e-15
G0_MIN = 0.0001
G0_MAX = 1e10

class Grain(object):
    '''Class for Grain'''
    def __init__(self, gtype, nsize, keywords, mass, rho, amin, amax, *args):
        self.gtype = gtype
        self.nsize = nsize
        self.keywords = keywords
        self.mass = mass
        self.rho = rho
        self.amin = amin
        self.amax = amax
        args = list(args)

        if 'logn' in keywords:
            self.a0 = args.pop(0)
            self.sigma = args.pop(0)
        elif 'plaw' in keywords:
            self.alpha = args.pop(0)
            if 'ed' in keywords:
                self.at = args.pop(0)
                self.ac = args.pop(0)
                self.gamma = args.pop(0)
            if 'cv' in keywords:
                self.au = args.pop(0)
                self.zeta = args.pop(0)
                self.eta = args.pop(0)
        else:
            raise ValueError

        if len(args) > 0: raise ValueError

    
    def __str__(self):
        '''make a line of formatted text like GRAIN.DAT'''
        if self.mass == 0.0: return ''
        if 'logn' in self.keywords:
            values = '%.4E %.4E'%(self.a0, self.sigma)
        elif 'plaw' in self.keywords:
            values = '%.4E'%self.alpha
            if 'ed' in self.keywords:
                values += ' %.4E %.4E %.4E'%(self.at, self.ac, self.gamma)
            if 'cv' in self.keywords:
                values += ' %.4E %.4E %.4E'%(self.au, self.zeta, self.eta)

        return '%s %d %s %.4E %.4E %.4E %.4E %s'%(self.gtype, self.nsize, '-'.join(self.keywords),
                                                  self.mass, self.rho, self.amin, self.amax, values)

    def correct_params(self):
        '''correct the mass if it's too small'''
        if self.mass < MASS_MIN: self.mass = MASS_MIN

    def copy(self):
        '''copy an instance of Grain'''
        return copy.deepcopy(self)
    
    def dup(self, **kwargs):
        '''copy an instance of Grain and change some parameters'''
        new = copy.deepcopy(self)
        for k,v in kwargs.items():
            setattr(new, k, v)
        return new


class GrainComposition(object):
    '''Class for Grain composition'''
    def __init__(self, run, g0, grains, umax=DEFAULT_UMAX, gamma=DEFAULT_GAMMA, alpha=DEFAULT_ALPHA):
        self.run = run
        self.g0 = g0
        self.gamma = gamma
        self.umax = umax
        self.alpha = alpha

        if isinstance(grains, dict):
            self.grains = OrderedDict(grains)
        else:
            self.grains = OrderedDict()
            for i, grain in enumerate(grains):
                self.grains[i] = grain

    def __mul__(self, num):
        new = self.copy()
        for grain in new.grains.keys():
            new.grains[grain].mass *= num
        return new
        
    def __getitem__(self, key):
        return self.grains[key]

    def __str__(self):
        '''make a formatted text like GRAIN.DAT'''
        string = '%s\n%.5f\n'%(self.run, self.g0)
        for grain in self.grains.values():
            s = str(grain)
            if not s=='': string += '%s\n'%s
        return string
    
    def correct_params(self):
        '''correct too small/large values'''
        if self.g0 < G0_MIN: self.g0 = G0_MIN
        if self.g0 > G0_MAX: self.g0 = G0_MAX
        for i in self.grains.keys():
            self.grains[i].correct_params()
        
    def copy(self):
        '''copy an instance of GrainComposition'''
        return copy.deepcopy(self)

    def dup(self, **kwargs):
        '''copy an instance of GrainComposition and change some parameters'''
        new = copy.deepcopy(self)
        for k,v in kwargs.items():
            setattr(new, k, v)
        return new
        
    def savefile(self, path):
        '''save the grain composition data like GRAIN.DAT'''
        f = open(path, 'w')
        f.write(str(self))
        f.close()

