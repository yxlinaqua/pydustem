import os

'''DustEM preferences'''
DUSTEM_DIR = '/home/tnakamura/workspace/dustem3.8_web/'
OUT_DIR = os.path.join(DUSTEM_DIR, 'out')
RES_FILE = os.path.join(OUT_DIR, 'SED.RES')
GRAIN = os.path.join(DUSTEM_DIR, 'data/GRAIN.DAT')
DUSTEM = os.path.join(DUSTEM_DIR, 'src/dustem')

''''DustEM output'''
RES_HEADER = 8 # lines

'''Filter data'''
FILTER_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)),'filterdata')
FILTER_PATH = {
                 'N2'  : 'n2.fad',  'N3'  : 'n3.fad',   'N4'  : 'n4.fad',    # AKARI/IRC NIR
                 'S7'  : 's7.fad',  'S9W' : 's9w.fad',  'S11' : 's11.fad',   # AKARI/IRC MIR-S
                 'L15' : 'l15.fad', 'L18W': 'l18w.fad', 'L24' : 'l24.fad',   # AKARI/IRC MIR-L
                 'M24' : 'm24.fad', 'M70' : 'm70.fad',  'M160': 'm160.fad',  # Spitzer/MIPS
              }

'''Grain parameters'''
DEFAULT_UMAX = 1e6
DEFAULT_GAMMA = 0.0
DEFAULT_ALPHA = 2.0
MASS_MIN = 1e-15
G0_MIN = 0.0001
G0_MAX = 1e10
