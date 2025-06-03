#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 15:17:07 2022

@author: piiamt
"""
from scipy.io import readsav
import numpy as np
import itertools
#import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os, os.path
from os.path import expanduser
home = '/lunarc/nobackup/projects/astro4/piia'
#'/projects/astro3/nobackup/piia' #expanduser('~')
folder = '/' + '1.0AU_1.5Mearth'#'1.5AU_1.0Mearth'  # + input('Data folder name:  ')
# 0.9AU_1.2Mearth
#1.2Mearth_300Myr_stable
# stable_1.1Mearth_300Myr
plt.rcParams.update({'font.size'           : 12, 
                     'mathtext.fontset'    : 'cm',
                     'font.family'         : 'serif',
                     'xtick.direction'     :'in',
                     'ytick.direction'     :'in',
                     #'axes.grid'           :True,
                     'grid.alpha'          : 0.7,
                     'grid.linestyle'      : 'dashed'
                    })

colors = np.array(['k', '#009e73', '#d55e00', '#cc79a7', '#0072b2', '#e69f00', '#56b4e9'])
RGBc = np.array([[0,0,0,0.2],[0,158/255,112/255,0.2],[213/255,94/255,0.0,0.2],[204/255,121/255,167/255,0.2],
                 [0,114/255,178/255,0.2],[230/255,159/255,0,0.2],[86/255,180/255,233/255,0.2]])
# aka           black   green     vermillion    pink        blue    orange      lskdfs
#################   CONSTANTS ############
mearth      = 5.9722e24     # kg
msun        = 1.98847e30    # kg
Myr_to_s    = 365.25*24*60*60*1000000
#%%
timefile =home+folder+ '/time_series.idl'
timef = readsav(timefile)
timekey = list(timef.values())[0].dtype
tser = list(timef.values())[0][0]




