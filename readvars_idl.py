#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 15:17:07 2022

@author: piiamt
"""
from scipy.io import readsav
import numpy as np
import itertools
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os, os.path
from os.path import expanduser
home = expanduser('~')
folder = '/' + '1.1Mearth_cdth0.01'  #input('Data folder name:  ')


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

varfile = home+folder+ '/variables.idl'
var = list((readsav(varfile)).values())[0][0]
nr = np.array(var.nr)
nt = np.array(var.nt)
t = np.array(var.t)
it = np.array(var.it)
imax = np.array(var.imax)
cname = np.array(var.cname)
r = np.array([var.r[0]])
iicomp = np.array([var.iicomp[0]])
tt = np.array([var.tt[0]])
rrho = np.array([var.rrho[0]])
pp = np.array([var.pp[0]])
ee = np.array([var.ee[0]])
ccp = np.array([var.ccp[0]])
kkeff = np.array([var.kkeff[0]])
pphi = np.array([var.pphi[0]])
aa = np.array([var.aa[0]])
ddm = np.array([var.ddm[0]])
ddv = np.array([var.ddv[0]])
aalpha = np.array([var.aalpha[0]])
ggammaad = np.array([var.ggammaad[0]])
eeta = np.array([var.eeta[0]])
rra = np.array([var.rra[0]])
nnu = np.array([var.nnu[0]])
ttmax = np.array([var.ttmax[0]])
ttsol = np.array([var.ttsol[0]])
ttliq = np.array([var.ttliq[0]])
tte = np.array([var.tte[0]])
ppe = np.array([var.ppe[0]])
rrhoe = np.array([var.rrhoe[0]])
wwe_vap = np.array([var.wwe_vap[0]])
wwe_liq = np.array([var.wwe_liq[0]])
ppe_sat = np.array([var.ppe_sat[0]])


filemax = len(os.listdir(home+folder+'/data'))-4
for no in np.arange(1, filemax+1):
    r = np.append(r, np.array([var.r[no]]), axis=0)
    iicomp = np.append(iicomp, np.array([var.iicomp[no]]), axis=0)
    tt = np.append(tt, np.array([var.tt[no]]), axis=0)
    rrho = np.append(rrho, np.array([var.rrho[no]]), axis=0)
    pp = np.append(pp, np.array([var.pp[no]]), axis=0)
    ee = np.append(ee, np.array([var.ee[no]]), axis=0)
    ccp = np.append(ccp, np.array([var.ccp[no]]), axis=0)
    kkeff = np.append(kkeff, np.array([var.kkeff[no]]), axis=0)
    pphi = np.append(pphi, np.array([var.pphi[no]]), axis=0)
    aa = np.append(aa, np.array([var.aa[no]]), axis=0)
    ddm = np.append(ddm, np.array([var.ddm[no]]), axis=0)
    ddv = np.append(ddv, np.array([var.ddv[no]]), axis=0)
    aalpha = np.append(aalpha, np.array([var.aalpha[no]]), axis=0)
    ggammaad = np.append(ggammaad, np.array([var.ggammaad[no]]), axis=0)
    eeta = np.append(eeta, np.array([var.eeta[no]]), axis=0)
    rra = np.append(rra, np.array([var.rra[no]]), axis=0)
    nnu = np.append(nnu, np.array([var.nnu[no]]), axis=0)
    ttmax = np.append(ttmax, np.array([var.ttmax[no]]), axis=0)
    ttsol = np.append(ttsol, np.array([var.ttsol[no]]), axis=0)
    ttliq = np.append(ttliq, np.array([var.ttliq[no]]), axis=0)
    tte = np.append(tte, np.array([var.tte[no]]), axis=0)
    ppe = np.append(ppe, np.array([var.ppe[no]]), axis=0)
    rrhoe = np.append(rrhoe, np.array([var.rrhoe[no]]), axis=0)
    wwe_vap = np.append(wwe_vap, np.array([var.wwe_vap[no]]), axis=0)
    wwe_liq = np.append(wwe_liq, np.array([var.wwe_liq[no]]), axis=0)
    ppe_sat = np.append(ppe_sat, np.array([var.ppe_sat[no]]), axis=0)




