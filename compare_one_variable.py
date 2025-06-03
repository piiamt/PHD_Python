'''
COMPARES SEVERAL SIMULATIONS' ONE VARIABLE OVER TIME
'''
import sys
import numpy as np
import matplotlib.pyplot as plt
import os, os.path
#from readts3 import tser3  # both combined
from matplotlib.colors import LogNorm
home = '/lustre/astro/piiamt/slow/'

#sims = np.array(['1.0AU_0.1Mearth','1.0AU_0.5Mearth', '1.0AU_1.0Mearth',
#                 '1.0AU_2.0Mearth','1.0AU_3.0Mearth', '1.0AU_4.0Mearth',
#                 '1.0AU_5.0Mearth'])

sims = np.loadtxt('../idl/simlist.txt', dtype='str')
folders=sims

plt.rcParams.update({'font.size'           : 7,#10, 
                     'mathtext.fontset'    : 'cm',
                     'font.family'         : 'serif',
                     'xtick.direction'     :'in',
                     'ytick.direction'     :'in',
                     #'axes.grid'           :True,
                     'grid.alpha'          : 0.7,
                     'grid.linestyle'      : 'dashed',
                     'axes.linewidth'      : 0.7,
                     'xtick.minor.width'   : 0.7,
                     'ytick.minor.width'   : 0.7,
                     'xtick.major.width'   : 0.7,
                     'ytick.major.width'   : 0.7
                    })

colours = np.array(['k', '#009e73', '#d55e00', '#cc79a7', '#0072b2', '#e69f00', '#56b4e9'])
RGBc = np.array([[0,0,0,0.2],[0,158/255,112/255,0.2],[213/255,94/255,0.0,0.2],[204/255,121/255,167/255,0.2],
                 [0,114/255,178/255,0.2],[230/255,159/255,0,0.2],[86/255,180/255,233/255,0.2]])

lines = np.array(['-','--','-.',':','-','--','-.',':'])
#################   CONSTANTS ############
mearth      = 5.9722e24     # kg
msun        = 1.98847e30    # kg
Myr_to_s    = 365.25*24*60*60*1000000

Aplas = np.array([])
Mplas = np.array([])
tenv90s = np.array([])
dragphis = np.array([])
ecolors = np.array([])
customs = np.array([])
fig, ax = plt.subplots(1,1,figsize=(3,3))
for folder in folders:
    timefile =home+folder+ '/time_series.idl'
    if not os.path.exists(timefile):
        folders = folders[np.where(folders!=folder)]
for folder in folders:
    f = folder
    if f[3]!='A':   # if Apla has 2 decimal places like 1.25AU
        Aplas = np.append(Aplas, float(f[:4]))
    elif f[3]=='A': # aka if Apla has 1 decimal place like 1.0AU
        Aplas = np.append(Aplas, float(f[:3]))
    Mplas = np.append(Mplas, np.load(home + f + '/pysave'+'/mpla.npy')[-1])
    ts = np.load(home + f + '/pysave/t.npy')
    Menv = np.load(home + f + '/pysave/menv.npy')
    Matm_O2 = np.load(home + f + '/pysave/matm_o2.npy')
    Matm_CO2 = np.load(home + f + '/pysave/matm_co2.npy')
    Matm = np.load(home + f + '/pysave/matm.npy')
    tmaxenv = ts[np.where(Menv==max(Menv))]
    if (Menv[-1]>(max(Menv)*0.1)):
        tenv90s = np.append(tenv90s, ts[-1])
    else:
        tenv90s = np.append(tenv90s, ts[np.where((ts>tmaxenv)&(Menv<(0.1*max(Menv))))][0])
    dragphi_O = np.load(home + f + '/pysave'+'/dragphi_O.npy')
#    dragphi_N = np.load(home + f + '/pysave'+'/dragphi_N.npy')
#    dragphi_C = np.load(home + f + '/pysave'+'/dragphi_C.npy')
    dragphis = np.append(dragphis, max(dragphi_O)/mearth)
    customs = np.append(customs,1)# max(Menv)/mearth)#max(Menv)/mearth)
    if max(dragphi_O)>=1:
        ecolors = np.append(ecolors, 'None')
    else: ecolors = np.append(ecolors, 'None')

at = ax.scatter(Aplas, Mplas/mearth, s=50, c=tenv90s/Myr_to_s, #customs
           edgecolor=ecolors, cmap='turbo')
ax.set_ylabel(r'Planet mass (M$_{\oplus}$)')
ax.set_xlabel('Orbital radius (AU)')
p = ax.get_position().get_points().flatten() #[[x0,y0],[x1,y1]]
at_cbar = fig.add_axes([p[0], p[3], p[2]-p[0], 0.015])
fig.colorbar(at, cax=at_cbar, orientation='horizontal',
            label=r'Time of envelope loss (Myr)')
at_cbar.xaxis.set_label_position('top')
at_cbar.xaxis.set_ticks_position('top')

plt.savefig(home+'/python/plots/notbugs_slow.pdf', bbox_inches='tight')

plt.show()

