import sys
"""
Created on Feb 15 2023
@author: piiamt
CREATES 2 CONTOURPLOTS ONE UNDER THE OTHER, MAKE SURE THAT THE SIMULATION PLOTTED
HAS HAS THE LINE 
adap_read_var, obj=var, /all, /lsave
RUN IN IDL BEFORE RUNNING THIS PROGRAM
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
home = '/lunarc/nobackup/projects/astro4/piia/'
folder = '1.0AU_2.0Mearth'  #input('Data folder name:  ')


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

timefile =home+folder+ '/time_series.idl'
timef = readsav(timefile)
timekey = list(timef.values())[0].dtype
tser = list(timef.values())[0][0]
varfile = home+folder+ '/variables.idl'
var = list((readsav(varfile)).values())[0][0]

from matplotlib.colors import LogNorm
home = '/lunarc/nobackup/projects/astro4/piia/'


############### CONTOURPLOT ################
fig, ax = plt.subplots(2,1, sharex=True, figsize=(10,10*5/6), gridspec_kw=dict(hspace=0.05))
from matplotlib.colors import TwoSlopeNorm

times = np.repeat(np.array([var.t]), len(var.r[0]), axis=0).T/Myr_to_s

colormap = 'inferno'#'jet' #'inferno'

ax[0].set_title('1.0 AU, 1.5 Mearth')
im1 = ax[0].contourf(times, var.r/1000, var.pphi, levels=50, cmap=colormap, vmin=0.0, vmax=1.0)
# ax[0].plot(timeseries.t/Myr_to_s, )
im2 = ax[1].contourf(times, var.r/1000, var.tt, levels=np.arange(-2,max(var.tt.flatten())//100*100+200, 50), cmap=colormap)
ax[1].plot(tser.t/Myr_to_s, tser.rsolid_cor/1000, ls='-', c='w', label='Solid core')
#ax[1].plot(tser.t/Myr_to_s, tser.rpla/1000, ls=':', c='w', label='Surface')
#ax[1].plot(tser.t/Myr_to_s, tser.rmagma_cor/1000, ls='-.', c='w', label='Magma core')
#ax[1].plot(tser.t/Myr_to_s, tser.rmagma_sur/1000, ls='', c='w', label='Magma surface')

ax[1].legend()
cbar1 = fig.colorbar(im1, ax=ax[0], pad=0.01)
cbar1.set_label('Melt fraction')
cbar2 = fig.colorbar(im2, ax=ax[1], pad=0.01)
cbar2.set_label('Temperature (K)')
# plt.colorbar()
ax[0].set_ylim(0,9500)
ax[1].set_ylim(0,9500)
ax[0].set_xlim(0,300)
ax[1].set_xlim(0,300)

ax[1].set_xlabel('Time (Myr)')
ax[0].set_ylabel('Radius (km)')
ax[1].set_ylabel('Radius (km)')

#plt.savefig(home+'python/plots/contour_magma_long.png', bbox_inches='tight')
plt.savefig(home+'python/plots/contour_magma_short.png', bbox_inches='tight')
#plt.savefig('1.1Mearth_300Myr_contourRHO_T.pdf', bbox_inches='tight')

plt.show()



