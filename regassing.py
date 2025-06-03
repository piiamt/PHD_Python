import sys
"""
Created on Feb 15 2023
@author: piiamt
PLOTS OLD ADAP SURFACE MAGMA OCEAN'S MELT FRACTION AND MATM UNDER THAT
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
home = '/lustre/astro/piiamt/dense'
folder1 = '/short_2.0Mearth'  

plt.rcParams.update({'font.size'           : 7,#10, 
                     'mathtext.fontset'    : 'cm',
                     'font.family'         : 'serif',
                     'xtick.direction'     :'out',
                     'ytick.direction'     :'out',
                     #'axes.grid'           :True,
                     'grid.alpha'          : 0.7,
                     'grid.linestyle'      : 'dashed',
                     'axes.linewidth'      : 0.7,
                     'xtick.minor.width'   : 0.7,
                     'ytick.minor.width'   : 0.7,
                     'xtick.major.width'   : 0.7,
                     'ytick.major.width'   : 0.7
                    })


colors = np.array(['k', '#009e73', '#d55e00', '#cc79a7', '#0072b2', '#e69f00', '#56b4e9'])
RGBc = np.array([[0,0,0,0.2],[0,158/255,112/255,0.2],[213/255,94/255,0.0,0.2],[204/255,121/255,167/255,0.2],
                 [0,114/255,178/255,0.2],[230/255,159/255,0,0.2],[86/255,180/255,233/255,0.2]])
# aka           black   green     vermillion    pink        blue    orange      lskdfs
#################   CONSTANTS ############
mearth      = 5.9722e24     # kg
msun        = 1.98847e30    # kg
Myr_to_s    = 365.25*24*60*60*1000000

varfile = home+folder1+ '/variables.idl'
var = list((readsav(varfile)).values())[0][0]
timef1 = readsav(home+folder1+'/time_series.idl')
tser1 = list(timef1.values())[0][0]

from matplotlib.colors import LogNorm


############### CONTOURPLOT ################
mm=1/25.4
from matplotlib.colors import TwoSlopeNorm

times1 = np.repeat(np.array([var.t]), len(var.r[0]), axis=0).T/Myr_to_s

rpla1 = tser1.rpla[-1]
rcor1 = tser1.rcor[-1]

colormap = 'turbo'#'jet' #'inferno' #'turbo'

#CHOOSE LOGARITHMIC INDEXES FOR FASTER PLOTTING
#timestotal = np.append(times[indexes], times2[indexes2], axis=0)
#rtotal = np.append(var.r[indexes], var2.r[indexes2], axis=0)
#pphitotal = np.append(var.pphi[indexes], var2.pphi[indexes2], axis=0)

fig, ax = plt.subplots(2,1, sharex=True, figsize=(88*mm,100*mm), 
        gridspec_kw=dict(hspace=0.00, wspace=0.05))
at1 = ax[0].contourf(times1, var.r/rpla1, var.pphi, 
        levels=100, cmap=colormap, vmin=0.0, vmax=1.0)

#print(var2.t/Myr_to_s)
#im1 = ax.contourf(times[indexes1], var.r[indexes1]/1000, var.pphi[indexes1], 
#        levels=100, interpolation='gaussian', cmap=colormap, vmin=0.0, vmax=1.0)
### FOR SMOOTH PDF CONTOUR EDGES, WARNING: SLOW
for c in at1.collections:
    c.set_edgecolor('face')

#ax[2].set_title(r'3.0 M$_{\oplus}$')

#ax[1].set_xticks([8.045,8.050,8.055,8.06])
#ax[1].set_xticklabels(['8.045','8.050','8.055','8.060'])

#ax[0].axhline(0.98, c='w', lw=0.7)
ax[0].axhline(rcor1/rpla1, c='w', lw=0.7)
#ax[1].axhline(0.98, c='w', lw=0.7)
#ax[2].axhline(0.98, c='w', lw=0.7)

#ax[2].text(31.329, 0.84, 'lid thickness limit\nto atmospheric\noutgassing', c='w')
#ax[2].text(31.329, 0.62, 'core-mantle-boundary', c='w')
ax[0].contour(times1, var.r/rpla1, var.pphi, 
              levels=[0.5], linewidths=0.5, colors='w')
ax[0].axvline(8.0665, c='w', lw=0.5)
ax[1].axvline(8.0665, c='k', lw=0.5, alpha=0.5)
ax[0].axvline(8.0672, c='w', lw=0.5)
ax[1].axvline(8.0672, c='k', lw=0.5, alpha=0.5)
#ax[0].axvline(5.0855, c='w', lw=0.5)
#ax[1].axvline(8.0604, c='w', lw=0.5)
#ax[2].axvline(31.3458, c='w', lw=0.5)

ax[1].plot(tser1.t/Myr_to_s, tser1.matm/mearth, c=colors[2],
        label='Atmosphere')
ax[1].plot(tser1.t/Myr_to_s, (tser1.msil_h2o+tser1.msil_co2+tser1.msil_n2)/mearth, 
        c=colors[0], label=r'Mantle volatiles')
ax[1].legend(frameon=False)

#ax.set_xscale('log')
#ax.legend()
cbar1 = fig.colorbar(at1, ax=ax, pad=0.00, ticks=np.array([0.0,0.2,0.4,0.6,0.8,1.0]))
cbar1.set_label('Melt fraction')
cbar1.ax.set_yticklabels(np.array(['0.0','0.2','0.4','0.6','0.8','1.0']))

#ax[0].set_xlim(5.6,5.8) # 0.1
#ax[2].set_xlim(15.403,15.426) # 3.0
#ax[0].set_xlim(5.07,5.092)
###ax[1].set_xlim(8.043, 8.085)
ax[1].set_xlim(8.046, 8.075)
#ax[2].set_xlim(31.328,31.35)
ax[0].set_ylim(0.6,1.02)
ax[1].set_ylim(1e-5,1.5e-3)
ax[1].set_yscale('log')

ax[0].set_ylabel('Radius (planet radii)')
ax[1].set_ylabel(r'Atmosphere mass (M$_{\oplus}$)')
ax[1].set_xlabel('Time (Myr)')
plt.savefig(home+'/python/plots/regassing.pdf', bbox_inches='tight')

plt.show()



