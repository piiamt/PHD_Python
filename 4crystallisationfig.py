import sys
"""
Created on Feb 15 2023
@author: piiamt
CREATES 3 CONTOURPLOTs, MAKE SURE THAT THE SIMULATION PLOTTED
HAS HAS THE LINE 
adap_read_var, obj=var, /all, /lsave
RUN IN IDL BEFORE RUNNING THIS PROGRAMi
WITH LOG TIME THE SIM MUST BE RUN WITH A LOT OF VAR FILES, 100 GIVES BAD PLOT
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
home = '/lustre/astro/piiamt/locean_outgas'
folder1 = '/1.0AU_0.1Mearth'  
folder2 = '/1.0AU_0.5Mearth'   
folder3 = '/1.0AU_1.0Mearth'
#folder1 = '/1.0AU_1.5Mearth'  
folder4 = '/1.0AU_2.0Mearth'   
#folder3 = '/1.0AU_3.0Mearth'

plt.rcParams.update({'font.size'           : 5.5,#10, 
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
varfile2 = home+folder2+'/variables.idl'
var2 = list((readsav(varfile2)).values())[0][0]
varfile3 = home+folder3+'/variables.idl'
var3 = list((readsav(varfile3)).values())[0][0]
varfile4 = home+folder4+'/variables.idl'
var4 = list((readsav(varfile4)).values())[0][0]
timef1 = readsav(home+folder1+'/time_series.idl')
tser1 = list(timef1.values())[0][0]
timef2 = readsav(home+folder2+'/time_series.idl')
tser2 = list(timef2.values())[0][0]
timef3 = readsav(home+folder3+'/time_series.idl')
tser3 = list(timef3.values())[0][0]
timef4 = readsav(home+folder4+'/time_series.idl')
tser4 = list(timef4.values())[0][0]



from matplotlib.colors import LogNorm


############### CONTOURPLOT ################
mm=1/25.4
from matplotlib.colors import TwoSlopeNorm

times1 = np.repeat(np.array([var.t]), len(var.r[0]), axis=0).T/Myr_to_s
times2 = np.repeat(np.array([var2.t]), len(var2.r[0]), axis=0).T/Myr_to_s
times3 = np.repeat(np.array([var3.t]), len(var3.r[0]), axis=0).T/Myr_to_s
times4 = np.repeat(np.array([var4.t]), len(var4.r[0]), axis=0).T/Myr_to_s

rpla1 = tser1.rpla[-1]
rcor1 = tser1.rcor[-1]
rpla2 = tser2.rpla[-1]
rcor2 = tser2.rcor[-1]
rpla3 = tser3.rpla[-1]
rcor3 = tser3.rcor[-1]
rpla4 = tser4.rpla[-1]
rcor4 = tser4.rcor[-1]

colormap = 'turbo'#'jet' #'inferno' #'turbo'

#CHOOSE LOGARITHMIC INDEXES FOR FASTER PLOTTING
#timestotal = np.append(times[indexes], times2[indexes2], axis=0)
#rtotal = np.append(var.r[indexes], var2.r[indexes2], axis=0)
#pphitotal = np.append(var.pphi[indexes], var2.pphi[indexes2], axis=0)

fig, ax = plt.subplots(1,4, sharey=True, figsize=(190*mm,65*mm), gridspec_kw=dict(hspace=0.24, wspace=0.05))
at1 = ax[0].contourf(times1, var.r/rpla1, var.pphi, 
        levels=100, cmap=colormap, vmin=0.0, vmax=1.0)
at2 = ax[1].contourf(times2, var2.r/rpla2, var2.pphi, 
        levels=100, cmap=colormap, vmin=0.0, vmax=1.0)
at3 = ax[2].contourf(times3, var3.r/rpla3, var3.pphi, 
        levels=100, cmap=colormap, vmin=0.0, vmax=1.0)
at4 = ax[3].contourf(times4, var4.r/rpla4, var4.pphi, 
        levels=100, cmap=colormap, vmin=0.0, vmax=1.0)

#print(var2.t/Myr_to_s)
#im1 = ax.contourf(times[indexes1], var.r[indexes1]/1000, var.pphi[indexes1], 
#        levels=100, interpolation='gaussian', cmap=colormap, vmin=0.0, vmax=1.0)
### FOR SMOOTH PDF CONTOUR EDGES, WARNING: SLOW
for c in at1.collections:
    c.set_edgecolor('face')
for c in at2.collections:
    c.set_edgecolor('face')
for c in at3.collections:
    c.set_edgecolor('face')
for c in at4.collections:
    c.set_edgecolor('face')

ax[0].set_title(r'0.1 M$_{\oplus}$')
ax[1].set_title(r'0.5 M$_{\oplus}$')
ax[2].set_title(r'1.0 M$_{\oplus}$')
ax[3].set_title(r'2.0 M$_{\oplus}$')
#ax[0].set_title(r'1.5 M$_{\oplus}$')
#ax[1].set_title(r'2.0 M$_{\oplus}$')
#ax[2].set_title(r'3.0 M$_{\oplus}$')

#ax[1].set_xticks([8.045,8.050,8.055,8.06])
#ax[1].set_xticklabels(['8.045','8.050','8.055','8.060'])

#ax[0].axhline(0.98, c='w', lw=0.7)
ax[0].axhline(rcor1/rpla1, c='w', lw=0.7)
#ax[1].axhline(0.98, c='w', lw=0.7)
ax[1].axhline(rcor2/rpla2, c='w', lw=0.7)
#ax[2].axhline(0.98, c='w', lw=0.7)
ax[2].axhline(rcor3/rpla3, c='w', lw=0.7)
ax[3].axhline(rcor4/rpla4, c='w', lw=0.7)

#ax[2].text(31.329, 0.84, 'lid thickness limit\nto atmospheric\noutgassing', c='w')
#ax[2].text(31.329, 0.62, 'core-mantle-boundary', c='w')
ax[0].contour(times1, var.r/rpla1, var.pphi, 
              levels=[0.5], linewidths=0.5, colors='w')
ax[1].contour(times2, var2.r/rpla2, var2.pphi, 
              levels=[0.5], linewidths=0.5, colors='w')
ax[2].contour(times3, var3.r/rpla3, var3.pphi, 
              levels=[0.5], linewidths=0.5, colors='w')
ax[3].contour(times4, var4.r/rpla4, var4.pphi, 
              levels=[0.5], linewidths=0.5, colors='w')

#ax[0].axvline(5.0855, c='w', lw=0.5)
#ax[1].axvline(8.0604, c='w', lw=0.5)
#ax[2].axvline(31.3458, c='w', lw=0.5)

#ax.set_xscale('log')
#ax.legend()
cbar1 = fig.colorbar(at3, ax=ax, pad=0.005, ticks=np.array([0.0,0.2,0.4,0.6,0.8,1.0]))
cbar1.set_label('Melt fraction')
cbar1.ax.set_yticklabels(np.array(['0.0','0.2','0.4','0.6','0.8','1.0']))

ax[0].set_xlim(5.6,5.8) # 0.1
ax[1].set_xlim(5.07,5.13) # 0.5
ax[2].set_xlim(5.435,5.455) # 1.0
#ax[0].set_xlim(6.326,6.346) # 1.5
ax[3].set_xlim(8.046,8.068) # 2.0
#ax[2].set_xlim(15.403,15.426) # 3.0
#ax[0].set_xlim(5.07,5.092)
#ax[1].set_xlim(8.043, 8.065)
#ax[2].set_xlim(31.328,31.35)
ax[0].set_ylim(0.6,1.05)

ax[0].set_ylabel('Radius (planet radii)')
ax[0].set_xlabel('Time (Myr)')
ax[1].set_xlabel('Time (Myr)')
ax[2].set_xlabel('Time (Myr)')
ax[3].set_xlabel('Time (Myr)')
plt.savefig(home+'/python/plots/4crystallisationfig_locean_outgas.pdf', bbox_inches='tight')

plt.show()



