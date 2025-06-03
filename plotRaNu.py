'''
plots a 2x3 figure of different lid sizes and planet masses affecting
convection timescale tau, Ra and Nu
'''

import numpy as np
from scipy.io import readsav
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os, os.path
home = '/lustre/astro/piiamt/locean_outgas/'
#folder = 'newnewnew_1.5Mearth_0.02'
folders = np.array(['pysave'])

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
colors = np.array(['k', '#009e73', '#d55e00', '#cc79a7', '#0072b2', '#e69f00', '#56b4e9'])
RGBc = np.array([[0,0,0,0.2],[0,158/255,112/255,0.2],[213/255,94/255,0.0,0.2],[204/255,121/255,167/255,0.2],
                 [0,114/255,178/255,0.2],[230/255,159/255,0,0.2],[86/255,180/255,233/255,0.2]])
### CONSTANTS ###
yr          = 365.25*24*60*60
Myr         = 365.25*24*60*60*1e6
AU          = 1.49e11 		# m
apla        = 1.5*AU
Myr         = 3.15576e13	# s
G           = 6.6743e-11	# N*m^2 / kg^2
kB          = 1.38e-23
rhopla      = 5343.7498 	# kg/m^3
me          = 5.9722e24		# kg

############### READING IN DATA
Mplas = np.load(home + folders[0] + '/Mplas.npy')
taus = np.load(home + folders[0] + '/taus.npy')
tausurs = np.load(home + folders[0] + '/tausurs.npy')
Ras = np.load(home + folders[0] + '/Ras.npy')
Nus = np.load(home + folders[0] + '/Nus.npy')
dTdrs = np.load(home + folders[0] + '/dTdrs.npy')
tsaves = np.load(home + folders[0] + '/tsaves.npy')
DTs = np.load(home + folders[0] + '/DTs.npy')

################# PLOTTING
from matplotlib.ticker import FormatStrFormatter
import matplotlib.ticker as tick
def y_fmt(x,y):
    if x==0:
        return 0
    return '{0:.1f}e{1:.0f}'.format(np.sign(x)*10**(-np.floor(np.log10(abs(x)))+np.log10(abs(x))), np.floor(np.log10(abs(x)))).replace('e',r'$\times 10^{')+'}$'

mm=1/25.4
fig, ax = plt.subplots(3,1, sharex=True, figsize=(77*mm,120*mm), gridspec_kw=dict(hspace=0.0, wspace=0.0))

ax[2].set_xlabel('Planet mass')
ax[0].set_ylabel(r'Magma Ra')
ax[1].set_ylabel(r'Magma Nu')
ax[2].set_ylabel(r'$\Delta T$')

ax[0].plot(Mplas/me, Ras, 
        c=colors[0], label=r'0.05R lid')

ax[1].plot(Mplas/me, Nus, 
        c=colors[1], label=r'0.05R lid')
ax[2].plot(Mplas/me, DTs, 
        c=colors[2], label=r'0.05R lid')

ax[0].set_yscale('log')
ax[1].set_yscale('log')
ax[0].set_ylim(1e20,1e23)
ax[1].set_ylim(1e5,1.6e6)
ax[2].set_ylim(0,0.17)
ax[0].yaxis.set_major_formatter(tick.FuncFormatter(y_fmt))
ax[1].yaxis.set_major_formatter(tick.FuncFormatter(y_fmt))

plt.savefig(home+'/python/plots/RaNu.pdf',bbox_inches='tight')
plt.show()
