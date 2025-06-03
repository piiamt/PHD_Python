'''
Uses the Mors star evolution tracks to retrieve and save the three XUV 
luminosities used in Johnstone+2021, also plotting them
'''
from scipy.io import readsav
import os, os.path
home = '/lustre/astro/piiamt/dense/'#dense_new/'#noatmloss_dense/'
import sys
sys.path.insert(0, '/groups/astro/piiamt/Mors')
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import Mors as mors
stardir = '/groups/astro/piiamt/packages/fs255_grid/'
#myStar = mors.StarEvo(starEvoDir='/groups/astro/piiamt/packages/fs255_grid/')

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
# aka           black   green     vermillion    pink        blue    orange      lskdfs
Myr_to_s    = 365.25*24*60*60*1000000


star03 = mors.Star(Mstar=1.0, Omega=1.5)
star55 = mors.Star(Mstar=1.0, Omega=5.5)
star32 = mors.Star(Mstar=1.0, Omega=32.0)
#print('Units of Lx is ',star.Units['Lx'],' and units of Age is ',star.Units['Age'])

### SAVING SLOW, MEDIUM AND FAST ROTATOR TRACKS
names = np.array(['slow','medium','fast'])
stars = np.array([star03, star55, star32])
for i in np.array([0,1,2]):
    filename = names[i]
    file = open('../adap/f90/XUV/'+filename+'.txt', 'w')
    for line in (stars[i].Tracks['Lx']+stars[i].Tracks['Leuv']):
        file.write(str(line) + '\n')
    file.close()
    file = open('../adap/f90/XUV/'+filename+'time.txt', 'w')
    for line in stars[i].Tracks['Age']:
        file.write(str(line) + '\n')
    file.close()

### RETRIEVING XUV OF ADAP SIMULATIONS
f = '1.0AU_2.0Mearth'
timefile =home+f+ '/time_series.idl'
timef = readsav(timefile)
tser = list(timef.values())[0][0]
LXUV1 = 5e23 # J/s aka W
adapXUV = LXUV1*(star03.Tracks['Age']/5)**(-0.75)  # W

#print(adapXUV[100], star03.Tracks['Lx'][10])
def LXtoLEUV(LX):
    return(10**4.8*LX**(0.86))

def LEUVtoLX(LEUV):
    return(10**(-4.8/0.86)*LEUV**(1/0.86))

def FXtoFXUV(FX):
    return()

def FXUVtoFX(FXUV):
    return()
#print(star03.PrintAvailableTracks())
#print(star03.Units['Feuv'], star03.Units['Fx'])

FXUV03 = (star03.Tracks['Leuv']+star03.Tracks['Lx'])*1e-7 #in W now
FXUV55 = (star55.Tracks['Leuv']+star55.Tracks['Lx'])*1e-7
FXUV32 = (star32.Tracks['Leuv']+star32.Tracks['Lx'])*1e-7

### PLOTTING XUV TRACKS
mm=1/25.4
fig,ax = plt.subplots(1,1,figsize=(88*mm,88*mm))
ax.plot(star03.Tracks['Age'], adapXUV, c=colors[0], label='Johansen et al. 2023')
ax.plot(star03.Tracks['Age'], FXUV03, c=colors[1], label=r'$\Omega_{\star}$=1.5 (slow)')
ax.plot(star55.Tracks['Age'], FXUV55, c=colors[2], label=r'$\Omega_{\star}$=5.5 (medium)')
ax.plot(star32.Tracks['Age'], FXUV32, c=colors[3], label=r'$\Omega_{\star}$=32 (fast)')
#secax = ax.secondary_yaxis('right', functions=(LXtoLEUV, LEUVtoLX))
#secax = ax.secondary_yaxis('right', functions=())
#secax.set_yscale('log')
#secax.set_ylabel('F_EUV (erg$\cdot$s$^{-1}\cdot$cm$^{-2}$)')
#secax.set_ylabel('L_EUV (erg/s)')
ax.set_xscale('log')
ax.set_yscale('log')
#ax.set_xlim(1,1000)
ax.set_xlabel('Age (Myr)')
ax.set_ylabel(r'L$_{\mathrm{XUV}}$ (W)')
#ax.set_ylabel(r'F$_{\mathrm{XUV}}$ (W/m$_{2}$)')
#ax.set_ylabel(r'Fxuv (erg$\cdot$s$^{-1}\cdot$cm$^{-2}$)')
#ax.set_ylabel('Lx (erg/s)')
ax.legend()
#plt.savefig('plots/LXUV_Johnstone.pdf', bbox_inches='tight')
plt.show()
