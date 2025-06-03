"""
CREATES HYDROSPHERE PLOT OF ALL SIMULATIONS AT 1 AU
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
OGhome = expanduser('~')
home0 = '/lustre/astro/piiamt/'
home = '/lustre/astro/piiamt/locean_outgas/'
from simlist import folders
folders = np.loadtxt('../idl/1AUsimlist.txt', dtype='str')
folders = np.array(['1.0AU_0.1Mearth', '1.0AU_0.5Mearth',
                    '1.0AU_1.0Mearth', '1.0AU_1.5Mearth',
                    '1.0AU_2.0Mearth', '1.0AU_2.5Mearth',
                    '1.0AU_3.0Mearth', '1.0AU_4.0Mearth'])
######################


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
# aka           black   green     vermillion    pink        blue    orange      lskdfs
#################   CONSTANTS ############
mearth      = 5.9722e24     # kg
msun        = 1.98847e30    # kg
Myr_to_s    = 365.25*24*60*60*1000000
meocean     = 1.4e21        # kg

### INITIATE EMPTY ARRAYS
Mplas = np.array([])
Aplas = np.array([]) #np.array([0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0])
Menvs = np.array([])
Psurs = np.array([])
Tsurs = np.array([])
losstimes = np.array([])
fCO2 = np.array([])
Matm_h2os = np.array([])
Matm_h2os0 = np.array([])
Mpla_h2os = np.array([])
Mpla_h2os0 = np.array([])
Msil_h2os = np.array([])
Msil_h2os0 = np.array([])
Msil_co2s = np.array([])
Msil_n2s = np.array([])
Mcor_co2s = np.array([])
Mcor_n2s = np.array([])
Mcor_h2os = np.array([])

### LOOPING THROUGH ALL THE FILES
for f in folders[:]:
    if f[3]!='A':	# if Apla has 2 decimal places like 1.25AU
        Apla = float(f[:4])        
    elif f[3]=='A':	# aka if Apla has 1 decimal place like 1.0AU
        Apla = float(f[:3])        
    timefile =home+ f + '/time_series.idl'
    timefile0 =home0+ f + '/time_series.idl'
    timef = readsav(timefile)
    timef0 = readsav(timefile0)
    timekey = list(timef.values())[0].dtype
    timekey0 = list(timef0.values())[0].dtype
    tser = list(timef.values())[0][0]
    tser0 = list(timef0.values())[0][0]
    Mpla = np.load(home + f + '/pysave'+'/mpla.npy')
    Menv = np.load(home + f + '/pysave'+'/menv.npy')
    ts = np.load(home + f + '/pysave'+'/t.npy')
    Msolid = np.load(home + f + '/pysave'+'/msolid.npy')
    Psur = np.load(home + f + '/pysave'+'/psur.npy')[-1]
    Tsur = np.load(home + f + '/pysave'+'/tsur.npy')[-1]
    Matm = np.load(home + f + '/pysave'+'/matm.npy')[-1]
    Matm_h2o = np.load(home + f + '/pysave/matm_h2o.npy')[-1]
    Mpla_h2o = np.load(home + f + '/pysave/mpla_h2o.npy')[-1]
    Matm_CO2 = np.load(home + f + '/pysave'+'/matm_co2.npy')[-1]
    Menvs = np.append(Menvs, Menv[-1])
    Msil_h2os = np.append(Msil_h2os, tser.msil_h2o[-1])
    Msil_h2os0 = np.append(Msil_h2os0, tser0.msil_h2o[-1])
    Msil_co2s = np.append(Msil_co2s, tser.msil_co2[-1])
    Msil_n2s = np.append(Msil_n2s, tser.msil_n2[-1])
    Mcor_h2os = np.append(Mcor_h2os, tser.mcor_h2o[-1])
    Mcor_n2s = np.append(Mcor_n2s, tser.mcor_n2[-1])
    Mcor_co2s = np.append(Mcor_co2s, tser.mcor_co2[-1])
    if Menv[-1]==0:
        losstime = ts[np.where((Msolid>0)&(ts>(5*Myr_to_s)))][0]
        losstimes = np.append(losstimes, losstime)
        Psurs = np.append(Psurs, Psur)
        Tsurs = np.append(Tsurs, Tsur)
        Aplas = np.append(Aplas, Apla)
        Mplas = np.append(Mplas, Mpla[-1])
        fCO2 = np.append(fCO2, Matm_CO2/Matm)
        Mpla_h2os = np.append(Mpla_h2os, Mpla_h2o)
        Matm_h2os = np.append(Matm_h2os, Matm_h2o)
        Mpla_h2os0 = np.append(Mpla_h2os0, tser0.mpla_h2o[-1])
        Matm_h2os0 = np.append(Matm_h2os0, tser0.matm_h2o[-1])
lostenv = np.where(Menvs==0)
envstay = np.where(Menvs>0)
liqgas = np.polyfit(([101325,22064000]), ([273.15,647.096]), 1)
solgas = np.polyfit(([0,0.006*101325]), ([0,275.15]), 1)

liq = np.where((Tsurs>273.15)&(Tsurs>(liqgas[0]*Psurs+liqgas[1]))&(Psurs<=(101325*218)))

#liq = np.where((Tsurs<647.096)&(Psurs<22064000)&(Tsurs>273.15)&(Psurs>101325))
crit = np.where((Tsurs>=647.096)|(Psurs>=22064000))
sol = np.where((Tsurs<=273.15)&(Tsurs>(solgas[0]*Psurs+solgas[1])))
gas = np.where((Tsurs<=(solgas[0]*Psurs+solgas[1]))&(Tsurs<=(liqgas[0]*Psurs+liqgas[1])))
#print('Tsurs',Tsurs)
#print('Mplas', Mplas/mearth)
#print('Psurs', Psurs)
co2 = np.array([[200,216.58,304.18], [2.0,5.185,73.8]])
co2_crit = np.array([[304.18, max(Tsurs)+10], [73.8, 73.8]])
#print('env lost w masses', Menvs[lostenv])
#print('anv stays w masses',Menvs[envstay])
#print('all Menvs are ',Menvs)
#print('Mplas',Mplas,Mplas[lostenv], 'and aplas',Aplas,Aplas[lostenv])
#%%
mm=1/25.4
### OCEAN MASS / PLANET MASS RATIO
fig, ax = plt.subplots(1,1, figsize=(88*mm,66*mm))
ax.set_xlabel('Planet mass')
ax.set_ylabel(r'Water mass (Earth ocean masses)')#('Water mass / planet mass')
#ax.scatter(Mplas/mearth, Msil_h2os/mearth, c=colors[2], label=r'H$_2$O')
#ax.scatter(Mplas/mearth, Msil_co2s/mearth, c=colors[4], label=r'CO$_2$')

ax.plot(Mplas/mearth, Msil_h2os/meocean, 
        c=colors[0], alpha=0.5, label=r'H$_2$O in mantle')

ax.plot(Mplas/mearth, (Mpla_h2os+Matm_h2os)/meocean, c=colors[4], label=r'Surface H$_2$O')
ax.plot(([0.1,4]), ([meocean/meocean, meocean/meocean]), 
        c=colors[0], alpha=0.4, ls='dashed', label='Earth ocean')


ax.plot(Mplas/mearth, Msil_h2os0/meocean, 
        c=colors[0], alpha=0.5, lw=0.7)

ax.plot(Mplas/mearth, (Mpla_h2os0+Matm_h2os0)/meocean, 
        c=colors[4], lw=0.7)

### MANUALLY SETTING ALL TO LIQ BC I CHECKED AND IT SHOULD ALL BE LIQ
liq = np.where(Tsurs>0)
#ax.scatter(Mplas[crit]/mearth, Mpla_h2os[crit]+Matm_h2os[crit]/meocean, 
#        facecolors=colors[0], edgecolor=colors[0], label='Supercritical')
#ax.scatter(Mplas[liq]/mearth, (Mpla_h2os[liq]+Matm_h2os[liq])/meocean,
#        s=12, facecolors='none', edgecolor=colors[6], label='Liquid')
#ax.scatter(Mplas[sol]/mearth, (Mpla_h2os[sol]+Matm_h2os[sol])/meocean,
#        s=12, facecolors='w', edgecolor='k', label='Solid')

#print('sol', Mplas[sol], 'liq', Mplas[liq], 'crit', Mplas[crit])

#ax.scatter(Mplas[gas]/mearth, (Mpla_h2os[gas]+Matm_h2os[gas])/Mplas[gas],
#        facecolors='w', edgecolor=colors[2], label='Gas')
#ax.set_ylim(0,6e-2)
ax.legend(frameon=False, bbox_to_anchor=(0.4,0.06,0.3,0.4))

#planets = ax.scatter(Mplas, Aplas, c=fCO2, cmap='inferno', s=2**7)
#p = ax.get_position().get_points().flatten()
#cbar=fig.add_axes([p[0],p[3]+0.005, p[2]-p[0],0.05])
#plt.colorbar(planets, cax=cbar, orientation='horizontal')
#cbar.xaxis.set_ticks_position('top')

#ax.set_xlabel('Temperature (K)')
#ax.set_ylabel(r'Pressure (bar)')
#ax.set_yscale('log')
#ax.plot(co2[0], co2[1], c=colors[0])
#ax.plot(co2_crit[0], co2_crit[1], c=colors[0], ls='--')
#planets = ax.scatter(Tsurs, Psurs*fCO2/1e5, s=2**8, 
#		c=Aplas, cmap='inferno')
#p = ax.get_position().get_points().flatten() # [[x0,y0], [x1,y1]]
#cbar = fig.add_axes([p[0], p[3]+0.008, p[2]-p[0], 0.05])
#plt.colorbar(planets, cax=cbar, orientation='horizontal',
#		label='Orbital radius (AU)')
#cbar.xaxis.set_label_position('top')
#cbar.xaxis.set_ticks_position('top')
#plt.savefig(home+'/python/plots/silreservoirs.pdf', bbox_inches='tight')
plt.savefig(home+'/python/plots/water_locean_outgas.pdf', bbox_inches='tight')
plt.show()
