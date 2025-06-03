import sys
"""
Created on Feb 27 2023
@author: piiamt
CREATES 3 PLOTS UNDER EACH OTHER WITH EACH PLOT SHOWING THE LOG T EVOLUTION OF
A DIFFERENT ATMOSPHERIC ELEMENT, COMPARING SIMULATIONS OF PLANET MASSES
0.1, 0.5, 2, 5 MEARTH IN DIFFERENT COLOURS.
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
#home = '/lustre/astro/piiamt/locean_outgas/long/'#dense_new/'#noatmloss_dense/'
home = '/lustre/astro/piiamt/delayed_XUV/'#dense_new/'#noatmloss_dense/'
home0 = '/lustre/astro/piiamt/'#log/testmedium/'
folders = np.array(['1.0AU_0.1Mearth','1.0AU_0.5Mearth',
                    '1.0AU_2.0Mearth','1.0AU_4.0Mearth'])
#                    '1.0AU_0.1Mearth_10Myr',
#                    '1.0AU_0.5Mearth_10Myr',
#                    '1.0AU_2.0Mearth_10Myr',
#                    '1.0AU_5.0Mearth_10Myr'])


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

colors = np.array(['k', '#009e73', '#d55e00', '#cc79a7',
                   'k', '#009e73', '#d55e00', '#cc79a7'])
    #'#0072b2', '#e69f00', '#56b4e9'])
RGBc = np.array([[0,0,0,0.2],[0,158/255,112/255,0.2],[213/255,94/255,0.0,0.2],[204/255,121/255,167/255,0.2],
                 [0,114/255,178/255,0.2],[230/255,159/255,0,0.2],[86/255,180/255,233/255,0.2]])
# aka           black   green     vermillion    pink        blue    orange      lskdfs
#################   CONSTANTS ############
mearth      = 5.9722e24     # kg
msun        = 1.98847e30    # kg
Myr_to_s    = 365.25*24*60*60*1000000
meatm       = 5.15e18       # kg
#timefile =home+folder+ '/time_series.idl'
#timef = readsav(timefile)
#timekey = list(timef.values())[0].dtype
#tser = list(timef.values())[0][0]
#varfile = home+folder+ '/variables.idl'
#var = list((readsav(varfile)).values())[0][0]

from matplotlib.colors import LogNorm
#home = '/lustre/astro/piiamt/'#noatmloss/'

Matm_CO2s = np.array([])#np.zeros(len(folders)-1)
Matm_O2s = np.array([])#np.zeros(len(folders)-1)
Matm_H2Os = np.array([])#np.zeros(len(folders)-1)
j = 0
for f in folders:
    #timefile =home+f+ '/time_series.idl'
    #timef = readsav(timefile)
    #tser = list(timef.values())[0][0]
    
    #Mpla = np.load(home + f + '/pysave'+'/mpla.npy')
    #A = Aplas[-1]
    #M = np.round(Mpla/mearth, 2)
    #dataAs = np.append(np.where(A==As), A)
    #dataMs = np.append(np.where(M==Ms), M)
    #Matm = np.load(home + f + '/pysave'+'/matm.npy')[-1]
    #Menv = np.load(home + f + '/pysave'+'/menv.npy')
    #Psur = np.load(home + f + '/pysave'+'/psur.npy')
    #Tsur = np.load(home + f + '/pysave'+'/tsur.npy')
    #ts = np.load(home + f + '/pysave'+'/t.npy')
    #Msolid = np.load(home + f + '/pysave'+'/msolid.npy')
    #MOBP = np.load(home + f + '/pysave/MOBP.npy')
    #IW = np.load(home + f + '/pysave/deltaIW.npy')
    #matm_CO2 = np.load(home + f + '/pysave'+'/matm_co2.npy')
    #matm_O2 = np.load(home + f + '/pysave'+'/matm_o2.npy')
    #matm_H2O = np.load(home + f + '/pysave'+'/matm_h2o.npy')
    #matm_CO = np.load(home + f + '/pysave'+'/matm_co.npy')
    #matm_N2 = np.load(home + f + '/pysave'+'/matm_n2.npy')
    #matm_H2 = np.load(home + f + '/pysave'+'/matm_h2.npy')
    #matm_SiO = np.load(home + f + '/pysave'+'/matm_sio.npy')
    #Matms = np.append(Matms, Matm)
    #Mplas = np.append(Mplas, Mpla[-1])
    #Menvs = np.append(Menvs, Menv[-1])
    #Psurs = np.append(Psurs, Psur[-1])
    #Tsurs = np.append(Tsurs, Tsur[-1])
    #MOBPs = np.append(MOBPs, MOBP[-1])
    #IWs = np.append(IWs, IW[-1])
    ### how many moles of each gas in atm:
    #mol_CO2 = matm_CO2 / muCO2
    #mol_CO = matm_CO / muCO
    #mol_N2 = matm_N2 / muN2
    #mol_H2O = matm_H2O / muH2O
    #mol_O2 = matm_O2 / muO2
    #mol_H2 = matm_H2 / muH2
    #mol_SiO = matm_SiO / muSiO
    #mol_atm = (matm_CO2/muCO2 + matm_CO/muCO + matm_N2/muN2 + matm_H2O/muH2O +
    #           matm_O2/muO2 + matm_H2/muH2 + matm_SiO/muSiO)
    #CO2s = np.append(CO2s, mol_CO2/mol_atm)
    #COs = np.append(COs, mol_CO/mol_atm)
    #N2s = np.append(N2s, mol_N2/mol_atm)
    #H2Os = np.append(H2Os, mol_H2O/mol_atm)
    #O2s = np.append(O2s, mol_O2/mol_atm)
    #H2s = np.append(H2s, mol_H2/mol_atm)
    #mCO2 = np.append(mCO2, matm_CO2/meatm)
    #mCO = np.append(mCO, matm_CO/meatm)
    #mO2 = np.append(mO2, matm_O2/meatm)
    #mH2O = np.append(mH2O, matm_H2O/meatm)
#    Matm_CO2s = np.append(Matm_CO2s, np.array([[matm_CO2]]))
#    Matm_O2s = np.append(Matm_O2s, np.array([[matm_O2]]))
#    Matm_H2Os = np.append(Matm_H2Os, np.array([[matm_H2O]]))
    j = j+1

############### 3PLOT  ################
mm=1/25.4
fig, ax = plt.subplots(3,2, sharex=True, figsize=(180*mm,120*mm), gridspec_kw=dict(hspace=0.0, wspace=0.0))

#ax[0].set_title(r'CO$_2$')
#ax[1].set_title(r'O$_2$')
#ax[2].set_title(r'H$_2$O')
lss = np.array(['-','--','-.',':'])
lss = np.array(['-', '-', '-', '-'])
labels = np.array([r'0.1 M$_{\oplus}$', r'0.5 M$_{\oplus}$', 
                   r'2.0 M$_{\oplus}$', r'4.0 M$_{\oplus}$'])

for i in np.arange(len(folders)):
    f = folders[i]
    timefile =home+f+ '/time_series.idl'
    timef = readsav(timefile)
    tser = list(timef.values())[0][0]
    ts = tser.t#np.load(home + f + '/pysave'+'/t.npy')
    if ts[-1]>(11*Myr_to_s):
        iplot = np.where(ts>(0*Myr_to_s))
        ts = ts[iplot]
        #envlosstime = ts[np.where((tser.menv<(max(tser.menv)*0.005))&(ts>(5.0*Myr_to_s)))][0]/Myr_to_s
        envlosstime = ts[np.where((tser.msolid>0)&(ts>(4.8*Myr_to_s)))][0]/Myr_to_s
    else:
        iplot = np.where(ts<=(10*Myr_to_s))
        ts = ts[iplot]
        try: 
            #envlosstime = ts[np.where((tser.menv<(max(tser.menv)*0.005))&(ts>(5.0*Myr_to_s)))][0]/Myr_to_s
            envlosstime = ts[np.where((tser.msolid>0)&(ts>(4.8*Myr_to_s)))][0]/Myr_to_s
        except:
            xxx=0
    #if f=='1.0AU_0.1Mearth': envlosstime=5.0*Myr_to_s
    #if f=='1.0AU_0.5Mearth': envlosstime=5.053*Myr_to_s
    Matm_CO2 = tser.matm_co2[iplot]#np.load(home + f + '/pysave'+'/matm_co2.npy')
    Matm_O2 = tser.matm_o2[iplot]#np.load(home + f + '/pysave'+'/matm_o2.npy')
    Matm_H2O = tser.matm_h2o[iplot]#np.load(home + f + '/pysave'+'/matm_h2o.npy')
    Matm_CO = tser.matm_co[iplot]#np.load(home + f + '/pysave/matm_co.npy')
    Matm_H2 = tser.matm_h2[iplot]#np.load(home + f + '/pysave/matm_h2.npy')
    Matm_N2 = tser.matm_n2[iplot]#np.load(home + f + '/pysave/matm_n2.npy')
    ax[0,0].vlines(envlosstime, 1e-9, 60, lw=0.6, ls='--', color=colors[i])
    ax[0,0].plot(ts/Myr_to_s, Matm_CO2/mearth, ls=lss[i], 
#    ax[0].plot(ts/Myr_to_s, tser.mcor_co2/meatm, ls=lss[i],
               c=colors[i], label=labels[i])
    ax[1,0].vlines(envlosstime, 1e-9, 60, lw=0.6, ls='--', color=colors[i])
    ax[1,0].plot(ts/Myr_to_s, Matm_O2/mearth, ls=lss[i], 
#    ax[1].plot(ts/Myr_to_s, tser.mcor_h2o/meatm, ls=lss[i],
               c=colors[i], label=labels[i])
    ax[2,0].vlines(envlosstime, 1e-9, 60, lw=0.6, ls='--', color=colors[i])
    ax[2,0].plot(ts/Myr_to_s, (tser.matm_h2o+tser.mpla_h2o)/mearth, ls=lss[i], 
#    ax[2].plot(ts/Myr_to_s, tser.msil_h2o/meatm, ls=lss[i],
               c=colors[i], label=labels[i])
    ax[0,1].vlines(envlosstime, 1e-9, 60, lw=0.6, ls='--', color=colors[i])
    ax[0,1].plot(ts/Myr_to_s, Matm_CO/mearth, ls=lss[i], 
               c=colors[i], label=labels[i])
    ax[1,1].vlines(envlosstime, 1e-9, 60, lw=0.6, ls='--', color=colors[i])
    ax[1,1].plot(ts/Myr_to_s, Matm_H2/mearth, ls=lss[i], 
               c=colors[i], label=labels[i])
    ax[2,1].vlines(envlosstime, 1e-9, 60, lw=0.6, ls='--', color=colors[i])
    ax[2,1].plot(ts/Myr_to_s, Matm_N2/mearth, ls=lss[i], 
               c=colors[i], label=labels[i])
#ax[1].plot(tser.t/Myr_to_s, tser.rpla/1000, ls=':', c='w', label='Surface')
#ax[1].plot(tser.t/Myr_to_s, tser.rmagma_cor/1000, ls='-.', c='w', label='Magma core')
#ax[1].plot(tser.t/Myr_to_s, tser.rmagma_sur/1000, ls='', c='w', label='Magma surface')

ax[2,1].legend(loc=1, frameon=False)
#ax[2,1].legend(loc='lower center', bbox_to_anchor=(0.6,0.33), frameon=False)
#ax[1].legend()
#ax[2].legend()
ax[0,0].text(1.4,1e-6,r'CO$_2$', fontsize=10)
ax[1,0].text(1.4,1e-6,r'O$_2$', fontsize=10)
ax[2,0].text(1.4,1e-6,r'H$_2$O', fontsize=10)
ax[0,1].text(1.4,1e-6, r'CO', fontsize=10)
ax[1,1].text(1.4,1e-6, r'H$_2$', fontsize=10)
ax[2,1].text(1.4,1e-6, r'N$_2$', fontsize=10)

for i in np.arange(len(folders)):
    f = folders[i]
    timefile =home0+f+ '/time_series.idl'
    timef = readsav(timefile)
    tser = list(timef.values())[0][0]
    ts = tser.t#np.load(home + f + '/pysave'+'/t.npy')
    if ts[-1]>(11*Myr_to_s):
        iplot = np.where(ts>(0*Myr_to_s))
        ts = ts[iplot]
        envlosstime = ts[np.where((tser.msolid>0)&(ts>(4.8*Myr_to_s)))][0]/Myr_to_s
    else:
        iplot = np.where(ts<=(10*Myr_to_s))
        ts = ts[iplot]
        try: 
            envlosstime = ts[np.where((tser.msolid>0)&(ts>(4.8*Myr_to_s)))][0]/Myr_to_s
        except:
            xxx=0
    Matm_CO2 = tser.matm_co2[iplot]#np.load(home + f + '/pysave'+'/matm_co2.npy')
    Matm_O2 = tser.matm_o2[iplot]#np.load(home + f + '/pysave'+'/matm_o2.npy')
    Matm_H2O = tser.matm_h2o[iplot]#np.load(home + f + '/pysave'+'/matm_h2o.npy')
    Matm_CO = tser.matm_co[iplot]#np.load(home + f + '/pysave/matm_co.npy')
    Matm_H2 = tser.matm_h2[iplot]#np.load(home + f + '/pysave/matm_h2.npy')
    Matm_N2 = tser.matm_n2[iplot]#np.load(home + f + '/pysave/matm_n2.npy')
    ax[0,0].plot(ts/Myr_to_s, Matm_CO2/mearth, ls=lss[i], 
               c=colors[i], lw=0.6)
    ax[1,0].plot(ts/Myr_to_s, Matm_O2/mearth, ls=lss[i], 
               c=colors[i], lw=0.6)
    ax[2,0].plot(ts/Myr_to_s, (tser.matm_h2o+tser.mpla_h2o)/mearth, ls=lss[i], 
               c=colors[i], lw=0.6)
    ax[0,1].plot(ts/Myr_to_s, Matm_CO/mearth, ls=lss[i], 
               c=colors[i], lw=0.6)
    ax[1,1].plot(ts/Myr_to_s, Matm_H2/mearth, ls=lss[i], 
               c=colors[i], lw=0.6)
    ax[2,1].plot(ts/Myr_to_s, Matm_N2/mearth, ls=lss[i], 
               c=colors[i], lw=0.6)

#ax[0,0].text(12,1e-3,'Envelope loss', fontsize=8)
#ax[0,0].arrow(11.7,2e-3,-2.5,0, color='k', width=5e-5, 
#        head_width=1.2e-3, head_length=0.6, head_starts_at_zero=False)
#ax[0,0].arrow(64,2e-3,20,0, color='k', width=5e-5, 
#        head_width=1.2e-3, head_length=10, head_starts_at_zero=False)

#ax[0].text(1.3,100, r'Core CO$_2$', fontsize=14)
#ax[1].text(1.3,1000, r'Core H$_2$O', fontsize=14)
#ax[2].text(1.3,100, r'Silicate H$_2$O', fontsize=14)

###ax[0,0].set_xlim(1e0,80)
###ax[0,1].set_xlim(1e0,80)
#ax[0,1].set_ylim(1e-9,8e-4)#(1e-4,1000)
#ax[1,1].set_ylim(1e-9,7e-3)#(1e-4,1000)
#ax[2,1].set_ylim(1e-10,2e-6)#(1e-4,5)
ax[0,0].set_ylim(1e-9,5e-2)#old max 5e-2
ax[1,0].set_ylim(1e-9,5e-2)#(1e-6,200000)
ax[2,0].set_ylim(1e-9,5e-2)#(1e-2,9000)
ax[0,1].set_ylim(1e-9,5e-2)#old max 5e-2
ax[1,1].set_ylim(1e-9,5e-2)#(1e-6,200000)
ax[2,1].set_ylim(1e-9,5e-2)#(1e-2,9000)

ax[0,0].set_yscale('log')
ax[1,0].set_yscale('log')
ax[2,0].set_yscale('log')
ax[0,1].set_yscale('log')
ax[1,1].set_yscale('log')
ax[2,1].set_yscale('log')
import matplotlib.ticker as ticker
#ax[0,0].yaxis.set_major_locator(ticker.FixedLocator([1e-9,1e-7,1e-5,1e-3]))
#ax[0,0].yaxis.set_minor_locator(ticker.FixedLocator([1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2]))

ax[0,0].yaxis.set_major_locator(ticker.LogLocator(base=10, numticks=8))
ax[0,1].yaxis.set_major_locator(ticker.LogLocator(base=10, numticks=8))
ax[1,0].yaxis.set_major_locator(ticker.LogLocator(base=10, numticks=8))
ax[1,1].yaxis.set_major_locator(ticker.LogLocator(base=10, numticks=8))
ax[2,0].yaxis.set_major_locator(ticker.LogLocator(base=10, numticks=8))
ax[2,1].yaxis.set_major_locator(ticker.LogLocator(base=10, numticks=8))

ax[0,0].tick_params(labelsize=6)
ax[0,1].tick_params(labelsize=6)
ax[1,0].tick_params(labelsize=6)
ax[1,1].tick_params(labelsize=6)
ax[2,0].tick_params(labelsize=6)
ax[2,1].tick_params(labelsize=6)
ax[0,0].tick_params(right=True)
ax[1,0].tick_params(right=True)
ax[2,0].tick_params(right=True)
#ax[0,0].set_yticks(np.array([1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2]))
#ax[0,0].set_yticklabels(np.array(['1e-9','1e-8','1e-7','1e-6','1e-5','1e-4','1e-3','1e-2']))
#ax[0,1].set_yticks([1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2])
ax[0,1].yaxis.set_label_position('right')
ax[0,1].yaxis.tick_right()
ax[0,1].tick_params(left=True)
#ax00 = ax[0,0].twinx()
#ax00.set_ylim(ax[0,0].get_ylim())
#ax00.set_yticks(ax[0,0].get_yticks())
#ax[1,0].set_yticks([1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2])
#ax[1,1].set_yticks([1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2])
ax[1,1].yaxis.set_label_position('right')
ax[1,1].yaxis.tick_right()
ax[1,1].tick_params(left=True)
#ax[2,0].set_yticks([1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2])
#ax[2,1].set_yticks([1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2])
ax[2,1].yaxis.set_label_position('right')
ax[2,1].yaxis.tick_right()
ax[2,1].tick_params(left=True)
#ax2 = ax[0,1].twinx()
#ax2.set_ylabels([1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2])
ax[0,0].set_xlim(1e0,300)#1e0,50)
ax[0,1].set_xlim(1e0,300)
ax[1,0].set_xlim(1e0,300)
ax[1,1].set_xlim(1e0,300)
ax[2,0].set_xlim(1e0,300)
ax[2,1].set_xlim(1e0,300)



ax[0,0].set_xscale('log')
ax[1,0].set_xscale('log')
ax[2,0].set_xscale('log')
ax[0,1].set_xscale('log')
ax[1,1].set_xscale('log')
ax[2,1].set_xscale('log')
ax[2,0].set_xlabel('Time (Myr)')
ax[2,1].set_xlabel('Time (Myr)')
ax[1,0].set_ylabel(r'Mass (M$_{\oplus}$)')

#plt.savefig(home+'python/plots/contour_magma_long.png', bbox_inches='tight')
plt.savefig(home+'/python/plots/atmcomparisons_delayed_XUV.pdf', bbox_inches='tight')
#plt.savefig('1.1Mearth_300Myr_contourRHO_T.pdf', bbox_inches='tight')

plt.show()



