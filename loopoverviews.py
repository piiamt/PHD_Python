###################################################################
### this loops through all finished simulations and saves data ####
### into npy files and overview plots as well
###################################################################

import sys
import os, os.path
from scipy.io import readsav
from os.path import expanduser
import numpy as np
import itertools
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import FormatStrFormatter
home = '/lustre/astro/piiamt/locean_outgas/'
#home = '/projects/astro3/nobackup/piia/'
folders = np.array(os.listdir(home))
simulations = folders[np.where((folders!='tester')&(folders!='overview')&(folders!='python')&(folders!='starter'))]
### simulations.txt must already be ready for idllooper
#folders = np.array(os.listdir(home))
#np.savetxt('../idl/simlist.txt', simulations, fmt='%s', newline='\n')
simulations = np.loadtxt('../idl/simlist_testmedium.txt', dtype='str')
simulations = np.array(['1.0AU_0.1Mearth', '1.0AU_0.5Mearth',
                    '1.0AU_1.0Mearth', '1.0AU_1.5Mearth',
                    '1.0AU_2.0Mearth', '1.0AU_2.5Mearth',
                    '1.0AU_3.0Mearth', '1.0AU_4.0Mearth'])

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
C = np.array([2.55893462e-9, -8.26050794e-7, 1.07405538e-4,
              -7.22411940e-3, 2.87447212e-1, -3.01444278])

G = 6.674e-11
rhocore = 5663.2660	
rhomagma = 3400.0

def IWfn(Ps):
    PsG = Ps/1e9
    IW = C[0]*PsG**5 + C[1]*PsG**4 + C[2]*PsG**3+C[3]*PsG**2 + C[4]*PsG + C[5]
    return(IW)

MOBPs = np.array([])
DeltaIWs = np.array([])

### DONT DO BAD SIMULATIONS
for folder in simulations:
    timefile =home+folder+ '/time_series.idl'
    if not os.path.exists(timefile):
        simulations = simulations[np.where(simulations!=folder)]

for folder in simulations[:]:
    timefile =home+folder+ '/time_series.idl'
    print('about to read '+home+folder+'/time_series.idl')
    timef = readsav(timefile)
    timekey = list(timef.values())[0].dtype
    tser = list(timef.values())[0][0]
        
    t5 = np.where(tser.t>(0*Myr_to_s)) # time filter
### ZOOMED IN TIME FILTER:
#    t5 = np.where((tser.t>(3*Myr_to_s))&(tser.t<15*Myr_to_s))
    ts = tser.t[t5]/Myr_to_s

### getting the needed constants from the idl time files
    rcor = tser.rcor[-1]
    rpla = tser.rpla[-1]
    mcore = tser.mcore[-1]
    psur = tser.psur[-1]
    g = G*mcore/(rcor**2)
    h = rpla - rcor
    P = rhomagma * g * h
    MOBPs = np.append(MOBPs, P+psur)

    if not os.path.exists(home+folder+'/pysave/t.npy'):
        fig, ax = plt.subplots(3,4, figsize=(12,9), sharex=True, 
                               gridspec_kw={'wspace':0.2,'hspace':0.2})
# plt.gca().set_yscale('log')
        fig.suptitle(folder, fontsize=15)

        ax[0,0].set_ylabel(r'Mass ($M_{\oplus}$)')
        ax[0,0].set_yscale('log')
       # plt.setp(ax, ylim=[xmin,xmax])    #sets same limits for all subplots
        ax[0,0].plot(ts, tser.matm[t5]/mearth, c=colors[1])
        ax[0,0].set_title('Atmosphere M', fontsize=10)
        ax[0,1].plot(ts, tser.mpla[t5]/mearth, c=colors[3])
        ax[0,1].set_ylim(np.mean(tser.mpla[t5])/mearth-0.005, 
        		np.mean(tser.mpla[t5]/mearth+0.005))
        #ax[0,1].set_yscale('log')
        ax[0,1].yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
        ax[0,1].set_title('Planet M', fontsize=10)
        ax[0,2].plot(ts, tser.matm_o2[t5]/tser.matm[t5], c=colors[0])
        ax[0,2].set_ylim(-0.05,1.05)
        ax[0,2].set_title(r'O$_2$ / Matm', fontsize=10)
        ax[0,3].plot(ts, tser.matm_h2o[t5]/tser.matm[t5], c=colors[5])
        ax[0,3].set_ylabel(r'Atm mass fraction')
        ax[0,3].yaxis.set_label_position('right')
        ax[0,3].set_ylim(-0.05,1.05)
        ax[0,3].set_title(r'H$_2$O / Matm', fontsize=10)


        ax[1,0].set_ylabel(r'Mass ($M_{\oplus}$)')
        ax[1,0].set_yscale('log')
        ax[1,0].plot(ts, tser.menv[t5]/mearth, c=colors[4])
        ax[1,0].set_title('Envelope M', fontsize=10)
        ax[1,1].plot(ts, tser.tsur[t5], c=colors[5])
        ax[1,1].set_title('Surface T (K)', fontsize=10)
        ax[1,2].plot(ts, tser.matm_co[t5]/tser.matm[t5], c=colors[2])
        ax[1,2].set_ylim(-0.05,1.05)
        ax[1,2].set_title(r'CO / Matm', fontsize=10)
        ax[1,3].plot(ts, tser.matm_co2[t5]/tser.matm[t5], c=colors[1])
        ax[1,3].set_ylabel(r'Atm mass fraction')
        ax[1,3].yaxis.set_label_position('right')
        ax[1,3].set_ylim(-0.05,1.05)
        ax[1,3].set_title(r'CO$_2$ / Matm', fontsize=10)


        ax[2,0].set_ylabel(r'Mass ($M_{\oplus}$)')
        ax[2,0].set_yscale('log')
        ax[2,0].plot(ts, tser.msolid[t5]/mearth, c=colors[0])
        ax[2,0].set_title('Solids M', fontsize=10)
        ax[2,1].plot(ts, tser.psur[t5]/1e5, c=colors[1])
        ax[2,1].set_yscale('log')
        ax[2,1].set_title('Surface pressure (bar)', fontsize=10)
        ax[2,2].plot(ts, tser.matm_n2[t5]/tser.matm[t5], c=colors[3])
        ax[2,2].set_ylim(-0.05,1.05)
        ax[2,2].set_title(r'N$_2$ / Matm', fontsize=10)
        ax[2,3].plot(ts, tser.matm_h2[t5]/tser.matm[t5], c=colors[4])
        ax[2,3].set_ylim(-0.05,1.05)
        ax[2,3].set_title(r'H$_2$ / Matm', fontsize=10)
        ax[2,3].set_ylabel(r'Atm mass fraction')
        ax[2,3].yaxis.set_label_position('right')

        ax[2,0].set_xlabel('Time (Myr)')
        ax[2,1].set_xlabel('Time (Myr)')
        ax[2,2].set_xlabel('Time (Myr)')
        ax[2,3].set_xlabel('Time (Myr)')

        #home = '/projects/astro3/nobackup/piia/'
        #OGhome = '/lunarc/nobackup/projects/astro4/piia'
        plt.savefig(home+'python/plots/'+folder+'.png', bbox_inches='tight')
        print('should have just saved '+folder)
        #plt.show()

    if not os.path.exists(home+folder+'/pysave/t.npy'):
        ### saving all the needed files to npy for later FAST plotting
        dirname = home + folder + '/pysave'
        os.mkdir(dirname)

        np.save(dirname+'/t', tser.t)
        np.save(dirname+'/mpla', tser.mpla)
        np.save(dirname+'/menv', tser.menv)
        np.save(dirname+'/matm', tser.matm)
        np.save(dirname+'/msolid', tser.msolid)
        np.save(dirname+'/psur', tser.psur)
        np.save(dirname+'/tsur', tser.tsur)
        np.save(dirname+'/matm_o2', tser.matm_o2)
        np.save(dirname+'/matm_co', tser.matm_co)
        np.save(dirname+'/matm_co2', tser.matm_co2)
        np.save(dirname+'/matm_h2o', tser.matm_h2o)
        np.save(dirname+'/matm_n2', tser.matm_n2)
        np.save(dirname+'/matm_h2', tser.matm_h2)
        np.save(dirname+'/MOBP', tser.mobp)
        np.save(dirname+'/deltaIW', IWfn(tser.mobp))
        np.save(dirname+'/matm_sio', tser.matm_sio)
        np.save(dirname+'/dragphi_O', tser.dragphi_o)
        np.save(dirname+'/dragphi_N', tser.dragphi_n)
        np.save(dirname+'/dragphi_C', tser.dragphi_c)
    dirname = home+folder+'/pysave'
    np.save(dirname+'/dragphi_O', tser.dragphi_o)
    np.save(dirname+'/dragphi_N', tser.dragphi_n)
    np.save(dirname+'/dragphi_C', tser.dragphi_c)
    np.save(dirname+'/rpla', tser.rpla)
    np.save(dirname+'/rmagma_cor', tser.rmagma_cor)
    np.save(dirname+'/mpla_h2o', tser.mpla_h2o)
#	timefile = home + sim + '/time_series.idl'
#	timef = readsav(timefile)
#	timekey = list(timef.values())[0].dtype
#	tser = list(timef.values())[0][0]
#	### getting the needed constants from the idl time files
    rcor = tser.rcor[-1]
    rpla = tser.rpla[-1]
    mcore = tser.mcore[-1]
    psur = tser.psur[-1]
    g = G*mcore/(rcor**23)
    h = rpla - rcor
    P = rhomagma * g * h
    MOBPs = np.append(MOBPs, tser.mobp[-1])
#    MOBPs = np.append(MOBPs, P+psur)
    DeltaIWs = np.append(DeltaIWs, IWfn(tser.mobp[-1]))
#    np.save(dirname+'/deltaIW', IWfn(tser.mobp))

#np.save(home+'/overview/MOBPs', MOBPs)
#mobps = np.load(home+'/overview/MOBPs.npy')
#np.save(home+'/overview/IW', IWfn(mobps))

