"""
SHOWS THE TIME OF ENVELOPE LOSS ON Y AXIS WITH THE PLANET MASS AND/OR 
ORBITAL RADIUS ON THE X AXIS.
@author: piiamt
"""
from scipy.io import readsav
import numpy as np
import itertools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os, os.path
from os.path import expanduser
#OGhome = expanduser('~')
#home = '/projects/astro3/nobackup/piia/'
#from simlist import folders
home = '/lunarc/nobackup/projects/astro4/piia/'
#simulations = np.array(os.listdir(home))
folders = np.loadtxt('../idl/simlist.txt', dtype='str')
######################


plt.rcParams.update({'font.size'           : 10, 
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
amu         = 1.660539040e-27
muO2        = 2*2.34*amu
muCO2       = 12.0107*amu + muO2
muCO        = (12.0107+2.34)*amu
muH2        = 2*1.00784*amu
muH2O       = muH2 + 2.34*amu
muN2        = 2*14.0067*amu
muSiO       = (28.0855+2.34)*amu
mearth      = 5.9722e24     # kg
msun        = 1.98847e30    # kg
Myr_to_s    = 365.25*24*60*60*1000000


#%% 
### INITIATE EMPTY ARRAYS
Mplas = np.array([])
Matms = np.array([])
Aplas = np.array([]) #np.array([0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0])
Menvs = np.array([])
Psurs = np.array([])
Tsurs = np.array([])
MOBPs = np.array([])
IWs = np.array([])
CO2s = np.array([])
O2s = np.array([])
H2Os = np.array([])
COs = np.array([])
N2s = np.array([])
H2s = np.array([])
dragphi_Os = np.array([])
dragphi_Cs = np.array([])
dragphi_Ns = np.array([])
losstimes = np.array([])

### DONT DO BAD SIMULATIONS
for folder in folders:
    timefile =home+folder+ '/time_series.idl'
    if not os.path.exists(timefile):
        folders = folders[np.where(folders!=folder)]


### LOOPING THROUGH ALL THE FILES
for f in folders[:]:
    if f[3]!='A':	# if Apla has 2 decimal places like 1.25AU
        Aplas = np.append(Aplas, float(f[:4]))
    elif f[3]=='A':	# aka if Apla has 1 decimal place like 1.0AU
        Aplas = np.append(Aplas, float(f[:3]))
    Mpla = np.load(home + f + '/pysave'+'/mpla.npy')
    Matm = np.load(home + f + '/pysave'+'/matm.npy')[-1]
    Menv = np.load(home + f + '/pysave'+'/menv.npy')
    Psur = np.load(home + f + '/pysave'+'/psur.npy')
    Tsur = np.load(home + f + '/pysave'+'/tsur.npy')
    ts = np.load(home + f + '/pysave'+'/t.npy')
    Msolid = np.load(home + f + '/pysave'+'/msolid.npy')
    MOBP = np.load(home + f + '/pysave/MOBP.npy')
    IW = np.load(home + f + '/pysave/deltaIW.npy')
    matm_CO2 = np.load(home + f + '/pysave'+'/matm_co2.npy')[-1]
    matm_O2 = np.load(home + f + '/pysave'+'/matm_o2.npy')[-1]
    matm_H2O = np.load(home + f + '/pysave'+'/matm_h2o.npy')[-1]
    matm_CO = np.load(home + f + '/pysave'+'/matm_co.npy')[-1]
    matm_N2 = np.load(home + f + '/pysave'+'/matm_n2.npy')[-1]
    matm_H2 = np.load(home + f + '/pysave'+'/matm_h2.npy')[-1]
    matm_SiO = np.load(home + f + '/pysave'+'/matm_sio.npy')[-1]
    dragphi_O = np.load(home + f + '/pysave'+'/dragphi_O.npy')
    dragphi_N = np.load(home + f + '/pysave'+'/dragphi_N.npy')
    dragphi_C = np.load(home + f + '/pysave'+'/dragphi_C.npy')
    Matms = np.append(Matms, Matm)
    Mplas = np.append(Mplas, Mpla[-1])
    Menvs = np.append(Menvs, Menv[-1])
    Psurs = np.append(Psurs, Psur[-1])
    Tsurs = np.append(Tsurs, Tsur[-1])
    MOBPs = np.append(MOBPs, MOBP[-1])
    IWs = np.append(IWs, IW[-1])
    ### how many moles of each gas in atm:
    mol_CO2 = matm_CO2 / muCO2
    mol_CO = matm_CO / muCO
    mol_N2 = matm_N2 / muN2
    mol_H2O = matm_H2O / muH2O
    mol_O2 = matm_O2 / muO2
    mol_H2 = matm_H2 / muH2
    mol_SiO = matm_SiO / muSiO
    mol_atm = (matm_CO2/muCO2 + matm_CO/muCO + matm_N2/muN2 + matm_H2O/muH2O +
               matm_O2/muO2 + matm_H2/muH2 + matm_SiO/muSiO)
    CO2s = np.append(CO2s, mol_CO2/mol_atm)
    COs = np.append(COs, mol_CO/mol_atm)
    N2s = np.append(N2s, mol_N2/mol_atm)
    H2Os = np.append(H2Os, mol_H2O/mol_atm)
    O2s = np.append(O2s, mol_O2/mol_atm)
    H2s = np.append(H2s, mol_H2/mol_atm)
    dragphi_Os = np.append(dragphi_Os, np.max(dragphi_O))
    dragphi_Ns = np.append(dragphi_Ns, np.max(dragphi_N))
    dragphi_Cs = np.append(dragphi_Cs, np.max(dragphi_C))
    if Menv[-1]==0:
        losstime = ts[np.where((Msolid>0)&(ts>(5*Myr_to_s)))][0]
        losstimes = np.append(losstimes, losstime)
lostenv = np.where(Menvs==0)
envstay = np.where(Menvs>0)
IW = IWs[lostenv]

smalls = np.where(dragphi_Os<1)
bigs = np.where(dragphi_Os>=1)

fig, ax = plt.subplots(1,1, figsize=(4,4))
ax.set_xlabel('Orbital radius (AU)')
ax.set_ylabel(r'Planet mass (M$_{\oplus}$)')
at = ax.scatter(Aplas, Mplas/mearth, s=2**6.8,
                #c=losstimes/Myr_to_s, cmap='inferno', edgecolor=None)
                c=dragphi_Os, cmap='turbo', edgecolor=None, norm=matplotlib.colors.LogNorm())
                #c=dragphi_Cs, cmap='turbo', edgecolor=None, norm=matplotlib.colors.LogNorm())
ax.scatter(Aplas[smalls], Mplas[smalls]/mearth, s=2**6.8,
           facecolors='none', edgecolors='k')
p = ax.get_position().get_points().flatten() #[[x0,y0],[x1,y1]]
at_cbar = fig.add_axes([p[0], p[3], p[2]-p[0], 0.02])
fig.colorbar(at, cax=at_cbar, orientation='horizontal', 
             #label='Time of envelope loss (Myr)')
             label='dragphi')
at_cbar.xaxis.set_label_position('top')
at_cbar.xaxis.set_ticks_position('top')


plt.savefig(home+'/python/plots/dragphi_Os.png', bbox_inches='tight')
plt.show()
