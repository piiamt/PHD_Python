"""
MAKES 6 DIFFERENT MEGACOMPAREPLOTS TO SHOW ATM COMPOSITION OF ALL SIMS
TO SHOW MAX ENV MASS AND ENV LOSS TIME IN ROW 0 AND 
IN ROW 1 SHOW THE MOLECULAR COMPOSITION OF THE 3 ROTATORS.
@author: piiamt
"""
from scipy.io import readsav
import numpy as np
import itertools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import ticker
import os, os.path
from os.path import expanduser
OGhome = '/groups/astro/piiamt/'
homes = np.array(['/lustre/astro/piiamt/log/slow/',
                  '/lustre/astro/piiamt/log/medium/',
                  '/lustre/astro/piiamt/log/fast/'])
folders = np.loadtxt('/groups/astro/piiamt/idl/simlist.txt', dtype='str')
fastslowmedium = np.array(['slow','medium','fast'])
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

colours = np.array(['k', '#009e73', '#d55e00', '#cc79a7', '#0072b2', '#e69f00', '#56b4e9'])
RGBc = np.array([[0,0,0,0.2],[0,158/255,112/255,0.2],[213/255,94/255,0.0,0.2],[204/255,121/255,167/255,0.2],
                 [0,114/255,178/255,0.2],[230/255,159/255,0,0.2],[86/255,180/255,233/255,0.2]])
# aka           black   green     vermillion    pink        blue    orange      lskdfs
#################   CONSTANTS ############
amu         = 1.660539040e-27
muO2        = 2*15.999*amu #CHECK
muCO2       = 12.0107*amu + muO2
muCO        = (12.0107+15.999)*amu
muH2        = 2*1.00784*amu
muH2O       = 2*muH2 + 15.999*amu
muN2        = 2*14.0067*amu
muSiO       = (28.0855+15.999)*amu
mearth      = 5.9722e24     # kg
msun        = 1.98847e30    # kg
Myr_to_s    = 365.25*24*60*60*1000000

nA          = 6.022140857e23
meatm       = 5.15e18       # kg

As = np.array([0.7, 0.9, 1.0, 1.1, 1.3, 1.4, 1.5, 1.6, 1.7])
Ms = np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
               1.1,1.3,1.5,2.0,2.5,3.0,3.2,3.5,3.7,4.0,4.2,4.5,4.7,5.0])
#datamap = np.zeros(len(As), len(Ms)) #here I will store the sorted data indexes

#%% 
### INITIATE EMPTY ARRAYS
Mplas1 = np.array([])#([[],[],[]])
Mplas2 = np.array([])
Mplas3 = np.array([])
Matms = np.array([])
Aplas1 = np.array([])#([[],[],[]]) #np.array([0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0])
Aplas2 = np.array([])
Aplas3 = np.array([])
Menvs = np.array([])
Psurs = np.array([])
Tsurs = np.array([])
MOBPs = np.array([])
IWs = np.array([])
mCO2 = np.array([])
CO2s = np.array([])
mO2 = np.array([])
O2s = np.array([])
mH2O = np.array([])
H2Os = np.array([])
mCO = np.array([])
COs = np.array([])
mN2 = np.array([])
N2s = np.array([])
mH2 = np.array([])
H2s = np.array([])
nHs = np.array([])
nOs1 = np.array([])#([[],[],[]])
nOs2 = np.array([])
nOs3 = np.array([])
nCs = np.array([])
Rs1 = np.array([])#([[],[],[]])
Rs2 = np.array([])
Rs3 = np.array([])
losstimes = np.array([])
dataAs = np.array([])
dataMs = np.array([])
maxMenvs1 = np.array([])#([[],[],[]])
maxMenvs2 = np.array([])
maxMenvs3 = np.array([])
tenv901 = np.array([])#([[],[],[]])
tenv902 = np.array([])
tenv903 = np.array([])
bugsslow_M = np.array([])
bugsslow_apla = np.array([])
bugsmedium_M = np.array([])
bugsmedium_apla = np.array([])
bugsfast_M = np.array([])
bugsfast_apla = np.array([])
P1000i = np.array([])
P10001 = np.array([])
P10002 = np.array([])
P10003 = np.array([])
ihome = 0
ifolder = 0

for home in homes:
    folders = np.loadtxt('/groups/astro/piiamt/idl/simlist_log'+fastslowmedium[ihome]+'.txt', dtype='str')
    ### DONT DO BAD SIMULATIONS
    Aplas = np.array([])
    Mplas = np.array([])
    nOs = np.array([])
    nCs = np.array([])
    nHs = np.array([])
    R = np.array([])
    maxMenvs = np.array([])
    tenv90 = np.array([])
    Mplas_bug = np.array([])
    aplas_bug = np.array([])
    for folder in folders:
        timefile =home+folder+ '/time_series.idl'
        if not os.path.exists(timefile):
            folders = folders[np.where(folders!=folder)]
### LOOPING THROUGH ALL THE FILES
    for f in folders[:]:
        if f[3]!='A':	# if Apla has 2 decimal places like 1.25AU
            Apla = float(f[:4])
        elif f[3]=='A':	# aka if Apla has 1 decimal place like 1.0AU
            Apla = float(f[:3])
        Mpla = np.load(home + f + '/pysave'+'/mpla.npy')
        A = Apla
        M = np.round(Mpla/mearth, 2)
        #dataAs = np.append(np.where(A==As), A)
        #dataMs = np.append(np.where(M==Ms), M)
        Matm = np.load(home + f + '/pysave'+'/matm.npy')
        Menv = np.load(home + f + '/pysave'+'/menv.npy')
        Psur = np.load(home + f + '/pysave'+'/psur.npy')
        Tsur = np.load(home + f + '/pysave'+'/tsur.npy')
        ts = np.load(home + f + '/pysave'+'/t.npy')
        Msolid = np.load(home + f + '/pysave'+'/msolid.npy')
        MOBP = np.load(home + f + '/pysave/MOBP.npy')
        IW = np.load(home + f + '/pysave/deltaIW.npy')
        tmaxenv = ts[np.where(Menv==max(Menv))]
        ### IS SIMULATION BUGGY?
        if (ts[-1]/Myr_to_s)<300:
            Mplas_bug = np.append(Mplas_bug, Mpla[-1])
            aplas_bug = np.append(aplas_bug, Apla)
        #ienv90 = np.where((ts>tmaxenv)&(Menv<(0.1*max(Menv))))[0][0]
        ##ienv90 = np.where((ts>tmaxenv)&(Menv==0))[0][0]
        ienv90 = [-1]#np.where(ts>(290*Myr_to_s))[0][0]#np.where((ts>tmaxenv)&(Menv==0))[0][0]
#        if (Menv[-1]>(max(Menv)*0.1)): #in case sim finished before env loss
#            tenv90 = np.append(tenv90, ts[-1])
#        else:
#            tenv90 = np.append(tenv90, ts[np.where((ts>tmaxenv)&(Menv<(0.1*max(Menv))))][0])
        matm_CO2 = np.load(home + f + '/pysave'+'/matm_co2.npy')[ienv90]#[-1]
        matm_O2 = np.load(home + f + '/pysave'+'/matm_o2.npy')[ienv90]#[-1]
        matm_H2O = np.load(home + f + '/pysave'+'/matm_h2o.npy')[ienv90]#[-1]
        matm_CO = np.load(home + f + '/pysave'+'/matm_co.npy')[ienv90]#[-1]
        matm_N2 = np.load(home + f + '/pysave'+'/matm_n2.npy')[ienv90]#[-1]
        matm_H2 = np.load(home + f + '/pysave'+'/matm_h2.npy')[ienv90]#[-1]
        matm_SiO = np.load(home + f + '/pysave'+'/matm_sio.npy')[ienv90]#[-1]
        #maxMenvs = np.append(maxMenvs, max(Menv))
        Matms = np.append(Matms, Matm[ienv90])
        #Mplas = np.append(Mplas, Mpla[-1])
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
        #mol_atm = (matm_CO2/muCO2 + matm_CO/muCO + matm_N2/muN2 + matm_H2O/muH2O +
        #           matm_O2/muO2 + matm_H2/muH2 + matm_SiO/muSiO)
        mol_atm = (matm_CO2/muCO2 + matm_CO/muCO + matm_H2O/muH2O +
                   matm_O2/muO2 + matm_H2/muH2)
        CO2s = np.append(CO2s, mol_CO2/mol_atm)
        COs = np.append(COs, mol_CO/mol_atm)
        N2s = np.append(N2s, mol_N2/mol_atm)
        H2Os = np.append(H2Os, mol_H2O/mol_atm)
        O2s = np.append(O2s, mol_O2/mol_atm)
        H2s = np.append(H2s, mol_H2/mol_atm)
        mCO2 = np.append(mCO2, matm_CO2/meatm)
        mCO = np.append(mCO, matm_CO/meatm)
        mO2 = np.append(mO2, matm_O2/meatm)
        mH2O = np.append(mH2O, matm_H2O/meatm)
        mN2 = np.append(mN2, matm_N2/meatm)
        mH2 = np.append(mH2, matm_H2/meatm)
        if Psur[-1]>1000:
            P1000i = np.append(P1000i, int(ifolder))
        Aplas = np.append(Aplas, Apla)
        Mplas = np.append(Mplas, Mpla[-1])
        maxMenvs = np.append(maxMenvs, max(Menv))
        nHs = np.append(nHs, 2*mol_H2*nA + 2*mol_H2O*nA)
        nOs = np.append(nOs, 2*mol_O2*nA + mol_H2O*nA + mol_CO*nA + 2*mol_CO2*nA)
        nCs = np.append(nCs, mol_CO*nA + mol_CO2*nA)
        if (Menv[-1]>(max(Menv)*0.1)): #in case sim finished before env loss
            tenv90 = np.append(tenv90, ts[-1])
        else:
            tenv90 = np.append(tenv90, ts[np.where((ts>tmaxenv)&(Menv<(0.1*max(Menv))))][0])
        if Menv[-1]==0:
            losstime = ts[np.where((Msolid>0)&(ts>(5*Myr_to_s)))][0]
            losstimes = np.append(losstimes, losstime)
        ifolder = ifolder + 1
    lostenv = np.where(Menvs==0)
    envstay = np.where(Menvs>0)
    IW = IWs[lostenv]
    R = 2*nCs + nHs/2
    if ihome==0:
        Mplas1=Mplas
        Aplas1=Aplas
        nOs1 = nOs
        maxMenvs1 = maxMenvs
        tenv90f1 = tenv90
        Rs1 = R
        bugsslow_M = Mplas_bug
        bugsslow_apla = aplas_bug
        P10001 = P1000i.astype(int)
    elif ihome==1:
        Mplas2=Mplas
        Aplas2=Aplas
        nOs2 = nOs
        maxMenvs2 = maxMenvs
        tenv90f2 = tenv90
        Rs2 = R
        bugsmedium_M = Mplas_bug
        bugsmedium_apla = aplas_bug
        P10002 = P1000i.astype(int)
    elif ihome==2:
        Mplas3=Mplas
        Aplas3=Aplas
        nOs3 = nOs
        maxMenvs3 = maxMenvs
        tenv90f3 = tenv90
        Rs3 = R
        bugsfast_M = Mplas_bug
        bugsfast_apla = aplas_bug
        P10003 = P1000i.astype(int)

    ihome = ihome + 1
    ifolder = 0
    P1000i = np.array([])

#%%
cmaps = np.array(['turbo', 'YlOrBr', 'Greens', 'Blues', 
		'inferno', 'Reds', 'Greys', 'Purples'])
cmaps = np.array(['turbo', 'turbo', 'turbo','turbo','turbo','turbo','turbo','turbo'])

#labels = np.array(['Surface pressure (bar)', r'X_CO$_2$', r'X_O$_2$', r'X_H$_2$O',
#		'Atmosphere mass (kg)', 'Surface temperature (K)', 'X_CO', r'X_H$_2$'])

mm=1/25.4
fig, ax = plt.subplots(2,3, sharex=True, sharey=True, figsize=(190*mm,120*mm), gridspec_kw=dict(hspace=0.24, wspace=0.05))
ax[1,1].set_xlabel('Orbital radius (AU)')
ax[0,0].set_ylabel(r'Planet mass (M$_{\oplus}$)')
ax[1,0].set_ylabel(r'Planet mass (M$_{\oplus}$)')
### PREPARING ARRAYS FOR CONTOUR PLOTS
Acontour,Mcontour = np.meshgrid(Aplas1, Mplas1) #basically tile(As,(len(Ms),1)) and 
from matplotlib import colors
ls = np.linspace(-4,3, 100)#np.logspace(-4, 3, 100)
#print(maxMenvs1, Aplas1,Mplas1)
#import matplotlib.tri as tri
#triang = tri.TriContourSet(Aplas1, Mplas1/mearth)
### FOLLOWING IS THE MAX ENV MASSES LOG PLOT
at1 = ax[0,0].tricontourf(Aplas1, Mplas1/mearth, 
                    maxMenvs1/Mplas1, #levels=100,
                    levels=np.logspace(-7,-3.9,100),#(-4,3.205,100),
                    norm=colors.LogNorm(), cmap='turbo')
at2 = ax[0,1].tricontourf(Aplas2, Mplas2/mearth, 
                    maxMenvs2/Mplas2, #levels=100,
                    levels=np.logspace(-7,-3.9,100),#(-4,3.205,100),
                    norm=colors.LogNorm(), cmap='turbo')
at3 = ax[0,2].tricontourf(Aplas3, Mplas3/mearth, 
                    maxMenvs3/Mplas3, #levels=100,
                    levels=np.logspace(-7,-3.9,100),#(-4,3.205,100),
                    norm=colors.LogNorm(), cmap='turbo')
print(min(maxMenvs1/Mplas1), max(maxMenvs1/Mplas1))
ax[0,0].set_title('Slow', pad=30)
ax[0,1].set_title('Medium', pad=30)
ax[0,2].set_title('Fast', pad=30)
levellist1 = np.array([5.1,10,20,100,300])
levellist2 = np.array([5.1,10,50,100])
levellist3 = np.array([5.1,10,30,60])
#levellist2 = np.array([8,12, 16, 25, 40, 60, 80, 100,110])#MEDIUM
#levellist1 = np.array([8, 20, 40, 70, 100, 150, 200, 250, 300]) #SLOW
#levellist3 = np.array([8, 12, 16, 20, 25, 30, 40, 50, 60])
#manual_locations1 = [(1.0,1.8),(1.1,2.5),(1.15,2.8),(1.2,3.0),(1.25,3.4),#SLOW
#                     (1.3,3.5),(1.35,3.8),(1.4,4.0),(1.6,4.7)]
#manual_locations2 = [(1.0,1.0), (1.05,2.6), (1.15,3.1), (1.2,3.6), #MEDIUM
#                    (1.3,4.0), (1.32,4.3),(1.36,4.4), (1.4,4.7), (1.5,4.95)]
#manual_locations3 = [(1.0,2.0), (1.1,2.8), (1.18,3.2), (1.22,3.5),#FAST
#                    (1.3,3.9), (1.38,4.3), (1.45,4.6), (1.54,4.7), (1.65,4.95)]

CS1 = ax[0,0].tricontour(Aplas1, Mplas1/mearth, 
              tenv90f1/Myr_to_s,
              levels=levellist1, linewidths=0.7, colors='w')
CS2 = ax[0,1].tricontour(Aplas2, Mplas2/mearth, 
              tenv90f2/Myr_to_s,
              levels=levellist2, linewidths=0.7, colors='w')
CS3 = ax[0,2].tricontour(Aplas3, Mplas3/mearth, 
              tenv90f3/Myr_to_s,
              levels=levellist3, linewidths=0.7, colors='w')
ax[0,0].clabel(CS1, fontsize=7,inline=True)#, manual=manual_locations1, inline=True)
ax[0,1].clabel(CS2, fontsize=7,inline=True)#, manual=manual_locations2, inline=True)
ax[0,2].clabel(CS3, fontsize=7,inline=True)#, manual=manual_locations3, inline=True)
### FOLLOWING IS THE QUANTIFYING ATOM NUMBERS PLOT
at4 = ax[1,0].tricontourf(Aplas1[P10001], Mplas1[P10001]/mearth,
                    Rs1[P10001]/nOs1[P10001], levels=np.linspace(0,2.3,100), cmap='turbo')
at5 = ax[1,1].tricontourf(Aplas2[P10002], Mplas2[P10002]/mearth,
                    Rs2[P10002]/nOs2[P10002], levels=np.linspace(0,2.3,100), cmap='turbo')
at6 = ax[1,2].tricontourf(Aplas3[P10003], Mplas3[P10003]/mearth,
                    Rs3[P10003]/nOs3[P10003], levels=np.linspace(0,2.3,100), cmap='turbo')

ax[1,0].tricontour(Aplas1[P10001], Mplas1[P10001]/mearth,
                Rs1[P10001]/nOs1[P10001], levels=[1.3],
                linewidths=1, colors='w')
ax[1,1].tricontour(Aplas2[P10002], Mplas2[P10002]/mearth,
                Rs2[P10002]/nOs2[P10002], levels=[1.3],
                linewidths=1, colors='w')
ax[1,2].tricontour(Aplas3[P10003], Mplas3[P10003]/mearth,
                Rs3[P10003]/nOs3[P10003], levels=[1.3],
                linewidths=1, colors='w')

ax[1,0].scatter(1,1, marker=r'$\oplus$', s=50, color='k')
ax[1,0].scatter(1.524, 0.107446849, marker=u'$\u2642$', s=50, color='k')#mars
ax[1,0].scatter(0.72333, 0.815, marker=u'$\u2640$', s=50, color='k')#venus
ax[1,1].scatter(1,1, marker=r'$\oplus$', s=50, color='k')
ax[1,1].scatter(1.524, 0.107446849, marker=u'$\u2642$', s=50, color='k')#mars
ax[1,1].scatter(0.72333, 0.815, marker=u'$\u2640$', s=50, color='k')#venus
ax[1,2].scatter(1,1, marker=r'$\oplus$', s=50, color='k')
ax[1,2].scatter(1.524, 0.107446849, marker=u'$\u2642$', s=50, color='k')#mars
ax[1,2].scatter(0.72333, 0.815, marker=u'$\u2640$', s=50, color='k')#venus

ax[1,2].text(0.9, 0.2, 'no remaining atmospheres')
### LOG Y AXES
ax[0,0].set_yscale("log")
ax[1,0].set_yscale("log")
ax[0,0].set_ylim(0.1,5)
ax[1,0].set_ylim(0.1,5)
### BUGGY SIMS AS X-ES
#ax[1,0].scatter(bugsslow_apla, bugsslow_M/mearth, color='k', marker='X')
#ax[1,1].scatter(bugsmedium_apla, bugsmedium_M/mearth, color='k', marker='X')
#ax[1,2].scatter(bugsfast_apla, bugsfast_M/mearth, color='k', marker='X')
# FOR SMOOTH PLOTS IN A PDF DO THIS AS WELL (SLOW):
ats = np.array([at1, at2, at3, at4, at5, at6])
for at in ats:
    for c in at.collections:
        c.set_edgecolor('face')
p1 = ax[0,0].get_position().get_points().flatten() #[[x0,y0],[x1,y1]]
p2 = ax[0,1].get_position().get_points().flatten() #[[x0,y0],[x1,y1]]
p3 = ax[0,2].get_position().get_points().flatten() #[[x0,y0],[x1,y1]]
p4 = ax[1,0].get_position().get_points().flatten() #[[x0,y0],[x1,y1]]
p5 = ax[1,1].get_position().get_points().flatten() #[[x0,y0],[x1,y1]]
p6 = ax[1,2].get_position().get_points().flatten() #[[x0,y0],[x1,y1]]
at1_cbar = fig.add_axes([p1[0], p1[3], p1[2]-p1[0], 0.015])
at2_cbar = fig.add_axes([p2[0], p2[3], p2[2]-p2[0], 0.015])
at3_cbar = fig.add_axes([p3[0], p3[3], p3[2]-p3[0], 0.015])
at4_cbar = fig.add_axes([p4[0], p4[3], p4[2]-p4[0], 0.015])
at5_cbar = fig.add_axes([p5[0], p5[3], p5[2]-p5[0], 0.015])
at6_cbar = fig.add_axes([p6[0], p6[3], p6[2]-p6[0], 0.015])
fig.colorbar(at1, cax=at1_cbar, orientation='horizontal', 
             ticks=([1e-7, 1e-6, 1e-5 ,1e-4]))#, format='%.2e')
at1_cbar.tick_params(labelsize=5)
fig.colorbar(at2, cax=at2_cbar, orientation='horizontal', 
             label=r'Maximum envelope mass (planet mass fraction)',
             ticks=([1e-7, 1e-6, 1e-5 ,1e-4]))
at2_cbar.tick_params(labelsize=5)
fig.colorbar(at3, cax=at3_cbar, orientation='horizontal', 
             ticks=([1e-7, 1e-6, 1e-5 ,1e-4]))#, format='%.2f')
at3_cbar.tick_params(labelsize=5)
fig.colorbar(at4, cax=at4_cbar, orientation='horizontal', 
             ticks=np.array([0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2]), format='%.1f')
at4_cbar.tick_params(labelsize=5)
fig.colorbar(at5, cax=at5_cbar, orientation='horizontal', 
             label=r'R/O',
             ticks=np.array([0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2]), format='%.1f')
at5_cbar.tick_params(labelsize=5)
fig.colorbar(at6, cax=at6_cbar, orientation='horizontal', 
             ticks=np.array([0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2]), format='%.1f')
at6_cbar.tick_params(labelsize=5)
at1_cbar.xaxis.set_label_position('top')
at1_cbar.xaxis.set_ticks_position('top')
at2_cbar.xaxis.set_label_position('top')
at2_cbar.xaxis.set_ticks_position('top')
at3_cbar.xaxis.set_label_position('top')
at3_cbar.xaxis.set_ticks_position('top')
at4_cbar.xaxis.set_label_position('top')
at4_cbar.xaxis.set_ticks_position('top')
at5_cbar.xaxis.set_label_position('top')
at5_cbar.xaxis.set_ticks_position('top')
at6_cbar.xaxis.set_label_position('top')
at6_cbar.xaxis.set_ticks_position('top')


plt.savefig(OGhome+'/python/plots/6contourfig_log.pdf', bbox_inches='tight')
plt.show()
