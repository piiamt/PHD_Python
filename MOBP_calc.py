'''
calculates the magma ocean base pressures of all simulations and 
saves them into a table in the /projects/.../piia/overview/ dir.

'''
from scipy.io import readsav
import numpy as np
import itertools
# import matplotlib
import matplotlib.pyplot as plt
# from matplotlib.colors import LogNorm
import os, os.path
from os.path import expanduser
OGhome = expanduser('~')
#home = '/projects/astro3/nobackup/piia/'
home = '/lunarc/nobackup/projects/astro4/piia/'
#from simlist import folders
folders = np.loadtxt('../idl/simlist.txt', dtype='str')
######################


plt.rcParams.update({'font.size'           : 13, 
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
#rcor = 4579023.0  	#
#rpla = 7104767.0 	#
#mcore = 2.2775820e24	#
G = 6.674e-11
rhocore = 5663.2660	
rhomagma = 3400.0

MOBPs = np.array([])

for sim in folders:
	print('Currently at folder ',sim)
	timefile = home + sim + '/time_series.idl'
	timef = readsav(timefile)
	timekey = list(timef.values())[0].dtype
	tser = list(timef.values())[0][0]
	### getting the needed constants from the idl time files
	rcor = tser.rcor[-1]
	rpla = tser.rpla[-1]
	mcore = tser.mcore[-1]
	psur = tser.psur[-1]
	g = G*mcore/(rcor**2)
	h = rpla - rcor
	P = rhomagma * g * h
	MOBPs = np.append(MOBPs, P+psur)

print('We have now gone through all the folders, saving')
np.save(home+'/overview/MOBPs', MOBPs)


