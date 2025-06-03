'''
Created 11.08.2022, author @piiamt
LOOPS THROUGH A GIVEN ARRAY OF ORBITAL RADII AND PLANET MASSES AND 
SUBMITS THEIR SIMULATIONS TO AURORA.
THE CODE THAT WILL LOOP THROUGH SAID SIMULATIONS AND RUN IDL AND OVERVIEW
PLOTTING ON THEM WILL BE IN ANOTHER FILE
'''

import numpy as np
import sys
import os, os.path
from os.path import expanduser
home = '/lustre/astro/piiamt/locean_outgas/slow/'
import shutil
#from simlist import folders
#simulations = np.loadtxt('simlist.txt', dtype='str')
folders = np.array(os.listdir(home))
simulations = folders[np.where((folders!='overview')&(folders!='python')&(folders!='starter'))]

##### CONSTANTS #####
apla_min = 0.7 
apla_max = 1.8 
apla_inc = 0.1 #0.05
mpla_min = 0.5 #0.9
mpla_max = 5.5 #1.05
mpla_inc = 0.5

myr         = 365.25*24*60*60*1e6   	# seconds
mearth      = 5.9722e24             	# kg
msun        = 1.98847e30            	# kg
AU          = 1.496e11              	# m
R0          = 250*1000              	# meters
rho0        = 1.421573328e20 / ((4/3)*np.pi*R0**3) 	 # kg/m^3
M0          = rho0 * (4/3)*np.pi*R0**3 	# kg (M of initial planetesimal)
myr5        = 5*myr
####################

def taccfn(Mpla):
	'''
	returns the accretion timescale for an end planet mass given in mearth
	with the starting planetesimal of R=250 km
	'''
	Mtot = Mpla * mearth		# kg
	t = myr5 / np.log(Mtot/M0)	# s
	return(t)

Aplas = np.array([1.0])#np.arange(apla_min, apla_max, apla_inc) # total 21
#Mplas = np.arange(mpla_min, mpla_max, mpla_inc)
#Mplas = np.insert(Mplas, 1, np.array([0.1,0.2,0.3,0.4,0.6,0.7,0.8,0.9]))
#Aplas = np.array([1.0])
#Mplas = np.array([0.1, 0.15, 0.24, 0.37, 0.57,
#                  0.88, 1.36, 2.09, 3.24, 5.0])
#Mplas = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 
#                  1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0])
Mplas = np.array([0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0])
n = 0

for aplaAU in Aplas:
    for mpla in Mplas:
        simname=str(np.round(aplaAU,2))+'AU_'+ str(np.round(mpla,2)) +'Mearth'
        if simname in simulations:
           print('the simulation '+simname+' already exists, next one')
        else:
            os.mkdir(home + simname)
            shutil.copy(home+'/starter/slurm_start_LU.sh',
                    home + simname + '/slurm_start_LU.sh')
            tacc = taccfn(mpla)
            apla = aplaAU * AU
            print('we have created the new directory')
            ### now we have the input constants and need to create new input.in
            with open(home+'/starter/input.in') as ogin:
                lines = ogin.readlines()
            with open(home+simname+'/input.in', 'a') as f:
                for line in lines:
                    consts = line.split()
                    if consts[0][:5]=='tacc=':
                        consts[0] = 'tacc='+'{:e}'.format(tacc)+','
                    if consts[0][:5]=='apla=':
                        consts[0] = 'apla='+'{:e}'.format(apla)+','
#                    if consts[0][:5]=='tmax=': # REMOVE LATER, THIS ONE JUST
#                    # TESTS WITH A SHORTER SIMULATION
#                        consts[0] = 'tmax='+'{:e}'.format(2*myr5)+','
                    f.write(''.join(consts)+'\n')
            ### set up and make ADAP:
            print('setting up ADAP')
            os.environ['SHELL']='tcsh'
            #os.system('tcsh ADAP_setupsrc\n')
            os.chdir(home+simname)
            print('entered tcsh, now the cwd is',os.getcwd())
            os.system('ADAP_setupsrc\n')
            print('ran adap setup')
           #os.system('cp ../../1.0AU_1.0Mearth/src/adap.x src')
           #os.system('make\n')
            if n==0:
                os.system('make\n')
                xlocation = home + simname + '/src/adap.x'
                print('made adap OG')
            else:
                os.system('cp ' + xlocation + ' src')
                print('made adap (by copying adap.x from'+xlocation+')')
            ### now have to submit the new input file with sbatch somehow
            os.system('sbatch slurm_start_LU.sh')
            print('submitted '+simname)
            simulations = np.append(simulations, simname)
            n = n+1

### adding new simulation folder name to simlist.txt:
#np.savetxt('../idl/simlist.txt', simulations, fmt='%s', newline='\n')
