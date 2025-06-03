'''
Makes txt file of all simulation folder names in given folder
'''
import numpy as np
import sys
import os, os.path
from os.path import expanduser
#home = '/projects/astro3/nobackup/piia'
home = '/lustre/astro/piiamt/delayed_XUV/'#noatmloss_dense/'
import shutil
folders = np.array(os.listdir(home))
simulations = folders[np.where((folders!='tester')&(folders!='overview')&(folders!='python')&(folders!='starter'))]
np.savetxt('../idl/simlist_delayed_XUV.txt', simulations, fmt='%s', newline='\n')
