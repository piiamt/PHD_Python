'''
Calculates pebble isolation mass at different orbital radii
'''
import numpy as np
import matplotlib.pyplot as plt

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
AU = 1.49e11
kB = 1.38e-23
mu0 = 2.34*1.660539040e-27
G = 6.67e-11
Ms = 1.99e30
alpha = 0.01
ME = 5.972e24
yr = 365.25*24*60*60

def T0(r): # r in AU
    return(150*r**(-3/7)*(0.7)**(2/7)*1**(-1/7))

def rho0(r):
    cs0 = (kB*T0(r)/mu0)**(1/2)
    Omega0 = (G*Ms/(r*AU)**3)**(1/2)
    H0 = cs0/Omega0
    Mdd = 1e-9*Ms*yr
    Sigma0 = Mdd/(3*np.pi*alpha*cs0*H0)
    return(Sigma0/((2*np.pi)**1/2*H0))

def dPdr(r):
    #T = T0(r)
    #rho = rho0(r)
    #return(rho/mu0*(-3*150/7)*r**(-10/7)*(0.7)**(2/7)*1**(-1/7))
    Mdd = 1e-9*Ms*yr
    C1 = Mdd*G*1*(mu0)**(1/2)/(3*np.pi*0.01*(2*np.pi*150*kB)**(1/2))*0.7**(-1/7)*1**(1/14)*(1.496e11)**(3/14)
    return(-39/14)

def Hdivr(r): # r in AU
    '''from Ida et al. 2016 "formation of ..."'''
    return(0.021*1**(1/7)*1**(-4/7)*r**(2/7))

def ffit(r):
    dP = dPdr(r)
    Hr = Hdivr(r)
    return((Hr/0.05)**3*(0.34*(np.log10(0.001)/np.log10(0.01))**4+0.66)*(1-(dP+2.5)/6))

def lambdaf(r):
    ff = ffit(r)
    return(0.00476/ff)

def Misostarf(r):
    ff = ffit(r)
    return(25*ff)

def Pif(r, mpla): # mpla in earth masses pls
    l = lambdaf(r)
    Misostar = Misostarf(r)
    return(l*(mpla-Misostar))

def Mmisof(r, mpla):
    Pi = Pif(r, mpla)
    Misostar = Misostarf(r)
    l = lambdaf(r)
    return(Misostar + Pi*ME/l)

radii = np.array([0.7,0.8,0.9,1.0,1.1,1.2,1.3])
print('from 0.7 to 1.3 AU Miso in ME is ',Misostarf(radii))

manyradii = np.linspace(0.7,1.7, 100)
mm=1/25.4
plt.figure(figsize=(88*mm,90*mm))
plt.plot(manyradii, Misostarf(manyradii), color=colors[1])
plt.xlabel('Orbital radius (AU)')
plt.ylabel(r'Pebble isolation mass (M$_{\oplus}$)')
plt.savefig('plots/pebbleiso.pdf', bbox_inches='tight')

plt.show()
