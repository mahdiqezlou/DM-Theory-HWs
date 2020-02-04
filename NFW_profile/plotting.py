import math
import matplotlib.pyplot as plt
import numpy as np
from unitsystem import UnitSystem as US
v_factor = 200
Kpc_to_cm = 3.086*10**21
G_const = 6.674*10**(-8.0)
# Sloar mass to g
M_sun = 1.989*10**33
halos = {'Mass':[10**10, 10**12, 10**14], 'c':[15.,10.,5.], 'label':['log(M) = 10, c = 15', 'log(M) = 12, c = 10', 'log(M) = 14, c = 5']}

us = US()


def plot_rho_to_rho_v(s = np.arange(0.02, 5., 0.02)):
    
    #halos = {'Mass':[10**10, 10**12, 10**14], 'c':[15.,10.,5.], 'label':['log(M) = 10, c = 15', 'log(M) = 12, c = 10', 'log(M) = 14, c = 5']}
    for i in range(3):
        x, y = rho_to_rho_v(M_vir=halos['Mass'][i], c=halos['c'][i], s=s)
        plt.loglog(x, y, label=halos['label'][i])
        
    plt.xlabel(r'$r/r_v$')
    plt.ylabel(r'$\rho / \rho_{vir}$')
    plt.legend()

def plot_mass_to_mass_v(s = np.arange(0.02, 5., 0.02)):
    
    #halos = {'Mass':[10**10, 10**12, 10**14], 'c':[15.,10.,5.], 'label':['log(M) = 10, c = 15', 'log(M) = 12, c = 10', 'log(M) = 14, c = 5']}
    for i in range(3):
        x, y = mass_to_mass_v(c=halos['c'][i], s=s)
        plt.plot(x, y, label=halos['label'][i])
        
    plt.xlabel(r'$r/r_v$')
    plt.ylabel(r'$M / M_{vir}$')
    plt.legend()

def plot_v_to_v_v(s = np.arange(0.02, 5., 0.02)):

    for i in range(3):
        x, y = v_to_v_v(c= halos['c'][i], s=s)
        plt.plot(x, y, label=halos['label'][i])

    plt.xlabel(r'$r/r_v$')
    plt.ylabel(r'$v_c / v_{c, vir}$')
    plt.legend()

def plot_rho(s = np.arange(0.02, 5., 0.02)):

    for  i in range(3):
        y = rho_to_rhoc(c=halos['c'][i], s=s)
        r_vir = get_r_vir(M_vir = halos['Mass'][i], c=halos['c'][i])
        x = s*r_vir/(Kpc_to_cm)

        y = y*rho_to_rhoc(c=halos['c'][i], s=s)*rho_crit()

        plt.loglog(x, y, label = halos['label'][i])

    plt.xlabel(r'$ r \ (\ Kpc \ )$')
    plt.ylabel(r'$\rho \ (\ g/cm^{-3}\ ) $')



def plot_mass(s = np.arange(0.02, 5., 0.02)):
    """ Plot mass profile in solar Mass vs Kpc """
    for i in range(3):
        x, y = mass_to_mass_v(c=halos['c'][i], s=s)
        r_vir = get_r_vir(M_vir = halos['Mass'][i], c=halos['c'][i])
        x = x*r_vir/(Kpc_to_cm)
        
        y = y*halos['Mass'][i]
        plt.semilogy(x, y, label=halos['label'][i])

    plt.xlabel(r'$r \ ( \ Kpc \ )$')
    plt.ylabel(r'$ M / M_{\odot} \ $')





def plot_v(s = np.arange(0.02, 5., 0.02)):
    """ plot circular velocity in physical units """
    for i in range(3):
        x , y = v_to_v_v(c = halos['c'][i], s=s)
        r_vir = get_r_vir(M_vir = halos['Mass'][i], c=halos['c'][i])
        x = x*r_vir/(Kpc_to_cm)
        ## velocity in km/s
        y = y*get_v_vir(M_vir=halos['Mass'][i], r_vir=r_vir)/(10**5)
        plt.semilogy(x, y, label=halos['label'][i])
    plt.xlabel(r'$r \ (\ kpc \ )$')
    plt.ylabel(r'$V_{circ} \ (\ km/s \ )$')
    plt.legend()



def rho_to_rho_v(M_vir, c, s= np.arange(0.02, 5., 0.02)):
    
    return (s, rho_to_rhoc(c, s)/rho_to_rhoc(c, 1))

def rho_to_rhoc(c, s):

    return (v_factor*g_func(c)*c**2)/(3*s*(1 + c*s)**2)


def rho_phys_units(c, s):

    return rho_crit()*rho_to_rhoc(c, s)

def g_func(c):
    g = 1./(np.log(1 + c) - c/(1 + c))
    return g

def rho_crit():
    return us.rho_crit(us.hubble(z=0., omegam0 = 0.3))


def get_r_vir(M_vir, c):
    """ virial radius in cm """
    r_vir = (3.*M_vir*M_sun/(4*math.pi*v_factor*rho_crit()))**(1/3)
    return r_vir

def get_v_vir(M_vir, r_vir):
    """ virial velocity in cm/s"""
    return np.sqrt(G_const*M_vir*M_sun/r_vir)

def mass_to_mass_v(c, s):
    # Eq 8
    return (s, g_func(c)*(np.log(1 + c*s) - (c*s/(1 + c*s))))

def v_to_v_v(c, s):
    """ Circulcar velocoty to Circular velocity at virial radius """

    return (s, (g_func(c)/s)*(np.log(1 + c*s) - (c*s/(1 + c*s))))
